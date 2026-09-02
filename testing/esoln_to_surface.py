#!/usr/bin/env python
"""
esoln_to_surface.py -- maps a ModEM binary .esoln E-field solution (e.g. the
TOTAL field written by ModEM's -E / secondary field formulation job) to a
surface-only ASCII file, in the SAME 14-column format FWD1D.f90's own
OUTPUT_SURFACE feature writes (and that plot_surface_output.py already knows
how to plot) -- so an existing, already-validated tool can be reused directly
for bulk plotting.

Uses read_esoln.py's reader for the file format (header, dx/dy/dz, per-
(period,mode) Ex/Ey/Ez blocks) -- see that module's docstring for the exact
binary layout, traced directly from ModEM's ioAscii.f90.

WHAT IS ACTUALLY COMPUTED
--------------------------
Only the HORIZONTAL E field (Ex, Ey) at the air/earth interface (the z-node
index = nzAir, i.e. the boundary between the air layers and the first earth
layer) is extracted. H is NOT available from a ModEM -E .esoln (that file
format only ever stores E -- see read_esoln.py's docstring) and the VERTICAL
E component (Ez) is skipped: Ez lives on z CELL centres, not the z NODE
exactly at the interface, so reporting it here would need an extra
across-the-interface interpolation that isn't needed for this horizontal-
field application (MT/GIC-style surface field maps).

To let plot_surface_output.py be reused UNMODIFIED, the output file keeps
that plotter's full 14-column layout (lat lon + Re/Im of Hr,Htheta,Hphi,
Er,Etheta,Ephi) -- but Hr,Htheta,Hphi and Er are ALWAYS written as 0+0j,
clearly labelled as unavailable/uncomputed in the header line, NOT as a
physical zero. Only Etheta and Ephi carry real data. (The resulting
*.surface.H.png from plot_surface_output.py is therefore a placeholder and
should be ignored/discarded when reviewing plots from this pipeline.)

SIGN CONVENTION (Ex,Ey -> Etheta,Ephi)
----------------------------------------
ModEM's grid axes are x=latitude (increasing NORTHWARD, confirmed directly:
read_esoln.py's ox = "south-edge latitude", so cumsum(dx) from ox increases
northward) and y=longitude (increasing EASTWARD). global1d's own native
components are Etheta (colatitude-increasing = SOUTHWARD-positive) and Ephi
(longitude-increasing = EASTWARD-positive, same sense as ModEM's y).

This sign relationship is not re-derived from scratch here -- it follows
directly from FWD1D.f90's own, already-tested MODEM-format writer
(write_modem_efield_block/write_modem_header): that writer requires
target_conv%theta_convention=LAT_S2N, under which apply_output_convention
BOTH reverses the theta array index AND NEGATES the theta component
(E%y = -E%y, see output_convention.f90) before the (sign-preserving, pure
axis-relabeling) transpose into ModEM's Ex slot. I.e. by construction,
ModEM's Ex (northward-positive) = -1 * (global1d's native Etheta,
southward-positive); ModEM's Ey (eastward-positive) = +1 * global1d's Ephi
(no sign flip in that leg of the writer). Reading a REAL ModEM esoln (as
this script does) is exactly the inverse mapping:
    Etheta = -Ex   Ephi = +Ey
Only used here to relabel components; the ACTUAL E values are ModEM's own
(the total, primary+secondary field), not converted between codes' unit
conventions in any other way (no additional physics rescale is applied).

GEOMETRY
--------
lat/lon come directly from the esoln's own header (ox,oy,dx,dy) -- no
separate .grd file is needed or read. Ex sits at (latitude CELL centres,
longitude NODES); Ey sits at (latitude NODES, longitude CELL centres)
(standard edge/Yee staggering -- see read_esoln.py's shape comments). Both
are averaged onto the common (latitude CELL centres, longitude CELL
centres) grid by a simple midpoint average of their two bracketing NODE
values along the mismatched axis -- exact, not approximate, since a cell
centre is by definition the midpoint of its two bounding nodes.

Usage:
  python esoln_to_surface.py --esoln FILE.esoln --out OUT_STEM
Writes one OUT_STEM.T<iFreq2>.mode<iMode2>.surface file per (period,mode)
block in the input (almost always exactly one for this pipeline's single-
mode, single-period sources).
"""
import argparse
import numpy as np

from read_esoln import read_esoln


def _cell_centers_from_origin(origin, widths):
    nodes = origin + np.concatenate(([0.0], np.cumsum(widths)))
    centers = 0.5 * (nodes[:-1] + nodes[1:])
    return nodes, centers


def surface_from_block(d, block):
    """Return (lat_centers, lon_centers, Etheta[nx,ny], Ephi[nx,ny])."""
    nzAir = d['nzAir']
    lat_nodes, lat_centers = _cell_centers_from_origin(d['ox'], d['dx'])
    lon_nodes, lon_centers = _cell_centers_from_origin(d['oy'], d['dy'])

    Ex = block['Ex'][:, :, nzAir]   # (nx, ny+1) -- lat centres, lon nodes
    Ey = block['Ey'][:, :, nzAir]   # (nx+1, ny) -- lat nodes,   lon centres

    Ex_c = 0.5 * (Ex[:, :-1] + Ex[:, 1:])   # -> (nx, ny) at lon centres
    Ey_c = 0.5 * (Ey[:-1, :] + Ey[1:, :])   # -> (nx, ny) at lat centres

    Etheta = -Ex_c
    Ephi = Ey_c
    return lat_centers, lon_centers, Etheta, Ephi


def write_surface_file(fname, esoln_fname, lat_centers, lon_centers, Etheta, Ephi, period, iMode):
    nx, ny = Etheta.shape
    zero = 0.0 + 0.0j
    with open(fname, 'w') as f:
        f.write(f"# esoln_to_surface.py: ModEM -E TOTAL E-field surface slice "
                 f"(air/earth interface) for period {period:.3f} secs, mode {iMode:02d}. "
                 f"H NOT AVAILABLE (fixed 0+0j); Er NOT COMPUTED (fixed 0+0j). "
                 f"SOURCE_ESOLN={esoln_fname}\n")
        f.write("# lat_deg lon_deg  Re(Hr) Im(Hr) Re(Hth) Im(Hth) Re(Hph) Im(Hph)"
                 "  Re(Er) Im(Er) Re(Eth) Im(Eth) Re(Eph) Im(Eph)\n")
        for i in range(nx):
            lat = lat_centers[i]
            for j in range(ny):
                lon = lon_centers[j]
                eth = Etheta[i, j]
                eph = Ephi[i, j]
                # column order: Hr,Hth,Hph,Er (all zero), then Eth,Eph (real)
                cols = []
                for v in [zero, zero, zero, zero]:
                    cols.append(v.real); cols.append(v.imag)
                cols.append(eth.real); cols.append(eth.imag)
                cols.append(eph.real); cols.append(eph.imag)
                f.write(f"{lat:11.4f}{lon:11.4f}" + "".join(f"{c:16.7e}" for c in cols) + "\n")
    print(f"wrote {fname}  ({nx}x{ny} = {nx*ny} surface points)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--esoln', required=True, help='ModEM binary .esoln file (e.g. -E output)')
    ap.add_argument('--out', required=True, help='output file stem (per-block suffix added)')
    args = ap.parse_args()

    d = read_esoln(args.esoln)
    print(f"{args.esoln}: nPer={d['nPer']} nMode={d['nMode']} nx={d['nx']} ny={d['ny']} "
          f"nz={d['nz']} nzAir={d['nzAir']}")

    for b in d['blocks']:
        lat_c, lon_c, Etheta, Ephi = surface_from_block(d, b)
        outname = f"{args.out}.T{b['iFreq']:02d}.mode{b['iMode']:02d}.surface"
        write_surface_file(outname, args.esoln, lat_c, lon_c, Etheta, Ephi,
                            b['period'], b['iMode'])


if __name__ == '__main__':
    main()
