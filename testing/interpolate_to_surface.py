#!/usr/bin/env python
"""
interpolate_to_surface.py -- INDEPENDENT cross-check of FWD1D.f90's
OUTPUT_SURFACE feature (subroutines surfaceField1d / surfaceField1d_s2).

Reads a global1d GLOBAL-format volume .efield and .hfield (the staggered
E/H cvectors written by write_cvector) plus the .grd grid file, then
INTERPOLATES every field component to the Earth's surface (r = R_earth)
at the lateral CELL CENTRES of the grid, and writes a *.surface file in
exactly the same format as FWD1D.f90's OUTPUT_SURFACE. Comparing this
file against the one FWD1D wrote directly validates that routine (the two
should agree to within lateral/radial interpolation error, since FWD1D
evaluates the analytic field exactly at the surface while this script
linearly interpolates the volume samples).

Any scalar normalization already baked into the volume (e.g. the ModEM
source normalization) is inherited automatically -- we interpolate the
STORED values, so the output is in whatever convention the volume is in.

CONVENTION: this script assumes the volume was written in the solver's
NATIVE convention (theta=COLAT_N2S, r=R_UP -- e.g. OUTPUT_CONVENTION=
'SUNEGBERT2012'), i.e. WITHOUT the theta/r array-index reversal that
apply_output_convention performs for LAT_S2N/R_DOWN targets. It reads the
convention metadata from the field-file header and refuses to run on a
reversed-convention file (the reversal + mid/node off-by-one bookkeeping
is deliberately out of scope here; regenerate the volume with
OUTPUT_CONVENTION='SUNEGBERT2012' for this cross-check).

Reuses read_field / read_grid_file / recenter_fake_pole / _KIND from
plot_global1d_output.py so the file/grid parsing and per-component
staggering are shared, not re-derived.

Usage:
  python interpolate_to_surface.py --efield E.efield --hfield H.hfield \
      --grid grid.grd --out out.surface \
      [--fake-center-lat 0 --fake-center-lon 90] [--R-earth 6371.0]
"""
import argparse
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from plot_global1d_output import read_field, read_grid_file, recenter_fake_pole, _KIND

R_EARTH_KM_DEFAULT = 6371.0   # matches earth%r0 = 6371.0e3 in FWD1D.f90


def _axis_positions(kind, node_vals):
    """Physical positions for one axis given its 'node'/'mid' staggering.
    node -> the node array itself (len N+1); mid -> cell centres (len N)."""
    if kind == 'mid':
        return 0.5 * (node_vals[:-1] + node_vals[1:])
    return node_vals


def _axis_count(kind, ncells):
    """Number of VALID samples along an axis (drops write_cvector's trailing
    zero-padding for 'mid' quantities): mid -> ncells, node -> ncells+1."""
    return ncells if kind == 'mid' else ncells + 1


def interp_component(arr, kinds, ph_deg, th_deg, r_m, tgt_ph, tgt_th, r_earth_m):
    """Interpolate one staggered component array to (tgt_ph[i], tgt_th[j],
    r_earth_m) for all lateral cell centres. Returns complex (nphi, ntheta).

    arr        : (nx+1, ny+1, nz+1) raw component (axes phi, theta, r)
    kinds      : (phi_kind, theta_kind, r_kind), each 'node' or 'mid'
    ph_deg     : phi node positions, deg (nx+1)
    th_deg     : theta node positions, deg (ny+1)  [colatitude]
    r_m        : r node positions, m (nz+1)  [top->bottom, decreasing]
    tgt_ph/th  : target cell-centre phi/theta, deg (nx / ny)
    """
    pk, tk, rk = kinds
    nx = len(ph_deg) - 1
    ny = len(th_deg) - 1
    nz = len(r_m) - 1

    phi_pos = _axis_positions(pk, ph_deg)
    th_pos = _axis_positions(tk, th_deg)
    r_pos = _axis_positions(rk, r_m)

    data = arr[:_axis_count(pk, nx), :_axis_count(tk, ny), :_axis_count(rk, nz)]

    # RegularGridInterpolator needs strictly increasing axes. phi and theta
    # (colatitude) are already increasing as read; r is decreasing (top->
    # bottom) so reverse both the r axis and the data along it.
    r_pos_inc = r_pos[::-1]
    data_inc = data[:, :, ::-1]

    # interpolate real and imaginary parts separately (robust across scipy
    # versions); allow tiny out-of-range at exact boundaries via extrapolation.
    def _rgi(vals):
        return RegularGridInterpolator((phi_pos, th_pos, r_pos_inc), vals,
                                       method='linear', bounds_error=False,
                                       fill_value=None)

    fre = _rgi(data_inc.real)
    fim = _rgi(data_inc.imag)

    PH, TH = np.meshgrid(tgt_ph, tgt_th, indexing='ij')   # (nx, ny)
    RR = np.full(PH.shape, r_earth_m)
    pts = np.stack([PH.ravel(), TH.ravel(), RR.ravel()], axis=-1)
    out = fre(pts) + 1j * fim(pts)
    return out.reshape(PH.shape)


def _check_native(conv, which):
    """Refuse reversed conventions (see module docstring)."""
    if conv is None:
        print(f"  ({which}: no convention metadata in header -- assuming NATIVE)")
        return
    theta = conv.get('theta')
    r = conv.get('r')
    if (theta is not None and theta != 'COLAT_N2S') or (r is not None and r != 'R_UP'):
        raise SystemExit(
            f"ERROR: {which} field is in a REVERSED convention "
            f"(theta={theta}, r={r}); this script only handles NATIVE "
            f"(theta=COLAT_N2S, r=R_UP). Regenerate the volume with "
            f"OUTPUT_CONVENTION='SUNEGBERT2012'.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--efield', required=True, help='GLOBAL-format volume .efield')
    ap.add_argument('--hfield', required=True, help='GLOBAL-format volume .hfield')
    ap.add_argument('--grid', required=True, help='.grd grid file')
    ap.add_argument('--out', required=True, help='output .surface file')
    ap.add_argument('--fake-center-lat', type=float, default=None)
    ap.add_argument('--fake-center-lon', type=float, default=None)
    ap.add_argument('--R-earth', type=float, default=R_EARTH_KM_DEFAULT,
                    help='surface radius in km (default 6371)')
    args = ap.parse_args()

    E = read_field(args.efield)
    H = read_field(args.hfield)
    _check_native(E['convention'], 'efield')
    _check_native(H['convention'], 'hfield')

    g = read_grid_file(args.grid)
    th_deg = g['th_deg'].copy()
    ph_deg = g['ph_deg'].copy()
    r_m = g['r_km'] * 1.0e3
    if (args.fake_center_lat is not None) != (args.fake_center_lon is not None):
        raise SystemExit('--fake-center-lat and --fake-center-lon must be given together')
    if args.fake_center_lat is not None:
        th_deg, ph_deg = recenter_fake_pole(th_deg, ph_deg,
                                            args.fake_center_lat, args.fake_center_lon)

    nx, ny = g['nx'], g['ny']
    tgt_ph = 0.5 * (ph_deg[:-1] + ph_deg[1:])   # cell-centre longitude (nx)
    tgt_th = 0.5 * (th_deg[:-1] + th_deg[1:])    # cell-centre colatitude (ny)
    r_earth_m = args.R_earth * 1.0e3

    # component -> (source field dict, raw slot). E is EDGE, H is FACE; the
    # per-(gridType,slot) staggering comes from _KIND (shared with the plotter).
    def comp(fld, slot):
        return interp_component(fld[slot], _KIND[fld['gridType']][slot],
                                ph_deg, th_deg, r_m, tgt_ph, tgt_th, r_earth_m)

    Hr = comp(H, 'r'); Hth = comp(H, 'theta'); Hph = comp(H, 'phi')
    Er = comp(E, 'r'); Eth = comp(E, 'theta'); Eph = comp(E, 'phi')

    with open(args.out, 'w') as f:
        f.write("# interpolate_to_surface.py: SURFACE fields (r=R_earth, cell centres) "
                "interpolated from volume E/H\n")
        f.write("# lat_deg lon_deg  Re(Hr) Im(Hr) Re(Hth) Im(Hth) Re(Hph) Im(Hph)"
                "  Re(Er) Im(Er) Re(Eth) Im(Eth) Re(Eph) Im(Eph)\n")
        for j in range(ny):
            lat = 90.0 - tgt_th[j]
            for i in range(nx):
                lon = tgt_ph[i]
                vals = [Hr[i, j], Hth[i, j], Hph[i, j], Er[i, j], Eth[i, j], Eph[i, j]]
                cols = []
                for v in vals:
                    cols.append(v.real); cols.append(v.imag)
                f.write(f"{lat:11.4f}{lon:11.4f}" + "".join(f"{c:16.7e}" for c in cols) + "\n")
    print(f"wrote {args.out}  ({nx}x{ny} = {nx*ny} surface points)")


if __name__ == '__main__':
    main()
