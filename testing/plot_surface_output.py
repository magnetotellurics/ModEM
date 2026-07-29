"""
Plot global1d's *.surface output -- the surface-only counterpart of
plot_global1d_output.py.

A *.surface file (see FWD1D.f90's OUTPUT_SURFACE / surfaceField1d) holds the
full H and E field at the Earth's surface (r = R_earth), evaluated at the
lateral CELL CENTRES of the grid, with all six components CO-LOCATED and an
explicit (lat, lon) for every point. Because of that, this plotter needs
NONE of the volume plotter's machinery -- no --grid, no fake-pole
recentering, no mid/node staggering handling, no convention index-reversal,
no r-index selection. The file is already at the surface, at cell centres,
in whatever OUTPUT_CONVENTION was requested; we just reshape and plot.

File format (FWD1D.f90 OUTPUT_SURFACE):
  line 1: '# ... for period P secs, mode M. OUTPUT_CONVENTION=... SOLVER=... time=...'
  line 2: '# lat_deg lon_deg  Re(Hr) Im(Hr) Re(Hth) Im(Hth) Re(Hph) Im(Hph)
                              Re(Er) Im(Er) Re(Eth) Im(Eth) Re(Eph) Im(Eph)'
  data:   lat lon + 12 reals, ordered theta(lat) OUTER, phi(lon) INNER
          (ny*nx rows).

Produces, per input file, an E figure (Etheta, Ephi) and an H figure
(Hr, Htheta, Hphi): Re/Im column pairs, RdBu_r pcolormesh, north at top,
longitude on [-180, 180). Style matches plot_global1d_output.py.

CLI:
  python plot_surface_output.py FILE.surface [FILE2.surface ...]
Output PNGs are written next to each input file (see out_png_names()), so
runs never overwrite each other.
"""
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# data-column layout of a .surface row (0-indexed): lat, lon, then Re/Im of
# Hr, Htheta, Hphi, Er, Etheta, Eph -- keyed by our short component names.
_COLS = {
    'Hr':     (2, 3),
    'Htheta': (4, 5),
    'Hphi':   (6, 7),
    'Er':     (8, 9),
    'Etheta': (10, 11),
    'Ephi':   (12, 13),
}


def parse_header(fname):
    """Read the first (comment) line and pull out period/mode/convention/
    solver tokens for the figure titles. Returns a dict (values may be '?')."""
    with open(fname) as f:
        line = f.readline()
    info = {'period': '?', 'mode': '?', 'OUTPUT_CONVENTION': '?', 'SOLVER': '?', 'time': '?'}
    # 'for period 3.981 secs, mode 01.'
    if 'period' in line:
        try:
            info['period'] = line.split('period', 1)[1].split('secs', 1)[0].strip()
        except IndexError:
            pass
    if 'mode' in line:
        try:
            info['mode'] = line.split('mode', 1)[1].strip().split()[0].rstrip('.')
        except IndexError:
            pass
    for tok in line.split():
        if '=' in tok:
            k, _, v = tok.partition('=')
            if k in info:
                info[k] = v.rstrip('.')
    return info


def read_surface(fname):
    """Return (lon_1d [nx, ascending, signed], lat_1d [ny], comps dict of
    (ny, nx) complex arrays). The file's rows are theta(lat)-outer,
    phi(lon)-inner, so the first nx rows share lat[0]."""
    data = np.loadtxt(fname, comments='#')
    if data.ndim != 2 or data.shape[1] < 14:
        raise ValueError(f"{fname}: expected >=14 columns (lat lon + 12 field reals), "
                         f"got shape {data.shape}")
    lat_col, lon_col = data[:, 0], data[:, 1]
    n = len(lat_col)

    # nx = length of the first constant-lat block (rows are lat-outer/lon-inner)
    changes = np.nonzero(lat_col != lat_col[0])[0]
    nx = int(changes[0]) if len(changes) else n
    if n % nx != 0:
        raise ValueError(f"{fname}: {n} rows not divisible by inferred nx={nx}")
    ny = n // nx

    lat_1d = lat_col[::nx]          # one lat per block (ny), native order (N->S)
    lon_1d = lon_col[:nx]           # nx lon values (first block)

    # longitude -> signed [-180, 180) and sort columns so the axis is
    # monotonic even for a full-globe (0..360) or wrap-around grid.
    lon_signed = ((lon_1d + 180.0) % 360.0) - 180.0
    order = np.argsort(lon_signed)
    lon_sorted = lon_signed[order]

    comps = {}
    for name, (rc, ic) in _COLS.items():
        c = (data[:, rc] + 1j * data[:, ic]).reshape(ny, nx)   # row=lat, col=lon
        comps[name] = c[:, order]                              # reorder columns to match lon_sorted
    return lon_sorted, lat_1d, comps


def out_png_names(fname):
    stem = fname[:-len('.surface')] if fname.endswith('.surface') else fname
    return stem + '.surface.E.png', stem + '.surface.H.png'


def _plot(lon, lat, panels, suptitle, fname):
    """panels: list of (data_2d, title). All panels share the one (lon, lat)
    cell-centre grid (unlike the volume plotter, where staggered components
    each carry their own grid)."""
    n = len(panels)
    fig, axes = plt.subplots(n, 2, figsize=(10, 4 * n), squeeze=False)
    for row, (data, title) in enumerate(panels):
        for col, (part, ri) in enumerate([(data.real, "Re"), (data.imag, "Im")]):
            ax = axes[row, col]
            # lat/lon are cell CENTRES -> shading='nearest'. matplotlib places
            # cells by coordinate value, so north (max lat) lands at the top
            # regardless of the array's native N->S row order.
            im = ax.pcolormesh(lon, lat, part, shading='nearest', cmap='RdBu_r')
            ax.set_title(f"{title} ({ri})", fontsize=10)
            plt.colorbar(im, ax=ax)
            ax.set_xlabel("longitude (deg)")
            ax.set_ylabel("latitude (deg)")
    plt.suptitle(suptitle, fontsize=10)
    plt.tight_layout()
    plt.savefig(fname, dpi=140)
    plt.close(fig)
    print(f"Saved -> {fname}")


def plot_file(fname):
    info = parse_header(fname)
    lon, lat, comps = read_surface(fname)
    title = (f"T={info['period']}s, mode {info['mode']} -- "
             f"{info['SOLVER']}, {info['OUTPUT_CONVENTION']}")
    print(f"{os.path.basename(fname)}: {len(lon)} lon x {len(lat)} lat  [{title}]")

    out_e, out_h = out_png_names(fname)
    _plot(lon, lat,
          [(comps['Etheta'], r"$E_\theta$"), (comps['Ephi'], r"$E_\phi$")],
          f"E surface field, {title}", out_e)
    _plot(lon, lat,
          [(comps['Hr'], r"$H_r$"), (comps['Htheta'], r"$H_\theta$"), (comps['Hphi'], r"$H_\phi$")],
          f"H surface field, {title}", out_h)

    for name in ('Etheta', 'Ephi', 'Hr', 'Htheta', 'Hphi'):
        print(f"  max|{name}|={np.max(np.abs(comps[name])):.3e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='+', help='one or more *.surface files')
    args = ap.parse_args()
    for f in args.files:
        plot_file(f)


if __name__ == '__main__':
    main()
