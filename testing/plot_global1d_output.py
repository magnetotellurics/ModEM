"""
Read global1d's actual Fortran .efield/.hfield output files and plot Ex/Ey
(and Hx/Hy/Hz) in the SAME style as demo_general_lm.py's MT-convention
figures (lat x [-180,180] longitude, Re/Im column pairs, RdBu_r pcolormesh).

Unlike demo_general_lm.py (which computes fields analytically via
spherical_em_induction.py), this script reads REAL Fortran output, so there
is no "source K" to plot -- only the field panels (Ex/Ey for E; Hx/Hy/Hz for
H) are shown.

--- File format (confirmed directly against field1d.f90/FWD1D.f90's
    write_cvector output and cross-checked against the MATLAB read_field.m
    reader -- see CLAUDE.md) ---

Line 1: comment header (skip)
Line 2: "  1" (unused)
Line 3: "  nx  ny  nz  gridType"   nx=#phi cells, ny=#theta cells, nz=#r cells
Lines 4..: "i j k  Re(x) Im(x) Re(y) Im(y) Re(z) Im(z)"  for
           (nx+1)*(ny+1)*(nz+1) rows, i (phi) fastest, then j (theta), then k (r)

Component convention (GridDef.f90): cvector%x = phi component, %y = theta
component, %z = r component -- NOT Cartesian x/y/z, and this x/phi,y/theta,
z/r slot assignment is FIXED regardless of grid staggering (see below); only
which physical field (H or E) occupies EDGE vs FACE changes.
(Confirmed against the actual file: MT.*.efield's "x" column is zero for
k=1 (top of the domain, zero tangential-phi E in air far from source) while
the "y" (theta) column is nonzero even near the top -- consistent with
x=phi, y=theta, z=r.)

--- Two different named "Hx/Hy/Hz" (or Ex/Ey/Ez) conventions -- do not
    confuse them; this script's PLOT_MT_CONVENTION targets the SECOND one
    (per user, 2026-07-22) ---

field1d.f90 was ORIGINALLY written for the global code following Uyeshima &
Schultz (2000):
  1) H on the primary grid (EDGEs), E on the dual grid (FACEs)
  2) Hx = H_phi, the *longitudinal* component (East-West axis)
  3) Hy = H_theta, the *latitudinal* component, North-to-South (theta = CO-latitude)
  4) Hz = H_r, points down
  5) global grid: includes the poles and zero-longitude wrapping

ModEM's regional convention is different:
  1) E on the primary grid (EDGEs), H on the dual grid (FACEs) -- OPPOSITE
     staggering from the global convention above
  2) Hy = H_phi, the *longitudinal* component (East-West axis)
  3) Hx = H_theta, the *latitudinal* component, South-to-North (theta =
     LATITUDE here, not co-latitude -- opposite sign sense from Uyeshima &
     Schultz's theta)
  4) Hz = H_r, points down
  5) regional grid: no special treatment for poles or zero-longitude wrap

FWD1D.f90 hardcodes primary_grid='E', i.e. E on EDGEs / H on FACEs -- this
matches the ModEM REGIONAL staggering, not the original global Uyeshima &
Schultz one, even though the underlying field1d.f90 physics/module is the
same code either way. (This is also why the l=1,m=-1 Ex/Ey/Hz-vs-Hx/Hy
diagnostic investigation in CLAUDE.md had to trace through Tnrp/Tns, the
FACE-branch potentials, rather than Tnr/Tnsp.)

PLOT_MT_CONVENTION=True below targets the ModEM regional naming (since
that's what PrimaryField.m and everything else in this comparison expects):
  Ex = -Etheta (South->North), Ey = +Ephi (East-West axis, positive=East --
       confirmed 2026-07-22: "East to West" in both convention descriptions
       is a descriptive axis label, NOT a sign statement; standard eastward
       phi-hat is positive, no extra flip), Ez/Hz = +Er/+Hr (down --
       confirmed 2026-07-22: global1d's Hr/Er ALREADY point down, i.e. no
       extra sign flip is needed to get Ez/Hz; H%z/E%z's leading minus sign
       in field1d.f90's formula is baked into the Mie-potential derivation
       itself, not a separate up->down relabeling); same pattern for H
       (Hx=-Htheta, Hy=+Hphi, Hz=+Hr).
PLOT_MT_CONVENTION=False shows the RAW physical Etheta/Ephi/Hr/Htheta/Hphi
components directly (no sign flips, no relabeling) -- these are the same
raw numbers either named convention above is built from; this script just
doesn't relabel them to either scheme in this mode.

The "mysterious factor" from the user's MATLAB snippet (an unresolved
overall sign/phase, still open per CLAUDE.md) is kept here as an explicit,
separately-toggleable FACTOR variable rather than baked in.

--- Regional / non-global grids (added 2026-07-25) ---

The .efield/.hfield files themselves do NOT store the actual (theta,phi) of
each node -- only i,j,k index triples plus field values (see the format
note above). read_field()/to_latlon_mt() historically ASSUMED a global grid
uniformly covering phi=[0,360], theta=[0,180] (np.linspace(0,360,nx+1) /
np.linspace(0,180,ny+1)), which is correct for e.g. global.1.0x1.0.grd but
WRONG for a regional grid (e.g. USA.0.25x0.25.grd, or any grid passed
through FWD1D.f90's optional "fake pole" recentering, see
recenter_grid_fake_pole in FWD1D.f90) -- regional grids have non-uniform
cell widths and do not span the full sphere.

If --grid is given (see parse_args()/CLI usage below), the actual node
positions are read directly from that .grd file via read_grid_file() (a
Python port of GridDef.f90's read_grid, same file format, same formulas)
instead of being assumed. If --fake-center-lat/--fake-center-lon are also
given, the SAME "fake pole" shift FWD1D.f90 applied at generation time
(recenter_grid_fake_pole) is reproduced here via recenter_fake_pole() --
these must match whatever values were actually passed to FWD1D when the
.efield/.hfield files were produced, or the plotted lat/lon axes will be
wrong. Omit --grid to keep the original global-uniform-grid behavior
unchanged.

CLI usage:
    python plot_global1d_output.py [efield_file] [hfield_file]
        [--grid GRID_FILE] [--fake-center-lat LAT_DEG --fake-center-lon LON_DEG]
Run with -h for the full option list. Both positional arguments are
optional and default to a built-in example (DEFAULT_EFIELD_FILE/
DEFAULT_HFIELD_FILE below). Output PNG filenames are derived from the
input file names (see output_png_names()), so different runs never
overwrite each other's plots.

--- Mid-theta/mid-phi zero-padding (found 2026-07-27) ---

write_cvector (sg_vector.f90) always emits a FIXED (nx+1)x(ny+1)x(nz+1)
grid per file, but each of cvector's %x/%y/%z is only physically valid over
an nx- or ny-sized subrange, depending on whether that component is
node-staggered or mid-staggered in that axis (create_cvector allocates e.g.
EDGE %x as (nx,ny+1,nz+1) -- mid phi, node theta -- vs FACE %x as
(nx+1,ny,nz) -- node phi, mid theta). write_cvector zero-fills the WHOLE
(nx+1,ny+1,nz+1) output array first, then copies only the valid subrange in
(sg_vector.f90 lines ~696-715) -- so the leftover row/column for a
mid-staggered component is a genuine zero in the FILE, not a "the source is
undefined there" physics statement. Confirmed against write_cvector's own
code:
    EDGE: %x mid-phi/node-theta, %y node-phi/mid-theta, %z node-phi/node-theta
    FACE: %x node-phi/mid-theta, %y mid-phi/node-theta, %z mid-phi/mid-theta
(see _KIND below). read_field()/to_latlon_mt() previously assumed every
component had full node-node data, so any mid-staggered component showed a
spurious zero row (mid-theta) or column (mid-phi) at the LAST index --
exactly the "Hy, Hz, Ex show zero at the South pole / minimum-latitude
boundary" symptom reported for primary_grid='E' (where MT Hy=+H%x [FACE,
mid-theta], Hz=+H%z [FACE, mid-AND-mid], Ex=-E%y [EDGE, mid-theta]).
to_latlon_mt now looks up each raw component's (phi_kind, theta_kind) from
_KIND and, for 'mid' axes, drops the padding row/column and positions the
remaining data at true CELL-CENTER coordinates instead of node coordinates
-- this changes to_latlon_mt's return value from one shared (lon,lat) pair
for all three components to three separate (lon,lat) pairs, since a
mid-staggered component's grid is genuinely offset from a node-staggered
one's. See the module-level _plot()/`__main__` usage for how the three are
now plotted on their own axes.
"""
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
# EFIELD_FILE/HFIELD_FILE/GRID_FILE/FAKE_CENTER_LAT/FAKE_CENTER_LON are now
# command-line arguments (see __main__ / parse_args() below) -- the values
# below are only used as their DEFAULTS when the script is run with no
# arguments, so `python plot_global1d_output.py` still works out of the box.
DEFAULT_EFIELD_FILE = r"C:\Users\Anna Kelbert\Developer\ModEM-global1d\MTsource\MT.1000sec.test_m-ve.sunegbert2012.E-grid.T01.efield"
DEFAULT_HFIELD_FILE = r"C:\Users\Anna Kelbert\Developer\ModEM-global1d\MTsource\MT.1000sec.test_m-ve.sunegbert2012.E-grid.T01.hfield"

PERIOD = 1000.0
SCALING = PERIOD / 5.0        # matches the user's MATLAB snippet exactly
FACTOR = 1.0                  # +0j  # "mysterious factor" placeholder -- still open, see CLAUDE.md
R_INDEX = 13                  # 1-indexed; k=13 = first layer at/below the surface
                               # (12 air layers above, per global.1.0x1.0.grd header "...12 0 47")

# If True, plot E as Ex (South->North, = -Etheta), Ey (West->East, = +Ephi),
# and H as Hx (=-Htheta), Hy (=+Hphi), Hz (down, =+Hr -- global1d's Hr/Er
# already point down, no extra flip needed), on a
# latitude x [-180,180] longitude grid -- this is the ModEM REGIONAL naming
# convention (Hx=Htheta-as-latitude, Hy=Hphi), matching PrimaryField.m
# exactly -- NOT the original global Uyeshima & Schultz (2000) convention
# field1d.f90 was originally written for, which swaps the Hx/Hy letters and
# uses co-latitude instead of latitude (see module docstring above). If
# False, plot Etheta/Ephi and Hr/Htheta/Hphi directly (raw, unrelabeled) on
# a colatitude x [0,360] grid.
PLOT_MT_CONVENTION = True


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot global1d .efield/.hfield output (Ex/Ey, Hx/Hy/Hz)."
    )
    parser.add_argument("efield_file", nargs="?", default=DEFAULT_EFIELD_FILE,
                         help="FWD1D .efield output file (default: the built-in example)")
    parser.add_argument("hfield_file", nargs="?", default=DEFAULT_HFIELD_FILE,
                         help="FWD1D .hfield output file (default: the built-in example)")
    parser.add_argument("--grid", default="", metavar="GRID_FILE",
                         help="Optional .grd file used to generate efield_file/hfield_file -- "
                              "plots on the ACTUAL (possibly regional, non-uniform) grid node "
                              "positions instead of assuming a global uniform grid. See "
                              "read_grid_file()/to_latlon_mt()'s docstrings.")
    parser.add_argument("--fake-center-lat", type=float, default=None, dest="fake_center_lat",
                         metavar="LAT_DEG",
                         help="If --grid was recentered via FWD1D's optional \"fake pole\" "
                              "arguments, the SAME latitude (deg) that was passed to FWD1D. "
                              "Must be given together with --fake-center-lon.")
    parser.add_argument("--fake-center-lon", type=float, default=None, dest="fake_center_lon",
                         metavar="LON_DEG",
                         help="See --fake-center-lat -- the matching longitude (deg).")
    return parser.parse_args()


def output_png_names(efield_file, hfield_file, mt_convention, r_index):
    """Derive descriptive output filenames from the INPUT file names, e.g.
    'MT.200ohmm.1000sec.Mode1.E-grid.T01.efield' + MT convention + r_index=13
    -> 'MT.200ohmm.1000sec.Mode1.E-grid.T01.Ex_Ey.r13.png' -- so the plot
    filename alone identifies which run and which components it came from."""
    def stem(path, suffix):
        base = os.path.basename(path)
        if base.lower().endswith(suffix):
            base = base[:-len(suffix)]
        return base

    comp_e = "Ex_Ey" if mt_convention else "Etheta_Ephi"
    comp_h = "Hx_Hy_Hz" if mt_convention else "Hr_Htheta_Hphi"
    out_e = f"{stem(efield_file, '.efield')}.{comp_e}.r{r_index}.png"
    out_h = f"{stem(hfield_file, '.hfield')}.{comp_h}.r{r_index}.png"
    return out_e, out_h


# ----------------------------------------------------------------------
# Fortran cvector-file reader (matches write_cvector / read_field.m exactly)
# ----------------------------------------------------------------------
def read_field(fname):
    with open(fname) as f:
        f.readline()  # header comment line
        f.readline()  # "  1"
        nx, ny, nz, gridType = f.readline().split()
        nx, ny, nz = int(nx), int(ny), int(nz)
        gridType = gridType.strip()

        n = (nx + 1) * (ny + 1) * (nz + 1)
        data = np.loadtxt(f, max_rows=n)

    # columns: i j k  xr xi yr yi zr zi  (1-indexed i,j,k; i=phi fastest)
    # NOTE: the file's own column labels are the generic cvector "x,y,z" slots,
    # but per GridDef.f90 these are NOT Cartesian -- they are phi,theta,r.
    # Named accordingly here (Cphi/Ctheta/Cr, not Cx/Cy/Cz) so nothing
    # downstream can mistake them for an already-MT-convention Cartesian
    # field -- global1d's raw output is phi/theta/r; MT convention (Ex/Ey/Ez)
    # is only applied later, explicitly, in to_latlon_mt() when requested.
    Cphi = data[:, 3] + 1j * data[:, 4]     # phi component
    Ctheta = data[:, 5] + 1j * data[:, 6]   # theta component
    Cr = data[:, 7] + 1j * data[:, 8]       # r component

    shape = (nx + 1, ny + 1, nz + 1)
    # Fortran/column-major reshape: i (first file column, fastest) -> axis 0
    Cphi = Cphi.reshape(shape, order='F')
    Ctheta = Ctheta.reshape(shape, order='F')
    Cr = Cr.reshape(shape, order='F')

    return dict(phi=Cphi, theta=Ctheta, r=Cr, nx=nx, ny=ny, nz=nz, gridType=gridType)


def read_grid_file(fname):
    """
    Python port of GridDef.f90's read_grid -- SAME .grd file format, SAME
    formulas, so it reproduces exactly the grid FWD1D.f90 actually used.

    File layout (whitespace-separated, one array per non-blank line):
      line 1: "nx ny nzAir nzCrust nzEarth [corner_lon]" -- the OPTIONAL
              trailing corner_lon (degrees) is the absolute longitude of the
              grid's lower-left corner (i=1, minimum longitude/westernmost).
              Co-latitude (line 3) and radius (line 4) are ALREADY absolute
              in the file, so no separate anchor is needed for those. If
              corner_lon is absent (older grid files), it defaults to 0,
              exactly reproducing the historical hardcoded ph[0]=0 behavior.
      line 2: nx values -- phi CELL WIDTHS, degrees (NOT absolute positions)
      line 3: ny+1 values -- theta NODE positions, degrees, ABSOLUTE (colatitude)
      line 4: nz+1 values -- radius NODE positions, km (top of air down to core)

    Returns dict(nx, ny, nz, th_deg (ny+1,), ph_deg (nx+1,), r_km (nz+1,),
    has_corner_lon (bool)).
    """
    with open(fname) as f:
        lines = [ln for ln in f if ln.strip()]   # drop blank separator lines
    header = lines[0].split()
    nx, ny, nzAir, nzCrust, nzEarth = (int(v) for v in header[:5])
    nz = nzAir + nzCrust + nzEarth
    has_corner_lon = len(header) > 5
    corner_lon = float(header[5]) if has_corner_lon else 0.0

    dp_deg = np.array(lines[1].split(), dtype=float)
    th_deg = np.array(lines[2].split(), dtype=float)    # already absolute co-latitude
    r_km = np.array(lines[3].split(), dtype=float)
    if len(dp_deg) != nx:
        raise ValueError(f"{fname}: expected {nx} phi cell widths, got {len(dp_deg)}")
    if len(th_deg) != ny + 1:
        raise ValueError(f"{fname}: expected {ny+1} theta nodes, got {len(th_deg)}")
    if len(r_km) != nz + 1:
        raise ValueError(f"{fname}: expected {nz+1} radius nodes, got {len(r_km)}")

    phi_shift = corner_lon % 360.0
    ph_deg = phi_shift + np.concatenate([[0.0], np.cumsum(dp_deg)])

    return dict(nx=nx, ny=ny, nz=nz, th_deg=th_deg, ph_deg=ph_deg, r_km=r_km,
                has_corner_lon=has_corner_lon)


def recenter_fake_pole(th_deg, ph_deg, center_lat_deg, center_lon_deg):
    """
    Python port of FWD1D.f90's recenter_grid_fake_pole -- reproduces the
    SAME uniform (theta,phi) node shift, so this must be called with the
    SAME center_lat_deg/center_lon_deg that were actually passed to FWD1D
    when generating the field file being plotted. See that subroutine's
    header comment (FWD1D.f90) for the full explanation of what this is
    (a translation, not a true spherical rotation) and why.

    Returns (th_deg_new, ph_deg_new), same shapes as the inputs.
    """
    theta_c_old = (th_deg.min() + th_deg.max()) / 2.0
    phi_c_old = (ph_deg.min() + ph_deg.max()) / 2.0

    theta_c_new = 90.0 - center_lat_deg
    phi_c_new = center_lon_deg % 360.0

    dtheta = theta_c_new - theta_c_old
    dphi = phi_c_new - phi_c_old

    th_deg_new = th_deg + dtheta
    ph_deg_new = (ph_deg + dphi) % 360.0

    if th_deg_new.min() <= 0.0 or th_deg_new.max() >= 180.0:
        raise ValueError(
            "recenter_fake_pole: requested center pushes the grid across a pole "
            "(theta out of the open interval (0,180) deg) -- choose a fake center "
            "farther from the poles, or a smaller/narrower grid extent"
        )

    return th_deg_new, ph_deg_new


# Per-(gridType, raw component) staggering, derived directly from
# create_cvector/write_cvector (sg_vector.f90) -- see module docstring,
# "Mid-theta/mid-phi zero-padding". 'node': full nx+1/ny+1 valid positions,
# matching the file's node grid exactly. 'mid': only nx/ny valid positions
# (a cell-center quantity); the file's leftover row/column is zero-padding,
# not data.
_KIND = {
    'EDGE': {'phi': ('mid', 'node'), 'theta': ('node', 'mid'), 'r': ('node', 'node')},
    'FACE': {'phi': ('node', 'mid'), 'theta': ('mid', 'node'), 'r': ('mid', 'mid')},
}


def _theta_positions(theta_deg_full, kind):
    """Returns (theta_positions, row_slice). 'mid' drops the last (padding)
    row and positions the remaining ny rows at true cell centers
    (average of adjacent node values) rather than node positions."""
    if kind == 'mid':
        return (theta_deg_full[:-1] + theta_deg_full[1:]) / 2.0, slice(0, -1)
    return theta_deg_full, slice(None)


def _phi_positions(phi_deg_full, kind, is_global):
    """Returns (phi_positions, col_slice), BEFORE signed-longitude
    conversion. 'mid': the last column is always zero-padding (never real
    data, regional or global) -- drop it, position the remaining nx columns
    at true cell centers. 'node', global: the last column is a TRUE
    360deg-wraparound DUPLICATE of the first (not padding) -- drop it for
    that different reason. 'node', regional: no duplicate: keep all nx+1
    columns."""
    if kind == 'mid':
        return (phi_deg_full[:-1] + phi_deg_full[1:]) / 2.0, slice(0, -1)
    if is_global:
        return phi_deg_full[:-1], slice(0, -1)
    return phi_deg_full, slice(None)


def _to_signed_lon_sorted(positions_deg, data_2d):
    """Convert ascending phi positions (any range) to signed longitude
    [-180,180) and sort for a monotonic pcolormesh axis (generalizes the
    old global-only searchsorted/circular-shift trick to also handle a
    regional patch straddling the +-180 meridian after recentering)."""
    lon = ((positions_deg + 180.0) % 360.0) - 180.0
    order = np.argsort(lon, kind='stable')
    return lon[order], data_2d[:, order]


def _component_latlon(raw2d, ph_deg_full, th_deg_full, phi_kind, theta_kind, is_global):
    """raw2d: shape (nx+1, ny+1), one component (phi, theta, or r), already
    sliced at a single r_index. Returns (lon_deg, lat_deg, data) with data
    shape (n_lat, n_lon), ready for pcolormesh(lon_deg, lat_deg, data)."""
    theta_pos, row_slice = _theta_positions(th_deg_full, theta_kind)
    phi_pos, col_slice = _phi_positions(ph_deg_full, phi_kind, is_global)
    data = raw2d[col_slice, row_slice].T   # -> (theta-like, phi-like)
    lat_deg = 90.0 - theta_pos
    lon_deg, data = _to_signed_lon_sorted(phi_pos, data)
    return lon_deg, lat_deg, data


def to_latlon_mt(field, r_index, scaling, factor, mt_convention=True,
                  th_deg=None, ph_deg=None):
    """
    field: dict from read_field() (phi/theta/r component arrays, shape
           (nx+1, ny+1, nz+1), axis0=phi ascending, axis1=theta
           (colatitude) ascending, axis2=r-index)

    th_deg/ph_deg (optional): the ACTUAL node positions (degrees) for this
    grid, e.g. from read_grid_file() (+ optionally recenter_fake_pole()).
    If omitted (both None), falls back to the ORIGINAL assumption of a
    global grid uniformly covering phi=[0,360], theta=[0,180]
    (np.linspace) -- correct for e.g. global.1.0x1.0.grd, WRONG for a
    regional grid. Pass them explicitly for any regional/non-global grid.

    Returns a dict of THREE independent (lon_deg, lat_deg, data) tuples --
    one per output field component -- NOT a single shared (lon,lat) pair,
    because a mid-staggered component's true grid is offset (by half a
    cell) from a node-staggered one's, and has one fewer valid row/column
    (see module docstring, "Mid-theta/mid-phi zero-padding"). Each data
    array has shape (n_lat, n_lon) for that component specifically, ready
    for pcolormesh(lon_deg, lat_deg, data).

    If mt_convention: keys are 'x','y','z' with Fx=South->North (=-Ftheta),
    Fy=West->East (=+Fphi), Fz=down (=+Fr -- global1d's Hr/Er already point
    down, no extra flip) -- the ModEM REGIONAL naming convention, matches
    PrimaryField.m exactly (NOT the original global Uyeshima & Schultz
    convention field1d.f90 was written for -- see module docstring). Note
    Fx's grid comes from the raw THETA component's staggering, Fy's from
    raw PHI, Fz's from raw R (the sign flip doesn't change staggering).

    If not mt_convention: keys are 'theta','phi','r', values unrelabeled,
    no sign flips.
    """
    nx, ny = field['nx'], field['ny']
    regional = th_deg is not None or ph_deg is not None
    if regional and (th_deg is None or ph_deg is None):
        raise ValueError("to_latlon_mt: th_deg and ph_deg must be given together")

    if not regional:
        ph_deg_full = np.linspace(0.0, 360.0, nx + 1)
        th_deg_full = np.linspace(0.0, 180.0, ny + 1)     # colatitude
    else:
        ph_deg_full = np.asarray(ph_deg, dtype=float)
        th_deg_full = np.asarray(th_deg, dtype=float)
        if len(ph_deg_full) != nx + 1 or len(th_deg_full) != ny + 1:
            raise ValueError(
                f"to_latlon_mt: grid size mismatch -- field file has nx={nx},ny={ny} "
                f"(expects {nx+1} phi / {ny+1} theta nodes), grid gave "
                f"{len(ph_deg_full)} phi / {len(th_deg_full)} theta nodes"
            )
    # "does phi genuinely wrap the full sphere" -- checked from the actual
    # extent rather than just "was --grid given", so an explicitly-supplied
    # GLOBAL .grd file (e.g. global.1.0x1.0.grd via --grid) is still handled
    # as global (drops the true 360deg-duplicate node-phi column) rather
    # than being misclassified as a non-wrapping regional patch.
    is_global = (ph_deg_full[-1] - ph_deg_full[0]) > 359.99

    gridType = field['gridType']
    if gridType not in _KIND:
        raise ValueError(f"to_latlon_mt: unknown gridType {gridType!r}")
    kind = _KIND[gridType]

    Ftheta = field['theta'][:, :, r_index - 1] * scaling * factor
    Fphi = field['phi'][:, :, r_index - 1] * scaling * factor
    Fr = field['r'][:, :, r_index - 1] * scaling * factor

    if not mt_convention:
        lon_t, lat_t, Dtheta = _component_latlon(Ftheta, ph_deg_full, th_deg_full, *kind['theta'], is_global)
        lon_p, lat_p, Dphi = _component_latlon(Fphi, ph_deg_full, th_deg_full, *kind['phi'], is_global)
        lon_r, lat_r, Dr = _component_latlon(Fr, ph_deg_full, th_deg_full, *kind['r'], is_global)
        return {'theta': (lon_t, lat_t, Dtheta), 'phi': (lon_p, lat_p, Dphi), 'r': (lon_r, lat_r, Dr)}

    lon_x, lat_x, Dx = _component_latlon(Ftheta, ph_deg_full, th_deg_full, *kind['theta'], is_global)
    lon_y, lat_y, Dy = _component_latlon(Fphi, ph_deg_full, th_deg_full, *kind['phi'], is_global)
    lon_z, lat_z, Dz = _component_latlon(Fr, ph_deg_full, th_deg_full, *kind['r'], is_global)
    return {'x': (lon_x, lat_x, -Dx),   # South->North
            'y': (lon_y, lat_y, Dy),    # West->East
            'z': (lon_z, lat_z, Dz)}    # down (already the case, see module docstring)


def _plot(panels, suptitle, fname, xlabel, ylabel, invert_y):
    """panels: list of (lon_or_phi_deg, lat_or_theta_deg, data, title) --
    EACH panel carries its OWN grid now (see to_latlon_mt's docstring: a
    mid-staggered component's grid is genuinely offset from a node-staggered
    one's, so they can no longer share one x_deg/y_deg pair)."""
    n = len(panels)
    fig, axes = plt.subplots(n, 2, figsize=(10, 4 * n))
    for row, (x_deg, y_deg, data, title) in enumerate(panels):
        for col, (part, ri) in enumerate([(data.real, "Re"), (data.imag, "Im")]):
            ax = axes[row, col]
            im = ax.pcolormesh(x_deg, y_deg, part, shading='auto', cmap='RdBu_r')
            ax.set_title(f"{title} ({ri})", fontsize=10)
            plt.colorbar(im, ax=ax)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            if invert_y:
                ax.invert_yaxis()   # colatitude: put theta=0 (N pole) at top
    plt.suptitle(suptitle, fontsize=10)
    plt.tight_layout()
    plt.savefig(fname, dpi=140)
    print(f"Saved -> {fname}")


if __name__ == "__main__":
    args = parse_args()

    Efield = read_field(args.efield_file)
    Hfield = read_field(args.hfield_file)
    print(f"E file: nx={Efield['nx']} ny={Efield['ny']} nz={Efield['nz']} gridType={Efield['gridType']}  ({args.efield_file})")
    print(f"H file: nx={Hfield['nx']} ny={Hfield['ny']} nz={Hfield['nz']} gridType={Hfield['gridType']}  ({args.hfield_file})")

    th_deg = ph_deg = None
    if args.grid:
        grid = read_grid_file(args.grid)
        th_deg, ph_deg = np.asarray(grid['th_deg']), np.asarray(grid['ph_deg'])
        print(f"Grid file: nx={grid['nx']} ny={grid['ny']} nz={grid['nz']}  (from {args.grid})")
        if grid['has_corner_lon']:
            print(f"  Grid has an explicit lower-left corner longitude (ph[0]={ph_deg[0]:.4f} deg)")
        else:
            print(f"  Grid has no corner longitude in its header -- using historical default (ph[0]=0)")
        if args.fake_center_lat is not None or args.fake_center_lon is not None:
            if args.fake_center_lat is None or args.fake_center_lon is None:
                raise ValueError("--fake-center-lat and --fake-center-lon must be given together")
            th_deg, ph_deg = recenter_fake_pole(th_deg, ph_deg, args.fake_center_lat, args.fake_center_lon)
            print(f"Applied fake-pole recentering to (lat,lon)=({args.fake_center_lat},{args.fake_center_lon})")
    elif args.fake_center_lat is not None or args.fake_center_lon is not None:
        raise ValueError("--fake-center-lat/--fake-center-lon require --grid to also be given")

    OUT_E_PNG, OUT_H_PNG = output_png_names(args.efield_file, args.hfield_file, PLOT_MT_CONVENTION, R_INDEX)

    E = to_latlon_mt(Efield, R_INDEX, SCALING, FACTOR, PLOT_MT_CONVENTION, th_deg=th_deg, ph_deg=ph_deg)
    H = to_latlon_mt(Hfield, R_INDEX, SCALING, FACTOR, PLOT_MT_CONVENTION, th_deg=th_deg, ph_deg=ph_deg)

    if PLOT_MT_CONVENTION:
        Ex_lon, Ex_lat, Ex = E['x']
        Ey_lon, Ey_lat, Ey = E['y']
        Hz_lon, Hz_lat, Hz = H['z']
        Hx_lon, Hx_lat, Hx = H['x']
        Hy_lon, Hy_lat, Hy = H['y']
        e_panels = [
            (Ex_lon, Ex_lat, Ex, r"$E_x$ (S$\to$N)"),
            (Ey_lon, Ey_lat, Ey, r"$E_y$ (W$\to$E)"),
        ]
        h_panels = [
            (Hz_lon, Hz_lat, Hz, r"$H_z$ (down)"),
            (Hx_lon, Hx_lat, Hx, r"$H_x$ (S$\to$N)"),
            (Hy_lon, Hy_lat, Hy, r"$H_y$ (W$\to$E)"),
        ]
        xlabel, ylabel, invert_y = "longitude (deg)", "latitude (deg)", False
        e_names, h_names = ("Ex", "Ey"), ("Hx", "Hy", "Hz")
    else:
        Et_lon, Et_lat, Etheta = E['theta']
        Ep_lon, Ep_lat, Ephi = E['phi']
        Ht_lon, Ht_lat, Ftheta = H['theta']
        Hp_lon, Hp_lat, Fphi = H['phi']
        Hr_lon, Hr_lat, Fr = H['r']
        e_panels = [
            (Et_lon, Et_lat, Etheta, r"$E_\theta$"),
            (Ep_lon, Ep_lat, Ephi, r"$E_\phi$"),
        ]
        h_panels = [
            (Hr_lon, Hr_lat, Fr, r"$H_r$"),
            (Ht_lon, Ht_lat, Ftheta, r"$H_\theta$"),
            (Hp_lon, Hp_lat, Fphi, r"$H_\phi$"),
        ]
        xlabel, ylabel, invert_y = "longitude / phi (deg)", "colatitude / theta (deg)", True
        e_names, h_names = ("Etheta", "Ephi"), ("Hr", "Htheta", "Hphi")

    _plot(e_panels, f"global1d output: E field, T={PERIOD}s, r-index={R_INDEX}",
          OUT_E_PNG, xlabel, ylabel, invert_y)
    _plot(h_panels, f"global1d output: H field, T={PERIOD}s, r-index={R_INDEX}",
          OUT_H_PNG, xlabel, ylabel, invert_y)

    for name, (_, _, data, _) in zip(e_names, e_panels):
        print(f"max|{name}|={np.max(np.abs(data)):.3e}", end="  ")
    print()
    for name, (_, _, data, _) in zip(h_names, h_panels):
        print(f"max|{name}|={np.max(np.abs(data)):.3e}", end="  ")
    print()
