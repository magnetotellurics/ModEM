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
    reader) ---

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

--- Output-convention-aware plotting (rewritten 2026-07-26 to match
    output_convention.f90's Part 2 rewrite) ---

field1d.f90/field1d_s2.f90's raw component slots are FIXED
regardless of any relabeling: cvector%x = phi component, %y = theta
component, %z = r component (see above) -- this never changes. What USED TO
change, via this script's own PLOT_MT_CONVENTION toggle, was an extra
Python-side relabeling step (Ex=-Etheta, Ey=+Ephi, Hz=+Hr, Hx=-Htheta,
Hy=+Hphi) applied AFTER reading the raw file, to reproduce the ModEM
regional "MT convention" naming by hand.

That relabeling is now done, correctly and completely, on the FORTRAN side
by EARTH/FWD/output_convention.f90's apply_output_convention -- FWD1D.f90's
output files are ALREADY in whatever convention was requested
(OUTPUT_CONVENTION), including the theta index-reversal + component
negation (colatitude N->S <-> latitude S->N) and r index-reversal + negation
(up/radius <-> down/depth) that PLOT_MT_CONVENTION used to approximate by
hand. Applying the OLD Python-side relabeling on top of an ALREADY-relabeled
file would double-apply it and silently produce a wrong (double-flipped)
result -- so this script no longer does ANY relabeling of its own. It plots
the file's raw phi/theta/r component slots directly, and the ONLY
convention-dependent adjustment it makes is to the GRID AXIS interpretation
(mapping array index -> physical latitude/depth correctly when the file's
own recorded theta_convention/r_convention indicate the array was
reversed relative to the .grd file's native colatitude/radius ordering --
see read_field()'s header-parsing and to_latlon_mt()'s th_deg_full/r_index
handling below). The active convention's name (and, if extracted, its five
dimension values) is read directly from the file's own header line and
printed in the figure title -- see read_field()'s CONVENTION_KEYS parsing.

Older field files, written before this convention metadata existed, have no
such header fields; read_field() falls back to convention=None in that case
and to_latlon_mt() makes no theta/r reversal adjustment (reproducing the
original, pre-Part-2 axis behavior) -- a printed note flags when this
fallback is in effect so a plot's axes are never silently mislabeled.

The "mysterious factor" from the user's MATLAB snippet (an unresolved
overall sign/phase, still open) is kept here as an explicit,
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
DEFAULT_EFIELD_FILE = r"C:\Users\Anna Kelbert\Developer\ModEM-global1d\MTsource\MT.1000sec.test_m-ve.s2.E-grid.T01.efield"
DEFAULT_HFIELD_FILE = r"C:\Users\Anna Kelbert\Developer\ModEM-global1d\MTsource\MT.1000sec.test_m-ve.s2.E-grid.T01.hfield"

PERIOD = 1000.0
SCALING = PERIOD / 5.0        # matches the user's MATLAB snippet exactly
FACTOR = 1.0                  # +0j  # "mysterious factor" placeholder -- still open
R_INDEX = 13                  # 1-indexed; k=13 = first layer at/below the surface
                               # (12 air layers above, per global.1.0x1.0.grd header "...12 0 47")
                               # -- this indexes into the field FILE's own k-axis; if the file's
                               # r_convention differs from the solvers' native R_UP, to_latlon()
                               # remaps it automatically so R_INDEX still means "the same physical
                               # layer counting from the top/air side", see NATIVE_R_CONVENTION below.

# Both field1d.f90 and field1d_s2.f90's native r_convention is
# R_UP (confirmed 2026-07-26, see output_convention.f90) -- used
# below to figure out whether a given output file's r-axis was reversed
# relative to the .grd file's own (always top-to-bottom) k-ordering.
NATIVE_R_CONVENTION = 'R_UP'


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
                              "read_grid_file()/to_latlon()'s docstrings.")
    parser.add_argument("--fake-center-lat", type=float, default=None, dest="fake_center_lat",
                         metavar="LAT_DEG",
                         help="If --grid was recentered via FWD1D's optional \"fake pole\" "
                              "arguments, the SAME latitude (deg) that was passed to FWD1D. "
                              "Must be given together with --fake-center-lon.")
    parser.add_argument("--fake-center-lon", type=float, default=None, dest="fake_center_lon",
                         metavar="LON_DEG",
                         help="See --fake-center-lat -- the matching longitude (deg).")
    return parser.parse_args()


def output_png_names(efield_file, hfield_file, convention_name, r_index):
    """Derive descriptive output filenames from the INPUT file names, e.g.
    'MT.200ohmm.1000sec.Mode1.E-grid.T01.efield' + convention 'EGBERTKELBERT2012'
    + r_index=13 -> 'MT.200ohmm.1000sec.Mode1.E-grid.T01.EGBERTKELBERT2012.Etheta_Ephi.r13.png'
    -- so the plot filename alone identifies which run, which convention, and
    which r-layer it came from. convention_name may be None (older files with
    no convention metadata) -- falls back to 'raw'. The efield_file/hfield_file
    stems are usually IDENTICAL except for their .efield/.hfield suffix (e.g.
    both from the same FWD1D run), so an explicit 'Etheta_Ephi'/'Hr_Htheta_Hphi'
    component tag is included as well -- WITHOUT it the two output filenames
    would collide and the H plot would silently overwrite the E plot (this bug
    existed transiently while output_png_names was being rewritten, caught by
    checking the actual saved filenames)."""
    def stem(path, suffix):
        base = os.path.basename(path)
        if base.lower().endswith(suffix):
            base = base[:-len(suffix)]
        return base

    tag = convention_name if convention_name else "raw"
    out_e = f"{stem(efield_file, '.efield')}.{tag}.Etheta_Ephi.r{r_index}.png"
    out_h = f"{stem(hfield_file, '.hfield')}.{tag}.Hr_Htheta_Hphi.r{r_index}.png"
    return out_e, out_h


# ----------------------------------------------------------------------
# Fortran cvector-file reader (matches write_cvector / read_field.m exactly)
# ----------------------------------------------------------------------

# Keys FWD1D.f90 writes into its header line, and how each parses. All are
# optional -- older field files (written before output_convention.f90
# existed) simply won't have any of these, and read_field() returns
# convention=None in that case. See FWD1D.f90's per-period loop, the two
# `hdr = "# FWD1D ..." // "OUTPUT_CONVENTION=" // ...` lines.
_CONVENTION_KEYS = ("OUTPUT_CONVENTION", "SOLVER", "time", "norm", "theta", "r")


def _parse_convention_header(line):
    """Extracts KEY=VALUE tokens (space-separated, no spaces within a value)
    for each of _CONVENTION_KEYS found in the header comment line. Returns
    None if none of the keys are present (older file, no convention
    metadata) -- otherwise a dict with whichever keys WERE found (all of
    them, for any file written by the current FWD1D.f90)."""
    found = {}
    for tok in line.split():
        if "=" in tok:
            key, _, val = tok.partition("=")
            if key in _CONVENTION_KEYS:
                found[key] = val
    return found if found else None


def read_field(fname):
    with open(fname) as f:
        header_line = f.readline()  # header comment line -- may carry convention metadata
        f.readline()  # "  1"
        nx, ny, nz, gridType = f.readline().split()
        nx, ny, nz = int(nx), int(ny), int(nz)
        gridType = gridType.strip()

        n = (nx + 1) * (ny + 1) * (nz + 1)
        data = np.loadtxt(f, max_rows=n)

    convention = _parse_convention_header(header_line)

    # columns: i j k  xr xi yr yi zr zi  (1-indexed i,j,k; i=phi fastest)
    # NOTE: the file's own column labels are the generic cvector "x,y,z" slots,
    # but per GridDef.f90 these are NOT Cartesian -- they are phi,theta,r.
    # Named accordingly here (Cphi/Ctheta/Cr, not Cx/Cy/Cz). These are the
    # RAW values as written by FWD1D.f90 -- already in whatever convention
    # is recorded in `convention` above (see output_convention.f90); this
    # script applies NO further relabeling to them, only grid-axis
    # adjustments (see to_latlon()).
    Cphi = data[:, 3] + 1j * data[:, 4]     # phi component
    Ctheta = data[:, 5] + 1j * data[:, 6]   # theta component
    Cr = data[:, 7] + 1j * data[:, 8]       # r component

    shape = (nx + 1, ny + 1, nz + 1)
    # Fortran/column-major reshape: i (first file column, fastest) -> axis 0
    Cphi = Cphi.reshape(shape, order='F')
    Ctheta = Ctheta.reshape(shape, order='F')
    Cr = Cr.reshape(shape, order='F')

    return dict(phi=Cphi, theta=Ctheta, r=Cr, nx=nx, ny=ny, nz=nz, gridType=gridType,
                convention=convention)


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
# "Mid-theta/mid-phi zero-padding". 'node': full nx+1/ny+1 (or nz+1 for r)
# valid positions, matching the file's node grid exactly. 'mid': only nx/ny
# (or nz for r) valid positions (a cell-center quantity); the file's leftover
# row/column/layer is zero-padding, not data. Three-tuples are
# (phi_kind, theta_kind, r_kind) -- added r_kind 2026-07-26 (see
# "r-axis mid/node off-by-one" fix in to_latlon()/_r_index_for_component()):
# EDGE %z and FACE %x/%y are 'mid'-r (size nz); everything else is 'node'-r
# (size nz+1).
_KIND = {
    'EDGE': {'phi': ('mid', 'node', 'node'), 'theta': ('node', 'mid', 'node'), 'r': ('node', 'node', 'mid')},
    'FACE': {'phi': ('node', 'mid', 'mid'), 'theta': ('mid', 'node', 'mid'), 'r': ('mid', 'mid', 'node')},
}


def _r_index_for_component(r_index, nz, r_kind, r_reversed):
    """Maps the user-facing r_index (meaning: 'the same physical layer,
    counting from the top/air side, in the solver's NATIVE r ordering') to
    the actual 1-indexed k-slot to read from the FILE's raw component array.

    If not r_reversed: file ordering already matches native -- no change.

    If r_reversed (apply_output_convention's r-transform ran): the reversal
    happens on the SOLVER's OWN cvector array, whose r-axis size is nz+1 for
    a 'node'-r component but only nz for a 'mid'-r component (see _KIND) --
    write_cvector then copies that array into file-slots 1..N (N=nz+1 or nz)
    UNCHANGED, always leaving file-slot N+1..nz+1 (if any) as zero-padding.
    So the correct remap is native_slot -> (N+1-native_slot) where N is THAT
    component's own r-axis size, NOT a single shared (nz+1) for every
    component -- using the wrong N here previously caused an off-by-one that,
    for a 'mid'-r component (e.g. E_phi/E%x on a FACE grid), could silently
    read an adjacent, wrong-signed layer instead of the intended one (found
    2026-07-26 while investigating a reported E_phi sign discrepancy under
    OUTPUT_CONVENTION=EGBERTKELBERT2012 on a regional grid)."""
    if not r_reversed:
        return r_index
    n = nz if r_kind == 'mid' else (nz + 1)
    return n + 1 - r_index


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


def to_latlon(field, r_index, scaling, factor, th_deg=None, ph_deg=None):
    """
    Maps field's raw phi/theta/r component arrays (from read_field(),
    ALREADY in whatever convention FWD1D.f90 produced -- no relabeling
    applied here or anywhere else in this script) to (lon_deg, lat_deg,
    data) triples ready for pcolormesh.

    The only convention-dependent adjustment made here is to the GRID AXIS
    interpretation: if field['convention'] records theta_convention=
    'LAT_S2N' (i.e. output_convention.f90's apply_output_convention
    reversed the theta index relative to the .grd file's native colatitude
    North->South ordering), th_deg_full is reversed here BEFORE computing
    node/mid-cell positions, so array index j still maps to the correct
    physical node -- verified algebraically to reproduce the exact node/
    mid-cell position a direct (un-reversed-array) derivation gives, for
    both node- and mid-theta components. Similarly, if r_convention differs
    from NATIVE_R_CONVENTION (both solvers' native value, 'R_UP'), r_index
    is remapped so it still selects the same physical layer counting from
    the top/air side. If field['convention'] is None (older file, no
    metadata), NO reversal is applied -- reproduces the original,
    pre-Part-2 axis behavior exactly.

    KNOWN GAP (pre-existing, not introduced here): r_index is applied
    uniformly to all three components' k-axis, but (mirroring the
    node/mid-theta distinction already handled for theta/phi) one of the
    three is actually a "mid-r" quantity on a given EDGE/FACE gridType,
    with its own zero-padding slot -- unlike theta/phi, this was never
    given the same mid/node-aware treatment, before or after this change.

    th_deg/ph_deg (optional): the ACTUAL, NATIVE node positions (degrees)
    for this grid (e.g. from read_grid_file()) -- always pass the grid
    file's own colatitude/top-to-bottom ordering here, NEVER pre-reversed;
    this function performs any needed reversal itself. If omitted (both
    None), falls back to the ORIGINAL assumption of a global grid uniformly
    covering phi=[0,360], theta=[0,180] (np.linspace) -- correct for e.g.
    global.1.0x1.0.grd, WRONG for a regional grid.

    Returns a dict with keys 'phi','theta','r' (the field's own raw
    component naming -- unrelabeled), each an independent (lon_deg,
    lat_deg, data) tuple (mid-staggered components have a genuinely offset
    grid from node-staggered ones, see module docstring).
    """
    nx, ny, nz = field['nx'], field['ny'], field['nz']
    conv = field.get('convention')
    regional = th_deg is not None or ph_deg is not None
    if regional and (th_deg is None or ph_deg is None):
        raise ValueError("to_latlon: th_deg and ph_deg must be given together")

    if not regional:
        ph_deg_full = np.linspace(0.0, 360.0, nx + 1)
        th_deg_full = np.linspace(0.0, 180.0, ny + 1)     # colatitude, native N->S
    else:
        ph_deg_full = np.asarray(ph_deg, dtype=float)
        th_deg_full = np.asarray(th_deg, dtype=float)
        if len(ph_deg_full) != nx + 1 or len(th_deg_full) != ny + 1:
            raise ValueError(
                f"to_latlon: grid size mismatch -- field file has nx={nx},ny={ny} "
                f"(expects {nx+1} phi / {ny+1} theta nodes), grid gave "
                f"{len(ph_deg_full)} phi / {len(th_deg_full)} theta nodes"
            )

    theta_reversed = bool(conv) and conv.get('theta') == 'LAT_S2N'
    if theta_reversed:
        th_deg_full = th_deg_full[::-1]

    r_reversed = bool(conv) and conv.get('r') is not None and conv.get('r') != NATIVE_R_CONVENTION

    # "does phi genuinely wrap the full sphere" -- checked from the actual
    # extent rather than just "was --grid given", so an explicitly-supplied
    # GLOBAL .grd file (e.g. global.1.0x1.0.grd via --grid) is still handled
    # as global (drops the true 360deg-duplicate node-phi column) rather
    # than being misclassified as a non-wrapping regional patch.
    is_global = (ph_deg_full[-1] - ph_deg_full[0]) > 359.99

    gridType = field['gridType']
    if gridType not in _KIND:
        raise ValueError(f"to_latlon: unknown gridType {gridType!r}")
    kind = _KIND[gridType]

    # Each component gets its OWN r-index -- theta/phi/r kinds can differ in
    # r-staggering too (see _KIND/_r_index_for_component), not just phi/theta.
    r_idx_theta = _r_index_for_component(r_index, nz, kind['theta'][2], r_reversed)
    r_idx_phi = _r_index_for_component(r_index, nz, kind['phi'][2], r_reversed)
    r_idx_r = _r_index_for_component(r_index, nz, kind['r'][2], r_reversed)

    Ftheta = field['theta'][:, :, r_idx_theta - 1] * scaling * factor
    Fphi = field['phi'][:, :, r_idx_phi - 1] * scaling * factor
    Fr = field['r'][:, :, r_idx_r - 1] * scaling * factor

    lon_t, lat_t, Dtheta = _component_latlon(Ftheta, ph_deg_full, th_deg_full, *kind['theta'][:2], is_global)
    lon_p, lat_p, Dphi = _component_latlon(Fphi, ph_deg_full, th_deg_full, *kind['phi'][:2], is_global)
    lon_r, lat_r, Dr = _component_latlon(Fr, ph_deg_full, th_deg_full, *kind['r'][:2], is_global)
    return {'theta': (lon_t, lat_t, Dtheta), 'phi': (lon_p, lat_p, Dphi), 'r': (lon_r, lat_r, Dr)}


def _plot(panels, suptitle, fname, xlabel, ylabel, invert_y):
    """panels: list of (lon_or_phi_deg, lat_or_theta_deg, data, title) --
    EACH panel carries its OWN grid now (see to_latlon's docstring: a
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

    # Convention metadata, extracted directly from the field files' own
    # header lines (read_field() -> _parse_convention_header()) -- no
    # command-line toggle needed any more, since the file itself records
    # what convention it's in. E and H files from the same FWD1D run always
    # carry identical convention metadata; fall back to whichever is present.
    conv = Efield['convention'] or Hfield['convention']
    if conv:
        conv_name = conv.get('OUTPUT_CONVENTION', '?')
        conv_desc = (f"{conv_name}  [solver={conv.get('SOLVER','?')} "
                     f"time={conv.get('time','?')} norm={conv.get('norm','?')} "
                     f"theta={conv.get('theta','?')} r={conv.get('r','?')}]")
        # Short form for the figure title (full conv_desc doesn't fit --
        # printed to the console instead, see below).
        conv_title = f"{conv.get('SOLVER','?')}, {conv_name}"
        theta_reversed = conv.get('theta') == 'LAT_S2N'
        print(f"Convention (from file header): {conv_desc}")
    else:
        conv_name = None
        conv_desc = "unknown (no convention metadata in file header -- older file?)"
        conv_title = "convention unknown"
        theta_reversed = False
        print(f"No convention metadata found in file header -- {conv_desc}; "
              f"axes assumed native (colatitude North->South, no reversal applied)")

    OUT_E_PNG, OUT_H_PNG = output_png_names(args.efield_file, args.hfield_file, conv_name, R_INDEX)

    E = to_latlon(Efield, R_INDEX, SCALING, FACTOR, th_deg=th_deg, ph_deg=ph_deg)
    H = to_latlon(Hfield, R_INDEX, SCALING, FACTOR, th_deg=th_deg, ph_deg=ph_deg)

    # Raw phi/theta/r component slots, plotted exactly as read -- NO
    # relabeling of any kind (see module docstring): the file is already in
    # whatever convention conv_desc above records.
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
    e_names, h_names = ("Etheta", "Ephi"), ("Hr", "Htheta", "Hphi")
    xlabel, ylabel = "longitude (deg)", "latitude (deg)"
    # invert_y: put north at the top of the plot. When theta was reversed
    # (LAT_S2N), latitude already ASCENDS with array index (south->north),
    # matplotlib's default already puts north at top -- no invert needed.
    # When not reversed (native colatitude N->S, or no metadata), latitude
    # DESCENDS with array index -- invert to put north at top, matching the
    # pre-Part-2 script's tested behavior for this case.
    invert_y = not theta_reversed

    # Short title only (solver + convention name) -- the full dimension
    # breakdown (conv_desc) is printed to the console above instead; it's too
    # long to fit in the figure title.
    _plot(e_panels, f"E field, T={PERIOD}s, r-index={R_INDEX} -- {conv_title}",
          OUT_E_PNG, xlabel, ylabel, invert_y)
    _plot(h_panels, f"H field, T={PERIOD}s, r-index={R_INDEX} -- {conv_title}",
          OUT_H_PNG, xlabel, ylabel, invert_y)

    for name, (_, _, data, _) in zip(e_names, e_panels):
        print(f"max|{name}|={np.max(np.abs(data)):.3e}", end="  ")
    print()
    for name, (_, _, data, _) in zip(h_names, h_panels):
        print(f"max|{name}|={np.max(np.abs(data)):.3e}", end="  ")
    print()
