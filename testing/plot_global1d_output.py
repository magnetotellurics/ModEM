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
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
#EFIELD_FILE = r"C:\Users\Anna Kelbert\Developer\ModEM-global1d\MTsource\MT.1000sec.Mode1.fix.E-grid.T01.efield"
#HFIELD_FILE = r"C:\Users\Anna Kelbert\Developer\ModEM-global1d\MTsource\MT.1000sec.Mode1.fix.E-grid.T01.hfield"
EFIELD_FILE = r"C:\Users\Anna Kelbert\Developer\ModEM-global1d\MTsource\MT.1000sec.Mode2.E-grid.T01.efield"
HFIELD_FILE = r"C:\Users\Anna Kelbert\Developer\ModEM-global1d\MTsource\MT.1000sec.Mode2.E-grid.T01.hfield"

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

OUT_E_PNG = "global1d_output_e.png"
OUT_H_PNG = "global1d_output_h.png"


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


def to_latlon_mt(field, r_index, scaling, factor, mt_convention=True):
    """
    field: dict from read_field() (phi/theta/r component arrays, shape
           (nx+1, ny+1, nz+1), axis0=phi(0..360 ascending), axis1=theta
           (colatitude, 0..180 ascending), axis2=r-index)

    If mt_convention: returns (lon_deg [-180,180], lat_deg [-90,90], Fx, Fy,
    Fz) where Fx=South->North (=-Ftheta), Fy=West->East (=+Fphi), Fz=down
    (=+Fr -- global1d's Hr/Er already point down, no extra flip) -- the
    ModEM REGIONAL naming convention, matches PrimaryField.m
    exactly (NOT the original global Uyeshima & Schultz convention
    field1d.f90 was written for -- see module docstring). The redundant
    phi=360deg node is dropped and phi is reindexed to [-180,180] (mirrors
    the MATLAB i180-circular-shift trick, generically for any nx) -- keeping
    it creates a zero-width duplicate column at longitude=0 once reindexed,
    which pcolormesh renders as a seam.

    If not mt_convention: returns (phi_deg [0,360], theta_deg [0,180],
    Ftheta, Fphi, Fr) directly, no sign flips, no reindexing.

    All field arrays shape (n_y, n_x), ready for pcolormesh(x, y, F).
    """
    nx, ny = field['nx'], field['ny']
    phi_deg = np.linspace(0.0, 360.0, nx + 1)
    theta_deg = np.linspace(0.0, 180.0, ny + 1)     # colatitude

    Ftheta = field['theta'][:, :, r_index - 1] * scaling * factor   # (phi, theta)
    Fphi = field['phi'][:, :, r_index - 1] * scaling * factor
    Fr = field['r'][:, :, r_index - 1] * scaling * factor

    if not mt_convention:
        return phi_deg, theta_deg, Ftheta.T, Fphi.T, Fr.T

    # drop the redundant phi=360deg node (physically == phi=0deg; confirmed
    # numerically identical to the phi=0 node to ~1e-14, i.e. floating-point
    # noise) -- keeping both creates a zero-width duplicate column at
    # longitude=0 once reindexed below, which pcolormesh renders as a seam.
    phi_deg = phi_deg[:-1]
    Ftheta = Ftheta[:-1, :]
    Fphi = Fphi[:-1, :]
    Fr = Fr[:-1, :]

    Fx = -Ftheta.T   # transpose -> (theta, phi) = (lat-like, lon-like); South->North
    Fy = Fphi.T       # West->East
    Fz = Fr.T         # down (confirmed 2026-07-22: global1d's Hz/Ez AND Hr/Er
                       # both already point down -- H%z/E%z's leading minus
                       # sign in field1d.f90 is baked into the Mie-potential
                       # derivation itself, not a separate up->down flip, so
                       # no additional sign change is needed here)

    lat_deg = 90.0 - theta_deg   # ascending, -90..90, N pole at theta=0 -> lat=+90

    # reindex phi -> longitude [-180,180], monotonic (mirrors the MATLAB
    # i180-circular-shift trick, generically for any nx)
    split = np.searchsorted(phi_deg, 180.0)
    lon_deg = np.concatenate([phi_deg[split:] - 360.0, phi_deg[:split]])
    Fx = np.concatenate([Fx[:, split:], Fx[:, :split]], axis=1)
    Fy = np.concatenate([Fy[:, split:], Fy[:, :split]], axis=1)
    Fz = np.concatenate([Fz[:, split:], Fz[:, :split]], axis=1)

    return lon_deg, lat_deg, Fx, Fy, Fz


def _plot(panels, suptitle, fname, x_deg, y_deg, xlabel, ylabel, invert_y):
    n = len(panels)
    fig, axes = plt.subplots(n, 2, figsize=(10, 4 * n))
    for row, (data, title) in enumerate(panels):
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
    Efield = read_field(EFIELD_FILE)
    Hfield = read_field(HFIELD_FILE)
    print(f"E file: nx={Efield['nx']} ny={Efield['ny']} nz={Efield['nz']} gridType={Efield['gridType']}")
    print(f"H file: nx={Hfield['nx']} ny={Hfield['ny']} nz={Hfield['nz']} gridType={Hfield['gridType']}")

    x_deg, y_deg, E1, E2, E3 = to_latlon_mt(Efield, R_INDEX, SCALING, FACTOR, PLOT_MT_CONVENTION)
    _, _, H1, H2, H3 = to_latlon_mt(Hfield, R_INDEX, SCALING, FACTOR, PLOT_MT_CONVENTION)

    if PLOT_MT_CONVENTION:
        Ex, Ey = E1, E2
        Hz, Hx, Hy = H3, H1, H2
        e_panels = [
            (Ex, r"$E_x$ (S$\to$N)"),
            (Ey, r"$E_y$ (W$\to$E)"),
        ]
        h_panels = [
            (Hz, r"$H_z$ (down)"),
            (Hx, r"$H_x$ (S$\to$N)"),
            (Hy, r"$H_y$ (W$\to$E)"),
        ]
        xlabel, ylabel, invert_y = "longitude (deg)", "latitude (deg)", False
        e_names, h_names = ("Ex", "Ey"), ("Hx", "Hy", "Hz")
    else:
        Etheta, Ephi = E1, E2
        Ftheta, Fphi, Fr = H1, H2, H3
        e_panels = [
            (Etheta, r"$E_\theta$"),
            (Ephi, r"$E_\phi$"),
        ]
        h_panels = [
            (Fr, r"$H_r$"),
            (Ftheta, r"$H_\theta$"),
            (Fphi, r"$H_\phi$"),
        ]
        xlabel, ylabel, invert_y = "longitude / phi (deg)", "colatitude / theta (deg)", True
        e_names, h_names = ("Etheta", "Ephi"), ("Hr", "Htheta", "Hphi")

    _plot(e_panels, f"global1d output: E field, T={PERIOD}s, r-index={R_INDEX}",
          OUT_E_PNG, x_deg, y_deg, xlabel, ylabel, invert_y)
    _plot(h_panels, f"global1d output: H field, T={PERIOD}s, r-index={R_INDEX}",
          OUT_H_PNG, x_deg, y_deg, xlabel, ylabel, invert_y)

    for name, (data, _) in zip(e_names, e_panels):
        print(f"max|{name}|={np.max(np.abs(data)):.3e}", end="  ")
    print()
    for name, (data, _) in zip(h_names, h_panels):
        print(f"max|{name}|={np.max(np.abs(data)):.3e}", end="  ")
    print()
