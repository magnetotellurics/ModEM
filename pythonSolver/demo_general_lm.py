"""
General (l, m) demo -- generalizes demo_l1_mneg1.py from a single (l=1, m=-1)
term to a superposition of any number of (K0, l, m) source terms, and
produces TWO figures from one shared setup: the induced E field
(general_lm_pattern_e.png) and the induced H field (general_lm_pattern_h.png).

KEY POINT (unchanged from demo_l1_mneg1.py): the radial solve (solve_layered /
closed_form_single_sphere) takes only `l` as an argument -- it is completely
insensitive to m, because l(l+1) is the eigenvalue of the angular Laplacian
for EVERY m. So each (l, m) term in the source list below gets its own radial
solve keyed only on its l (solves are cached so repeated l's aren't re-solved),
and the different m's only change the ANGULAR shape of the source current and
of the resulting fields, via C_lm/Y_lm -- exactly as in demo_l1_mneg1.py, just
summed over multiple terms the way demo_superposition.py sums over multiple l
(m=0) terms.

Fields (and the source pattern) from each (K0, l, m) term are computed at the
SAME evaluation radius and summed coherently (complex superposition) before
plotting -- this reproduces demo_l1_mneg1.py exactly when
SOURCE_TERMS = [(1.0, 1, -1)].
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from spherical_em_induction import solve_layered, C_lm, fields_from_R_general

# ----------------------------------------------------------------------
# Configuration: Earth model, evaluation radius, and the (K0, l, m) source
# terms to superpose. Add/remove/edit entries in SOURCE_TERMS freely --
# everything below (both figures) adapts automatically.
# ----------------------------------------------------------------------

mu0 = 4 * np.pi * 1e-7
mu_earth = mu0
sigma = 1e-2                 # S/m, uniform mantle-like conductor
a = 6.371e6                  # Earth radius, m
c = a + 1.1e5                # ~110 km altitude source shell
T_period = 1000
omega = 2 * np.pi / T_period
r_eval = a                   # ground level

# If True, plot E as Ex (South->North, = -Etheta), Ey (West->East, = +Ephi),
# and H as Hx (=-Htheta), Hy (=+Hphi), Hz (down, =-Hr), on a
# latitude x [-180,180] longitude grid -- MT convention, matches
# PrimaryField.m's sign convention exactly. If False, plot Etheta/Ephi and
# Hr/Htheta/Hphi directly on a colatitude x [0,360] grid.
PLOT_MT_CONVENTION = True

SOURCE_TERMS = [
    # (K0, l, m)
    (1.0, 1, 0),
    #(-0.5, 1, 1),
    #(0.5, 1, -1),
    #(0.4, 2, 0),
    #(0.3, 2, 2),
]

# ----------------------------------------------------------------------
# Angular grid. theta (colatitude, 0..pi) is always used for the physics
# (Y_lm/S_l/C_lm need colatitude). phi spans [-pi,pi) under MT convention
# (-> longitude [-180,180] directly, no wraparound reordering needed) or
# [0,2pi) otherwise.
# ----------------------------------------------------------------------
theta = np.linspace(0.02, np.pi - 0.02, 90)
if PLOT_MT_CONVENTION:
    phi = np.linspace(-np.pi, np.pi, 180)
else:
    phi = np.linspace(0, 2 * np.pi, 180)
TH, PH = np.meshgrid(theta, phi, indexing='ij')

if PLOT_MT_CONVENTION:
    X = np.degrees(PH)                    # longitude, [-180,180]
    Y = 90.0 - np.degrees(TH)              # latitude, [-90,90], +90=N pole
    xlabel, ylabel = "longitude (deg)", "latitude (deg)"
else:
    X = np.degrees(PH)                    # longitude, [0,360]
    Y = np.degrees(TH)                    # colatitude, [0,180], 0=N pole
    xlabel, ylabel = "longitude / phi (deg)", "colatitude / theta (deg)"

# ----------------------------------------------------------------------
# Radial solves, cached by l (m does not affect the radial problem at all --
# see module docstring above / spherical_em_induction.py's own comments)
# ----------------------------------------------------------------------
degrees = sorted({l for _, l, _ in SOURCE_TERMS})
radial_solutions = {}
for l in degrees:
    sol = solve_layered(l, [a], [sigma], [mu0], mu0, c, 1.0, omega)  # K0=1 here;
    Rval = sol['R'](np.array([r_eval]))[0]                            # actual K0 applied
    Rpval = sol['Rp'](np.array([r_eval]))[0]                          # per-term below
    radial_solutions[l] = (Rval, Rpval)

# ----------------------------------------------------------------------
# Superpose source pattern and induced E, H fields over all (K0, l, m) terms
# ----------------------------------------------------------------------
Ctheta_total = np.zeros_like(TH, dtype=complex)
Cphi_total = np.zeros_like(TH, dtype=complex)
Etheta_total = np.zeros_like(TH, dtype=complex)
Ephi_total = np.zeros_like(TH, dtype=complex)
Hr_total = np.zeros_like(TH, dtype=complex)
Htheta_total = np.zeros_like(TH, dtype=complex)
Hphi_total = np.zeros_like(TH, dtype=complex)

for K0, l, m in SOURCE_TERMS:
    Rval, Rpval = radial_solutions[l]

    Ctheta, Cphi = C_lm(l, m, TH, PH)
    Ctheta_total += K0 * Ctheta
    Cphi_total += K0 * Cphi

    fields = fields_from_R_general(l, m, r_eval, TH, PH, K0 * Rval, K0 * Rpval, mu0, omega)
    Etheta_total += fields['Etheta']
    Ephi_total += fields['Ephi']
    Hr_total += fields['Hr']
    Htheta_total += fields['Htheta']
    Hphi_total += fields['Hphi']

# MT-convention fields: Ex=-Etheta, Ey=+Ephi, Hx=-Htheta, Hy=+Hphi, Hz=-Hr
# (matches PrimaryField.m's sign convention exactly).
Ex_total = -Etheta_total
Ey_total = Ephi_total
Hx_total = -Htheta_total
Hy_total = Hphi_total
Hz_total = -Hr_total

terms_str = ", ".join(f"(K0={K0},l={l},m={m})" for K0, l, m in SOURCE_TERMS)


def _plot(panels, suptitle, fname):
    n = len(panels)
    fig, axes = plt.subplots(n, 2, figsize=(10, 4 * n))
    for row, (data_re, data_im, title) in enumerate(panels):
        for col, (data, ri) in enumerate([(data_re, "Re"), (data_im, "Im")]):
            ax = axes[row, col]
            im = ax.pcolormesh(X, Y, data, shading='auto', cmap='RdBu_r')
            ax.set_title(f"{title} ({ri})", fontsize=10)
            plt.colorbar(im, ax=ax)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            if not PLOT_MT_CONVENTION:
                ax.invert_yaxis()   # colatitude: put theta=0 (N pole) at top
            # latitude already increases upward -> no inversion needed for MT convention
    plt.suptitle(f"{suptitle}\n{terms_str}", fontsize=10)
    plt.tight_layout()
    plt.savefig(fname, dpi=140)
    print(f"Saved -> {fname}")


# ----------------------------------------------------------------------
# FIGURE 1: electric field (4 rows x 2 cols: Ktheta, Etheta/Ex, Kphi, Ephi/Ey)
# ----------------------------------------------------------------------
if PLOT_MT_CONVENTION:
    e_panels = [
        (Ctheta_total.real, Ctheta_total.imag, r"source $K_\theta$, superposed"),
        (Ex_total.real, Ex_total.imag, rf"induced $E_x$ (S$\to$N) at $r$={r_eval/a:.2g}$a$"),
        (Cphi_total.real, Cphi_total.imag, r"source $K_\phi$, superposed"),
        (Ey_total.real, Ey_total.imag, rf"induced $E_y$ (W$\to$E) at $r$={r_eval/a:.2g}$a$"),
    ]
else:
    e_panels = [
        (Ctheta_total.real, Ctheta_total.imag, r"source $K_\theta$, superposed"),
        (Etheta_total.real, Etheta_total.imag, rf"induced $E_\theta$ at $r$={r_eval/a:.2g}$a$"),
        (Cphi_total.real, Cphi_total.imag, r"source $K_\phi$, superposed"),
        (Ephi_total.real, Ephi_total.imag, rf"induced $E_\phi$ at $r$={r_eval/a:.2g}$a$"),
    ]
_plot(e_panels, "Superposed toroidal source and induced E field pattern", "general_lm_pattern_e.png")

# ----------------------------------------------------------------------
# FIGURE 2: magnetic field (3 rows x 2 cols: Hr/Hz, Htheta/Hx, Hphi/Hy)
# ----------------------------------------------------------------------
if PLOT_MT_CONVENTION:
    h_panels = [
        (Hz_total.real, Hz_total.imag, rf"induced $H_z$ (down) at $r$={r_eval/a:.2g}$a$"),
        (Hx_total.real, Hx_total.imag, rf"induced $H_x$ (S$\to$N) at $r$={r_eval/a:.2g}$a$"),
        (Hy_total.real, Hy_total.imag, rf"induced $H_y$ (W$\to$E) at $r$={r_eval/a:.2g}$a$"),
    ]
else:
    h_panels = [
        (Hr_total.real, Hr_total.imag, rf"induced $H_r$ at $r$={r_eval/a:.2g}$a$"),
        (Htheta_total.real, Htheta_total.imag, rf"induced $H_\theta$ at $r$={r_eval/a:.2g}$a$"),
        (Hphi_total.real, Hphi_total.imag, rf"induced $H_\phi$ at $r$={r_eval/a:.2g}$a$"),
    ]
_plot(h_panels, "Superposed induced H field pattern", "general_lm_pattern_h.png")

print(f"Source terms: {SOURCE_TERMS}")
if PLOT_MT_CONVENTION:
    print(f"max|Ex| = {np.max(np.abs(Ex_total)):.3e}   max|Ey| = {np.max(np.abs(Ey_total)):.3e}")
    print(f"max|Hx| = {np.max(np.abs(Hx_total)):.3e}   max|Hy| = {np.max(np.abs(Hy_total)):.3e}   max|Hz| = {np.max(np.abs(Hz_total)):.3e}")
else:
    print(f"max|Etheta| = {np.max(np.abs(Etheta_total)):.3e}   max|Ephi| = {np.max(np.abs(Ephi_total)):.3e}")
    print(f"max|Hr| = {np.max(np.abs(Hr_total)):.3e}   max|Htheta| = {np.max(np.abs(Htheta_total)):.3e}   max|Hphi| = {np.max(np.abs(Hphi_total)):.3e}")
