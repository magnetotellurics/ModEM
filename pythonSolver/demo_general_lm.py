"""
General (l, m) demo -- generalizes demo_l1_mneg1.py from a single (l=1, m=-1)
term to a superposition of any number of source terms, and produces TWO
figures from one shared setup: the induced E field (general_lm_pattern_e.png)
and the induced H field (general_lm_pattern_h.png).

INPUT FORMAT: source coefficients are given in GLOBAL1D_COEFF, the SAME
packed format global1d's .prm files / field1d.f90's and field1d_sunegbert2012.f90's
`coeff(:)` arrays use: for degree l=0..LMAX, a block of (2l+1) complex
values ordered m=0,+1,-1,+2,-2,...,+l,-l (l=0's single entry is always
unused -- no monopole term). Paste directly from a .prm file's real/imag
columns, or from a Fortran `coeff(:)` array literal, in this same order --
see unpack_global1d_coeff() below for the exact index formula.

OUTPUT CONVENTION: all computed fields (and the plotted source pattern) are
in field1d.f90/field1d_sunegbert2012.f90's NATIVE e^{-iwt} time convention, directly
comparable to their raw .hfield/.efield output with no further conjugation
needed. This requires THREE corrections on top of pythonSolver's own native
e^{+iwt} machinery (solve_layered / fields_from_R_general / C_lm), all
found and verified 2026-07-23 while cross-validating field1d_sunegbert2012.f90 (see
CLAUDE.md, "Cross-convention comparison rule", and
pythonSolver/reference_earth_l1mneg1.py for the original derivation/check):
  1. normalize by sqrt(l(l+1)) -- pythonSolver's own convention is
     T=R(r)*Y_l^m/sqrt(l(l+1)), field1d_sunegbert2012's is T_l(r)*Y_l^m directly (no
     norm division), so this factor is needed on TOP of solve_layered's own
     "A" (external-amplitude) normalization;
  2. flip the sign of m -- conjugating a field also conjugates its e^{imphi}
     angular dependence, turning an (l,+m) pattern into an (l,-m)-like one,
     so field1d's (l,+m) term corresponds to conj(pythonSolver's (l,-m)
     term), NOT conj(pythonSolver's (l,+m) term);
  3. H picks up an extra overall -1 that E (and, by the same "polar vector"
     reasoning, the source current K) does not -- H is a pseudovector
     (reverses under true time-reversal, being sourced by currents), E is a
     polar vector (does not). NOTE: the -1-for-H-only rule was empirically
     verified for E/H fields specifically (l=1,m=-1, Earth-scale, see
     CLAUDE.md); applying the SAME (E-type, +1) rule to the source pattern
     Ktheta/Kphi below is a reasonable extension (same "polar vector" type
     as E) but has not been independently re-verified on its own.

COEFFICIENT SCALING: none -- each raw GLOBAL1D_COEFF value is used as-is, no
extra rescaling. (A Period/5 rescaling, matching MATLAB's PrimaryField.m
"shc" scaling before TSModel.ShcInc, was tried 2026-07-23 and removed
2026-07-25 to keep this demo consistent with FWD1D.f90/field1d_sunegbert2012.f90's
own driver-level behaviour, which applies no such rescaling either; see
CLAUDE.md, "Matching MATLAB/SIEM's practical coefficient scaling" for the
historical note.)

KEY POINT (unchanged from demo_l1_mneg1.py): the radial solve (solve_layered)
takes only `l` as an argument -- it is completely insensitive to m, because
l(l+1) is the eigenvalue of the angular Laplacian for EVERY m. So each
(l, m) term gets its own radial solve keyed only on its l (cached, so
repeated l's aren't re-solved), and the different m's only change the
ANGULAR shape of the source current and of the resulting fields.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from spherical_em_induction import solve_layered, C_lm, fields_from_R_general

# ----------------------------------------------------------------------
# Configuration: Earth model and evaluation radius.
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

# ----------------------------------------------------------------------
# Source coefficients, global1d packed format -- see module docstring.
# Default below reproduces the old demo's default, SOURCE_TERMS=[(0.5,1,-1)]:
# LMAX=1, l=1 block (indices 1..3) = [m=0, m=+1, m=-1] = [0, 0, 0.5].
# ----------------------------------------------------------------------
LMAX = 1
GLOBAL1D_COEFF = [
    0.0 + 0.0j,   # l=0, m=0  (unused -- no monopole)
    0.0 + 0.0j,   # l=1, m=0
    0.0 + 0.0j,   # l=1, m=+1
    0.5 + 0.0j,   # l=1, m=-1
]


def unpack_global1d_coeff(coeff, lmax):
    """Unpack a global1d-format packed coeff(:) array into a list of
    (coeff, l, m) tuples (m=-l..l, zero entries dropped). Index formula
    (0-indexed here; matches field1d_sunegbert2012.f90's 1-indexed
    coeff(l*l+1), coeff(l*l+2*m), coeff(l*l+2*m+1) exactly, offset by 1):
        coeff2(l, 0)  = coeff[l*l]
        coeff2(l, +m) = coeff[l*l + 2*m - 1]   for m=1..l
        coeff2(l, -m) = coeff[l*l + 2*m]       for m=1..l
    """
    expected_len = (lmax + 1) ** 2
    if len(coeff) != expected_len:
        raise ValueError(f"GLOBAL1D_COEFF has length {len(coeff)}, expected "
                          f"(LMAX+1)^2={expected_len} for LMAX={lmax}")
    terms = []
    for l in range(1, lmax + 1):
        base = l * l
        if coeff[base] != 0:
            terms.append((coeff[base], l, 0))
        for m in range(1, l + 1):
            if coeff[base + 2 * m - 1] != 0:
                terms.append((coeff[base + 2 * m - 1], l, m))
            if coeff[base + 2 * m] != 0:
                terms.append((coeff[base + 2 * m], l, -m))
    return terms


SOURCE_TERMS = unpack_global1d_coeff(GLOBAL1D_COEFF, LMAX)

# ----------------------------------------------------------------------
# Angular grid. theta (colatitude, 0..pi) is always used for the physics
# (Y_lm/C_lm need colatitude). phi spans [-pi,pi) under MT convention
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
# see module docstring above / spherical_em_induction.py's own comments).
# Solved once at K0=1 (pythonSolver's own convention); A is solve_layered's
# own external-amplitude readback for that K0=1 solve, used below to
# normalize to unit external amplitude before rescaling by the user's actual
# (global1d-convention) coefficient.
# ----------------------------------------------------------------------
degrees = sorted({l for _, l, _ in SOURCE_TERMS})
radial_solutions = {}
for l in degrees:
    sol = solve_layered(l, [a], [sigma], [mu0], mu0, c, 1.0, omega)
    Rval1 = sol['R'](np.array([r_eval]))[0]
    Rpval1 = sol['Rp'](np.array([r_eval]))[0]
    A1 = sol['A']
    radial_solutions[l] = (Rval1, Rpval1, A1)

# ----------------------------------------------------------------------
# Superpose source pattern and induced E, H fields over all source terms,
# converting each term to field1d_sunegbert2012's native e^{-iwt} convention via the
# norm + m-flip + H-sign correction described in the module docstring.
# ----------------------------------------------------------------------
Ctheta_total = np.zeros_like(TH, dtype=complex)
Cphi_total = np.zeros_like(TH, dtype=complex)
Etheta_total = np.zeros_like(TH, dtype=complex)
Ephi_total = np.zeros_like(TH, dtype=complex)
Hr_total = np.zeros_like(TH, dtype=complex)
Htheta_total = np.zeros_like(TH, dtype=complex)
Hphi_total = np.zeros_like(TH, dtype=complex)

for coeff_raw, l, m in SOURCE_TERMS:
    Rval1, Rpval1, A1 = radial_solutions[l]
    norm = np.sqrt(l * (l + 1))
    m_python = -m   # correction 2: m-flip

    coeff = coeff_raw   # no extra rescaling -- see module docstring, "COEFFICIENT SCALING"

    fields_py = fields_from_R_general(l, m_python, r_eval, TH, PH, Rval1, Rpval1, mu0, omega)
    Ctheta_py, Cphi_py = C_lm(l, m_python, TH, PH)

    scale_E = coeff * norm / A1          # correction 1 (norm) folded in; H sign is +1 for E/K
    scale_H = -scale_E                    # correction 3: extra -1 for H only

    Ctheta_total += scale_E * np.conj(Ctheta_py)
    Cphi_total += scale_E * np.conj(Cphi_py)
    Etheta_total += scale_E * np.conj(fields_py['Etheta'])
    Ephi_total += scale_E * np.conj(fields_py['Ephi'])
    Hr_total += scale_H * np.conj(fields_py['Hr'])
    Htheta_total += scale_H * np.conj(fields_py['Htheta'])
    Hphi_total += scale_H * np.conj(fields_py['Hphi'])

# MT-convention fields: Ex=-Etheta, Ey=+Ephi, Hx=-Htheta, Hy=+Hphi, Hz=+Hr
# (matches PrimaryField.m's sign convention exactly; these relationships are
# convention-independent, unaffected by the e^{-iwt}/e^{+iwt} correction above.
# NOTE: Hz=+Hr, NOT -Hr -- "down" is already the natural sign of Hr in
# global1d's own convention, confirmed directly by the user and applied the
# same way in plot_global1d_output.py's Fz=Fr.T; an earlier version of this
# line had the extra flip, which showed up as a sign mismatch on Hz alone
# once compared directly against field1d_sunegbert2012 output, 2026-07-23).
Ex_total = -Etheta_total
Ey_total = Ephi_total
Hx_total = -Htheta_total
Hy_total = Hphi_total
Hz_total = Hr_total

terms_str = ", ".join(f"(coeff={coeff},l={l},m={m})" for coeff, l, m in SOURCE_TERMS)


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
    plt.suptitle(f"{suptitle} (native e^-iwt, field1d convention)\n{terms_str}", fontsize=10)
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

print(f"Source terms (global1d/coeff convention, LMAX={LMAX}): {SOURCE_TERMS}")
if PLOT_MT_CONVENTION:
    print(f"max|Ex| = {np.max(np.abs(Ex_total)):.3e}   max|Ey| = {np.max(np.abs(Ey_total)):.3e}")
    print(f"max|Hx| = {np.max(np.abs(Hx_total)):.3e}   max|Hy| = {np.max(np.abs(Hy_total)):.3e}   max|Hz| = {np.max(np.abs(Hz_total)):.3e}")
else:
    print(f"max|Etheta| = {np.max(np.abs(Etheta_total)):.3e}   max|Ephi| = {np.max(np.abs(Ephi_total)):.3e}")
    print(f"max|Hr| = {np.max(np.abs(Hr_total)):.3e}   max|Htheta| = {np.max(np.abs(Htheta_total)):.3e}   max|Hphi| = {np.max(np.abs(Hphi_total)):.3e}")
