"""
Reference values for test_unit_sphere.f90.

Single homogeneous conducting sphere, radius r0, conductivity sigma, in
vacuum, driven by a UNIT external (l,m) multipole (Theta_ext(r)=r^(l+1),
amplitude 1). Computes the field components using Sun & Egbert (2012)
eq.(5)-(6), independently derived (see field1d_sunegbert2012.f90's module
docstring for the derivation).

Physical parameters match test_unit_sphere.f90 EXACTLY:
  r0 = 1000 m, sigma = 1 S/m, omega = 1 rad/s, mu0 = 1.256637e-6 (field1d.f90's
  literal), pi = 3.14159265357898 (field1d.f90's literal), l=2, m=1.

Merged 2026-07-25 from the former separate reference_unit_sphere.py (which
compared, indirectly via a s=-1/conjugate dance, against field1d.f90 only) and
reference_unit_sphere_sunegbert2012.py (which compared, directly and exactly,
against field1d_sunegbert2012.f90) into ONE script with TWO clearly labeled
sections, since both use the identical closed-form machinery below -- only
the pass/fail criterion differs between the two solvers (see each section's
own docstring for why).

IMPORTANT: field1d.f90/field1d_sunegbert2012.f90 are BOTH native e^{-iwt}
(kl=sqrt(i*omega*mu0*sigma)), so s=+1 below is directly comparable to either
solver's raw output with no extra time-convention conjugation.
"""
import numpy as np
from scipy.special import jv, lpmv, factorial

# ---- exact literals from field1d.f90 ----
mu0 = 1.256637e-6
pi_f = 3.14159265357898

r0 = 1000.0
sigma = 1.0
omega = 1.0
l, m = 2, 1

pi = pi_f
d2r = pi / 180.0


def sph_jn_complex(l, z):
    return np.sqrt(pi / (2 * z)) * jv(l + 0.5, z)


def rj(l, k, r):
    if k == 0:
        return r ** (l + 1)
    return r * sph_jn_complex(l, k * r)


def rj_prime(l, k, r, h=1e-7):
    if k == 0:
        return (l + 1) * r ** l
    hr = h * r
    return (rj(l, k, r + hr) - rj(l, k, r - hr)) / (2 * hr)


def solve_unit_external(l, k1, a):
    """Match interior C*psi_l(k1 r) to exterior r^(l+1) + B*r^-l with unit
    external amplitude (coefficient of r^(l+1) = 1), via continuity of
    Theta and Theta' at r=a (eq 12; no thin sheet, tau=0).

    Derivation: exterior Text(r)=r^(l+1)+B*r^-l, Text'(r)=(l+1)*r^l-l*B*r^-(l+1);
    interior Tint(r)=C*psi_l(k1 r), Tint'(r)=C*Fap (Fap=d/dr[psi_l(k1 r)] at r=a).
    Continuity of T at a:  -B*a^-l + C*Fa  = a^(l+1)
    Continuity of T' at a: (l+1)*a^l - l*B*a^-(l+1) = C*Fap
                            => l*a^-(l+1)*B + C*Fap = (l+1)*a^l   (note: +Fap,
    not -Fap -- an earlier version of this function had a sign error here,
    found 2026-07-23 by direct comparison against field1d_sunegbert2012.f90's own
    eq(23)-based radial solve, which -- independently -- matches to ~1e-9
    once this sign is corrected)."""
    Fa, Fap = rj(l, k1, a), rj_prime(l, k1, a)
    M = np.array([[-a ** (-l), Fa],
                  [l * a ** (-l - 1), Fap]], dtype=complex)
    rhs = np.array([a ** (l + 1), (l + 1) * a ** l], dtype=complex)
    B, C = np.linalg.solve(M, rhs)
    return B, C


def Theta_ext_region(r, l, B):
    Th = r ** (l + 1) + B * r ** (-l)
    Thp = (l + 1) * r ** l - l * B * r ** (-l - 1)
    return Th, Thp


def Y_lm(l, m, theta, phi):
    norm = np.sqrt((2 * l + 1) / (4 * pi) * factorial(l - m) / factorial(l + m))
    return norm * lpmv(m, l, np.cos(theta)) * np.exp(1j * m * phi)


def dY_dtheta(l, m, theta, phi, h=1e-7):
    return (Y_lm(l, m, theta + h, phi) - Y_lm(l, m, theta - h, phi)) / (2 * h)


def fields(Th, Thp, r, theta, phi, l, m, s):
    """s=+1: e^{-iwt} (field1d.f90/field1d_sunegbert2012.f90 native). s=-1: e^{+iwt}."""
    Yv = Y_lm(l, m, theta, phi)
    Yt = dY_dtheta(l, m, theta, phi)
    Yp = (1j * m / np.sin(theta)) * Yv
    Hr = (l * (l + 1) / r ** 2) * Th * Yv
    Htheta = (1 / r) * Thp * Yt
    Hphi = (1 / r) * Thp * Yp
    Etheta = s * (1j * omega * mu0 / r) * Th * Yp
    Ephi = -s * (1j * omega * mu0 / r) * Th * Yt
    return dict(Hr=Hr, Htheta=Htheta, Hphi=Hphi, Etheta=Etheta, Ephi=Ephi)


# ---- exact grid geometry from test_unit_sphere.f90 ----
th_node5 = (60.0 + 4 * (150.0 - 60.0) / 8) * d2r      # grid%th(5)
dt5 = (150.0 - 60.0) / 8 * d2r
th_mid5 = th_node5 + dt5 / 2                            # grid%th(5)+dt(5)/2
ph_node1 = 0.0                                           # grid%ph(1)
dp1 = 360.0 / 8 * d2r
ph_mid1 = ph_node1 + dp1 / 2                             # grid%ph(1)+dp(1)/2

Rr2 = 1525.0   # m
Rs2 = 2000.0   # m

points = {
    "Hr (H%z)":     (Rr2, th_node5, ph_node1, "Hr"),
    "Hphi (H%x)":   (Rs2, th_node5, ph_mid1,  "Hphi"),
    "Htheta (H%y)": (Rs2, th_mid5,  ph_node1, "Htheta"),
    "Etheta (E%y)": (Rr2, th_node5, ph_mid1,  "Etheta"),
    "Ephi (E%x)":   (Rr2, th_mid5,  ph_node1, "Ephi"),
}

k_ft = np.sqrt(1j * omega * mu0 * sigma)     # native e^{-iwt} k, both solvers
k_py = np.conj(k_ft)                          # e^{+iwt} convention (Section 1 only)
B_ft, C_ft = solve_unit_external(l, k_ft, r0)
B_py, C_py = solve_unit_external(l, k_py, r0)

print(f"k (native e^-iwt) = {k_ft:.8e},  |k*r0| = {abs(k_ft)*r0:.6f}")
print(f"B (external-region beta/alpha, native e^-iwt) = {B_ft:.10e}")

# ============================================================================
# Section 1: field1d.f90 (KELBERT2014)
#
# KELBERT2014 carries its own internal REAL normalization constants per
# component (R0^2/r^n factors, an extra 1/(l(l+1)) on the tangential H and E
# components that H_r does not have) AND applies a final conjg() to every
# assembled component intended to compensate for a "+m,-m conjugate-pairing"
# reconstruction trick that this test (sourcing only m=+1, not the pair) never
# actually exercises as intended. As analyzed 2026-07-25 (see CLAUDE.md), the
# result is NOT a clean scalar/phase discrepancy for m!=0 -- so do NOT expect
# ratio = 1+0i, or even a single common phase across all five components, when
# comparing these numbers against test_unit_sphere.f90's KELBERT2014 output.
# Printed for visual/diagnostic inspection only (kept in its original s=-1
# form: what MUST match, if KELBERT2014's own bookkeeping were internally
# consistent for this term, is the PHASE of field1d_component / conj(this
# s=-1 value) -- that ratio should come out real; per the analysis above and
# in CLAUDE.md, it generally will NOT for m!=0, which is exactly the point).
# ============================================================================
print()
print("=== Section 1: field1d.f90 (KELBERT2014) -- diagnostic only, no clean match expected ===")
print(f"{'component':16s} {'(r,theta,phi)':45s} {'value @ s=-1 (+iwt)':>45s}")
for name, (r, th, ph, which) in points.items():
    Th_py, Thp_py = Theta_ext_region(r, l, B_py)
    F = fields(Th_py, Thp_py, r, th, ph, l, m, s=-1)
    val = F[which]
    print(f"{name:16s} r={r:8.2f} th={th:.10f} ph={ph:.10f}   {val.real:+.10e} {val.imag:+.10e}j")

# ============================================================================
# Section 2: field1d_sunegbert2012.f90 (SUNEGBERT2012)
#
# SUNEGBERT2012 uses this exact formula and this exact "alpha_l^T=1" (unit
# external amplitude) normalization directly -- so the expected match against
# test_unit_sphere.f90's SUNEGBERT2012 output is EXACT (ratio = 1+0i for all
# five components), not just phase-equal.
# ============================================================================
print()
print("=== Section 2: field1d_sunegbert2012.f90 (SUNEGBERT2012) -- expect EXACT match, ratio=1+0i ===")
print(f"{'component':16s} {'(r,theta,phi)':45s} {'value @ native e^-iwt':>45s}")
for name, (r, th, ph, which) in points.items():
    Th, Thp = Theta_ext_region(r, l, B_ft)
    F = fields(Th, Thp, r, th, ph, l, m, s=+1)
    val = F[which]
    print(f"{name:16s} r={r:8.2f} th={th:.10f} ph={ph:.10f}   {val.real:+.10e} {val.imag:+.10e}j")

print()
print("Compare Section 1 against test_unit_sphere.f90's 'KELBERT2014 (field1d.f90)' block")
print("(diagnostic only) and Section 2 against its 'SUNEGBERT2012 (field1d_sunegbert2012.f90)'")
print("block (expect ratio FIELD1D_SUNEGBERT2012/reference = 1+0i).")
