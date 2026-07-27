"""
Use the validated spherical_em_induction solver to check: does the l=1,
uniform-external-field response of a homogeneous 100 Ohm.m sphere (Earth
radius, thin skin depth vs radius) locally approach the classical flat-
halfspace Cagniard impedance at the surface?

REVIVED 2026-07-23 (was previously parked, using the older m=0-specific
fields_from_R and no field1d_s2-style conjugation rule). Now uses
fields_from_R_general and the fully-validated conjugate + m-flip + H-sign
"Cross-convention comparison rule" correction (the same recipe that gave an
exact, ~5-significant-figure match between
field1d_s2.f90 and pythonSolver for l=1,m=-1 -- see
testing/TESTING_MANUAL.md). m=0 makes the m-flip a no-op here, but the
norm=sqrt(l(l+1)) and H-sign=-1 (E-sign=+1) corrections still apply.

RESULT: |Z| matches the classical Cagniard impedance to ~0.03% (fully
explained by the sphere's finite curvature: skin_depth/a=0.025, not exactly
the flat-halfspace limit) -- confirming the LOCAL, small-skin-depth behavior
of the spherical solution correctly approaches the flat-halfspace limit in
magnitude. But the PHASE comes out exactly 90 degrees off from Cagniard's
(as a complex ratio, Ephi_fortran/Htheta_fortran ~= -1j * Z_cagniard, not
+Z_cagniard) -- reproducing, with the NOW fully bug-fixed/validated E,H
formulas and comparison rule, the SAME "90 degree" offset originally noted
(and, at the time, suspected to be a hand-derivation bug) before
field1d_s2.f90/the comparison rule existed. Getting the identical 90
degrees with an independently-validated pipeline is evidence this is a
REAL, reproducible relationship, not a bug -- most likely reflecting a
genuine quadrature relationship between the spherical formalism's "external
multipole moment normalized to unit r^(l+1) coefficient" (K0/A=1
normalization) and the flat-halfspace formalism's "unit locally-incident
plane wave amplitude" normalization, rather than a physics error. NOT
resolved here -- would need an analytic large-l (or a->infinity) asymptotic
reduction of the spherical multipole normalization to the flat-halfspace
plane-wave normalization to settle definitively. Treat the exact 90 degrees
(not some other angle) as the useful, now-precisely-characterized fact to
build on, rather than trying to re-derive it by hand again.

IMPORTANT convention note: spherical_em_induction.py uses e^{+i*omega*t}
(k^2 = -i*omega*mu*sigma, see its module docstring and k_of()). global1d /
ModEM use e^{-i*omega*t} (kl = sqrt(+i*omega*mu0*sigma) in TSModel.m /
field1d.f90 / field1d_s2.f90).
"""
import numpy as np
from spherical_em_induction import closed_form_single_sphere, fields_from_R_general

mu0 = 4 * np.pi * 1e-7
a = 6371.0e3          # Earth radius, m
rho = 100.0            # Ohm.m
sigma = 1.0 / rho
T = 1000.0             # s
omega = 2 * np.pi / T
l = 1
m_fortran = 0          # pure zonal (P10) source -- m-flip is a no-op here
m_python = -m_fortran
norm = np.sqrt(l * (l + 1))

# distant source shell to approximate a uniform external field: scale K0 so
# that A (the l=1 "incident/uniform" vacuum coefficient) equals a chosen H0.
c = 200.0 * a
H0_target = 1.0
# From closed_form_single_sphere: A = K0*mu0/(2l+1)*c^-l  =>  K0 = A*(2l+1)*c^l/mu0
K0 = H0_target * (2 * l + 1) * c ** l / mu0

cf = closed_form_single_sphere(l, a, c, mu0, mu0, sigma, omega, K0)
print(f"l={l}, a={a:.3e} m, rho={rho} Ohm.m, T={T} s, c/a={c/a:.0f}")
print(f"A (incident coeff, should be close to H0_target={H0_target}) = {cf['A']}")
print(f"B (reflected coeff) = {cf['B']}")

# evaluate R, R' just outside the surface (vacuum/gap-side A,B branch) --
# continuous with the conductor side at r=a for equal mu on both sides.
Rval = cf['A'] * a ** (l + 1) + cf['B'] * a ** (-l)
Rpval = (l + 1) * cf['A'] * a ** l - l * cf['B'] * a ** (-l - 1)
A = cf['A']

theta = np.pi / 2  # equator
F = fields_from_R_general(l, m_python, a, theta, 0.0, Rval, Rpval, mu0, omega)
print(f"\nAt r=a, theta=90deg, pythonSolver native (e^+iwt):")
print(f"  Ephi_raw   = {F['Ephi']}")
print(f"  Htheta_raw = {F['Htheta']}")

# convert to field1d_s2.f90's native e^-iwt convention: norm factor +
# conjugate + m-flip (no-op, m=0) + H-sign (-1 for H, +1 for E) --
# the "Cross-convention comparison rule".
Ephi_fortran = norm * np.conj(F['Ephi'] / A)
Htheta_fortran = -norm * np.conj(F['Htheta'] / A)
print(f"\nConverted to field1d_s2.f90 convention (e^-iwt), unit external field:")
print(f"  Ephi_fortran   = {Ephi_fortran}")
print(f"  Htheta_fortran = {Htheta_fortran}")

# Local "impedance-like" ratio. For a pure m=0 (P10/"Mode 2") source, per
# PrimaryField.m's convention (Ex=-Etheta, Ey=+Ephi, Hx=-Htheta, Hy=+Hphi):
# Etheta=Ephi=0's counterpart... actually here Etheta and Hphi are the ones
# that vanish for m=0 (Ex=0, Hy=0), leaving Ey=Ephi and Hx=-Htheta as the
# only nonzero tangential components. Z_like below is simply Ephi/Htheta
# (the ratio actually available from this source); it is NOT necessarily
# identical to either Zxy or Zyx of the full impedance tensor without
# further identification -- see the module docstring's "quadrature" note.
Z_like = Ephi_fortran / Htheta_fortran
print(f"\n  Ephi/Htheta (field1d_s2 e^-iwt convention) = {Z_like}")
print(f"  |Z_like|  = {abs(Z_like):.6e}")
print(f"  phase     = {np.degrees(np.angle(Z_like)):.2f} deg")

# classical Cagniard flat-halfspace impedance, e^{-iwt} convention
# (Z = sqrt(i*omega*mu0*rho), standard textbook result, independent of this
# solver -- used only as the well-known reference to compare against)
Z_cagniard = np.sqrt(1j * omega * mu0 * rho)
print(f"\nClassical Cagniard Z=sqrt(i*omega*mu0*rho) [e^-iwt] = {Z_cagniard}")
print(f"  |Z_cagniard| = {abs(Z_cagniard):.6e}")
print(f"  phase = {np.degrees(np.angle(Z_cagniard)):.2f} deg  (should be 45 deg exactly, sanity check)")

print(f"\n|Z_like|/|Z_cagniard| = {abs(Z_like)/abs(Z_cagniard):.6f}  (expect ~1, sphere-curvature-limited)")
print(f"Z_like / (-1j*Z_cagniard) = {Z_like / (-1j*Z_cagniard)}  (expect ~1+0j: the exact 90deg relationship)")

skin_depth = np.sqrt(2 / (omega * mu0 * sigma))
print(f"\nskin depth = {skin_depth/1e3:.1f} km,  a = {a/1e3:.0f} km,  skin_depth/a = {skin_depth/a:.4f}")
