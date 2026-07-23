"""
Independent verification of spherical_em_induction.py's OWN Rval/Rpval
(from solve_layered, K0-normalized to a unit external field via a distant
source shell) against the same closed-form psi_1 formulas already validated
for Fortran's Tnr/Tnsp/Tnrp/Tns in uniform_sphere_Tnr_predict.py -- but this
time using pythonSolver's NATIVE e^{+i*omega*t} convention (k=sqrt(-i*omega*
mu*sigma), via spherical_em_induction.py's own k_of()), not global1d's.

Motivation (see CLAUDE.md, 2026-07-22 entry): comparing global1d against
pythonSolver for l=1,m=-1 shows Hz,Ex,Ey (all built from T(r)/Rval) picking
up an extra factor of (-1i) relative to Hx,Hy (built from T'(r)/Rpval),
after accounting for the e^{+iwt}<->e^{-iwt} conjugation. Fortran's Tnr/Tns
(=T(r)) and Tnsp/Tnrp (=T'(r)) have ALREADY been independently verified
correct (test_Tnr_uniform_sphere.f90 vs uniform_sphere_Tnr_predict.py). This
script runs the SAME kind of check on the pythonSolver side: is Rval (not
Rpval) where the discrepancy actually lives?

Derivation reused (see uniform_sphere_Tnr_predict.py for the full trace):
for a single uniform layer down to r=0, in EITHER time convention,
  A1 = (1 + x*psi_1'(x)/psi_1(x)) / 3,   x = k*r0
  T(r0)/A1 -> the "per unit external field" (A=1 normalized) potential value
  T'(r0)/A1 = (3 - T(r0)/A1) / r0
This algebra never assumed a particular sign of k, so it applies unchanged
to pythonSolver's k = sqrt(-i*omega*mu*sigma).

Normalization: solve_layered's own coefficient "A" (vacuum r^(l+1) term) is
the DIRECT analog of Fortran's A1-normalization reference (unit external
field). Using a distant source shell (c=200*a) with K0 chosen so A~=1 (as
already confirmed numerically in compare_global1d_vs_cagniard.py: "A ~= 1.0"
to high precision for this c/a), Rval/Rpval come out ALREADY in the
A=1-normalized basis -- directly comparable to the psi_1 predictions below,
no extra rescaling needed.
"""
import cmath
import numpy as np
from spherical_em_induction import solve_layered, k_of

# --- same uniform-sphere test config used throughout this session's
#     l=1,m=-1 investigation (matches demo_general_lm.py) ---
mu0 = 4 * np.pi * 1e-7
a = 6.371e6          # Earth radius, m
rho = 100.0          # Ohm.m
sigma = 1.0 / rho
T = 1000.0           # s
omega = 2 * np.pi / T
l = 1

c = 200.0 * a        # distant source shell -> approximates uniform ext. field
H0_target = 1.0
K0 = H0_target * (2 * l + 1) * c ** l / mu0   # -> A ~= H0_target = 1

# --- pythonSolver's own Rval, Rpval at r=a (matches demo_general_lm.py's
#     r_eval = a exactly) ---
sol = solve_layered(l, [a], [sigma], [mu0], mu0, c, K0, omega)
Rval = sol['R'](np.array([a]))[0]
Rpval = sol['Rp'](np.array([a]))[0]
A_check = sol['A']
print(f"solve_layered: A = {A_check}  (should be close to 1.0 -- confirms")
print(f"               the distant-shell normalization is working as expected)")
print(f"Rval(a)  (pythonSolver, e^+iwt) = {Rval}")
print(f"Rpval(a) (pythonSolver, e^+iwt) = {Rpval}")

# --- independent closed-form prediction, SAME formulas as
#     uniform_sphere_Tnr_predict.py, but with pythonSolver's OWN k
#     convention (k_of imported directly from spherical_em_induction.py --
#     not reimplemented -- to guarantee an apples-to-apples k) ---
def psi1(z):
    return cmath.sin(z) / z - cmath.cos(z)

def psi1p(z):
    return cmath.cos(z) / z - cmath.sin(z) / z ** 2 + cmath.sin(z)

k_py = k_of(omega, mu0, sigma)     # pythonSolver's own e^{+iwt} k = sqrt(-i*omega*mu*sigma)
x = k_py * a
ratio = psi1p(x) / psi1(x)
A1 = (1 + x * ratio) / 3
Tnr_pred = 1 / A1
Tnsp_pred = (3 - Tnr_pred) / a

print(f"\nk_py = {k_py}   x = k_py*a = {x}")
print(f"Predicted Rval(a)  (closed-form, pythonSolver convention) = {Tnr_pred}")
print(f"Predicted Rpval(a) (closed-form, pythonSolver convention) = {Tnsp_pred}")

# --- comparison ---
def relerr(got, expected):
    return abs(got - expected) / abs(expected)

print(f"\n--- Comparison (raw, before accounting for pythonSolver's R0^2 units convention) ---")
ratio_val_raw = Rval / Tnr_pred
ratio_pval_raw = Rpval / Tnsp_pred
print(f"Rval / predicted  = {ratio_val_raw}")
print(f"Rpval / predicted = {ratio_pval_raw}")
print(f"a^2 = {a**2:.6e}  (ratios above ~= a^2: Fortran's Tnr is non-dimensionalized")
print(f"with an explicit R0^2 reinserted in H%z's formula, Tnr*R0^2/Rr^2; pythonSolver's")
print(f"Hr = Rval*Slp/(mu*r^2*sinth) has no such R0^2 -- expected units difference, not a bug)")

print(f"\n--- Comparison (R0^2-corrected -- apples to apples) ---")
Rval_corr = Rval / a**2
Rpval_corr = Rpval / a**2
print(f"Rval/a^2:  pythonSolver={Rval_corr:.6e}  predicted={Tnr_pred:.6e}  rel.err={relerr(Rval_corr, Tnr_pred):.3e}")
print(f"Rpval/a^2: pythonSolver={Rpval_corr:.6e}  predicted={Tnsp_pred:.6e}  rel.err={relerr(Rpval_corr, Tnsp_pred):.3e}")
print(f"\nRval/a^2  / predicted  = {Rval_corr / Tnr_pred}   (expect ~1.0 if correct)")
print(f"Rpval/a^2 / predicted = {Rpval_corr / Tnsp_pred}   (expect ~1.0 if correct)")
