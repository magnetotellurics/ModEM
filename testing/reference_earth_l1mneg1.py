"""
Reference values for testing/test_earth_l1mneg1.f90.

Uses spherical_em_induction.py's OWN solve_layered/fields_from_R_general
directly (both independently validated/fixed earlier in this session -- see
pythonSolver/test_pythonsolver_Rval_Rpval.py for solve_layered's Rval/Rpval,
and pythonSolver/test_pythonsolver_faraday.py for fields_from_R_general's
Faraday-law self-consistency, both passing to ~machine precision after the
2026-07-22 fixes). Unlike reference_unit_sphere_sunegbert2012.py (which used a
separately re-derived closed form), this script cross-checks
field1d_sunegbert2012.f90 against the SAME validated pythonSolver machinery used
throughout this project.

Physical setup matches test_earth_l1mneg1.f90 EXACTLY: r0=6.371e6 m, uniform
100 Ohm.m (sigma=0.01 S/m) from centre to surface, T=1000 s, pure (l=1,m=-1)
unit external multipole.

Normalization: rather than relying on the "distant source shell" (c=200*a)
trick to be an exact unit external field (it's only an approximation, though
a very good one -- see CLAUDE.md, "A ~= 1.0 to high precision"), this script
reads solve_layered's own returned external amplitude 'A' and divides by it.
pythonSolver's OWN convention is T(r,theta,phi)=R(r)*Y_l^m/sqrt(l(l+1)) (see
fields_from_R_general's docstring) -- an EXTRA factor of sqrt(l(l+1)) beyond
1/A is needed to reach field1d_sunegbert2012.f90's convention (T_l(r) multiplies
Y_l^m directly, no norm division). Found 2026-07-23 by direct comparison
against field1d_sunegbert2012.f90's actual output (l=1,m=-1): dividing by A alone
left a clean, uniform-across-all-5-components factor of sqrt(l(l+1))=sqrt(2)
in the magnitude.

Time convention and m: pythonSolver is native e^{+iwt}; field1d_sunegbert2012.f90 is
native e^{-iwt} (same convention as field1d.f90). For the SAME physical
source, F_fortran(l,+m) = conj(F_python(l,-m)) -- NOTE THE SIGN FLIP ON m --
because conjugating a field also conjugates its e^{i*m*phi} angular
dependence, turning an (l,+m) pattern into an (l,-m)-like one; conjugating
python's (l,+m) output directly (no m-flip) does NOT correspond to the same
source as fortran's (l,+m). Found 2026-07-23 the same way: conjugating
python's (l,+1) output instead of (l,-1) turned the leftover discrepancy
from (magnitude sqrt(2), varying phase per component) into (magnitude
sqrt(2), phase EXACTLY 0 or 180 degrees, split cleanly by field type).

H-vs-E sign: on top of the m-flip, H-components pick up an extra overall -1
that E-components do not, under this same conjugate+m-flip operation -- H is
a pseudovector (reverses under true time-reversal, since it is sourced by
currents) while E is a polar vector (does not) -- so this is expected, not a
bug. Confirmed empirically 2026-07-23: after the norm and m-flip fixes, all
5 components' ratios (field1d_sunegbert2012 / this-reference) were real, with H's
exactly -1 and E's exactly +1 (before including this sign); baking the sign
in here makes the expected ratio exactly +1 for every component.

Set SHOW_NAIVE_COMPARISON=True below to also print the WRONG (naive
conjugate-only, no norm/m-flip/H-sign correction) reference values, for a
side-by-side illustration of exactly what discrepancy pattern each missing
correction produces -- see testing/TESTING_MANUAL.md for the worked example.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'pythonSolver'))

import numpy as np
from spherical_em_induction import solve_layered, fields_from_R_general

SHOW_NAIVE_COMPARISON = False

mu0 = 4 * np.pi * 1e-7
a = 6.371e6
rho = 100.0
sigma = 1.0 / rho
T = 1000.0
omega = 2 * np.pi / T
l = 1
m_fortran = -1
m_python = -m_fortran   # sign-flip rule, see module docstring
norm = np.sqrt(l * (l + 1))

c = 200.0 * a
K0 = 1.0 * (2 * l + 1) * c ** l / mu0   # distant-shell trick, A ~= 1 (not relied on exactly)

sol = solve_layered(l, [a], [sigma], [mu0], mu0, c, K0, omega)
A = sol['A']
print(f"solve_layered: A (raw external amplitude for this K0) = {A}")

# ---- exact grid geometry from test_earth_l1mneg1.f90 ----
d2r = np.pi / 180.0
th_node5 = (60.0 + 4 * (150.0 - 60.0) / 8) * d2r      # grid%th(5)
dt5 = (150.0 - 60.0) / 8 * d2r
th_mid5 = th_node5 + dt5 / 2                            # grid%th(5)+dt(5)/2
ph_node1 = 0.0                                           # grid%ph(1)
dp1 = 360.0 / 8 * d2r
ph_mid1 = ph_node1 + dp1 / 2                             # grid%ph(1)+dp(1)/2

Rr2 = 6373000.0 - (6373000.0 - 6371500.0) / 2   # Rr(2), m
Rs2 = 6373000.0                                  # Rs(2), m

# (..., H_or_E sign: -1 for H components, +1 for E components -- see docstring)
points = {
    "H%z [Hr]":     (Rr2, th_node5, ph_node1, "Hr",     -1),
    "H%x [Hphi]":   (Rs2, th_node5, ph_mid1,  "Hphi",   -1),
    "H%y [Htheta]": (Rs2, th_mid5,  ph_node1, "Htheta", -1),
    "E%y [Etheta]": (Rr2, th_node5, ph_mid1,  "Etheta", +1),
    "E%x [Ephi]":   (Rr2, th_mid5,  ph_node1, "Ephi",   +1),
}

print(f"\n{'component':16s} {'(r,theta,phi)':45s} {'value @ native e^-iwt (fully corrected)':>50s}")
for name, (r, th, ph, which, sign) in points.items():
    Rval = sol['R'](np.array([r]))[0]
    Rpval = sol['Rp'](np.array([r]))[0]
    F = fields_from_R_general(l, m_python, r, th, ph, Rval, Rpval, mu0, omega)
    val_fortran_convention = sign * norm * np.conj(F[which] / A)
    print(f"{name:16s} r={r:12.3f} th={th:.10f} ph={ph:.10f}   "
          f"{val_fortran_convention.real:+.10e} {val_fortran_convention.imag:+.10e}j")

print()
print("Compare each value directly against test_earth_l1mneg1.f90's printed")
print("output for the same component (expect ratio FIELD1D_SE12/reference ~= 1+0i).")

if SHOW_NAIVE_COMPARISON:
    print()
    print("=" * 78)
    print("NAIVE (WRONG) comparison for illustration: conj(F_python(l,+m)/A),")
    print("no norm factor, no m-flip, no H-sign -- see TESTING_MANUAL.md")
    print("=" * 78)
    for name, (r, th, ph, which, sign) in points.items():
        Rval = sol['R'](np.array([r]))[0]
        Rpval = sol['Rp'](np.array([r]))[0]
        F_naive = fields_from_R_general(l, m_fortran, r, th, ph, Rval, Rpval, mu0, omega)
        val_naive = np.conj(F_naive[which] / A)
        print(f"{name:16s} {val_naive.real:+.10e} {val_naive.imag:+.10e}j")
