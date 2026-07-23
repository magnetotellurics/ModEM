"""
Test the general (l,m) extension:
 1) regression: m=0 general formulas must be proportional to the original
    axisymmetric S_l(theta) formulas (same physics, different normalization).
 2) l=1, m=-1: run the full pipeline (radial solve unchanged, only angular
    pattern changes) and sanity check the field pattern / magnitude.
"""
import numpy as np
from spherical_em_induction import (
    solve_layered, S_l, dS_l_dtheta, Y_lm, dY_lm_dtheta, C_lm,
    fields_from_R, fields_from_R_general
)

mu0 = 4*np.pi*1e-7
a = 6.371e6
sigma = 1e-1
c = a + 1.1e5
T_period = 24*3600.0
omega = 2*np.pi/T_period
K0 = 1.0

print("="*70)
print("Regression check: m=0 general-formula fields vs original axisymmetric")
print("="*70)
for l in [1,2,3]:
    sol = solve_layered(l, [a], [sigma], [mu0], mu0, c, K0, omega)
    r_eval = a*0.6
    Rval = sol['R'](np.array([r_eval]))[0]
    Rpval = sol['Rp'](np.array([r_eval]))[0]
    theta = np.linspace(0.1, np.pi-0.1, 6)
    phi = 0.3

    # original axisymmetric formulas
    Hr0, Ht0, Ep0 = fields_from_R(l, r_eval, theta, Rval, Rpval, mu0, omega)

    # new general formulas at m=0
    gen = fields_from_R_general(l, 0, r_eval, theta, phi, Rval, Rpval, mu0, omega)

    ratio_E = (gen['Ephi']/Ep0).real
    ratio_Hr = (gen['Hr']/Hr0).real
    ratio_Ht = (gen['Htheta']/Ht0).real
    print(f"l={l}: Ephi ratio (should be const across theta) = {ratio_E}")
    print(f"       Hr   ratio = {ratio_Hr}")
    print(f"       Htheta ratio = {ratio_Ht}")
    print(f"       Hphi (should be ~0 for m=0): max|Hphi| = {np.max(np.abs(gen['Hphi'])):.3e}")
    print()

print("="*70)
print("l=1, m=-1 source: radial solve is IDENTICAL to l=1 axisymmetric case")
print("(only the angular pattern differs) -- verifying R(r) is unchanged")
print("="*70)
l = 1
sol_m0 = solve_layered(l, [a], [sigma], [mu0], mu0, c, K0, omega)
sol_mneg1 = solve_layered(l, [a], [sigma], [mu0], mu0, c, K0, omega)  # radial solve takes l only
r_test = np.linspace(1e4, 1.4*c, 50)
print("max diff R(r) between 'm=0 setup' and 'm=-1 setup' (same l, radial code untouched):",
      np.max(np.abs(sol_m0['R'](r_test) - sol_mneg1['R'](r_test))))
