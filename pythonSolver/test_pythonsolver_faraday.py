"""
Self-contained check: does spherical_em_induction.py's own field reconstruction
(fields_from_R_general) satisfy Faraday's law, curl(E) = i*omega*mu*H, at a
fixed radius (radial component only -- needs no radial derivatives, just
numerical theta/phi derivatives of Etheta, Ephi at fixed r)?

    (curl E)_r = (1/(r sin(theta))) * [ d(sin(theta)*Ephi)/dtheta - dEtheta/dphi ]

This isolates whether Hr and Etheta/Ephi -- BOTH built from the
now-independently-verified-correct Rval (see test_pythonsolver_Rval_Rpval.py)
-- are combined with the angular functions in a MUTUALLY CONSISTENT way
inside pythonSolver itself, with no reference to Fortran at all. If this
check fails, the bug is inside spherical_em_induction.py's own formulas
(likely a missing/extra factor of i in Etheta/Ephi or Hr specifically). If
it holds, that clears pythonSolver's field reconstruction too, and points
squarely at Fortran's H%z/E%x/E%y assembly in sourceField1d.
"""
import numpy as np
from spherical_em_induction import solve_layered, fields_from_R_general

mu0 = 4 * np.pi * 1e-7
a = 6.371e6
rho = 100.0
sigma = 1.0 / rho
T = 1000.0
omega = 2 * np.pi / T
l, m = 1, -1

c = 200.0 * a
K0 = 1.0 * (2 * l + 1) * c ** l / mu0   # A ~= 1 (unit external field), as before

sol = solve_layered(l, [a], [sigma], [mu0], mu0, c, K0, omega)
r_eval = a
Rval = sol['R'](np.array([r_eval]))[0]
Rpval = sol['Rp'](np.array([r_eval]))[0]

theta0, phi0 = 1.1, 0.7   # generic test point, away from poles

def fields(theta, phi):
    return fields_from_R_general(l, m, r_eval, theta, phi, Rval, Rpval, mu0, omega)

f0 = fields(theta0, phi0)
Hr0 = f0['Hr']

# numerical derivatives of Etheta, Ephi at fixed r=r_eval
h = 1e-6
def Etheta_of(theta, phi):
    return fields_from_R_general(l, m, r_eval, theta, phi, Rval, Rpval, mu0, omega)['Etheta']
def Ephi_of(theta, phi):
    return fields_from_R_general(l, m, r_eval, theta, phi, Rval, Rpval, mu0, omega)['Ephi']

dEphi_dtheta = (np.sin(theta0 + h) * Ephi_of(theta0 + h, phi0)
                - np.sin(theta0 - h) * Ephi_of(theta0 - h, phi0)) / (2 * h)
dEtheta_dphi = (Etheta_of(theta0, phi0 + h) - Etheta_of(theta0, phi0 - h)) / (2 * h)

curlE_r = (dEphi_dtheta - dEtheta_dphi) / (r_eval * np.sin(theta0))
rhs = 1j * omega * mu0 * Hr0

print(f"(curl E)_r (numerical, from Etheta/Ephi) = {curlE_r}")
print(f"i*omega*mu0*Hr (from Hr formula)         = {rhs}")
print(f"ratio (curl E)_r / (i*omega*mu0*Hr)       = {curlE_r / rhs}")
print(f"  (expect ~1.0 if Etheta/Ephi/Hr are mutually consistent per Faraday's law;")
print(f"   ~ -1 or ~ +-1j would indicate a real bug INSIDE spherical_em_induction.py)")

relerr = abs(curlE_r - rhs) / abs(rhs)
print(f"\nrelative error = {relerr:.3e}")
