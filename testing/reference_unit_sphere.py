"""
Reference values for test_unit_sphere.f90.

Single homogeneous conducting sphere, radius r0, conductivity sigma, in
vacuum, driven by a UNIT external (l,m) multipole (Theta_ext(r)=r^(l+1),
amplitude 1). Computes the field components using Sun & Egbert (2012)
eq.(5)-(6), independently derived and numerically validated (see
validate_conjugation.py) against spherical_em_induction.py.

Physical parameters match test_unit_sphere.f90 EXACTLY:
  r0 = 1000 m, sigma = 1 S/m, omega = 1 rad/s, mu0 = 1.256637e-6 (field1d.f90's
  literal), pi = 3.14159265357898 (field1d.f90's literal), l=2, m=1.

Reports both time conventions (s=+1 native to field1d.f90's e^{-iwt}; s=-1 is
the e^{+iwt} convention used in spherical_em_induction.py) at the 5 exact
(r,theta,phi) points field1d.f90 will evaluate for H%z, H%x, H%y, E%y, E%x
given the test's hand-built grid (i=1, j=5, k=2 throughout).

IMPORTANT: field1d.f90 carries its own internal REAL normalization constants
per component (R0^2/r^n factors, and an extra 1/(l(l+1)) on the tangential
H and E components that H_r does not have) that this script does not
reproduce -- so absolute magnitudes will NOT match field1d.f90's raw output.
What MUST match, if field1d.f90 is bug-free, is the PHASE (complex argument)
of the ratio  field1d_component / conj(reference_component_at_s=-1)  --
that ratio must come out REAL (phase 0 or 180 deg) for all five components.
A phase of +/-90 degrees on H%z, E%y, E%x specifically (and not on H%x, H%y)
is the exact signature of the bug reported earlier in this conversation.
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
    Fa, Fap = rj(l, k1, a), rj_prime(l, k1, a)
    M = np.array([[-a ** (-l), Fa],
                  [l * a ** (-l - 1), -Fap]], dtype=complex)
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
    """s=+1: e^{-iwt} (field1d.f90 native). s=-1: e^{+iwt}."""
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
    "Hr (H%z)":     (Rr2, th_node5, ph_node1),
    "Hphi (H%x)":   (Rs2, th_node5, ph_mid1),
    "Htheta (H%y)": (Rs2, th_mid5,  ph_node1),
    "Etheta (E%y)": (Rr2, th_node5, ph_mid1),
    "Ephi (E%x)":   (Rr2, th_mid5,  ph_node1),
}
which_field = {"Hr (H%z)": "Hr", "Hphi (H%x)": "Hphi", "Htheta (H%y)": "Htheta",
               "Etheta (E%y)": "Etheta", "Ephi (E%x)": "Ephi"}

k_ft = np.sqrt(1j * omega * mu0 * sigma)     # field1d.f90's k (e^{-iwt})
k_py = np.conj(k_ft)                          # e^{+iwt} convention

B_ft, C_ft = solve_unit_external(l, k_ft, r0)
B_py, C_py = solve_unit_external(l, k_py, r0)

print(f"k (field1d convention) = {k_ft:.8e},  |k*r0| = {abs(k_ft)*r0:.6f}")
print()
print(f"{'component':16s} {'(r,theta,phi)':45s} {'value @ s=-1 (+iwt)':>45s}")
for name, (r, th, ph) in points.items():
    Th_py, Thp_py = Theta_ext_region(r, l, B_py)
    F_py = fields(Th_py, Thp_py, r, th, ph, l, m, s=-1)
    val = F_py[which_field[name]]
    print(f"{name:16s} r={r:8.2f} th={th:.10f} ph={ph:.10f}   {val.real:+.10e} {val.imag:+.10e}j")
