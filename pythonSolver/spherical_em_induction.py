"""
Spherical-harmonic solver for the EM induction problem in a radially
layered, spherically symmetric Earth, driven by an external ionospheric-type
current source expanded in spherical harmonics.

Physics / conventions
----------------------
Time dependence e^{i w t}.  Axisymmetric (m=0), purely toroidal-source
("poloidal-field") problem: the source is an azimuthal surface current
  K_phi(theta) = sum_l K0_l * sin(theta) * P_l'(cos theta)
on a shell of radius c (the "ionosphere"), c > a = Earth radius.

For each degree l this reduces to a 1-D radial problem for the flux
function  psi(r,theta) = R_l(r) * S_l(theta),  S_l = sin(theta) P_l'(cos theta):

    R'' - l(l+1)/r^2 R + k^2 R = 0,      k^2 = -i*omega*mu*sigma

with general solution
    R(r) = alpha * r*j_l(k r) + beta * r*y_l(k r)      (k != 0, conducting)
    R(r) = alpha * r^(l+1)    + beta * r^(-l)          (k  = 0, vacuum / DC)

Two solvers are provided for each l:

  * closed_form_single_sphere(...)   -- the exact 4-unknown solution for a
        SINGLE homogeneous conducting sphere (0<r<a) + vacuum gap + source
        shell at c + vacuum beyond. This reproduces (and was checked
        against) the closed-form benchmark formulas.

  * solve_layered(...)               -- a general transfer-matrix solver for
        an arbitrary number of concentric homogeneous conducting shells,
        built directly from the same physics (continuity of R and of R'/mu
        at every interface, jump in R'/mu at the source shell r=c).
        Setting a single shell reproduces closed_form_single_sphere to
        machine precision (this is the validation test in test_validate.py).

Fields are recovered from R_l(r) via
    H_r(r,theta)     =  R_l(r)  * S_l'(theta) / (mu r^2 sin theta)
    H_theta(r,theta) = -R_l'(r) * S_l(theta)  / (mu r   sin theta)
    E_phi(r,theta)   = -i*omega * R_l(r) * S_l(theta)/sin(theta)     (conductor only, psi=r sin(theta) A_phi)
"""

import numpy as np
from scipy.special import spherical_jn, spherical_yn, lpmv


# ----------------------------------------------------------------------
# angular function S_l(theta) = sin(theta) P_l'(cos(theta)), and its
# theta-derivative, needed for H_r, H_theta reconstruction.
# We get P_l'(cos theta) from the associated Legendre function of order 1:
#   P_l^1(x) = -sqrt(1-x^2) P_l'(x)  =>  P_l'(x) = -P_l^1(x)/sin(theta)
# so S_l(theta) = sin(theta)*P_l'(cos theta) = -P_l^1(cos theta)
# ----------------------------------------------------------------------

def S_l(l, theta):
    """
    CORRECTED angular eigenfunction of Document 1's Eq (9).

    Document 1 Eq (7) defines Sl(theta) = sin(theta)*Pl'(cos theta) and claims
    (Eq 9) that this satisfies  sin(theta) d/dtheta[(1/sin theta) dSl/dtheta]
    = -l(l+1) Sl(theta).  This is checked here to be FALSE (nonzero residual,
    confirmed both symbolically and numerically for l=1..4). The function that
    actually satisfies Eq (9) with eigenvalue -l(l+1) is

        Theta_l(theta) = sin(theta) * [sin(theta) Pl'(cos theta)]
                        = -sin(theta) * P_l^1(cos theta)

    i.e. Document 1's Sl is missing one factor of sin(theta). This is what
    S_l() below returns. (The radial closed-form solutions R(r), A,B,C,D in
    this module never depended on Sl's specific functional form -- only on
    l(l+1) as the separation constant -- so they are unaffected by this fix;
    only the angular source shape / Hr,Htheta,Ephi reconstruction below use it.)
    """
    x = np.cos(theta)
    return -np.sin(theta) * lpmv(1, l, x)


def dS_l_dtheta(l, theta, h=1e-6):
    # robust central-difference derivative (S_l has no elementary closed form
    # that's simpler than this without extra recursion bookkeeping)
    return (S_l(l, theta + h) - S_l(l, theta - h)) / (2 * h)


# ----------------------------------------------------------------------
# Radial basis functions and derivatives
# ----------------------------------------------------------------------

def k_of(omega, mu, sigma):
    """k^2 = -i*omega*mu*sigma  (document's e^{+iwt} convention, Remark 1)."""
    return np.sqrt(-1j * omega * mu * sigma + 0j)


def rj(l, k, r):
    """r * j_l(k r), valid for k=0 too (returns r^(l+1)/(2l+1)!! * ... );
    here we just special-case k=0 to the exact power-law limit."""
    if k == 0:
        return r ** (l + 1)
    return r * spherical_jn(l, k * r)


def rj_prime(l, k, r, h=None):
    """ d/dr [ r j_l(k r) ] = j_l(kr) + k r j_l'(kr) """
    if k == 0:
        return (l + 1) * r ** l
    return spherical_jn(l, k * r) + k * r * spherical_jn(l, k * r, derivative=True)


def ry(l, k, r):
    if k == 0:
        return r ** (-l)
    return r * spherical_yn(l, k * r)


def ry_prime(l, k, r):
    if k == 0:
        return -l * r ** (-l - 1)
    return spherical_yn(l, k * r) + k * r * spherical_yn(l, k * r, derivative=True)


# ----------------------------------------------------------------------
# 1) Closed-form single-homogeneous-sphere benchmark
#    (this is Document 1, Section 5.5 / 5.6, reproduced & validated)
# ----------------------------------------------------------------------

def closed_form_single_sphere(l, a, c, mu0, mu1, sigma, omega, K0):
    """Returns dict with A,B,C,D coefficients for:
       R1(r) = C r j_l(k1 r),      0<r<a
       R2(r) = A r^(l+1) + B r^-l, a<r<c
       R3(r) = D r^-l,             r>c
    """
    if omega == 0:
        Delta = l * (mu0 + mu1) + mu0
        A = K0 * mu0 / (2 * l + 1) * c ** (-l)
        C = K0 * mu0 * mu1 / Delta * c ** (-l)
        B = K0 * mu0 * (l + 1) * (mu1 - mu0) / ((2 * l + 1) * Delta) * a ** (2 * l + 1) * c ** (-l)
        D = A * c ** (2 * l + 1) + B
        return dict(A=A, B=B, C=C, D=D)
    k1 = k_of(omega, mu1, sigma)
    Fa = rj(l, k1, a)
    Fap = rj_prime(l, k1, a)
    denom = l * mu1 * Fa + mu0 * a * Fap
    A = K0 * mu0 / (2 * l + 1) * c ** (-l)
    C = K0 * mu0 * mu1 / denom * a ** (l + 1) * c ** (-l)
    B = K0 * mu0 * ((l + 1) * mu1 * Fa - mu0 * a * Fap) / ((2 * l + 1) * denom) * a ** (2 * l + 1) * c ** (-l)
    D = A * c ** (2 * l + 1) + B
    return dict(A=A, B=B, C=C, D=D)


# ----------------------------------------------------------------------
# 2) General transfer-matrix layered-Earth solver (arbitrary # of shells)
#    "reuses the approach": spherical-harmonic separation + matching R and
#    R'/mu at every interface, exactly like Sun & Egbert Section 2, but
#    written as one global linear solve per l instead of a recursion.
# ----------------------------------------------------------------------

def solve_layered(l, radii, sigmas, mus, mu0, c, K0, omega):
    """
    radii  : increasing list [r1, r2, ..., rN=a] of shell OUTER radii
             (shell j occupies r_{j-1} < r < r_j, with r_0 = 0)
    sigmas : conductivity of each shell, length N
    mus    : permeability of each shell, length N
    mu0    : vacuum permeability
    c      : source-shell radius (c > radii[-1])
    K0     : source amplitude
    omega  : angular frequency (omega=0 -> magnetostatic)

    Returns: list of per-shell coefficients [(C1,), (Cj,Ej) for j=2..N],
             plus (A,B,D) for the two vacuum regions, and callables to
             evaluate R(r) and its derivative anywhere.
    """
    N = len(radii)
    assert len(sigmas) == N and len(mus) == N
    assert c > radii[-1]

    ks = [k_of(omega, mus[j], sigmas[j]) if (omega != 0 and sigmas[j] > 0) else 0.0
          for j in range(N)]

    # unknown vector: [C1, C2, E2, C3, E3, ..., CN, EN, A, B, D]
    n_shell_unknowns = 1 + 2 * (N - 1)
    n_unknowns = n_shell_unknowns + 3
    M = np.zeros((n_unknowns, n_unknowns), dtype=complex)
    rhs = np.zeros(n_unknowns, dtype=complex)

    def shell_R(j, r, coeffs):
        # coeffs: (Cj,) for j==0 (innermost), else (Cj,Ej)
        if j == 0:
            return coeffs[0] * rj(l, ks[0], r)
        return coeffs[0] * rj(l, ks[j], r) + coeffs[1] * ry(l, ks[j], r)

    def shell_Rp(j, r, coeffs):
        if j == 0:
            return coeffs[0] * rj_prime(l, ks[0], r)
        return coeffs[0] * rj_prime(l, ks[j], r) + coeffs[1] * ry_prime(l, ks[j], r)

    # index bookkeeping for unknown vector slots
    idx = []  # idx[j] = tuple of column indices for shell j's coeffs
    pos = 0
    idx.append((pos,)); pos += 1
    for j in range(1, N):
        idx.append((pos, pos + 1)); pos += 2
    iA, iB, iD = pos, pos + 1, pos + 2

    row = 0
    # interfaces between shells j-1 and j, at r = radii[j-1], for j=1..N-1
    for j in range(1, N):
        r_if = radii[j - 1]
        # continuity of R
        for col, sign in zip(idx[j - 1], (1, 1)[:len(idx[j - 1])]):
            pass
        # build row via basis evaluation at r_if for shell j-1 (outer side= inner shell's outer radius)
        # shell j-1 evaluated at its own outer boundary:
        if j - 1 == 0:
            M[row, idx[0][0]] = rj(l, ks[0], r_if)
        else:
            M[row, idx[j - 1][0]] = rj(l, ks[j - 1], r_if)
            M[row, idx[j - 1][1]] = ry(l, ks[j - 1], r_if)
        M[row, idx[j][0]] = -rj(l, ks[j], r_if)
        M[row, idx[j][1]] = -ry(l, ks[j], r_if)
        rhs[row] = 0
        row += 1

        # continuity of R'/mu
        if j - 1 == 0:
            M[row, idx[0][0]] = rj_prime(l, ks[0], r_if) / mus[0]
        else:
            M[row, idx[j - 1][0]] = rj_prime(l, ks[j - 1], r_if) / mus[j - 1]
            M[row, idx[j - 1][1]] = ry_prime(l, ks[j - 1], r_if) / mus[j - 1]
        M[row, idx[j][0]] = -rj_prime(l, ks[j], r_if) / mus[j]
        M[row, idx[j][1]] = -ry_prime(l, ks[j], r_if) / mus[j]
        rhs[row] = 0
        row += 1

    # interface at r = a = radii[-1] between last shell (N-1) and vacuum (A,B)
    a = radii[-1]
    if N - 1 == 0:
        M[row, idx[0][0]] = rj(l, ks[0], a)
    else:
        M[row, idx[N - 1][0]] = rj(l, ks[N - 1], a)
        M[row, idx[N - 1][1]] = ry(l, ks[N - 1], a)
    M[row, iA] = -a ** (l + 1)
    M[row, iB] = -a ** (-l)
    rhs[row] = 0
    row += 1

    if N - 1 == 0:
        M[row, idx[0][0]] = rj_prime(l, ks[0], a) / mus[0]
    else:
        M[row, idx[N - 1][0]] = rj_prime(l, ks[N - 1], a) / mus[N - 1]
        M[row, idx[N - 1][1]] = ry_prime(l, ks[N - 1], a) / mus[N - 1]
    M[row, iA] = -(l + 1) * a ** l / mu0
    M[row, iB] = -(-l) * a ** (-l - 1) / mu0
    rhs[row] = 0
    row += 1

    # interface at r = c: continuity of R, jump in R'/mu0 = K0
    M[row, iA] = c ** (l + 1)
    M[row, iB] = c ** (-l)
    M[row, iD] = -c ** (-l)
    rhs[row] = 0
    row += 1

    M[row, iA] = (l + 1) * c ** l / mu0
    M[row, iB] = (-l) * c ** (-l - 1) / mu0
    M[row, iD] = -(-l) * c ** (-l - 1) / mu0
    rhs[row] = K0
    row += 1

    sol = np.linalg.solve(M, rhs)

    shell_coeffs = [tuple(sol[i] for i in idx[j]) for j in range(N)]
    A, B, D = sol[iA], sol[iB], sol[iD]

    def R(r):
        r = np.atleast_1d(r).astype(float)
        out = np.zeros_like(r, dtype=complex)
        for i, rv in enumerate(r):
            if rv <= radii[-1]:
                j = next(jj for jj in range(N) if rv <= radii[jj])
                out[i] = shell_R(j, rv, shell_coeffs[j])
            elif rv <= c:
                out[i] = A * rv ** (l + 1) + B * rv ** (-l)
            else:
                out[i] = D * rv ** (-l)
        return out

    def Rp(r):
        r = np.atleast_1d(r).astype(float)
        out = np.zeros_like(r, dtype=complex)
        for i, rv in enumerate(r):
            if rv <= radii[-1]:
                j = next(jj for jj in range(N) if rv <= radii[jj])
                out[i] = shell_Rp(j, rv, shell_coeffs[j])
            elif rv <= c:
                out[i] = (l + 1) * A * rv ** l - l * B * rv ** (-l - 1)
            else:
                out[i] = -l * D * rv ** (-l - 1)
        return out

    def mu_at(r):
        for j in range(N):
            if r <= radii[j]:
                return mus[j]
        return mu0

    return dict(shell_coeffs=shell_coeffs, A=A, B=B, D=D, R=R, Rp=Rp, mu_at=mu_at)


# ----------------------------------------------------------------------
# Field reconstruction from a solved R(r), R'(r) for a single degree l
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# 3) General (l, m) toroidal-mode vector spherical harmonics.
#
# The radial physics (rj, ry, k_of, closed_form_single_sphere,
# solve_layered) does NOT depend on m at all -- l(l+1) is the eigenvalue
# of the angular Laplacian for every m, so the same R(r) solved above for
# a given l applies unchanged to any m. Only the ANGULAR pattern of the
# source current and of the resulting fields changes with m.
#
# The axisymmetric (m=0) source used above, K_phi(theta) = K0 sin(theta)
# P_l'(cos theta), is the m=0 special case of a toroidal-mode vector
# spherical harmonic source
#     K(theta,phi) = K0 * [ Ctheta_lm(theta,phi) theta_hat
#                          + Cphi_lm(theta,phi)   phi_hat ]
# with (Morse & Feshbach convention, matching Sun & Egbert Eq. 26)
#     C_l^m = -(1/sqrt(l(l+1))) r_hat x grad_a Y_l^m
# which works out to
#     Ctheta_lm =  (i*m)/(sin(theta)*sqrt(l(l+1))) * Y_l^m(theta,phi)
#     Cphi_lm   = -1/sqrt(l(l+1)) * dY_l^m/dtheta
# For m=0 this reduces to Ctheta=0, Cphi_lm = -Y_l^0'(theta)/sqrt(l(l+1))
# which is proportional to S_l(theta) used above (regression-tested below).
# ----------------------------------------------------------------------

from scipy.special import factorial as _factorial


def Y_lm(l, m, theta, phi):
    """Fully normalized complex spherical harmonic (orthonormal on unit sphere)."""
    norm = np.sqrt((2 * l + 1) / (4 * np.pi) * _factorial(l - m) / _factorial(l + m))
    P = lpmv(m, l, np.cos(theta))
    return norm * P * np.exp(1j * m * phi)


def dY_lm_dtheta(l, m, theta, phi, h=1e-6):
    return (Y_lm(l, m, theta + h, phi) - Y_lm(l, m, theta - h, phi)) / (2 * h)


def C_lm(l, m, theta, phi):
    """Toroidal-mode vector spherical harmonic components (Ctheta, Cphi)."""
    norm = np.sqrt(l * (l + 1))
    Y = Y_lm(l, m, theta, phi)
    dYdt = dY_lm_dtheta(l, m, theta, phi)
    Ctheta = (1j * m) / (np.sin(theta) * norm) * Y
    Cphi = -dYdt / norm
    return Ctheta, Cphi


def fields_from_R_general(l, m, r, theta, phi, Rval, Rpval, mu, omega):
    """
    General (l,m) field reconstruction for a pure-toroidal-potential
    (T-only, i.e. P=0) solution T(r,theta,phi) = R(r) Y_l^m(theta,phi)/sqrt(l(l+1)):

        H = (1/r) grad_a(dT/dr) - (r_hat/r^2) laplacian_a(T)
        E = -i*omega*mu * r_hat x grad_a(T)      (source-free / conductor interior)

    Reduces to fields_from_R() above when m=0 (regression-tested).
    """
    Ctheta, Cphi = C_lm(l, m, theta, phi)
    Y = Y_lm(l, m, theta, phi)
    dYdt = dY_lm_dtheta(l, m, theta, phi)
    norm = np.sqrt(l * (l + 1))

    # E = -i*omega*mu * r_hat x grad_a(T);  r_hat x grad_a(Y) = -sqrt(l(l+1)) * (Ctheta,Cphi)
    Etheta = 1j * omega * mu * Rval * norm * Ctheta / norm  # = i*omega*mu*Rval*Ctheta
    Ephi = 1j * omega * mu * Rval * Cphi

    # H components from T = R(r) Y/sqrt(l(l+1))
    Hr = -l * (l + 1) * Rval * Y / (norm * r ** 2)
    Htheta = Rpval / (r * norm) * dYdt
    Hphi = Rpval / (r * norm) * (1j * m / np.sin(theta)) * Y

    return dict(Etheta=Etheta, Ephi=Ephi, Hr=Hr, Htheta=Htheta, Hphi=Hphi)


def fields_from_R(l, r, theta, Rval, Rpval, mu, omega):
    Sl = S_l(l, theta)
    Slp = dS_l_dtheta(l, theta)
    sin_t = np.sin(theta)
    Hr = Rval * Slp / (mu * r ** 2 * sin_t)
    Htheta = -Rpval * Sl / (mu * r * sin_t)
    Ephi = -1j * omega * Rval * Sl / (r * sin_t) if omega != 0 else np.zeros_like(Rval)
    return Hr, Htheta, Ephi
