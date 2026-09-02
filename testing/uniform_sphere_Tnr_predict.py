"""
Independent prediction of sourcePotential's Tnr(surface, l=1) for a uniform
200 Ohm.m sphere at T=1000s, r0=6371km -- to compare directly against the
value printed by test_Tnr_uniform_sphere.f90.

This is NOT a re-derivation of the E-field/Faraday-law physics (that attempt
had an unresolved sign error). It is the algebraic trace of
what sourcePotential ITSELF computes (field1d.f90:472-731), for the special
case of a single uniform layer down to r=0 (lmax=1):

  sourcePotential does two normalizing divisions:
    (1) Tnr /= tn        where tn = psi_1(k*r0)      [rbsls0 type 1, at r=r0]
    (2) Tnr /= tni       where tni = A1  (the "incident" airprop decomposition
                                          coefficient, from boundary values
                                          T(r0)=1, T'(r0)=k*psi_1'(k*r0)/psi_1(k*r0))

  Because both rbsls0 calls in step (1) share the identical z0=k*r0 and l=1,
  the internal numerical-stability scale factor (exp(-|Im z0|), z0^-l, etc.)
  cancels EXACTLY in the ratio, leaving:

      Tnr_code(r, l=1) = psi_1(k*r) / [psi_1(k*r0) * A1]

  At r=r0 (surface): Tnr_code(r0,1) = 1/A1  exactly (psi_1(kr0) cancels top
  and bottom). This is the quantity predicted below.

  A1 = ( 1 + (k*r0)*psi_1'(k*r0)/psi_1(k*r0) ) / 3

This algebraic reduction is checked against two independent, unambiguous
physical limits (not just internal self-consistency):
  - sigma -> 0 (no induction): A1 -> 1 (rn0=1, rnp0->0 in this limit since
    tnp->0 as k->0, so A1=(1+0)/3... wait see note below), i.e. Tnr->1.
    (Tnr=1 at the surface in the sigma->0 limit is the correct trivial
    "total = incident" result, consistent with no induction happening.)
  - sigma -> infinity (perfect conductor): psi_1'(x)/psi_1(x) -> i (shown
    via the same asymptotic argument as the analytic Hr check), so
    A1 -> (1+i*k*r0)/3 -> grows without bound -> Tnr=1/A1 -> 0, matching
    total-field flux exclusion at the surface (no field can be "total" if
    it's all excluded -- consistent with the earlier Hr->0 check).

--- Tnsp addition (2026-07-22) ---
Motivated by: comparing global1d against pythonSolver for l=1,m=-1 shows
Hx,Hy (built from Tnsp) conjugate cleanly as expected, but Ex,Ey,Hz (all
built from Tnr) come out with Re/Im SWAPPED instead -- i.e. an extra factor
of i relative to simple conjugation. Since Tnr and Tnsp are the only two
things that structurally differ between those two groups, and it can be
shown that R(r) and R'(r) BOTH transform via simple conjugation under the
time-convention switch k->conj(k) (real-power-series argument applied
term-by-term to r*j_l(kr) and its derivative -- no inherent asymmetry from
that switch alone), any discrepancy has to be Fortran-specific, inside
sourcePotential's Tnr computation specifically (not Tnsp, not the E/H
wrapper formulas). This addition predicts Tnsp(surface,l=1) too, from the
SAME already-derived quantities (no new independent derivation, low risk):

  Tnsp(r) = d(Tnr)/dr = k*psi_1'(kr) / [psi_1(kr0)*A1]
  At r=r0: Tnsp(r0) = [k*psi_1'(kr0)/psi_1(kr0)] / A1

  Using 3*A1 = 1 + (k*r0)*psi_1'(kr0)/psi_1(kr0)  (already established above):
    k*psi_1'(kr0)/psi_1(kr0) = (3*A1 - 1)/r0
  So:
    Tnsp(r0) = (3*A1 - 1) / (r0*A1) = (3 - 1/A1)/r0 = (3 - Tnr(r0))/r0

Compare BOTH Tnr(surface,1) and Tnsp(surface,1) printed here against
test_Tnr_uniform_sphere.f90's output for the same parameters, to see
directly whether Tnr specifically (not Tnsp) is where the extra i enters.
"""
import cmath, math

mu0 = 1.256637e-6

def psi1(z):
    return cmath.sin(z)/z - cmath.cos(z)

def psi1p(z):
    return cmath.cos(z)/z - cmath.sin(z)/z**2 + cmath.sin(z)

def predict_Tnr_Tnsp_surface(rho, T, r0):
    sigma = 1.0/rho
    omega = 2*math.pi/T
    k = cmath.sqrt(1j*omega*mu0*sigma)
    x = k*r0
    ratio = psi1p(x)/psi1(x)          # psi_1'(x)/psi_1(x)
    A1 = (1 + x*ratio)/3
    Tnr = 1/A1
    Tnsp = (3 - Tnr)/r0
    return Tnr, Tnsp, x, A1

if __name__ == "__main__":
    # matches demo_general_lm.py's current test config (sigma=1e-2 S/m -> 100
    # Ohm.m, T=1000s, a=6.371e6m) so this compares directly against the l=1,
    # m=-1 swap investigation
    r0 = 6.371e6
    rho = 100.0
    T = 1000.0
    Tnr_pred, Tnsp_pred, x, A1 = predict_Tnr_Tnsp_surface(rho, T, r0)
    print(f"x = k*r0 = {x}")
    print(f"A1 = {A1}")
    print(f"Predicted Tnr(surface, l=1)  = {Tnr_pred}")
    print(f"  Re = {Tnr_pred.real:.12e}   Im = {Tnr_pred.imag:.12e}")
    print(f"Predicted Tnsp(surface, l=1) = {Tnsp_pred}")
    print(f"  Re = {Tnsp_pred.real:.12e}   Im = {Tnsp_pred.imag:.12e}")

    # sanity checks against the two unambiguous physical limits
    Tnr_dc, Tnsp_dc, _, A1_dc = predict_Tnr_Tnsp_surface(1e30, T, r0)   # rho->inf, sigma->0
    print(f"\nsigma->0 check: A1={A1_dc}, Tnr={Tnr_dc}  (expect Tnr->1)")

    # (rho=1e-12 overflows cmath.sin/cos at these frequencies -- the Fortran
    # code handles this via its exp(-|Im z0|) scaling trick; for this sanity
    # check we just use a milder but still much-more-conductive rho and show
    # the trend, which is enough given the DC limit above is already exact)
    for rho_pc in (10.0, 1.0):
        Tnr_pc, Tnsp_pc, x_pc, A1_pc = predict_Tnr_Tnsp_surface(rho_pc, T, r0)
        print(f"sigma->inf trend: rho={rho_pc:<6} |x|={abs(x_pc):9.3e}  |Tnr|={abs(Tnr_pc):.3e}  (should shrink toward 0)")

    print(f"\nCompare 'Predicted Tnr/Tnsp(surface, l=1)' above directly against the")
    print(f"values printed by test_Tnr_uniform_sphere.f90 for the same parameters")
    print(f"(rho={rho} Ohm.m, T={T}s, r0={r0:.3e}m).")
