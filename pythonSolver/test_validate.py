"""
Validate the general transfer-matrix layered-Earth solver against the
closed-form single-sphere benchmark (Document 1, Section 5), for both the
magnetostatic (omega=0) and diffusive (omega!=0) cases, across several l.

Also checks the omega->0 limit of the diffusive general solver approaches
the magnetostatic closed form (internal consistency check suggested in
the source document).
"""
import numpy as np
from spherical_em_induction import closed_form_single_sphere, solve_layered

mu0 = 4 * np.pi * 1e-7
mu1 = mu0 * 3.0     # give it a nontrivial permeability contrast too
sigma = 2.5
a, c = 0.85, 1.2
K0 = 1.7

print("=" * 78)
print("TEST 1: general N=1-shell solver vs closed-form benchmark")
print("=" * 78)
print(f"{'l':>3} {'omega':>10} {'max rel err in R(r) on grid':>30}")

r_test = np.linspace(0.05, 3.0, 400)

worst = 0.0
for l in [1, 2, 3, 5, 8]:
    for omega in [0.0, 1e-4, 3e-3, 5e-2]:
        cf = closed_form_single_sphere(l, a, c, mu0, mu1, sigma, omega, K0)
        gen = solve_layered(l, [a], [sigma], [mu1], mu0, c, K0, omega)

        # closed-form R(r) piecewise
        def R_cf(rv):
            if rv <= a:
                if omega == 0:
                    return cf['C'] * rv ** (l + 1)
                from spherical_em_induction import k_of, rj
                k1 = k_of(omega, mu1, sigma)
                return cf['C'] * rj(l, k1, rv)
            elif rv <= c:
                return cf['A'] * rv ** (l + 1) + cf['B'] * rv ** (-l)
            else:
                return cf['D'] * rv ** (-l)

        R_cf_vals = np.array([R_cf(rv) for rv in r_test])
        R_gen_vals = gen['R'](r_test)
        err = np.max(np.abs(R_cf_vals - R_gen_vals)) / np.max(np.abs(R_cf_vals))
        worst = max(worst, err)
        print(f"{l:>3} {omega:>10.1e} {err:>30.3e}")

print(f"\nWorst relative error across all (l, omega) tested: {worst:.3e}")
assert worst < 1e-9, "General solver does NOT match closed-form benchmark!"
print("PASSED: general layered solver reproduces the exact closed-form solution.")

print()
print("=" * 78)
print("TEST 2: diffusive solver -> magnetostatic limit as omega -> 0")
print("=" * 78)
for l in [1, 3]:
    cf0 = closed_form_single_sphere(l, a, c, mu0, mu1, sigma, 0.0, K0)
    for omega in [1e-2, 1e-4, 1e-6, 1e-8]:
        cf_w = closed_form_single_sphere(l, a, c, mu0, mu1, sigma, omega, K0)
        relA = abs(cf_w['A'] - cf0['A']) / abs(cf0['A'])
        relC = abs(cf_w['C'] - cf0['C']) / abs(cf0['C'])
        print(f"l={l}  omega={omega:.0e}   rel err A={relA:.2e}  rel err C={relC:.2e}")
    print()

print("=" * 78)
print("TEST 3: multi-shell (N=2) sanity check -- splitting one shell into")
print("two identical sub-shells at the same total sigma/mu must reproduce")
print("the single-shell result exactly (this exercises the new N>1 code")
print("path that has no closed-form counterpart to check against).")
print("=" * 78)
for l in [1, 4]:
    for omega in [0.0, 2e-3]:
        gen1 = solve_layered(l, [a], [sigma], [mu1], mu0, c, K0, omega)
        r_mid = a * 0.5
        gen2 = solve_layered(l, [r_mid, a], [sigma, sigma], [mu1, mu1], mu0, c, K0, omega)
        R1 = gen1['R'](r_test)
        R2 = gen2['R'](r_test)
        err = np.max(np.abs(R1 - R2)) / np.max(np.abs(R1))
        print(f"l={l} omega={omega:.1e}  max rel diff (1-shell vs 2-identical-shells) = {err:.3e}")
        assert err < 1e-9
print("PASSED: splitting a shell into sub-shells of identical properties")
print("reproduces the single-shell solution (internal consistency of the")
print("general N-layer code, independent of the closed-form benchmark).")
