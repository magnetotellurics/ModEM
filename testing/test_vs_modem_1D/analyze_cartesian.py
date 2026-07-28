"""
Cartesian-ModEM sanity check of the G / EGBERTKELBERT2012_MODEM normalization
(see cartesian_sanity_check.md). Cartesian ModEM is trusted/accurate, unlike
the experimental spherical build (which needs -DFORCE_SPHERICAL).

Two checks on a rho=100 halfspace (both codes, ~1000 km air):
  1. air-column mechanism: does Cartesian ModEM's surface E equal the analytic
     Z_e/(Z_e - i*w*mu0*d_air) that underlies G?
  2. full preset: r = E_global1d_normalized(fake equator) / E_ModEM_cartesian.

Inputs (both generated LOCALLY in this directory by run_cartesian_sanity.sh):
  out_cart.esoln                          -- Cartesian ModEM, Halfspace.rho
  hs.s1.EGBERTKELBERT2012_MODEM.esoln     -- global1d halfspace (default preset)

Usage:
  bash run_cartesian_sanity.sh    # generates both .esoln (needs ModEM+FWD1D) then runs this
  python analyze_cartesian.py     # re-run the analysis on already-generated .esoln
The ONLY external dependency is the compiled Cartesian ModEM (Mod3DMT_SP2) and the
global1d FWD1D binary -- their paths are set at the top of run_cartesian_sanity.sh.
All model/template/control files (Halfspace.rho, Halfspace.fwd, halfspace_template.dat,
halfspace.prm, p10_halfspace.prm) live in this directory.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + os.sep + '..')
import numpy as np
from read_esoln import read_esoln

MU0 = 4e-7 * np.pi
RHO = 100.0


def main():
    dc = read_esoln('out_cart.esoln')            # generated locally by run_cartesian_sanity.sh
    dz = dc['dz']
    d_air = dz[:dc['nzAir']].sum()
    kc = dc['nzAir']
    sigma = 1.0 / RHO

    print(f'd_air = {d_air/1e3:.2f} km')
    print('\n=== Check 1: Cartesian ModEM surface E  vs  analytic Z_e/(Z_e - i*w*mu0*d_air) ===')
    print(' T(s)      |E/E_pred|   phase')
    for b in sorted(dc['blocks'], key=lambda x: x['omega']):
        if b['iMode'] != 2:
            continue
        om = b['omega']; T = 2 * np.pi / om
        Em = b['Ex'][dc['nx'] // 2, dc['ny'] // 2, kc]
        g = np.sqrt(-1j * om * MU0 * sigma)
        g = g if g.real > 0 else -g
        Ze = -1j * om * MU0 / g
        Ep = Ze / (Ze - 1j * om * MU0 * d_air)
        r = Em / Ep
        print(f'{T:8.2f}   {abs(r):.3f}      {np.degrees(np.angle(r)):+.1f}')

    hs = 'hs.s1.EGBERTKELBERT2012_MODEM.esoln'    # generated locally by run_cartesian_sanity.sh
    if not os.path.exists(hs):
        print('\n(global1d halfspace .esoln not present; skipping Check 2)')
        return
    dg = read_esoln(hs)

    def byom(d, im):
        return {round(b['omega'], 6): b for b in d['blocks'] if b['iMode'] == im}
    cc, gg = byom(dc, 1), byom(dg, 1)          # iMode=1 = Ey/east in both
    kg = 12
    print('\n=== Check 2: surface Ey (V/m) at the equator, and r = E_global1d_norm / E_ModEM_cart ===')
    print('   E_global1d_normalized -- global1d at its sin(theta)=1 fake equator, default preset')
    print('   E_ModEM_cartesian     -- Cartesian ModEM surface Ey (laterally-uniform mean)')
    print(f"  {'T(s)':>8}  {'E_global1d_normalized (Re,Im)':>32}  {'E_ModEM_cartesian (Re,Im)':>32}  {'|r|':>6} {'arg':>6}")
    for om in sorted(cc.keys(), reverse=True):
        omg = min(gg.keys(), key=lambda o: abs(o - om))
        if abs(omg - om) / om > 0.02:
            continue
        T = 2 * np.pi / om
        Ec, Eg = cc[om]['Ey'], gg[omg]['Ey']
        cv = np.array([Ec[i, j, kc] for i in range(2, dc['nx'] - 1)
                       for j in range(2, dc['ny'] - 2)])
        cmean = np.mean(cv)
        j = Eg.shape[1] // 2
        ge = Eg[np.argmax(np.abs(Eg[:, j, kg])), j, kg]
        r = ge / cmean
        print(f"  {T:8.2f}  {ge.real:+.5e}{ge.imag:+.5e}j  {cmean.real:+.5e}{cmean.imag:+.5e}j  "
              f"{abs(r):.3f} {np.degrees(np.angle(r)):+5.1f}")


if __name__ == '__main__':
    main()
