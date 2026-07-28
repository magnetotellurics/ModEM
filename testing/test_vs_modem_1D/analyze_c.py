"""
Characterize the global1d-vs-ModEM E-field relationship c = E_global1d / E_ModEM.

Establishes (see zonal_sign_test_report.md) that c is a single multiplicative
complex constant per period, c = i*omega*mu0*G -- the potential-vs-field
source-normalization difference -- NOT a sign flip or conjugation.

This measures the RAW ratio, so it needs global1d built with the RAW
OUTPUT_CONVENTION='EGBERTKELBERT2012' (NOT the default _MODEM, which already
divides the output by i*omega*mu0*G and would make c trivially ~1). The
driver run_spherical_comparison.sh handles that recompile automatically for
its optional analyze_c step.

Inputs (both LOCAL to this directory):
  out.esoln                 -- spherical ModEM (USA_small), from run_spherical_comparison.sh
  <global1d.esoln>          -- global1d raw (EGBERTKELBERT2012), passed as argv[1]
Only external dependency: the compiled ModEM (Mod3DMT_SP2_SPH) + FWD1D -- paths
set at the top of run_spherical_comparison.sh.

Usage: python analyze_c.py <global1d.esoln> [iMode]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + os.sep + '..')
import numpy as np
from read_esoln import read_esoln


def profiles(dm, dg, imode):
    bm = [b for b in dm['blocks'] if b['iMode'] == imode]
    bg = [b for b in dg['blocks'] if b['iMode'] == imode]
    comp = 'Ey' if imode == 1 else 'Ex'   # ModEM iMode=1 Ey-pol, iMode=2 Ex-pol
    im, jm = dm['nx'] // 2, dm['ny'] // 2
    ig, jg = dg['nx'] // 2, dg['ny'] // 2
    for bmk in bm:
        om = bmk['omega']
        bgk = min(bg, key=lambda b: abs(b['omega'] - om))
        T = 2 * np.pi / om
        Em = bmk[comp][im, jm, :]
        Eg = bgk[comp][ig, jg, :]
        yield T, om, Em, Eg


if __name__ == '__main__':
    g1d = sys.argv[1]
    imode = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    dm = read_esoln('out.esoln')
    dg = read_esoln(g1d)
    for T, om, Em, Eg in profiles(dm, dg, imode):
        k = 12  # surface (k=13, 1-indexed)
        c = Eg[k] / Em[k]
        print(f"T={T:8.1f}s  omega={om:.4e}  surface c: |c|={abs(c):.4e} "
              f"arg={np.degrees(np.angle(c)):+7.2f}deg  |c|/omega={abs(c)/om:.4e}")
