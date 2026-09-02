"""
Plot numerical (general transfer-matrix solver) vs analytical closed-form
solution for a single spherical-harmonic degree l=10 (P10) source shell.

Uses Earth-scale parameters (single homogeneous conducting sphere, so the
closed-form benchmark applies directly) and a diurnal excitation period.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from spherical_em_induction import (
    solve_layered, closed_form_single_sphere, k_of, rj, S_l
)

mu0 = 4 * np.pi * 1e-7
mu_earth = mu0
sigma = 1e-2                 # S/m, uniform mantle-like conductor
a = 6.371e6                  # Earth radius, m
c = a + 1.1e5                # ~110 km altitude source shell
T_period = 1000
omega = 2 * np.pi / T_period
K0 = 1.0
l = 1

# --- solve both ways ---
cf = closed_form_single_sphere(l, a, c, mu0, mu_earth, sigma, omega, K0)
gen = solve_layered(l, [a], [sigma], [mu_earth], mu0, c, K0, omega)

def R_closed_form(rv):
    if rv <= a:
        k1 = k_of(omega, mu_earth, sigma)
        return cf['C'] * rj(l, k1, rv)
    elif rv <= c:
        return cf['A'] * rv**(l+1) + cf['B'] * rv**(-l)
    else:
        return cf['D'] * rv**(-l)

r_grid = np.linspace(1e4, 1.5 * c, 2000)
R_cf = np.array([R_closed_form(rv) for rv in r_grid])
R_num = gen['R'](r_grid)

rel_err = np.abs(R_num - R_cf) / np.max(np.abs(R_cf))

fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                          gridspec_kw={'height_ratios': [3, 1.4]})

ax = axes[0]
ax.plot(r_grid/1e3, R_cf.real, 'C0-', lw=2.5, label='analytical, Re(R)')
ax.plot(r_grid/1e3, R_cf.imag, 'C1-', lw=2.5, label='analytical, Im(R)')
# sparse markers for the numerical solver so both curves are visible
idx = np.arange(0, len(r_grid), 40)
ax.plot(r_grid[idx]/1e3, R_num[idx].real, 'ko', ms=4, mfc='none', label='numerical, Re(R)')
ax.plot(r_grid[idx]/1e3, R_num[idx].imag, 'k^', ms=4, mfc='none', label='numerical, Im(R)')
ax.axvline(a/1e3, color='gray', ls=':', lw=1)
ax.axvline(c/1e3, color='gray', ls=':', lw=1)
ax.text(a/1e3, ax.get_ylim()[1]*0.9, ' Earth surface (a)', rotation=90, va='top', fontsize=8, color='gray')
ax.text(c/1e3, ax.get_ylim()[1]*0.9, ' source shell (c)', rotation=90, va='top', fontsize=8, color='gray')
ax.set_ylabel(r"$R_{\ell}(r)$")
ax.set_title(rf"$\ell$={l} (P$_{{{l}}}$) source: numerical vs. closed-form analytical solution")
ax.legend(fontsize=9, ncol=2)

ax2 = axes[1]
ax2.semilogy(r_grid/1e3, rel_err + 1e-20, 'r-')
ax2.axvline(a/1e3, color='gray', ls=':', lw=1)
ax2.axvline(c/1e3, color='gray', ls=':', lw=1)
ax2.set_ylabel("relative error")
ax2.set_xlabel("radius r (km)")
ax2.set_ylim(1e-17, 1e-13)

plt.tight_layout()
plt.savefig("l10_numerical_vs_analytical.png", dpi=140)

print(f"l = {l}")
print(f"max |R| over grid          : {np.max(np.abs(R_cf)):.6e}")
print(f"max relative error (num-cf): {np.max(rel_err):.3e}")
