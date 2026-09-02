# The toroidal source potential, the uniform inducing field, and the global1d↔ModEM normalization factor

*Author: Anna Kelbert with Claude, 2026-07. Companion PDF: `docs/source_normalization.pdf`.*

## Summary

global1d (`field1d.f90`/`field1d_s2.f90`) drives its solution with an **external toroidal potential
amplitude** `coeff(l,m)`. ModEM's regional 3-D solver drives its solution with a **boundary
electric field** of unit amplitude. Because a *potential* and a *field* are related by Faraday's law
`∇×E = iωμ₀H`, the two "unit sources" differ, and the raw electric fields the two codes produce for
their respective unit sources differ by a single multiplicative complex constant per period,

```
c = E_global1d / E_ModEM = i · ω · μ₀ · G ,
```

with `G` a real, geometric constant set by the **air-layer thickness** and the surface magnetic
response. This document explains the physics (the toroidal potential defines a quasi-uniform
*magnetic* field, not an electric one), derives `G` from the air-layer geometry, and defines a
normalized output preset that removes the factor and matches ModEM's normalization.

There is **no bug and no time-convention/sign error** in either code: `c` cancels identically in the
impedance `Z = E/H` (which is why global1d and ModEM agree on `Z`, apparent resistivity, and
transfer functions to <1%). The apparent "`−Re,+Im` vs `+Re,+Im` sign flip" seen when comparing raw
fields is simply the `≈90°` phase of the `iωμ₀` factor.

---

## 1. The toroidal potential and the fields it produces

Under the quasi-static (pre-Maxwell) approximation, a divergence-free electromagnetic field on a
sphere admits a Mie / poloidal–toroidal decomposition. In the insulating air between the ionospheric
source and the ground there are no currents, so the **poloidal** electric potential `P` vanishes
identically (Sun & Egbert 2012, eqs 21–22), and an external MT/ionospheric source drives only the
**toroidal** potential `T`. With `P≡0`, Sun & Egbert's eqs (5)–(6) reduce to

$$
\mathbf{E} = -i\omega\mu_0\,\frac{\hat{\mathbf r}}{r}\times\nabla_a T,
\qquad
\mathbf{H} = \frac{1}{r}\nabla_a\!\frac{\partial T}{\partial r} - \frac{\hat{\mathbf r}}{r^2}\nabla_a^2 T .
$$

Writing `T(r,θ,φ) = Σ coeff(l,m) T_l(r) Y_l^m(θ,φ)` and using
`r̂×θ̂ = φ̂`, `r̂×φ̂ = −θ̂`, `∇²_a Y_l^m = −l(l+1)Y_l^m`, the components are
(`field1d_s2.f90`, lines 42–47):

| | |
|---|---|
| `E_θ = +(iωμ₀/r)·T_l(r)·Y_p` | `Y_p = (1/sinθ)∂Y/∂φ` |
| `E_φ = −(iωμ₀/r)·T_l(r)·Y_t` | `Y_t = ∂Y/∂θ` |
| `E_r = 0` | ← **E is purely tangential (toroidal)** |
| `H_r = +(l(l+1)/r²)·T_l(r)·Y` | ← **H has a radial part (poloidal)** |
| `H_θ = +(1/r)·T_l'(r)·Y_t` | |
| `H_φ = +(1/r)·T_l'(r)·Y_p` | |

So the single scalar `T` generates a **toroidal electric field** (`E_r=0`, purely horizontal — the
name "toroidal potential" refers to *this*) and a **poloidal magnetic field** (`H_r≠0` — the
vertical/Z magnetic component measured at the ground). This is the standard deep-Earth induction
picture (Schuster–Lamb, Banks): the source is a poloidal magnetic field and E is the toroidal
induced response.

## 2. The (external) toroidal potential is a quasi-uniform *magnetic* field

In the insulating air, `T` satisfies Laplace's equation radially, so
`T_l(r) = α r^{l+1} + β r^{-l}` — an **external** part (`α`, the source) plus an **internal**
part (`β`, the induced response). global1d normalizes so that `coeff(l,m)=1` sets the coefficient of
`r^{l+1}` to exactly 1 (`field1d_s2.f90`, lines 437–458), i.e. `α = coeff`.

Take the external part alone for the P10 (l=1, m=0) case, with the fully-normalized
`Y_1^0 = \sqrt{3/4\pi}\,\cosθ`, so `T_1^{ext} = α r^2`, `T_1'^{ext} = 2αr`:

$$
H_r = \frac{l(l+1)}{r^2}T_1 Y = \frac{2}{r^2}(\alpha r^2)\sqrt{\tfrac{3}{4\pi}}\cosθ = 2\alpha\sqrt{\tfrac{3}{4\pi}}\cosθ,
\qquad
H_θ = \frac{1}{r}T_1' Y_t = -2\alpha\sqrt{\tfrac{3}{4\pi}}\sinθ .
$$

Since `cosθ\,\hat{\mathbf r} − \sinθ\,\hat{\boldsymbol θ} = \hat{\mathbf z}` exactly,

$$
\boxed{\;\mathbf H^{ext} = 2\alpha\sqrt{\tfrac{3}{4\pi}}\,\hat{\mathbf z}\;}
\qquad\text{— a spatially uniform magnetic field along the axis.}
$$

The electric field it induces is **not** uniform, and it carries the explicit `iωμ₀`:

$$
E_φ = -\frac{iωμ_0}{r}T_1 Y_t = +iωμ_0\,\alpha\, r\sqrt{\tfrac{3}{4\pi}}\sinθ \quad (\propto \sinθ,\ \text{and}\ \propto i\omega\mu_0).
$$

**So yes: the (external, l=1) toroidal potential defines a quasi-uniform inducing *magnetic* field
H, and the electric field is the `iωμ₀`-scaled induced response.** "Quasi-uniform" because the l=1
*external* part is exactly uniform; the total field also contains the induced internal part and (for
a realistic source) higher degrees. This is precisely the MT plane-wave source approximation.

## 3. global1d vs ModEM: `c = iωμ₀ G`

global1d's unit source (`coeff=1`, a potential ∝ uniform H) and ModEM's unit source (a boundary
E-field of amplitude 1) therefore differ by the potential-vs-field relationship. Comparing the two
codes' raw surface `E` for their respective unit sources (`USA.0.25×0.25` grid, fake pole `(0,90)`,
`layered.prm`), the ratio

```
c = E_global1d / E_ModEM
```

is (verified numerically):

- **constant over depth** (`|c|`, `arg c` steady from the top of the air column down through the
  upper mantle) — a linear (multiplicative) relationship, not a conjugation;
- **scaling as ω**: `|c|/ω` constant and `arg c → 90°` at short period ⇒ `c = i·ω·μ₀·G`;
- **mode-independent**: Mode 1 (zonal) and Mode 2 (non-zonal) give identical `arg c`.

The `i` is Faraday's law; the small residual (`arg c` drifting `90°→95°`, `|c|/ω` +9% over
`T=10→1000 s`) is the earth's induction response (§4) plus the flat-earth-vs-spherical difference in
ModEM's boundary condition.

## 4. Derivation of `G` from the air-layer geometry

### 4.1 ModEM's surface field: the air-thickness factor

ModEM computes its boundary condition with the Weerachai/REBOCC "quasi-1D" scheme: for each edge
column it solves a 1-D layered problem driven by a unit E-field at the **top of the domain**, a
height `d_air` above the surface (`Mod2d/WSfwd1Dmod.f90`; `d_air ≈ 1000 km` here). In the insulating
air the 1-D diffusion equation reduces to `d²E/dz² = 0`, so **E is linear in height and H is uniform**
through the air (consistent with §2). With `z` downward, `E(0)=1` at the top and the earth presenting
a surface impedance `Z_e = E_s/H_s`,

$$
E(z) = 1 + iωμ_0 H\,z,\quad H=\text{const},\qquad
E_{\text{ModEM}}(\text{surface}) = \frac{Z_e}{Z_e - iωμ_0 d_{air}} .
$$

Because `|Z_e| \ll ωμ_0 d_{air}` for a `1000\,km` air column (by a factor ~10 at `1000 s`, ~300 at
`10 s`), `E_{\text{ModEM}}(\text{surface}) \approx -Z_e/(iωμ_0 d_{air})` — **ModEM's surface E is set
by the air thickness.**

### 4.2 The ratio, and the exact flat-air form of `G`

global1d's surface field obeys the same physics, `E_φ = Z_e H_θ` at the surface (same `Z_e`, since
both codes solve the same layered earth). Taking the ratio, `Z_e` cancels:

$$
c = \frac{E_{\text{global1d}}}{E_{\text{ModEM}}}
  = H_θ\,(Z_e - iωμ_0 d_{air})
  = E_φ - iωμ_0 d_{air} H_θ
  = -iωμ_0\,H_θ\,(d_{air} + C),
$$

where `C \equiv Z_e/(iωμ_0)` is the **C-response** (a complex length, the induction depth,
`≈ 4, 28, 132 km` at `T = 10, 100, 1000 s` here). Hence

$$
\boxed{\;G = -H_θ(\text{surface})\,(d_{air} + C) \;\approx\; -H_θ\,d_{air}\;}
$$

The leading term is `|H_θ|·d_{air}` — the surface horizontal magnetic field times the air-layer
thickness. `C` is the period-dependent residual (it produces the observed `arg c` drift, growing
with period as the induction reaches deeper).

### 4.3 The geometric value of `H_θ(surface)`

From `T_1(r) = r^2 + β/r` (unit external coefficient) at the surface `r_0`,
`T_1'(r_0) = 2r_0 - β/r_0^2`, so

$$
H_θ(\text{surface}) = \frac{1}{r_0}T_1'(r_0)\,Y_t = -\Big(2 - \tfrac{β}{r_0^3}\Big)\sqrt{\tfrac{3}{4\pi}}\sinθ
= -(2 - Q_1)\sqrt{\tfrac{3}{4\pi}}\sinθ,
$$

where `Q_1 = β/r_0^{3}` is the l=1 **Q-response** (induced/external ratio). Its limits are `Q_1=0`
(insulator: `H_θ = -2\sqrt{3/4\pi}\sinθ`, purely the external uniform field) and `Q_1 = l/(l+1) = 1/2`
(perfect conductor). At MT/geomagnetic periods the mantle is a good conductor, `Q_1 \to 1/2`, giving

$$
H_θ(\text{surface}) \to -\tfrac{3}{2}\sqrt{\tfrac{3}{4\pi}}\sinθ ,
$$

and therefore, at the equator (`sinθ=1`, where the fake-pole comparison is evaluated),

$$
\boxed{\;G = \tfrac{3}{2}\sqrt{\dfrac{3}{4\pi}}\;d_{air} \;\approx\; 0.733\,d_{air}\;}
$$

The factor `3/2 = (2 - Q_1)` combines the uniform-external (2) and induced (`−Q_1`) parts; `√(3/4π)`
is the fully-normalized `Y_1^0` amplitude; `d_air` is the air-layer thickness (`grid%r(1)·1000 − r_0`
in metres, `≈ 1.000×10^6 m`).

### 4.4 Numerical verification

Measured `|H_θ|`, `|G|` (global1d, `coeff=1`, surface):

| `T` (s) | `|H_θ|` | `|G|/d_air` (measured) | `arg G` |
|---|---|---|---|
| 10 | 0.718 | 0.719 | +0.65° |
| 100 | 0.721 | 0.730 | +1.59° |
| 1000 | 0.720 | 0.783 | +4.73° |

Derived value `(3/2)\sqrt{3/4\pi} = 0.733`. Measured `|G|/d_air` straddles it (best at `100 s`); the
upward drift and the growing `arg G` are the `C`-response `(d_air+C)` correction — exactly the
`≈5°` residual noted in §3. `H_θ` is nearly period-independent (`0.718–0.721`), confirming it is a
geometric quantity in the deep-induction limit, `Q_1≈0.5`.

## 5. The normalized output preset

To express global1d's fields in ModEM's boundary-field normalization, multiply the output by `1/c`:

$$
(\mathbf E,\mathbf H)_{\text{normalized}} = \frac{1}{i\omega\mu_0 G}\,(\mathbf E,\mathbf H),
\qquad G = \tfrac{3}{2}\sqrt{\tfrac{3}{4\pi}}\,d_{air},
$$

with `d_air` taken from the grid and `ω` from the period. Applying the *same* complex factor to both
`E` and `H` (equivalently, rotating the source `coeff` by `1/i = -i` and scaling) **preserves the
impedance `Z=E/H`** — the physics is untouched; only the arbitrary source-normalization convention
is changed. This is implemented as the `EGBERTKELBERT2012_MODEM` output-convention preset
(`output_convention.f90`).

**Accuracy.** The factor `1/(iωμ₀G)` removes the exact `iωμ₀` and the air-thickness scaling. The
residual — the `C`-response `(d_air+C)` term, the finite-`Q` departure from `3/2`, and the
flat-earth-vs-spherical difference — leaves a mismatch of a few percent in amplitude and up to `~5°`
in phase at long period (→ 0 at short period). These residuals are genuine model/physics differences
(global1d being the more accurate, fully-spherical solution), not convention artifacts.

## 6. Validation against Cartesian ModEM

The spherical ModEM runs used to *measure* `c` above (`Mod3DMT_SP2_SPH`) require building with
`-DBUILD_SPHERICAL -DFORCE_SPHERICAL`, which bypasses a guard warning that the spherical version is
*experimental and "will not work with traditional MT source setup"* — and the internal `COMPUTE_BC`
boundary condition IS that traditional (WS/REBOCC per-column) setup. So the derivation was cross-
checked against **Cartesian ModEM**, which is trusted and accurate, on a ρ=100 Ω·m halfspace (same
1000 km / 12-layer air; `ModEM/debug-serial/CylinderModel/Halfspace.rho`).

*Air-column mechanism* (the heart of `G`): Cartesian ModEM's surface E equals the analytic
`Z_e/(Z_e − iωμ₀ d_air)` (with the `e^{−iωt}` halfspace impedance `Z_e = −iωμ₀/\sqrt{−iωμ₀σ}`) to
**0.2–0.6 %, ±0.1° at every period** from 4 to 4000 s.

*Full preset* — `r = E_{\text{global1d,norm}}(\text{equator}) / E_{\text{ModEM,cartesian}}`:

| `T` (s) | `|r|` | `arg r` |
|---|---|---|
| 3.98 | 1.006 | +0.2° |
| 63 | 1.019 | +0.9° |
| 1000 | 1.073 | +3.6° |
| 3981 | 1.145 | +6.4° |

The match is **essentially perfect in the flat-earth limit** (`|r|=1.006`, `+0.2°` at 4 s — where `G`
is exact and `Q_1→½`), and the residual grows **monotonically and cleanly** with period as Earth's
curvature enters. Cartesian ModEM is laterally uniform to `0.00 %`. This confirms the `iωμ₀G`
framework and its air-geometry derivation on trusted code, and shows the residual is the **genuine
flat-vs-spherical** physics — global1d being the more accurate spherical solution — not an artifact
of the experimental spherical BC. (Full detail: `testing/test_vs_modem_1D/cartesian_sanity_check.md`.)

### Generality

The derivation above is for the l=1 (P10 / MT-mode) source at the equator, in the deep-induction
limit. For a general multi-degree source or off-equator evaluation, `H_θ(surface)` — and hence `G` —
depends on `(l, θ)` via `H_θ = -[(l+1) - l\,Q_l]\,r_0^{\,l-1}\,∂Y_l^0/∂θ`, so a single scalar `G` is
exact only for the l=1 equatorial case that the MT modes reduce to. For the intended MT-interop use
case this is the relevant limit.
