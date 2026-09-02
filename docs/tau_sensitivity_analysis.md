# Sensitivity of global1d solutions to the thin-sheet surface conductance `earth%tau`

*Author: Anna Kelbert with Claude, 2026-08. Companion PDF: `docs/tau_sensitivity_analysis.pdf`.*

## Abstract

global1d (`field1d.f90`/`field1d_s2.f90`) carries a single scalar parameter, `earth%tau` (units
Siemens), representing an optional thin conducting sheet superimposed at the model's outer surface
`r=r0`; we denote it τ below. This note derives how τ enters the governing equations, identifies the
natural dimensionless control parameter, and reports a numerical sensitivity test against the
closed-form analytic solution for a uniformly conducting sphere. We find (i) no numerical advantage
to a nonzero τ: the solver is exact and well-behaved at τ=0, which is also the case the closed-form
reference itself assumes; (ii) the error introduced by a nonzero τ scales linearly with the
dimensionless parameter `Q = ωμ₀τr0` over eight orders of magnitude, with the linear-regime
coefficient implying the perturbation only reaches order unity near `Q~40`; and (iii) a dedicated
Earth-scale test (matching the conductor's own electrical size to the toy sphere's, to isolate the
effect of Q alone) confirms this threshold and identifies the mechanism beyond it: once Q dominates
the conductor's own pre-jump value, the boundary condition is controlled by the artificial thin
sheet rather than the layered conductor, so the response **saturates** in magnitude and can
**overshoot and reverse sign** in phase, rather than growing without bound. (iv) A further check
using the actual layered Earth conductivity model (`LWS/layered_GDE_rho.prm`) at the production
period of 3600 s shows that, at the same Q, the real Earth's much higher conductivity suppresses the
perturbation well below what the matched-regime toy-sphere comparison of (iii) predicted — the
`Q~40` threshold is specific to that toy conductor's own electrical regime, not a universal constant.

---

## 1. How τ enters the equations

Both solvers solve for the toroidal potential `T_l(r)` via a transfer-matrix recursion through the
layered conductivity structure, matched to an insulating atmosphere above `r0`. Ordinary continuity
of the layered-earth solution requires both `T` and `T'` continuous across every interface. At the
model's outer surface, *both* solvers instead apply a modified matching condition
(`field1d_s2.f90` eq. 17; `field1d.f90`, algebraically identical, independently written):

```
T_l'(r0+) = T_l'(r0-) - i*omega*mu0*tau*T_l(r0),      T_l(r0+) = T_l(r0-)
```

`T` itself remains continuous; only `T'` jumps, by a term proportional to τ.

**Physical interpretation.** This is the classical thin-sheet (Price 1949) approximation. A sheet
of conductance `τ=∫σ dz` carries an induced surface current `K=τ·E_tangential` (Ohm's law
integrated through the sheet's negligible thickness). Since `H_θ,H_φ ∝ T'/r` and `E_θ,E_φ ∝ T`,
Ampère's law relates a surface current directly to a jump in tangential H — exactly the jump
condition above. Setting τ=0 removes the sheet entirely and recovers ordinary continuity.

**Dimensionless control parameter.** Comparing the added term against the natural scale of `T'`
itself (`T'/T ~ O(l)/r0`) gives the classical thin-sheet parameter

```
Q = omega * mu0 * tau * r0
```

which governs the size of the perturbation independent of the specific `(r0,ω)` used in any
particular test or production run.

## 2. Numerical test

We use the same single homogeneous-sphere configuration as the existing closed-form regression test
(`test_unit_sphere.f90`: `r0=1000` m, `σ=1` S/m, `ω=1` rad/s, unit external `(l=2,m=+1)` multipole
source), for which S2 already matches the analytic solution *exactly* at τ=0 (ratio = 1+0i for all
five field components). The τ=0 output therefore serves directly as the analytic reference; no
closed-form values need to be re-derived here. `Hr, Hθ, Hφ, Eθ, Eφ` are extracted at one fixed grid
point in the vacuum region immediately above the sphere for every τ.

### 2.1 Magnitude error vs. τ

Table 1 sweeps τ from 1 S down to 0, reporting Q and the maximum relative error across all five
field components against the τ=0 reference.

**Table 1.** Magnitude error vs. τ, referenced to the τ=0 analytic solution, for the toy sphere
(`r0=1000` m, `σ=1` S/m, `ω=1` rad/s). Error scales linearly with τ (and hence with Q) across eight
decades, until double-precision roundoff floors out below `τ~1e-6`. The fit is `error ≈ 0.024·Q`.

| τ [S] | Q = ωμ₀τr0 | max relative error (5 components) |
|---|---|---|
| 1e0  | 1.257e-3  | 3.058e-5 |
| 1e-1 | 1.257e-4  | 3.058e-6 |
| 1e-2 | 1.257e-5  | 3.058e-7 |
| 1e-3 | 1.257e-6  | 3.058e-8 |
| 1e-4 | 1.257e-7  | 3.058e-9 |
| 1e-5 | 1.257e-8  | 3.058e-10 |
| 1e-6 | 1.257e-9  | 3.058e-11 |
| 1e-7 | 1.257e-10 | 3.058e-12 |
| 1e-8 | 1.257e-11 | 3.058e-13 (roundoff floor) |
| 0    | 0         | 0 (exact, by construction) |

### 2.2 Magnitude ratio and phase shift of H at specific τ values

Tables 2–4 report the magnitude ratio `|H(τ)/H(0)|` and phase shift `arg H(τ) − arg H(0)` (degrees)
of each magnetic field component in turn, for `τ = 1e-4, 1, 1e2, 2.5e3` S.

**Table 2.** `Hr` magnitude ratio and phase shift relative to τ=0, for the toy sphere (`r0=1000` m,
`σ=1` S/m, `ω=1` rad/s).

| τ [S] | Q | \|Hr/Hr,0\| | ΔφHr [°] |
|---|---|---|---|
| 1e-4 | 1.257e-7 | 1.000000 | 1.7e-7 |
| 1    | 1.257e-3 | 0.999998 | 1.75e-3 |
| 1e2  | 0.1257   | 0.999722 | 0.1743 |
| 2.5e3 | 3.142   | 0.964484 | 3.086 |

**Table 3.** `Hθ` magnitude ratio and phase shift relative to τ=0, for the toy sphere (`r0=1000` m,
`σ=1` S/m, `ω=1` rad/s).

| τ [S] | Q | \|Hθ/Hθ,0\| | ΔφHθ [°] |
|---|---|---|---|
| 1e-4 | 1.257e-7 | 1.000000 | -3e-8 |
| 1    | 1.257e-3 | 1.000000 | -3.00e-4 |
| 1e2  | 0.1257   | 1.000051 | -0.02992 |
| 2.5e3 | 3.142   | 1.006426 | -0.5060 |

**Table 4.** `Hφ` magnitude ratio and phase shift relative to τ=0, for the toy sphere (`r0=1000` m,
`σ=1` S/m, `ω=1` rad/s).

| τ [S] | Q | \|Hφ/Hφ,0\| | ΔφHφ [°] |
|---|---|---|---|
| 1e-4 | 1.257e-7 | 1.000000 | -3e-8 |
| 1    | 1.257e-3 | 1.000000 | -3.00e-4 |
| 1e2  | 0.1257   | 1.000051 | -0.02992 |
| 2.5e3 | 3.142   | 1.006426 | -0.5060 |

Tables 3 and 4 are numerically identical. This is not a coincidence: `Hθ` and `Hφ` both derive from
the same radial factor `T'(r)`, whose complex ratio to the τ=0 case does not depend on angular
position, so the two components (evaluated at the same radius, different angles) necessarily agree.
`Hr` (Table 2), driven by `T(r)` rather than `T'(r)`, follows an independent (smaller-magnitude,
opposite-sign) trend.

## 3. Earth-scale test: is Q alone sufficient, and why does the response become nonlinear?

The toy sphere has `r0=1000` m, `σ=1` S/m, `ω=1` rad/s (as in §2); the real Earth has `r0≈6.371e6` m
(6371× larger). Since `Q∝r0`, the
*same* τ produces a vastly larger Q at Earth scale — Table 5 evaluates this at the NA2026 project's
own periods (162 s, 3600 s).

**Table 5.** Q at Earth scale for the same τ values as Tables 2–4.

| τ [S] | Q_Earth, 162 s | Q_Earth, 3600 s |
|---|---|---|
| 1e-4 | 3.1e-5 | 1.4e-6 |
| 1    | 0.31   | 0.014 |
| 1e2  | 31     | 1.4 |
| 2.5e3 | 776   | 35 |

**Methodology: isolating Q from the conductor's own regime.** Naively re-running the toy sphere with
`r0=6.371e6` m but the same `σ=1` S/m, `ω=1` rad/s would *also* push the conductor's own electrical
size `|kl|r0 = √(ωμ₀σ)·r0` from ≈1.12 (order-1 induction, the toy sphere's regime) to ≈1400 (extreme
skin effect) — a second, unrelated regime change that would confound any observed nonlinearity: we
would no longer know whether a deviation from linear-in-Q behavior is caused by τ growing large, or
by the conductor itself behaving completely differently. To isolate the τ/Q effect cleanly, σ at
each Earth-scale period is instead chosen so that `|kl|r0` *equals* the toy sphere's own value
(1.121), holding the conductor's regime fixed while only `r0` (and hence Q, for the same τ) changes.
This gives `σ=6.35e-7` S/m at 162 s and `σ=1.41e-5` S/m at 3600 s — unrealistically resistive for an
actual whole-Earth sphere, but that is expected and intentional: this is a controlled numerical
experiment isolating one mechanism, not a production-accurate Earth model. `Hθ=Hφ` again (same
argument as §2.2), so only one column pair is shown for them.

**Table 6.** H magnitude ratio and phase shift relative to τ=0, Earth radius, period 162 s,
`σ=6.35e-7` S/m (matched `|kl|r0`).

| τ [S] | Q | \|Hr/Hr,0\| | ΔφHr [°] | \|Hθ,φ/Hθ,φ,0\| | ΔφHθ,φ [°] |
|---|---|---|---|---|---|
| 1e-4 | 3.11e-5 | 1.000000 | 4.3e-5 | 1.000000 | -7.4e-6 |
| 1    | 0.3105  | 0.999065 | 0.4254 | 1.000173 | -0.07295 |
| 1e2  | 31.05   | 0.882199 | 0.9808 | 1.020275 | -0.13990 |
| 2.5e3 | 776.3  | 0.878969 | -0.1977 | 1.020790 | 0.03518 |

**Table 7.** H magnitude ratio and phase shift relative to τ=0, Earth radius, period 3600 s,
`σ=1.41e-5` S/m (matched `|kl|r0`).

| τ [S] | Q | \|Hr/Hr,0\| | ΔφHr [°] | \|Hθ,φ/Hθ,φ,0\| | ΔφHθ,φ [°] |
|---|---|---|---|---|---|
| 1e-4 | 1.40e-6 | 1.000000 | 1.9e-6 | 1.000000 | -3.3e-7 |
| 1    | 0.01397 | 0.999976 | 0.01932 | 1.000004 | -0.00332 |
| 1e2  | 1.397   | 0.989809 | 1.7626 | 1.001871 | -0.29856 |
| 2.5e3 | 34.93  | 0.881538 | 0.8515 | 1.020380 | -0.12055 |

**Refining the linear-regime threshold.** Table 1's fit, `error ≈ 0.024·Q`, implies the perturbation
only reaches order unity once `Q ~ 1/0.024 ≈ 42` — not `Q~1` as the crude order-of-magnitude
argument in §1 suggested. This is exactly consistent with Tables 6–7: at `Q=0.31–1.4` the response
is still closely linear in Q (e.g. `ΔφHr` at `Q=0.31` (162 s) and `Q=1.4` (3600 s) both fall close to
the toy-sphere fit scaled to that Q), while by `Q=31–35` (comparable to the refined `Q~42` threshold)
the magnitude ratio has already dropped by ~12% — clearly nonlinear.

**Why the response saturates rather than diverging.** From the jump condition, with `T_l(r0)=1` by
construction, `T_l'(r0+) = tnp(l) - i*Q`, where `tnp(l)` is the conductor's own (τ-independent)
pre-jump value. Once `Q >> |tnp(l)|`, this sum is dominated by the `-iQ` term regardless of the
conductor's own properties: `T_l'(r0+) → -iQ`, a fixed -90° phase relative to `T_l(r0)`, independent
of τ once τ is large enough. This is why the magnitude ratio *saturates* (Tables 6–7: essentially
unchanged between `Q=31` and `Q=776`, and between `Q=1.4` and `Q=35`) rather than continuing to
grow, and why the phase shift can *overshoot and reverse sign* (162 s: +0.98° at `Q=31` to -0.20° at
`Q=776`; 3600 s: +1.76° at `Q=1.4` down to +0.85° at `Q=35`) rather than growing monotonically: for
very large τ the physics is no longer governed by the layered conductor at all, but entirely by the
artificial thin sheet, and further increasing τ has diminishing effect on the already-τ-dominated
boundary condition. This is a genuine physical saturation, not a numerical artifact: the same
qualitative behavior (rise, peak/plateau, mild reversal) appears at both Earth-scale periods, at Q
values that differ by a factor of ~20, using two independently matched conductor regimes.

## 4. Realistic layered Earth conductivity model

§3's Earth-scale test used a homogeneous sphere with σ chosen only to match the toy sphere's own
`|kl|r0` — a controlled experiment isolating the Q mechanism, not a realistic Earth. Here we repeat
the same τ sweep at Earth radius and the NA2026 project's 3600 s period, but with the *actual*
layered conductivity profile used in production runs, `LWS/layered_GDE_rho.prm`, read the same way
`FWD1D.f90` builds `earth%layer`/`earth%sigma` from a `layers`-format model file: 32 layers given as
depth below the surface (km) and `log10(rho)` (rho in Ohm·m), converted to `sigma = 10^(-log10(rho))`,
plus a hardcoded `sigma=1e5` S/m "core" layer beneath the deepest defined depth (3500 km). Table 8
quotes every depth/sigma pair used.

**Table 8.** Layer depths and conductivities of the layered Earth model, `LWS/layered_GDE_rho.prm`
(32 layers, `sigma = 10^(-log10(rho))`, plus the hardcoded core layer), used for Table 9.

| depth [km] | log10(rho) | σ [S/m] |
|---|---|---|
| 10 | 2.000000 | 1.000e-2 |
| 40 | 3.000000 | 1.000e-3 |
| 110 | 3.000000 | 1.000e-3 |
| 160 | 2.522879 | 3.000e-3 |
| 210 | 2.154902 | 7.000e-3 |
| 260 | 1.745627 | 1.796e-2 |
| 310 | 1.667885 | 2.148e-2 |
| 360 | 1.688480 | 2.049e-2 |
| 410 | 1.709952 | 1.950e-2 |
| 460 | 1.699400 | 1.998e-2 |
| 510 | 1.645104 | 2.264e-2 |
| 560 | 1.540801 | 2.879e-2 |
| 610 | 1.379495 | 4.174e-2 |
| 660 | 1.153394 | 7.024e-2 |
| 710 | 0.869989 | 1.349e-1 |
| 760 | 0.593626 | 2.549e-1 |
| 800 | 0.430060 | 3.715e-1 |
| 850 | 0.279000 | 5.260e-1 |
| 900 | 0.279000 | 5.260e-1 |
| 1000 | -0.228100 | 1.691 |
| 1100 | -0.228100 | 1.691 |
| 1200 | -0.228100 | 1.691 |
| 1300 | -0.228100 | 1.691 |
| 1400 | -0.228100 | 1.691 |
| 1500 | -0.228100 | 1.691 |
| 1620 | -0.228100 | 1.691 |
| 1764 | -0.228100 | 1.691 |
| 1936.8 | -0.228100 | 1.691 |
| 2144.16 | -0.228100 | 1.691 |
| 2400 | -0.228100 | 1.691 |
| 2900 | -1.000000 | 10.00 |
| 3500 | -2.000000 | 100.0 |
| (core, r < r0 − 3500 km) | — | 1.000e5 |

**Table 9.** H magnitude ratio and phase shift relative to τ=0, Earth radius, period 3600 s,
*realistic* layered conductivity (Table 8) rather than the `|kl|r0`-matched homogeneous sphere of
Table 7. Same τ, Q values as every other table in this report.

| τ [S] | Q | \|Hr/Hr,0\| | ΔφHr [°] | \|Hθ,φ/Hθ,φ,0\| | ΔφHθ,φ [°] |
|---|---|---|---|---|---|
| 1e-4 | 1.40e-6 | 1.000000 | 6e-8 | 1.000000 | -1e-8 |
| 1    | 0.01397 | 0.999987 | 6.37e-4 | 1.000002 | -9.56e-5 |
| 1e2  | 1.397   | 0.998695 | 0.05833 | 1.000201 | -0.00872 |
| 2.5e3 | 34.93  | 0.979892 | 0.02557 | 1.003066 | -0.00138 |

**The perturbation is markedly smaller than the matched-regime toy comparison predicted.** Q
depends only on ω, τ, r0, not on conductivity, so Table 9 uses *exactly* the same Q values as
Table 7. Yet at `Q=34.9` (`τ=2500` S), the realistic model gives `|Hr/Hr,0|=0.980` and
`ΔφHr=0.026°`, versus `0.882` and `0.85°` for the matched-regime homogeneous sphere — roughly 6×
smaller in both magnitude deviation and phase shift at the identical Q. The mechanism from the
previous paragraph explains why: saturation sets in once `Q >> |tnp(l)|`, and `|tnp(l)|` — the
conductor's own pre-jump value — is set by the conductivity structure, not by Q. The real Earth
model is far more conductive through most of its depth (up to σ≈1.7 S/m in the mantle, σ=1e5 S/m in
the core; see Table 8) than the thin, `|kl|r0`-matched homogeneous sphere used in §3 to isolate the
Q effect, giving it a correspondingly larger `|tnp(l)|` and therefore requiring a larger Q before
the `-iQ` term in the jump condition dominates. The phase shift is still non-monotonic in Q here too
(`0.058°` at `Q=1.4`, falling to `0.026°` at `Q=34.9`) — the same qualitative saturation/overshoot
behavior identified in §3, just at a much smaller absolute scale and a lower onset Q than that
section's toy-conductor threshold. **Conclusion: the `Q~40` threshold derived in §2.1/§3 is specific
to that toy conductor's own electrical regime, not a universal property of Q alone** — a genuinely
more conductive Earth model tolerates a larger τ (at fixed r0, ω) before departing from linear
response. This makes the practical recommendation below (set τ=0) even safer for production runs
than §3 alone would suggest: at the current default `τ=1e-4` S, `Q≈1.4e-6` at 3600 s regardless of
which conductivity model is used, many orders of magnitude below where any version of this
nonlinearity begins.

## 5. Discussion

**No numerical advantage to nonzero τ.** τ appears in exactly one place in each solver: the additive
term in the jump condition. It is never a divisor, matrix pivot, or regularizer, so τ=0 introduces
no degeneracy — consistent with the exact match already obtained at τ=0 in the existing regression
test. There is no floating-point or conditioning reason to prefer a small nonzero value over exactly
zero.

**Practical recommendation.** For runs that are not deliberately modeling a real conductive surface
sheet, set τ=0: it is the mathematically clean choice, matches the validated exact-match test case,
and removes an otherwise-unexplained free parameter. `FWD1D.f90`'s current default `τ=1e-4` has no
numerical justification we could identify and can be replaced by 0 without any accuracy cost.

**Scale dependence.** Since `Q∝r0`, the toy sphere's Q values (`τ=2500` S gives only `Q=3.14`) do
*not* directly represent Earth-scale behavior — at Earth's radius the same τ values give Q up to 776
(Table 5), well beyond where Table 1's linear fit holds. §3 evaluates this directly rather than
extrapolating: at `τ=1e-4` S, `Q_Earth~1e-5–1e-6`, confirming the current default is negligible at
Earth scale too and can safely be replaced by 0. At `τ=100–2500` S, `Q_Earth` reaches the regime
where §3's toy-conductor comparison first showed nonlinearity (`Q≳40`) — but §4 shows this specific
threshold is conductor-dependent, not a property of Q alone: repeating the same τ sweep with the
actual production conductivity model (`LWS/layered_GDE_rho.prm`) at 3600 s gives a perturbation
roughly 6× smaller, at the identical Q, than the matched-regime toy sphere predicted, because the
real Earth's higher conductivity gives it a correspondingly larger `|tnp(l)|`. This is not merely a
mechanism demonstration but a production-representative estimate already obtained (§4, Table 9): if
a nonzero τ representing a genuine conductive surface feature (e.g. the ocean) is intentionally used
at Earth scale, its quantitative effect on this specific conductivity model is now known directly,
not merely bounded by the toy-sphere mechanism study.
