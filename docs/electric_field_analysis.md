# Electric field computation in `sourceField1d`

*Author: Anna Kelbert with Claude, 2026-09.*

## Summary

`sourceField1d` (`EARTH/FWD/field1d.f90`) computes both the magnetic field **H** (poloidal) and
the electric field **E** (toroidal) of the 1D layered-sphere solution, from a single scalar
toroidal potential `T_l(r)`. **H** and **E** are built from opposite vector spherical harmonic
(VSH) patterns of the same potential, related by Faraday's law. This note derives that relation,
gives the exact formulas and grid staggering `sourceField1d` implements, and records the one
physical convention choice (radial direction) that distinguishes this solver from its companion,
`field1d_s2.f90`.

---

## 1. Poloidal and toroidal decomposition

Any smooth vector field **F** on a spherical surface splits into three orthogonal parts:

```
F = F_r r_hat  +  grad_tang P  +  r_hat x grad_tang T
```

`grad_tang P` (the tangential gradient of a scalar `P`) is the **poloidal** part; `r_hat x
grad_tang T` (the curl of a radial vector) is the **toroidal** part. For a magnetotelluric source,
**H** is poloidal (`H_r ≠ 0`) and **E** is toroidal (`E_r = 0`, purely tangential) — the standard
TE-mode picture of deep-Earth induction.

In vector spherical harmonics (degree `l`, order `m`), with `vsharm` computing
`Y(l,m+1)=Y_l^m`, `Yt(l,m+1)=∂Y/∂θ`, `Yp(l,m+1)=(∂Y/∂φ)/sinθ`:

```
Psi_lm = grad_tang Y_lm         ->  (Yt theta_hat + Yp phi_hat)     -- POLOIDAL pattern
Phi_lm = r_hat x grad_tang Y_lm ->  (-Yp theta_hat + Yt phi_hat)    -- TOROIDAL pattern
```

**H** uses `Psi_lm` (theta component `Yt`, phi component `+Yp`); **E** uses `Phi_lm` (theta
component `-Yp`, phi component `+Yt`) — the same pair `(Yt,Yp)`, swapped and one sign flipped.
This is the fundamental Maxwell duality for the TE mode:

```
curl E = i*omega*mu0*H     (Faraday)   -- toroidal curl gives a poloidal field
curl H = sigma*E           (Ampere, quasi-static)
```

---

## 2. The potentials, and the H field

`sourcePotential` returns, at every radius and degree:

| Array | Value | Radii | Drives |
|---|---|---|---|
| `Tnr(k,l)` | `T(r)` | `Rr(1:Nr)`, cell centers | `H_r`, both `E` components |
| `Tnsp(k,l)` | `dT/dr` | `Rs(1:Nr+1)`, cell faces | `H_theta`, `H_phi` |

`Tnsp = d(Tnr)/dr` — the physical radial derivative of `Tnr`.

`sourceField1d` assembles, per component (Fortran coefficient ordering `m=0,+1,-1,+2,-2,...`):

```
H%x(i,j,k) [phi]   = sum_l  Yp(l,.) . coefl . Tnsp(k,l) . R0^2/Rs(k) / l(l+1)
H%y(i,j,k) [theta] = sum_l  Yt(l,.) . coefl . Tnsp(k,l) . R0^2/Rs(k) / l(l+1)
H%z(i,j,k) [r]     = sum_l  Y(l,.)  . coefl . Tnr(k,l)  . R0^2/Rr(k)^2
```

(`H_r` carries no `l(l+1)` denominator — consistent with how `Tnr`/`Tnsp` are normalized in
`sourcePotential`.) Each `l`-sum is built directly from the `m=0` term plus each `+m/-m` pair
(§4) — no post-hoc conjugation of any kind.

---

## 3. Electric field: derivation and formulas

Start from `curl E = i*omega*mu0*H`, `E_r = 0`. The theta-component of Faraday's law gives:

```
-(1/r) d(r E_phi)/dr = i*omega*mu0 * H_theta = i*omega*mu0 * Tnsp(r) * R0^2/(r*l(l+1)) * Yt
```

Integrating (`Tnsp = dTnr/dr`):

```
E_phi = -i*omega*mu0 * Tnr(r) * R0^2 / (r*l(l+1)) * Yt
```

The phi-component gives, by the same steps:

```
E_theta = +i*omega*mu0 * Tnr(r) * R0^2 / (r*l(l+1)) * Yp
```

`E_r=0` throughout (TE mode). The `E_r` (radial) component of Faraday's law is automatically
satisfied — not an independent constraint — since `E_theta`/`E_phi` were built from the SAME
potential `T_l(r)` whose value (`Tnr`) and derivative (`Tnsp`) already solve one consistent
radial ODE; the poloidal/toroidal formalism guarantees the two are compatible by construction.

**Result:**

```
E_theta  =  +i*omega*mu0 * Tnr(r) * R0^2 / (r * l(l+1))  *  Yp
E_phi    =  -i*omega*mu0 * Tnr(r) * R0^2 / (r * l(l+1))  *  Yt
E_r      =   0
```

`E` uses the existing `Tnr` array directly — no separate potential is needed.

### Angular patterns compared

| Component | H pattern | H potential | E pattern | E potential |
|---|---|---|---|---|
| theta (H%y / E%y) | `+Yt` | `Tnsp` at `Rs` | `+Yp` | `i*omega*mu0 * Tnr` at `Rr` |
| phi (H%x / E%x) | `+Yp` | `Tnsp` at `Rs` | `-Yt` | `i*omega*mu0 * Tnr` at `Rr` |
| radial (H%z / E%z) | `+Y` | `Tnr` at `Rr` | `0` | — |

`E_tang` pairs with `H_r` (same potential `T(r)`, different angular pattern); `H_tang` pairs with
`dT/dr` (`Tnsp`). This is the Faraday-Ampere duality.

---

## 4. Coefficient assembly: the `±m` pairing

Fortran coefficient ordering differs from MATLAB's:

| Physical `m` | Fortran `coefl` index | MATLAB `amp` index |
|---|---|---|
| `m=0` | `1` | `1` |
| `m=+m` | `2*m` | `m+1` |
| `m=-m` | `2*m+1` | `l+1+m` |

Each degree-`l` sum adds the `m=0` term directly, then each `+m/-m` pair via the
Condon-Shortley identity `Y_l^{-m} = (-1)^m·conjg(Y_l^m)`:

```fortran
C = (VSH(l,m+1)*coefl(2*m) + (-1)**m*conjg(VSH(l,m+1))*coefl(2*m+1)) / (l*(l+1))
```

for the `l(l+1)`-normalized components (`H_theta`, `H_phi`, `E_theta`, `E_phi`); `H_r` uses the
same pairing without the `l(l+1)` divisor. `E%x`'s accumulated pair carries an overall leading
minus sign, matching the `-Yt` pattern above. The `(-1)^m` factor is required: omitting it
reproduces the physical sum only for even `m`.

No component is conjugated after assembly — each `l`-sum, once the `±m` pairing above is applied
term-by-term, is already the correct physical field.

---

## 5. Radial (source-amplitude) normalization

The formulas in §§2–3 are written directly in terms of `coefl` — the coefficient array as
consumed natively by `sourceField1d` — and `Tnr`/`Tnsp`, the radial potential/derivative solved
once per degree `l`, independent of `m` and of the coefficient values. Because `Tnr`/`Tnsp` never
depend on `coeff`, every formula above (`H_r`, `H_theta`, `H_phi`, `E_theta`, `E_phi`) is exactly
**linear in `coefl`** — a coefficient's absolute scale enters every field component identically.

`field1d.f90` natively expects `coeff(l,m)=1` to mean a **unit external surface radial field**:
`H_r(R0) = Y_l^m(θ,φ)` — the `RADIAL_SURFACE` convention in `EARTH/FWD/output_convention.f90`.
Three other conventions are supported for source coefficients authored differently: an amplitude
of the external potential term `T_l^ext(r) = α_l^T·(r/R0)^(l+1)` (`RADIAL_MULTIPOLE`,
`field1d_s2.f90`'s own native convention, unit `α_l^T=1`); a literal `r^(l+1)` coefficient
(`RADIAL_DIMENSIONAL`, Sun & Egbert's own text convention, `r` in metres); or the classical
external geomagnetic potential's Gauss coefficient (`RADIAL_POTENTIAL`,
`V = R0·Σ_l (r/R0)^l·ε_l^m·Y_l^m`, used for ionospheric-current/Sq sources).

To use coefficients authored in one of these other conventions with the formulas in §§2–3,
rescale them **per degree `l`** (the same factor for every `m` of that degree, before
substituting as `coefl`):

```
coefl_for_formula(l,m) = coeff_in_other_convention(l,m) * A(SURFACE,l) / A(other_convention,l)
```

with the per-degree amplitude `A(l)`, relative to the `MULTIPOLE` baseline (identical to
`output_convention.f90`'s `radial_amplitude()`):

| `radial_norm` | `A(l)` |
|---|---|
| `MULTIPOLE` | `1` |
| `SURFACE` | `R0^2 / (l(l+1))` |
| `DIMENSIONAL` | `R0^(l+1)` (overflows double precision for `l ≳ 44`) |
| `POTENTIAL` | `R0^2 / (l+1)` |

Because this rescale is a single real scalar per degree `l`, applied identically wherever `coefl`
appears in §§2–3, it multiplies `H_r`, `H_theta`, `H_phi`, `E_theta`, `E_phi` all by the **same**
factor for that degree — an overall amplitude change only. It never touches the angular pattern
(`Y`/`Yt`/`Yp`), the relative `H`/`E` component ratios, or the impedance `Z = E/H`
(amplitude-independent by construction). This linearity is exactly why rescaling the source
coefficients before the solve (`rescale_source_radial`) and rescaling the finished `H`/`E`
output by the same factor afterward are equivalent operations.

---

## 6. Grid staggering

`H_r`/`E`'s `Tnr` live at cell-center radii `Rr`; `H_theta`/`H_phi`'s `Tnsp` at cell-face radii
`Rs` (`k` loops over `Nrr=Nr` vs. `Nrs=Nr+1` respectively).

| Component | theta eval | phi eval | radii | j range |
|---|---|---|---|---|
| `H%x` (phi) | node `th(j)` | mid `ph+dp/2` | `Rs` | `2..Nt` |
| `H%y` (theta) | mid `th+dt/2` | node `ph(i)` | `Rs` | `1..Nt` |
| `H%z` (r) | node `th(j)` | node `ph(i)` | `Rr` | `j1..j2` (pole-aware) |
| `E%y` (theta) | node `th(j)` | mid `ph+dp/2` | `Rr` | `j1..j2` (pole-aware) |
| `E%x` (phi) | mid `th+dt/2` | node `ph(i)` | `Rr` | `1..Nt` |
| `E%z` (r) | — | — | — | zero |

Each `E` component shares its angular stagger with the `H` component using the **same VSH
factor** — `E%y` (`Yp`) with `H%x` (also `Yp`, node theta/mid phi), `E%x` (`Yt`) with `H%y` (also
`Yt`, mid theta/node phi) — not with the `H` component it pairs with physically via Faraday's
law. `E%y`, not `E%x`, is therefore the pole-adjacent (node-theta) component and carries the
pole-aware `j1,j2` range: `j1=1` unless `grid%th(1)` is a true pole (`j1=2` then); `j2=Nt+1`
unless `grid%th(Nt+1)` is a true pole (`j2=Nt` then). This range is only ever narrowed for a
genuinely global grid reaching both poles — a regional grid keeps its full node range.
`node_P_lm` (shared with `H%z`) is reused for `E%y`; `edge_P_lm` (shared with `H%y`) is reused
for `E%x` — no separate Legendre computation is needed for `E`.

---

## 7. Scope and the radial-direction convention

`E_r=0` is exact for the TE (poloidal-`H`) source this solver targets; a source with TM-mode
content would need separate treatment (not needed for standard Sq/ring-current MT sources).

`field1d.f90`'s raw output above is natively `e^{-iωt}`, with `r̂` pointing away from Earth's
center (increasing radius) — identical, in every bookkeeping respect but one, to the companion
solver `field1d_s2.f90`. This was a deliberate choice, not the solver's original one:
`field1d.f90` was written in 2011 as a faithful implementation of Kelbert (2006)'s own
conventions (the Global3D Fortran code in ModEM), which use `r̂` pointing toward Earth's center
and `e^{+iωt}` output — genuinely different, equally valid physical conventions. `field1d.f90`
was later re-pointed to match `field1d_s2.f90` so the two independently-derived solvers agree
natively with no compensating sign in comparison code; Kelbert (2006)'s original convention
remains available as the named `KELBERT2006` preset in `EARTH/FWD/output_convention.f90`. See
`docs/output_conventions.md` for the full convention system and the derivation of exactly which
field components that radial-direction choice affects.
