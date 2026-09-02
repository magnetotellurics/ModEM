# Testing manual: `field1d_s2.f90` validation

This is a step-by-step walkthrough of the test suite that validates
`EARTH/FWD/field1d_s2.f90` (the independent Sun & Egbert 2012 solver) —
first against an independently-derived closed-form solution, then against
`pythonSolver/spherical_em_induction.py`. This is the practical
"how do I run it and what should I see" reference.

All commands below assume `gfortran` and a Python environment with
`numpy`/`scipy` are on your PATH, and that you run from the `testing/`
directory unless stated otherwise.

## 0. Build

```
make -f Makefile.test1d
```

builds all four test executables: `test_unit_sphere` (§1, l=2,m=1 closed-form
check for BOTH solvers), `test_earth_l1mneg1` (§2, l=1,m=-1 vs pythonSolver
check), `test_Tnr_uniform_sphere` (§3, `field1d.f90`'s `sourcePotential`
check), and `test_s1_vs_s2_l1m0` (direct `field1d.f90` vs
`field1d_s2.f90` side-by-side, l=1,m=0 — see its own header comment; both
solvers agree up to a single real normalization constant, per §1's finding).
`test_unit_sphere`/
`test_earth_l1mneg1`/`test_s1_vs_s2_l1m0` link against both
`field1d_s2.f90` and `field1d.f90`; `test_Tnr_uniform_sphere` links
against `field1d.f90` only. (`test_unit_sphere_s1.f90` and
`test_unit_sphere_s2.f90` — which built the identical model/grid/
source — were merged into the single `test_unit_sphere.f90` above, 2026-07-25;
consolidated 2026-07-23 from a former separate `build_test_Tnr.sh` script.)

---

## 1. Closed-form check: l=2, m=1, single 1000 m sphere — BOTH solvers

Physical setup: uniform 1 Ω·m sphere, radius 1000 m, ω=1 rad/s, unit
external (l=2, m=+1) multipole. `test_unit_sphere.f90` runs BOTH solvers
against this same setup; `reference_unit_sphere.py` prints the matching
closed-form reference in two labeled sections, since the two solvers use
different internal normalizations:

- **`field1d_s2.f90` (S2)** uses the paper's own literal normalization
  ("coefficient of `r^(l+1)`=1"), so this comparison matches the closed
  form **exactly** (not just in phase).
- **`field1d.f90` (S1)** carries its own internal normalization
  (`tni`-referenced, `R0^2/l(l+1)` factors) and, separately, natively
  outputs `H_r`/`E_theta`/`E_phi` with the opposite sign convention from S2
  (S1's `r̂` points toward Earth's center, S2's points away — two different
  but equally valid physical conventions). Its ratio against the closed
  form is a single **real** constant — **≈+1.6667e5 for the three H
  components, ≈−1.6667e5 for the two E components** (the sign split is the
  ordinary H-vs-E pseudovector/polar-vector distinction) — not `1+0i`,
  since S1's normalization was never calibrated to match the closed form's
  amplitude.

```
python reference_unit_sphere.py
./test_unit_sphere
```

### Expected output (regenerated 2026-09)

`reference_unit_sphere.py`:
```
=== Section 1: field1d.f90 (S1) ===
component        (r,theta,phi)                                                           value @ s=-1 (+iwt)
Hr (H%z)         r= 1525.00 th=1.8325957146 ph=0.0000000000   +1.7614866793e-06 -7.6439111111e-09j
Hphi (H%x)       r= 2000.00 th=1.8325957146 ph=0.3926990817   -4.5857415588e-07 +1.1047658692e-06j
Htheta (H%y)     r= 2000.00 th=1.9307704850 ph=0.0000000000   +3.4747063632e-06 +2.5902344387e-09j
Etheta (E%y)     r= 1525.00 th=1.8325957146 ph=0.3926990817   +5.3908754226e-10 +2.2056156756e-10j
Ephi (E%x)       r= 1525.00 th=1.9307704850 ph=0.0000000000   +7.3422534503e-12 +1.6919717486e-09j

=== Section 2: field1d_s2.f90 (S2) -- expect EXACT match, ratio=1+0i ===
component        (r,theta,phi)                                                         value @ native e^-iwt
Hr (H%z)         r= 1525.00 th=1.8325957146 ph=0.0000000000   +1.7614866793e-06 +7.6439111111e-09j
Hphi (H%x)       r= 2000.00 th=1.8325957146 ph=0.3926990817   -4.5692654243e-07 +1.1054483330e-06j
Htheta (H%y)     r= 2000.00 th=1.9307704850 ph=0.0000000000   +3.4747063632e-06 -2.5902344387e-09j
Etheta (E%y)     r= 1525.00 th=1.8325957146 ph=0.3926990817   -5.3715303688e-10 -2.2523187669e-10j
Ephi (E%x)       r= 1525.00 th=1.9307704850 ph=0.0000000000   +7.3422534503e-12 -1.6919717486e-09j
```

`./test_unit_sphere`:
```
 S1 (field1d.f90) output (real, imag) -- compare vs reference_unit_sphere.py section 1:
  H%z(i=1,j=5,k=Rr=2) [Hr]     =    2.935804715788943E-01   1.282912668115673E-03
  H%x(i=1,j=5,k=Rs=2) [Hphi]   =   -7.615349024114372E-02   1.842418565234570E-01
  H%y(i=1,j=5,k=Rs=2) [Htheta] =    5.791179446223432E-01  -4.347309285245833E-04
  E%y(i=1,j=5,k=Rr=2) [Etheta] =   -8.952418337538885E-05  -3.754129342278147E-05
  E%x(i=1,j=5,k=Rr=2) [Ephi]   =    1.232284119260964E-06  -2.819946843171999E-04

 S2 (field1d_s2.f90) output (real, imag) -- compare vs
 reference_unit_sphere.py section 2 ("field1d_s2.f90 (S2)");
 expect EXACT match, ratio = 1+0i for all five components:
  H%z(i=1,j=5,k=Rr=2) [Hr]     =    1.761482829459403E-06   7.697476133112672E-09
  H%x(i=1,j=5,k=Rs=2) [Hphi]   =   -4.569209414340766E-07   1.105451139147801E-06
  H%y(i=1,j=5,k=Rs=2) [Htheta] =    3.474707667738790E-06  -2.608385613308313E-09
  E%y(i=1,j=5,k=Rr=2) [Etheta] =   -5.371451089734440E-10  -2.252477642384539E-10
  E%x(i=1,j=5,k=Rr=2) [Ephi]   =    7.393704955394321E-12  -1.691968133423683E-09
```

**Pass criterion, S2 block**: ratio (Fortran / reference) ≈ 1+0i for all 5
components, to ~5 significant figures (residual ~1e-5, from
`field1d_s2.f90`'s deliberate `r0+1m` epsilon shift and the reference's
finite-difference `ψ_l'` approximation — both harmless).

**Pass criterion, S1 block**: ratio (Fortran / reference Section 1) is a
single real constant, no phase residual beyond a fraction of a degree —
**≈+1.6667e5 for `Hr`/`Hphi`/`Htheta`, ≈−1.6667e5 for `Etheta`/`Ephi`**.
Not `1+0i` (S1's normalization isn't calibrated to the closed form's
amplitude) and not yet an automated pass/fail check in the driver — verify
by eye, or compute the ratio directly from the two blocks above.

---

## 2. Cross-check against pythonSolver: l=1, m=-1, Earth scale

Physical setup: uniform 100 Ω·m Earth (r0=6.371e6 m), T=1000 s, unit
external (l=1, m=-1) multipole — the historically problematic case from
earlier in this project's debugging history.

```
python reference_earth_l1mneg1.py
./test_earth_l1mneg1
```

### Expected (correct) output

`reference_earth_l1mneg1.py`:
```
solve_layered: A (raw external amplitude for this K0) = (1+1.1747755004857644e-24j)

component        (r,theta,phi)                                            value @ native e^-iwt (fully corrected)
H%z [Hr]         r= 6372250.000 th=1.8325957146 ph=0.0000000000   +2.5388253699e-02 +2.4371129852e-02j
H%x [Hphi]       r= 6373000.000 th=1.8325957146 ph=0.3926990817   -4.0322158396e-01 -9.4050895938e-01j
H%y [Htheta]     r= 6373000.000 th=1.9307704850 ph=0.0000000000   -3.6043042032e-01 +4.4422153984e-03j
E%y [Etheta]     r= 6372250.000 th=1.8325957146 ph=0.3926990817   +8.5377953726e-04 +3.3337246535e-04j
E%x [Ephi]       r= 6372250.000 th=1.9307704850 ph=0.0000000000   -2.2358130622e-04 +2.3291242380e-04j
```

`./test_earth_l1mneg1`:
```
 FIELD1D_SE12 output (real, imag):
  H%z(i=1,j=5,k=Rr=2) [Hr]     =    2.538794744247850E-02   2.437113760124664E-02
  H%x(i=1,j=5,k=Rs=2) [Hphi]   =   -4.032216483037955E-01  -9.405091042585518E-01
  H%y(i=1,j=5,k=Rs=2) [Htheta] =   -3.604304761582562E-01   4.442216810940462E-03
  E%y(i=1,j=5,k=Rr=2) [Etheta] =    8.537722454616832E-04   3.333757041459036E-04
  E%x(i=1,j=5,k=Rr=2) [Ephi]   =   -2.235813773162212E-04   2.329096142040912E-04
```

**Pass criterion**: same as step 1, ratio ≈ 1+0i for all 5 components to
~5 significant figures.

### What goes wrong if you get the sign convention wrong

`reference_earth_l1mneg1.py` has `SHOW_NAIVE_COMPARISON` at the top — set it
`True` to also print the WRONG reference values (naive
`conj(pythonSolver(l,+m)/A)`, with none of the norm/m-flip/H-sign
corrections listed below). This is exactly the mistake that produced the
original, long-unresolved "extra factor" mystery earlier in this project.
Comparing the Fortran output above against that naive reference gives:

| component | \|ratio\| | phase (deg) | after applying only the √(l(l+1)) norm fix |
|---|---|---|---|
| Hr     | 1.414205 | 0°    | 0.999994 |
| Hphi   | 1.414214 | 135°  | 1.000000 |
| Htheta | 1.414214 | -0°   | 1.000000 |
| Etheta | 1.414205 | -45°  | 0.999994 |
| Ephi   | 1.414205 | -180° | 0.999994 |

Two tells that something's wrong in the comparison methodology (not the
solver itself), if you ever see this pattern again:
- **A uniform `√(l(l+1))` magnitude offset across every component** →
  missing pythonSolver's own `T=R(r)·Y_l^m/√(l(l+1))` normalization (fix:
  multiply the reference by `√(l(l+1))`, on top of dividing by
  `solve_layered`'s `A`).
- **A phase offset that depends on which angular point/component you're
  looking at, clustering at 0°/±45°/±90°/135°/180°** (never a clean, single
  real or imaginary ratio across all 5 components at once) → conjugating
  pythonSolver's `(l,+m)` output instead of its `(l,-m)` output (fix: flip
  the sign of `m` before calling `fields_from_R_general`/`C_lm`, in addition
  to conjugating).
- Once both of the above are fixed, a **clean, real ±1 split between H and E
  components** (not yet ±1 exactly, but real and consistent within each
  field type) is the remaining, expected, physical H-vs-E pseudovector sign
  — bake it in (sign=-1 for H, +1 for E) to reach ratio=+1 for everything.

---

## 3. `field1d.f90`'s own Tnr/Tnsp check (uniform sphere, m=0)

```
make -f Makefile.test1d test_Tnr_uniform_sphere   # already built by step 0's plain `make -f Makefile.test1d`
python uniform_sphere_Tnr_predict.py
./test_Tnr_uniform_sphere
```

Validates `field1d.f90`'s `sourcePotential`'s `Tnr`/`Tnsp`/`Tnrp`/`Tns`
(l=1, uniform 100 Ω·m sphere, T=1000s) against an independent closed-form
ψ₁-based prediction. Last confirmed-passing values (from this project's
history): `Tnr=(3.747126300279E-02, 3.653564680352E-02)`,
`Tnsp=(4.649991236165E-07, -5.731721451757E-09)`, with `Tnrp≡Tnsp` and
`Tns≡Tnr` (expected, since the single test radius sits exactly on both the
cell-centre and face grids).

---

## 4. Faraday's-law self-consistency, real `.hfield`/`.efield` output

```
python test_hfield_efield_faraday.py
```

Checks `(curl E)_r = i·ω·μ0·H_r` via staggering-aware finite differences on
ACTUAL Fortran-written `.hfield`/`.efield` files (edit the `EFIELD_FILE`/
`HFIELD_FILE`/`R_INDEX`/`I_H`/`J_H` constants at the top to point at
whichever output you want to check). Note: this only validates the
*differential* consistency between E and H at the SAME (l,m) — it does not
catch an overall wrong normalization or sign shared by both (that needs
steps 1–2's closed-form/cross-solver checks instead).

---

## 5. `earth%tau` (thin-sheet surface conductance) sensitivity

```
make -f Makefile.test1d test_tau_sweep
./test_tau_sweep
```

Sweeps `earth%tau` from 1 S down to 0 and reports the magnitude error and phase shift of the
S2-solver output relative to the `tau=0` case (which is itself exact, since S2 matches the
closed-form solution exactly at `tau=0` — see §1). Also reports the same comparison at Earth's
actual radius, at two production periods, both for a homogeneous sphere matched to the toy
sphere's own electrical regime and for the actual layered conductivity model in
`LWS/layered_GDE_rho.prm`. No pass/fail criterion — this is a sensitivity study, not a
regression check; every number it prints is quoted directly in
`docs/tau_sensitivity_analysis.md`/`.pdf`, which explains what they mean.

---

## Summary of what's validated where

| What | How | Where |
|---|---|---|
| `field1d_s2.f90` radial + angular + assembly, l=2,m=1 | closed form, exact | §1 |
| `field1d.f90` radial + angular + assembly, l=2,m=1 | closed form, clean real ratio (S1's own normalization, not calibrated to 1) | §1 |
| `field1d_s2.f90` radial + angular + assembly, l=1,m=-1, Earth scale | pythonSolver cross-check, exact | §2 |
| `field1d.f90`'s `sourcePotential` radial functions only | closed form, exact | §3 |
| Any `.hfield`/`.efield` output's E/H differential consistency | Faraday's law, numerical | §4 |
| `earth%tau` thin-sheet sensitivity, toy sphere + Earth scale (matched-regime + realistic model) | magnitude/phase sweep, sensitivity study | §5 |
| `field1d.f90` vs `field1d_s2.f90` direct comparison, l=1,m=0 | side-by-side, real single-constant ratio (S1's normalization, `H`/`E` opposite sign) | `test_s1_vs_s2_l1m0` |
| `spherical_em_induction.py`'s own `solve_layered` | closed form, exact | `pythonSolver/test_pythonsolver_Rval_Rpval.py` |
| `spherical_em_induction.py`'s own `fields_from_R_general` | Faraday's law, self-consistency | `pythonSolver/test_pythonsolver_faraday.py` |
| `spherical_em_induction.py`'s general-`m` vs `solve_layered`/closed-form regression | internal regression | `pythonSolver/test_validate.py`, `test_general_lm.py` |
