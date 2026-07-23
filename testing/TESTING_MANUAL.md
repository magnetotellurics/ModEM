# Testing manual: `field1d_sunegbert2012.f90` validation

This is a step-by-step walkthrough of the test suite that validates
`EARTH/FWD/field1d_sunegbert2012.f90` (the independent Sun & Egbert 2012 solver) —
first against an independently-derived closed-form solution, then against
`pythonSolver/spherical_em_induction.py`. See `CLAUDE.md` for the full
narrative (bugs found/fixed, derivations) — this file is the practical
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
check), and `test_kelbert2014_vs_sunegbert2012_l1m0` (direct `field1d.f90` vs
`field1d_sunegbert2012.f90` side-by-side, l=1,m=0 — see its own header
comment and `CLAUDE.md` for what this one is investigating). `test_unit_sphere`/
`test_earth_l1mneg1`/`test_kelbert2014_vs_sunegbert2012_l1m0` link against both
`field1d_sunegbert2012.f90` and `field1d.f90`; `test_Tnr_uniform_sphere` links
against `field1d.f90` only. (`test_unit_sphere_kelbert2014.f90` and
`test_unit_sphere_sunegbert2012.f90` — which built the identical model/grid/
source — were merged into the single `test_unit_sphere.f90` above, 2026-07-25;
consolidated 2026-07-23 from a former separate `build_test_Tnr.sh` script.)

---

## 1. Closed-form check: l=2, m=1, single 1000 m sphere — BOTH solvers

Physical setup: uniform 1 Ω·m sphere, radius 1000 m, ω=1 rad/s, unit
external (l=2, m=+1) multipole. `test_unit_sphere.f90` runs BOTH solvers
against this same setup; `reference_unit_sphere.py` prints the matching
closed-form reference in two labeled sections, since the two solvers are
NOT expected to match the closed form (or each other) the same way:

- **`field1d_sunegbert2012.f90` (SUNEGBERT2012)** uses the paper's own literal
  normalization ("coefficient of `r^(l+1)`=1"), so this comparison is
  expected to match the closed form **exactly** (not just in phase).
- **`field1d.f90` (KELBERT2014)** carries its own internal normalization
  (`tni`-referenced, `R0^2/l(l+1)` factors) AND applies a final `conjg()` to
  every assembled component that was designed to compensate for a `+m,-m`
  conjugate-pairing reconstruction trick — but this test, like every other
  single-`(l,m)` test in this suite, sources only ONE of the `+m/-m` pair,
  so that pairing-compensation is never exercised as intended. As analyzed
  2026-07-25 (see `CLAUDE.md`), the result for `m≠0` is NOT a clean scalar/
  phase discrepancy — it reconstructs a genuinely different angular pattern
  (the `m`-flipped one) combined with a conjugated radial function. Printed
  for diagnostic inspection only; do **not** expect a real, let alone
  `1+0i`, ratio here. This is the open "absolute sign"/`conjg()` issue,
  not a test bug — see §`test_kelbert2014_vs_sunegbert2012_l1m0` and
  `CLAUDE.md` for the ongoing investigation.

```
python reference_unit_sphere.py
./test_unit_sphere
```

### Expected (correct) output

`reference_unit_sphere.py`, Section 2 (SUNEGBERT2012):
```
component        (r,theta,phi)                                                       value @ native e^-iwt
Hr (H%z)         r= 1525.00 th=1.8325957146 ph=0.0000000000   +1.7667764256e+03 +7.6668657838e+00j
Hphi (H%x)       r= 2000.00 th=1.8325957146 ph=0.3926990817   -4.5829869330e+02 +1.1087679955e+03j
Htheta (H%y)     r= 2000.00 th=1.9307704850 ph=0.0000000000   +3.4851409098e+03 -2.5980129153e+00j
Etheta (E%y)     r= 1525.00 th=1.8325957146 ph=0.3926990817   -5.3876610799e-01 -2.2590824824e-01j
Ephi (E%x)       r= 1525.00 th=1.9307704850 ph=0.0000000000   +7.3643022448e-03 -1.6970527415e+00j
```

`./test_unit_sphere`, SUNEGBERT2012 block:
```
 SUNEGBERT2012 (field1d_sunegbert2012.f90) output (real, imag) -- compare vs
 reference_unit_sphere.py section 2 ("field1d_sunegbert2012.f90 (SUNEGBERT2012)");
 expect EXACT match, ratio = 1+0i for all five components:
  H%z(i=1,j=5,k=Rr=2) [Hr]     =    1.766772564157752E+03   7.720591661637887E+00
  H%x(i=1,j=5,k=Rs=2) [Hphi]   =   -4.582930754781241E+02   1.108770810024113E+03
  H%y(i=1,j=5,k=Rs=2) [Htheta] =    3.485142218339718E+03  -2.616218597913464E+00
  E%y(i=1,j=5,k=Rr=2) [Etheta] =   -5.387581562728363E-01  -2.259241834997097E-01
  E%x(i=1,j=5,k=Rr=2) [Ephi]   =    7.415908258769074E-03  -1.697049115420322E+00
```

**Pass criterion (SUNEGBERT2012 block only)**: ratio (Fortran / reference)
≈ 1+0i for all 5 components, to ~5 significant figures (residual ~1e-5,
from `field1d_sunegbert2012.f90`'s deliberate `r0+1m` epsilon shift and the
reference's finite-difference `ψ_l'` approximation — both harmless, see
`CLAUDE.md`).

The KELBERT2014 block (Section 1 of the reference script, `./test_unit_sphere`'s
first output block) is diagnostic only — see the note above; do not expect a
passing ratio there.

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
corrections from `CLAUDE.md`'s "Cross-convention comparison rule"). This is
exactly the mistake that produced the original, long-unresolved "extra
factor" mystery earlier in this project. Comparing the Fortran output above
against that naive reference gives:

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

## Summary of what's validated where

| What | How | Where |
|---|---|---|
| `field1d_sunegbert2012.f90` radial + angular + assembly, l=2,m=1 | closed form, exact | §1 |
| `field1d.f90` radial + angular + assembly, l=2,m=1 | closed form, diagnostic only (no clean match expected) | §1 |
| `field1d_sunegbert2012.f90` radial + angular + assembly, l=1,m=-1, Earth scale | pythonSolver cross-check, exact | §2 |
| `field1d.f90`'s `sourcePotential` radial functions only | closed form, exact | §3 |
| Any `.hfield`/`.efield` output's E/H differential consistency | Faraday's law, numerical | §4 |
| `field1d.f90` vs `field1d_sunegbert2012.f90` direct comparison, l=1,m=0 | side-by-side, investigates `conjg()`/absolute-sign issue | `test_kelbert2014_vs_sunegbert2012_l1m0` (see `CLAUDE.md`) |
| `spherical_em_induction.py`'s own `solve_layered` | closed form, exact | `pythonSolver/test_pythonsolver_Rval_Rpval.py` |
| `spherical_em_induction.py`'s own `fields_from_R_general` | Faraday's law, self-consistency | `pythonSolver/test_pythonsolver_faraday.py` |
| `spherical_em_induction.py`'s general-`m` vs `solve_layered`/closed-form regression | internal regression | `pythonSolver/test_validate.py`, `test_general_lm.py` |
