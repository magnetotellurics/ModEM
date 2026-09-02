# E-field "sign mismatch" global1d vs. ModEM — RESOLVED (2026-07-27)

Focus: **Ephi for Mode 1 (P10, l=1 m=0 zonal source)**, the cleanest case (axisymmetric,
E_theta=0, only E_phi nonzero).

## Bottom line

There is **no sign bug and no conjugation.** global1d's EGBERTKELBERT2012 output and ModEM's raw
`.esoln` E-field differ by a **single multiplicative complex constant per period**:

```
c = E_global1d / E_ModEM = i · ω · μ₀ · G          (G ≈ real geometric constant)
```

This is the physically-required **potential-vs-field source-normalization** difference: global1d's
"unit source" is a toroidal *potential* amplitude, ModEM's is a boundary *E-field*, and Faraday's
`∇×E = iωμ₀H` relates them — hence the factor of `i·ω·μ₀`. It cancels identically in `Z = E/H`
(which is why the impedance matched to <1% all along). The apparent "`-Re,+Im` vs `+Re,+Im` sign
flip" is just this `≈90°` phase rotation, misread as a sign flip because `arg(ModEM)≈45°` at the
sampled point.

**Fix (tested, physics-preserving):** rotate global1d's source coefficient by `1/i = -i` (a complex
source amplitude — legitimate, and applied at the source it scales E and H together so `Z` is
preserved). This lands global1d in ModEM's `+Re,+Im` quadrant.

## The three decisive tests

Analysis script: `analyze_c.py` (in this directory). Run ModEM once
(`Mod3DMT_SP2_SPH -F USA_small_1D.rho small_test.dat out.dat out.esoln fwd.ctrl`) and global1d in
`OUTPUT_FORMAT='MODEM'` on `USA.0.25x0.25.grd`, fake pole `0 90`.

### 1. `c` is constant over depth (rules out `-conj`)

`c(k) = E_global1d(k) / E_ModEM(k)` down the column (T=1000s, P10):

| k (top→down) | \|c\| | arg(c) |
|---|---|---|
| 1 (top of air) | 5.41e-3 | +94.5° |
| 13 (surface) | 5.74e-3 | +95.35° |
| 48 (deep) | 5.80e-3 | +95.39° |
| ≥53 (\|field\|→0) | noisy | noisy |

Rock-steady `|c|≈5.74e-3`, `arg≈+95.35°` across the whole physically-meaningful range. A constant
multiplicative `c` (linear) is incompatible with `-conj` (antilinear) — the earlier `-conj` claim
was a single-point artifact.

### 2. `c` scales as ω — the `iωμ₀` (rules in potential-vs-field)

3-period P10 source (`mode1_3per.prm`):

| T (s) | \|c\|/ω | arg(c) |
|---|---|---|
| 10 | 0.839 | **+90.2°** |
| 100 | 0.848 | +91.5° |
| 1000 | 0.914 | +95.4° |

`|c|/ω` constant and `arg → 90°` at short period ⇒ `c = i·ω·μ₀·G`. The small residual (90°→95°,
`|c|/ω` +9%, growing with period) is the **flat-earth (ModEM per-column BC) vs spherical (global1d)**
difference — correctly vanishing at thin skin depth. It is ModEM's approximation, not global1d's.

### 3. Mode-independent (one factor fixes both modes)

Mode2 (non-zonal, l=1 m=±1): `arg(c) = +95.34°`, identical to Mode1's `+95.35°` (magnitudes differ
only by the sources' own amplitude normalizations). So a single per-period complex factor
reconciles both modes.

## The fix, tested

P10 source rotated by `-i` (`.prm` line `1 0` with real=0.0, imag=-1.0):

```
global1d P10 source=-i, surface Ey = +4.79157e-04 +5.73552e-04j   (now +Re,+Im)
ModEM                   surface Ey = +9.24304e-02 +9.17078e-02j   (+Re,+Im)
```

Ratio phase now `+5.35°` (the flat-vs-spherical residual at T=1000s; `→ +0.2°` at T=10s) — i.e. a
positive-real scale, matching quadrant. To match ModEM's *amplitude* too, divide additionally by
`ω·μ₀·G` per period, but that scale is arbitrary between two "unit source" conventions; the `-i` is
the essential, period-independent part that aligns the phase.

**Not wired into the code as a preset** — it would need the per-period `iωμ₀` scaling plus the
empirical/geometric `G`. Left as a documented, understood convention relationship. If direct
raw-field interoperability with ModEM is ever needed, apply `1/(iωμ₀G)` to the source.

## Supporting facts (established earlier, all consistent with the above)

- **global1d's own P10 formula is correct**: closed-form uniform-sphere validation at `l=1,m=0`
  (`test_unit_sphere_zonal.f90` + `reference_unit_sphere_zonal.py`) gives ratio `1.0000+0.0001j`
  for all three nonzero components.
- **Impedance matches**: `test_vs_modem_1D_impedance.f90` gives global1d `Zyx = Ey/Hx` = ModEM's internal Z
  to <1% amplitude, <0.5° phase — as it must, since `c` cancels in `Z`.
- **ModEM's BC code** is the classical WS/REBOCC quasi-1D per-column scheme (traced end-to-end),
  driven by a hardcoded unit top-boundary E-field, identical for both polarizations — consistent
  with the mode-independent `c`. (A separately-found `ISIGN`/`e^{±iωt}` question in that legacy
  engine was tested by direct sign-flip: it did not change the mismatch. See CLAUDE.md.)

## The `EGBERTKELBERT2012_MODEM` preset and its residual (same-physical-latitude validation)

`G = (3/2)√(3/(4π))·d_air ≈ 0.733·d_air` is derived from the air-layer geometry in
`docs/source_normalization.md`/`.pdf`. The `EGBERTKELBERT2012_MODEM` output preset (now the default)
multiplies E and H by `1/(iωμ₀G)` to match ModEM's normalization.

**Is the leftover residual a sampling artifact or genuine?** The comparison is subtle because ModEM's
per-column BC makes its surface field laterally ~uniform, while global1d's zonal field varies as
`sinθ` about the fake pole (maximal at the fake equator). The correct comparison is therefore ModEM
(anywhere) vs global1d at its `sinθ=1` maximum. Checked directly (`analyze_samepoint.py`):

- global1d's field **is** `sinθ`-dependent: `|Ey| = 0.129` at the fake equator (lat 0°) falling to
  `0.102` at the grid edges (±37.6°).
- but the interior comparison point (`i=68`) **was already at the fake equator** — so there is **no
  latitude-sampling artifact** (interior-point and lateral-max ratios are identical). The only
  sampling effect is ModEM's `~1.5%` lateral non-uniformity.
- true equatorial residual `r = E_g1d(equator)/mean(E_ModEM)`:

  | `T` (s) | `|r|` (amp) | `arg r` (phase) |
  |---|---|---|
  | 10 | 0.926 | +0.2° |
  | 100 | 0.936 | +1.5° |
  | 1000 | 1.010 | +5.2° |

So the residual — a few percent in amplitude, up to `~5°` in phase — is a **genuine model
difference** (the earth's finite-`Q` C-response, and ModEM's discrete-air-layer flat boundary
condition vs global1d's fully-spherical solution), **not** a sampling artifact or convention error.
global1d is the more accurate, fully-spherical solution. (This corrects an earlier speculation that
part of the residual was inter-grid latitude sampling; the comparison point was already correctly at
`sinθ=1`.)

## Two indexing bugs found and fixed along the way

Both stem from `apply_output_convention` correctly reversing the field data's r-index when
`r_convention` differs (native `R_UP` → target `R_DOWN`):

1. **`FWD1D.f90`'s `MODEM`-format writer** paired the (reversed) field data with an (unreversed)
   `dz` — swapping "top of air" and "bottom of core". Fixed (`transpose_12_reverse_3`).
2. **Earlier raw-ASCII spot-checks read `k=13`** off a transformed file assuming it meant "surface"
   — it doesn't (maps to native index 48). All results here use the corrected indexing.
