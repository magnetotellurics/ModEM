#!/bin/bash
# =============================================================================
# Self-contained spherical-ModEM comparison for the EGBERTKELBERT2012_MODEM
# normalization (same-physical-latitude check; see zonal_sign_test_report.md).
# NOTE: the spherical ModEM build (Mod3DMT_SP2_SPH) is EXPERIMENTAL -- it must be
# built with -DBUILD_SPHERICAL -DFORCE_SPHERICAL and warns it "will not work with
# traditional MT source setup"; the trusted cross-check is run_cartesian_sanity.sh.
#
# Run in WSL / Linux:
#     wsl bash run_spherical_comparison.sh      # from Windows PowerShell, or
#     bash run_spherical_comparison.sh          # from a WSL shell in this directory
#
# Only external dependencies are the two compiled binaries below. The spherical
# ModEM inputs (USA_small_1D.rho, small_test.dat, fwd.ctrl) are in ../ (testing/);
# the global1d grid/model are in the repo; mode1_3per.prm is in this directory.
# =============================================================================
set -e
cd "$(dirname "$0")"

# -------------------- external binaries: EDIT THESE PATHS --------------------
MODEM_SPH=../../../ModEM/build-serial/Mod3DMT_SP2_SPH  # spherical ModEM (experimental)
FWD1D=../../FWD1D                                      # global1d forward solver (default preset)
# -----------------------------------------------------------------------------
GRID=../../MTsource/USA.0.25x0.25.grd
LAYERED=../../MTsource/layered.prm

echo "[1/3] spherical ModEM forward on USA_small -> out.esoln"
"$MODEM_SPH" -F ../USA_small_1D.rho ../small_test.dat out.dat out.esoln ../fwd.ctrl

echo "[2/3] global1d 3-period, default EGBERTKELBERT2012_MODEM preset, fake pole (0,90) -> normchk.*.esoln"
"$FWD1D" "$LAYERED" mode1_3per.prm "$GRID" normchk 0 90

echo "[3/3] same-physical-latitude analysis (needs a python with numpy+scipy):"
python analyze_samepoint.py \
  || echo ">> .esoln generated. Now run:  python analyze_samepoint.py  (with your scipy-enabled python)"

# =============================================================================
# OPTIONAL: raw c = i*omega*mu0*G measurement (analyze_c.py). This needs global1d
# built with the RAW EGBERTKELBERT2012 convention (the default _MODEM already
# divides by i*omega*mu0*G). Uncomment to run -- it recompiles FWD1D, then restores:
# -----------------------------------------------------------------------------
# sed -i "s/'EGBERTKELBERT2012_MODEM'  ! DEFAULT/'EGBERTKELBERT2012'/" ../../FWD1D.f90 || true
# ( cd ../.. && make -f Makefile1d FWD1D >/dev/null )
# "$FWD1D" "$LAYERED" mode1_3per.prm "$GRID" g1d3_raw 0 90
# python analyze_c.py g1d3_raw.s1.EGBERTKELBERT2012.esoln 1 || echo "run: python analyze_c.py g1d3_raw.s1.EGBERTKELBERT2012.esoln 1"
# git -C ../.. checkout -- FWD1D.f90 ; ( cd ../.. && make -f Makefile1d FWD1D >/dev/null )
# =============================================================================
