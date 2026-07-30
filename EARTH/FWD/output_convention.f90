! *****************************************************************************
module output_convention
    ! Composable output-convention layer for the 1D layered-sphere solvers
    ! (field1d.f90/S1, field1d_s2.f90/S2).
    ! Neither solver's internal numerics (Legendre/vsharm recursions, radial
    ! solve) are touched anywhere in this module. All convention handling is
    ! either (a) a rescale of the source spherical-harmonic coefficients,
    ! applied BEFORE the solver runs (rescale_source_coeffs), or (b) a
    ! transform of the finished output cvectors, applied AFTER the solver
    ! runs (apply_output_convention). This module does not `use` either
    ! solver module -- solvers/drivers depend on this module, never the
    ! reverse.
    !
    ! Five field-bookkeeping dimensions are tracked (grid scope -- regional
    ! vs global, poles, zero-longitude wrap -- is explicitly out of scope;
    ! that is already handled by the existing primary_grid flag and the
    ! separate "fake pole" grid-recentering option in FWD1D.f90):
    !   1) time convention        (e^{+iwt} vs e^{-iwt})
    !   2) harmonic normalization (Schmidt semi-norm vs fully-normalized)
    !   3) Condon-Shortley phase  (included vs excluded)
    !   4) theta convention       (colatitude N->S vs latitude S->N)
    !   5) r convention           (down/depth-increasing vs up/radius-increasing)
    ! plus primary_grid ('H' or 'E'), recorded for bookkeeping/validation only
    ! -- it is not touched by either transform below, since both solvers
    ! already accept primary_grid as an independent flag.
    !
    ! ==========================================================================
    ! NATIVE CONVENTIONS, as of 2026-07-26 (field1d.f90 fixed to match
    ! Sun & Egbert (2012) convention exactly):
    !
    ! Both field1d.f90 (post the 2026-07-26 paired-conjg() removal and
    ! T(r)-formula sign fix) and field1d_s2.f90 now produce
    ! IDENTICAL raw output (time=e^{-iwt}, norm=fully-normalized [vsharm's Y
    ! is independently sympy-verified to match the standard orthonormal
    ! Y_l^m], Condon-Shortley INCLUDED [vsharm also verified to match the
    ! CS-phase Y_l^m exactly],
    ! theta=colatitude North->South) -- confirmed via
    ! testing/test_s1_vs_s2_l1m0.f90 giving ratio=+1,
    ! angle=0deg for every nonzero component. MT-relabeling-wise,
    ! Hz(down)=-H%z(raw) for BOTH solvers now (confirmed 2026-07-26; this was
    ! Hz(down)=+H%z(raw) before the T(r)-formula sign fix -- the relationship
    ! flipped along with the formula).
    !
    ! primary_grid (2026-07-27, per user direction): S1's history predates S2
    ! and originally defaulted to 'H' staggering, but field1d.f90 has always
    ! worked correctly with EITHER primary_grid (both EDGE and FACE branches
    ! are fully implemented and validated this session) -- there is no
    ! technical reason for S1 and S2 to default to different staggering, so
    ! S1's native primary_grid is now also 'E', matching S2 exactly. The two
    ! solvers' native conventions are therefore now IDENTICAL in every
    ! dimension including primary_grid -- native_convention() below simply
    ! returns the SUNEGBERT2012 preset for either solver name, with no
    ! separate S1-specific construction (the earlier version of this module
    ! built S1's native convention by hand with primary_grid='H'; that is no
    ! longer correct and has been removed, not just left stale).
    !
    ! IMPORTANT: per explicit user direction (2026-07-26), the KELBERT2006
    ! preset below is kept EXACTLY as originally specified (time=+iwt,
    ! norm=Schmidt-seminorm, Condon-Shortley=off, r=down) even though this NO
    ! LONGER describes field1d.f90's actual native output -- it is a fixed,
    ! deliberately-named target convention a user may request, not a claim
    ! about what any solver produces raw. Consequently native_convention()
    ! below does NOT alias S1's native convention to the KELBERT2006
    ! preset (unlike the convenience shown in the original implementation
    ! spec, which predates the T(r)-formula sign fix and assumed r=down was
    ! still field1d.f90's genuine native behavior) -- it constructs
    ! field1d.f90's TRUE native convention explicitly instead. Producing
    ! genuine KELBERT2006-convention output from EITHER solver now requires
    ! applying the full time+r (and, per rescale_source_coeffs, norm+CS)
    ! transform, same as any other non-native target.
    ! ==========================================================================

    use math_constants, only: prec, PI, MU_0
    use sg_vector, only: cvector

    implicit none

    private
    public :: output_convention_t
    public :: TIME_POSITIVE, TIME_NEGATIVE, NORM_SCHMIDT, NORM_FULL, &
              THETA_COLAT, THETA_LAT, R_DOWN, R_UP
    public :: KELBERT2006, SUNEGBERT2012, EGBERTKELBERT2012, EGBERTKELBERT2012_MODEM, LWS
    public :: RADIAL_SURFACE, RADIAL_MULTIPOLE, RADIAL_DIMENSIONAL, RADIAL_POTENTIAL
    public :: get_convention, native_convention, native_radial
    public :: rescale_source_coeffs, apply_output_convention, apply_modem_normalization
    public :: rescale_source_radial, radial_amplitude

    ! ------- enumerated tokens (character params keep the .prm/config human-readable) -------
    character(len=10), parameter :: TIME_POSITIVE = 'PLUS_IWT'
    character(len=10), parameter :: TIME_NEGATIVE = 'MINUS_IWT'
    character(len=16), parameter :: NORM_SCHMIDT   = 'SCHMIDT_SEMINORM'
    character(len=16), parameter :: NORM_FULL      = 'FULLY_NORM'
    character(len=10), parameter :: THETA_COLAT = 'COLAT_N2S'  ! colatitude, North->South
    character(len=10), parameter :: THETA_LAT   = 'LAT_S2N'    ! latitude,   South->North
    character(len=6),  parameter :: R_DOWN = 'R_DOWN'          ! r increases downward (depth)
    character(len=6),  parameter :: R_UP   = 'R_UP'            ! r increases upward (radius)
    ! Per-degree RADIAL (source-amplitude) normalization -- i.e. what
    ! coeff(l,m)=1 means for the radial part of the field. This is the ONE
    ! convention dimension where the two solvers genuinely differ natively
    ! (S1=SURFACE, S2=MULTIPOLE, see native_radial); reconciling it is what
    ! makes S1 and S2 produce identical output for the same OUTPUT_CONVENTION.
    ! See radial_amplitude() for the exact per-degree factors -- that single
    ! function is the place to adjust or add a per-degree radial normalization.
    !
    ! Two distinct FORMULATIONS are represented here, not just four scale
    ! choices:
    !   TOROIDAL POTENTIAL (T) family -- RADIAL_MULTIPOLE, RADIAL_SURFACE,
    !     RADIAL_DIMENSIONAL. All three parameterize coeff(l,m) as an
    !     amplitude of the solvers' own toroidal potential T_l(r), whose
    !     source (regular-at-origin) term in the insulating atmosphere is
    !     T_l^ext(r) = alpha_l^T * (r/R0)^(l+1) (Sun & Egbert eq 6) -- they
    !     differ only in what numeric scale alpha_l^T=1 corresponds to.
    !   EXTERNAL (SCALAR) POTENTIAL (V) family -- RADIAL_POTENTIAL. A
    !     genuinely different formulation: coeff(l,m) parameterizes the
    !     classical external/inducing geomagnetic potential V instead of T,
    !     i.e. the Gauss-coefficient convention V = R0*sum_l (r/R0)^l *
    !     epsilon_l^m * Y_l^m used throughout classical geomagnetism and
    !     ionospheric-current (Sq/LWS-style) source modeling. See
    !     radial_amplitude() for the T<->V derivation and A(l).
    character(len=12), parameter :: RADIAL_MULTIPOLE   = 'MULTIPOLE'   ! S2 native; Sun & Egbert eq (6), unit (r/a)^(l+1) potential coeff; A = 1
    character(len=12), parameter :: RADIAL_SURFACE     = 'SURFACE'     ! S1 native; unit surface radial field; A = R0^2/(l(l+1))
    character(len=12), parameter :: RADIAL_DIMENSIONAL = 'DIMENSIONAL' ! Sun & Egbert text, unit r^(l+1) (r in metres); A = R0^(l+1) -- OVERFLOWS for l >~ 44
    character(len=12), parameter :: RADIAL_POTENTIAL   = 'V_POTENTIAL' ! external scalar potential V (Gauss coefficient) parameterization; A = R0^2/(l+1)

    ! ------- the convention type: five field-bookkeeping dimensions + primary_grid -------
    type :: output_convention_t
        character(len=24) :: name               ! 'KELBERT2006' etc, for logging
        character(len=10) :: time_convention    ! TIME_POSITIVE | TIME_NEGATIVE
        character(len=16) :: harmonic_norm      ! NORM_SCHMIDT | NORM_FULL
        logical            :: condon_shortley   ! .true. = Condon-Shortley phase included
        character(len=10) :: theta_convention   ! THETA_COLAT | THETA_LAT
        character(len=6)  :: r_convention       ! R_DOWN | R_UP
        character(len=1)  :: primary_grid       ! 'H' | 'E'  (records intended staggering)
        ! .true. = additionally rescale the OUTPUT fields by 1/(i*omega*mu0*G)
        ! to match ModEM's boundary-E-field source normalization instead of
        ! global1d's toroidal-potential normalization (see apply_modem_norm-
        ! alization and docs/source_normalization.pdf). Purely a field-value
        ! rescale of the finished output; does NOT change any of the seven
        ! bookkeeping dimensions above. Default .false. for every preset
        ! except EGBERTKELBERT2012_MODEM.
        logical            :: modem_normalize = .false.
        ! Per-degree RADIAL source-amplitude normalization (RADIAL_MULTIPOLE |
        ! RADIAL_SURFACE | RADIAL_DIMENSIONAL). Reconciles the two solvers'
        ! different native radial normalizations (S1=SURFACE, S2=MULTIPOLE) so
        ! that BOTH produce identical output at the chosen scale. Default
        ! SURFACE (O(1) fields, S1-native, ModEM-compatible). Applied to the
        ! source coeff per degree before the solve (rescale_source_radial).
        character(len=12)  :: radial_norm = RADIAL_SURFACE
    end type output_convention_t

    ! ------- the predefined conventions (fixed, named targets) -------
    !
    ! Per-preset radial_norm choices:
    !   SUNEGBERT2012          RADIAL_MULTIPOLE -- so this preset IS S2's EXACT
    !                          native convention (S2 + SUNEGBERT2012 = identity;
    !                          see native_convention). Sun & Egbert eq-(6),
    !                          unit (r/a)^(l+1) potential-coefficient scale.
    !   KELBERT2006,           RADIAL_SURFACE   -- S1-native / O(1) / ModEM-
    !   EGBERTKELBERT2012,     compatible surface-field scale. (Under
    !   EGBERTKELBERT2012_MODEM  modem_normalize, the radial target is forced to
    !                          SURFACE anyway -- see FWD1D.f90 -- so the value
    !                          here is what a NON-modem variant would use.)
    !   LWS                    RADIAL_POTENTIAL (external scalar potential V /
    !                          Gauss-coefficient parameterization -- a genuinely
    !                          DIFFERENT formulation than the other presets'
    !                          toroidal-potential T family, see the RADIAL_*
    !                          declarations above), AND condon_shortley=.false.
    !                          (see its own note below).
    type(output_convention_t), parameter :: KELBERT2006 = output_convention_t( &
        'KELBERT2006', TIME_POSITIVE, NORM_SCHMIDT, .false., THETA_COLAT, R_DOWN, 'H', .false., RADIAL_SURFACE)

    ! SUNEGBERT2012 = S2's EXACT native convention (see native_convention):
    ! MINUS_IWT, FULLY_NORM, Condon-Shortley included, colatitude N->S, r up,
    ! E-primary, no ModEM rescale, MULTIPOLE radial. Selecting this preset with
    ! SOLVER='S2' is a mechanical identity (all convention transforms are no-ops).
    type(output_convention_t), parameter :: SUNEGBERT2012 = output_convention_t( &
        'SUNEGBERT2012', TIME_NEGATIVE, NORM_FULL, .true., THETA_COLAT, R_UP, 'E', .false., RADIAL_MULTIPOLE)

    ! condon_shortley=.true. here (2026-07-26, per direct user instruction, for
    ! validation purposes): confirmed empirically (test_diagnose_egbertkelbert.f90)
    ! that ModEM's own .prm source-coefficient files are ALREADY written in
    ! vsharm's native CS-included basis (they reproduce ModEM correctly when fed
    ! UNSCALED into either solver, see test_vs_modem_1D_impedance.f90) -- .false. here
    ! (the original spec's assumption) spuriously flips every odd-m source
    ! (e.g. Mode2, l=1 m=+-1) via rescale_source_coeffs's (-1)^m term while
    ! leaving m=0 sources (e.g. Mode1) untouched, producing a mode-dependent
    ! sign inconsistency. NOTE (2026-07-26, A. Kelbert): the historical ModEM
    ! solver this convention is named for never actually used spherical
    ! harmonics at all; a future CUSTOM/named preset with condon_shortley=
    ! .false. may still be needed for that use case -- track separately,
    ! do not silently revert this field back to .false.
    type(output_convention_t), parameter :: EGBERTKELBERT2012 = output_convention_t( &
        'EGBERTKELBERT2012', TIME_NEGATIVE, NORM_FULL, .true., THETA_LAT, R_DOWN, 'E', .false., RADIAL_SURFACE)

    ! EGBERTKELBERT2012_MODEM: identical field bookkeeping to EGBERTKELBERT2012,
    ! but additionally rescales the OUTPUT E and H by 1/(i*omega*mu0*G) so the
    ! amplitude/phase match ModEM's boundary-E-field source normalization
    ! rather than global1d's toroidal-potential normalization. The two differ
    ! by exactly c = i*omega*mu0*G (G a real geometric constant set by the
    ! air-layer thickness) -- the potential-vs-field relationship of Faraday's
    ! law; see apply_modem_normalization and docs/source_normalization.pdf for
    ! the full derivation. Applying the SAME complex factor to both E and H
    ! preserves Z=E/H, so no physics is changed -- only the (arbitrary) source-
    ! normalization convention. Intended for direct raw-field interoperability
    ! with ModEM (e.g. feeding global1d fields into ModEM as boundary values);
    ! the match is exact for the l=1/P10-at-equator MT modes in the deep-
    ! induction limit, with a few-% amplitude / <~5deg phase residual (the
    ! earth's C-response and the flat-vs-spherical difference) -- both genuine
    ! model differences, not convention artifacts. VALIDATED against trusted
    ! CARTESIAN ModEM (rho=100 halfspace): the air-column formula matches to
    ! 0.2-0.6% at all periods, and this preset gives |ratio|=1.006 (perfect)
    ! at short period rising cleanly to 1.145 at 3981s -- confirming the
    ! residual is genuine flat-vs-spherical physics, not an artifact of the
    ! experimental spherical ModEM BC. See docs/source_normalization.md/.pdf
    ! sec.6 and testing/test_vs_modem_1D/cartesian_sanity_check.md.
    type(output_convention_t), parameter :: EGBERTKELBERT2012_MODEM = output_convention_t( &
        'EGBERTKELBERT2012_MODEM', TIME_NEGATIVE, NORM_FULL, .true., THETA_LAT, R_DOWN, 'E', .true., RADIAL_SURFACE)

    ! LWS (NASA "Living With a Star" project): S1's native convention (see
    ! native_convention('S1'): MINUS_IWT, FULLY_NORM, colatitude N->S, r up,
    ! E-primary) but with TWO deviations from native, both testing hypotheses
    ! about how the LWS-style (ionospheric-current / Sq) source coefficients
    ! this project uses are actually authored:
    !
    !   (1) radial_norm = RADIAL_POTENTIAL, not RADIAL_SURFACE. This is a
    !       genuinely different FORMULATION, not just a different scale: LWS/
    !       Sq-current source coefficients are conventionally given as the
    !       classical external geomagnetic potential's Gauss coefficients
    !       (coeff(l,m) = epsilon_l^m, parameterizing V = R0*sum_l (r/R0)^l *
    !       epsilon_l^m * Y_l^m), NOT as an amplitude of the solvers' own
    !       toroidal potential T. RADIAL_POTENTIAL performs exactly that
    !       reinterpretation (A(l) = R0^2/(l+1), see radial_amplitude's
    !       derivation) so a coeff(:) authored as external-potential Gauss
    !       coefficients drives the solver correctly.
    !   (2) condon_shortley = .false. It is NOT a solver-native output for
    !       either solver (both solvers are natively CS-included, via
    !       vsharm) -- the output layer flips the (-1)^m Condon-Shortley phase
    !       (rescale_source_coeffs / apply_norm_cs) so the fields come out in
    !       the non-CS basis, testing the hypothesis that the LWS / MATLAB-SIEM
    !       (@TSModel, tanField.m / radField.m) .prm source coefficients are
    !       authored in the non-CS convention.
    !
    ! Run OUTPUT_CONVENTION='LWS' and compare against the MATLAB/SIEM output to
    ! test both hypotheses. As with every preset, the output is defined by the
    ! CONVENTION alone: S1+LWS and S2+LWS produce identical fields (both
    ! transforms are applied identically for both solvers). The solver
    ! internals are unchanged and remain CS-included / T-parameterized.
    type(output_convention_t), parameter :: LWS = output_convention_t( &
        'LWS', TIME_NEGATIVE, NORM_FULL, .false., THETA_COLAT, R_UP, 'E', .false., RADIAL_POTENTIAL)

contains

    ! resolve a name string -> convention (for a config-file / CUSTOM interface later)
    function get_convention(name) result(conv)
        character(*), intent(in)   :: name
        type(output_convention_t)  :: conv
        select case (trim(name))
            case ('KELBERT2006');           conv = KELBERT2006
            case ('SUNEGBERT2012');         conv = SUNEGBERT2012
            case ('EGBERTKELBERT2012');     conv = EGBERTKELBERT2012
            case ('EGBERTKELBERT2012_MODEM'); conv = EGBERTKELBERT2012_MODEM
            case ('LWS');         conv = LWS
            ! NOTE: 'NATIVE' is NOT resolved here -- it is solver-specific (S1
            ! and S2 have different native conventions), so FWD1D.f90 resolves
            ! OUTPUT_CONVENTION='NATIVE' to native_convention(SOLVER) directly.
            case default
                write(0,*) 'output_convention: unknown convention name: ', trim(name)
                stop 1
        end select
    end function get_convention

    ! Native (raw, un-transformed) convention actually produced by each solver.
    ! This is the SINGLE SOURCE OF TRUTH for "native": every convention transform
    ! (rescale_source_coeffs, rescale_source_radial, apply_output_convention)
    ! converts native_convention(solver) -> target, and OUTPUT_CONVENTION='NATIVE'
    ! (resolved in FWD1D.f90) simply sets target = native_convention(SOLVER),
    ! which makes ALL of those transforms mechanical no-ops -> exact raw output.
    !
    ! S1 and S2 share the same native time/norm/CS/theta/r/primary bookkeeping
    ! (the SUNEGBERT2012 preset's values -- MINUS_IWT, FULLY_NORM, CS-included,
    ! colatitude N->S, r up, E-primary), and differ natively in exactly ONE
    ! dimension: the RADIAL source-amplitude normalization (S1=SURFACE,
    ! S2=MULTIPOLE, see native_radial). Because SUNEGBERT2012 itself now carries
    ! RADIAL_MULTIPOLE, native_convention('S2') == SUNEGBERT2012 exactly (S2 +
    ! SUNEGBERT2012 is the identity); native_convention('S1') is the same but
    ! with radial_norm = SURFACE.
    function native_convention(solver) result(conv)
        character(*), intent(in)   :: solver
        type(output_convention_t)  :: conv
        select case (trim(solver))
            case ('S1', 'S2')
                conv = SUNEGBERT2012                     ! shared native bookkeeping
                conv%radial_norm = native_radial(solver) ! the one solver-specific dimension
            case default
                write(0,*) 'output_convention: unknown solver: ', trim(solver)
                stop 1
        end select
    end function native_convention

    ! Native per-degree RADIAL (source-amplitude) normalization of each solver.
    ! S1 (field1d.f90) assembles Hr with an R0^2/r^2 factor, so its external
    ! radial field at the surface is exactly Y_l^m for coeff=1 -> RADIAL_SURFACE.
    ! S2 (field1d_s2.f90) uses Sun & Egbert eq (6) (l(l+1)/r^2) with the
    ! dimensionless (r/a)^(l+1)=1 potential -> RADIAL_MULTIPOLE.
    function native_radial(solver) result(rn)
        character(*), intent(in) :: solver
        character(len=12)        :: rn
        select case (trim(solver))
            case ('S1'); rn = RADIAL_SURFACE
            case ('S2'); rn = RADIAL_MULTIPOLE
            case default
                write(0,*) 'output_convention: unknown solver: ', trim(solver)
                stop 1
        end select
    end function native_radial

    !===========================================================================
    ! (a) Source-coefficient rescale -- applied BEFORE the solve.
    !===========================================================================

    subroutine rescale_source_coeffs(coeff, lmax, target_conv)
        ! Rescales coeff(:) (in the solvers' own m=0,1,-1,2,-2,... flat
        ! ordering) so that, once combined with the solver's FIXED internal
        ! basis functions (vsharm's Y_l^m -- independently verified
        ! to be fully-normalized AND Condon-Shortley-included,
        ! i.e. exactly the SUNEGBERT2012 preset's norm/CS values), the resulting
        ! field represents what target_conv's own norm/CS basis would give.
        !
        ! Concretely: for basis functions Y_target(l,m) = ratio(l,m) *
        ! Y_vsharm(l,m), feeding rescaled_coeff = coeff * ratio into the
        ! (unchanged) solver gives rescaled_coeff * Y_vsharm = coeff *
        ! Y_target -- the same physical field as if Y_target had been used
        ! directly. ratio(l,m) is therefore computed relative to vsharm's OWN
        ! fixed convention (norm=FULL, CS=true), NOT relative to some other
        ! nominal "source file" convention -- this is what makes
        ! target_conv=SUNEGBERT2012 (matching vsharm exactly) a
        ! mechanically verifiable no-op (ratio=1 for every l,m), independent
        ! of the harder-to-externally-verify Schmidt-normalization constant
        ! below (see apply_norm_cs).
        !
        ! CAVEAT (flagged per the spec's own warning about factor-of-sqrt(2)
        ! slips): the Schmidt-vs-fully-normalized ratio's absolute value has
        ! NOT been independently cross-checked against an external reference
        ! (neither Sun & Egbert (2012) Section 2 nor pythonSolver/
        ! spherical_em_induction.py define a Schmidt-seminormalized basis at
        ! all -- both use the fully-normalized convention only, which IS
        ! independently verified). The Condon-Shortley (-1)^m factor and the
        ! overall STRUCTURE of the transform are on solid footing (see
        ! apply_norm_cs); only the exact Schmidt/full normalization constant
        ! is a documented, flagged assumption -- see testing/test_norm_condon_shortley.f90
        ! for the self-consistency checks that ARE possible without an
        ! external reference (round-trip identity, m=0 special case).
        complex(prec), dimension(:), intent(inout) :: coeff  ! in m=0,1,-1,2,-2,... order
        integer, intent(in)                        :: lmax
        type(output_convention_t), intent(in)      :: target_conv
        integer :: l, m, idx

        idx = 1                                        ! coeff(1) is l=0 (unused)
        do l = 1, lmax
            idx = idx + 1
            call apply_norm_cs(coeff(idx), l, 0, target_conv)      ! m=0 term
            do m = 1, l
                call apply_norm_cs(coeff(idx+1), l,  m, target_conv)   ! +m
                call apply_norm_cs(coeff(idx+2), l, -m, target_conv)   ! -m
                idx = idx + 2
            end do
        end do

    end subroutine rescale_source_coeffs

    subroutine apply_norm_cs(c, l, m, target_conv)
        ! Rescales one coefficient c (degree l, order m) by ratio(l,m) =
        ! [target_conv's basis] / [vsharm's fixed FULL/CS-included basis].
        complex(prec), intent(inout)          :: c
        integer, intent(in)                   :: l, m
        type(output_convention_t), intent(in) :: target_conv
        real(prec) :: delta_m0

        delta_m0 = 0.0_prec
        if (m == 0) delta_m0 = 1.0_prec

        ! Condon-Shortley: vsharm is fixed CS-included (.true.). Flip by
        ! (-1)^m iff target_conv wants CS excluded -- see the derivation in
        ! rescale_source_coeffs's header comment (this is NOT conditioned on
        ! any nominal "source file convention", only on vsharm's own fixed,
        ! verified behavior).
        if (.not. target_conv%condon_shortley) then
            c = c * ((-1.0_prec)**m)
        end if

        ! Harmonic normalization: vsharm is fixed FULLY-NORMALIZED. Rescale
        ! by the Schmidt/full ratio iff target_conv wants Schmidt semi-norm.
        ! Ratio direction: Schmidt seminorm S_l^m = sqrt((2-delta_m0)*(l-m)!/(l+m)!) * P_l^m
        ! (no (2l+1)/(4*pi) growth factor); fully-normalized Y_l^m = sqrt((2l+1)/(4*pi)
        ! * (l-m)!/(l+m)!) * P_l^m * exp(i*m*phi) (vsharm's own, independently
        ! verified convention) -- so Y_schmidt/Y_full = sqrt(4*pi*(2-delta_m0)/(2l+1)),
        ! the reciprocal of the ratio literally written in the implementation
        ! spec (flagged there as needing independent verification, which was
        ! not possible in this offline session -- see the CAVEAT above).
        if (target_conv%harmonic_norm == NORM_SCHMIDT) then
            c = c * sqrt(4.0_prec * PI * (2.0_prec - delta_m0) / real(2*l+1, prec))
        end if

    end subroutine apply_norm_cs

    !===========================================================================
    ! (a2) Radial (source-amplitude) renormalization -- applied BEFORE the solve,
    !      right after rescale_source_coeffs. Converts the source coeff from the
    !      SOLVER's native radial normalization (native_radial(solver)) to the
    !      TARGET radial normalization requested by the OUTPUT_CONVENTION
    !      (target_radial), per DEGREE l. Because BOTH solvers are then expressed
    !      in the SAME target radial scale, S1 and S2 produce IDENTICAL output
    !      (up to each solver's independent numerical accuracy) for the same
    !      coeff and the same OUTPUT_CONVENTION.
    !
    !      Why this is needed and is NOT covered by rescale_source_coeffs above:
    !      that routine reconciles the ANGULAR (spherical-harmonic norm + CS
    !      phase) convention only. The RADIAL / source-amplitude normalization --
    !      what "coeff(l,m)=1" means for the radial part of the field -- is a
    !      genuinely separate convention dimension, and the ONE dimension where
    !      S1 and S2 differ natively (S1=SURFACE, S2=MULTIPOLE, see
    !      native_radial). The per-degree amplitudes of the named conventions are
    !      defined in radial_amplitude(); the rescale factor is simply
    !      A(target,l) / A(native(solver),l). Applied identically to E and H, so
    !      Z = E/H is unchanged.
    !===========================================================================
    subroutine rescale_source_radial(coeff, lmax, solver, target_radial, r0)
        complex(prec), dimension(:), intent(inout) :: coeff  ! m=0,1,-1,2,-2,... order
        integer, intent(in)                        :: lmax
        character(len=*), intent(in)               :: solver        ! 'S1' | 'S2'
        character(len=*), intent(in)               :: target_radial ! RADIAL_* token
        real(prec), intent(in)                     :: r0            ! Earth radius, metres
        integer           :: l, m, idx
        real(prec)        :: fac
        character(len=12) :: native

        native = native_radial(solver)
        idx = 1                                        ! coeff(1) is l=0 (unused)
        do l = 1, lmax
            fac = radial_amplitude(target_radial, l, r0) / radial_amplitude(native, l, r0)
            idx = idx + 1
            coeff(idx) = coeff(idx) * fac                    ! m=0
            do m = 1, l
                coeff(idx+1) = coeff(idx+1) * fac            ! +m
                coeff(idx+2) = coeff(idx+2) * fac            ! -m
                idx = idx + 2
            end do
        end do
    end subroutine rescale_source_radial

    real(prec) function radial_amplitude(radial_norm, l, r0) result(A)
        ! Per-degree amplitude A(l) of a named RADIAL normalization, relative to
        ! the MULTIPOLE baseline (Sun & Egbert eq (6) with the dimensionless
        ! (r/a)^(l+1)=1 potential, A=1). A field expressed in convention X equals
        ! (the same field in MULTIPOLE) * A(X,l). THIS is the single place to
        ! adjust or add a per-degree radial normalization: add a case with your
        ! own A(l) formula.
        !   RADIAL_MULTIPOLE   : A = 1               -- S2 native; coeff of (r/a)^(l+1)=1
        !   RADIAL_SURFACE     : A = R0^2/(l(l+1))   -- S1 native; unit surface radial field
        !   RADIAL_DIMENSIONAL : A = R0^(l+1)        -- Sun & Egbert text, coeff of r^(l+1)=1
        !                                               (r in metres) -- OVERFLOWS for l >~ 44
        !   RADIAL_POTENTIAL   : A = R0^2/(l+1)      -- external scalar potential V
        !                                               (Gauss coefficient) parameterization.
        !     Derivation: in the insulating atmosphere (sigma=0, curl H=0), H=-grad(V)
        !     is valid. Comparing that identity to the solvers' own H formulas
        !     (H_theta=(1/r)*T_l'(r)*Yt, H_r=(l(l+1)/r^2)*T_l(r)*Y -- see
        !     field1d_s2.f90's header) gives V_l(r) = -T_l'(r) (and T_l''=
        !     l(l+1)/r^2*T_l confirms both solve the same source-free radial ODE).
        !     The classical external/inducing potential expansion is
        !     V = R0*sum_l (r/R0)^l * epsilon_l^m * Y_l^m (the regular-at-origin,
        !     r-growing term -- the standard Gauss-coefficient convention used in
        !     geomagnetism and ionospheric-current/Sq source modeling). Substituting
        !     T_l^ext(r) = alpha_l^T*(r/R0)^(l+1) (MULTIPOLE's own source term) and
        !     differentiating: V_l^ext(r) = -alpha_l^T*(l+1)/R0*(r/R0)^l, so matching
        !     powers gives epsilon_l^m proportional to alpha_l^T/(l+1), i.e.
        !     A(l) = R0^2/(l+1) -- the RADIAL_SURFACE formula with the extra 1/l
        !     factor removed (as for RADIAL_SURFACE, the sign is absorbed into the
        !     convention -- A(l) is a plain positive real scale, matching every
        !     other radial_norm case here).
        character(len=*), intent(in) :: radial_norm
        integer, intent(in)          :: l
        real(prec), intent(in)       :: r0
        select case (trim(radial_norm))
            case (RADIAL_MULTIPOLE)
                A = 1.0_prec
            case (RADIAL_SURFACE)
                A = r0*r0 / real(l*(l+1), prec)
            case (RADIAL_DIMENSIONAL)
                if (l+1 > 44) then
                    write(0,*) 'output_convention: RADIAL_DIMENSIONAL overflows double precision at l=',l, &
                               ' (R0^(l+1) > ~1e308). Use RADIAL_SURFACE or RADIAL_MULTIPOLE for high degree.'
                    stop 1
                end if
                A = r0**(l+1)
            case (RADIAL_POTENTIAL)
                A = r0*r0 / real(l+1, prec)
            case default
                write(0,*) 'output_convention: unknown radial_norm: ', trim(radial_norm)
                stop 1
        end select
    end function radial_amplitude

    !===========================================================================
    ! (b) Output-array transform -- applied AFTER the solve.
    !===========================================================================

    subroutine apply_output_convention(H, E, native, target_conv)
        type(cvector), intent(inout)          :: H, E
        type(output_convention_t), intent(in) :: native, target_conv

        ! 1) time convention: conjugate everything iff they differ
        if (native%time_convention /= target_conv%time_convention) then
            H%x = conjg(H%x); H%y = conjg(H%y); H%z = conjg(H%z)
            E%x = conjg(E%x); E%y = conjg(E%y); E%z = conjg(E%z)
        end if

        ! 4) theta: colatitude(N->S) <-> latitude(S->N) is a relabeling of the
        ! SAME physical field under theta -> pi-theta: reverse the theta (j,
        ! 2nd) index of every component array (each component's own extent,
        ! see reverse_theta_index), AND negate the theta component
        ! (Hy=H_theta, Ey=E_theta) to account for theta-hat's direction
        ! flipping under the relabeling. Index reversal and component
        ! negation are a matched pair -- see reverse_theta_index's header.
        if (native%theta_convention /= target_conv%theta_convention) then
            call reverse_theta_index(H); call reverse_theta_index(E)
            H%y = -H%y ;  E%y = -E%y
        end if

        ! 5) r: down(depth) <-> up(radius) is a reversal of the r (k, 3rd)
        ! index of every component array, AND negation of the r component
        ! (Hz=H_r, Ez=E_r).
        if (native%r_convention /= target_conv%r_convention) then
            call reverse_r_index(H); call reverse_r_index(E)
            H%z = -H%z ;  E%z = -E%z
        end if

    end subroutine apply_output_convention

    subroutine reverse_theta_index(V)
        ! Reverses the theta (2nd) index of each of V%x,V%y,V%z IN PLACE,
        ! using EACH component's OWN extent along that axis (NOT a shared N)
        ! -- on EDGE, V%x/V%z have ny+1 node-theta slots while V%y has ny
        ! mid-theta slots; on FACE the roles swap (see sg_vector.f90's
        ! create_cvector). Reversing "size(V%comp,2):1:-1" per component
        ! handles both node- and mid-theta arrays correctly without needing
        ! to know which is which -- reversing the index order is exactly
        ! what turns a "colatitude, North->South" traversal into a
        ! "latitude, South->North" one (or the reverse), for both mid-cell
        ! and node arrays alike; it is the reversal itself (not any special
        ! handling of the +1 staggered extent) that makes this correct.
        type(cvector), intent(inout) :: V
        V%x = V%x(:, size(V%x,2):1:-1, :)
        V%y = V%y(:, size(V%y,2):1:-1, :)
        V%z = V%z(:, size(V%z,2):1:-1, :)
    end subroutine reverse_theta_index

    subroutine reverse_r_index(V)
        ! Reverses the r (3rd) index of each of V%x,V%y,V%z IN PLACE, same
        ! per-component-extent handling as reverse_theta_index.
        type(cvector), intent(inout) :: V
        V%x = V%x(:, :, size(V%x,3):1:-1)
        V%y = V%y(:, :, size(V%y,3):1:-1)
        V%z = V%z(:, :, size(V%z,3):1:-1)
    end subroutine reverse_r_index

    !===========================================================================
    ! (c) ModEM source-normalization rescale -- applied AFTER apply_output_-
    !     convention, only for target conventions with modem_normalize=.true.
    !===========================================================================

    subroutine apply_modem_normalization(H, E, omega, d_air)
        ! Rescales BOTH H and E (every component) by the single complex factor
        !     1 / (i * omega * mu0 * G),   G = (3/2) * sqrt(3/(4*pi)) * d_air,
        ! converting global1d's toroidal-potential ("unit external multipole")
        ! source normalization into ModEM's boundary-E-field ("unit surface E")
        ! normalization. The two differ by exactly c = i*omega*mu0*G -- the
        ! potential-vs-field factor of Faraday's law (see the header comment on
        ! the EGBERTKELBERT2012_MODEM preset and, for the full derivation,
        ! docs/source_normalization.md / .pdf).
        !
        ! G is derived from the AIR-LAYER GEOMETRY: d_air is the air-column
        ! thickness (top-of-grid radius minus Earth radius, in metres); the
        ! (3/2)*sqrt(3/(4*pi)) factor is |H_theta(surface)| per unit source for
        ! the l=1 (P10/MT) mode at the equator in the deep-induction limit
        ! (Q_1 -> 1/2): (2 - Q_1) = 3/2 combines the uniform-external (2) and
        ! induced (-Q_1) parts, sqrt(3/(4*pi)) is the fully-normalized Y_1^0
        ! amplitude. This is EXACT for l=1 at the equator in that limit; a
        ! few-% amplitude / <~5deg phase residual remains (the earth's finite-Q
        ! C-response and ModEM's flat-vs-spherical BC), both genuine model
        ! differences. Applying the same factor to E and H preserves Z=E/H.
        !
        ! d_air must be passed by the caller (from the grid: grid%r(1)*1e3 minus
        ! Earth radius r0), and omega=2*pi/period for the current period, since
        ! this module does not own the grid or the period loop.
        type(cvector), intent(inout) :: H, E
        real(prec), intent(in)       :: omega, d_air
        real(prec)                   :: G
        complex(prec)                :: factor

        G = 1.5_prec * sqrt(3.0_prec/(4.0_prec*PI)) * d_air
        factor = 1.0_prec / dcmplx(0.0_prec, omega*MU_0*G)   ! 1/(i*omega*mu0*G)

        H%x = factor*H%x ; H%y = factor*H%y ; H%z = factor*H%z
        E%x = factor*E%x ; E%y = factor*E%y ; E%z = factor*E%z
    end subroutine apply_modem_normalization

end module output_convention
