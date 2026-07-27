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

    use math_constants, only: prec, PI
    use sg_vector, only: cvector

    implicit none

    private
    public :: output_convention_t
    public :: TIME_POSITIVE, TIME_NEGATIVE, NORM_SCHMIDT, NORM_FULL, &
              THETA_COLAT, THETA_LAT, R_DOWN, R_UP
    public :: KELBERT2006, SUNEGBERT2012, EGBERTKELBERT2012
    public :: get_convention, native_convention
    public :: rescale_source_coeffs, apply_output_convention

    ! ------- enumerated tokens (character params keep the .prm/config human-readable) -------
    character(len=10), parameter :: TIME_POSITIVE = 'PLUS_IWT'
    character(len=10), parameter :: TIME_NEGATIVE = 'MINUS_IWT'
    character(len=16), parameter :: NORM_SCHMIDT   = 'SCHMIDT_SEMINORM'
    character(len=16), parameter :: NORM_FULL      = 'FULLY_NORM'
    character(len=10), parameter :: THETA_COLAT = 'COLAT_N2S'  ! colatitude, North->South
    character(len=10), parameter :: THETA_LAT   = 'LAT_S2N'    ! latitude,   South->North
    character(len=6),  parameter :: R_DOWN = 'R_DOWN'          ! r increases downward (depth)
    character(len=6),  parameter :: R_UP   = 'R_UP'            ! r increases upward (radius)

    ! ------- the convention type: five field-bookkeeping dimensions + primary_grid -------
    type :: output_convention_t
        character(len=24) :: name              ! 'KELBERT2006' etc, for logging
        character(len=10) :: time_convention    ! TIME_POSITIVE | TIME_NEGATIVE
        character(len=16) :: harmonic_norm      ! NORM_SCHMIDT | NORM_FULL
        logical            :: condon_shortley   ! .true. = Condon-Shortley phase included
        character(len=10) :: theta_convention   ! THETA_COLAT | THETA_LAT
        character(len=6)  :: r_convention       ! R_DOWN | R_UP
        character(len=1)  :: primary_grid       ! 'H' | 'E'  (records intended staggering)
    end type output_convention_t

    ! ------- the three predefined conventions (kept exactly as specified --
    ! see the NATIVE CONVENTIONS note above for why these are fixed, named
    ! targets, not necessarily any solver's raw native output) -------
    type(output_convention_t), parameter :: KELBERT2006 = output_convention_t( &
        'KELBERT2006', TIME_POSITIVE, NORM_SCHMIDT, .false., THETA_COLAT, R_DOWN, 'H')

    type(output_convention_t), parameter :: SUNEGBERT2012 = output_convention_t( &
        'SUNEGBERT2012', TIME_NEGATIVE, NORM_FULL, .true., THETA_COLAT, R_UP, 'E')

    ! condon_shortley=.true. here (2026-07-26, per direct user instruction, for
    ! validation purposes): confirmed empirically (test_diagnose_egbertkelbert.f90)
    ! that ModEM's own .prm source-coefficient files are ALREADY written in
    ! vsharm's native CS-included basis (they reproduce ModEM correctly when fed
    ! UNSCALED into either solver, see test_vs_modem_1D.f90) -- .false. here
    ! (the original spec's assumption) spuriously flips every odd-m source
    ! (e.g. Mode2, l=1 m=+-1) via rescale_source_coeffs's (-1)^m term while
    ! leaving m=0 sources (e.g. Mode1) untouched, producing a mode-dependent
    ! sign inconsistency. NOTE (2026-07-26, per user): the historical ModEM
    ! solver this convention is named for never actually used spherical
    ! harmonics at all; a future CUSTOM/named preset with condon_shortley=
    ! .false. may still be needed for that use case -- track separately,
    ! do not silently revert this field back to .false.
    type(output_convention_t), parameter :: EGBERTKELBERT2012 = output_convention_t( &
        'EGBERTKELBERT2012', TIME_NEGATIVE, NORM_FULL, .true., THETA_LAT, R_DOWN, 'E')

contains

    ! resolve a name string -> convention (for a config-file / CUSTOM interface later)
    function get_convention(name) result(conv)
        character(*), intent(in)   :: name
        type(output_convention_t)  :: conv
        select case (trim(name))
            case ('KELBERT2006');        conv = KELBERT2006
            case ('SUNEGBERT2012');      conv = SUNEGBERT2012
            case ('EGBERTKELBERT2012');  conv = EGBERTKELBERT2012
            case default
                write(0,*) 'output_convention: unknown convention name: ', trim(name)
                stop 1
        end select
    end function get_convention

    ! Native (raw, un-transformed) convention actually produced by each
    ! solver, empirically verified -- see the NATIVE CONVENTIONS note above.
    ! As of 2026-07-27, S1 and S2 have IDENTICAL native conventions in every
    ! dimension (including primary_grid, now 'E' for both) -- both cases
    ! simply return the SUNEGBERT2012 preset (S2's native convention). Kept
    ! as a solver-name lookup (rather than having callers use the
    ! SUNEGBERT2012 preset directly) so
    ! FWD1D.f90 stays agnostic to the fact that the two happen to coincide
    ! today -- if a future solver's native output genuinely differs, only
    ! this function needs a new case, not every call site.
    function native_convention(solver) result(conv)
        character(*), intent(in)   :: solver
        type(output_convention_t)  :: conv
        select case (trim(solver))
            case ('S1', 'S2')
                conv = SUNEGBERT2012
            case default
                write(0,*) 'output_convention: unknown solver: ', trim(solver)
                stop 1
        end select
    end function native_convention

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

end module output_convention
