program test_norm_condon_shortley
! ****************************************************************************
! Unit test for output_convention.f90's rescale_source_coeffs (and its
! per-coefficient helper apply_norm_cs), per FWD1D_output_convention_
! implementation_spec.md section 2.5.
!
! Checks hand-computed norm/Condon-Shortley ratios against what the code
! actually produces, for several (l,m) including the m=0 special case, for
! all three named target conventions:
!
!   - SUNEGBERT2012 (norm=FULL, CS=true): matches vsharm's own fixed
!     basis EXACTLY (see output_convention.f90's header note) -- this is a
!     MECHANICALLY VERIFIABLE no-op, ratio=1 for every (l,m), independent of
!     any external Schmidt-normalization reference.
!   - EGBERTKELBERT2012 (norm=FULL, CS=true as of 2026-07-26 -- see
!     output_convention.f90's own note on this preset): also matches vsharm's
!     basis exactly, so ALSO a no-op (ratio=1 for every l,m), same as SUNEGBERT2012 --
!     this was NOT the case before the CS=true fix (it used to expect
!     ratio=(-1)^m); this test's own expectation must track that fix, not
!     the original spec's since-corrected assumption.
!   - KELBERT2006 (norm=SCHMIDT, CS=false): ratio = (-1)^m * sqrt(4*pi*(2-
!     delta_m0)/(2l+1)) -- exercises the Schmidt-norm formula this project
!     chose to implement (see output_convention.f90's CAVEAT comment: this
!     formula's absolute value has NOT been independently cross-checked
!     against an external reference, since neither Sun & Egbert (2012) nor
!     pythonSolver define a Schmidt-seminormalized basis at all -- this test
!     confirms the CODE matches its OWN documented formula exactly, not that
!     the formula itself matches some external Schmidt-norm literature
!     definition).
! ****************************************************************************
    use output_convention
    use math_constants, only: prec, PI
    implicit none

    integer, parameter :: lmax = 2
    complex(prec) :: coeff(9)   ! (lmax+1)^2 = 9, ordering: l=0(unused),l=1(m=0,1,-1),l=2(m=0,1,-1,2,-2)
    complex(prec) :: orig(9)
    real(prec) :: expect_re, tol
    integer :: l, m, idx
    logical :: all_pass

    tol = 1.0e-10_prec
    all_pass = .true.

    write(*,*) '=== test_norm_condon_shortley: rescale_source_coeffs vs hand-computed ratios ==='
    write(*,*)

    ! ---- SUNEGBERT2012: expect ratio=1 for every (l,m) (mechanically verifiable no-op) ----
    call fill_coeff(orig)
    coeff = orig
    call rescale_source_coeffs(coeff, lmax, SUNEGBERT2012)
    call check_block('SUNEGBERT2012', orig, coeff, target_is_sunegbert2012=.true.)

    ! ---- EGBERTKELBERT2012 (CS=true, norm=FULL as of 2026-07-26): expect ratio=1, same as SUNEGBERT2012 ----
    call fill_coeff(orig)
    coeff = orig
    call rescale_source_coeffs(coeff, lmax, EGBERTKELBERT2012)
    call check_block('EGBERTKELBERT2012', orig, coeff, target_is_sunegbert2012=.true.)

    ! ---- KELBERT2006 (CS=false, norm=SCHMIDT): expect ratio=(-1)^m*sqrt(4pi(2-d_m0)/(2l+1)) ----
    call fill_coeff(orig)
    coeff = orig
    call rescale_source_coeffs(coeff, lmax, KELBERT2006)
    call check_block('KELBERT2006', orig, coeff, target_is_sunegbert2012=.false., norm_full=.false.)

    write(*,*)
    if (all_pass) then
        write(*,*) 'ALL CHECKS PASSED'
    else
        write(*,*) 'AT LEAST ONE CHECK FAILED -- see above'
        stop 1
    end if

contains

    subroutine fill_coeff(c)
        complex(prec), intent(out) :: c(9)
        c = dcmplx(1.0d0, 0.0d0)   ! unit coefficient everywhere, easiest to hand-check ratios directly
    end subroutine fill_coeff

    ! index layout: idx=1 (l=0,unused), idx=2 (l=1,m=0), idx=3(l=1,m=+1), idx=4(l=1,m=-1),
    ! idx=5(l=2,m=0), idx=6(l=2,m=+1), idx=7(l=2,m=-1), idx=8(l=2,m=+2), idx=9(l=2,m=-2)
    subroutine check_block(label, orig_c, new_c, target_is_sunegbert2012, norm_full)
        character(*), intent(in)  :: label
        complex(prec), intent(in) :: orig_c(9), new_c(9)
        logical, intent(in)       :: target_is_sunegbert2012
        logical, intent(in), optional :: norm_full

        write(*,'(a,a)') 'Target: ', label
        call check_one(2, 1, 0, orig_c, new_c, target_is_sunegbert2012, norm_full)
        call check_one(3, 1, 1, orig_c, new_c, target_is_sunegbert2012, norm_full)
        call check_one(4, 1,-1, orig_c, new_c, target_is_sunegbert2012, norm_full)
        call check_one(5, 2, 0, orig_c, new_c, target_is_sunegbert2012, norm_full)
        call check_one(6, 2, 1, orig_c, new_c, target_is_sunegbert2012, norm_full)
        call check_one(7, 2,-1, orig_c, new_c, target_is_sunegbert2012, norm_full)
        call check_one(8, 2, 2, orig_c, new_c, target_is_sunegbert2012, norm_full)
        call check_one(9, 2,-2, orig_c, new_c, target_is_sunegbert2012, norm_full)
        write(*,*)
    end subroutine check_block

    subroutine check_one(idx, l, m, orig_c, new_c, target_is_sunegbert2012, norm_full)
        integer, intent(in) :: idx, l, m
        complex(prec), intent(in) :: orig_c(9), new_c(9)
        logical, intent(in) :: target_is_sunegbert2012
        logical, intent(in), optional :: norm_full
        real(prec) :: expect_ratio, actual_ratio, delta_m0, diff

        if (target_is_sunegbert2012) then
            expect_ratio = 1.0_prec
        else
            delta_m0 = 0.0_prec
            if (m == 0) delta_m0 = 1.0_prec
            expect_ratio = (-1.0_prec)**m
            if (.not. norm_full) then
                expect_ratio = expect_ratio * sqrt(4.0_prec*PI*(2.0_prec-delta_m0)/real(2*l+1,prec))
            end if
        end if

        actual_ratio = dble(new_c(idx)) / dble(orig_c(idx))
        diff = abs(actual_ratio - expect_ratio)

        if (diff < tol) then
            write(*,'(a,i2,a,i2,a,f14.8,a,f14.8,a)') '  l=',l,' m=',m,':  ratio=',actual_ratio, &
                '  expect=',expect_ratio,'  PASS'
        else
            write(*,'(a,i2,a,i2,a,f14.8,a,f14.8,a)') '  l=',l,' m=',m,':  ratio=',actual_ratio, &
                '  expect=',expect_ratio,'  FAIL'
            all_pass = .false.
        end if

    end subroutine check_one

end program test_norm_condon_shortley
