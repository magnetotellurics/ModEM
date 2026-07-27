program test_reversals
! ****************************************************************************
! Unit test for output_convention.f90's apply_output_convention theta/r
! index-reversal + component-negation logic, per FWD1D_output_convention_
! implementation_spec.md section 2.5.
!
! Two checks, both off-equator/off-center (not just at a point where a
! reversal might look like a no-op):
!
! (1) PURE INDEX-REVERSAL CORRECTNESS. Fills every component of an EDGE H /
!     FACE E cvector pair with a synthetic value encoding its OWN (i,j,k)
!     array position (val = i*10000+j*100+k), separately for EACH of the six
!     arrays (H%x,H%y,H%z,E%x,E%y,E%z) -- these have DIFFERENT extents along
!     the theta/r axes on EDGE vs FACE (see sg_vector.f90's create_cvector).
!     After calling apply_output_convention with native/target conventions
!     that differ ONLY in theta_convention (time/r/norm/CS held fixed),
!     asserts new(i,j,k) == -old(i, N_theta_component+1-j, k) for EACH
!     component's OWN N_theta (its own size along dim 2) -- this is exactly
!     the "off-by-one silently shifts the field by half a cell" trap the
!     spec warns about: using the WRONG N (e.g. a shared N instead of each
!     component's own) would fail this check for at least one component.
!     Repeated for r_convention (dim 3).
!
! (2) PHYSICAL SANITY CHECK. Fills H%y (H_theta) with the analytic pattern
!     H_theta = cos(theta), evaluated at each node's OWN grid%th(j) (mid-cell
!     average for EDGE H%y, which is mid-theta -- see field1d.f90's own
!     staggering table). After the theta reversal+negation, checks that the
!     new value at an OFF-EQUATOR point equals the field "re-evaluated under
!     theta -> pi-theta", i.e. -cos(theta_old_at_the_mirrored_node) --
!     confirming the negation direction (not just index bookkeeping) is
!     correct, per the spec's explicit warning that a pure equatorial check
!     cannot distinguish "correct reversal+negate" from "a lucky point".
! ****************************************************************************
    use output_convention
    use griddef
    use sg_vector
    use math_constants, only: prec, PI
    implicit none

    type(grid_t)   :: grid
    type(cvector)  :: h1d, e1d, h1d_orig, e1d_orig
    integer, parameter :: nx = 4, ny = 6, nz = 3
    integer :: ii, jj, kk
    type(output_convention_t) :: native, target_theta, target_r
    logical :: all_pass

    all_pass = .true.

    write(*,*) '=== test_reversals: apply_output_convention theta/r index reversal ==='
    write(*,*)

    ! ---- build a simple, non-degenerate global-ish grid (not symmetric in
    ! nx/ny/nz, so a transposed-axis bug would not accidentally pass) ----
    call create_grid(nx, ny, nz, grid)
    grid%nx = nx; grid%ny = ny; grid%nz = nz
    grid%nzAir = nz; grid%nzCrust = 0; grid%nzEarth = 0
    do ii = 1, nx+1
        grid%ph(ii) = (ii-1)*(360.0d0/nx) * d2r
    end do
    do ii = 1, nx
        grid%dp(ii) = grid%ph(ii+1) - grid%ph(ii)
    end do
    do jj = 1, ny+1
        grid%th(jj) = (5.0d0 + (jj-1)*(170.0d0/ny)) * d2r   ! avoid exact poles
    end do
    do jj = 1, ny
        grid%dt(jj) = grid%th(jj+1) - grid%th(jj)
    end do
    do kk = 1, nz+1
        grid%r(kk) = 6371.0d0 - (kk-1)*10.0d0
    end do
    do kk = 1, nz
        grid%dr(kk) = grid%r(kk) - grid%r(kk+1)
    end do
    grid%allocated = .true.

    call create_cvector(grid, h1d, EDGE)
    call create_cvector(grid, e1d, FACE)
    call create_cvector(grid, h1d_orig, EDGE)
    call create_cvector(grid, e1d_orig, FACE)

    ! ---- fill with synthetic index-encoding pattern ----
    call fill_index_pattern(h1d)
    call fill_index_pattern(e1d)
    h1d_orig = h1d
    e1d_orig = e1d

    ! ---- (1a) theta reversal check ----
    native       = output_convention_t('NATIVE', TIME_NEGATIVE, NORM_FULL, .true., THETA_COLAT, R_UP, 'H')
    target_theta = output_convention_t('TARGET', TIME_NEGATIVE, NORM_FULL, .true., THETA_LAT,   R_UP, 'H')
    call apply_output_convention(h1d, e1d, native, target_theta)
    write(*,*) '--- (1a) theta index-reversal, synthetic pattern ---'
    ! Only the theta-DIRECTED component (%y, since y=theta in this codebase's
    ! x=phi/y=theta/z=r convention) gets NEGATED -- phi-hat/r-hat do not flip
    ! under a pure theta relabeling, so %x/%z are index-reversed only.
    call check_theta_reversal('H%x', h1d%x, h1d_orig%x, negate=.false.)
    call check_theta_reversal('H%y', h1d%y, h1d_orig%y, negate=.true.)
    call check_theta_reversal('H%z', h1d%z, h1d_orig%z, negate=.false.)
    call check_theta_reversal('E%x', e1d%x, e1d_orig%x, negate=.false.)
    call check_theta_reversal('E%y', e1d%y, e1d_orig%y, negate=.true.)
    call check_theta_reversal('E%z', e1d%z, e1d_orig%z, negate=.false.)
    write(*,*)

    ! ---- (1b) r reversal check (reset to fresh pattern first) ----
    call fill_index_pattern(h1d)
    call fill_index_pattern(e1d)
    h1d_orig = h1d
    e1d_orig = e1d
    target_r = output_convention_t('TARGET', TIME_NEGATIVE, NORM_FULL, .true., THETA_COLAT, R_DOWN, 'H')
    call apply_output_convention(h1d, e1d, native, target_r)
    write(*,*) '--- (1b) r index-reversal, synthetic pattern ---'
    ! Only the r-DIRECTED component (%z) gets NEGATED; phi-hat/theta-hat do
    ! not flip under a pure r (up/down) relabeling.
    call check_r_reversal('H%x', h1d%x, h1d_orig%x, negate=.false.)
    call check_r_reversal('H%y', h1d%y, h1d_orig%y, negate=.false.)
    call check_r_reversal('H%z', h1d%z, h1d_orig%z, negate=.true.)
    call check_r_reversal('E%x', e1d%x, e1d_orig%x, negate=.false.)
    call check_r_reversal('E%y', e1d%y, e1d_orig%y, negate=.false.)
    call check_r_reversal('E%z', e1d%z, e1d_orig%z, negate=.true.)
    write(*,*)

    ! ---- (2) physical sanity check: H_theta = cos(theta), off-equator ----
    call fill_costheta_pattern(h1d, grid)
    h1d_orig = h1d
    call apply_output_convention(h1d, e1d, native, target_theta)
    write(*,*) '--- (2) physical check: H%y = cos(theta), off-equator, post reversal+negate ---'
    call check_costheta(h1d%y, h1d_orig%y, grid)

    write(*,*)
    if (all_pass) then
        write(*,*) 'ALL CHECKS PASSED'
    else
        write(*,*) 'AT LEAST ONE CHECK FAILED -- see above'
        stop 1
    end if

contains

    subroutine fill_index_pattern(V)
        type(cvector), intent(inout) :: V
        integer :: i,j,k
        do k = 1,size(V%x,3); do j = 1,size(V%x,2); do i = 1,size(V%x,1)
            V%x(i,j,k) = dcmplx(dble(i*10000+j*100+k), 0.0d0)
        end do; end do; end do
        do k = 1,size(V%y,3); do j = 1,size(V%y,2); do i = 1,size(V%y,1)
            V%y(i,j,k) = dcmplx(dble(i*10000+j*100+k), 0.0d0)
        end do; end do; end do
        do k = 1,size(V%z,3); do j = 1,size(V%z,2); do i = 1,size(V%z,1)
            V%z(i,j,k) = dcmplx(dble(i*10000+j*100+k), 0.0d0)
        end do; end do; end do
    end subroutine fill_index_pattern

    subroutine fill_costheta_pattern(V, g)
        ! H%y (EDGE, mid-theta) -- fill with cos(mid-theta), independent of i,k
        type(cvector), intent(inout) :: V
        type(grid_t), intent(in) :: g
        integer :: i,j,k
        real(8) :: th_mid
        do j = 1,size(V%y,2)
            th_mid = g%th(j) + g%dt(j)/2.0d0
            do k = 1,size(V%y,3); do i = 1,size(V%y,1)
                V%y(i,j,k) = dcmplx(cos(th_mid), 0.0d0)
            end do; end do
        end do
    end subroutine fill_costheta_pattern

    subroutine check_theta_reversal(label, new_arr, old_arr, negate)
        character(*), intent(in) :: label
        complex(prec), dimension(:,:,:), intent(in) :: new_arr, old_arr
        logical, intent(in) :: negate
        integer :: i,j,k,n2,jold
        real(prec) :: diff, sgn
        logical :: ok
        n2 = size(new_arr,2)
        sgn = 1.0_prec
        if (negate) sgn = -1.0_prec
        ok = .true.
        do k = 1,size(new_arr,3); do j = 1,n2; do i = 1,size(new_arr,1)
            jold = n2+1-j
            diff = abs(new_arr(i,j,k) - sgn*old_arr(i,jold,k))
            if (diff > 1.0e-10_prec) ok = .false.
        end do; end do; end do
        if (ok) then
            write(*,'(a,a,a,i3)') '  ', label, ':  PASS  (N_theta=',n2
        else
            write(*,'(a,a,a)') '  ', label, ':  FAIL'
            all_pass = .false.
        end if
    end subroutine check_theta_reversal

    subroutine check_r_reversal(label, new_arr, old_arr, negate)
        character(*), intent(in) :: label
        complex(prec), dimension(:,:,:), intent(in) :: new_arr, old_arr
        logical, intent(in) :: negate
        integer :: i,j,k,n3,kold
        real(prec) :: diff, sgn
        logical :: ok
        n3 = size(new_arr,3)
        sgn = 1.0_prec
        if (negate) sgn = -1.0_prec
        ok = .true.
        do k = 1,n3; do j = 1,size(new_arr,2); do i = 1,size(new_arr,1)
            kold = n3+1-k
            diff = abs(new_arr(i,j,k) - sgn*old_arr(i,j,kold))
            if (diff > 1.0e-10_prec) ok = .false.
        end do; end do; end do
        if (ok) then
            write(*,'(a,a,a,i3)') '  ', label, ':  PASS  (N_r=',n3
        else
            write(*,'(a,a,a)') '  ', label, ':  FAIL'
            all_pass = .false.
        end if
    end subroutine check_r_reversal

    subroutine check_costheta(new_y, old_y, g)
        complex(prec), dimension(:,:,:), intent(in) :: new_y, old_y
        type(grid_t), intent(in) :: g
        integer :: j, n2, jold, jcheck
        real(prec) :: expect, actual, diff
        logical :: ok
        n2 = size(new_y,2)
        ok = .true.
        ! pick an off-equator node (not the middle) to avoid a degenerate check
        jcheck = 2
        jold = n2+1-jcheck
        expect = -cos(g%th(jold) + g%dt(jold)/2.0d0)
        actual = dble(new_y(1,jcheck,1))
        diff = abs(actual-expect)
        if (diff > 1.0e-10_prec) ok = .false.
        write(*,'(a,i3,a,f12.8,a,f12.8)') '  at j=',jcheck,':  actual=',actual,'  expect=-cos(theta_mirrored)=',expect
        if (ok) then
            write(*,*) '  H%y physical reversal check:  PASS'
        else
            write(*,*) '  H%y physical reversal check:  FAIL'
            all_pass = .false.
        end if
    end subroutine check_costheta

end program test_reversals
