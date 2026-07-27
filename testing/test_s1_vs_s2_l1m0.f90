program test_s1_vs_s2_l1m0
! ****************************************************************************
! Direct side-by-side comparison of field1d.f90's sourceField1d (the ORIGINAL
! solver, Uyeshima & Schultz-style normalization) against field1d_s2.f90's
! sourceField1d_s2 (the independently-derived, independently-validated --
! see TESTING_MANUAL.md sections 1-2 -- Sun & Egbert 2012 solver), for the
! SAME model, grid, and pure (l=1, m=0 = "P10"/zonal) unit source.
!
! Motivation: field1d.f90 had a long-unresolved "absolute sign" question for
! exactly this l=1,m=0 (E-field phase/sign convention for the P10 term) --
! the ModEM ground-truth comparison found Mode1's E_phi imaginary part had
! the OPPOSITE sign from what MATLAB/Fortran produced, never conclusively
! resolved at the time. field1d_s2.f90 is now independently validated in an
! ABSOLUTE sense (closed form + pythonSolver cross-check, both to ~5
! significant figures), so it can serve as ground truth: the two solvers'
! outputs are expected to differ by an overall complex constant per
! component-type (different normalization conventions), but that constant
! should be the SAME sign/phase across all 5 components if field1d.f90's own
! internal sign bookkeeping is consistent for this term; a component whose
! ratio breaks that pattern points at exactly where the absolute sign issue
! lives.
!
! Model/grid/staggering are IDENTICAL to test_earth_l1mneg1.f90 (uniform 100
! Ohm.m sphere, r0=6.371e6 m, T=1000s, H on EDGE / E on FACE), so the two
! solvers are evaluated at EXACTLY the same (r,theta,phi) points -- only the
! source m changes (0 instead of -1) and both solvers run on the same model.
! ****************************************************************************

    use field1d, only: conf1d_t, sourceField1d
    use field1d_s2, only: sourceField1d_s2
    use griddef
    use sg_vector
    implicit none

    type(conf1d_t)                 :: earth
    type(grid_t)                   :: grid
    type(cvector)                  :: h1d_old, e1d_old, h1d_s2, e1d_s2
    complex(8), allocatable        :: coeff(:)
    integer, parameter             :: lmax = 1
    integer, parameter             :: nx = 8, ny = 8, nz = 2
    integer                        :: ncoeff, istat, ii, jj
    real(8)                        :: period, omega
    real(8)                        :: th_lo_deg, th_hi_deg
    integer                        :: i0, j0, kRr, kRs

    ! ---- physical test parameters (identical to test_earth_l1mneg1.f90) ----
    real(8), parameter             :: r0_m    = 6.371d6      ! Earth radius, m
    real(8), parameter             :: sigma1  = 1.0d-2       ! S/m (100 Ohm.m)
    real(8), parameter             :: period0 = 1000.0d0     ! s

    write(*,*) '=== test_s1_vs_s2_l1m0: field1d.f90 vs field1d_s2.f90 (l=1,m=0, P10) ==='

    ! ---- build the earth model: ONE homogeneous layer = the whole sphere ----
    earth%r0   = r0_m
    earth%tau  = 0.0d0
    earth%tol  = 1.0d-9
    allocate(earth%layer(1), earth%sigma(1), STAT=istat)
    earth%layer(1) = 0.0d0
    earth%sigma(1) = sigma1
    earth%allocated = .true.

    ! ---- build a small hand-made grid entirely in the vacuum region above r0 ----
    call create_grid(nx, ny, nz, grid)
    grid%nx = nx; grid%ny = ny; grid%nz = nz
    grid%nzAir = nz; grid%nzCrust = 0; grid%nzEarth = 0

    do ii = 1, nx+1
        grid%ph(ii) = (ii-1) * (360.0d0/nx) * d2r
    end do
    do ii = 1, nx
        grid%dp(ii) = grid%ph(ii+1) - grid%ph(ii)
    end do

    th_lo_deg = 60.0d0; th_hi_deg = 150.0d0
    do jj = 1, ny+1
        grid%th(jj) = (th_lo_deg + (jj-1)*(th_hi_deg-th_lo_deg)/ny) * d2r
    end do
    do jj = 1, ny
        grid%dt(jj) = grid%th(jj+1) - grid%th(jj)
    end do

    grid%r(1) = 6374.00d0
    grid%r(2) = 6373.00d0
    grid%r(3) = 6371.50d0
    grid%dr(1) = grid%r(1) - grid%r(2)
    grid%dr(2) = grid%r(2) - grid%r(3)
    grid%allocated = .true.

    earth%rmax = 1.0d3 * grid%r(1)

    ! ---- unit source: pure (l=1, m=0) multipole ("P10"), amplitude 1+0i ----
    ncoeff = (lmax+1)**2
    allocate(coeff(ncoeff), STAT=istat)
    coeff = dcmplx(0.0d0, 0.0d0)
    ! ordering within each l block is m=0,+1,-1,+2,-2,...
    ! l=0: index 1 (unused, no monopole); l=1: indices 2,3,4 (m=0,+1,-1)
    coeff(2) = dcmplx(1.0d0, 0.0d0)   ! l=1, m=0

    omega  = 2.0d0*pi/period0
    period = period0

    ! ---- run BOTH solvers (H on EDGE, E on FACE -- native staggering for both) ----
    call create_cvector(grid, h1d_old, EDGE)
    call create_cvector(grid, e1d_old, FACE)
    call sourceField1d(earth, lmax, coeff, period, grid, h1d_old, e1d_old)

    call create_cvector(grid, h1d_s2, EDGE)
    call create_cvector(grid, e1d_s2, FACE)
    call sourceField1d_s2(earth, lmax, coeff, period, grid, h1d_s2, e1d_s2)

    ! ---- extract the 5 field components at one shared (i=1, j=5) angular index ----
    ! (identical staggering/indices to test_earth_l1mneg1.f90):
    !   H%x (phi,   EDGE x): node_th  x mid_ph  x Rs
    !   H%y (theta, EDGE y): mid_th   x node_ph x Rs
    !   H%z (r,     EDGE z): node_th  x node_ph x Rr
    !   E%y (theta, FACE y): node_th  x mid_ph  x Rr
    !   E%x (phi,   FACE x): mid_th   x node_ph x Rr
    i0 = 1
    j0 = 5
    kRs = 2
    kRr = 2

    write(*,*)
    write(*,'(a,3(a24))') '  component  ', 'field1d.f90 (old)', 'field1d_s2.f90', 'ratio old/s2'
    call report('Hr    ', h1d_old%z(i0,j0,kRr), h1d_s2%z(i0,j0,kRr))
    call report('Hphi  ', h1d_old%x(i0,j0,kRs), h1d_s2%x(i0,j0,kRs))
    call report('Htheta', h1d_old%y(i0,j0,kRs), h1d_s2%y(i0,j0,kRs))
    call report('Etheta', e1d_old%y(i0,j0,kRr), e1d_s2%y(i0,j0,kRr))
    call report('Ephi  ', e1d_old%x(i0,j0,kRr), e1d_s2%x(i0,j0,kRr))

    write(*,*)
    write(*,*) 'If field1d.f90 is internally sign-consistent for this term, all 5 ratios'
    write(*,*) 'above should share the same phase (angle) -- a common REAL positive or'
    write(*,*) 'negative multiple is fine (different normalization convention), but any'
    write(*,*) 'component whose ratio phase differs from the others pinpoints exactly'
    write(*,*) 'where the absolute-sign issue lives.'

    deallocate(coeff, STAT=istat)
    deallocate(earth%layer, earth%sigma, STAT=istat)
    call deall_cvector(h1d_old)
    call deall_cvector(e1d_old)
    call deall_cvector(h1d_s2)
    call deall_cvector(e1d_s2)
    call deall_grid(grid)

contains

    subroutine report(name, v_old, v_s2)
        character(*), intent(in) :: name
        complex(8), intent(in)   :: v_old, v_s2
        complex(8)               :: ratio
        real(8)                  :: ang

        ratio = v_old / v_s2
        ang = atan2(aimag(ratio), dble(ratio)) * 180.0d0 / pi
        write(*,'(2x,a,2x,2es12.4,2x,2es12.4,2x,2es12.4,a,f7.2,a)') &
            name, v_old, v_s2, ratio, '  (angle=', ang, ' deg)'
    end subroutine report

end program test_s1_vs_s2_l1m0
