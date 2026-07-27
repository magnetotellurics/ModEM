program test_earth_l1mneg1
! ****************************************************************************
! Cross-check of field1d_s2.f90 against pythonSolver/spherical_em_induction.py
! (fields_from_R_general, now fixed -- see 2026-07-22/23 fixes) at Earth scale
! and physical parameters already used repeatedly earlier in this project:
! r0=6.371e6 m, uniform 100 Ohm.m from centre to surface, T=1000 s. Source is
! a PURE (l=1, m=-1) unit multipole -- the original, historically unresolved
! comparison case from earlier in this conversation ("extra (-1i)-type factor
! on Hz,Ex,Ey").
!
! Grid is hand-built (no grid file), entirely in the vacuum region above r0,
! same pattern as test_unit_sphere_s2.f90 (which already validated
! field1d_s2 for l=2,m=1 against a closed-form single-sphere solution to
! ~5 significant figures).
!
! Companion script: reference_earth_l1mneg1.py -- uses spherical_em_induction's
! own solve_layered/fields_from_R_general directly (not a separately re-derived
! closed form this time, since those are now independently validated), with
! the SAME uniform 100 Ohm.m model, normalized via solve_layered's own
! returned "A" (external amplitude) so the K0/source-shell trick's accuracy
! doesn't matter -- and the e^{+iwt}->e^{-iwt} conjugation needed to compare
! against field1d_s2.f90's native convention.
! ****************************************************************************

    use field1d_s2
    use field1d, only: conf1d_t
    use griddef
    use sg_vector
    implicit none

    type(conf1d_t)                 :: earth
    type(grid_t)                   :: grid
    type(cvector)                  :: h1d, e1d
    complex(8), allocatable        :: coeff(:)
    integer, parameter             :: lmax = 1
    integer, parameter             :: nx = 8, ny = 8, nz = 2
    integer                        :: ncoeff, istat, ii, jj
    real(8)                        :: period, omega
    real(8)                        :: th_lo_deg, th_hi_deg
    integer                        :: i0, j0, kRr, kRs

    ! ---- physical test parameters (must match reference_earth_l1mneg1.py) ----
    real(8), parameter             :: r0_m    = 6.371d6      ! Earth radius, m
    real(8), parameter             :: sigma1  = 1.0d-2       ! S/m (100 Ohm.m)
    real(8), parameter             :: period0 = 1000.0d0     ! s

    write(*,*) '=== test_earth_l1mneg1: field1d_s2.f90 vs spherical_em_induction.py (l=1,m=-1) ==='

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

    ! longitude: 0..360 deg, 8 cells
    do ii = 1, nx+1
        grid%ph(ii) = (ii-1) * (360.0d0/nx) * d2r
    end do
    do ii = 1, nx
        grid%dp(ii) = grid%ph(ii+1) - grid%ph(ii)
    end do

    ! co-latitude: 60..150 deg, 8 cells (well away from both poles)
    th_lo_deg = 60.0d0; th_hi_deg = 150.0d0
    do jj = 1, ny+1
        grid%th(jj) = (th_lo_deg + (jj-1)*(th_hi_deg-th_lo_deg)/ny) * d2r
    end do
    do jj = 1, ny
        grid%dt(jj) = grid%th(jj+1) - grid%th(jj)
    end do

    ! radius, in km, decreasing top->bottom, entirely above r0 = 6371.0 km
    grid%r(1) = 6374.00d0
    grid%r(2) = 6373.00d0
    grid%r(3) = 6371.50d0
    grid%dr(1) = grid%r(1) - grid%r(2)
    grid%dr(2) = grid%r(2) - grid%r(3)
    grid%allocated = .true.

    earth%rmax = 1.0d3 * grid%r(1)

    ! ---- unit source: pure (l=1, m=-1) multipole, amplitude 1+0i ----
    ncoeff = (lmax+1)**2
    allocate(coeff(ncoeff), STAT=istat)
    coeff = dcmplx(0.0d0, 0.0d0)
    ! ordering within each l block is m=0,+1,-1,+2,-2,...
    ! l=0: index 1 (unused, no monopole)
    ! l=1: indices 2,3,4      (m=0,+1,-1)
    coeff(4) = dcmplx(1.0d0, 0.0d0)   ! l=1, m=-1

    omega  = 2.0d0*pi/period0
    period = period0

    ! ---- run the solver (H on EDGE, E on FACE -- field1d_s2's native staggering) ----
    call create_cvector(grid, h1d, EDGE)
    call create_cvector(grid, e1d, FACE)
    call sourceField1d_s2(earth, lmax, coeff, period, grid, h1d, e1d)

    ! ---- extract the 5 field components at one shared (i=1, j=5) angular index ----
    ! (same staggering/indices as test_unit_sphere_s2.f90):
    !   H%x (phi,   EDGE x): node_th  x mid_ph  x Rs
    !   H%y (theta, EDGE y): mid_th   x node_ph x Rs
    !   H%z (r,     EDGE z): node_th  x node_ph x Rr
    !   E%y (theta, FACE y): node_th  x mid_ph  x Rr
    !   E%x (phi,   FACE x): mid_th   x node_ph x Rr
    i0 = 1
    j0 = 5
    kRs = 2   ! Rs(2) = grid%r(2)*1e3 m
    kRr = 2   ! Rr(2) = Rs(2) - dr(2)/2*1e3 m

    write(*,*)
    write(*,*) 'Grid check:'
    write(*,'(a,f18.12)') '  grid%th(5)      [rad] = ', grid%th(5)
    write(*,'(a,f18.12)') '  grid%th(5)+dt/2 [rad] = ', grid%th(5)+grid%dt(5)/2.0d0
    write(*,'(a,f18.12)') '  grid%ph(1)             = ', grid%ph(1)
    write(*,'(a,f18.12)') '  grid%ph(1)+dp/2        = ', grid%ph(1)+grid%dp(1)/2.0d0
    write(*,'(a,f18.3)')  '  Rs(2)  [m]             = ', 1.0d3*grid%r(2)
    write(*,'(a,f18.3)')  '  Rr(2)  [m]             = ', 1.0d3*grid%r(2) - 1.0d3*grid%dr(2)/2.0d0

    write(*,*)
    write(*,*) 'FIELD1D_SE12 output (real, imag):'
    write(*,'(a,2es24.15)') '  H%z(i=1,j=5,k=Rr=2) [Hr]     = ', h1d%z(i0,j0,kRr)
    write(*,'(a,2es24.15)') '  H%x(i=1,j=5,k=Rs=2) [Hphi]   = ', h1d%x(i0,j0,kRs)
    write(*,'(a,2es24.15)') '  H%y(i=1,j=5,k=Rs=2) [Htheta] = ', h1d%y(i0,j0,kRs)
    write(*,'(a,2es24.15)') '  E%y(i=1,j=5,k=Rr=2) [Etheta] = ', e1d%y(i0,j0,kRr)
    write(*,'(a,2es24.15)') '  E%x(i=1,j=5,k=Rr=2) [Ephi]   = ', e1d%x(i0,j0,kRr)

    write(*,*)
    write(*,*) 'Compare against: python reference_earth_l1mneg1.py'

    deallocate(coeff, STAT=istat)
    deallocate(earth%layer, earth%sigma, STAT=istat)
    call deall_cvector(h1d)
    call deall_cvector(e1d)
    call deall_grid(grid)

end program test_earth_l1mneg1
