program test_unit_sphere
! ****************************************************************************
! Validates sourceField1d (field1d.f90) against the closed-form analytic
! solution for a single homogeneous conducting sphere of radius r0=earth%r0,
! conductivity sigma, in vacuum, driven by a UNIT external multipole source
! of a single (l,m) degree (l=2, m=1). No thin sheet (tau=0), no separate
! core layer -- the single earth%layer/earth%sigma entry IS the whole sphere.
!
! This isolates exactly the same physics validated independently in Python
! (Sun & Egbert 2012 eq. 5-6, single conducting sphere, unit external
! r^(l+1) inducing field) -- see companion script reference_unit_sphere.py,
! which must be run with IDENTICAL physical parameters (r0, sigma, omega,
! mu0=1.256637e-6 as hardcoded in field1d.f90) to produce the expected
! values printed below.
!
! Grid is built by hand (no grid file) with a small 8x8x2 patch entirely
! above the sphere (vacuum region, r > r0), so every extracted field sample
! sits in the region covered by the closed-form solution.
! ****************************************************************************

    use field1d
    use griddef
    use sg_vector
    implicit none

    type(conf1d_t)                 :: earth
    type(grid_t)                   :: grid
    type(cvector)                  :: h1d, e1d
    complex(8), allocatable        :: coeff(:)
    integer, parameter             :: lmax = 2
    integer, parameter             :: nx = 8, ny = 8, nz = 2
    integer                        :: ncoeff, istat, ii, jj
    real(8)                        :: period, omega
    real(8)                        :: th_lo_deg, th_hi_deg
    integer                        :: i0, j0, kRr, kRs

    ! ---- physical test parameters (must match reference_unit_sphere.py) ----
    real(8), parameter             :: r0_m    = 1000.0d0     ! sphere radius, m
    real(8), parameter             :: sigma1  = 1.0d0        ! S/m
    real(8), parameter             :: omega0  = 1.0d0        ! rad/s

    write(*,*) '=== test_unit_sphere: field1d.f90 vs analytic single-sphere solution ==='

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

    ! radius, in km, decreasing top->bottom, entirely above r0=1 km
    grid%r(1) = 3.00d0
    grid%r(2) = 2.00d0
    grid%r(3) = 1.05d0
    grid%dr(1) = grid%r(1) - grid%r(2)
    grid%dr(2) = grid%r(2) - grid%r(3)
    grid%allocated = .true.

    earth%rmax = 1.0d3 * grid%r(1)

    ! ---- unit source: pure (l=2, m=+1) multipole, amplitude 1+0i ----
    ncoeff = (lmax+1)**2
    allocate(coeff(ncoeff), STAT=istat)
    coeff = dcmplx(0.0d0, 0.0d0)
    ! ordering within each l block is m=0,+1,-1,+2,-2,...
    ! l=0: index 1 (unused, no monopole)
    ! l=1: indices 2,3,4      (m=0,+1,-1)
    ! l=2: indices 5,6,7,8,9  (m=0,+1,-1,+2,-2)
    coeff(6) = dcmplx(1.0d0, 0.0d0)   ! l=2, m=+1

    omega  = omega0
    period = 2.0d0*pi/omega

    ! ---- run the solver (H on EDGE, E on FACE -- field1d.f90's native staggering) ----
    call create_cvector(grid, h1d, EDGE)
    call create_cvector(grid, e1d, FACE)
    call sourceField1d(earth, lmax, coeff, period, grid, h1d, e1d)

    ! ---- extract the 5 field components at one shared (i=1, j=5) angular index ----
    ! matching the staggering documented in field1d.f90's EDGE-branch header comment:
    !   H%x (phi,   EDGE x): node_th  x mid_ph  x Rs
    !   H%y (theta, EDGE y): mid_th   x node_ph x Rs
    !   H%z (r,     EDGE z): node_th  x node_ph x Rr
    !   E%y (theta, FACE y): node_th  x mid_ph  x Rr
    !   E%x (phi,   FACE x): mid_th   x node_ph x Rr
    i0 = 1
    j0 = 5
    kRs = 2   ! Rs(2) = grid%r(2)*1e3 = 2000 m
    kRr = 2   ! Rr(2) = Rs(2) - dr(2)/2*1e3 = 1525 m

    write(*,*)
    write(*,*) 'Grid check:'
    write(*,'(a,f18.12)') '  grid%th(5)      [rad, deg] = ', grid%th(5)
    write(*,'(a,f18.12)') '  grid%th(5)+dt/2 [rad]      = ', grid%th(5)+grid%dt(5)/2.0d0
    write(*,'(a,f18.12)') '  grid%ph(1)                 = ', grid%ph(1)
    write(*,'(a,f18.12)') '  grid%ph(1)+dp/2            = ', grid%ph(1)+grid%dp(1)/2.0d0
    write(*,'(a,f18.6)')  '  Rs(2)  [m]                 = ', 1.0d3*grid%r(2)
    write(*,'(a,f18.6)')  '  Rr(2)  [m]                 = ', 1.0d3*grid%r(2) - 1.0d3*grid%dr(2)/2.0d0

    write(*,*)
    write(*,*) 'FIELD1D output (real, imag):'
    write(*,'(a,2es24.15)') '  H%z(i=1,j=5,k=Rr=2) = ', h1d%z(i0,j0,kRr)
    write(*,'(a,2es24.15)') '  H%x(i=1,j=5,k=Rs=2) = ', h1d%x(i0,j0,kRs)
    write(*,'(a,2es24.15)') '  H%y(i=1,j=5,k=Rs=2) = ', h1d%y(i0,j0,kRs)
    write(*,'(a,2es24.15)') '  E%y(i=1,j=5,k=Rr=2) = ', e1d%y(i0,j0,kRr)
    write(*,'(a,2es24.15)') '  E%x(i=1,j=5,k=Rr=2) = ', e1d%x(i0,j0,kRr)

    write(*,*)
    write(*,*) 'Compare the five complex numbers above against reference_unit_sphere.py'
    write(*,*) '(run with r0=1000, sigma=1, omega=1, l=2, m=1 -- see header of that script).'
    write(*,*) 'Report the ratio FIELD1D / reference for each component.'

    deallocate(coeff, STAT=istat)
    deallocate(earth%layer, earth%sigma, STAT=istat)
    call deall_cvector(h1d)
    call deall_cvector(e1d)
    call deall_grid(grid)

end program test_unit_sphere
