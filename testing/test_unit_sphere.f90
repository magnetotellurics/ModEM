program test_unit_sphere
! ****************************************************************************
! Validates BOTH sourceField1d (field1d.f90, S1) and
! sourceField1d_s2 (field1d_s2.f90, S2) against the
! SAME closed-form analytic solution for a single homogeneous conducting
! sphere of radius r0=earth%r0, conductivity sigma, in vacuum, driven by a
! UNIT external multipole source of a single (l,m) degree (l=2, m=+1). No
! thin sheet (tau=0), no separate core layer -- the single earth%layer/
! earth%sigma entry IS the whole sphere. (Merged 2026-07-25 from the former
! separate test_unit_sphere_s1.f90 and test_unit_sphere_s2.f90 --
! both built the identical model/grid/source, so there was no reason to keep
! them apart. See testing/reference_unit_sphere.py for the merged closed-form
! reference.)
!
! IMPORTANT -- the two solvers are NOT expected to match this closed form (or
! each other) the same way:
!   - S2 uses the paper's own literal "alpha_l^T=1" (unit external
!     r^(l+1) amplitude) normalization directly, with no extra bookkeeping,
!     so the match is EXACT: ratio = 1+0i for all five components.
!   - S1 carries its own internal normalization (tni-referenced,
!     R0^2/l(l+1) factors) AND applies a final conjg() to every assembled
!     component that was designed to compensate for a "+m,-m conjugate-
!     pairing" reconstruction trick -- but this test (like all the other
!     single-(l,m) tests in this suite) sources only ONE of the +m/-m pair,
!     so that pairing-compensation is never actually exercised as intended.
!     As analyzed 2026-07-25, the result is NOT a clean
!     scalar/phase discrepancy for m!=0: it reconstructs a DIFFERENT angular
!     pattern (the m-flipped one) combined with a conjugated radial function,
!     not just a scaled/rotated version of the correct field. So S1's
!     ratio against this closed form is expected to be neither 1+0i nor a
!     single common phase across all five components -- this is the open
!     "absolute sign"/conjg() issue under investigation, NOT a test bug.
!     Report the raw S1 numbers for visual inspection only.
!
! Grid is built by hand (no grid file) with a small 8x8x2 patch entirely
! above the sphere (vacuum region, r > r0), so every extracted field sample
! sits in the region covered by the closed-form solution.
! ****************************************************************************

    use field1d, only: conf1d_t, sourceField1d
    use field1d_s2, only: sourceField1d_s2
    use griddef
    use sg_vector
    implicit none

    type(conf1d_t)                 :: earth
    type(grid_t)                   :: grid
    type(cvector)                  :: h1d_s1, e1d_s1, h1d_s2, e1d_s2
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

    write(*,*) '=== test_unit_sphere: field1d.f90 (S1) and field1d_s2.f90 (S2)'
    write(*,*) '    vs analytic single-sphere solution ==='

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
    ! ordering within each l block is m=0,+1,-1,+2,-2
    ! l=0: index 1 (unused, no monopole)
    ! l=1: indices 2,3,4      (m=0,+1,-1)
    ! l=2: indices 5,6,7,8,9  (m=0,+1,-1,+2,-2)
    coeff(6) = dcmplx(1.0d0, 0.0d0)   ! l=2, m=+1

    omega  = omega0
    period = 2.0d0*pi/omega

    ! ---- run BOTH solvers (H on EDGE, E on FACE -- native staggering for both) ----
    call create_cvector(grid, h1d_s1, EDGE)
    call create_cvector(grid, e1d_s1, FACE)
    call sourceField1d(earth, lmax, coeff, period, grid, h1d_s1, e1d_s1)

    call create_cvector(grid, h1d_s2, EDGE)
    call create_cvector(grid, e1d_s2, FACE)
    call sourceField1d_s2(earth, lmax, coeff, period, grid, h1d_s2, e1d_s2)

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
    write(*,*) 'S1 (field1d.f90) output (real, imag) -- compare vs reference_unit_sphere.py'
    write(*,*) 'section 1 ("field1d.f90 (S1)"); NOT expected to match exactly, see header note:'
    write(*,'(a,2es24.15)') '  H%z(i=1,j=5,k=Rr=2) [Hr]     = ', h1d_s1%z(i0,j0,kRr)
    write(*,'(a,2es24.15)') '  H%x(i=1,j=5,k=Rs=2) [Hphi]   = ', h1d_s1%x(i0,j0,kRs)
    write(*,'(a,2es24.15)') '  H%y(i=1,j=5,k=Rs=2) [Htheta] = ', h1d_s1%y(i0,j0,kRs)
    write(*,'(a,2es24.15)') '  E%y(i=1,j=5,k=Rr=2) [Etheta] = ', e1d_s1%y(i0,j0,kRr)
    write(*,'(a,2es24.15)') '  E%x(i=1,j=5,k=Rr=2) [Ephi]   = ', e1d_s1%x(i0,j0,kRr)

    write(*,*)
    write(*,*) 'S2 (field1d_s2.f90) output (real, imag) -- compare vs'
    write(*,*) 'reference_unit_sphere.py section 2 ("field1d_s2.f90 (S2)");'
    write(*,*) 'expect EXACT match, ratio = 1+0i for all five components:'
    write(*,'(a,2es24.15)') '  H%z(i=1,j=5,k=Rr=2) [Hr]     = ', h1d_s2%z(i0,j0,kRr)
    write(*,'(a,2es24.15)') '  H%x(i=1,j=5,k=Rs=2) [Hphi]   = ', h1d_s2%x(i0,j0,kRs)
    write(*,'(a,2es24.15)') '  H%y(i=1,j=5,k=Rs=2) [Htheta] = ', h1d_s2%y(i0,j0,kRs)
    write(*,'(a,2es24.15)') '  E%y(i=1,j=5,k=Rr=2) [Etheta] = ', e1d_s2%y(i0,j0,kRr)
    write(*,'(a,2es24.15)') '  E%x(i=1,j=5,k=Rr=2) [Ephi]   = ', e1d_s2%x(i0,j0,kRr)

    deallocate(coeff, STAT=istat)
    deallocate(earth%layer, earth%sigma, STAT=istat)
    call deall_cvector(h1d_s1)
    call deall_cvector(e1d_s1)
    call deall_cvector(h1d_s2)
    call deall_cvector(e1d_s2)
    call deall_grid(grid)

end program test_unit_sphere
