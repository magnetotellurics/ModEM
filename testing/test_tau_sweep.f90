program test_tau_sweep
! Sweeps earth%tau (thin-sheet surface conductance -- see field1d_s2.f90's
! eq(17) jump condition T'(r0+) = T'(r0-) - i*omega*mu0*tau*T(r0)) from a
! representative range down toward zero, and measures the relative error
! AND phase shift of the magnetic field against the tau=0 case, which IS
! the analytic closed-form solution: S2 matches the closed-form single-
! sphere reference EXACTLY at tau=0 (see test_unit_sphere.f90 /
! reference_unit_sphere.py, ratio=1+0i for all 5 field components). So
! tau=0's own S2 output can be used directly as the "analytic" baseline
! here, with no need to separately re-derive/hardcode the closed-form
! numbers in this driver.
!
! SAME sphere/grid/source setup as test_unit_sphere.f90 (l=2,m=+1 unit
! source, r0=1000 m, sigma=1 S/m, omega=1 rad/s), reused verbatim so the
! tau=0 baseline is identical to that already-validated exact match.
!
! Also reports the SAME magnitude/phase tables at Earth's actual radius
! (r0=6.371e6 m), at the NA2026 project's own periods (162s, 3600s), to
! test whether the toy sphere's own Q=omega*mu0*tau*r0 ~ O(1) nonlinearity
! threshold and error~0.024*Q linear coefficient generalize, or are
! specific to that sphere's own conductor regime. sigma at each Earth-
! scale period is chosen so |kl|*r0 = sqrt(omega*mu0*sigma)*r0 EQUALS the
! toy sphere's own value: this isolates the tau/Q effect cleanly, since
! naively reusing sigma=1 S/m at Earth's radius would ALSO push the
! conductor's own electrical size |kl|*r0 from ~1.12 (order-1 induction,
! same regime as the toy test) to ~1400 (extreme skin effect) -- a second,
! unrelated regime change that would confound any observed nonlinearity.
!
! Also reports the same table using a REALISTIC layered Earth conductivity
! profile (not a matched-regime homogeneous sphere): the depths/log10(rho)
! values from LWS/layered_GDE_rho.prm, read the same way FWD1D.f90 builds
! earth%layer/earth%sigma from a "layers"-format model file (surface-down
! depths in km, sigma=10**(-log10rho) per layer, plus a hardcoded
! sigma=1e5 S/m "core" layer below the deepest defined depth), at r0 =
! R0_EARTH and period = 3600 s.
!
! Backs docs/tau_sensitivity_analysis.tex/.md/.pdf (2026-08): confirms (i)
! no numerical advantage to a nonzero tau -- the solver is exact and
! well-behaved at tau=0 -- and (ii) the error/phase perturbation scales
! linearly with the dimensionless thin-sheet parameter Q=omega*mu0*tau*r0
! up to Q~O(1). Table 1's output matches that report's Table 1 exactly;
! the toy-sphere block of Table 2 output matches that report's Tables 2-4;
! the matched-|kl|r0 Earth-scale blocks match its Tables 6-7; the GDE
! layered-model depth/sigma printout matches its Table 8 and the results
! block matches its Table 9 (sigma=1, omega=1 elsewhere unless otherwise
! stated).

    use field1d, only: conf1d_t
    use field1d_s2, only: sourceField1d_s2
    use griddef
    use sg_vector
    use math_constants, only: MU_0
    implicit none

    type(conf1d_t)                 :: earth
    type(grid_t)                   :: grid
    type(cvector)                  :: h1d, e1d, h1d_ref, e1d_ref
    complex(8), allocatable        :: coeff(:)
    integer, parameter             :: lmax = 2
    integer, parameter             :: nx = 8, ny = 8, nz = 2
    integer                        :: ncoeff, istat, ii, jj, itau
    real(8)                        :: period, omega
    real(8)                        :: th_lo_deg, th_hi_deg
    integer                        :: i0, j0, kRr, kRs

    ! ---- physical test parameters (must match test_unit_sphere.f90) ----
    real(8), parameter             :: r0_m    = 1000.0d0     ! sphere radius, m
    real(8), parameter             :: sigma1  = 1.0d0        ! S/m
    real(8), parameter             :: omega0  = 1.0d0        ! rad/s
    real(8), parameter             :: RAD2DEG = 180.0d0 / 3.14159265358979324d0
    real(8), parameter             :: TWOPI   = 2.0d0 * 3.14159265358979324d0
    real(8), parameter             :: R0_EARTH = 6.371d6     ! Earth radius, m

    integer, parameter             :: ntau = 10
    real(8), parameter             :: tau_sweep(ntau) = &
        (/ 1.0d0, 1.0d-1, 1.0d-2, 1.0d-3, 1.0d-4, 1.0d-5, 1.0d-6, 1.0d-7, 1.0d-8, 0.0d0 /)
    real(8)                        :: Q, errHr, errHth, errHph, errEth, errEph, maxerr
    complex(8)                     :: refHr, refHth, refHph, refEth, refEph
    complex(8)                     :: curHr, curHth, curHph, curEth, curEph

    integer, parameter             :: ntau2 = 4
    real(8), parameter             :: tau_report(ntau2) = &
        (/ 1.0d-4, 1.0d0, 1.0d2, 2.5d3 /)

    real(8)                        :: klr0_toy, sigma_earth

    ! ---- realistic layered Earth conductivity model (LWS/layered_GDE_rho.prm) ----
    ! Layer depths [km] (top of layer, measured from the surface) and
    ! log10(rho[Ohm-m]) values, transcribed directly from the "log degree 0
    ! layer <depth>" / "0 0 <value>" blocks of that file (32 layers). FWD1D.f90
    ! reads a file of this format and builds earth%sigma = 10**(-logrho), with
    ! an additional sigma=1e5 S/m "core" layer appended below the deepest
    ! defined depth (3500 km) -- reproduced identically below.
    integer, parameter             :: nL_gde = 32
    real(8), parameter             :: depth_gde(nL_gde) = &
        (/   10.00d0,   40.00d0,  110.00d0,  160.00d0,  210.00d0,  260.00d0, &
            310.00d0,  360.00d0,  410.00d0,  460.00d0,  510.00d0,  560.00d0, &
            610.00d0,  660.00d0,  710.00d0,  760.00d0,  800.00d0,  850.00d0, &
            900.00d0, 1000.00d0, 1100.00d0, 1200.00d0, 1300.00d0, 1400.00d0, &
           1500.00d0, 1620.00d0, 1764.00d0, 1936.80d0, 2144.16d0, 2400.00d0, &
           2900.00d0, 3500.00d0 /)
    real(8), parameter             :: logrho_gde(nL_gde) = &
        (/  2.000000d0,  3.000000d0,  3.000000d0,  2.522879d0,  2.154902d0, &
            1.745627d0,  1.667885d0,  1.688480d0,  1.709952d0,  1.699400d0, &
            1.645104d0,  1.540801d0,  1.379495d0,  1.153394d0,  0.869989d0, &
            0.593626d0,  0.430060d0,  0.279000d0,  0.279000d0, -0.228100d0, &
           -0.228100d0, -0.228100d0, -0.228100d0, -0.228100d0, -0.228100d0, &
           -0.228100d0, -0.228100d0, -0.228100d0, -0.228100d0, -0.228100d0, &
           -1.000000d0, -2.000000d0 /)

    write(*,*) '=== test_tau_sweep: earth%tau sensitivity vs analytic (tau=0) reference ==='
    write(*,*) '    sphere r0=',r0_m,' m, sigma=',sigma1,' S/m, omega=',omega0,' rad/s'
    write(*,*)

    ! ---- build the earth model: ONE homogeneous layer = the whole sphere ----
    earth%r0   = r0_m
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
    coeff(6) = dcmplx(1.0d0, 0.0d0)   ! l=2, m=+1

    omega  = omega0
    period = 2.0d0*pi/omega

    i0 = 1
    j0 = 5
    kRs = 2
    kRr = 2

    ! ---- reference: tau=0 (the analytic/exact case) ----
    earth%tau = 0.0d0
    call create_cvector(grid, h1d_ref, EDGE)
    call create_cvector(grid, e1d_ref, FACE)
    call sourceField1d_s2(earth, lmax, coeff, period, grid, h1d_ref, e1d_ref)
    refHr  = h1d_ref%z(i0,j0,kRr)
    refHph = h1d_ref%x(i0,j0,kRs)
    refHth = h1d_ref%y(i0,j0,kRs)
    refEth = e1d_ref%y(i0,j0,kRr)
    refEph = e1d_ref%x(i0,j0,kRr)

    write(*,'(a)') '--- Table 1: magnitude-error sweep (all 5 field components) ---'
    write(*,'(a)') '  tau [S]        Q=omega*mu0*tau*r0   maxrelerr(5 comps)'
    write(*,'(a)') '  -------------  -------------------  -------------------'

    do itau = 1,ntau
        earth%tau = tau_sweep(itau)
        Q = omega * MU_0 * earth%tau * r0_m

        call create_cvector(grid, h1d, EDGE)
        call create_cvector(grid, e1d, FACE)
        call sourceField1d_s2(earth, lmax, coeff, period, grid, h1d, e1d)

        curHr  = h1d%z(i0,j0,kRr)
        curHph = h1d%x(i0,j0,kRs)
        curHth = h1d%y(i0,j0,kRs)
        curEth = e1d%y(i0,j0,kRr)
        curEph = e1d%x(i0,j0,kRr)

        errHr  = abs(curHr  - refHr ) / abs(refHr)
        errHth = abs(curHth - refHth) / abs(refHth)
        errHph = abs(curHph - refHph) / abs(refHph)
        errEth = abs(curEth - refEth) / abs(refEth)
        errEph = abs(curEph - refEph) / abs(refEph)
        maxerr = max(errHr, errHth, errHph, errEth, errEph)

        write(*,'(2x,es13.6,2x,es19.6,2x,es19.6)') earth%tau, Q, maxerr

        call deall_cvector(h1d)
        call deall_cvector(e1d)
    end do

    write(*,*)
    write(*,'(a)') '--- Table 2 (toy sphere, r0=1000m): H magnitude ratio and phase shift vs tau=0 ---'
    call tau_report_at_scale(r0_m, sigma1, omega0, refHr, refHth, refHph, 'toy sphere (r0=1000 m, omega=1 rad/s)')

    ! ---- Earth-scale scenarios: r0=R0_EARTH, sigma chosen to preserve the
    ! toy sphere's own |kl|*r0 (isolates the Q/tau effect from the
    ! conductor's own electrical-size regime -- see header comment) ----
    klr0_toy = sqrt(omega0 * MU_0 * sigma1) * r0_m

    write(*,*)
    write(*,'(a,es12.5)') 'Earth-scale scenarios: preserving toy sphere |kl|*r0 = ', klr0_toy

    omega = TWOPI / 162.0d0
    sigma_earth = (klr0_toy/R0_EARTH)**2 / (omega * MU_0)
    write(*,*)
    write(*,'(a,es12.5,a)') '--- Table (Earth, period=162s): sigma=',sigma_earth,' S/m, H ratio/phase vs tau=0 ---'
    call tau_report_earth_scale(R0_EARTH, sigma_earth, omega, 'Earth r0, period=162s')

    omega = TWOPI / 3600.0d0
    sigma_earth = (klr0_toy/R0_EARTH)**2 / (omega * MU_0)
    write(*,*)
    write(*,'(a,es12.5,a)') '--- Table (Earth, period=3600s): sigma=',sigma_earth,' S/m, H ratio/phase vs tau=0 ---'
    call tau_report_earth_scale(R0_EARTH, sigma_earth, omega, 'Earth r0, period=3600s')

    ! ---- Earth scale, REALISTIC layered GDE conductivity model, period=3600s ----
    write(*,*)
    write(*,'(a)') '--- Table (Earth, layered GDE model, period=3600s): H ratio/phase vs tau=0 ---'
    call tau_report_earth_GDE(3600.0d0, 'Earth r0, GDE layered model, period=3600s')

    deallocate(coeff, STAT=istat)
    deallocate(earth%layer, earth%sigma, STAT=istat)
    call deall_cvector(h1d_ref)
    call deall_cvector(e1d_ref)
    call deall_grid(grid)

contains

    real(8) function wrap_deg(x) result(w)
        real(8), intent(in) :: x
        w = x
        do while (w > 180.0d0)
            w = w - 360.0d0
        end do
        do while (w <= -180.0d0)
            w = w + 360.0d0
        end do
    end function wrap_deg

    ! Prints the 4-tau-value H magnitude-ratio/phase-shift table for the
    ! MODULE-level grid/coeff/period (toy sphere scale), given the already-
    ! computed tau=0 reference values ref{Hr,Hth,Hph}.
    subroutine tau_report_at_scale(r0_val, sigma_val, omega_val, refHr_val, refHth_val, refHph_val, label)
        real(8), intent(in)      :: r0_val, sigma_val, omega_val
        complex(8), intent(in)   :: refHr_val, refHth_val, refHph_val
        character(*), intent(in) :: label
        real(8)                  :: Qloc, dphHr, dphHth, dphHph, magHr, magHth, magHph
        complex(8)               :: cHr, cHth, cHph
        type(cvector)             :: hloc, eloc
        integer                   :: k

        write(*,'(a)') '  tau[S]  Q  |Hr/Hr0| dPh_Hr  |Hth/Hth0| dPh_Hth  |Hph/Hph0| dPh_Hph'
        do k = 1,ntau2
            earth%tau = tau_report(k)
            Qloc = omega_val * MU_0 * earth%tau * r0_val

            call create_cvector(grid, hloc, EDGE)
            call create_cvector(grid, eloc, FACE)
            call sourceField1d_s2(earth, lmax, coeff, 2.0d0*pi/omega_val, grid, hloc, eloc)

            cHr  = hloc%z(i0,j0,kRr)
            cHph = hloc%x(i0,j0,kRs)
            cHth = hloc%y(i0,j0,kRs)

            magHr  = abs(cHr) / abs(refHr_val)
            magHth = abs(cHth) / abs(refHth_val)
            magHph = abs(cHph) / abs(refHph_val)

            dphHr  = wrap_deg( (atan2(aimag(cHr), real(cHr))   - atan2(aimag(refHr_val), real(refHr_val)))   * RAD2DEG )
            dphHth = wrap_deg( (atan2(aimag(cHth),real(cHth))  - atan2(aimag(refHth_val),real(refHth_val)))  * RAD2DEG )
            dphHph = wrap_deg( (atan2(aimag(cHph),real(cHph))  - atan2(aimag(refHph_val),real(refHph_val)))  * RAD2DEG )

            write(*,'(2x,es11.4,2x,es11.4,2x,f10.6,2x,f15.8,2x,f11.6,2x,f16.8,2x,f11.6,2x,f16.8)') &
                earth%tau, Qloc, magHr, dphHr, magHth, dphHth, magHph, dphHph

            call deall_cvector(hloc)
            call deall_cvector(eloc)
        end do
    end subroutine tau_report_at_scale

    ! Same as tau_report_at_scale, but for an INDEPENDENT (r0,sigma,omega)
    ! scenario: builds its own rescaled grid (proportional to r0_val/r0_m,
    ! preserving the same relative Rs/Rr extraction points above the
    ! surface) and its own tau=0 reference, then reuses the same reporting
    ! loop/format.
    subroutine tau_report_earth_scale(r0_val, sigma_val, omega_val, label)
        real(8), intent(in)      :: r0_val, sigma_val, omega_val
        character(*), intent(in) :: label
        type(conf1d_t)            :: earth_loc
        type(grid_t)               :: grid_loc
        type(cvector)               :: h0, e0
        real(8)                    :: scale_factor
        complex(8)                 :: refHr_loc, refHth_loc, refHph_loc
        integer                     :: kk

        scale_factor = r0_val / r0_m

        earth_loc%r0  = r0_val
        earth_loc%tol = 1.0d-9
        allocate(earth_loc%layer(1), earth_loc%sigma(1))
        earth_loc%layer(1) = 0.0d0
        earth_loc%sigma(1) = sigma_val
        earth_loc%allocated = .true.

        call create_grid(nx, ny, nz, grid_loc)
        grid_loc%nx = nx; grid_loc%ny = ny; grid_loc%nz = nz
        grid_loc%nzAir = nz; grid_loc%nzCrust = 0; grid_loc%nzEarth = 0
        do kk = 1, nx+1
            grid_loc%ph(kk) = grid%ph(kk)
        end do
        do kk = 1, nx
            grid_loc%dp(kk) = grid%dp(kk)
        end do
        do kk = 1, ny+1
            grid_loc%th(kk) = grid%th(kk)
        end do
        do kk = 1, ny
            grid_loc%dt(kk) = grid%dt(kk)
        end do
        grid_loc%r(1) = grid%r(1) * scale_factor
        grid_loc%r(2) = grid%r(2) * scale_factor
        grid_loc%r(3) = grid%r(3) * scale_factor
        grid_loc%dr(1) = grid_loc%r(1) - grid_loc%r(2)
        grid_loc%dr(2) = grid_loc%r(2) - grid_loc%r(3)
        grid_loc%allocated = .true.
        earth_loc%rmax = 1.0d3 * grid_loc%r(1)

        earth_loc%tau = 0.0d0
        call create_cvector(grid_loc, h0, EDGE)
        call create_cvector(grid_loc, e0, FACE)
        call sourceField1d_s2(earth_loc, lmax, coeff, 2.0d0*pi/omega_val, grid_loc, h0, e0)
        refHr_loc  = h0%z(i0,j0,kRr)
        refHph_loc = h0%x(i0,j0,kRs)
        refHth_loc = h0%y(i0,j0,kRs)
        call deall_cvector(h0)
        call deall_cvector(e0)

        write(*,'(a)') '  tau[S]  Q  |Hr/Hr0| dPh_Hr  |Hth/Hth0| dPh_Hth  |Hph/Hph0| dPh_Hph'
        do kk = 1,ntau2
            earth_loc%tau = tau_report(kk)
            call run_one(earth_loc, grid_loc, refHr_loc, refHth_loc, refHph_loc, omega_val, r0_val)
        end do

        deallocate(earth_loc%layer, earth_loc%sigma)
        call deall_grid(grid_loc)
    end subroutine tau_report_earth_scale

    ! Same reporting loop as tau_report_earth_scale, but builds the REAL
    ! layered Earth conductivity model (nL_gde depth/log10rho pairs above,
    ! plus the hardcoded sigma=1e5 S/m core layer) at r0=R0_EARTH, instead of
    ! a homogeneous sphere with sigma matched to preserve |kl|*r0. The
    ! evaluation grid is the same toy grid rescaled by R0_EARTH/r0_m, so the
    ! extraction point sits at the same relative altitude above the surface
    ! as every other table in this program.
    subroutine tau_report_earth_GDE(period_val, label)
        real(8), intent(in)      :: period_val
        character(*), intent(in) :: label
        type(conf1d_t)            :: earth_gde
        type(grid_t)               :: grid_gde
        type(cvector)               :: h0, e0
        real(8)                    :: scale_factor, omega_val
        complex(8)                 :: refHr_loc, refHth_loc, refHph_loc
        integer                     :: kk

        omega_val = TWOPI / period_val
        scale_factor = R0_EARTH / r0_m

        earth_gde%r0  = R0_EARTH
        earth_gde%tol = 1.0d-9
        allocate(earth_gde%layer(nL_gde+1), earth_gde%sigma(nL_gde+1))
        earth_gde%layer(1) = 0.0d0
        earth_gde%layer(2:nL_gde+1) = 1.0d3 * depth_gde(1:nL_gde)
        earth_gde%sigma(1:nL_gde) = 10.0d0**( -logrho_gde(1:nL_gde) )
        earth_gde%sigma(nL_gde+1) = 10.0d0**5.0d0   ! core conductivity, as in FWD1D.f90
        earth_gde%allocated = .true.

        write(*,'(a)') 'Layer depths [km] (top of layer, from surface), LWS/layered_GDE_rho.prm:'
        write(*,'(8f10.2)') depth_gde
        write(*,'(a)') 'log10(rho [Ohm-m]) per layer:'
        write(*,'(8f10.6)') logrho_gde
        write(*,'(a)') 'sigma [S/m] per layer (layer 33 = hardcoded core, sigma=1e5 S/m):'
        write(*,'(8es12.4)') earth_gde%sigma

        call create_grid(nx, ny, nz, grid_gde)
        grid_gde%nx = nx; grid_gde%ny = ny; grid_gde%nz = nz
        grid_gde%nzAir = nz; grid_gde%nzCrust = 0; grid_gde%nzEarth = 0
        do kk = 1, nx+1
            grid_gde%ph(kk) = grid%ph(kk)
        end do
        do kk = 1, nx
            grid_gde%dp(kk) = grid%dp(kk)
        end do
        do kk = 1, ny+1
            grid_gde%th(kk) = grid%th(kk)
        end do
        do kk = 1, ny
            grid_gde%dt(kk) = grid%dt(kk)
        end do
        grid_gde%r(1) = grid%r(1) * scale_factor
        grid_gde%r(2) = grid%r(2) * scale_factor
        grid_gde%r(3) = grid%r(3) * scale_factor
        grid_gde%dr(1) = grid_gde%r(1) - grid_gde%r(2)
        grid_gde%dr(2) = grid_gde%r(2) - grid_gde%r(3)
        grid_gde%allocated = .true.
        earth_gde%rmax = 1.0d3 * grid_gde%r(1)

        earth_gde%tau = 0.0d0
        call create_cvector(grid_gde, h0, EDGE)
        call create_cvector(grid_gde, e0, FACE)
        call sourceField1d_s2(earth_gde, lmax, coeff, period_val, grid_gde, h0, e0)
        refHr_loc  = h0%z(i0,j0,kRr)
        refHph_loc = h0%x(i0,j0,kRs)
        refHth_loc = h0%y(i0,j0,kRs)
        call deall_cvector(h0)
        call deall_cvector(e0)

        write(*,'(a)') '  tau[S]  Q  |Hr/Hr0| dPh_Hr  |Hth/Hth0| dPh_Hth  |Hph/Hph0| dPh_Hph'
        do kk = 1,ntau2
            earth_gde%tau = tau_report(kk)
            call run_one(earth_gde, grid_gde, refHr_loc, refHth_loc, refHph_loc, omega_val, R0_EARTH)
        end do

        deallocate(earth_gde%layer, earth_gde%sigma)
        call deall_grid(grid_gde)
    end subroutine tau_report_earth_GDE

    subroutine run_one(earth_arg, grid_arg, refHr_val, refHth_val, refHph_val, omega_val, r0_val)
        type(conf1d_t), intent(in) :: earth_arg
        type(grid_t), intent(in)   :: grid_arg
        complex(8), intent(in)     :: refHr_val, refHth_val, refHph_val
        real(8), intent(in)        :: omega_val, r0_val
        type(cvector)               :: hloc, eloc
        real(8)                     :: Qloc, dphHr, dphHth, dphHph, magHr, magHth, magHph
        complex(8)                  :: cHr, cHth, cHph

        Qloc = omega_val * MU_0 * earth_arg%tau * r0_val

        call create_cvector(grid_arg, hloc, EDGE)
        call create_cvector(grid_arg, eloc, FACE)
        call sourceField1d_s2(earth_arg, lmax, coeff, 2.0d0*pi/omega_val, grid_arg, hloc, eloc)

        cHr  = hloc%z(i0,j0,kRr)
        cHph = hloc%x(i0,j0,kRs)
        cHth = hloc%y(i0,j0,kRs)

        magHr  = abs(cHr) / abs(refHr_val)
        magHth = abs(cHth) / abs(refHth_val)
        magHph = abs(cHph) / abs(refHph_val)

        dphHr  = wrap_deg( (atan2(aimag(cHr), real(cHr))   - atan2(aimag(refHr_val), real(refHr_val)))   * RAD2DEG )
        dphHth = wrap_deg( (atan2(aimag(cHth),real(cHth))  - atan2(aimag(refHth_val),real(refHth_val)))  * RAD2DEG )
        dphHph = wrap_deg( (atan2(aimag(cHph),real(cHph))  - atan2(aimag(refHph_val),real(refHph_val)))  * RAD2DEG )

        write(*,'(2x,es11.4,2x,es11.4,2x,f10.6,2x,f15.8,2x,f11.6,2x,f16.8,2x,f11.6,2x,f16.8)') &
            earth_arg%tau, Qloc, magHr, dphHr, magHth, dphHth, magHph, dphHph

        call deall_cvector(hloc)
        call deall_cvector(eloc)
    end subroutine run_one

end program test_tau_sweep
