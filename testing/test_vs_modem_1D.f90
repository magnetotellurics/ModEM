program test_vs_modem_1D
! ****************************************************************************
! Merges test_kelbert2014_vs_modem_1D.f90 and test_sunegbert2012_vs_modem_1D.f90
! into one driver (same pattern as test_unit_sphere.f90, which runs both
! solvers side by side against one reference) -- runs BOTH field1d.f90
! (KELBERT2014) and field1d_sunegbert2012.f90 (SUNEGBERT2012) on the identical
! earth model/grid/coefficients and checks both against the same independent
! ModEM 3D-code benchmark (testing/small_predicted.dat, computed on
! testing/USA_small_1D.rho/.prm's 47-layer model with traditional ModEM
! sources in spherical coordinates -- imperfect since it doesn't account for
! Earth's sphericity in its BCs, but a good-enough benchmark per the user).
!
! Both solvers are now expected to pass CLEANLY, with no compensating sign of
! any kind: field1d.f90's three T(r)-valued formulas (Hr, Etheta, Ephi) were
! fixed (2026-07-26) to carry the same explicit sign as Sun & Egbert (2012)
! eq (5)-(6) -- see CLAUDE.md and the NOTE above sourceField1d in field1d.f90
! -- after which field1d.f90 and field1d_sunegbert2012.f90 agree exactly (up
! to a real normalization constant) in testing/test_kelbert2014_vs_sunegbert2012_l1m0.f90,
! and field1d_sunegbert2012.f90 was already independently confirmed to match
! ModEM directly with no sign correction needed at all.
!
! Two source patterns checked, both from this project's existing
! MTsource.1000sec.Mode{1,2}.prm coefficient values (hardcoded here directly,
! matching those files exactly -- l=1 only):
!   Mode2 (l=1,m=0, zonal): coefl(m=0)=1, coefl(m=+-1)=0 -- Zyx = Ey/Hx.
!   Mode1 (l=1,m=+-1): coefl(m=+1)=-0.5, coefl(m=-1)=+0.5 -- Zxy = Ex/Hy,
!     evaluated at phi=90deg (Ex,Hy are ~sin(phi), zero at phi=0; the exact
!     phi chosen doesn't affect the ratio -- verified analytically, see
!     CLAUDE.md).
!
! Units: field1d.f90/field1d_sunegbert2012.f90 output raw SI (E in V/m, H in
! A/m), giving Z=E/H in Ohms; small_predicted.dat reports Z in practical
! [mV/km]/[nT] (uses B in nT, not H in A/m): Z[mV/km/nT] = (1e-3/mu0) *
! Z[Ohms].
!
! Time convention: ModEM's INTERNAL convention is e^{-i*omega*t}, same as
! both solvers here -- raw fields need NO conjugation. Only
! small_predicted.dat's own impedance OUTPUT is exported/labeled
! e^{+i*omega*t} (its header says so explicitly), so it is the reference Z
! values (hardcoded below) that get conjugated once before comparison.
! ****************************************************************************

    use field1d, only: conf1d_t, sourceField1d
    use field1d_sunegbert2012, only: sourceField1d_sunegbert2012
    use griddef
    use sg_vector
    implicit none

    type(conf1d_t)                 :: earth
    type(grid_t)                   :: grid
    type(cvector)                  :: h1d, e1d
    complex(8), allocatable        :: coeff(:)
    integer, parameter             :: lmax = 1
    integer, parameter             :: nx = 4, ny = 2, nz = 2  ! nz=2 matches the THREE r-values (r(1:3)/dr(1:2)) set below
    integer                        :: istat, ii, jj, ip
    real(8)                        :: period
    integer                        :: i_midph, i_nodeph, j_node, j_mid, kRr, kRs
    real(8), allocatable           :: depths_km(:), logrho(:)
    integer                        :: nL

    real(8), parameter             :: r0_m = 6.371d6
    integer, parameter             :: nper = 3
    real(8), parameter             :: periods(nper) = (/ 10.0d0, 100.0d0, 1000.0d0 /)

    ! ModEM reference (testing/small_predicted.dat, e^{+iwt} per its own header),
    ! transcribed directly, station S001, T=10,100,1000s in that order.
    complex(8), parameter :: ZXY_modem(nper) = (/ &
        dcmplx(2.337319d0, 0.9753857d0), &
        dcmplx(1.641010d0, 0.6809239d0), &
        dcmplx(0.6345598d0, 0.5235962d0) /)
    complex(8), parameter :: ZYX_modem(nper) = (/ &
        dcmplx(-2.337399d0, -0.9754012d0), &
        dcmplx(-1.641054d0, -0.6809613d0), &
        dcmplx(-0.6345333d0, -0.5236308d0) /)

    complex(8) :: Ex, Ey, Hx, Hy, Zxy_calc, Zyx_calc, Zxy_ref, Zyx_ref
    real(8) :: pct_amp, pct_phase

    real(8), parameter :: mu0_local = 1.256637d-6   ! matches field1d.f90's own literal
    real(8), parameter :: Z_SI_to_practical = 1.0d-3 / mu0_local

    write(*,*) '=== test_vs_modem_1D: field1d.f90 AND field1d_sunegbert2012.f90 vs ModEM small_predicted.dat ==='
    write(*,*)

    ! ---- earth model: USA_small_1D.prm's own 47 layers + FWD1D.f90's hardcoded core ----
    call read_layered_prm('USA_small_1D.prm', nL, depths_km, logrho)
    write(*,'(a,i3,a)') ' Read ', nL, ' layers from USA_small_1D.prm'

    earth%r0   = r0_m
    earth%tau  = 0.0d0     ! no extra thin-sheet conductance -- comparing against a plain ModEM 1D benchmark
    earth%tol  = 1.0d-9
    allocate(earth%layer(nL+1), earth%sigma(nL+1), STAT=istat)
    earth%layer(1) = 0.0d0
    earth%layer(2:nL+1) = 1.0d3 * depths_km(1:nL)
    earth%sigma(1:nL) = 10.0d0**(-logrho(1:nL))
    earth%sigma(nL+1) = 10.0d0**(5.0d0)  ! core conductivity, matching FWD1D.f90's own hardcoded value
    earth%allocated = .true.

    ! ---- hand-built grid: phi mid-cell at 90deg (Mode1's Zxy needs this --
    ! see header note), theta node at 90deg exactly (equator), thin shell
    ! just above the surface r0 ----
    call create_grid(nx, ny, nz, grid)
    grid%nx = nx; grid%ny = ny; grid%nz = nz
    grid%nzAir = nz; grid%nzCrust = 0; grid%nzEarth = 0

    do ii = 1, nx+1
        grid%ph(ii) = (45.0d0 + (ii-1)*90.0d0) * d2r   ! nodes 45,135,225,315,405 -> mid(i=1)=90deg
    end do
    do ii = 1, nx
        grid%dp(ii) = grid%ph(ii+1) - grid%ph(ii)
    end do

    do jj = 1, ny+1
        grid%th(jj) = (88.0d0 + (jj-1)*2.0d0) * d2r    ! nodes 88,90,92 -> node(j=2)=90deg exactly
    end do
    do jj = 1, ny
        grid%dt(jj) = grid%th(jj+1) - grid%th(jj)
    end do

    grid%r(1) = 6371.10d0
    grid%r(2) = 6371.05d0
    grid%r(3) = 6371.01d0   ! ~29m above the true surface r0=6371.000km after Rr's half-cell offset below
    grid%dr(1) = grid%r(1) - grid%r(2)
    grid%dr(2) = grid%r(2) - grid%r(3)
    grid%allocated = .true.

    earth%rmax = 1.0d3 * grid%r(1)

    i_midph  = 1   ! mid-cell(i=1) = 90deg  (H%x, E%y -- "Yp-based", node-theta)
    i_nodeph = 1   ! node(i=1)     = 45deg  (H%y, E%x -- "Yt-based", mid-theta; phi-independent for Mode2)
    j_node   = 2   ! node(j=2)     = 90deg exactly (H%x, E%y)
    j_mid    = 2   ! mid-cell(j=2) = 91deg  (H%y, E%x; 1deg off equator)
    kRs = 2   ! Rs(2) = grid%r(2)*1e3
    kRr = 2   ! Rr(2) = Rs(2) - dr(2)/2*1e3, ~29m above r0

    call create_cvector(grid, h1d, EDGE)
    call create_cvector(grid, e1d, FACE)

    write(*,'(a,f10.4,a)') ' Rs(2) = ', 1.0d3*grid%r(2), ' m'
    write(*,'(a,f10.4,a)') ' Rr(2) = ', 1.0d3*grid%r(2) - 1.0d3*grid%dr(2)/2.0d0, ' m  (r0 = 6371000 m)'
    write(*,*)

    allocate(coeff((lmax+1)**2), STAT=istat)

    ! ================= Mode2 (l=1,m=0, zonal) -- Zyx = Ey/Hx =================
    coeff = dcmplx(0.0d0,0.0d0)
    coeff(2) = dcmplx(1.0d0, 0.0d0)   ! l=1,m=0 -- matches MTsource.1000sec.Mode2.prm

    write(*,*) '--- Mode2 (zonal, l=1 m=0): Zyx = Ey/Hx ---'
    write(*,'(a8,a14,4a24,2a12)') 'T(s)','solver','Zyx_calc','Zyx_ModEM(conj)','','','%|Z|diff','phase diff(deg)'
    do ip = 1,nper
        period = periods(ip)
        Zyx_ref = conjg(ZYX_modem(ip))         ! file is e^{+iwt}; solvers are e^{-iwt}

        call sourceField1d(earth, lmax, coeff, period, grid, h1d, e1d)
        Ey =  e1d%x(i_nodeph, j_mid, kRr)     ! Ey = +Ephi = +E%x
        Hx = -h1d%y(i_nodeph, j_mid, kRs)     ! Hx = -Htheta = -H%y
        Zyx_calc = (Ey/Hx) * Z_SI_to_practical
        pct_amp   = 100.0d0*abs(abs(Zyx_calc)-abs(Zyx_ref))/abs(Zyx_ref)
        pct_phase = (180.0d0/pi_local())*abs(atan2(aimag(Zyx_calc),dble(Zyx_calc)) &
                                            - atan2(aimag(Zyx_ref),dble(Zyx_ref)))
        write(*,'(f8.1,a14,2(2es12.4,2x),2f12.4)') period, 'KELBERT2014', Zyx_calc, Zyx_ref, pct_amp, pct_phase

        call sourceField1d_sunegbert2012(earth, lmax, coeff, period, grid, h1d, e1d)
        Ey =  e1d%x(i_nodeph, j_mid, kRr)
        Hx = -h1d%y(i_nodeph, j_mid, kRs)
        Zyx_calc = (Ey/Hx) * Z_SI_to_practical
        pct_amp   = 100.0d0*abs(abs(Zyx_calc)-abs(Zyx_ref))/abs(Zyx_ref)
        pct_phase = (180.0d0/pi_local())*abs(atan2(aimag(Zyx_calc),dble(Zyx_calc)) &
                                            - atan2(aimag(Zyx_ref),dble(Zyx_ref)))
        write(*,'(f8.1,a14,2(2es12.4,2x),2f12.4)') period, 'SUNEGBERT2012', Zyx_calc, Zyx_ref, pct_amp, pct_phase
    end do
    write(*,*)

    ! ================= Mode1 (l=1,m=+-1) -- Zxy = Ex/Hy, evaluated at phi=90deg =================
    coeff = dcmplx(0.0d0,0.0d0)
    coeff(3) = dcmplx(-0.5d0, 0.0d0)  ! l=1,m=+1 -- matches MTsource.1000sec.Mode1.prm
    coeff(4) = dcmplx( 0.5d0, 0.0d0)  ! l=1,m=-1

    write(*,*) '--- Mode1 (l=1 m=+-1): Zxy = Ex/Hy, evaluated at phi=90deg ---'
    write(*,'(a8,a14,4a24,2a12)') 'T(s)','solver','Zxy_calc','Zxy_ModEM(conj)','','','%|Z|diff','phase diff(deg)'
    do ip = 1,nper
        period = periods(ip)
        Zxy_ref = conjg(ZXY_modem(ip))

        call sourceField1d(earth, lmax, coeff, period, grid, h1d, e1d)
        Ex = -e1d%y(i_midph, j_node, kRr)     ! Ex = -Etheta = -E%y
        Hy =  h1d%x(i_midph, j_node, kRs)     ! Hy = +Hphi = +H%x
        Zxy_calc = (Ex/Hy) * Z_SI_to_practical
        pct_amp   = 100.0d0*abs(abs(Zxy_calc)-abs(Zxy_ref))/abs(Zxy_ref)
        pct_phase = (180.0d0/pi_local())*abs(atan2(aimag(Zxy_calc),dble(Zxy_calc)) &
                                            - atan2(aimag(Zxy_ref),dble(Zxy_ref)))
        write(*,'(f8.1,a14,2(2es12.4,2x),2f12.4)') period, 'KELBERT2014', Zxy_calc, Zxy_ref, pct_amp, pct_phase

        call sourceField1d_sunegbert2012(earth, lmax, coeff, period, grid, h1d, e1d)
        Ex = -e1d%y(i_midph, j_node, kRr)
        Hy =  h1d%x(i_midph, j_node, kRs)
        Zxy_calc = (Ex/Hy) * Z_SI_to_practical
        pct_amp   = 100.0d0*abs(abs(Zxy_calc)-abs(Zxy_ref))/abs(Zxy_ref)
        pct_phase = (180.0d0/pi_local())*abs(atan2(aimag(Zxy_calc),dble(Zxy_calc)) &
                                            - atan2(aimag(Zxy_ref),dble(Zxy_ref)))
        write(*,'(f8.1,a14,2(2es12.4,2x),2f12.4)') period, 'SUNEGBERT2012', Zxy_calc, Zxy_ref, pct_amp, pct_phase
    end do

    deallocate(coeff, depths_km, logrho, STAT=istat)
    deallocate(earth%layer, earth%sigma, STAT=istat)
    call deall_cvector(h1d)
    call deall_cvector(e1d)
    call deall_grid(grid)

contains

    real(8) function pi_local()
        pi_local = 3.14159265357898d0
    end function pi_local

    subroutine read_layered_prm(fname, nL, depths_km, logrho)
        ! Minimal, self-contained reader for field1d.f90's layered .prm format.
        character(*), intent(in)                        :: fname
        integer, intent(out)                             :: nL
        real(8), allocatable, intent(out)                :: depths_km(:), logrho(:)
        integer, parameter  :: iun = 79
        character(200)      :: line
        integer              :: ios, l_dum, m_dum
        real(8)              :: depth_val, logrho_val

        open(iun, file=fname, status='old', action='read')
        read(iun,'(A)') line   ! "Format: ..." header line

        ! first pass: count layers
        nL = 0
        do
            read(iun,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line,'log degree') > 0) nL = nL + 1
        end do
        rewind(iun)
        allocate(depths_km(nL), logrho(nL), STAT=ios)

        ! second pass: extract depth + logrho per layer
        read(iun,'(A)') line   ! skip "Format: ..." again
        nL = 0
        do
            read(iun,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line,'log degree') > 0) then
                nL = nL + 1
                read(line(index(line,'layer')+5:),*) depth_val
                depths_km(nL) = depth_val
                read(iun,'(A)') line          ! "l  m  value  min  max" header -- skip
                read(iun,*) l_dum, m_dum, logrho_val
                logrho(nL) = logrho_val
            end if
        end do
        close(iun)

    end subroutine read_layered_prm

end program test_vs_modem_1D
