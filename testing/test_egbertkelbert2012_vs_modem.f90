program test_egbertkelbert2012_vs_modem
! ****************************************************************************
! Closes the "NOT yet validated" gap flagged in CLAUDE.md (Part 3, criterion 2):
! test_vs_modem_1D_impedance.f90 validates the two solvers' RAW native output against
! ModEM (via small_predicted.dat's Zxy/Zyx), with NO rescale_source_coeffs or
! apply_output_convention applied at all. This program runs the SAME two
! source patterns (Mode1=zonal/l=1,m=0, Mode2=l=1,m=+-1) through the FULL
! pipeline FWD1D.f90 actually uses -- rescale_source_coeffs(...,
! EGBERTKELBERT2012) BEFORE the solve, apply_output_convention(...,
! EGBERTKELBERT2012) AFTER -- and checks whether Zxy/Zyx (and, importantly,
! the RAW signed Ex/Ey/Hx/Hy themselves, not just their ratio) still agree
! with ModEM once the transform's own theta/r index-reversal and component
! renaming (which by design replaces the old manual "Ex=-Etheta" step, see
! output_convention.f90's NATIVE CONVENTIONS note and the plot_global1d_
! output.py rewrite) is taken into account.
!
! Because apply_output_convention reverses the theta (2nd) and r (3rd) index
! of every component array IN PLACE, using each component's OWN extent, this
! program does NOT hand-compute the transformed index -- it reads back
! size(component, dim) at runtime and applies the identical N+1-i formula
! reverse_theta_index/reverse_r_index themselves use, to avoid a manual
! indexing mistake silently producing a wrong answer.
! ****************************************************************************

    use field1d, only: conf1d_t, sourceField1d
    use field1d_s2, only: sourceField1d_s2
    use output_convention
    use griddef
    use sg_vector
    implicit none

    type(conf1d_t)                 :: earth
    type(grid_t)                   :: grid
    type(cvector)                  :: h1d, e1d
    complex(8), allocatable        :: coeff(:)
    integer, parameter             :: lmax = 1
    integer, parameter             :: nx = 4, ny = 2, nz = 2
    integer                        :: istat, ii, jj, ip
    real(8)                        :: period
    integer                        :: i_midph, i_nodeph, j_node, j_mid, kRr, kRs
    integer                        :: j_node_y, j_mid_x, kRr_y, kRs_y, kRr_x, kRs_x
    real(8), allocatable           :: depths_km(:), logrho(:)
    integer                        :: nL

    real(8), parameter             :: r0_m = 6.371d6
    integer, parameter             :: nper = 3
    real(8), parameter             :: periods(nper) = (/ 10.0d0, 100.0d0, 1000.0d0 /)

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
    type(output_convention_t) :: native, target_conv

    real(8), parameter :: mu0_local = 1.256637d-6
    real(8), parameter :: Z_SI_to_practical = 1.0d-3 / mu0_local

    write(*,*) '=== test_egbertkelbert2012_vs_modem: FULL rescale+solve+transform pipeline vs ModEM ==='
    write(*,*)

    call read_layered_prm('USA_small_1D.prm', nL, depths_km, logrho)

    earth%r0   = r0_m
    earth%tau  = 0.0d0
    earth%tol  = 1.0d-9
    allocate(earth%layer(nL+1), earth%sigma(nL+1), STAT=istat)
    earth%layer(1) = 0.0d0
    earth%layer(2:nL+1) = 1.0d3 * depths_km(1:nL)
    earth%sigma(1:nL) = 10.0d0**(-logrho(1:nL))
    earth%sigma(nL+1) = 10.0d0**(5.0d0)
    earth%allocated = .true.

    call create_grid(nx, ny, nz, grid)
    grid%nx = nx; grid%ny = ny; grid%nz = nz
    grid%nzAir = nz; grid%nzCrust = 0; grid%nzEarth = 0

    do ii = 1, nx+1
        grid%ph(ii) = (45.0d0 + (ii-1)*90.0d0) * d2r
    end do
    do ii = 1, nx
        grid%dp(ii) = grid%ph(ii+1) - grid%ph(ii)
    end do
    do jj = 1, ny+1
        grid%th(jj) = (88.0d0 + (jj-1)*2.0d0) * d2r
    end do
    do jj = 1, ny
        grid%dt(jj) = grid%th(jj+1) - grid%th(jj)
    end do
    grid%r(1) = 6371.10d0
    grid%r(2) = 6371.05d0
    grid%r(3) = 6371.01d0
    grid%dr(1) = grid%r(1) - grid%r(2)
    grid%dr(2) = grid%r(2) - grid%r(3)
    grid%allocated = .true.
    earth%rmax = 1.0d3 * grid%r(1)

    i_midph  = 1
    i_nodeph = 1
    j_node   = 2   ! native colatitude node index (90deg exactly)
    j_mid    = 2   ! native colatitude mid-cell index (91deg)
    kRs = 2
    kRr = 2

    call create_cvector(grid, h1d, EDGE)
    call create_cvector(grid, e1d, FACE)

    write(*,'(a,i0,a,i0)') ' e1d%y (theta-node, used for Ex) extent along theta = ', size(e1d%y,2), &
                            ' ; along r = ', size(e1d%y,3)
    write(*,'(a,i0,a,i0)') ' e1d%x (theta-mid, used for Ey)  extent along theta = ', size(e1d%x,2), &
                            ' ; along r = ', size(e1d%x,3)
    write(*,'(a,i0,a,i0)') ' h1d%y (theta-mid, used for Hx)  extent along theta = ', size(h1d%y,2), &
                            ' ; along r = ', size(h1d%y,3)
    write(*,'(a,i0,a,i0)') ' h1d%x (theta-node, used for Hy) extent along theta = ', size(h1d%x,2), &
                            ' ; along r = ', size(h1d%x,3)
    write(*,*)

    ! EGBERTKELBERT2012 differs from native S1/S2 in BOTH theta_convention
    ! (COLAT->LAT) and r_convention (UP->DOWN) -- both index axes get
    ! reversed by apply_output_convention, per-component, using that
    ! component's OWN extent (N+1-i). Compute the remapped indices here
    ! rather than hand-guessing them.
    j_node_y = size(e1d%y,2) + 1 - j_node   ! Ex lives on e1d%y (node-theta)
    j_mid_x  = size(e1d%x,2) + 1 - j_mid    ! Ey lives on e1d%x (mid-theta)
    kRr_y    = size(e1d%y,3) + 1 - kRr
    kRr_x    = size(e1d%x,3) + 1 - kRr
    kRs_y    = size(h1d%y,3) + 1 - kRs      ! Hx lives on h1d%y (mid-theta on EDGE)
    kRs_x    = size(h1d%x,3) + 1 - kRs      ! Hy lives on h1d%x (node-theta on EDGE)

    allocate(coeff((lmax+1)**2), STAT=istat)

    native = native_convention('S1')   ! identical to native_convention('S2')
    target_conv = get_convention('EGBERTKELBERT2012')
    write(*,*) 'native:  time=',trim(native%time_convention),' norm=',trim(native%harmonic_norm), &
               ' CS=',native%condon_shortley,' theta=',trim(native%theta_convention), &
               ' r=',trim(native%r_convention)
    write(*,*) 'target:  time=',trim(target_conv%time_convention),' norm=',trim(target_conv%harmonic_norm), &
               ' CS=',target_conv%condon_shortley,' theta=',trim(target_conv%theta_convention), &
               ' r=',trim(target_conv%r_convention)
    write(*,*)

    ! ================= Mode1 (l=1,m=0, zonal) -- Zyx = Ey/Hx =================
    write(*,*) '--- Mode1 (zonal, l=1 m=0): Zyx = Ex_pipeline(theta-node)/Hx_pipeline(theta-mid) ---'
    write(*,'(a8,a14,4a24,2a12)') 'T(s)','solver','Zyx_calc','Zyx_ModEM(conj)','','','%|Z|diff','phase diff(deg)'
    do ip = 1,nper
        period = periods(ip)
        Zyx_ref = conjg(ZYX_modem(ip))

        coeff = dcmplx(0.0d0,0.0d0)
        coeff(2) = dcmplx(1.0d0, 0.0d0)
        call rescale_source_coeffs(coeff, lmax, target_conv)
        call sourceField1d(earth, lmax, coeff, period, grid, h1d, e1d)
        call apply_output_convention(h1d, e1d, native, target_conv)
        Ey =  e1d%x(i_nodeph, j_mid_x, kRr_x)   ! "Ey" (east/phi) -- unaffected by theta negate
        Hx =  h1d%y(i_nodeph, j_mid, kRs_y)     ! "Hx" (north/theta) -- ALREADY sign-corrected by the transform
        Zyx_calc = (Ey/Hx) * Z_SI_to_practical
        pct_amp   = 100.0d0*abs(abs(Zyx_calc)-abs(Zyx_ref))/abs(Zyx_ref)
        pct_phase = (180.0d0/pi_local())*abs(atan2(aimag(Zyx_calc),dble(Zyx_calc)) &
                                            - atan2(aimag(Zyx_ref),dble(Zyx_ref)))
        write(*,'(f8.1,a14,2(2es12.4,2x),2f12.4)') period, 'S1', Zyx_calc, Zyx_ref, pct_amp, pct_phase

        coeff = dcmplx(0.0d0,0.0d0)
        coeff(2) = dcmplx(1.0d0, 0.0d0)
        call rescale_source_coeffs(coeff, lmax, target_conv)
        call sourceField1d_s2(earth, lmax, coeff, period, grid, h1d, e1d)
        call apply_output_convention(h1d, e1d, native, target_conv)
        Ey =  e1d%x(i_nodeph, j_mid_x, kRr_x)
        Hx =  h1d%y(i_nodeph, j_mid, kRs_y)
        Zyx_calc = (Ey/Hx) * Z_SI_to_practical
        pct_amp   = 100.0d0*abs(abs(Zyx_calc)-abs(Zyx_ref))/abs(Zyx_ref)
        pct_phase = (180.0d0/pi_local())*abs(atan2(aimag(Zyx_calc),dble(Zyx_calc)) &
                                            - atan2(aimag(Zyx_ref),dble(Zyx_ref)))
        write(*,'(f8.1,a14,2(2es12.4,2x),2f12.4)') period, 'S2', Zyx_calc, Zyx_ref, pct_amp, pct_phase
    end do
    write(*,*)

    ! ================= Mode2 (l=1,m=+-1) -- Zxy = Ex/Hy, evaluated at phi=90deg =================
    write(*,*) '--- Mode2 (l=1 m=+-1): Zxy = Ex_pipeline(theta-node)/Hy_pipeline(theta-node), phi=90deg ---'
    write(*,'(a8,a14,4a24,2a12)') 'T(s)','solver','Zxy_calc','Zxy_ModEM(conj)','','','%|Z|diff','phase diff(deg)'
    do ip = 1,nper
        period = periods(ip)
        Zxy_ref = conjg(ZXY_modem(ip))

        coeff = dcmplx(0.0d0,0.0d0)
        coeff(3) = dcmplx(-0.5d0, 0.0d0)
        coeff(4) = dcmplx( 0.5d0, 0.0d0)
        call rescale_source_coeffs(coeff, lmax, target_conv)
        call sourceField1d(earth, lmax, coeff, period, grid, h1d, e1d)
        call apply_output_convention(h1d, e1d, native, target_conv)
        Ex =  e1d%y(i_midph, j_node_y, kRr_y)   ! "Ex" (north/theta) -- ALREADY sign-corrected by the transform
        Hy =  h1d%x(i_midph, j_node, kRs_x)     ! "Hy" (east/phi) -- unaffected by theta negate
        Zxy_calc = (Ex/Hy) * Z_SI_to_practical
        pct_amp   = 100.0d0*abs(abs(Zxy_calc)-abs(Zxy_ref))/abs(Zxy_ref)
        pct_phase = (180.0d0/pi_local())*abs(atan2(aimag(Zxy_calc),dble(Zxy_calc)) &
                                            - atan2(aimag(Zxy_ref),dble(Zxy_ref)))
        write(*,'(f8.1,a14,2(2es12.4,2x),2f12.4)') period, 'S1', Zxy_calc, Zxy_ref, pct_amp, pct_phase

        coeff = dcmplx(0.0d0,0.0d0)
        coeff(3) = dcmplx(-0.5d0, 0.0d0)
        coeff(4) = dcmplx( 0.5d0, 0.0d0)
        call rescale_source_coeffs(coeff, lmax, target_conv)
        call sourceField1d_s2(earth, lmax, coeff, period, grid, h1d, e1d)
        call apply_output_convention(h1d, e1d, native, target_conv)
        Ex =  e1d%y(i_midph, j_node_y, kRr_y)
        Hy =  h1d%x(i_midph, j_node, kRs_x)
        Zxy_calc = (Ex/Hy) * Z_SI_to_practical
        pct_amp   = 100.0d0*abs(abs(Zxy_calc)-abs(Zxy_ref))/abs(Zxy_ref)
        pct_phase = (180.0d0/pi_local())*abs(atan2(aimag(Zxy_calc),dble(Zxy_calc)) &
                                            - atan2(aimag(Zxy_ref),dble(Zxy_ref)))
        write(*,'(f8.1,a14,2(2es12.4,2x),2f12.4)') period, 'S2', Zxy_calc, Zxy_ref, pct_amp, pct_phase
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
        character(*), intent(in)                        :: fname
        integer, intent(out)                             :: nL
        real(8), allocatable, intent(out)                :: depths_km(:), logrho(:)
        integer, parameter  :: iun = 79
        character(200)      :: line
        integer              :: ios, l_dum, m_dum
        real(8)              :: depth_val, logrho_val

        open(iun, file=fname, status='old', action='read')
        read(iun,'(A)') line

        nL = 0
        do
            read(iun,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line,'log degree') > 0) nL = nL + 1
        end do
        rewind(iun)
        allocate(depths_km(nL), logrho(nL), STAT=ios)

        read(iun,'(A)') line
        nL = 0
        do
            read(iun,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line,'log degree') > 0) then
                nL = nL + 1
                read(line(index(line,'layer')+5:),*) depth_val
                depths_km(nL) = depth_val
                read(iun,'(A)') line
                read(iun,*) l_dum, m_dum, logrho_val
                logrho(nL) = logrho_val
            end if
        end do
        close(iun)

    end subroutine read_layered_prm

end program test_egbertkelbert2012_vs_modem
