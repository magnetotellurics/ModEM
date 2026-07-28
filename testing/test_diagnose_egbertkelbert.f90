program test_diagnose_egbertkelbert
! ****************************************************************************
! DIAGNOSTIC (not a permanent regression test): compares (a) the raw native
! E-field with the OLD-STYLE manual MT relabeling (Ex=-Etheta, Ey=+Ephi --
! this is what test_vs_modem_1D_impedance.f90 already validated against ModEM via
! Zxy/Zyx, <1% agreement) against (b) the SAME physical quantity computed by
! running the FULL new pipeline: rescale_source_coeffs(...,EGBERTKELBERT2012)
! BEFORE the solve, then apply_output_convention(...,EGBERTKELBERT2012) AFTER
! -- exactly what FWD1D.f90 now does. If (a) and (b) disagree, that pinpoints
! whether the bug is in rescale_source_coeffs, apply_output_convention, or
! their interaction.
! ****************************************************************************

    use field1d, only: conf1d_t, sourceField1d
    use output_convention
    use griddef
    use sg_vector
    implicit none

    type(conf1d_t)                 :: earth
    type(grid_t)                   :: grid
    type(cvector)                  :: h1d, e1d
    complex(8), allocatable        :: coeff(:), coeff_native(:)
    integer, parameter             :: lmax = 1
    integer, parameter             :: nx = 4, ny = 2, nz = 2
    integer                        :: istat, ii, jj
    real(8)                        :: period
    integer                        :: i_midph, i_nodeph, j_node, j_mid, kRr, kRs
    real(8), allocatable           :: depths_km(:), logrho(:)
    integer                        :: nL
    real(8), parameter             :: r0_m = 6.371d6
    type(output_convention_t)      :: native, target_conv
    complex(8) :: Ex_native, Ey_native, Ex_new, Ey_new

    period = 1000.0d0

    write(*,*) '=== test_diagnose_egbertkelbert: native-manual vs full-pipeline Ex/Ey ==='
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
    j_node   = 2
    j_mid    = 2
    kRs = 2
    kRr = 2

    call create_cvector(grid, h1d, EDGE)
    call create_cvector(grid, e1d, FACE)

    native = native_convention('S1')
    target_conv = get_convention('EGBERTKELBERT2012')

    allocate(coeff((lmax+1)**2), coeff_native((lmax+1)**2), STAT=istat)

    ! ================= Mode2 (l=1,m=+-1) =================
    coeff_native = dcmplx(0.0d0,0.0d0)
    coeff_native(3) = dcmplx(-0.5d0, 0.0d0)
    coeff_native(4) = dcmplx( 0.5d0, 0.0d0)

    ! (a) native, manual relabeling -- Ex = -Etheta = -E%y
    coeff = coeff_native
    call sourceField1d(earth, lmax, coeff, period, grid, h1d, e1d)
    Ex_native = -e1d%y(i_midph, j_node, kRr)

    ! (b) full pipeline -- rescale THEN solve THEN transform
    coeff = coeff_native
    call rescale_source_coeffs(coeff, lmax, target_conv)
    write(*,'(a,4f10.5)') 'Mode2 coeff before/after rescale (m=+1,-1): ', &
        dble(coeff_native(3)), dble(coeff_native(4)), dble(coeff(3)), dble(coeff(4))
    call sourceField1d(earth, lmax, coeff, period, grid, h1d, e1d)
    call apply_output_convention(h1d, e1d, native, target_conv)
    ! under EGBERTKELBERT2012 (theta reversed), the physical point at j_node
    ! (native) is now at index ny+2-j_node in the transformed array -- and
    ! %y should ALREADY represent Ex directly (no further manual relabeling)
    Ex_new = e1d%y(i_midph, (ny+2)-j_node, kRr)

    write(*,'(a,2es14.6,a,2es14.6)') 'Mode2 Ex_native(manual)= ', Ex_native, '   Ex_new(pipeline)= ', Ex_new
    write(*,'(a,2es14.6)') '  ratio Ex_new/Ex_native = ', Ex_new/Ex_native
    write(*,*)

    ! ================= Mode1 (l=1,m=0) =================
    coeff_native = dcmplx(0.0d0,0.0d0)
    coeff_native(2) = dcmplx(1.0d0, 0.0d0)

    coeff = coeff_native
    call sourceField1d(earth, lmax, coeff, period, grid, h1d, e1d)
    Ey_native = e1d%x(i_nodeph, j_mid, kRr)

    coeff = coeff_native
    call rescale_source_coeffs(coeff, lmax, target_conv)
    write(*,'(a,2f10.5)') 'Mode1 coeff before/after rescale (m=0): ', dble(coeff_native(2)), dble(coeff(2))
    call sourceField1d(earth, lmax, coeff, period, grid, h1d, e1d)
    call apply_output_convention(h1d, e1d, native, target_conv)
    Ey_new = e1d%x(i_nodeph, j_mid, kRr)   ! Ephi/%x untouched by theta transform's index range here (mid-theta)

    write(*,'(a,2es14.6,a,2es14.6)') 'Mode1 Ey_native(manual)= ', Ey_native, '   Ey_new(pipeline)= ', Ey_new
    write(*,'(a,2es14.6)') '  ratio Ey_new/Ey_native = ', Ey_new/Ey_native

    deallocate(coeff, coeff_native, depths_km, logrho, STAT=istat)
    deallocate(earth%layer, earth%sigma, STAT=istat)
    call deall_cvector(h1d)
    call deall_cvector(e1d)
    call deall_grid(grid)

contains

    subroutine read_layered_prm(fname, nL, depths_km, logrho)
        character(*), intent(in)                        :: fname
        integer, intent(out)                             :: nL
        real(8), allocatable, intent(out)                :: depths_km(:), logrho(:)
        integer, parameter  :: iun = 80
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

end program test_diagnose_egbertkelbert
