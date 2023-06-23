!
!> Derived class to define a CSEM Source with E_p computed using EM1D
!
module SourceCSEM_EM1D
    !
    use SourceCSEM
    use TransmitterCSEM
    use TransmitterArray
    use ModelParameterCell_SG
    use ModelParameterCell_SG_VTI
    !
    type, extends( SourceCSEM_t ) :: SourceCSEM_EM1D_t
        !
        integer :: i_tx
        !
        real( kind=prec ) :: location(3)
        !
        class( Vector_t ), allocatable :: cond_anomaly_h, cond_anomaly_v
        !
        contains
            !
            final :: SourceCSEM_EM1D_dtor
            !
            procedure, public :: createE => createE_SourceCSEM_EM1D
            !
            procedure, public :: createRHS => createRHS_SourceCSEM_EM1D
            !
            procedure, public :: set1DModel => set1DModel_SourceCSEM_EM1D
            !
            procedure, private :: setTemp_SourceCSEM_EM1D, create_Ep_from_EM1D
            !
            procedure, private :: create_background_data, create_source_data
            !
    end type SourceCSEM_EM1D_T
    !
    integer :: nlay1D_temp
    real( kind=prec ), dimension(:), allocatable, private :: zlay1D_temp     !> (m)   Depth to top of each layer, first layer ignored 
    real( kind=prec ), dimension(:), allocatable, private :: sig1D_temp_h    !> (S/m) Layer conductivities in the horizontal direction used in VTI
    real( kind=prec ), dimension(:), allocatable, private :: sig1D_temp_v    !> (S/m) Layer conductivities in the vertical direction used in VTI   
    !
    interface SourceCSEM_EM1D_t
        module procedure SourceCSEM_EM1D_ctor
    end interface SourceCSEM_EM1D_t
    !
contains
    !
    !> SourceCSEM_EM1D constructor
    function SourceCSEM_EM1D_ctor( model_operator, sigma, period, location, i_tx ) result( self )
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        class( ModelParameter_t ), target, intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: period, location(3)
        integer, intent( in ) :: i_tx
        !
        type( SourceCSEM_EM1D_t ) :: self
        !
        !write( *, * ) "Constructor SourceCSEM_EM1D_t"
        !
        call self%init
        !
        self%model_operator => model_operator
        !
        self%sigma => sigma
        !
        self%period = period
        !
        self%location = location
        !
        self%i_tx = i_tx
        !
        self%non_zero_source = .TRUE.
        !
        self%non_zero_bc = .FALSE.
        !
    end function SourceCSEM_EM1D_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    subroutine SourceCSEM_EM1D_dtor( self )
        implicit none
        !
        type( SourceCSEM_EM1D_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor SourceCSEM_EM1D_t"
        !
        call self%dealloc
        !
    end subroutine SourceCSEM_EM1D_dtor
    !
    !> Set self%E from forward modeling 1D
    !
    subroutine createE_SourceCSEM_EM1D( self )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( inout ) :: self
        !
        type( backgrounddata ) :: bgdat    !cond_cell description, coordinates, data
        type( sorec ) :: src               !source specification
        type( freqdata ) :: freqdat        !frequency dependent specifications
        type( refl_struct ) :: refl_var    !all variables that have to be remembered while computing 1D fields
        !
        integer :: ifreq, icur, comm
        complex( kind=prec ) :: i_omega_mu
        !
        !> Verbose...
        write( *, * ) "          - Extract CSEM Source from EM1D"
        !
        call self%set1DModel
        !call self%set1DModel( self%location(1), self%location(2) )
        !
        bgdat%omega = ( 2.0 * PI / self%period )
        bgdat%dowhat = 1
        !
        call self%create_background_data( bgdat )
        !
        call self%create_source_data( src, freqdat )
        !
        ifreq = 1
        icur = 1
        comm = 1
        refl_var%nzrecHxy = 0
        refl_var%nzrecHz = 0
        call reflectivity_unified( src, bgdat, refl_var, ifreq, icur, comm ) !> Output field will be saved in bgdat
        !
        call self%create_Ep_from_EM1D( self%sigma%metric%grid, bgdat )       !> Put the 1D field into E_p
        !
        i_omega_mu = cmplx( 0., real( -1.0d0 * isign * mu_0 * ( 2.0 * PI / self%period ), kind=prec ), kind=prec )
        !
        !> Construct E from E_p
        allocate( self%E(1), source = E_p )
        !
        select type( E => self%E(1) )
            !
            class is( cVector3D_SG_t )
                !
                select type( E_P )
                    !
                    class is( cVector3D_SG_t )
                        !
                        E%x = self%cond_anomaly_h%getAxis("x") * E_P%x
                        E%y = self%cond_anomaly_h%getAxis("y") * E_P%y
                        E%z = self%cond_anomaly_v%getAxis("z") * E_P%z
                        !
                    class default
                        stop "createE_SourceCSEM_EM1D > Unclassified E_P"
                end select
                !
            class default
                stop "createE_SourceCSEM_EM1D > Unclassified E"
        end select
        !
        call self%E(1)%mult( i_omega_mu )
        !
        deallocate( bgdat%sigv, bgdat%sigh )
        deallocate( bgdat%epsrv, bgdat%epsrh )
        deallocate( bgdat%zbound )
        deallocate( bgdat%Expos, bgdat%Eypos, bgdat%Ezpos )
        deallocate( bgdat%Ex, bgdat%Ey, bgdat%Ez  )
        !
        deallocate( src%nelem  )
        deallocate( src%pos  )
        deallocate( src%ljx, src%ljy, src%ljz )
        deallocate( src%akx, src%aky, src%akz )
        !
        call self%createRHS
        !
    end subroutine createE_SourceCSEM_EM1D
    !
    !> No subroutine briefing
    !
    subroutine create_Ep_from_EM1D( self, grid, bgdat )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( inout ) :: self
        class( Grid_t ), intent( in ) :: grid
        type( backgrounddata ), intent( in ) :: bgdat
        !
        integer ix, iy, iz, counter
        !
        if( allocated( E_p ) ) deallocate( E_p )
        allocate( E_p, source = cVector3D_SG_t( grid, EDGE ) )
        !
        !> Fill e_vector (cVector3D_SG) from E2D (Esoln2DTM_t)
        select type( E_P )
            !
            class is( cVector3D_SG_t )
                !
                counter = 1
                !> E-field corresponding to these nodes is Ex
                do iz = 1,grid%Nz+1    !Edge Z
                    do iy = 1,grid%Ny+1     !Edge Y
                        do ix = 1,grid%Nx       !Center X
                            E_p%x(ix,iy,iz) = bgdat%Ex(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
                counter = 1
                !> E-field corresponing to these nodes is Ey
                do iz = 1, grid%Nz+1    !Edge Z
                    do iy = 1, grid%Ny      !Center y
                        do ix = 1, grid%Nx+1    !Edge x
                            E_p%y(ix,iy,iz) = bgdat%Ey(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
                counter = 1
                !> E-field corresponing to these nodes is Ez
                do iz = 1,grid%Nz !Center Z
                    do iy = 1,grid%Ny+1 !Edge y
                        do ix = 1,grid%Nx+1 !Edge x
                            E_p%z(ix,iy,iz) = bgdat%Ez(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
            class default
                stop "create_Ep_from_EM1D > Unclassified E_P"
        end select
        !
    end subroutine create_Ep_from_EM1D
    !
    !> Set RHS from self%E
    !
    subroutine createRHS_SourceCSEM_EM1D( self )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( inout ) :: self
        !
        if( allocated( self%rhs ) ) deallocate( self%rhs )
        allocate( cVector3D_SG_t :: self%rhs(1) )
        !
        self%rhs(1) = self%E(1)
        !
        call self%rhs(1)%mult( self%model_operator%metric%Vedge )
        !
    end subroutine createRHS_SourceCSEM_EM1D
    !
    !> No subroutine briefing
    !
    subroutine setTemp_SourceCSEM_EM1D( self, sigma_cell_h, sigma_cell_v )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: sigma_cell_h, sigma_cell_v
        !
        integer :: i, j, k, nzAir
        real( kind=prec ) :: wt, temp_sigma_value
        complex( kind=prec ), allocatable :: v(:, :, :)
        !
        ! for layer boundaries use z-edges of 3D grid
        if( allocated( zlay1D_temp ) ) then
            deallocate( zlay1D_temp, sig1D_temp_h, sig1D_temp_v )
        end if
        !
        nlay1D_temp = sigma_cell_h%grid%Nz
        !
        allocate( zlay1D_temp( nlay1D_temp ) )
        allocate( sig1D_temp_h( nlay1D_temp ) )
        allocate( sig1D_temp_v( nlay1D_temp ) )
        !
        do i = 1, nlay1D_temp
            zlay1D_temp(i) = sigma_cell_h%grid%zEdge(i)
        end do
        !
        ! For create sig1D, we divide this process into two parts (1) for air layers and 
        !    (2) for earth layers
        ! For air layer, sig1D equal to air layer conductivity
        ! For earth layer, The Geometric mean is be used to create sig1D
        !
        nzAir = sigma_cell_h%grid%nzAir
        !
        nlay1D = sigma_cell_h%nz + nzAir
        !
        v = sigma_cell_h%getV()
        sig1D_temp_h(1:nzAir) = SIGMA_AIR !v(1,1,1:nzAir)
        !
        v = sigma_cell_v%getV()
        sig1D_temp_v(1:nzAir) = SIGMA_AIR !v(1,1,1:nzAir)
        !
        if( trim( get_1D_from ) == "Geometric_mean" ) then
            !
            v = sigma_cell_h%getV()
            !
            do k = nzAir+1, nlay1D_temp
                !
                wt = R_ZERO
                temp_sigma_value = R_ZERO
                !
                do i = 1, sigma_cell_h%grid%Nx
                    do j = 1, sigma_cell_h%grid%Ny
                        !
                        wt = wt + sigma_cell_h%grid%dx(i) * sigma_cell_h%grid%dy(j)
                        !
                        temp_sigma_value = temp_sigma_value + v(i,j,k-nzAir) * &
                        sigma_cell_h%grid%dx(i) * sigma_cell_h%grid%dy(j)
                        !
                    enddo
                enddo
                !
                sig1D_temp_h(k) = exp( temp_sigma_value / wt )
                !
            enddo
            !
            v = sigma_cell_v%getV()
            !
            do k = nzAir+1, nlay1D_temp
                !
                wt = R_ZERO
                temp_sigma_value = R_ZERO
                !
                do i = 1, sigma_cell_v%grid%Nx
                    do j = 1, sigma_cell_v%grid%Ny
                        !
                        wt = wt + sigma_cell_v%grid%dx(i) * sigma_cell_v%grid%dy(j)
                        !
                        temp_sigma_value = temp_sigma_value + v(i,j,k-nzAir) * &
                        sigma_cell_v%grid%dx(i) * sigma_cell_v%grid%dy(j)
                        !
                    enddo
                enddo
                !
                sig1D_temp_v(k) = exp( temp_sigma_value / wt )
                !
            enddo
            !
        elseif( trim( get_1D_from ) == "At_Tx_Position" ) then
            !
            stop "Error: setTemp_SourceCSEM_Dipole1D > At_Tx_Position not implemented yet"
            !
        elseif( trim(get_1d_from) == "Geometric_mean_around_Tx" ) then
            !
            stop "Error: setTemp_SourceCSEM_Dipole1D > Geometric_mean_around_Tx not implemented yet"
            !
        elseif( trim(get_1d_from) == "Full_Geometric_mean" ) then
            !
            stop "Error: setTemp_SourceCSEM_Dipole1D > Full_Geometric_mean not implemented yet"
            !
        elseif( trim( get_1d_from ) == "Fixed_Value" ) then
            !
            stop "Error: setTemp_SourceCSEM_Dipole1D > Fixed_Value not implemented yet"
            !
        else
            !
            stop "Error: setTemp_SourceCSEM_Dipole1D > Unknown get_1d_from"
            !
        endif
        !
    end subroutine setTemp_SourceCSEM_EM1D
    !
    !> this is a private routine, used to extract layer averages from
    !> a 3D conductivity parameter (sigma) and set up
    !> (1) nlay1D    !> Number of layers
    !> (2) sig1D => !> (S/m) Layer conductivities 
    !> (3) zlay1D => !> (m)   Depth to top of each layer, first layer ignored the 1D Model  z_P, sigma_P
    !
    subroutine set1DModel_SourceCSEM_EM1D( self )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( inout ) :: self
        !
        class( Scalar_t ), allocatable :: sigma_cell_h, sigma_cell_v
        !
        select type( sigma => self%sigma )
            !
            class is( ModelParameterCell_SG_t )
                !
                allocate( sigma_cell_h, source = sigma%cell_cond )
                !
                call self%setTemp_SourceCSEM_EM1D( sigma_cell_h, sigma_cell_h )
                !
                sig1D = sig1D_temp_h
                !
                call self%setCondAnomally( sigma_cell_h, self%cond_anomaly_h )
                !
                sig1D = sig1D_temp_v
                !
                call self%setCondAnomally( sigma_cell_h, self%cond_anomaly_v )
                !
                deallocate( sigma_cell_h )
                !
            class is( ModelParameterCell_SG_VTI_t )
                !
                allocate( sigma_cell_h, source = sigma%cell_cond_h )
                allocate( sigma_cell_v, source = sigma%cell_cond_v )
                !
                call self%setTemp_SourceCSEM_EM1D( sigma_cell_h, sigma_cell_v )
                !
                sig1D = sig1D_temp_h
                !
                call self%setCondAnomally( sigma_cell_h, self%cond_anomaly_h )
                !
                deallocate( sigma_cell_h )
                !
                sig1D = sig1D_temp_v
                !
                call self%setCondAnomally( sigma_cell_v, self%cond_anomaly_v )
                !
                deallocate( sigma_cell_v )
                !
            class default
                stop "Error: set1DModel_SourceCSEM_EM1D > Unclassified sigma"
            !
        end select
        !
    end subroutine set1DModel_SourceCSEM_EM1D
    !
    !> No Subroutine briefing
    !
    subroutine create_background_data( self, bgdat )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( in ) :: self
        type( backgrounddata ), intent( inout ) :: bgdat 
        !
        integer :: counter, ilay, ix, iy, iz, ierr
        integer( kind=int32 ) :: nx1, ny1, nz1    !nr of points in my domain for which fields are computed
        !
        bgdat%nlay= nlay1D_temp
        !allocate vectors for medium properties
        allocate(bgdat%sigv(bgdat%nlay),bgdat%sigh(bgdat%nlay),bgdat%epsrv(bgdat%nlay),bgdat%epsrh(bgdat%nlay), stat=ierr)
        !if(ierr.ne.0) call alloc_error(pid,'readinput','sig, epsr',ierr)
        !allocate vector for layer boundary depths: 1 element less than nr of layers
        allocate(bgdat%zbound(bgdat%nlay-1),stat=ierr)
        !if(ierr.ne.0) call alloc_error(pid,'readinput','zbound',ierr)
        !
        bgdat%rsplmin = 50.0
        bgdat%aniso = vti !default vti, change to iso if any layer is isotropic!!!> OR keep it VTI in all cases: if isotropic case, sig1D_temp_h==sig1D_temp_v 
        !
        !write(*,*) 'bgdat%nlay', bgdat%nlay
        do ilay = 1, bgdat%nlay
            bgdat%sigh(ilay) = sig1D_temp_h(ilay)
            bgdat%sigv(ilay) = sig1D_temp_v(ilay)
            bgdat%epsrh(ilay) = 1.0
            bgdat%epsrv(ilay) = 1.0
            !write(65,*) ilay, sig1D_temp_h(ilay), sig1D_temp_v(ilay)
        enddo 
        !
        do ilay = 1, bgdat%nlay - 1
            bgdat%zbound(ilay) = zlay1D_temp(ilay)
        enddo
        !
        nx1 = (self%sigma%metric%grid%Nx) * (self%sigma%metric%grid%Ny+1) * (self%sigma%metric%grid%Nz+1)
        ny1 = (self%sigma%metric%grid%Nx+1) * (self%sigma%metric%grid%Ny) * (self%sigma%metric%grid%Nz+1)
        nz1 = (self%sigma%metric%grid%Nx+1) * (self%sigma%metric%grid%Ny+1) * (self%sigma%metric%grid%Nz)
        !write(*,*)"nx1,ny1,nz1: ", nx1,ny1,nz1
        bgdat%nExy = 0
        bgdat%nEx = nx1
        bgdat%nEy = ny1
        bgdat%nEz = nz1
        !
        bgdat%nHxy = 0
        bgdat%nHx = 0
        bgdat%nHy = 0
        bgdat%nHz = 0
        !
        bgdat%allcomp_samecoord = .FALSE.
        bgdat%allzrec_samecoord = .TRUE.
        !
        allocate(bgdat%Expos(bgdat%nEx,3),bgdat%Eypos(bgdat%nEy,3),bgdat%Ezpos(bgdat%nEz,3), stat=ierr)
        if( ierr .NE. 0 ) call alloc_error(pid,'In-backgroundfield','Epos',ierr)
        !
        !> allocate(bgdat%Hxpos(bgdat%nHx,3),bgdat%Hypos(bgdat%nHy,3),bgdat%Hzpos(bgdat%nHz,3), stat=ierr)
        !> if(ierr.ne.0) call alloc_error(pid,'In-backgroundfield','Hpos',ierr)
        !
        ix = bgdat%nEx*1.5
        iy = bgdat%nEy*1.5
        iz = bgdat%nEz*1.5
        allocate(bgdat%Ex(ix),bgdat%Ey(iy),bgdat%Ez(iz), stat=ierr)
        if( ierr .NE. 0 ) call alloc_error(pid,'Out-backgroundfield','E fields',ierr)    
        !
        !> allocate(bgdat%Hx(nxyz),bgdat%Hy(nxyz),bgdat%Hz(nxyz), stat=ierr)
        !> if(ierr.ne.0) call alloc_error(pid,'Out-backgroundfield','H fields',ierr)  
        !
        bgdat%Expos=0
        bgdat%Eypos=0
        bgdat%Ezpos=0
        !
        bgdat%Ex = 0._real64
        bgdat%Ey = 0._real64
        bgdat%Ez = 0._real64
        !
        counter = 1
        !> E-field corresponding to these nodes is Ex
        do iz = 1, self%sigma%metric%grid%Nz+1 !Edge Z
            do iy = 1, self%sigma%metric%grid%Ny+1 !Edge Y
                do ix = 1,self%sigma%metric%grid%Nx !Center X
                    bgdat%Expos(counter,1) = self%sigma%metric%grid%xCenter(ix)
                    bgdat%Expos(counter,2) = self%sigma%metric%grid%yEdge(iy)
                    bgdat%Expos(counter,3) = self%sigma%metric%grid%zEdge(iz)
                    counter = counter + 1
                enddo
            enddo
        enddo
        !
        counter = 1
        !> E-field corresponing to these nodes is Ey
        do iz = 1,self%sigma%metric%grid%Nz+1 !Edge Z
            do iy = 1,self%sigma%metric%grid%Ny !Center y
                do ix = 1,self%sigma%metric%grid%Nx+1 !Edge x
                    bgdat%Eypos(counter,1) = self%sigma%metric%grid%xEdge(ix)
                    bgdat%Eypos(counter,2) = self%sigma%metric%grid%yCenter(iy)
                    bgdat%Eypos(counter,3) = self%sigma%metric%grid%zEdge(iz)
                    counter = counter + 1
                enddo
            enddo
        enddo
        !
        counter = 1
        !> E-field corresponing to these nodes is Ez
        do iz = 1,self%sigma%metric%grid%Nz !Center Z
            do iy = 1,self%sigma%metric%grid%Ny+1 !Edge y
                do ix = 1,self%sigma%metric%grid%Nx+1 !Edge x
                    bgdat%Ezpos(counter,1)= self%sigma%metric%grid%xEdge(ix)
                    bgdat%Ezpos(counter,2) = self%sigma%metric%grid%yEdge(iy)
                    bgdat%Ezpos(counter,3) = self%sigma%metric%grid%zCenter(iz)
                    counter = counter + 1
                enddo
            enddo
        enddo
        !
    end subroutine create_background_data
    !
    !> No Subroutine briefing
    !
    subroutine create_source_data( self, src, freqdat )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( in ) :: self
        type( sorec ), intent( inout ) :: src !source specification
        type( freqdata ), intent( inout ) :: freqdat !frequency dependent specifications
        !
        integer :: nelem, ierr
        !
        !> fill the Freq object with the required information
        freqdat%nfreq = 1
        !allocate frequency vector
        allocate( freqdat%omega(1), stat = ierr )
        !
        !> Instantiate Transmitter's Source - According to transmitter type and chosen via control file
        select type( Tx => transmitters( self%i_tx )%Tx )
            !
            class is( TransmitterCSEM_t )
                !
                freqdat%omega(1)= 2.0 * PI / Tx%period     !angular frequencies (2*pi*f, f in Hertz)
                !
                !> fill the sources object with the required information    
                src%type=1   !receiver = 0, dipole = 1, wire = 2, star = 3
                src%srcname="Tx"  !> dummy name
                !
                nelem=1
                allocate(src%nelem(nelem), stat=ierr)
                allocate(src%pos(3,nelem), stat=ierr)
                allocate(src%ljx(nelem), stat=ierr)
                allocate(src%ljy(nelem), stat=ierr)
                allocate(src%ljz(nelem), stat=ierr)
                !
                allocate(src%akx(nelem), stat=ierr)
                allocate(src%aky(nelem), stat=ierr)
                allocate(src%akz(nelem), stat=ierr)
                !
                src%nelem(nelem)= nelem
                src%pos(1,nelem)=Tx%location(1)
                src%pos(2,nelem)=Tx%location(2)
                src%pos(3,nelem)=Tx%location(3)
                !
                src%ljx(nelem)=cos( D2R * Tx%azimuth )
                src%ljy(nelem)=sin( D2R * Tx%azimuth )
                src%ljz(nelem)=0.0
                !
                src%akx(nelem)=0.0
                src%aky(nelem)=0.0
                src%akz(nelem)=0.0
                !
                src%elsrc = .TRUE. 
                !
                allocate(src%cur(src%nelem(1),1),stat=ierr)
                src%cur(1,1)=cmplx(1.0,0.0)
                !
            class default
                stop "Error: create_source_data > Not a CSEM Transmitter"
            !
        end select
        !
    end subroutine create_source_data
    !
end module SourceCSEM_EM1D
