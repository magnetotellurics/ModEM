!
!> Derived class to define a CSEM Source with E_p computed using EM1D
!
module SourceCSEM_EM1D
    !
    use EM1D
    !
    use Source
    use Constants
    use cVector3D_SG
    use rVector3D_SG
    use Grid3D_SG
    use ModelOperator
    use ModelParameterCell_SG
    use TransmitterArray
    use TransmitterCSEM
    !
    type, extends( Source_t ) :: SourceCSEM_EM1D_t
        !
        integer :: i_tx
        !
        real( kind=prec ) :: location(3)
        !
        type( rVector3D_SG_t ) :: cond_anomaly_h, cond_anomaly_v
        type( rVector3D_SG_t ) :: cond_nomaly_h, cond_nomaly_v
        !
        contains
            !
            final :: SourceCSEM_EM1D_dtor
            !
            procedure, public :: createE => createE_SourceCSEM_EM1D
            !
            procedure, public :: createRHS => createRHS_SourceCSEM_EM1D
            !
            procedure, private :: create_Ep_from_EM1D, set1DModel_VTI, setAnomConductivity_VTI
            !
            procedure, private :: create_background_data, create_source_data
            !
    end type SourceCSEM_EM1D_T
    !
    integer, private :: nlay1D_temp
    !
    real( kind=prec ), dimension(:), allocatable, private :: zlay1D_temp     ! (m)   Depth to top of each layer, first layer ignored 
    real( kind=prec ), dimension(:), allocatable, private :: sig1D_temp_h    ! (S/m) Layer conductivities in the horizontal direction used in VTI
    real( kind=prec ), dimension(:), allocatable, private :: sig1D_temp_v    ! (S/m) Layer conductivities in the vertical direction used in VTI   
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
        type( backgrounddata ) :: bgdat    !model description, coordinates, data
        type( sorec ) :: src               !source specification
        type( freqdata ) :: freqdat        !frequency dependent specifications
        type( refl_struct ) :: refl_var    !all variables that have to be remembered while computing 1D fields
        !
        integer :: ifreq, icur, comm, istat
        complex( kind=prec ) :: i_omega_mu 
        !
        !> Verbose...
        write( *, * ) "          - Extract CSEM Source from EM1D"
        !
        call self%set1DModel_VTI( self%sigma, self%location(1), self%location(2) )
        !
        call self%setAnomConductivity_VTI()
        !
        bgdat%omega = ( 2.0 * PI / self%period )
        bgdat%dowhat = 1
        !
        call self%create_background_data( bgdat )
        !
        call self%create_source_data( src, freqdat )
        !
        ifreq=1
        icur=1
        comm=1
        refl_var%nzrecHxy=0
        refl_var%nzrecHz=0
        call reflectivity_unified( src, bgdat, refl_var, ifreq, icur, comm ) ! Output field will be saved in bgdat
        !
        call self%create_Ep_from_EM1D( self%sigma%metric%grid, bgdat )       ! Put the 1D field into E_p
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
                        E%x = self%cond_anomaly_h%x * E_P%x
                        E%y = self%cond_anomaly_h%y * E_P%y
                        E%z = self%cond_anomaly_v%z * E_P%z
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
        deallocate( bgdat%sigv, stat = istat )
        deallocate( bgdat%sigh, stat = istat )
        deallocate( bgdat%epsrv, stat = istat )
        deallocate( bgdat%epsrh, stat = istat )
        deallocate( bgdat%zbound, stat = istat )
        deallocate( bgdat%Expos, stat = istat )
        deallocate( bgdat%Eypos, stat = istat )
        deallocate( bgdat%Ezpos, stat = istat )
        deallocate( bgdat%Ex, stat = istat )
        deallocate( bgdat%Ey, stat = istat )
        deallocate( bgdat%Ez, stat = istat )
        !
        deallocate( src%nelem, stat = istat )
        deallocate( src%pos, stat = istat )
        deallocate( src%ljx, stat = istat )
        deallocate( src%ljy, stat = istat )
        deallocate( src%ljz, stat = istat )
        deallocate( src%akx, stat = istat )
        deallocate( src%aky, stat = istat )
        deallocate( src%akz, stat = istat )
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
                ! E-field corresponding to these nodes is Ex
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
                ! E-field corresponing to these nodes is Ey
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
                ! E-field corresponing to these nodes is Ez
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
    !> this is a private routine, used to extract layer averages from
    !> a 3D conductivity parameter (sigma) and set up
    !> (1) nlay1D    ! Number of layers
    !> (2) sig1D => ! (S/m) Layer conductivities 
    !> (3) zlay1D => ! (m)   Depth to top of each layer, first layer ignored the 1D Model  z_P, sigma_P
    !
    subroutine set1DModel_VTI( self, sigma, xTx1D, yTx1D, FromFile )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: sigma 
        real( kind=prec ), intent( in ) :: xTx1D, yTx1D 
        logical, intent( in ), optional :: FromFile
        !
        !   local variables ... this is an easy, but not necessarily most efficient
        !   way to get an average background layered conductivity ...
        !    could add routines to modelParameter module to do this more directly
        character(:), allocatable :: param_type
        type( rScalar3D_SG_t ) :: model, sigmaCell_h, sigmaCell_v 
        class( ModelParameter_t ), allocatable :: aModel, Anomalous_model
        integer :: nzEarth, i, j, k, ixTx, iyTx, counter, Nx, Ny, Nz, nzAir
        real( kind=prec ) :: wt, vAir, asigma, temp_sigma_value
        !
        ! first define conductivity on cells: in h and v
        ! call get_vti(sigma)
        ! write(*,*) 'In get 1D'
        call sigma%modelParamToCell( sigmaCell_h, param_type, cCond_v=sigmaCell_v )
        !
        ! Get the grid spec from either h or v sigmaCell
        Nx = sigmaCell_h%grid%Nx
        Ny = sigmaCell_h%grid%Ny
        Nz = sigmaCell_h%grid%Nz
        !
        nzEarth = sigmaCell_h%grid%nzEarth
        nzAir = sigmaCell_h%grid%nzAir
        !  write(55,*)Nx,Ny,Nz,nzAir
        nlay1D_temp = Nz
        !
        ! Get the Tx position (cell #) in X and Y
        ixTx = minNode( xTx1D, sigmaCell_h%grid%xEdge )
        iyTx = minNode( yTx1D, sigmaCell_h%grid%yEdge )
        !
        ! for layer boundaries use z-edges of 3D grid
        if( allocated( zlay1D_temp ) ) then
            deallocate( zlay1D_temp, sig1D_temp_h, sig1D_temp_v )
        endif
        !
        allocate( zlay1D_temp( nlay1D_temp ) )
        allocate( sig1D_temp_h( nlay1D_temp ) )
        allocate( sig1D_temp_v( nlay1D_temp ) )
        do k=1, nlay1D_temp
            zlay1D_temp(k) = sigmaCell_h%grid%zEdge(k)
        enddo
        !
        ! For create sig1D, we divide this process into two parts (1) for air layers and 
        !    (2) for earth layers
        ! For air layer, sig1D equal to air layer conductivity
        ! For earth layer, The Geometric mean is be used to create sig1D
        sig1D_temp_h(1:nzAir) = sigmaCell_h%v(1,1,1:nzAir)
        sig1D_temp_v(1:nzAir) = sigmaCell_v%v(1,1,1:nzAir)
        !
        !> Verbose
        write( *, * ) "          - Get 1D according to: ", trim( get_1d_from )
        !
        if(trim(get_1d_from) == "Geometric_mean") then
            !
            do k = nzAir+1,nlay1D_temp
                wt = R_ZERO
                temp_sigma_value=R_ZERO
                do i = 1,Nx
                    do j = 1,Ny
                        if(log(sigmaCell_h%v(i,j,k)) .gt. -20.0 ) then
                            wt = wt + sigmaCell_h%grid%dx(i)*sigmaCell_h%grid%dy(j)
                            temp_sigma_value = temp_sigma_value + log(sigmaCell_h%v(i,j,k))* &
                            sigmaCell_h%grid%dx(i)*sigmaCell_h%grid%dy(j)
                        endif
                    enddo
                enddo
                sig1D_temp_h(k) = exp(temp_sigma_value/wt)
            enddo
            !
            do k = nzAir+1,nlay1D_temp
                wt = R_ZERO
                temp_sigma_value=R_ZERO
                do i = 1,Nx
                    do j = 1,Ny
                        if(log(sigmaCell_v%v(i,j,k)) .gt. -20.0  ) then
                            wt = wt + sigmaCell_v%grid%dx(i)*sigmaCell_v%grid%dy(j)
                            temp_sigma_value = temp_sigma_value + log(sigmaCell_v%v(i,j,k))* &
                            sigmaCell_v%grid%dx(i)*sigmaCell_v%grid%dy(j)
                        endif
                    enddo
                enddo
                sig1D_temp_v(k) = exp(temp_sigma_value/wt)
            enddo
        !
        elseif( trim(get_1d_from) == "Tx_Position" ) then
            !
            do k = nzAir+1,nlay1D_temp
                sig1D_temp_h(k)=sigmaCell_h%v(ixTx,iyTx,k)
            enddo
            do k = nzAir+1,nlay1D_temp
                sig1D_temp_v(k)=sigmaCell_v%v(ixTx,iyTx,k)
            enddo  
            do k = nzAir+1,nlay1D_temp
                write(70,*) k, sig1D_temp_h(k),sig1D_temp_v(k)
            enddo
            !
        elseif( trim(get_1d_from) == "Mean_around_Tx" ) then
            !
            do k = nzAir+1,nlay1D_temp
                wt = R_ZERO
                do i = ixTx-5,ixTx+5
                    do j = iyTx-5,iyTx+5
                        if(log(sigmaCell_h%v(i,j,k)) .gt. -20.0  ) then
                            wt = wt + sigmaCell_h%grid%dx(i)*sigmaCell_h%grid%dy(j)
                            sig1D_temp_h(k) = sig1D_temp_h(k) + log(sigmaCell_h%v(i,j,k))* &
                            sigmaCell_h%grid%dx(i)*sigmaCell_h%grid%dy(j)
                        endif
                    enddo
                enddo
                sig1D_temp_h(k) = exp(sig1D_temp_h(k)/wt)
            enddo
            !
            do k = nzAir+1,nlay1D_temp
                wt = R_ZERO
                do i = ixTx-5,ixTx+5
                    do j = iyTx-5,iyTx+5
                        if(log(sigmaCell_v%v(i,j,k)) .gt. -20.0  ) then
                            wt = wt + sigmaCell_v%grid%dx(i)*sigmaCell_v%grid%dy(j)
                            sig1D_temp_v(k) = sig1D_temp_v(k) + log(sigmaCell_v%v(i,j,k))* &
                            sigmaCell_v%grid%dx(i)*sigmaCell_v%grid%dy(j)
                        endif
                    enddo
                enddo
                sig1D_temp_v(k) = exp(sig1D_temp_v(k)/wt)
            enddo
            !
        elseif( trim( get_1d_from ) == "Geometric_mean" ) then
            wt = R_ZERO
            temp_sigma_value=R_ZERO
            counter=0
            do k = nzAir+1,nlay1D_temp
                do i = 1,Nx
                    do j = 1,Ny
                        if(log(sigmaCell_h%v(i,j,k)) .gt. -20.0  ) then
                            counter=counter+1
                            wt = wt + sigmaCell_h%grid%dx(i)*sigmaCell_h%grid%dy(j)*sigmaCell_h%grid%dz(k)
                            temp_sigma_value = temp_sigma_value + log(sigmaCell_h%v(i,j,k))
                        endif
                    enddo
                enddo
            enddo
            do k = nzAir+1,nlay1D_temp
                sig1D_temp_h(k) = exp(temp_sigma_value/counter)    
            enddo
            !
            wt = R_ZERO
            temp_sigma_value=R_ZERO
            counter=0
            do k = nzAir+1,nlay1D_temp
                do i = 1,Nx
                    do j = 1,Ny
                        if(log(sigmaCell_v%v(i,j,k)) .gt. -20.0  ) then
                            counter=counter+1
                            wt = wt + sigmaCell_v%grid%dx(i)*sigmaCell_v%grid%dy(j)*sigmaCell_v%grid%dz(k)
                            temp_sigma_value = temp_sigma_value + log(sigmaCell_v%v(i,j,k))
                        endif
                    enddo
                enddo
            enddo
            !
            do k = nzAir+1,nlay1D_temp
                sig1D_temp_v(k) = exp( temp_sigma_value / counter )
            enddo
            !
        elseif( trim( get_1d_from ) == "Fixed" ) then
            !
            temp_sigma_value=sigmaCell_h%v(ixTx,iyTx,nzAir+1) !the value exactly below the Tx
            do k = nzAir+1,nlay1D_temp
                sig1D_temp_h(k) = temp_sigma_value
            enddo
            !
            temp_sigma_value=sigmaCell_v%v(ixTx,iyTx,nzAir+1) !the value exactly below the Tx
            do k = nzAir+1,nlay1D_temp
                sig1D_temp_v(k) = temp_sigma_value
            enddo
            !
        endif
        !
        !call sigma%getValue( param_type, model, vAir )  !It is just to get model structure (place holder)
        !
        ! Put the background (Primary) h "condNomaly" conductivities in ModEM model format
        model%v=R_ZERO
        do k = 1,nzEarth
            asigma = sig1D_temp_h(k+nzAir)
            if( trim( param_type ) == LOGE ) asigma = log( asigma )
            do i = 1, Nx
                do j = 1, Ny
                    model%v(i,j,k) = asigma
                enddo
            enddo
        enddo
        !
        allocate( amodel, source = sigma )
        !
        call amodel%setType( param_type )
        !
        !call amodel%setValue( param_type, model, vAir )
        !
        call amodel%PDEmapping( self%cond_nomaly_h )
        !
        ! Put the background (Primary) v "condNomaly" conductivities in ModEM model format
        model%v=R_ZERO
        do k = 1,nzEarth
            asigma = sig1D_temp_v(k+nzAir)
            if( trim( param_type ) == LOGE ) asigma = log( asigma )
            do i = 1,Nx
                do j = 1,Ny
                    model%v(i,j,k) = asigma
                    !
                enddo
            enddo
        enddo
        !
        amodel = sigma
        call amodel%setType( param_type )
        !
        !call amodel%setValue( param_type, model, vAir )
        !
        call amodel%PDEmapping( self%cond_nomaly_v )
        !
        deallocate( amodel )
        !
    end subroutine set1DModel_VTI
    !
    !> This is a private routine that sets anomalous conductivity
    !> in module variable condAnomaly using input model parameter sigma,
    !> and layered background conductivity (already set in module 
    !> variables z_P and sigma_P by a call to setPrimaryCond)
    !
    subroutine setAnomConductivity_VTI( self )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( inout ) :: self
        !
        class( Vector_t ), allocatable :: ccond
        !
        allocate( ccond, source = rVector3D_SG_t( self%sigma%metric%grid, EDGE ) )
        !
        ! map conductivity onto edges
        call self%sigma%PDEMapping( ccond )
        !
        self%cond_anomaly_h = ccond
        call self%cond_anomaly_h%sub( self%cond_nomaly_h )
        !
        self%cond_anomaly_v = ccond
        call self%cond_anomaly_v%sub( self%cond_nomaly_v )
        !
        deallocate( ccond )
        !
    end subroutine setAnomConductivity_VTI
    !
    ! No Subroutine briefing
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
        bgdat%aniso = vti !default vti, change to iso if any layer is isotropic!!! OR keep it VTI in all cases: if isotropic case, sig1D_temp_h==sig1D_temp_v 
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
        ! allocate(bgdat%Hxpos(bgdat%nHx,3),bgdat%Hypos(bgdat%nHy,3),bgdat%Hzpos(bgdat%nHz,3), stat=ierr)
        ! if(ierr.ne.0) call alloc_error(pid,'In-backgroundfield','Hpos',ierr)
        !
        ix = bgdat%nEx*1.5
        iy = bgdat%nEy*1.5
        iz = bgdat%nEz*1.5
        allocate(bgdat%Ex(ix),bgdat%Ey(iy),bgdat%Ez(iz), stat=ierr)
        if( ierr .NE. 0 ) call alloc_error(pid,'Out-backgroundfield','E fields',ierr)    
        !
        ! allocate(bgdat%Hx(nxyz),bgdat%Hy(nxyz),bgdat%Hz(nxyz), stat=ierr)
        ! if(ierr.ne.0) call alloc_error(pid,'Out-backgroundfield','H fields',ierr)  
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
        ! E-field corresponding to these nodes is Ex
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
        ! E-field corresponing to these nodes is Ey
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
        ! E-field corresponing to these nodes is Ez
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
    ! No Subroutine briefing
    !
    subroutine create_source_data( self, src, freqdat )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( in ) :: self
        type( sorec ), intent( inout ) :: src          !source specification
        type( freqdata ), intent( inout ) :: freqdat    !frequency dependent specifications
        !
        integer :: nelem, ierr
        class( Transmitter_t ), pointer :: Tx
        !
        ! fill the Freq object with the required information
        freqdat%nfreq = 1
        !allocate frequency vector
        allocate( freqdat%omega(1), stat = ierr )
        !
        Tx => getTransmitter( self%i_tx )
        !
        !> Instantiate Transmitter's Source - According to transmitter type and chosen via control file
        select type( Tx )
            !
            class is( TransmitterCSEM_t )
                !
                freqdat%omega(1)= 2.0 * PI / Tx%period     !angular frequencies (2*pi*f, f in Hertz)
                !
                ! fill the sources object with the required information    
                src%type=1   !receiver = 0, dipole = 1, wire = 2, star = 3
                src%srcname="Tx"  ! dummy name
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
    !>    This is a utility routine, used by several data functional
    !>    set up routines, and for other interpolation functions
    !>    Returns index ix such that    xNode(ix) <= x < xNode(ix+1)
    !>    If x is out of range:
    !>    x < xNode(1) returns 0; if x> xNode(nx) returns nx
    !>    Assumes xNode is strictly increasing; does not check this
    !>    NOTE: as presently coded, when xNode is called with center
    !>    (face) node positions, this routine will return zero for
    !>    the coordinates in the outer half cell nearest the boundary
    !>    If evaluation over the complete model domain is to be allowed
    !>    a more general interpolation rule will be required.
    !>    A.K.: modified to allow input of any size, nx = size(xNode).
    !
    function minNode( x, xNode ) result( ix )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ), dimension(:), intent( in ) :: xNode
        !
        integer :: ix, i
        !
        do i = 1, size( xNode )
            if( clean( xNode(i) ) .GT. clean(x) ) then
                ix = i-1
                exit
            endif
        enddo
        !
    end function minNode
    !
    !>    This is a utility routine, used by several data functional
    !>    set up routines, and for other interpolation functions
    !>    Returns index ix such that    xNode(ix) <= x < xNode(ix+1)
    !>    If x is out of range:
    !>    x > xNode(1) returns 0; if x< xNode(nx) returns nx
    !>    Assumes xNode is strictly decreasing; does not check this
    !>    NOTE: as presently coded, when xNode is called with center
    !>    (face) node positions, this routine will return zero for
    !>    the coordinates in the outer half cell nearest the boundary
    !>    If evaluation over the complete model domain is to be allowed
    !>    a more general interpolation rule will be required.
    !>    A.K.: modified to allow input of any size, nx = size(xNode).
    !
    function maxNode(x, xNode) result(ix)
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ), dimension(:), intent( in ) :: xNode
        !
        integer :: ix, i
        !
        do i = 1, size(xNode)
           if( clean( xNode(i)) .LT. clean(x) ) then
                ix = i-1
                exit
           endif
        enddo
        !
    end function maxNode
    !
    !> This is a utility routine that provides an expression used to battle
    !> against machine error problems. It returns the same real or real(8)
    !> as the input, but without the extra digits at the end that are often
    !> a cause of wrong comparisons in the if statements. ALWAYS use clean(x)
    !> instead of x in an inequality!!!
    !> R_LARGE is defined in the module math_constants
    !> A.K.
    !
    function clean( x )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ) :: clean
        !
        clean = dnint(x*R_LARGE)/R_LARGE
        !
    end function clean
    !
end module SourceCSEM_EM1D
