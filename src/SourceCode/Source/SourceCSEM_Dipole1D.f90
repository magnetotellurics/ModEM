!
!> Derived class to define a CSEM Source with E_p computed using Dipole1D
!
module SourceCSEM_Dipole1D
    !
    use SourceCSEM
    !
    type, extends( SourceCSEM_t ) :: SourceCSEM_Dipole1D_t
        !
        real( kind=prec ) :: azimuth, dip, moment
        !
        class( Vector_t ), allocatable :: cond_anomaly
        !
        contains
            !
            final :: SourceCSEM_Dipole1D_dtor
            !
            procedure, public :: createE => createE_SourceCSEM_Dipole1D
            !
            procedure, public :: set1DModel => set1DModel_SourceCSEM_Dipole1D
            !
            procedure, private :: create_Ep_from_Dipole1D
            !
    end type SourceCSEM_Dipole1D_T
    !
    interface SourceCSEM_Dipole1D_t
        module procedure SourceCSEM_Dipole1D_ctor
    end interface SourceCSEM_Dipole1D_t
    !
contains
    !
    !> SourceCSEM_Dipole1D constructor
    !
    function SourceCSEM_Dipole1D_ctor( model_operator, sigma, period, location, dip, azimuth, moment ) result( self )
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        class( ModelParameter_t ), target, intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: period, azimuth, dip, moment, location(3)
        !
        type( SourceCSEM_Dipole1D_t ) :: self
        !
        !write( *, * ) "Constructor SourceCSEM_Dipole1D_t"
        !
        call self%baseInit
        !
        self%model_operator => model_operator
        !
        self%sigma => sigma
        !
        self%period = period
        !
        self%location = location
        !
        self%dip = dip
        !
        self%azimuth = azimuth
        !
        self%moment = moment
        !
        self%non_zero_source = .TRUE.
        !
        self%non_zero_bc = .FALSE.
        !
    end function SourceCSEM_Dipole1D_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine baseDealloc().
    !
    subroutine SourceCSEM_Dipole1D_dtor( self )
        implicit none
        !
        type( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor SourceCSEM_Dipole1D_t"
        !
        call self%baseDealloc
        !
        if( allocated( sig1D ) ) deallocate( sig1D )
        !
        if( allocated( zlay1D ) ) deallocate( zlay1D )
        !
    end subroutine SourceCSEM_Dipole1D_dtor
    !
    !> Set self%E from forward modeling 1D
    !
    subroutine createE_SourceCSEM_Dipole1D( self )
        implicit none
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        !
        type( Grid3D_SG_t ) :: grid_sg
        type( rVector3D_SG_t ) :: cond_anomaly_sg
        complex( kind=prec ) :: i_omega_mu
        integer :: ix, iy, iz
        !
        !> Get the Transmitter setting:
        xTx1D = self%location(1)
        yTx1D = self%location(2)
        zTx1D = self%location(3)
        !
        ftx1D = 1.0d0/self%period
        sdm1D = self%moment        !> (Am), dipole moment. Normalize to unit source moment
        azimuthTx1D = self%azimuth !> (degrees) 
        dipTx1D = self%dip
        !
        HTmethod1D = "kk_ht_201"   !> Use 201 point HT digital filters.
        outputdomain1D = "spatial" !> Assume spatial domain comps
        lbcomp = .FALSE.           !> This is changed to true if magnetics in data file
        lUseSpline1D = .TRUE.      !> Use spline interpolation for faster 1D computations
        linversion = .FALSE.       !> Compute derivatives with respect to self%sigma(layers)
        !
        phaseConvention = "lag"    !> The usual default is lag, where phase becomes larger 
                                   !> positive values with increasing range.
        lenTx1D = 00.d0            !> (m) Dipole length 0 = point dipole
        numIntegPts = 0            !> Number of points to use for Gauss quadrature integration for finite dipole
        !
        !> Verbose
        write( *, "( a52, a14 )" ) "- SourceCSEM_Dipole1D according to: ", trim( get_1d_from )
        !
        call self%set1DModel
        !
        !> Initialize the 1D vectors where to compute the e_field field
        select type( grid => self%sigma%metric%grid )
            !
            class is( Grid3D_SG_t )
                !
                call initilize_1d_vectors( grid )
            !
            class is( Grid3D_MR_t )
                !
                grid_sg = param_grid
                !
                call grid_sg%setAirLayers
                !
                call initilize_1d_vectors( grid_sg )
                !
            class default
                call errStop( "createE_SourceCSEM_Dipole1D > grid must be Grid3D_SG_t" )
            !
        end select
        !
        call comp_dipole1D !> Calculate e_field-Field by Key's code
        !
        call self%create_Ep_from_Dipole1D( self%sigma%metric%grid )
        !
        deallocate( zlay1D )
        !
        !>
        allocate( self%E( 1 ) )
        !
        self%E(1) = self%E_p
        !
        select type( cond_anomaly => self%cond_anomaly )
            !
            class is( rVector3D_SG_t )
                !
                call self%E(1)%mult( cond_anomaly )
                !
            class is( rVector3D_mr_t )
                !
                call cond_anomaly%toSG( cond_anomaly_sg )
                !
                call self%E(1)%mult( cond_anomaly_sg )
                !
            class default
                call errStop( "createE_SourceCSEM_Dipole1D > grid must be Grid3D_SG_t" )
            !
        end select
        !
        i_omega_mu = cmplx( 0., real( -1.0d0 * isign * mu_0 * ( 2.0 * PI / self%period ), kind=prec ), kind=prec )
        !
        call self%E(1)%mult( i_omega_mu )
        !
        call self%createRHS
        !
    end subroutine createE_SourceCSEM_Dipole1D
    !
    !> No subroutine briefing
    !
    subroutine set1DModel_SourceCSEM_Dipole1D( self )
        implicit none
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        !
        integer :: i, k
        real( kind=prec ), dimension( nlay1D ):: sig, zlay0
        !
        call self%setCondAnomally( self%cond_anomaly, 1 )
        !
        ! WORKS FINE WITH EM1D, BUT NOT HERE
        !
        ! !
        ! !> Merge Layers (Michael Commer)
        ! i = nlay1D
        ! !
        ! k = 1
        ! sig(1) = SIGMA_AIR
        ! zlay0(1) = 0d0
        ! !
        ! ! air layer
        ! do i = self%sigma%metric%grid%nzAir + 1, nlay1D
            ! ! if either sig_H or sig_V change from layer k to k+1, add new layer
            ! if( abs( sig1D(i) - sig( k ) ) > SIGMA_MIN ) then
                ! !
                ! k = k + 1
                ! sig( k ) = sig1D(i)
                ! zlay0( k ) = zlay1D(i)
                ! !
            ! endif
        ! enddo
        ! !
        ! ! reset temp. 1D-model arrays
        ! do i = 1, k ! new layers
            ! !
            ! sig1D(i) = sig(i)
            ! zlay1D(i) = zlay0(i)
            ! !
        ! enddo
        ! !
        ! nlay1D = k
        ! !
    end subroutine set1DModel_SourceCSEM_Dipole1D
    !
    !> No subroutine briefing
    !
    subroutine initilize_1d_vectors( grid )
        implicit none
        !
        type( Grid3D_SG_t ), intent( in ) :: grid 
        !
        integer counter, ix, iy, iz
        !
        n1D = ( grid%Nx ) * ( grid%Ny+1 ) * ( grid%Nz+1 )
        n1D = n1D + ( grid%Nx+1 ) * ( grid%Ny ) * ( grid%Nz+1 )
        n1D = n1D + ( grid%Nx+1 ) * ( grid%Ny+1 ) * ( grid%Nz )
        !
        if( allocated( x1D ) ) deallocate( x1D )
        if( allocated( y1D ) ) deallocate( y1D )
        if( allocated( z1D ) ) deallocate( z1D )
        !
        allocate ( x1D(n1D), y1D(n1D), z1D(n1D) )
        !
        if( allocated( ex1D ) ) deallocate( ex1D )
        if( allocated( ey1D ) ) deallocate( ey1D )
        if( allocated( jz1D ) ) deallocate( jz1D )
        !
        allocate ( ex1D(n1D), ey1D(n1D), jz1D(n1D) )
        !
        if( allocated( bx1D ) ) deallocate( bx1D )
        if( allocated( by1D ) ) deallocate( by1D )
        if( allocated( bz1D ) ) deallocate( bz1D )
        !
        allocate ( bx1D(n1D), by1D(n1D), bz1D(n1D) )
        !
        !> Create position vector that the primary field has to be calculated
        !
        counter = 1
        !
        !> e_field-field corresponding to these nodes is Ex
        do iz = 1, grid%Nz+1 !Edge Z
            do iy = 1, grid%Ny+1 !Edge Y
                do ix = 1, grid%Nx !Center X
                    !
                    x1D(counter) = grid%x_center(ix)
                    y1D(counter) = grid%y_edge(iy)
                    z1D(counter) = grid%z_edge(iz)
                    !
                    counter = counter + 1
                    !
                enddo
            enddo
        enddo
        !
        !> e_field-field corresponding to these nodes is Ey
        do iz = 1, grid%Nz+1 !Edge Z
            do iy = 1, grid%Ny !Center y
                do ix = 1, grid%Nx+1 !Edge x
                    !
                    x1D(counter) = grid%x_edge(ix)
                    y1D(counter) = grid%y_center(iy)
                    z1D(counter) = grid%z_edge(iz)
                    !
                    counter = counter + 1
                    !
                enddo
            enddo
        enddo
        !
        !> e_field-field corresponding to these nodes is Ez
        do iz = 1, grid%Nz !Center Z
            do iy = 1, grid%Ny+1 !Edge y
                do ix = 1, grid%Nx+1 !Edge x
                    !
                    x1D(counter) = grid%x_edge(ix)
                    y1D(counter) = grid%y_edge(iy)
                    z1D(counter) = grid%z_center(iz)
                    !
                    counter = counter + 1
                    !
                enddo
            enddo
        enddo
        !
    end subroutine initilize_1d_vectors 
    !
    !> No subroutine briefing
    !
    subroutine create_Ep_from_Dipole1D( self, grid )
        implicit none
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        class( Grid_t ), intent( in ) :: grid
        !
        integer ix, iy, iz, counter
        !
        self%E_p = cVector3D_SG_t( grid, EDGE )
        !
        counter = 1
        !
        !> e_field-field corresponding to these nodes is Ex
        do iz = 1, grid%Nz+1 !Edge Z
            do iy = 1, grid%Ny+1 !Edge Y
                do ix = 1, grid%Nx !Center X
                    !
                    self%E_p%x(ix,iy,iz) = ex1D(counter)
                    counter = counter + 1
                    !
                enddo
            enddo
        enddo
        !
        !> e_field-field corresponding to these nodes is Ey
        do iz = 1, grid%Nz+1 !Edge Z
            do iy = 1, grid%Ny !Center y
                do ix = 1, grid%Nx+1 !Edge x
                    self%E_p%y(ix,iy,iz) = ey1D(counter)
                    counter = counter + 1
                enddo
            enddo
        enddo
        !
        !> e_field-field corresponding to these nodes is Ez
        do iz = 1, grid%Nz !Center Z
            do iy = 1, grid%Ny+1 !Edge y
                do ix = 1, grid%Nx+1 !Edge x
                    self%E_p%z(ix,iy,iz) = jz1D(counter)
                    counter = counter + 1
                enddo
            enddo
        enddo
        !
        deallocate( x1D, y1D, z1D )
        deallocate( ex1D, ey1D, jz1D )
        deallocate( bx1D, by1D, bz1D )
        !
    end subroutine create_Ep_from_Dipole1D
    !
end module SourceCSEM_Dipole1D
!