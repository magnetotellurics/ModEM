!
!> Derived class to define a CSEM Source with E_p computed using Dipole1D
!
module SourceCSEM_Dipole1D
    !
    use SourceCSEM
    use ModelParameterCell_SG
    !
    type, extends( SourceCSEM_t ) :: SourceCSEM_Dipole1D_t
        !
        real( kind=prec ) :: azimuth, dip, moment
        !
        type( rVector3D_SG_t ) :: cond_anomaly
        !
        contains
            !
            final :: SourceCSEM_Dipole1D_dtor
            !
            procedure, public :: createE => createE_SourceCSEM_Dipole1D
            !
            procedure, public :: createRHS => createRHS_SourceCSEM_Dipole1D
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
    !>     Calls the base routine dealloc().
    !
    subroutine SourceCSEM_Dipole1D_dtor( self )
        implicit none
        !
        type( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor SourceCSEM_Dipole1D_t"
        !
        call self%dealloc
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
        complex( kind=prec ) :: i_omega_mu
        !
        !> Get the Transmitter setting:
        xTx1D = self%location(1)
        yTx1D = self%location(2)
        zTx1D = self%location(3)
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
        !> Verbose...
        write( *, * ) "          - Extract CSEM Source from Dipole 1D"
        !
        call self%set1DModel
        !
        call initilize_1d_vectors( self%sigma%metric%grid ) !> Initilize the 1D vectors where to compupte the e_field field
        !
        call comp_dipole1D !> Calculate e_field-Field by Key's code
        !
        call self%create_Ep_from_Dipole1D( self%sigma%metric%grid )
        !
        deallocate( zlay1D )
        !
        allocate( self%E(1), source = E_p )
        !
        call self%E(1)%mult( self%cond_anomaly )
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
        call self%setCondAnomally( self%cond_anomaly, 1 )
        !
    end subroutine set1DModel_SourceCSEM_Dipole1D
    !
    !> No subroutine briefing
    !
    subroutine initilize_1d_vectors( grid )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid 
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
        !====================================================================
        !> Create position vector that the primary field has to be calculated
        !====================================================================
        counter = 1
        !
        !> e_field-field corresponding to these nodes is Ex
        do iz = 1,grid%Nz+1 !Edge Z
            do iy = 1,grid%Ny+1 !Edge Y
                do ix = 1,grid%Nx !Center X
                    x1D(counter) = grid%x_center(ix)
                    y1D(counter) = grid%y_edge(iy)
                    z1D(counter) = grid%z_edge(iz)
                    counter = counter + 1
                enddo
            enddo
        enddo
        !
        !> e_field-field corresponding to these nodes is Ey
        do iz = 1,grid%Nz+1 !Edge Z
            do iy = 1,grid%Ny !Center y
                do ix = 1,grid%Nx+1 !Edge x
                    x1D(counter) = grid%x_edge(ix)
                    y1D(counter) = grid%y_center(iy)
                    z1D(counter) = grid%z_edge(iz)
                    counter = counter + 1
                enddo
            enddo
        enddo
        !
        !> e_field-field corresponding to these nodes is Ez
        do iz = 1,grid%Nz !Center Z
            do iy = 1,grid%Ny+1 !Edge y
                do ix = 1,grid%Nx+1 !Edge x
                    x1D(counter) = grid%x_edge(ix)
                    y1D(counter) = grid%y_edge(iy)
                    z1D(counter) = grid%z_center(iz)
                    counter = counter + 1
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
        if( allocated( E_p ) ) deallocate( E_p )
        allocate( cVector3D_SG_t :: E_p )
        !
        E_p = cVector3D_SG_t( grid, EDGE )
        !
        select type( E_p )
            !
            class is( cVector3D_SG_t )
                !
                counter = 1
                !
                !> e_field-field corresponding to these nodes is Ex
                do iz = 1,grid%Nz+1 !Edge Z
                    do iy = 1,grid%Ny+1 !Edge Y
                        do ix = 1,grid%Nx !Center X
                            E_p%x(ix,iy,iz) = ex1D(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
                !> e_field-field corresponding to these nodes is Ey
                do iz = 1,grid%Nz+1 !Edge Z
                    do iy = 1,grid%Ny !Center y
                        do ix = 1,grid%Nx+1 !Edge x
                            E_p%y(ix,iy,iz) = ey1D(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
                !> e_field-field corresponding to these nodes is Ez
                do iz = 1,grid%Nz !Center Z
                    do iy = 1,grid%Ny+1 !Edge y
                        do ix = 1,grid%Nx+1 !Edge x
                            E_p%z(ix,iy,iz) = jz1D(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
        end select
        !
        deallocate( x1D, y1D, z1D )
        deallocate( ex1D, ey1D, jz1D )
        deallocate( bx1D, by1D, bz1D )
        !
    end subroutine create_Ep_from_Dipole1D
    !
    !> Set RHS from self%E
    !
    subroutine createRHS_SourceCSEM_Dipole1D( self )
        implicit none
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        !
        if( allocated( self%rhs ) ) deallocate( self%rhs )
        allocate( cVector3D_SG_t :: self%rhs(1) )
        !
        self%rhs(1) = self%E(1)
        !
        call self%rhs(1)%mult( self%model_operator%metric%Vedge )
        !
    end subroutine createRHS_SourceCSEM_Dipole1D
    !
end module SourceCSEM_Dipole1D
