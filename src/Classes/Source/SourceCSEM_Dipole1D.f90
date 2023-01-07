!
!> Derived class to define a MT Source with boundary data computed by 1D solutions
!
module SourceCSEM_Dipole1D
    !
    use dipole1d
    !
    use Constants
    use cVector3D_SG
    use rVector3D_SG
    use Grid3D_SG
    use Source
    use ModelOperator
    use ModelParameterCell_SG
    !
    type, extends( Source_t ) :: SourceCSEM_Dipole1D_t
        !
        real( kind=prec ) :: azimuth, dip, moment, location(3)
        !
        type( rVector3D_SG_t ) :: cond_anomaly_h
        !
        contains
            !
            final :: SourceCSEM_Dipole1D_dtor
            !
            procedure, public :: createE => createE_SourceCSEM_Dipole1D
            procedure, public :: createRHS => createRHS_SourceCSEM_Dipole1D
            procedure, public :: create_Ep_from_Dipole1D
            procedure, public :: set1DModel
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
        call self%init()
        !
        self%model_operator => model_operator
        self%sigma => sigma
        !
        self%period   = period
        self%location = location
        self%dip      = dip
        self%azimuth  = azimuth
        self%moment   = moment
        !
        self%non_zero_source = .TRUE.
        !
    end function SourceCSEM_Dipole1D_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    subroutine SourceCSEM_Dipole1D_dtor( self )
        implicit none
        !
        type( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor SourceCSEM_Dipole1D_t"
        !
        call self%dealloc()
        !
    end subroutine SourceCSEM_Dipole1D_dtor
    !
    !> Set self%E from forward modeling 1D
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
        sdm1D = self%moment          !> (Am), dipole moment. Normalize to unit source moment
        azimuthTx1D = self%azimuth   !> (degrees) 
        dipTx1D     = self%dip
        !
        HTmethod1D     = "kk_ht_201" !> Use 201 point HT digital filters.
        outputdomain1D = "spatial"   !> Assume spatial domain comps
        lbcomp         = .FALSE.     !> This is changed to true if magnetics in data file
        lUseSpline1D   = .TRUE.      !> Use spline interpolation for faster 1D computations
        linversion     = .FALSE.     !> Compute derivatives with respect to self%sigma(layers)
        !
        phaseConvention = "lag"      !> The usual default is lag, where phase becomes larger 
        !> positive values with increasing range.
        lenTx1D         = 00.d0      !> (m) Dipole length 0 = point dipole
        numIntegPts     = 0          !> Number of points to use for Gauss quadrature integration for finite dipole
        !
        !> Verbose...
        write( *, * ) "          - Extract CSEM Source from Dipole 1D"
        !
        call self%set1DModel( xTx1D, yTx1D )
        !
        call initilize_1d_vectors( self%sigma%metric%grid ) !> Initilize the 1D vectors where to compupte the e_field field
        !
        call comp_dipole1D !> Calculate e_field-Field by Key"s code
        !
        call self%create_Ep_from_Dipole1D( self%sigma%metric%grid )
        !
        deallocate( zlay1D )
        !
        call self%E(1)%mult( self%cond_anomaly_h )
        !
        i_omega_mu = cmplx( 0., real( 1.0d0 * isign * MU_0 * ( 2.0 * PI / self%period ), kind=prec ), kind=prec )
        !
        call self%E(1)%mult( i_omega_mu )
        !
        call self%createRHS
        !
    end subroutine createE_SourceCSEM_Dipole1D
    !
    !> No subroutine briefing
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
        if( allocated( x1D ) ) then  
            deallocate( x1D, y1D, z1D )
            deallocate( ex1D, ey1D, jz1D )
            deallocate( bx1D, by1D, bz1D )
        endif
        !
        allocate ( x1D(n1D), y1D(n1D), z1D(n1D) )
        allocate ( ex1D(n1D), ey1D(n1D), jz1D(n1D) )
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
                    x1D(counter) = grid%xCenter(ix)
                    y1D(counter) = grid%yEdge(iy)
                    z1D(counter) = grid%zEdge(iz)
                    counter = counter + 1
                enddo
            enddo
        enddo
        !
        !> e_field-field corresponing to these nodes is Ey
        do iz = 1,grid%Nz+1 !Edge Z
            do iy = 1,grid%Ny !Center y
                do ix = 1,grid%Nx+1 !Edge x
                    x1D(counter) = grid%xEdge(ix)
                    y1D(counter) = grid%yCenter(iy)
                    z1D(counter) = grid%zEdge(iz)
                    counter = counter + 1
                enddo
            enddo
        enddo
        !
        !> e_field-field corresponing to these nodes is Ez
        do iz = 1,grid%Nz !Center Z
            do iy = 1,grid%Ny+1 !Edge y
                do ix = 1,grid%Nx+1 !Edge x
                    x1D(counter) = grid%xEdge(ix)
                    y1D(counter) = grid%yEdge(iy)
                    z1D(counter) = grid%zCenter(iz)
                    counter = counter + 1
                enddo
            enddo
        enddo
        !
    end subroutine initilize_1d_vectors 
    !
    !> No subroutine briefing
    subroutine create_Ep_from_Dipole1D( self, grid )
        implicit none
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        class( Grid_t ), intent( in ) :: grid
        !
        integer ix, iy, iz, counter
        !
        allocate( cVector3D_SG_t :: self%E(1) )
        !
        self%E(1) = cVector3D_SG_t( grid, EDGE )
        !
        select type( E => self%E(1) )
            !
            class is( cVector3D_SG_t )
                !
                counter = 1
                !
                !> e_field-field corresponding to these nodes is Ex
                do iz = 1,grid%Nz+1 !Edge Z
                    do iy = 1,grid%Ny+1 !Edge Y
                        do ix = 1,grid%Nx !Center X
                            E%x(ix,iy,iz) = ex1D(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
                !> e_field-field corresponding to these nodes is Ey
                do iz = 1,grid%Nz+1 !Edge Z
                    do iy = 1,grid%Ny !Center y
                        do ix = 1,grid%Nx+1 !Edge x
                            E%y(ix,iy,iz) = ey1D(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
                !> e_field-field corresponding to these nodes is Ez
                do iz = 1,grid%Nz !Center Z
                    do iy = 1,grid%Ny+1 !Edge y
                        do ix = 1,grid%Nx+1 !Edge x
                            E%z(ix,iy,iz) = jz1D(counter)
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
    !> No subroutine briefing
    subroutine set1DModel( self, xTx1D, yTx1D )
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        real( kind=prec ),intent( in ) :: xTx1D, yTx1D 
        !
        type( rScalar3D_SG_t ) :: sigma_cell, model
        character( len=80 ) :: param_type
        type( ModelParameterCell_SG_t ) :: aModel
        !
        integer :: nzEarth, nzAir, i, j, k, ixTx, iyTx, counter
        real( kind=prec ) :: wt, asigma, temp_sigma_value
        !
        class( Vector_t ), allocatable :: model_param_map, amodel_map
        !
        !>   first define conductivity on cells  
        !>   (extract into variable which is public)
        !call modelParamToCell(sigma, sigma_cell, param_type)
        !
        select type( sigma => self%sigma )
            class is( ModelParameterCell_SG_t )
                !
                sigma_cell = sigma%cell_cond
                nlay1D = sigma_cell%nz + sigma_cell%grid%nzAir
                nzEarth = sigma_cell%grid%nzEarth
                nzAir = sigma_cell%grid%nzAir
                !
                ixTx = minNode( xTx1D, sigma_cell%grid%xEdge )
                iyTx = minNode( yTx1D, sigma_cell%grid%yEdge )
                !
                if( allocated( zlay1D ) ) deallocate( zlay1D, sig1D )
                !
                allocate( zlay1D(nlay1D) )
                allocate( sig1D(nlay1D) )
                !
                do k=1,nlay1D
                    zlay1D(k) = sigma_cell%grid%zEdge(k)
                enddo
                !
                !> For create sig1D, we divide this process into two parts (1) for air layers and 
                !>    (2) for earth layers
                !> For air layer, sig1D equal to air layer conductivity
                !> For earth layer, The Geometric mean is be used to create sig1D
                !
                sig1D(1:nzAir) = SIGMA_AIR !sigma_cell%v(1,1,1:nzAir)
                !
                !> Verbose
                write( *, * ) "          - Get 1D according to: ", trim(get_1D_from)
                !
                if( trim(get_1D_from) =="Geometric_mean" ) then
                    !
                    do k = nzAir+1,nlay1D
                        wt = R_ZERO
                        temp_sigma_value=R_ZERO
                        !
                        do i = 1,sigma_cell%grid%Nx
                            do j = 1,sigma_cell%grid%Ny
                                wt = wt + sigma_cell%grid%dx(i) * sigma_cell%grid%dy(j)
                                !
                                temp_sigma_value = temp_sigma_value + (sigma_cell%v(i,j,k-nzAir))* &
                                sigma_cell%grid%dx(i) * sigma_cell%grid%dy(j)
                            enddo
                        enddo
                        !
                        sig1D(k) = exp( temp_sigma_value / wt )
                        !
                   enddo
                   !
                else if( trim( get_1D_from ) =="At_Tx_Position" ) then
                    !
                    do k = nzAir+1,nlay1D
                        sig1D(k)=sigma_cell%v(ixTx,iyTx,k-nzAir)
                    enddo
                    !
                else if( trim( get_1d_from ) == "Geometric_mean_around_Tx" ) then
                    do k = nzAir+1,nlay1D
                        !
                        wt = R_ZERO
                        !
                        do i = ixTx-5,ixTx+5
                            do j = iyTx-5,iyTx+5
                                !
                                wt = wt + sigma_cell%grid%dx(i)*sigma_cell%grid%dy(j)
                                !
                                sig1D(k) = sig1D(k) + sigma_cell%v(i,j,k-nzAir) * &
                                sigma_cell%grid%dx(i)*sigma_cell%grid%dy(j)
                                !
                            enddo
                        enddo
                        !
                        sig1D(k) = exp(sig1D(k)/wt)
                        !
                    enddo
                    !
                else if( trim( get_1d_from ) == "Full_Geometric_mean" ) then
                    !
                    wt = R_ZERO
                    temp_sigma_value=R_ZERO
                    counter=0
                    !
                    do k = nzAir+1,nlay1D
                        do i = 1,sigma_cell%grid%Nx
                            do j = 1,sigma_cell%grid%Ny
                                !
                                counter=counter+1
                                !
                                wt = wt + sigma_cell%grid%dx(i)*sigma_cell%grid%dy(j)*sigma_cell%grid%dz(k)
                                !
                                temp_sigma_value = temp_sigma_value + sigma_cell%v(i,j,k-nzAir)
                                !
                            enddo
                        enddo
                    enddo
                    !
                    do k = nzAir+1,nlay1D
                        !
                        sig1D(k) = exp(temp_sigma_value/counter)
                        !
                    enddo
                    !
                else if( trim( get_1d_from ) == "Fixed_Value" ) then
                    !
                    temp_sigma_value = sigma_cell%v( ixTx, iyTx, k-nzAir ) !the value exactly below the Tx
                    !
                    do k = nzAir+1,nlay1D
                        !
                        sig1D(k) = temp_sigma_value
                        !
                    enddo
                    !
                else
                    !
                    stop "Error: set1DModel > Unknow get_1d_from"
                    !
                endif
                !
                model = sigma_cell
                !
                !> Put the background (Primary) "condNomaly" conductivities in ModEM model format
                model%v = R_ZERO
                do k = nzAir+1, nlay1D
                    asigma = sig1D(k)
                    !
                    if( trim( sigma%param_type ) == LOGE ) asigma = log( asigma )
                    !
                    do i = 1,sigma_cell%grid%Nx
                        do j = 1,sigma_cell%grid%Ny
                            model%v( i, j, k-nzAir ) = ( asigma )
                        enddo
                    enddo
                enddo   
                !
                amodel = sigma
                !
                amodel%cell_cond = model
                !
                call sigma%PDEmapping( model_param_map )
                call amodel%PDEmapping( amodel_map )
                !
                self%cond_anomaly_h = model_param_map
                call self%cond_anomaly_h%sub( amodel_map )
                !
        end select
        !
    end subroutine set1DModel
    !
    !> Set RHS from self%E
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
    !> LARGE_REAL is defined in the module math_constants
    !> A.K.
    function clean(x)
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ) :: clean
        !
        clean = dnint(x*LARGE_REAL)/LARGE_REAL
        !
    end function clean
    !
end module SourceCSEM_Dipole1D
