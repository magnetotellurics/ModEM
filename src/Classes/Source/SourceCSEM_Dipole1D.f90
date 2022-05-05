! *************
! 
! Derived class to define a MT Source with boundary data computed by 1D solutions
! 
! Last modified at 10/11/2021 by Paulo Werdt
! 
! *************
! 
!**
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
        type( rVector3D_SG_t ) :: CondAnomaly_h
        !
        contains
            !
            final :: SourceCSEM_Dipole1D_dtor
            !
            procedure, public :: setRHS => setRHS_CSEM_Dipole1D
            procedure, public :: setE   => setE_CSEM_Dipole1D
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
    ! SourceCSEM_Dipole1D constructor
    function SourceCSEM_Dipole1D_ctor( model_operator, model_parameter, period, location, dip, azimuth, moment, E ) result( self )
        implicit none
        !
        class( ModelOperator_t ), target, intent( in )  :: model_operator
        class( ModelParameter_t ), target, intent( in ) :: model_parameter
        real( kind=prec ), intent( in )                 :: period, azimuth, dip, moment, location(3)
        class( cVector_t ), intent( in ), optional      :: E
        !
        type( SourceCSEM_Dipole1D_t ) :: self
        !
        !write( *, * ) "Constructor SourceCSEM_Dipole1D_t"
        !
        call self%init()
        !
        self%non_zero_source = .TRUE.
        !
        self%model_operator  => model_operator
        self%model_parameter => model_parameter
        !
        self%period   = period
        self%location = location
        self%dip      = dip
        self%azimuth  = azimuth
        self%moment   = moment
        !
        if( present( E ) ) then
            !
            self%E = E
            !
            call self%setRHS()
            !
        endif
        !
    end function SourceCSEM_Dipole1D_ctor
    !
    ! SourceCSEM_Dipole1D destructor
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
    ! Set self%E from forward modelling 1D
    subroutine setE_CSEM_Dipole1D( self, polarization )
        implicit none
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        integer, intent( in )                           :: polarization
        !
        ! Get the Transmitter setting:
        xTx1D = self%location(1)
        yTx1D = self%location(2)
        zTx1D = self%location(3)
        ftx1D = 1.0d0/self%period
        sdm1D = self%moment          ! (Am), dipole moment. Normalize to unit source moment
        azimuthTx1D = self%azimuth   ! (degrees) 
        dipTx1D     = self%dip
        !
        HTmethod1D     = "kk_ht_201" ! Use 201 point HT digital filters.
        outputdomain1D = "spatial"   ! Assume spatial domain comps
        lbcomp         = .FALSE.     ! This is changed to true if magnetics in data file
        lUseSpline1D   = .TRUE.      ! Use spline interpolation for faster 1D computations
        linversion     = .FALSE.     ! Compute derivatives with respect to self%model_parameter(layers)
        !
        phaseConvention = "lag"      ! The usual default is lag, where phase becomes larger 
        ! positive values with increasing range.
        lenTx1D         = 00.d0      ! (m) Dipole length 0 = point dipole
        numIntegPts     = 0          ! Number of points to use for Gauss quadrature integration for finite dipole
        !
        ! Verbosis...
        write( *, * ) "    -> Extracting CSEM Source from Dipole 1D"
        !
        call self%set1DModel( xTx1D, yTx1D )
        !
        call initilize_1d_vectors( self%model_parameter%grid ) ! Initilize the 1D vectors where to compupte the E field
        !
        call comp_dipole1D ! Calculate E-Field by Key"s code
        !
        call self%create_Ep_from_Dipole1D( self%model_parameter%grid )
        !
        call self%setRHS
        !
    end subroutine setE_CSEM_Dipole1D
    !
    ! Set RHS from self%E
    subroutine setRHS_CSEM_Dipole1D( self )
        implicit none
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        !
        complex( kind=prec ) :: i_omega_mu
        !
        !
        i_omega_mu = cmplx( 0., real( 1.0d0 * ISIGN * MU_0 * ( 2.0 * PI / self%period ) ), kind=prec )
        !
        select type( grid => self%model_parameter%grid )
            !
            class is( Grid3D_SG_t )
                !
                if( allocated( self%rhs ) ) deallocate( self%rhs )
                allocate( self%rhs, source = cVector3D_SG_t( grid, EDGE ) )
                    !
                    self%rhs = self%E * self%CondAnomaly_h
                    !
                    self%rhs = self%rhs * i_omega_mu
                    !
        end select
        !
    end subroutine setRHS_CSEM_Dipole1D
    !
    subroutine initilize_1d_vectors( grid )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid 
        !
        integer counter, ix, iy, iz
        !
        !
        n1D = ( grid%Nx ) * ( grid%Ny+1 ) * ( grid%Nz+1 )
        n1D = n1D + ( grid%Nx+1 ) * ( grid%Ny ) * ( grid%Nz+1 )
        n1D = n1D + ( grid%Nx+1 ) * ( grid%Ny+1 ) * ( grid%Nz )
        !
        if( allocated( x1D ) ) then  
            deallocate( x1D, y1D, z1D )
            deallocate( ex1D, ey1D, jz1D )
            deallocate( bx1D, by1D, bz1D )
        end if
        !
        allocate ( x1D(n1D), y1D(n1D), z1D(n1D) )
        allocate ( ex1D(n1D), ey1D(n1D), jz1D(n1D) )
        allocate ( bx1D(n1D), by1D(n1D), bz1D(n1D) )
        !
        !====================================================================
        ! Create position vector that the primary field has to be calculated
        !====================================================================
        counter = 1
        !
        ! E-field corresponing to these nodes is Ex
        do iz = 1,grid%Nz+1 !Edge Z
            do iy = 1,grid%Ny+1 !Edge Y
                do ix = 1,grid%Nx !Center X
                    x1D(counter) = grid%xCenter(ix)
                    y1D(counter) = grid%yEdge(iy)
                    z1D(counter) = grid%zEdge(iz)
                    counter = counter + 1
                end do
            end do
        end do
        !
        ! E-field corresponing to these nodes is Ey
        do iz = 1,grid%Nz+1 !Edge Z
            do iy = 1,grid%Ny !Center y
                do ix = 1,grid%Nx+1 !Edge x
                    x1D(counter) = grid%xEdge(ix)
                    y1D(counter) = grid%yCenter(iy)
                    z1D(counter) = grid%zEdge(iz)
                    counter = counter + 1
                end do
            end do
        end do
        !
        ! E-field corresponing to these nodes is Ez
        do iz = 1,grid%Nz !Center Z
            do iy = 1,grid%Ny+1 !Edge y
                do ix = 1,grid%Nx+1 !Edge x
                    x1D(counter) = grid%xEdge(ix)
                    y1D(counter) = grid%yEdge(iy)
                    z1D(counter) = grid%zCenter(iz)
                    counter = counter + 1
                end do
            end do
        end do
        !
    end subroutine initilize_1d_vectors 

    subroutine create_Ep_from_Dipole1D( self, grid )
        implicit none
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        class( Grid_t ), intent(in)                     :: grid
        !
        integer ix, iy, iz, counter
        !
        select type( grid  )
            class is( Grid3D_SG_t )
                !
                if( allocated( self%E ) ) deallocate( self%E )
                allocate( self%E, source = cVector3D_SG_t( grid, EDGE ) )
                !
                select type( E => self%E )
                    !
                    class is( cVector3D_SG_t )
                        !
                        counter = 1
                        !
                        ! E-field corresponing to these nodes is Ex
                        do iz = 1,grid%Nz+1 !Edge Z
                            do iy = 1,grid%Ny+1 !Edge Y
                                do ix = 1,grid%Nx !Center X
                                    E%x(ix,iy,iz) = ex1D(counter)
                                    counter = counter + 1
                                end do
                            end do
                        end do
                        !
                        ! E-field corresponing to these nodes is Ey
                        do iz = 1,grid%Nz+1 !Edge Z
                            do iy = 1,grid%Ny !Center y
                                do ix = 1,grid%Nx+1 !Edge x
                                    E%y(ix,iy,iz) = ey1D(counter)
                                    counter = counter + 1
                                end do
                            end do
                        end do
                        !
                        ! E-field corresponing to these nodes is Ez
                        do iz = 1,grid%Nz !Center Z
                            do iy = 1,grid%Ny+1 !Edge y
                                do ix = 1,grid%Nx+1 !Edge x
                                    E%z(ix,iy,iz) = jz1D(counter)
                                    counter = counter + 1
                                end do
                            end do
                        end do
                        !
                        deallocate( x1D, y1D, z1D )
                        deallocate( ex1D, ey1D, jz1D )
                        deallocate( bx1D, by1D, bz1D )
                        !
                end select
                !
            end select
            !
    end subroutine create_Ep_from_Dipole1D
    !
    subroutine set1DModel( self, xTx1D, yTx1D )
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        real( kind=prec ),intent(in)                    :: xTx1D, yTx1D 
        !
        !
        type( rScalar3D_SG_t )          :: sigma_cell, model
        character( len=80 )             :: paramtype
        type( ModelParameterCell_SG_t ) :: aModel, Anomalous_model
        !
        integer :: nzEarth, nzAir, i, j, k, ixTx, iyTx, counter
        real( kind=prec ) :: wt, asigma, temp_sigma_value
        !
        !   first define conductivity on cells  
        !   (extract into variable which is public)
        !call modelParamToCell(model_parameter, sigma_cell, paramtype)
        !
        select type( model_parameter => self%model_parameter )
            class is( ModelParameterCell_SG_t )
                !
                !
                sigma_cell = model_parameter%cellCond
                nlay1D = sigma_cell%nz+sigma_cell%grid%nzAir
                nzEarth = sigma_cell%grid%nzEarth
                nzAir = sigma_cell%grid%nzAir

                ixTx= minNode(xTx1D, sigma_cell%grid%xEdge)  
                iyTx= minNode(yTx1D, sigma_cell%grid%yEdge)
                !
                if(allocated(zlay1D)) deallocate(zlay1D, sig1D)
                !
                allocate( zlay1D(nlay1D) )
                allocate( sig1D(nlay1D) )
                !
                do k=1,nlay1D
                    zlay1D(k) = sigma_cell%grid%zEdge(k)
                end do
                !
                ! For create sig1D, we divide this process into two parts (1) for air layers and 
                !    (2) for earth layers
                ! For air layer, sig1D equal to air layer conductivity
                ! For earth layer, The Geometric mean is be used to create sig1D
                !
                sig1D(1:nzAir) = SIGMA_AIR !sigma_cell%v(1,1,1:nzAir)
                !
                if( trim(get_1D_from) =="Geometric_mean" ) then
                    !
                    do k = nzAir+1,nlay1D
                        wt = R_ZERO
                        temp_sigma_value=R_ZERO
                        !
                        do i = 1,sigma_cell%grid%Nx
                            do j = 1,sigma_cell%grid%Ny
                                wt = wt + sigma_cell%grid%dx(i)*sigma_cell%grid%dy(j)
                                !
                                temp_sigma_value = temp_sigma_value + (sigma_cell%v(i,j,k-nzAir))* &
                                sigma_cell%grid%dx(i)*sigma_cell%grid%dy(j)
                            end do
                        end do
                        !
                        sig1D(k) = exp(temp_sigma_value/wt)
                        !
                   end do
                   !
                else if( trim( get_1D_from ) =="At_Tx_Position" ) then
                    !
                    do k = nzAir+1,nlay1D
                        sig1D(k)=sigma_cell%v(ixTx,iyTx,k-nzAir)
                    end do
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
                            end do
                        end do
                        !
                        sig1D(k) = exp(sig1D(k)/wt)
                        !
                    end do
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
                            end do
                        end do
                    end do
                    !
                    do k = nzAir+1,nlay1D
                        !
                        sig1D(k) = exp(temp_sigma_value/counter)
                        !
                    end do
                    !
                else if( trim( get_1d_from )== "Fixed_Value" ) then
                    !
                    temp_sigma_value = sigma_cell%v( ixTx, iyTx, k-nzAir ) !the value exactly below the Tx
                    !
                    do k = nzAir+1,nlay1D
                        !
                        sig1D(k) = temp_sigma_value
                        !
                    end do
                    !
                end if
                !
                model = sigma_cell
                !
                ! Put the background (Primary) "condNomaly" conductivities in ModEM model format
                model%v = R_ZERO
                do k = nzAir+1, nlay1D
                    asigma = sig1D(k)
                    !
                    if( trim( model_parameter%ParamType ) == LOGE ) asigma = log( asigma )
                    !
                    do i = 1,sigma_cell%grid%Nx
                        do j = 1,sigma_cell%grid%Ny
                            model%v( i, j, k-nzAir ) = ( asigma )
                        end do
                    end do
                end do   
                !
                amodel = model_parameter
                !
                amodel%cellCond = model
                !
                self%CondAnomaly_h = model_parameter%PDEmapping() - model_parameter%dPDEmapping( amodel )
                !
        end select
        !
    end subroutine set1DModel
    !
    function minNode( x, xNode ) result( ix )
        !    This is a utility routine, used by several data functional
        !    set up routines, and for other interpolation functions
        !    Returns index ix such that    xNode(ix) <= x < xNode(ix+1)
        !    If x is out of range:
        !    x < xNode(1) returns 0; if x> xNode(nx) returns nx
        !    Assumes xNode is strictly increasing; does not check this
        !    NOTE: as presently coded, when xNode is called with center
        !    (face) node positions, this routine will return zero for
        !    the coordinates in the outer half cell nearest the boundary
        !    If evaluation over the complete model domain is to be allowed
        !    a more general interpolation rule will be required.
        !    A.K.: modified to allow input of any size, nx = size(xNode).
        !
        implicit none
        !
        real( kind=prec ), intent( in )               :: x
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
    ! **************************************************************************
    function maxNode(x, xNode) result(ix)
        !    This is a utility routine, used by several data functional
        !    set up routines, and for other interpolation functions
        !    Returns index ix such that    xNode(ix) <= x < xNode(ix+1)
        !    If x is out of range:
        !    x > xNode(1) returns 0; if x< xNode(nx) returns nx
        !    Assumes xNode is strictly decreasing; does not check this
        !    NOTE: as presently coded, when xNode is called with center
        !    (face) node positions, this routine will return zero for
        !    the coordinates in the outer half cell nearest the boundary
        !    If evaluation over the complete model domain is to be allowed
        !    a more general interpolation rule will be required.
        !    A.K.: modified to allow input of any size, nx = size(xNode).
        !
        implicit none
        !
        real( kind=prec ), intent( in )               :: x
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
    function clean(x)
        ! This is a utility routine that provides an expression used to battle
        ! against machine error problems. It returns the same real or real(8)
        ! as the input, but without the extra digits at the end that are often
        ! a cause of wrong comparisons in the if statements. ALWAYS use clean(x)
        ! instead of x in an inequality!!!
        ! LARGE_REAL is defined in the module math_constants
        ! A.K.
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec )               :: clean
        !
        clean = dnint(x*LARGE_REAL)/LARGE_REAL
        !
    end function clean
    !
end module SourceCSEM_Dipole1D
