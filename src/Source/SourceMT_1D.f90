!
!> Derived class to define a MT Source with boundary data computed by 1D solutions
!>
!> *************
!>
module SourceMT_1D
    !
    use Constants
    use cVector3D_SG
    use Source
    use ModelOperator
    use ModelParameter1D
    use Forward1D
    !
    type, extends( Source_t ) :: SourceMT_1D_t
        !
        !> No derived properties
        !
        contains
            !
            final :: SourceMT_1D_dtor
            !
            procedure, public :: createE => createESourceMT_1D
            procedure, public :: createRHS => createRHSSourceMT_1D
            !
    end type SourceMT_1D_T
    !
    interface SourceMT_1D_t
        module procedure SourceMT_1D_ctor
    end interface SourceMT_1D_t
    !
contains
    !
    !> SourceMT_1D constructor
    function SourceMT_1D_ctor( model_operator, sigma, period ) result( self )
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        class( ModelParameter_t ), target, intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: period
        !
        type( SourceMT_1D_t ) :: self
        !
        !write( *, * ) "Constructor SourceMT_1D_t"
        !
        call self%init
        !
        self%model_operator => model_operator
        !
        self%sigma => sigma
        !
        self%period = period
        !
        self%non_zero_bc = .TRUE.
        !
    end function SourceMT_1D_ctor
    !
    !> Deconstructor routine:
    !>     Call the base routine dealloc().
    subroutine SourceMT_1D_dtor( self )
        implicit none
        !
        type( SourceMT_1D_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor SourceMT_1D_t"
        !
        call self%dealloc
        !
    end subroutine SourceMT_1D_dtor
    !
    !> Set self%E from forward modeling 1D
    subroutine createESourceMT_1D( self )
        implicit none
        !
        class( SourceMT_1D_t ), intent( inout ) :: self
        !
        type( ModelParameter1D_t ) :: model_parameter_1D
        type( Forward1D_t ) :: forward_1D
        complex( kind=prec ), allocatable, dimension(:) :: E1D
        !
        integer :: ix, iy, pol
        !
        !> Get Model1D from average conductivity 3D
        model_parameter_1D = self%sigma%slice1D( 1, 1 )
        !
        forward_1D = Forward1D_t( model_parameter_1D )
        !
        call forward_1D%setFreq( 2.0 * PI / self%period )
        !
        !> NOTE: E1D is defined at layer interfaces -- dzEdge(nz+1) 
        allocate( E1D( self%model_operator%metric%grid%nz + 1 ) )
        !
        !> Solve 1D and store the result in E1D structure
        call forward_1D%solve( E1D )
        !
        !> 
        allocate( cVector3D_SG_t :: self%E( 2 ) )
        !
        do pol = 1, 2
            !
            self%E( pol ) = cVector3D_SG_t( self%model_operator%metric%grid, EDGE )
            !
            !> Fill e_vector (cVector3D_SG) from E1D (Esoln1DTM_t)
            !>     Note that Ez components are all left set to 9
            select type( E => self%E( pol ) )
                !
                class is( cVector3D_SG_t )
                    !
                    !> 1st polarization case: Only y components are non-zero
                    if( pol == 1 ) then
                        do ix = 1, self%model_operator%metric%grid%nx+1
                            do iy = 1, self%model_operator%metric%grid%ny
                                E%y( ix, iy, : ) = E1D
                            enddo
                        enddo
                    !
                    !> 2nd polarization case: Only x components are non-zero
                    else
                        do ix = 1, self%model_operator%metric%grid%nx
                            do iy = 1, self%model_operator%metric%grid%ny+1
                                E%x( ix, iy, : ) = E1D
                            enddo
                        enddo
                        !
                    endif
                    !
                class default
                    stop "Error: createESourceMT_1D: Unclassified Vector"
            end select
            !
        enddo
        !
        deallocate( E1D )
        !
        call self%createRHS()
        !
    end subroutine createESourceMT_1D
    !
    !> Set RHS from self%E
    subroutine createRHSSourceMT_1D( self )
        implicit none
        !
        class( SourceMT_1D_t ), intent( inout ) :: self
        !
        class( Vector_t ), allocatable :: e_boundary
        !
        integer :: pol
        !
        if( allocated( self%rhs ) ) deallocate( self%rhs )
        allocate( cVector3D_SG_t :: self%rhs( 2 ) )
        !
        do pol = 1, 2
            !
            call self%E( pol )%boundary( e_boundary )
            !
            self%rhs( pol ) = cVector3D_SG_t( self%model_operator%metric%grid, EDGE )
            !
            call self%model_operator%MultAib( e_boundary, self%rhs( pol ) )
            !
            deallocate( e_boundary )
            !
            call self%rhs( pol )%mult( C_MinusOne )
            !
        enddo
        !
    end subroutine createRHSSourceMT_1D
    !
end module SourceMT_1D
