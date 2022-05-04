! *************
! 
! Base class to define a Source
! 
! Last modified at 10/11/2021 by Paulo Werdt
! 
! *************
! 
module Source
    !
    use cVector
    use ModelOperator
    use ModelParameter
    !
    character(:), allocatable :: source_type
    character ( len=11 ), parameter :: SRC_MT_1D = "SourceMT_1D"
    character ( len=11 ), parameter :: SRC_MT_2D = "SourceMT_2D"
	!
    character(:), allocatable :: get_1D_from
	!
    type, abstract :: Source_t
        !
        real( kind=prec ) :: period
        integer :: polarization
        !
        class( ModelOperator_t ), pointer  :: model_operator
        class( ModelParameter_t ), pointer :: model_parameter
		!
        class( cVector_t ), allocatable    :: rhs, E
        !
        logical                            :: non_zero_source, adjt
        !
        contains
            !
            procedure, public :: init    => initializeSource
            procedure, public :: dealloc => deallocateSource
            !
            procedure( interface_set_rhs ), deferred, public :: setRHS
            procedure( interface_set_e ), deferred, public   :: setE
            !
    end type Source_t
    !
    abstract interface
        !
        subroutine interface_set_rhs( self )
           !
           import :: Source_t
           class( Source_t ), intent( inout ) :: self
           !
        end subroutine interface_set_rhs
        !
        subroutine interface_set_e( self, polarization )
           !
           import :: Source_t, prec
           class( Source_t ), intent( inout ) :: self
           integer, intent( in )              :: polarization
           !
        end subroutine interface_set_e
        !
    end interface
    !
    contains
    !
    subroutine initializeSource( self )
        implicit none
        !
        class( Source_t ), intent( inout ) :: self
        !
        self%period = 0.0
        self%polarization = 0
        !
        self%non_zero_source = .FALSE.
        self%adjt            = .FALSE.
        !
    end subroutine initializeSource
    !
    subroutine deallocateSource( self )
        implicit none
        !
        class( Source_t ), intent( inout ) :: self
        !
        if( allocated( self%rhs ) ) deallocate( self%rhs )
        if( allocated( self%E ) ) deallocate( self%E )
        !
    end subroutine deallocateSource
    !
end module Source
