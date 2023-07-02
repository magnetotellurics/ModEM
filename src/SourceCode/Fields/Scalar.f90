!
!> Abstract class to define a Scalar 
!
module Scalar
    !
    use Field
    !
    type, abstract, extends( Field_t ) :: Scalar_t
        !
        integer, dimension(3) :: NdV
        !
        integer :: Nxyz
        !
    contains
        !
        !> Scalar Interfaces
        procedure( interface_get_v_scalar ), deferred, public :: getV
        procedure( interface_set_v_scalar ), deferred, public :: setV
        !
    end type Scalar_t
    !
    !>
    abstract interface
        !
        !> No interface function briefing
        !
        function interface_get_v_scalar( self ) result( v )
            import :: Scalar_t, prec
            !
            class( Scalar_t ), intent( in ) :: self
            !
            complex( kind=prec ), allocatable :: v(:, :, :)
            !
        end function interface_get_v_scalar
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_v_scalar( self, v )
            import :: Scalar_t, prec
            !
            class( Scalar_t ), intent( inout ) :: self
            complex( kind=prec ), allocatable, intent( in ) :: v(:, :, :)
            !
        end subroutine interface_set_v_scalar
        !
    end interface
    !
end module Scalar
