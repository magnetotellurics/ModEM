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
        !> No derived methods
        !
    end type Scalar_t
    !
end module Scalar
