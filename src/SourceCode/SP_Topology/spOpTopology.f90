!
!> Abstract Base class to define a SpOpTopology
!
module SpOpTopology
    !
    use SpOpTools
    !
    !> Sparse grad topology : maps from all nodes to all edges
    type( spMatCSR_Real ) :: G
    !
    !> Sparse curl topology : maps from all edges to all faces
    type( spMatCSR_Real ) :: T
    !
    type, abstract :: SpOpTopology_t
        !
        !> No base properties
        !
        contains
            !
            procedure( interface_curl_SpOpTopology ), deferred, public :: curl
            procedure( interface_grad_SpOpTopology ), deferred, public :: grad
            !
    end type SpOpTopology_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_curl_SpOpTopology( self, T ) 
            import :: SpOpTopology_t, spMatCSR_Real
            !
            class( SpOpTopology_t ), intent( in ) :: self
            type( spMatCSR_Real ), intent( inout ) :: T
        end subroutine interface_curl_SpOpTopology
        !
        !> No interface subroutine briefing
        subroutine interface_grad_SpOpTopology( self, G ) 
            import :: SpOpTopology_t, spMatCSR_Real
            !
            class( SpOpTopology_t ), intent( in ) :: self
            type( spMatCSR_Real ), intent( inout ) :: G
        end subroutine interface_grad_SpOpTopology
        !
    end interface
    !
end module SpOpTopology
