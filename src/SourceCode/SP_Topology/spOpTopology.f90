!
!> Abstract Base class to define a spOpTopology
!
module spOpTopology
    !
    use SpOpTools
    !
    !> Sparse grad topology : maps from all nodes to all edges
    type( spMatCSR_Real ) :: G
    !
    !> Sparse curl topology : maps from all edges to all faces
    type( spMatCSR_Real ) :: T
    !
    type, abstract :: spOpTopology_t
        !
        !> No base properties
        !
        contains
            !
            procedure( interface_curl_spoptopology ), deferred, public :: curl
            procedure( interface_grad_spoptopology ), deferred, public :: grad
            !
    end type spOpTopology_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_curl_spoptopology( self, T ) 
            import :: spOpTopology_t, spMatCSR_Real
            !
            class( spOpTopology_t ), intent( inout ) :: self
            type( spMatCSR_Real ), intent( inout ) :: T
        end subroutine interface_curl_spoptopology
        !
        !> No interface subroutine briefing
        subroutine interface_grad_spoptopology( self, G ) 
            import :: spOpTopology_t, spMatCSR_Real
            !
            class( spOpTopology_t ), intent( inout ) :: self
            type( spMatCSR_Real ), intent( inout ) :: G
        end subroutine interface_grad_spoptopology
        !
    end interface
    !
end module spOpTopology
