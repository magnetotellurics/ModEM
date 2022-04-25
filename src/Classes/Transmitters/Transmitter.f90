! *************
! 
! Base class to define a Transmitter
! 
! Last modified at 10/11/2021 by Paulo Werdt
! 
! *************
! 
module Transmitter
    ! 
    use Constants
    use Source
    use ForwardSolver
    use cVector
    !
    type, abstract :: Transmitter_t
        !
        integer :: id, n_pol, fwd_key(8)
        !
        real( kind=prec ) :: period
        !
        class( ForwardSolver_t ), pointer :: forward_solver 
        class( Source_t ), pointer        :: source
        !
        class( cVector_t ), allocatable    :: e_all(:)
        integer, allocatable, dimension(:) :: receiver_indexes
        !
    contains
        !
        procedure, public :: init     => initializeTx
        procedure, public :: dealloc  => deallocateTx
        !
        procedure, public :: setSource => setSourceTx
        procedure, public :: setForwardSolver => setForwardSolverTx
        !
        procedure, public :: updateFwdKey
        !
        procedure, public :: isEqual => isEqualTransmitter
        !
        procedure, public :: updateReceiverIndexesArray
        !
        procedure( interface_solve_fwd_tx ), deferred, public :: solveFWD
        !
        procedure( interface_write_tx ), deferred, public     :: write
        !
    end type Transmitter_t
    !
    abstract interface
        !
        subroutine interface_solve_fwd_tx( self )
            import :: Transmitter_t
            class( Transmitter_t ), intent( inout ) :: self
        end subroutine interface_solve_fwd_tx
        !
        function interface_is_equal_tx( self, other ) result( equal )
            import :: Transmitter_t
            class( Transmitter_t ), intent( in ) :: self, other
            logical                                        :: equal
        end function interface_is_equal_tx
        !
        subroutine interface_write_tx( self )
            import :: Transmitter_t
            class( Transmitter_t ), intent(in) :: self
        end subroutine interface_write_tx
        !
    end interface
    !
    contains
        !
        ! Compare two transmitters
        function isEqualTransmitter( self, other ) result( equal )
            implicit none
            !
            class( Transmitter_t ), intent( in ) :: self
            class( Transmitter_t ), intent( in ) :: other
            logical                              :: equal
            !
            equal = .FALSE.
            !
            if( ABS( self%period - other%period ) < TOL6 ) then
                equal = .TRUE.
            endif
            !
        end function isEqualTransmitter
        !
        subroutine initializeTx( self )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            self%id = 0
            self%n_pol = 0
            call self%updateFwdKey()
            !
            self%period = 0.0
            !
        end subroutine initializeTx
        !
        subroutine deallocateTx( self )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            if( allocated( self%e_all ) ) deallocate( self%e_all )
            !
            if( allocated( self%receiver_indexes ) ) deallocate( self%receiver_indexes )
            !
        end subroutine deallocateTx
        !
        subroutine updateFwdKey( self )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            call date_and_time( values=self%fwd_key )
            !
        end subroutine updateFwdKey
        !
        subroutine setForwardSolverTx( self, forward_solver )
            !
            class( Transmitter_t ), intent( inout )        :: self
            class( ForwardSolver_t ), target, intent( in ) :: forward_solver
            !
            self%forward_solver => forward_solver
            !
        end subroutine setForwardSolverTx
        !
        subroutine setSourceTx( self, source )
            !
            class( Transmitter_t ), intent( inout ) :: self
            class( Source_t ), target, intent( in ) :: source
            !
            self%source => source
            !
        end subroutine setSourceTx
        !
        subroutine updateReceiverIndexesArray( self, new_int )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            integer, intent( in )                   :: new_int
            !
            integer, allocatable, dimension(:)      :: temp_array
            integer                                 :: idx, istat
            !
            if( .NOT. allocated( self%receiver_indexes ) ) then
                allocate( self%receiver_indexes(1) )
                self%receiver_indexes(1) = new_int
            else
                !
                do idx = 1, size( self%receiver_indexes )
                    if ( new_int == self%receiver_indexes( idx ) ) then
                        return
                    end if
                end do
                !
                allocate( temp_array( size( self%receiver_indexes ) + 1 ), STAT=istat )
                temp_array( 1 : size( self%receiver_indexes ) ) = self%receiver_indexes
                temp_array( size( self%receiver_indexes ) + 1 ) = new_int
                self%receiver_indexes = temp_array
                !
                deallocate( temp_array )
                !
            endif
            !
         end subroutine updateReceiverIndexesArray
         !
end module Transmitter
