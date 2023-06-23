!
!> Abstract base class to define a Transmitter
!
module Transmitter
    !
    use FileUnits
    use ForwardSolver
    use ModelParameter
    use ModelOperator
    use SourceCSEM
    use SourceInteriorForce
    use rVector3D_SG
    !
    type, abstract :: Transmitter_t
        !
        integer :: i_tx, n_pol, i_sol, fwd_key(8)
        !
        real( kind=prec ) :: period
        !
        class( Source_t ), allocatable :: source
        !
        class( ForwardSolver_t ), pointer :: forward_solver
        !
        class( Vector_t ), allocatable, dimension(:) :: e_sol_0, e_sol_1, e_sens
        !
        integer, allocatable, dimension(:) :: receiver_indexes
        !
        contains
            !
            procedure( interface_solve_tx ), deferred, public :: solve
            !
            procedure( interface_is_equal_tx ), deferred, public :: isEqual
            !
            procedure( interface_print_tx ), deferred, public :: print
            !
            procedure, public :: init => initializeTx
            !
            procedure, public :: dealloc => deallocateTx
            !
            procedure, public :: updateFwdKey => updateFwdKeyTx
            !
            procedure, public :: updateReceiverIndexesArray
            !
            procedure, public :: setSource => setSourceTx
            !
            procedure, public :: getSolutionVector => getSolutionVectorTx
            !
            procedure, public :: PMult => PMult_Tx
            !
            procedure, public :: PMult_t => PMult_t_Tx
            !
            procedure, public :: writeESolution
            !
    end type Transmitter_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_solve_tx( self )
            import :: Transmitter_t
            class( Transmitter_t ), intent( inout ) :: self
        end subroutine interface_solve_tx
        !
        !> No interface function briefing
        function interface_is_equal_tx( self, other ) result( equal )
            import :: Transmitter_t
            class( Transmitter_t ), intent( in ) :: self, other
            logical :: equal
        end function interface_is_equal_tx
        !
        !> No interface subroutine briefing
        subroutine interface_print_tx( self )
            import :: Transmitter_t
            class( Transmitter_t ), intent( in ) :: self
        end subroutine interface_print_tx
        !
    end interface
    !
    contains
        !
        !> Initialize transmitter base variables (Avoid initialization on declaration).
        subroutine initializeTx( self )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            self%period = R_ZERO
            !
            self%forward_solver => null()
            !
            self%i_tx = 0
            self%n_pol = 0
            self%i_sol = 0
            !
            call self%updateFwdKey()
            !
        end subroutine initializeTx
        !
        !> Free the memory used by all allocatable variables belonging to this transmitter.
        !> Called before anything in the destructor of derived classes.
        !
        subroutine deallocateTx( self )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            if( allocated( self%source ) ) deallocate( self%source )
            !
            if( allocated( E_p ) ) deallocate( E_p )
            !
            if( allocated( self%e_sol_0 ) ) deallocate( self%e_sol_0 )
            !
            if( allocated( self%e_sol_1 ) ) deallocate( self%e_sol_1 )
            !
            if( allocated( self%e_sens ) ) deallocate( self%e_sens )
            !
            if( allocated( self%receiver_indexes ) ) deallocate( self%receiver_indexes )
            !
        end subroutine deallocateTx
        !
        !> No procedure briefing
        !
        subroutine updateFwdKeyTx( self )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            call date_and_time( values=self%fwd_key )
            !
        end subroutine updateFwdKeyTx
        !
        !> Add a receiver index to the integer array (receiver_indexes).
        !> Increasing the size of the array, if the index does not already exist.
        !
        subroutine updateReceiverIndexesArray( self, new_int )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            integer, intent( in ) :: new_int
            !
            integer :: idx, ni
            integer, allocatable, dimension(:) :: temp_array
            !
            if( .NOT. allocated( self%receiver_indexes ) ) then
                !
                allocate( self%receiver_indexes(1) )
                !
                self%receiver_indexes(1) = new_int
            else
                !
                ni = size( self%receiver_indexes )
                !
                do idx = 1, ni
                    if( new_int == self%receiver_indexes( idx ) ) then
                        return
                    endif
                enddo
                !
                allocate( temp_array( ni + 1 ) )
                !
                temp_array( 1 : ni ) = self%receiver_indexes
                !
                temp_array( ni + 1 ) = new_int
                !
                deallocate( self%receiver_indexes )
                !
                allocate( self%receiver_indexes, source = temp_array )
                !
                deallocate( temp_array )
                !
            endif
            !
        end subroutine updateReceiverIndexesArray
        !
        !> Allocate the source of this transmitter if it is allocated.
        !> And define a new source for this transmitter, sent as an argument.
        !
        subroutine setSourceTx( self, source )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            class( Source_t ), intent( in ) :: source
            !
            if( allocated( self%source ) ) deallocate( self%source )
            allocate( self%source, source = source )
            !
        end subroutine setSourceTx
        !
        !> Allocate the source of this transmitter if it is allocated.
        !> And define a new source for this transmitter, sent as an argument.
        !
        subroutine getSolutionVectorTx( self, pol, solution )
            implicit none
            !
            class( Transmitter_t ), intent( in ) :: self
            integer, intent( in ) :: pol
            class( Vector_t ), pointer, intent( out ) :: solution
            !
            if( self%i_sol == 0 ) then
                allocate( solution, source = self%e_sol_0( pol ) )
            else
                allocate( solution, source = self%e_sol_1( pol ) )
            endif
            !
        end subroutine getSolutionVectorTx
        !
        !> Returns a SourceInteriorForce from two distinct models, with the same ModelOperator.
        !
        function PMult_Tx( self, sigma, dsigma, model_operator ) result( source_int_force )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in ) :: sigma, dsigma
            class( ModelOperator_t ), intent( in ) :: model_operator
            !
            type( SourceInteriorForce_t ) :: source_int_force
            !
            class( Vector_t ), allocatable, dimension(:) :: bSrc
            class( Vector_t ), pointer :: solution
            type( rVector3D_SG_t ) :: map_e_vector
            complex( kind=prec ) :: minus_i_omega_mu
            integer :: pol
            !
            ! Verbose
            !write( *, * ) "               - Start PMult"
            !
            !> Get map_e_vector from dPDEmapping
            call sigma%dPDEmapping( dsigma, map_e_vector )
            !
            !> ON WORKING
            minus_i_omega_mu = -isign * mu_0 * cmplx( 0., ( 2.0 * PI / self%period ), kind=prec )
            !
            !> Initialize and fill bSrc
            allocate( cVector3D_SG_t :: bSrc( self%n_pol ) )
            !
            do pol = 1, self%n_pol
                !
                call self%getSolutionVector( pol, solution )
                !
                bSrc( pol ) = solution
                !
                deallocate( solution )
                !
                call bSrc( pol )%mult( map_e_vector )
                !
                call bSrc( pol )%mult( minus_i_omega_mu )
                !
            enddo
            !
            !> Instantiates the source and sets its E, creating Rhs from it.
            source_int_force = SourceInteriorForce_t( model_operator, sigma, self%period )
            !
            call source_int_force%setE( bSrc )
            !
            !> Free up local memory
            deallocate( bSrc )
            !
            ! Verbose
            !write( *, * ) "               - Finish PMult"
            !
        end function PMult_Tx
        !
        !> Defines a new model (dsigma) from a previous model and e_sens for this transmitter.
        !
        subroutine PMult_t_Tx( self, sigma, dsigma )
            implicit none
            !
            class( Transmitter_t ), intent( in ) :: self
            class( ModelParameter_t ), intent( in ) :: sigma
            class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
            !
            class( Vector_t ), allocatable, dimension(:) :: eSens
            class( Vector_t ), pointer :: solution
            class( Vector_t ), allocatable :: real_sens
            complex( kind=prec ) :: minus_i_omega_mu
            integer :: pol
            !
            ! Verbose
            !write( *, * ) "               - Start PMult_t"
            !
            if( .NOT. allocated( self%e_sens ) ) then
                stop "Error: PMult_t_Tx > eSens not allocated on the Tx"
            endif
            !
            !> Copy e_sens to a local variable to keep its original value.
            allocate( eSens, source = self%e_sens )
            !
            call self%getSolutionVector( 1, solution )
            !
            call eSens(1)%mult( solution )
            !
            deallocate( solution )
            !
            !> Loop over all other polarizations, adding them to the first position
            do pol = 2, self%n_pol
                !
                call self%getSolutionVector( pol, solution )
                !
                call eSens( pol )%mult( solution )
                !
                deallocate( solution )
                !
                call eSens(1)%add( eSens( pol ) )
                !
            enddo
            !
            minus_i_omega_mu = -isign * mu_0 * cmplx( 0., ( 2.0 * PI / self%period ), kind=prec )
            !
            call eSens(1)%mult( minus_i_omega_mu )
            !
            call eSens(1)%getReal( real_sens )
            !
            !> Free up local memory
            deallocate( eSens )
            !
            !> Get dsigma from dPDEmapping_T, using first position of eSens
            call sigma%dPDEmapping_T( real_sens, dsigma )
            !
            deallocate( real_sens )
            !
        end subroutine PMult_t_Tx
        !
        !> No subroutine briefing
        !
        subroutine writeESolution( self )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            integer :: i_pol, ios
            !
            character( len=20 ) :: ModeName
            character :: char_i_pol
            real( kind=prec ) :: omega
            !
            omega = ( 2.0 * PI / self%period )
            !
            !> Loop over all polarizations (MT n_pol = 2)
            do i_pol = 1, self%n_pol
                !
                write( char_i_pol, "(i1)") i_pol
                ModeName = "E"//char_i_pol
                !
                open( ioESolution, action = "write", position = "append", form = "unformatted", iostat = ios )
                !
                if( ios /= 0 ) then
                    stop "Error opening file in solveTransmitterMT: e_solution"
                else
                    !
                    !> write the frequency header - 1 record
                    write( ioESolution ) omega, self%i_tx, i_pol, ModeName
                    !
                    call self%e_sol_0( i_pol )%write( ioESolution )
                    !
                    close( ioESolution )
                    !
                endif
                !
            enddo
            !
        end subroutine writeESolution
        !
end module Transmitter
