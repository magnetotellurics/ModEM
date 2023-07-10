!
!> Derived class to define a MT Transmitter
!
module TransmitterMT
    !
    use Transmitter
    !
    type, extends( Transmitter_t ), public :: TransmitterMT_t
        !
        !> No derived properties
        !
        contains
            !
            final :: TransmitterMT_dtor
            !
            procedure, public :: solve => solveTransmitterMT
            !
            procedure, public :: isEqual => isEqualTransmitterMT
            !
            procedure, public :: print => printTransmitterMT
            !
    end type TransmitterMT_t
    !
    interface TransmitterMT_t
        module procedure TransmitterMT_ctor
    end interface TransmitterMT_t
    !
    contains
    !
    !> TransmitterMT constructor
    !
    function TransmitterMT_ctor( period ) result ( self )
        implicit none
        !
        type( TransmitterMT_t ) :: self
        !
        real( kind=prec ), intent( in ) :: period
        !
        !write( *, * ) "Constructor TransmitterMT_t"
        !
        call self%baseInit
        !
        self%n_pol = 2
        !
        self%period = period
        !
        !> self%PMult_ptr => PMult_E
        !
        !> self%PMult_t_ptr => PMult_t_E
        !
    end function TransmitterMT_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine baseDealloc().
    !
    subroutine TransmitterMT_dtor( self )
        implicit none
        !
        type( TransmitterMT_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor TransmitterMT_t:", self%id
        !
        call self%baseDealloc
        !
    end subroutine TransmitterMT_dtor
    !
    !> Calculate e_sol_0, e_sol_1 or e_sens from with ForwardSolver
    !> Depending of the Source%adjoint
    !
    subroutine solveTransmitterMT( self )
        implicit none
        !
        class( TransmitterMT_t ), intent( inout ) :: self
        !
        integer :: i_pol
        !
        if( .NOT. allocated( self%source ) ) then
            stop "Error: solveTransmitterMT > source not allocated!"
        endif
        !
        !> First allocate e_sol_0 or e_sens, according to the Source case
        if( self%source%calc_sens ) then
            !
            if( allocated( self%e_sens ) ) deallocate( self%e_sens )
            allocate( cVector3D_SG_t :: self%e_sens(2) )
            !
        else
            !
            !> 
            if( self%i_sol == 0 ) then
                !
                if( allocated( self%e_sol_0 ) ) deallocate( self%e_sol_0 )
                allocate( cVector3D_SG_t :: self%e_sol_0(2) )
                !
            else
                !
                if( allocated( self%e_sol_1 ) ) deallocate( self%e_sol_1 )
                allocate( cVector3D_SG_t :: self%e_sol_1(2) )
                !
            endif
            !
        endif
        !
        !> Calculate e_sol_0 or e_sens through ForwardSolver
        !> For all polarizations (MT n_pol = 2)
        do i_pol = 1, self%n_pol
            !
            !> Verbose
            if( self%source%calc_sens ) then
                !
                write( *, "( a44, es10.2, a6, i2 )" ) "- Solving MT e_sens Tx for period=", self%period, ", pol=", i_pol
                !
                call self%forward_solver%createESolution( i_pol, self%source, self%e_sens( i_pol ) )
                !
            else
                !
                if( self%i_sol == 0 ) then
                    !
                    write( *, "( a42, es10.2, a6, i2 )" ) "- Solving MT e_sol_0 for period=", self%period, ", pol=", i_pol
                    !
                    call self%forward_solver%createESolution( i_pol, self%source, self%e_sol_0( i_pol ) )
                    !
                else
                    !
                    write( *, "( a42, es10.2, a6, i2 )" ) "- Solving MT e_sol_1 for period=", self%period, ", pol=", i_pol
                    !
                    call self%forward_solver%createESolution( i_pol, self%source, self%e_sol_1( i_pol ) )
                    !
                endif
                !
            endif
            !
        enddo
        !
    end subroutine solveTransmitterMT
    !
    !> No subroutine briefing
    !
    function isEqualTransmitterMT( self, other ) result( equal )
        implicit none
        !
        class( TransmitterMT_t ), intent( in ) :: self
        class( Transmitter_t ), intent( in ) :: other
        !
        logical :: equal
        !
        equal = .FALSE.
        !
        select type( other )
            !
            class is( TransmitterMT_t )
                !
                if( ABS( self%period - other%period ) < TOL6 ) then
                    equal = .TRUE.
                endif
                !
            class default
                equal = .FALSE.
            !
        end select
        !
    end function isEqualTransmitterMT
    !
    !> No subroutine briefing
    !
    subroutine printTransmitterMT( self )
        implicit none
        !
        class( TransmitterMT_t ), intent( in ) :: self
        !
        integer :: iRx
        !
        write( *, "( A30, I8, A9, es16.5, A6, I8)" ) &
        "TransmitterMT", self%i_tx, &
        ", Period=",    self%period, &
        ", NRx=", size( self%receiver_indexes )
        !
    end subroutine printTransmitterMT
    !
end module TransmitterMT
