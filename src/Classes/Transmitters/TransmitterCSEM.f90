!
!> Derived class to define a CSEM Transmitter
!
module TransmitterCSEM
    !
    use Transmitter 
    !
    type, extends( Transmitter_t ), public :: TransmitterCSEM_t
        !
        real( kind=prec ) :: location(3), azimuth, dip, moment
        character(:), allocatable :: dipole
        !
        contains
            !
            final :: TransmitterCSEM_dtor
            !
            procedure, public :: solve => solveTransmitterCSEM
            !
            procedure, public :: isEqual => isEqualTransmitterCSEM
            !
            procedure, public :: print => printTransmitterCSEM
            !
    end type TransmitterCSEM_t
    !
    interface TransmitterCSEM_t
        module procedure TransmitterCSEM_ctor
    end interface TransmitterCSEM_t
    !
contains
    !
    !> Parametrized constructor
    function TransmitterCSEM_ctor( period, location, azimuth, dip, moment, dipole ) result ( self )
        !
        type( TransmitterCSEM_t ) :: self
        !
        real( kind=prec ), intent( in ) :: period, azimuth, dip, moment, location(3)
        character(:), allocatable, intent( in ) :: dipole
        !
        !> write( *, * ) "Constructor TransmitterCSEM_t"
        !
        call self%init()
        !
        self%n_pol = 1
        !
        self%period = period
        !
        self%omega = ( 2.0 * PI / self%period )
        !
        self%location = location
        !
        self%azimuth = azimuth
        !
        self%dip = dip
        !
        self%moment = moment
        !
        self%dipole = dipole
        !
        !self%PMult_ptr => PMult_E
        !
        !self%PMult_t_ptr => PMult_t_E
        !
    end function TransmitterCSEM_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    subroutine TransmitterCSEM_dtor( self )
        implicit none
        !
        type( TransmitterCSEM_t ) :: self
        !
        !> write( *, * ) "Destructor TransmitterCSEM_t"
        !
        call self%dealloc()
        !
        deallocate( self%dipole )
        !
    end subroutine TransmitterCSEM_dtor
    !
    !> No function briefing
    function isEqualTransmitterCSEM( self, other ) result( equal )
        implicit none
        !
        class( TransmitterCSEM_t ), intent( in ) :: self
        class( Transmitter_t ), intent( in ) :: other
        !
        logical :: equal
        !
        equal = .FALSE.
        !
        select type( other )
            !
            class is( TransmitterCSEM_t )
                !
                if( ABS( self%period - other%period ) < TOL6 .AND.   &
                         self%location(1) == other%location(1) .AND. &
                         self%location(2) == other%location(2) .AND. &
                         self%location(3) == other%location(3) ) then
                    equal = .TRUE.
                endif
                !
            class default
                equal = .FALSE.
            !
        end select
        !
    end function isEqualTransmitterCSEM
    !
    !> No subroutine briefing
    subroutine solveTransmitterCSEM( self )
        implicit none
        !
        class( TransmitterCSEM_t ), intent( inout ) :: self
        !
        integer :: ios
        !
        character( len=20 ) :: ModeName
        !
        if( .NOT. allocated( self%source ) ) then
            stop "Error: solveTransmitterCSEM > source not allocated!"
        endif
        !
        !> First allocate e_sol or e_sens, according to the Source case
        if( self%source%sens ) then
            !
            if( allocated( self%e_sens ) ) deallocate( self%e_sens )
            allocate( cVector3D_SG_t :: self%e_sens(1) )
            !
        else
            !
            if( allocated( self%e_sol ) ) deallocate( self%e_sol )
            allocate( cVector3D_SG_t :: self%e_sol(1) )
            !
        endif
        !
        !> Calculate e_sol or e_sens through ForwardSolver
        !> For one polarization (CSEM n_pol = 1)
        if( self%source%sens ) then
            !
            call self%forward_solver%createESolution( 1, self%source, self%e_sens(1) )
            !
            call self%e_sens(1)%add( self%source%E(1) )
            !
        else
            !
            call self%forward_solver%createESolution( 1, self%source, self%e_sol(1) )
            !
            call self%e_sol(1)%add( self%source%E(1) )
            !
        endif
        !
    end subroutine solveTransmitterCSEM
    !
    !> No subroutine briefing
    subroutine printTransmitterCSEM( self )
        implicit none
        !
        class( TransmitterCSEM_t ), intent( in ) :: self
        !
        integer :: iRx
        !
        write( *, "( A30, I8, A13, f16.3, f16.3, f16.3, A10, es16.5, A6, I8)" ) &
        "TransmitterCSEM", self%i_tx, &
        ", Location= [", self%location(1), self%location(2), self%location(3), &
        "], Period=",    self%period, &
        ", NRx=", size( self%receiver_indexes )
        !
    end subroutine printTransmitterCSEM
    !
end module TransmitterCSEM
