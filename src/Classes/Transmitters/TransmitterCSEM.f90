! *************
! 
! Derived class to define a CSEM Transmitter
!
! *************
! 
module TransmitterCSEM
    ! 
    use FileUnits
    use Transmitter 
    use cVector3D_SG
    !
    type, extends( Transmitter_t ), public :: TransmitterCSEM_t
        !
        real( kind=prec )         :: location(3), azimuth, dip, moment
        character(:), allocatable :: dipole
        !
        contains
            !
            final    :: TransmitterCSEM_dtor
            !
            procedure, public :: solveFWD => solveFWDTransmitterCSEM
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
    ! Parametrized constructor
    function TransmitterCSEM_ctor( period, location, azimuth, dip, moment, dipole ) result ( self )
        !
        type( TransmitterCSEM_t ) :: self
        !
        real( kind=prec ), intent( in )         :: period, azimuth, dip, moment, location(3)
        character(:), allocatable, intent( in ) :: dipole
        !
        ! write(*,*) "Constructor TransmitterCSEM_t"
        !
        call self%init()
        !
        self%n_pol = 1
        self%period = period
        self%location = location
        self%azimuth = azimuth
        self%dip = dip
        self%moment = moment
        self%dipole = dipole
        !
        !self%pMult_ptr => pMult_E
        !
        !self%pMult_t_ptr => pMult_t_E
        !
    end function TransmitterCSEM_ctor
    !
    ! Destructor
    subroutine TransmitterCSEM_dtor( self )
        implicit none
        !
        type( TransmitterCSEM_t )    :: self
        !
        ! write(*,*) "Destructor TransmitterCSEM_t"
        !
        call self%dealloc()
        !
    end subroutine TransmitterCSEM_dtor
    !
    function isEqualTransmitterCSEM( self, other ) result( equal )
        implicit none
        !
        class( TransmitterCSEM_t ), intent( in ) :: self
        class( Transmitter_t ), intent( in )     :: other
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
    subroutine solveFWDTransmitterCSEM( self )
        implicit none
        !
        class( TransmitterCSEM_t ), intent( inout ) :: self
        !
        integer           :: ios
        real( kind=prec ) :: omega
        !
        character( len=20 ) :: ModeName
        !
        !
        omega = 2.0 * PI / self%period
        !
        allocate( cVector3D_SG_t :: self%e_all( self%n_pol ) )
        !
        ! Verbose...
        write( *, * ) "          SolveFWD for CSEM Tx:", self%id, " -> Period:", self%period
        !
        call self%source%setE( 1 )
        !
        select type( mgrid => self%source%model_operator%metric%grid )
            class is( Grid3D_SG_t )
                !
                self%e_all( 1 ) = cVector3D_SG_t( mgrid, EDGE )
                !
        end select
        !
        call self%forward_solver%getESolution( self%source, self%e_all( 1 ) )
        call self%e_all( 1 )%add( self%source%E )
        !
        ModeName = "Ex"
        !
        open( ioESolution, file = e_solution_file_name, action = "write", position = "append", form = "unformatted", iostat = ios )
        !
        if( ios /= 0 ) then
            stop "Error opening file in solveFWDTransmitterCSEM: e_solution"
        else
            !
            ! write the frequency header - 1 record
            write( ioESolution ) omega, self%id, 1, ModeName
            !
            call self%e_all( 1 )%write( ioESolution )
            !
            close( ioESolution )
            !
        endif
        !
    end subroutine solveFWDTransmitterCSEM
    !
    subroutine printTransmitterCSEM( self )
        implicit none
        !
        class( TransmitterCSEM_t ), intent( in ) :: self
        !
        integer :: iRx
        !
        write( *, "( A30, I5, A12, f8.2, f8.2, f8.2, A10, es10.2, A7, I5)" ) &
        "               TransmitterCSEM", self%id, &
        ": Location[", self%location(1), self%location(2), self%location(3), &
        "], Period: ",    self%period, &
        ", NRx: ", size( self%receiver_indexes )
        !
    end subroutine printTransmitterCSEM
    !
end module TransmitterCSEM
