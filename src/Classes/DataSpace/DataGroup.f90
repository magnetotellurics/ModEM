! *************
! 
! Base class to define a Data Group
!
! *************
! 
module DataGroup
    !
    use String
    !
    use Constants
    !
    type :: DataGroup_t
        !
        integer :: id, n_data, id_rx, id_tx
        type( String_t ), allocatable, dimension(:) :: components
        real( kind=prec ), dimension(:), allocatable :: reals, imaginaries, errors ! All Together
        !
        integer, private :: counter
        !
    contains
        !
        final :: DataGroup_dtor
        !
        procedure, public :: add => addDataDg
        !
        procedure, public :: isEqual => isEqualDg
        !
        procedure, public :: print => printDataGroup
        !
    end type DataGroup_t
    !
    interface DataGroup_t
        module procedure DataGroup_ctor
    end interface DataGroup_t
    !
contains
    !
    ! Parametrized constructor
    function DataGroup_ctor( id_rx, id_tx, n_data ) result( self )
        implicit none
        !
        type( DataGroup_t ) :: self
        !
        integer, intent( in ) :: id_rx, id_tx, n_data
        !
        self%id_rx = id_rx
        self%id_tx = id_tx
        self%n_data = n_data
        !
        self%counter = 1
        !
        allocate( self%components( n_data ) )
        !
        allocate( self%reals( n_data ) )
        allocate( self%imaginaries( n_data ) )
        allocate( self%errors( n_data ) )
        !
    end function DataGroup_ctor
    !
    subroutine DataGroup_dtor( self )
        implicit none
        !
        type( DataGroup_t ), intent( in out ) :: self
        integer :: i, asize
        !
        !write( *, * ) "Destructor DataGroup_t: ", self%id
        !
        !asize = size( self%components )
        !do i = asize, 1, -(1)
            !deallocate( self%components(i)%str )
        !enddo
        !deallocate( self%components )
        !
        deallocate( self%reals )
        deallocate( self%imaginaries )
        deallocate( self%errors )
        !
    end subroutine DataGroup_dtor
    !
    subroutine addDataDg( self, component, rvalue, imaginary, error )
        implicit none
        !
        class( DataGroup_t ), intent( inout )   :: self
        character(:), allocatable, intent( in ) :: component
        real( kind=prec ), intent( in )         :: rvalue, imaginary, error
        !
        self%components( self%counter )%str = component
        !
        self%reals( self%counter ) = rvalue
        !
        self%imaginaries( self%counter ) = imaginary
        !
        self%errors( self%counter ) = error
        !
        self%counter = self%counter + 1
        !
    end subroutine addDataDg
    !
    function isEqualDg( self, other ) result ( equal )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: self
        class( DataGroup_t ), intent( in ) :: other
        logical :: equal
        !
        integer :: i_data
        !
        equal = .TRUE.
        !
        if( self%id_tx /= other%id_tx .OR. self%id_rx /= other%id_rx ) then
            equal = .FALSE.
            return
        end if
        !
        do i_data = 1, self%n_data
            !
            if( self%components( i_data )%str /= other%components( i_data )%str .OR. &
                ABS( self%reals( i_data ) - other%reals( i_data ) ) >= TOL6 .OR. &
                ABS( self%imaginaries( i_data ) - other%imaginaries( i_data ) ) >= TOL6 .OR. &
                ABS( self%errors( i_data ) - other%errors( i_data ) ) >= TOL6 ) then
                    equal = .FALSE.
                    exit
            end if
            !
        enddo
        !
    end function isEqualDg
    !
    subroutine printDataGroup( self )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: self
        integer                            :: i_data
        !
        write( *, * ) "    Write DataGroup_t Id: ", self%id
        write( *, * ) "             Receiver Id: ", self%id_rx
        write( *, * ) "          Transmitter Id: ", self%id_tx
        write( *, * ) self%n_data, " data_rows:"
        !
        do i_data = 1, self%n_data
            !
            write( *, * ) i_data, ":", self%components( i_data )%str, self%reals( i_data ), self%imaginaries( i_data ), self%errors( i_data )
            !
        enddo
        !
    end subroutine printDataGroup
    !
end module DataGroup
