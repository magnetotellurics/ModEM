! *************
! 
! Derived class to define a MT Transmitter
!
! *************
! 
module TransmitterMT
    ! 
    use FileUnits
    use cVector3D_SG
    use Transmitter
    !
    type, extends( Transmitter_t ), public :: TransmitterMT_t
        !
        ! PROPERTIES HERE
        !
        contains
            !
            final :: TransmitterMT_dtor
            !
            procedure, public :: solveFWD => solveFWDTransmitterMT
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
    ! TransmitterMT constructor
    function TransmitterMT_ctor( period ) result ( self )
        implicit none
        !
        type( TransmitterMT_t ) :: self
        !
        real( kind=prec ), intent( in ) :: period
        !
        !write(*,*) "Constructor TransmitterMT_t"
        !
        call self%init()
        !
        self%n_pol = 2
        !
        self%period = period
        !
        ! self%pMult_ptr => pMult_E
        !
        ! self%pMult_t_ptr => pMult_t_E
        !
    end function TransmitterMT_ctor
    !
    ! TransmitterMT destructor
    subroutine TransmitterMT_dtor( self )
        implicit none
        !
        type( TransmitterMT_t ), intent( inout ) :: self
        !
        !write(*,*) "Destructor TransmitterMT_t:", self%id
        !
        call self%dealloc()
        !
    end subroutine TransmitterMT_dtor
    !
    ! Set self%e_all from forward modeling solver
    subroutine solveFWDTransmitterMT( self )
        implicit none
        !
        class( TransmitterMT_t ), intent( inout ) :: self
        !
        integer           :: i_pol, ios
        real( kind=prec ) :: omega
        !
        character( len=20 ) :: ModeName
            !
            character(:), allocatable :: title
            !
        !
        omega = 2.0 * PI / self%period
        !
        allocate( cVector3D_SG_t :: self%e_all( self%n_pol ) )
        !
        ! Loop over all polarizations (MT n_pol = 2)
        do i_pol = 1, self%n_pol
            !
            ! Verbose
            write( *, * ) "          SolveFWD for MT Tx:", self%id, " -> Period:", self%period, " - Polarization:", i_pol
            !
            call self%source%setE( i_pol )
            !
            select type( grid => self%source%model_operator%metric%grid )
                class is( Grid3D_SG_t )
                    !
                    self%e_all( i_pol ) = cVector3D_SG_t( grid, EDGE )
                    !
                class default
                    stop "Error: solveFWDTransmitterMT: undefined grid"
                    !
            end select
            !
            call self%forward_solver%getESolution( self%source, self%e_all( i_pol ) )
            !
            if( i_pol == 1 ) then
                ModeName = "Ey"
            else
                ModeName = "Ex"
            endif
            !
            open( ioESolution, file = e_solution_file_name, action = "write", position = "append", form = "unformatted", iostat = ios )
            !
            if( ios /= 0 ) then
                stop "Error opening file in solveFWDTransmitterMT: e_solution"
            else
                !
                ! write the frequency header - 1 record
                write( ioESolution ) omega, self%id, i_pol, ModeName
                !
                call self%e_all( i_pol )%write( ioESolution )
                !
                close( ioESolution )
                !
            endif
        !
        enddo
        !
    end subroutine solveFWDTransmitterMT
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
    subroutine printTransmitterMT( self )
        implicit none
        !
        class( TransmitterMT_t ), intent( in ) :: self
        !
        integer :: iRx
        !
        write( *, "( A29, I5, A10, es10.2, A7, I5)" ) &
        "               TransmitterMT:", self%id, &
        ", Period: ",    self%period, &
        ", NRx: ", size( self%receiver_indexes )
        !
    end subroutine printTransmitterMT
    !
end module TransmitterMT
