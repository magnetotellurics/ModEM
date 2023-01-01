!*************
!>
!> Abstract Base class to define a Receiver
!
module Receiver
    !
    use String
    !
    use Transmitter
    use cVector3D_SG
    use cSparseVector3D_SG
    use ModelOperator
    use DataGroup
    !
    type, abstract :: Receiver_t
        !
        integer :: id, rx_type, n_comp
        !
        character(:), allocatable :: code
        !
        real( kind=prec ), dimension(3) :: location
        !
        logical :: is_complex, interpolation_set
        !
        complex( kind=prec ), allocatable, dimension(:) :: response
        !
        complex( kind=prec ), allocatable, dimension(:,:) :: I_BB
        !
        type( cSparseVector3D_SG_t ) :: Lex, Ley, Lez, Lbx, Lby, Lbz
        !
        type( DataGroup_t ) :: data_group
        !
        class( Vector_t ), allocatable, dimension(:,:) :: lrows
        !
        type( String_t ), allocatable, dimension(:) :: EHxy, comp_names
        !
        contains
            !
            procedure( interface_set_lrows_receiver ), deferred, public :: setLRows
            !
            procedure( interface_predicted_data_receiver ), deferred, public :: predictedData
            !
            procedure( interface_is_equal_receiver ), deferred, public :: isEqualRx
            !
            procedure( interface_print_receiver ), deferred, public :: print
            !
            procedure, public :: evaluationFunction => evaluationFunctionRx
            !
            procedure, public :: savePredictedData => savePredictedDataRx
            !
            procedure, public :: init    => initializeRx
            !
            procedure, public :: dealloc => deallocateRx
            !
    end type Receiver_t
    !
    public :: getStringReceiverType, getIntReceiverType
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_set_lrows_receiver( self, transmitter )
            !
            import :: Receiver_t, Transmitter_t, Vector_t
            !
            class( Receiver_t ), intent( inout ) :: self
            class( Transmitter_t ), intent( in ) :: transmitter
        end subroutine interface_set_lrows_receiver
        !
        !> No interface subroutine briefing
        subroutine interface_predicted_data_receiver( self, transmitter )
            !
            import :: Receiver_t, Transmitter_t
            !
            class( Receiver_t ), intent( inout ) :: self
            class( Transmitter_t ), intent( in ) :: transmitter
            !
        end subroutine interface_predicted_data_receiver
        !
        !> No interface subroutine briefing
        subroutine interface_save_receiver_receiver( self, tx )
            !
            import :: Receiver_t, Transmitter_t
            !
            class( Receiver_t ), intent( in ) :: self
            class( Transmitter_t ), intent( in ) :: tx
            !
        end subroutine interface_save_receiver_receiver
        !
        !> No interface subroutine briefing
        subroutine interface_save_receiver( self, tx )
            !
            import :: Receiver_t, Transmitter_t
            !
            class( Receiver_t ), intent( inout ) :: self
            class( Transmitter_t ), intent( in ) :: tx
            !
        end subroutine interface_save_receiver
        !
        !> No interface function briefing
        function interface_is_equal_receiver( self, other ) result( equal )
            !
            import :: Receiver_t
            !
            class( Receiver_t ), intent( in ) :: self, other
            logical :: equal
            !
        end function interface_is_equal_receiver
        !
        !> No interface subroutine briefing
        subroutine interface_print_receiver( self )
            !
            import :: Receiver_t
            !
            class( Receiver_t ), intent( in ) :: self
            !
        end subroutine interface_print_receiver
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine initializeRx( self )
        implicit none
        !
        class( Receiver_t ), intent( inout ) :: self
        !
        self%id = 0
        !
        self%n_comp = 0
        !
        self%code = ""
        !
        self%rx_type = 0
        !
        self%location = R_ZERO
        !
        self%is_complex = .FALSE.
        !
        self%interpolation_set = .FALSE.
        !
    end subroutine initializeRx
    !
    !> No subroutine briefing
    subroutine deallocateRx( self )
        implicit none
        !
        class( Receiver_t ), intent( inout ) :: self
        !
        integer :: i, asize
        !
        asize = size( self%comp_names )
        do i = asize, 1, -(1)
            deallocate( self%comp_names(i)%str )
        enddo
        deallocate( self%comp_names )
        !
        if( allocated( self%EHxy ) ) then
            asize = size( self%EHxy )
            do i = asize, 1, -(1)
                if( allocated( self%EHxy(i)%str ) ) deallocate( self%EHxy(i)%str )
            enddo
            deallocate( self%EHxy )
        endif
        !
        if( allocated( self%code ) ) deallocate( self%code )
        !
        if( allocated( self%I_BB ) ) deallocate( self%I_BB )
        !
        if( allocated( self%response ) ) deallocate( self%response )
        !
        if( allocated( self%lrows ) ) deallocate( self%lrows )
        !
    end subroutine deallocateRx
    !
    !> No subroutine briefing
    subroutine evaluationFunctionRx( self, model_operator )
        implicit none
        !
        class( Receiver_t ), intent( inout ) :: self
        class( ModelOperator_t ), intent( in ) :: model_operator
        class( Vector_t ), allocatable :: temp_full_vec
        !
        integer :: k
        !
        class( Vector_t ), allocatable :: e_h, lh
        !
        do k = 1, size( self%EHxy )
            !
            select case( self%EHxy(k)%str )
                !
                case( "Ex", "Ey", "Ez" )
                    allocate( e_h, source = cVector3D_SG_t( model_operator%metric%grid, EDGE ) )
                    !
                case( "Bx", "By", "Bz" )
                    allocate( e_h, source = cVector3D_SG_t( model_operator%metric%grid, FACE ) )
                    !
                case default
                    stop "Error: evaluationFunctionRx: Unknown EHxy"
            end select
            !
            select case( self%EHxy(k)%str )
                !
                case( "Ex" )
                    !
                    call e_h%interpFunc( self%location, "x", temp_full_vec )
                    !
                    call self%Lex%fromFullVector( temp_full_vec )
                    !
                case( "Ey" )
                    !
                    call e_h%interpFunc( self%location, "y", temp_full_vec )
                    !
                    call self%Ley%fromFullVector( temp_full_vec )
                    !
                case( "Ez" )
                    !
                    call e_h%interpFunc( self%location, "z", temp_full_vec )
                    !
                    call self%Lez%fromFullVector( temp_full_vec )
                    !
                case( "Bx" )
                    !
                    call e_h%interpFunc( self%location, "x", lh )
                    !
                    allocate( temp_full_vec, source = cVector3D_SG_t( model_operator%metric%grid, EDGE ) )
                    !
                    call model_operator%multCurlT( lh, temp_full_vec )
                    !
                    deallocate( lh )
                    !
                    call self%Lbx%fromFullVector( temp_full_vec )
                    !
                case( "By" )
                    !
                    call e_h%interpFunc( self%location, "y", lh )
                    !
                    allocate( temp_full_vec, source = cVector3D_SG_t( model_operator%metric%grid, EDGE ) )
                    !
                    call model_operator%multCurlT( lh, temp_full_vec )
                    !
                    deallocate( lh )
                    !
                    call self%Lby%fromFullVector( temp_full_vec )
                    !
                case( "Bz" )
                    !
                    call e_h%interpFunc( self%location, "z", lh )
                    !
                    allocate( temp_full_vec, source = cVector3D_SG_t( model_operator%metric%grid, EDGE ) )
                    !
                    call model_operator%multCurlT( lh, temp_full_vec )
                    !
                    deallocate( lh )
                    !
                    call self%Lbz%fromFullVector( temp_full_vec )
                    !
            end select
            !
            deallocate( e_h, temp_full_vec )
            !
        enddo
        !
    end subroutine evaluationFunctionRx
    !
    !> No subroutine briefing
    subroutine savePredictedDataRx( self, transmitter )
        implicit none
        !
        class( Receiver_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        !
        character(:), allocatable :: component
        real( kind=prec ) :: real_part, imaginary, error
        !
        integer :: i
        !
        self%data_group = DataGroup_t( self%id, transmitter%id, self%n_comp )
        !
        do i = 1, self%n_comp
            !
            component = trim( self%comp_names(i)%str )
            real_part = real( self%response(i), kind=prec )
            imaginary = real( aimag( self%response(i) ), kind=prec )
            error = 1.0
            !
            call self%data_group%put( component, real_part, imaginary, error )
            !
        enddo
        !
    end subroutine savePredictedDataRx
    !
    !> No function briefing
    function getStringReceiverType( int_receiver_type ) result( str_receiver_type )
    !
        integer, intent( in ) :: int_receiver_type
        character(:), allocatable :: str_receiver_type
        !
        select case( int_receiver_type )
            !
            case( 1 )
                str_receiver_type = "Full_Impedance"
            case( 2 )
                str_receiver_type = "Full_Interstation_TF"
            case( 3 )
                str_receiver_type = "Off_Diagonal_Rho_Phase"
            case( 4 )
                str_receiver_type = "Phase_Tensor"
            case( 5 )
                str_receiver_type = "Off_Diagonal_Impedance"
            case( 6 )
                str_receiver_type = "Ex_Field"
            case( 7 )
                str_receiver_type = "Ey_Field"
            case( 8 )
                str_receiver_type = "Bx_Field"
            case( 9 )
                str_receiver_type = "By_Field"
            case( 10 )
                str_receiver_type = "Bz_Field"
            case( 11 )
                str_receiver_type = "Full_Vertical_Components"
            case( 12 )
                str_receiver_type = "Full_Vertical_Magnetic"
            case default
                write( *, * ) "Error: Unknown receiver type :[", int_receiver_type, "]"
                stop "Receiver.f08: getStringReceiverType()"
            !
        end select
        !
    end function getStringReceiverType
    !
    !> No function briefing
    function getIntReceiverType( str_receiver_type ) result( int_receiver_type )
        !
        character(:), allocatable, intent( in ) :: str_receiver_type
        integer :: int_receiver_type
        !
        int_receiver_type = 0
        !
        select case( str_receiver_type )
            !
            case( "Full_Impedance" )
                int_receiver_type = 1
            case( "Full_Interstation_TF" )
                int_receiver_type = 2
            case( "Off_Diagonal_Rho_Phase" )
                int_receiver_type = 3
            case( "Phase_Tensor" )
                int_receiver_type = 4
            case( "Off_Diagonal_Impedance" )
                int_receiver_type = 5
            case( "Ex_Field" )
                int_receiver_type = 6
            case( "Ey_Field" )
                int_receiver_type = 7
            case( "Bx_Field" )
                int_receiver_type = 8
            case( "By_Field" )
                int_receiver_type = 9
            case( "Bz_Field" )
                int_receiver_type = 10
            case( "Full_Vertical_Components" )
                int_receiver_type = 11
            case( "Full_Vertical_Magnetic" )
                int_receiver_type = 12
            case default
                write( *, * ) "Error: Unknown receiver type :[", str_receiver_type, "]"
                stop "Receiver.f08: getIntReceiverType()"
            !
        end select
        !
    end function getIntReceiverType
    !
end module Receiver
