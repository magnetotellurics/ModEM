!*************
!
! Base class to define a Receiver
!
!*************
!
module Receiver
    !
    use String
    !
    use Transmitter
    use cVector3D_SG
    use cSparseVector3D_SG
    use ModelOperator
    use DataHandleFArray
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
        type( String_t ), allocatable, dimension(:) :: EHxy, comp_names
        !
        complex( kind=prec ), allocatable, dimension(:) :: response
        !
        type( cSparseVector3D_SG_t ), allocatable :: Lex, Ley, Lez, Lbx, Lby, Lbz
        !
        type( Dh_t ), allocatable, dimension(:) :: predicted_data
        !
        contains
            !
            ! Base interfaces
            procedure( interface_is_equal_rx ), deferred, public :: isEqualRx
            !
            procedure( interface_predicted_data ), deferred, public :: predictedData
            !
            procedure( interface_save_predicted_data ), deferred, public :: savePredictedData
            !
            procedure( interface_write_rx ), deferred, public :: write
            !
            ! Class procedures
            procedure, public :: evaluationFunction => evaluationFunctionRx
            !
            procedure, public :: init    => initializeRx
            procedure, public :: dealloc => deallocateRx
            !
    end type Receiver_t
    !
    public :: getStringReceiverType, getIntReceiverType
    !
    abstract interface
        !
        subroutine interface_predicted_data( self, transmitter )
            !
            import :: Receiver_t, Transmitter_t
            !
            class( Receiver_t ), intent( inout )  :: self
            class( Transmitter_t ), intent( in )  :: transmitter
            !
        end subroutine interface_predicted_data
        !
        subroutine interface_save_predicted_data_rx( self, tx )
            !
            import :: Receiver_t, Transmitter_t
            !
            class( Receiver_t ), intent( in )    :: self
            class( Transmitter_t ), intent( in ) :: tx
            !
        end subroutine interface_save_predicted_data_rx
        !
        function interface_is_equal_rx( self, other ) result( equal )
            !
            import :: Receiver_t
            !
            class( Receiver_t ), intent( in ) :: self, other
            logical :: equal
            !
        end function interface_is_equal_rx
        !
        subroutine interface_write_predicted_data_rx( self )
            !
            import :: Receiver_t
            !
            class( Receiver_t ), intent( in ) :: self
            !
        end subroutine interface_write_predicted_data_rx
        !
        subroutine interface_save_predicted_data( self, tx )
            !
            import :: Receiver_t, Transmitter_t
            !
            class( Receiver_t ), intent( inout ) :: self
            class( Transmitter_t ), intent( in ) :: tx
            !
        end subroutine interface_save_predicted_data
        !
        subroutine interface_write_rx( self )
            !
            import :: Receiver_t
            !
            class( Receiver_t ), intent(in) :: self
            !
        end subroutine interface_write_rx
        !
    end interface
    !
contains
    !
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
        self%location = 0.0
        !
        self%is_complex = .FALSE.
        !
        self%interpolation_set = .FALSE.
        !
        self%Lex = cSparsevector3D_SG_t()
        self%Ley = cSparsevector3D_SG_t()
        self%Lez = cSparsevector3D_SG_t()
        self%Lbx = cSparsevector3D_SG_t()
        self%Lby = cSparsevector3D_SG_t()
        self%Lbz = cSparsevector3D_SG_t()
        !
    end subroutine initializeRx
    !
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
        asize = size( self%EHxy )
        do i = asize, 1, -(1)
            deallocate( self%EHxy(i)%str )
        enddo
        deallocate( self%EHxy )
        !
        if( allocated( self%predicted_data ) ) call deallocateDataHandleArray( self%predicted_data )
        !
        if( allocated( self%Lex ) ) deallocate( self%Lex )
		if( allocated( self%Ley ) ) deallocate( self%Ley )
		if( allocated( self%Lez ) ) deallocate( self%Lez )
		!
		if( allocated( self%Lbx ) ) deallocate( self%Lbx )
		if( allocated( self%Lbx ) ) deallocate( self%Lbx )
		if( allocated( self%Lbz ) ) deallocate( self%Lbz )
		!
    end subroutine deallocateRx
    !
    subroutine evaluationFunctionRx( self, model_operator )
        implicit none
        !
        class( Receiver_t ), intent( inout )   :: self
        class( ModelOperator_t ), intent( in ) :: model_operator
        class( cVector_t ), allocatable        :: temp_full_vec
        !
        integer              :: k
        !
        class( cVector_t ), allocatable :: e, h, lh
        !
        !
        do k = 1, size( self%EHxy )
            !
            select case( self%EHxy(k)%str )
                !
                case( "Ex" )
                    !
                    select type( grid => model_operator%metric%grid )
                        class is( Grid3D_SG_t )
                            allocate( e, source = cVector3D_SG_t( grid, EDGE ) )
                        class default
                            stop "evaluationFunctionRx: Unclassified grid for ex"
                    end select
                    !
                    call e%interpFunc( self%location, "x", temp_full_vec )
                    !
                    select type( temp_full_vec )
                        class is( cVector3D_SG_t )
                            !
                            call full2Sparse( self%Lex, temp_full_vec )
                            !
                        class default
                            stop "evaluationFunctionRx: Unclassified temp_full_vec_ex"
                    end select
                    !
                    deallocate( e )
                    !
                case( "Ey" )
                    !
                    select type( grid => model_operator%metric%grid )
                        class is( Grid3D_SG_t )
                            allocate( e, source = cVector3D_SG_t( grid, EDGE ) )
                        class default
                            stop "evaluationFunctionRx: Unclassified grid for ey"
                    end select
                    !
                    call e%interpFunc( self%location, "y", temp_full_vec )
                    !
                    select type( temp_full_vec )
                        class is( cVector3D_SG_t )
                            !
                            call full2Sparse( self%ley, temp_full_vec )
                            !
                        class default
                            stop "evaluationFunctionRx: Unclassified temp_full_vec_ey"
                    end select
                    !
                    deallocate( e )
                    !
                case( "Ez" )
                    !
                    select type( grid => model_operator%metric%grid )
                        class is( Grid3D_SG_t )
                            allocate( e, source = cVector3D_SG_t( grid, EDGE ) )
                        class default
                            stop "evaluationFunctionRx: Unclassified grid for ez"
                    end select
                    !
                    call e%interpFunc( self%location, "z", temp_full_vec )
                    !
                    select type( temp_full_vec )
                        class is( cVector3D_SG_t )
                            !
                            call full2Sparse( self%Lez, temp_full_vec )
                            !
                        class default
                            stop "evaluationFunctionRx: Unclassified temp_full_vec_ez"
                    end select
                    !
                    deallocate( e )
                    !
                case( "Bx" )
                    !
                     select type( grid => model_operator%metric%grid )
                          class is( Grid3D_SG_t )
                             allocate( h, source = cVector3D_SG_t( grid, FACE ) )
                          class default
                             stop "evaluationFunctionRx: Unclassified grid for hx"
                     end select
                    !
                    call h%interpFunc( self%location, "x", lh )
                    !
                    deallocate( h )
                    !
                    select type( lh )
                        class is( cVector3D_SG_t )
                            if( allocated( temp_full_vec ) ) deallocate( temp_full_vec )
                            allocate( temp_full_vec, source = cVector3D_SG_t( lh%grid, EDGE ) )
                            !
                        class default
                            write( *, * ) "ERROR:Receiver::evaluationFunction:"
                            stop          "            Unknow lh type"
                    end select
                    !
                    call model_operator%multCurlT( lh, temp_full_vec )
                    !
                    deallocate( lh )
                    !
                    !call temp_full_vec%mults( isign * comega )
                    !
                    select type( temp_full_vec )
                        class is( cVector3D_SG_t )
                            !
                            call full2Sparse( self%Lbx, temp_full_vec )
                            !
                        class default
                            stop "evaluationFunctionRx: Unclassified temp_full_vec_bx"
                    end select
                    !
                case( "By" )
                    !
                     select type( grid => model_operator%metric%grid )
                          class is( Grid3D_SG_t )
                             allocate( h, source = cVector3D_SG_t( grid, FACE ) )
                          class default
                             stop "evaluationFunctionRx: Unclassified grid for hy"
                     end select
                    ! 
                    call h%interpFunc( self%location, "y", lh )
                    !
                    deallocate( h )
                    !
                    select type( lh )
                        class is( cVector3D_SG_t )
                            if( allocated( temp_full_vec ) ) deallocate( temp_full_vec )
                            allocate( temp_full_vec, source = cVector3D_SG_t( lh%grid, EDGE ) )
                            !
                        class default
                            write( *, * ) "ERROR:Receiver::evaluationFunction:"
                            stop          "            Unknow lh type"
                    end select
                    !
                    call model_operator%multCurlT( lh, temp_full_vec )
                    !
                    deallocate( lh )
                    !
                    !call temp_full_vec%mults( isign * comega )
                    !
                    select type( temp_full_vec )
                        class is( cVector3D_SG_t )
                            !
                            call full2Sparse( self%Lby, temp_full_vec )
                            !
                        class default
                            stop "evaluationFunctionRx: Unclassified temp_full_vec_by"
                    end select
                    !
                case( "Bz" )
                    !
                     select type( grid => model_operator%metric%grid )
                          class is( Grid3D_SG_t )
                             allocate( h, source = cVector3D_SG_t( grid, FACE ) )
                          class default
                             stop "evaluationFunctionRx: Unclassified grid for hz"
                     end select
                    !
                    call h%interpFunc( self%location, "z", lh )
                    !
                    deallocate( h )
                    !
                    select type( lh )
                        class is( cVector3D_SG_t )
                            if( allocated( temp_full_vec ) ) deallocate( temp_full_vec )
                            allocate( temp_full_vec, source = cVector3D_SG_t( lh%grid, EDGE ) )
                            !
                        class default
                            write( *, * ) "ERROR:Receiver::evaluationFunction:"
                            stop          "            Unknow lh type"
                    end select
                    !
                    call model_operator%multCurlT( lh, temp_full_vec )
                    !
                    deallocate( lh )
                    !
                    !call temp_full_vec%mults( isign * comega )
                    !
                    select type( temp_full_vec )
                        class is( cVector3D_SG_t )
                            !
                            call full2Sparse( self%Lbz, temp_full_vec )
                            !
                        class default
                            stop "evaluationFunctionRx: Unclassified temp_full_vec_bz"
                    end select
                    !
            end select
            !
            deallocate( temp_full_vec )
            !
        end do
        !
    end subroutine evaluationFunctionRx
    !
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
                write( *, * ) "unknow receiver type :[", int_receiver_type, "]"
                STOP "Receiver.f08: getStringReceiverType()"
            !
        end select
        !
    end function getStringReceiverType
    !
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
                write( *, * ) "unknow receiver type :[", str_receiver_type, "]"
                STOP "Receiver.f08: getIntReceiverType()"
            !
        end select
        !
    end function getIntReceiverType
    !
end module Receiver
