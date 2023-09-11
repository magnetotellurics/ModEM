!
!> Abstract Base class to define a Receiver
!
module Receiver
    !
    use String
    use DataGroup
    use Transmitter
    use cVectorSparse3D_SG
    !
    type, abstract :: Receiver_t
        !
        integer :: i_rx, rx_type, n_comp
        !
        character(:), allocatable :: code, units
        !
        real( kind=prec ), dimension(3) :: location
        !
        logical :: is_complex, interpolation_set
        !
        complex( kind=prec ), allocatable, dimension(:) :: response
        !
        complex( kind=prec ), allocatable, dimension(:,:) :: I_BB
        !
        type( cVectorSparse3D_SG_t ) :: Lex, Ley, Lez, Lbx, Lby, Lbz
        !
        type( GenVector_t ), allocatable, dimension(:,:) :: lrows
        !
        type( String_t ), allocatable, dimension(:) :: EHxy, comp_names
        !
        contains
            !
            !> Base interfaces
            procedure( interface_set_lrows_receiver ), deferred, public :: setLRows
            !
            procedure( interface_predicted_data_receiver ), deferred, public :: predictedData
            !
            procedure( interface_is_equal_receiver ), deferred, public :: isEqualRx
            !
            procedure( interface_print_receiver ), deferred, public :: print
            !
            !> Base routines
            procedure, public :: evaluationFunction => evaluationFunction_Receiver
            !
            procedure, public :: savePredictedData => savePredictedData_Receiver
            !
            procedure, public :: baseInit => initialize_Receiver
            !
            procedure, public :: baseDealloc => deallocate_Receiver
            !
    end type Receiver_t
    !
    !> Module routines
    public :: getStringReceiverType, getIntReceiverType
    public :: ImpUnits
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_set_lrows_receiver( self, transmitter )
            import :: Receiver_t, Transmitter_t
            !
            class( Receiver_t ), intent( inout ) :: self
            class( Transmitter_t ), intent( in ) :: transmitter
        end subroutine interface_set_lrows_receiver
        !
        !> No interface subroutine briefing
        subroutine interface_predicted_data_receiver( self, transmitter, data_group )
            import :: Receiver_t, Transmitter_t, DataGroup_t
            !
            class( Receiver_t ), intent( inout ) :: self
            class( Transmitter_t ), intent( in ) :: transmitter
            type( DataGroup_t ), intent( out ), optional :: data_group
        end subroutine interface_predicted_data_receiver
        !
        !> No interface subroutine briefing
        subroutine interface_save_receiver_receiver( self, tx )
            import :: Receiver_t, Transmitter_t
            !
            class( Receiver_t ), intent( in ) :: self
            class( Transmitter_t ), intent( in ) :: tx
        end subroutine interface_save_receiver_receiver
        !
        !> No interface subroutine briefing
        subroutine interface_save_receiver( self, tx )
            import :: Receiver_t, Transmitter_t
            !
            class( Receiver_t ), intent( inout ) :: self
            class( Transmitter_t ), intent( in ) :: tx
        end subroutine interface_save_receiver
        !
        !> No interface function briefing
        function interface_is_equal_receiver( self, other ) result( equal )
            import :: Receiver_t
            !
            class( Receiver_t ), intent( in ) :: self, other
            logical :: equal
        end function interface_is_equal_receiver
        !
        !> No interface subroutine briefing
        subroutine interface_print_receiver( self )
            import :: Receiver_t
            !
            class( Receiver_t ), intent( in ) :: self
        end subroutine interface_print_receiver
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine initialize_Receiver( self )
        implicit none
        !
        class( Receiver_t ), intent( inout ) :: self
        !
        self%i_rx = 0
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
    end subroutine initialize_Receiver
    !
    !> No subroutine briefing
    subroutine deallocate_Receiver( self )
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
    end subroutine deallocate_Receiver
    ! !
    ! !> No subroutine briefing
    ! subroutine deallocateLRows( transmitter )
        ! implicit none
        ! !
        ! class( Transmitter_t ), intent( in ) :: transmitter
        ! integer :: i_pol, i_comp, nrx, irx
        ! !
        ! if( allocated( self%lrows ) ) then
            ! !
            ! asize = size( self%EHxy )
            ! do i_comp = asize, 1, -(1)
                ! if( allocated( self%lrows( i_pol, i_comp )%v ) ) deallocate( self%lrows( i_pol, i_comp )%v )
            ! enddo
            ! !
            ! deallocate( self%EHxy )
            ! !
        ! endif
        ! !
    ! end subroutine deallocateLRows
    ! !
    !> No subroutine briefing
    subroutine evaluationFunction_Receiver( self, model_operator )
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
            !> Create e_h vector with proper type
            !> Eletrical(E) - EDGE and Magnetic(B) - FACE
            select case( self%EHxy(k)%str )
                !
                case( "Ex", "Ey", "Ez" )
                    !
                    call model_operator%metric%createVector( complex_t, EDGE, e_h )
                    !
                case( "Bx", "By", "Bz" )
                    !
                    call model_operator%metric%createVector( complex_t, FACE, e_h )
                    !
                case default
                    call errStop( "evaluationFunction_Receiver > Unknown EHxy" )
                !
            end select
            !
            !> Different behaviour for E and B components
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
                    call model_operator%metric%createVector( complex_t, EDGE, temp_full_vec )
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
                    call model_operator%metric%createVector( complex_t, EDGE, temp_full_vec )
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
                    call model_operator%metric%createVector( complex_t, EDGE, temp_full_vec )
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
    end subroutine evaluationFunction_Receiver
    !
    !> No subroutine briefing
    subroutine savePredictedData_Receiver( self, transmitter, data_group )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        type( DataGroup_t ), intent( out ) :: data_group
        !
        real( kind=prec ) :: real_part, imaginary, error
        !
        integer :: i
        !
        data_group = DataGroup_t( self%i_rx, transmitter%i_tx, self%n_comp, .FALSE. )
        !
        data_group%is_complex = self%is_complex
        !
        do i = 1, self%n_comp
            !
            real_part = real( self%response(i), kind=prec )
            imaginary = real( aimag( self%response(i) ), kind=prec )
            error = R_ZERO
            !
            call data_group%put( real_part, imaginary, error )
            !
        enddo
        !
    end subroutine savePredictedData_Receiver
    !
    !> No subroutine briefing
    !
    function getStringReceiverType( int_receiver_type ) result( str_receiver_type )
    !
        integer, intent( in ) :: int_receiver_type
        character(:), allocatable :: str_receiver_type
        !
        select case( int_receiver_type )
            !
            case(1)
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
                call errStop( "getStringReceiverType > Unknown receiver type" )
            !
        end select
        !
    end function getStringReceiverType
    !
    !> No subroutine briefing
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
                call errStop( "getIntReceiverType > Unknown receiver type :["//str_receiver_type//"]" )
            !
        end select
        !
    end function getIntReceiverType
    !
    !> No subroutine briefing
    !
    function ImpUnits( oldUnits, newUnits ) result( SI_factor )
        implicit none
        !
        character(*), intent( in ) :: oldUnits, newUnits
        !
        real( kind=prec ) :: SI_factor
        !
        real( kind=prec ) :: factor1, factor2
        !
        ! if the quantity is dimensionless, do nothing
        if( index( oldUnits, "[]"  ) > 0 .OR. index( newUnits, "[]" ) > 0 ) then
            SI_factor = ONE
            return
        endif
        !
        ! first convert the old units to [V/m]/[T]
        if( index( oldUnits, "[V/m]/[T]" ) > 0 ) then
            ! SI units for E/B
            factor1 = ONE
        elseif( index( oldUnits, "[mV/km]/[nT]" ) > 0 ) then
            ! practical units for E/B
            factor1 = ONE * 1000.0
        elseif( index( oldUnits, "[V/m]/[A/m]" ) > 0 .OR. index( oldUnits, "Ohm" ) > 0 ) then
            ! SI units for E/H
            factor1 = ONE * 1000.0 * 10000.0 / ( 4 * PI ) ! approx. 796000.0
        elseif( index( oldUnits, "[V/m]" ) > 0 ) then
            ! SI units for E
            factor1 = ONE
        elseif( index( oldUnits, "[T]" ) > 0 ) then
            ! SI units for B
            factor1 = ONE
        else
            call errStop( "ImpUnits > Unknown input units ["//trim( oldUnits )//"]" )
        endif
        !
        ! now convert [V/m]/[T] to the new units
        if( index( newUnits, "[V/m]/[T]" ) > 0 ) then
            ! SI units for E/B
            factor2 = ONE
        elseif( index( newUnits, "[mV/km]/[nT]" ) > 0 ) then
            ! practical units for E/B
            factor2 = ONE / 1000.0
        elseif( index( newUnits, "[V/m]/[A/m]") > 0 .OR. index( newUnits, "Ohm" ) > 0 ) then
            ! SI units for E/H
            factor2 = ONE / ( 1000.0 * 10000.0 / ( 4 * PI ) )
        elseif(index( newUnits, "[V/m]") > 0 ) then
            ! SI units for E
            factor2 = ONE
        elseif( index( newUnits, "[T]") > 0 ) then
            ! SI units for B
            factor2 = ONE
        else
            call errStop( "ImpUnits > Unknown output units ["//trim( newUnits )//"]" )
        endif
        !
        SI_factor = factor1 * factor2
        !
    end function ImpUnits
    !
end module Receiver
