!
!> Abstract Base class to define a data file
!
module DataFile
    !
    use Constants
    use String
    !
    use DataEntryArray
    use DataEntryMT
    use DataEntryMT_REF
    use DataEntryCSEM
    !
    use ReceiverFullImpedance
    use ReceiverFullVerticalMagnetic
    use ReceiverOffDiagonalImpedance
    use ReceiverSingleField
    use ReceiverArray
    !
    use TransmitterMT
    use TransmitterCSEM
    use TransmitterArray
    !
    use DataGroupTxArray
    !
    type, abstract :: DataFile_t
        !
        integer :: n_tx, n_rx
        !
        character(:), allocatable :: file_name
        !
        logical :: conjugate
        !
        real( kind=prec ) :: SI_factor
        !
        character(20) :: units, units_in_file
        !
        type( De_t ), allocatable, dimension(:) :: data_entries
        !
        type( DataGroup_t ), allocatable, dimension(:) :: measured_data
        !
        contains
            !
            procedure, public :: init => initializeDataFile
            !
            procedure, public :: dealloc => deallocateDataFile
            !
            procedure, public :: loadReceiversAndTransmitters
            !
            procedure, public :: contructMeasuredDataGroupTxArray
            !
    end type DataFile_t
    !
contains
    !
    !> No subroutine briefing
    subroutine initializeDataFile( self )
        implicit none
        !
        class( DataFile_t ), intent( inout ) :: self
        !
        self%n_tx = 0
        !
        self%n_rx = 0
        !
        self%file_name = ""
        !
        self%conjugate = .FALSE.
        !
        self%SI_factor = ONE
        !
        self%units = "[]"
        !
        self%units_in_file = "[]"
        !
    end subroutine initializeDataFile
    !
    !> No subroutine briefing
    subroutine deallocateDataFile( self )
        implicit none
        !
        class( DataFile_t ), intent( inout ) :: self
        !
        deallocate( self%data_entries )
        !
        deallocate( self%measured_data )
        !
    end subroutine deallocateDataFile
    !
    !> Procedure loadReceiversAndTransmitters
    !> Load all Receivers (Based on Location and Component type) and all Transmitters (Based on Period)
    subroutine loadReceiversAndTransmitters( self, data_entry )
        implicit none
        !
        class( DataFile_t ), intent( inout ) :: self
        class( DataEntry_t ), intent( in ) :: data_entry
        !
        class( Receiver_t ), allocatable :: receiver
        class( Transmitter_t ), pointer :: transmitter
        class( DataGroup_t ), pointer :: data_group
        integer :: i_tx, n_tx, i_rx, rx_type, i_dg
        real( kind=prec ) :: azimuth
        complex( kind=prec ) :: c_value
        !
        call updateDataEntryArray( self%data_entries, data_entry )
        !
        ! TRANSMITTERS
        select type ( data_entry )
            !
            class is ( DataEntryMT_t )
                !
                i_tx = updateTransmitterArray( TransmitterMT_t( data_entry%period ) )
                !
                self%units = "[V/m]/[T]"
                !
            class is ( DataEntryMT_REF_t )
                !
                i_tx = updateTransmitterArray( TransmitterMT_t( data_entry%period ) )
                !
                self%units = "[V/m]/[T]"
                !
            class is ( DataEntryCSEM_t )
                !
                i_tx = updateTransmitterArray( TransmitterCSEM_t( data_entry%period, data_entry%tx_location, data_entry%azimuth, data_entry%dip, data_entry%moment, data_entry%dipole ) )
                !
                self%units = "[V/m]"
                !
        end select
        !
        rx_type = getIntReceiverType( data_entry%type )
        !
        select case( data_entry%type )
            !
            case( "Ex_Field" )
                !
                azimuth = 1.0
                allocate( receiver, source = ReceiverSingleField_t( data_entry%location, azimuth, rx_type ) )
                !
            case( "Ey_Field" )
                !
                azimuth = 2.0
                allocate( receiver, source = ReceiverSingleField_t( data_entry%location, azimuth, rx_type ) )
                !
            case( "Bx_Field" )
                !
                azimuth = 3.0
                allocate( receiver, source = ReceiverSingleField_t( data_entry%location, azimuth, rx_type ) )
                !
            case( "By_Field" )
                !
                azimuth = 4.0
                allocate( receiver, source = ReceiverSingleField_t( data_entry%location, azimuth, rx_type ) )
                !
            case( "Bz_Field" )
                !
                azimuth = 5.0
                allocate( receiver, source = ReceiverSingleField_t( data_entry%location, azimuth, rx_type ) )
                !
            case( "Full_Impedance" )
                !
                allocate( receiver, source = ReceiverFullImpedance_t( data_entry%location, rx_type ) )
                !
            case( "Full_Interstation_TF" )
                !
                stop "Error: DataManager.f08: loadReceiversAndTransmitters(): To implement Full_Interstation_TF !!!!"
                !
            case( "Off_Diagonal_Rho_Phase" )
                !
                stop "Error: DataManager.f08: loadReceiversAndTransmitters(): To implement Off_Diagonal_Rho_Phase !!!!"
                !
            case( "Phase_Tensor" )
                !
                stop "Error: DataManager.f08: loadReceiversAndTransmitters(): To implement Phase_Tensor !!!!"
                !
            case( "Off_Diagonal_Impedance" )
                !
                allocate( receiver, source = ReceiverOffDiagonalImpedance_t( data_entry%location, rx_type ) )
                !
            case( "Full_Vertical_Components", "Full_Vertical_Magnetic" )
                !
                allocate( receiver, source = ReceiverFullVerticalMagnetic_t( data_entry%location, rx_type ) )
                !
            case default
                write( *, * ) "Unknown component type :[", data_entry%type, "]"
                stop "Error: DataManager.f08: loadReceiversAndTransmitters()"
            !
        end select
        !
        receiver%is_complex = data_entry%isComplex()
        !
        receiver%code = data_entry%code
        !
        i_rx = updateReceiverArray( receiver )
        !
        data_group => null()
        !
        if( allocated( self%measured_data ) ) then
            !
            do i_dg = 1, size( self%measured_data )
                !
                if( self%measured_data( i_dg )%i_rx == i_rx .AND. self%measured_data( i_dg )%i_tx == i_tx ) then
                    !
                    data_group => getDataByIndex( self%measured_data, i_dg )
                    !
                    exit
                    !
                endif
                !
            enddo
            !
        endif
        !
        !> 
        if( self%conjugate ) then
            !
            c_value = cmplx( data_entry%real, -data_entry%imaginary, kind=prec )
        else
            !
            c_value = cmplx( data_entry%real, data_entry%imaginary, kind=prec )
        endif
        !
        self%SI_factor = ImpUnits( self%units_in_file, self%units )
        !
        c_value = c_value * self%SI_factor
        !
        if( associated( data_group ) ) then
            !
            call data_group%put( data_entry%component, real( c_value, kind=prec ), real( aimag( c_value ), kind=prec ), data_entry%error )
            !
        else
            !
            allocate( data_group, source = DataGroup_t( i_rx, i_tx, receiver%n_comp, .TRUE. ) )
            !
            call data_group%put( data_entry%component, real( c_value, kind=prec ), real( aimag( c_value ), kind=prec ), data_entry%error )
            !
            call updateDataGroupArray( self%measured_data, data_group )
            !
            deallocate( data_group )
            !
        endif
        !
        deallocate( receiver )
        !
        n_tx = size( transmitters )
        !
        !> Loop over transmitters
        do i_tx = 1, n_tx
            !
            transmitter => getTransmitter( i_tx )
            !
            select type( transmitter )
                !
                class is( TransmitterMT_t )
                    !
                    select type( data_entry )
                        !
                        class is( DataEntryMT_t )
                            !
                            if( ABS( transmitter%period - data_entry%period ) < TOL6 ) then
                                !
                                call transmitter%updateReceiverIndexesArray( i_rx )
                                !
                                exit
                                !
                            endif
                            !
                    end select
                    !
                class is( TransmitterCSEM_t )
                    !
                    select type( data_entry )
                        !
                        class is( DataEntryCSEM_t )
                            !
                            if( ABS( transmitter%period - data_entry%period ) < TOL6      .AND.   &
                                     transmitter%location(1) == data_entry%tx_location(1) .AND.   &
                                     transmitter%location(2) == data_entry%tx_location(2) .AND.   &
                                     transmitter%location(3) == data_entry%tx_location(3) ) then
                                !
                                call transmitter%updateReceiverIndexesArray( i_rx )
                                !
                                exit
                                !
                            endif
                            !
                    end select
                !
                class default
                    stop "Error: DataFile.f90 > loadReceiversAndTransmitters > unclassified Transmitter"
                !
            end select
            !
        enddo
        !
    end subroutine loadReceiversAndTransmitters
    !
    !> Procedure getLineNumber
    !> Return the number of lines of a given file
    function getLineNumber( funit ) result( line_counter )
        implicit none
        !
        integer, intent( in ) :: funit
        !
        integer :: line_counter
        character(80) :: line_text
        !
        line_counter = 0
        !
        rewind( funit )
        !
        do
            read( funit, "(a)", END = 10 ) line_text
            line_text = adjustl( line_text )
            if( index( line_text, "#" ) .EQ. 0 ) then
                !
                line_counter = line_counter + 1
            endif
        enddo
        !
10     return
        !
    end function getLineNumber
    !
    !> Returns a pointer, allowing modifications directly to a DataGroup at a given index
    function getDataByIndex( data_array, i_data ) result( data_group )
        implicit none
        !
        type( DataGroup_t ), dimension(:), target, intent( in ) :: data_array
        integer, intent( in ) :: i_data
        !
        type( DataGroup_t ), pointer :: data_group
        !
        data_group => data_array( i_data )
        !
    end function getDataByIndex
    !
    !> Dynamically add a new DataGroup to the array, always via reallocation.
    subroutine updateDataGroupArray( data_array, data_group )
        implicit none
        !
        type( DataGroup_t ), allocatable, dimension(:), intent( inout ) :: data_array
        type( DataGroup_t ), intent( inout ) :: data_group
        !
        integer :: n_data
        !
        type( DataGroup_t ), allocatable, dimension(:) :: temp_array
        !
        if( .NOT. allocated( data_array ) ) then
            !
            data_group%i_dg = 1
            !
            allocate( data_array(1) )
            !
            data_array(1) = data_group
            !
        else
            !
            n_data = size( data_array )
            !
            allocate( temp_array( n_data + 1 ) )
            !
            temp_array( 1 : n_data ) = data_array(:)
            !
            data_group%i_dg = n_data + 1
            !
            temp_array( n_data + 1 ) = data_group
            !
            deallocate( data_array )
            !
            allocate( data_array, source = temp_array )
            !
            deallocate( temp_array )
            !
        endif
        !
    end subroutine updateDataGroupArray
    !
    !> Allocate and Initialize the predicted data array, 
    !> according to the arrangement of the Transmitter-Receiver pairs of the input
    subroutine contructMeasuredDataGroupTxArray( self )
        implicit none
        !
        class( DataFile_t ), intent( inout ) :: self
        !
        !> Auxiliary variable to group data under a single transmitter index
        type( DataGroupTx_t ) :: tx_data
        !
        !> Local indexes
        integer :: i_tx, i_data
        !
        !> If is the first time, create the predicted data array, 
        !> according to the arrangement of the Transmitter-Receiver pairs of the input
        if( .NOT. allocated( all_measured_data ) ) then
            !
            write( *, * ) "     - Create All Measured Data"
            !
            !> Create an array of DataGroupTx to store the predicted data in the same format as the measured data
            !> Enabling the use of grouped predicted data in future jobs (all_measured_data)
            do i_tx = 1, size( transmitters )
                !
                tx_data = DataGroupTx_t( i_tx )
                !
                do i_data = 1, size( self%measured_data )
                    !
                    if( self%measured_data( i_data )%i_tx == i_tx ) then
                        !
                        call tx_data%put( self%measured_data( i_data ) )
                        !
                    endif
                enddo
                !
                call updateDataGroupTxArray( all_measured_data, tx_data )
                !
            enddo
            !
        else
            stop "Error: contructMeasuredDataGroupTxArray > Unnecessary recreation of predicted data array"
        endif
        !
    end subroutine contructMeasuredDataGroupTxArray
    !
    !> ????
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
        end if
        !
        ! first convert the old units to [V/m]/[T]
        if( index( oldUnits, "[V/m]/[T]" ) > 0 ) then
            ! SI units for E/B
            factor1 = ONE
        else if( index( oldUnits, "[mV/km]/[nT]" ) > 0 ) then
            ! practical units for E/B
            factor1 = ONE * 1000.0
        else if( index( oldUnits, "[V/m]/[A/m]" ) > 0 .OR. index( oldUnits, "Ohm" ) > 0 ) then
            ! SI units for E/H
            factor1 = ONE * 1000.0 * 10000.0 / ( 4 * PI ) ! approx. 796000.0
        else if( index( oldUnits, "[V/m]" ) > 0 ) then
            ! SI units for E
            factor1 = ONE
        else if( index( oldUnits, "[T]" ) > 0 ) then
            ! SI units for B
            factor1 = ONE
        else
            stop "Error: Unknown input units in ImpUnits: "//trim( oldUnits )
        end if
        !
        ! now convert [V/m]/[T] to the new units
        if( index( newUnits, "[V/m]/[T]" ) > 0 ) then
            ! SI units for E/B
            factor2 = ONE
        else if( index( newUnits, "[mV/km]/[nT]" ) > 0 ) then
            ! practical units for E/B
            factor2 = ONE / 1000.0
        else if( index( newUnits, "[V/m]/[A/m]") > 0 .OR. index( newUnits, "Ohm" ) > 0 ) then
            ! SI units for E/H
            factor2 = ONE / ( 1000.0 * 10000.0 / ( 4 * PI ) )
        else if(index( newUnits, "[V/m]") > 0 ) then
            ! SI units for E
            factor2 = ONE
        else if( index( newUnits, "[T]") > 0 ) then
            ! SI units for B
            factor2 = ONE
        else
            stop "Error: Unknown output units in ImpUnits: "//trim( newUnits )
        end if
        !
        SI_factor = factor1 * factor2
        !
    end function ImpUnits
    !
end module DataFile
