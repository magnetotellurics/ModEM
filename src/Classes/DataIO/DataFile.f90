!*************
!
! Base class to read a data file
!
!*************
!
module DataFile
    !
    use Constants
    use String
    !
    use DataEntryFArray
    use DataEntryMT
    use DataEntryMT_REF
    use DataEntryCSEM
    !
    use ReceiverFullImpedance
    use ReceiverFullVerticalMagnetic
    use ReceiverOffDiagonalImpedance
    use ReceiverSingleField
    use ReceiverFArray
    use TransmitterMT
    use TransmitterCSEM
    use TransmitterFArray
    !
    use DataGroupFArray
    !
    type, abstract :: DataFile_t
        !
        integer                   :: nTx, nRx
        character(:), allocatable :: file_name
        !
        type( De_t ), allocatable, dimension(:) :: data_entries
        !
        contains
            !
            procedure, public  :: init    => initializeDataFile
            procedure, public  :: dealloc => deallocateDataFile
            procedure, public  :: loadReceiversAndTransmitters
            !
    end type DataFile_t
    !
contains
    !
    subroutine initializeDataFile( self )
        class( DataFile_t ), intent( inout ) :: self
        !
        self%nTx = 0
        self%nRx = 0
        !
    end subroutine initializeDataFile
    !
    subroutine deallocateDataFile( self )
        class( DataFile_t ), intent( inout ) :: self
        !
        deallocate( self%data_entries )
        !
    end subroutine deallocateDataFile
    !
    ! Load all Receivers (Based on Location and Component type) and all Transmitters (Based on Period)
    subroutine loadReceiversAndTransmitters( self, data_entry )
        implicit none
        !
        class( DataFile_t ), intent( inout ) :: self
        class( DataEntry_t ), intent( in )   :: data_entry
        !
        class( Receiver_t ), allocatable :: receiver
        class( Transmitter_t ), pointer  :: transmitter
        class( DataGroup_t ), pointer    :: data_group
        integer                          :: iTx, nTx, rx_id, rx_type, iDg, dg_index
        real( kind=prec )                :: azimuth
        !
        call updateDataEntryArray( self%data_entries, data_entry )
        !
        ! TRANSMITTERS
        !
        select type ( data_entry )
            !
            class is ( DataEntryMT_t )
                !
                iTx = updateTransmitterArray( TransmitterMT_t( data_entry%period ) )
                !
            class is ( DataEntryMT_REF_t )
                !
                iTx = updateTransmitterArray( TransmitterMT_t( data_entry%period ) )
                !
            class is ( DataEntryCSEM_t )
                !
                iTx = updateTransmitterArray( TransmitterCSEM_t( data_entry%period, data_entry%tx_location, data_entry%azimuth, data_entry%dip, data_entry%moment, data_entry%dipole ) )
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
                receiver = ReceiverSingleField_t( data_entry%location, azimuth, rx_type )
                !
            case( "Ey_Field" )
                !
                azimuth = 2.0
                receiver = ReceiverSingleField_t( data_entry%location, azimuth, rx_type )
                !
            case( "Bx_Field" )
                !
                azimuth = 3.0
                receiver = ReceiverSingleField_t( data_entry%location, azimuth, rx_type )
                !
            case( "By_Field" )
                !
                azimuth = 4.0
                receiver = ReceiverSingleField_t( data_entry%location, azimuth, rx_type )
                !
            case( "Bz_Field" )
                !
                azimuth = 5.0
                receiver = ReceiverSingleField_t( data_entry%location, azimuth, rx_type )
                !
            case( "Full_Impedance" )
                !
                receiver = ReceiverFullImpedance_t( data_entry%location, rx_type )
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
                receiver = ReceiverOffDiagonalImpedance_t( data_entry%location, rx_type )
                !
            case( "Full_Vertical_Components", "Full_Vertical_Magnetic" )
                !
                receiver = ReceiverFullVerticalMagnetic_t( data_entry%location, rx_type )
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
        rx_id = updateReceiverArray( receiver )
        !
        ! Set data groups with input data
        data_group => null()
        !
        if( allocated( original_data ) ) then
            !
            do iDg = 1, size( original_data )
                !
                if( original_data( iDg )%id_rx == rx_id .AND. original_data( iDg )%id_tx == iTx ) then
                    !
                    data_group => getDataGroupByIndex( original_data, iDg )
                    !
                    exit
                    !
                endif
                !
            enddo
            !
        endif
        !
        if( associated( data_group ) ) then
            !
            call data_group%add( data_entry%component, data_entry%real, data_entry%imaginary, data_entry%error )
            !
        else
            !
            allocate( data_group, source = DataGroup_t( rx_id, iTx, receiver%n_comp ) )
            !
            call data_group%add( data_entry%component, data_entry%real, data_entry%imaginary, data_entry%error )
            !
            call updateDataGroupArray( original_data, data_group )
            !
            deallocate( data_group )
            !
        endif
        !
        deallocate( receiver )
        !
        nTx = size( transmitters )
        !
        ! Loop over transmitters
        do iTx = 1, nTx
            !
            transmitter => getTransmitter( iTx )
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
                                call transmitter%updateReceiverIndexesArray( rx_id )
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
                                call transmitter%updateReceiverIndexesArray( rx_id )
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
    ! Return the number of lines of a given file
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
            end if
        end do
        !
10     return
        !
    end function getLineNumber
    !
end module DataFile
