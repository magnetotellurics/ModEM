!*************
!
! Base class to read a data file
!
! Last modified at 05/2021 by Paulo Werdt
!
!*************
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
    use ReceiverFArray
    use TransmitterMT
    use TransmitterCSEM
    use TransmitterFArray
    !
    type, abstract :: DataFile_t
        !
        integer                   :: nTx, nRx
        character(:), allocatable :: fine_name
        !
        class( DataEntryArray_t ), pointer :: data_entries
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
        allocate( self%data_entries, source = DataEntryArray_t() )
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
        class( Receiver_t ), pointer    :: receiver
        class( Transmitter_t ), pointer :: transmitter
        integer                         :: iTx, nTx, rx_id, rx_type
        real ( kind=prec )              :: azimuth
        !
        call self%data_entries%add( data_entry )
        !
        ! TRANSMITTERS
        !
        select type ( data_entry )
            !
            class is ( DataEntryMT_t )
                !
                allocate( transmitter, source = TransmitterMT_t( data_entry%period ) )
                !
            class is ( DataEntryMT_REF_t )
                !
                allocate( transmitter, source = TransmitterMT_t( data_entry%period ) )
                !
            class is ( DataEntryCSEM_t )
                !
                allocate( transmitter, source = TransmitterCSEM_t( data_entry%period, data_entry%tx_location, data_entry%azimuth, data_entry%dip, data_entry%moment, data_entry%dipole ) )
                !
        end select
        !
        if( updateTransmitterArray( transmitter ) == 0 ) deallocate( transmitter )
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
                stop "DataManager.f08: loadReceiversAndTransmitters(): To implement Full_Interstation_TF !!!!"
                !
            case( "Off_Diagonal_Rho_Phase" )
                !
                stop "DataManager.f08: loadReceiversAndTransmitters(): To implement Off_Diagonal_Rho_Phase !!!!"
                !
            case( "Phase_Tensor" )
                !
                stop "DataManager.f08: loadReceiversAndTransmitters(): To implement Phase_Tensor !!!!"
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
                write( *, * ) "unknow component type :[", data_entry%type, "]"
                stop "DataManager.f08: loadReceiversAndTransmitters()"
            !
        end select
        !
        receiver%is_complex = data_entry%isComplex()
        !
        receiver%code = data_entry%code
        !
        rx_id = updateReceiverArray( receiver )
        !
        nTx = size( transmitters )
        !
        ! LOOP OVER TRANSMITTERS
        do iTx = 1, nTx
            !
            transmitter => getTransmitter( iTx )
            !
            if( ABS( transmitter%period - data_entry%period ) < TOL6 ) then
                !
                call transmitter%updateReceiverIndexesArray( rx_id )
                !
                exit
                !
            endif
            !
        enddo
        !
    end subroutine loadReceiversAndTransmitters
    !
end module DataFile
