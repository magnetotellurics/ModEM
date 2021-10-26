!*************
!
! Class to create an array of transmitters related with their receivers,
! from an instance of DataFile_t object readed from a data file
!
! Last modified at 08/06/2021 by Paulo Werdt
!
!*************
!
module DataManager
   !
   use FileUnits
   !
   use StandardDataFile
   use DataGroupArray
   use ReceiverFullImpedance
   use ReceiverFullVerticalMagnetic
   use ReceiverOffDiagonalImpedance
   use ReceiverSingleField
   use ReceiverArray
   use TransmitterMT
   use TransmitterCSEM
   use TransmitterArray
   !
   ! Global Array of Transmitters
   class( TransmitterArray_t ), pointer, save, public :: transmitters
   !
   ! Global Array of Receivers
   class( ReceiverArray_t ), pointer, save, public :: receivers
   !
   type :: DataManager_t 
      !
      class( DataFile_t ), pointer :: data_file
      !
   contains
      !
      final :: DataManager_dtor
      !
      procedure, private :: loadReceiversAndTransmitters, relateData
      !
   end type DataManager_t
   !
   interface DataManager_t
      module procedure DataManager_ctor
   end interface DataManager_t
   !
contains
   !
   ! Create Data File object, create Receivers, create transmitters and create
   function DataManager_ctor( file_name) result( self )
      !
      character(:), allocatable, intent( in ) :: file_name
	  type( DataManager_t ), pointer          :: self
      !
      ! write(*,*) "Constructor DataManager_t"
      !
      allocate( DataManager_t :: self )
      !
      self%data_file => StandardDataFile_t( ioStartup, file_name )
      !
      call self%loadReceiversAndTransmitters()
      !
      call self%relateData()
      !
      write( *, * ) receivers%size(), " Receivers"
      !
      write( *, * ) transmitters%size(), " Transmitters"
      !
   end function DataManager_ctor 
   !
   subroutine DataManager_dtor( self )
      implicit none
      !
      type( DataManager_t ), intent( in out ) :: self
      !
      !write(*,*) "Destructor DataManager_t"
      !
      deallocate( self%data_file )
      !
   end subroutine DataManager_dtor
   !
   ! Load all Receivers (Based on Location and Component type) and all Transmitters (Based on Period)
   subroutine loadReceiversAndTransmitters( self )
      implicit none
      !
      class( DataManager_t ), intent( inout ) :: self
      !
      class( DataEntry_t ), pointer   :: data_entry
      class( Receiver_t ), pointer    :: receiver
      class( Transmitter_t ), pointer :: transmitter
      integer                         :: iDe, nDe, iRx, iTx
      real ( kind=prec )              :: azimuth
      !
      receivers => ReceiverArray_t()
      !
      transmitters => TransmitterArray_t()
      !
      ! Loop over all Data Entries...
      nDe = self%data_file%data_entries%size()
      !
      do iDe = 1, nDe
         !
         data_entry => self%data_file%get( iDe )
         !
         ! RECEIVERS
         !
         iRx = receivers%size() + 1
         !
         selectcase( data_entry%type )
            !
            case( "Ex_Field" )
               !
               azimuth = 1.0
               receiver => ReceiverSingleField_t( iRx, data_entry%xyz, azimuth )
               !
            case( "Ey_Field" )
               !
               azimuth = 2.0
               receiver => ReceiverSingleField_t( iRx, data_entry%xyz, azimuth )
               !
            case( "Bx_Field" )
               !
               azimuth = 3.0
               receiver => ReceiverSingleField_t( iRx, data_entry%xyz, azimuth )
               !
            case( "By_Field" )
               !
               azimuth = 4.0
               receiver => ReceiverSingleField_t( iRx, data_entry%xyz, azimuth )
               !
            case( "Bz_Field" )
               !
               azimuth = 5.0
               receiver => ReceiverSingleField_t( iRx, data_entry%xyz, azimuth )
               !
            case( "Full_Impedance" )
               !
               receiver => ReceiverFullImpedance_t( iRx, data_entry%xyz )
               !
            case( "Full_Interstation_TF" )
               !
               STOP "DataManager.f08: loadReceiversAndTransmitters(): To implement Full_Interstation_TF !!!!"
               !
            case( "Off_Diagonal_Rho_Phase" )
               !
               STOP "DataManager.f08: loadReceiversAndTransmitters(): To implement Off_Diagonal_Rho_Phase !!!!"
               !
            case( "Phase_Tensor" )
               !
               STOP "DataManager.f08: loadReceiversAndTransmitters(): To implement Phase_Tensor !!!!"
               !
            case( "Off_Diagonal_Impedance" )
               !
               receiver => ReceiverOffDiagonalImpedance_t( iRx, data_entry%xyz )
               !
            case( "Full_Vertical_Components", "Full_Vertical_Magnetic" )
               !
               receiver => ReceiverFullVerticalMagnetic_t( iRx, data_entry%xyz )
               !
            case default
               write(*,*) "unknow component type :[", data_entry%type, "]"
               STOP "DataManager.f08: loadReceiversAndTransmitters()"
            !
         end select
         !
         receiver%is_complex = data_entry%isComplex()
         !
		 receiver%code = data_entry%code
		 !
         if( .NOT. receivers%has( receiver ) ) then 
            call receivers%add( receiver )
         end if
         !
         ! TRANSMITTERS
         !
         iTx = transmitters%size() + 1
         !
         selecttype ( data_entry )
            !
            class is ( DataEntryMT_t )
               !
               transmitter => Transmitter_t( iTx, data_entry%period )
               !
            class is ( DataEntryMT_REF_t )
               !
               transmitter => Transmitter_t( iTx, data_entry%period )
               !
            class is ( DataEntryCSEM_t )
               !
               transmitter => Transmitter_t( iTx, data_entry%period, data_entry%tx_xyz )
               !
         end select
         !
         if( .NOT. transmitters%has( transmitter ) ) then 
            call transmitters%add( transmitter )
         end if
         !
      enddo
      !
   end subroutine loadReceiversAndTransmitters
   !
   ! create DataGroups and relate with Receivers and Transmitters
   subroutine relateData( self )
      implicit none
      !
      class( DataManager_t ), intent( in out ) :: self
	  !
      class( DataEntry_t ), pointer   :: data_entry
      class( Transmitter_t ), pointer :: transmitter
      class( Receiver_t ), pointer    :: receiver
      class( DataGroup_t ), pointer   :: data_group
      integer                         :: iDg, iDe, nDe, iRx, nRx, iTx, nTx, counter
      !
      iDg = 1
      !
      ! Loop over all Data Entries... 
      nDe = self%data_file%data_entries%size()
      !
      iDe = 1
      do while ( iDe <= nDe ) 
         !
         data_entry => self%data_file%get( iDe )
         !
         selectcase( data_entry%type )
            !
            case( "Ex_Field", "Ey_Field", "Bx_Field", "By_Field","Bz_Field" )
               !
               data_group => DataGroup_t( iDg, 1 )
               !
            case( "Full_Impedance", "Full_Interstation_TF", "Off_Diagonal_Rho_Phase", "Phase_Tensor" )
               !
               data_group => DataGroup_t( iDg, 4 )
               !
            case( "Off_Diagonal_Impedance", "Full_Vertical_Components" )
               !
               data_group => DataGroup_t( iDg, 2 )
               !
            case default
               write(*,*) "unknow type :[", data_entry%type, "]"
               STOP "DataManager.f08: relateData()"
            !
         end select
         !
         call data_group%add( data_entry%component, data_entry%real, data_entry%imaginary, data_entry%error )
         !
         do counter = 1, data_group%n_data - 1
            data_entry => self%data_file%get( iDe + counter )
            call data_group%add( data_entry%component, data_entry%real, data_entry%imaginary, data_entry%error )
         end do
         !
         ! LOOP OVER RECEIVERS
         nRx = receivers%size()
         do iRx = 1, nRx
            !
            receiver => receivers%get( iRx )
            !
            if( receiver%location(1) == data_entry%xyz(1) .AND.   &
               receiver%location(2) == data_entry%xyz(2) .AND.   &
               receiver%location(3) == data_entry%xyz(3) ) then
               !
               data_group%id_rx = receiver%id
               !
               if( .NOT. receiver%has( data_group ) ) then
                  call receiver%add( data_group )
               end if
               !
               exit
               !
            endif
            !
         enddo
         !
         ! LOOP OVER TRANSMITTERS
         nTx = transmitters%size()
         do iTx = 1, nTx
            !
            transmitter => transmitters%get( iTx )
            !
            if( ABS( transmitter%period - data_entry%period ) < TOL6 ) then
               !
               data_group%id_tx = transmitter%id
               !
               if( .NOT. transmitter%has( receiver%id ) ) then
                  call transmitter%add( receiver%id )
               end if
               !
               exit
               !
            endif
            !
         enddo
         !
         iDe = iDe + data_group%n_data
         !
         if( iDe < nDe ) iDg = iDg + 1
         !
      enddo
      !
   end subroutine relateData
   !
end module DataManager
