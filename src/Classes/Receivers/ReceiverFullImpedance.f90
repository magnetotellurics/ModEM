! *************
! 
! Derived class to define a Full_Impedance Receiver
! 
! Last modified at 08/06/2021 by Paulo Werdt
! 
! *************
! 
module ReceiverFullImpedance
   !
   use FileUnits
   use DataEntryMT
   use Receiver
   !
   type, extends( Receiver_t ), public :: ReceiverFullImpedance_t
      !
      ! SPECIFIC PROPERTIES HERE
      !
      contains
         !
         final :: ReceiverFullImpedance_dtor
         !
         procedure, public :: predictedData => predictedDataFullImpedance
         !
         procedure, public :: savePredictedData => savePredictedDataFullImpedance
         procedure, public :: writePredictedData => writePredictedDataFullImpedance
         procedure, public :: write => writeReceiverFullImpedance
         !
   end type ReceiverFullImpedance_t
   !
   interface ReceiverFullImpedance_t
      module procedure ReceiverFullImpedance_ctor
   end interface ReceiverFullImpedance_t
   !
contains
   !
   function ReceiverFullImpedance_ctor( id, location ) result( self )
      implicit none
      !
      integer, intent( in )           :: id
      real( kind=prec ), intent( in ) :: location(3)
      type( ReceiverFullImpedance_t ) :: self
      !
      !write(*,*) "Constructor ReceiverFullImpedance_t"
      !
      call self%init()
      !
      self%id = id
      self%location = location
      self%n_comp = 4
      self%is_complex = .TRUE.
      !
      !
      allocate( character(2) :: self%EHxy( 4 ) )
      !
      self%EHxy(1)="Ex"
      self%EHxy(2)="Ey"
      self%EHxy(3)="Bx"
      self%EHxy(4)="By"
      !
      ! components required to get the full impdedance tensor Z [Zxx, Zxy, Zyx, Zyy]
      allocate( character(3) :: self%comp_names( 4 ) )
      !
     self%comp_names(1) = "ZXX"
     self%comp_names(2) = "ZXY"
     self%comp_names(3) = "ZYX"
     self%comp_names(4) = "ZYY"
     !
   end function ReceiverFullImpedance_ctor
   !
   subroutine ReceiverFullImpedance_dtor( self )
      implicit none
      !
      type( ReceiverFullImpedance_t ), intent( in out ) :: self
      !
      !write(*,*) "Destructor ReceiverFullImpedance_t:", self%id
      !
      call self%dealloc()
      !
   end subroutine ReceiverFullImpedance_dtor
   !
   subroutine predictedDataFullImpedance( self, model_operator, transmitter )
      implicit none
      !
      class( ReceiverFullImpedance_t ), intent( inout ) :: self
      class( ModelOperator_t ), intent( in )            :: model_operator
      class( Transmitter_t ), intent( in )              :: transmitter
      !
      class( cVector_t ), allocatable   :: e_tx_pol_1, e_tx_pol_2
      complex( kind=prec ), allocatable :: BB(:,:), det
      real( kind=prec )                 :: omega
      integer                           :: i, j, ij
      !
      omega = ( 2.0 * PI / transmitter%period )
      !
      ! Set Vectors Lex, Ley, Lbx, Lby
      call self%evaluationFunction( model_operator, omega )
      !
      ! get e_all from the Tx 1st polarization
      allocate( e_tx_pol_1, source = transmitter%e_all( 1 ) )
      !
      ! get e_all from the Tx 2nd polarization
      allocate( e_tx_pol_2, source = transmitter%e_all( 2 ) )
      !
      allocate( complex( kind=prec ) :: self%EE( 2, 2 ) )
      !
      self%EE( 1, 1 ) = self%Lex .dot. e_tx_pol_1
      self%EE( 2, 1 ) = self%Ley .dot. e_tx_pol_1
      self%EE( 1, 2 ) = self%Lex .dot. e_tx_pol_2
      self%EE( 2, 2 ) = self%Ley .dot. e_tx_pol_2
      !
      !write(*,*) "EE:"
      !write(*,*) self%EE( 1, 1 ), self%EE( 1, 2 )
      !write(*,*) self%EE( 2, 1 ), self%EE( 2, 2 )
      !
      allocate( complex( kind=prec ) :: BB( 2, 2 ) )
      !
      BB( 1, 1 ) = self%Lbx .dot. e_tx_pol_1
      BB( 2, 1 ) = self%Lby .dot. e_tx_pol_1
      BB( 1, 2 ) = self%Lbx .dot. e_tx_pol_2
      BB( 2, 2 ) = self%Lby .dot. e_tx_pol_2
      !
      deallocate( e_tx_pol_1 )
      deallocate( e_tx_pol_2 )
      !
      !write(*,*) "BB:"
      !write(*,*) BB( 1, 1 ), BB( 1, 2 )
      !write(*,*) BB( 2, 1 ), BB( 2, 2 )
      !
      !invert horizontal B matrix using Kramer"s rule.
      det = BB( 1, 1 ) * BB( 2, 2 ) - BB( 1, 2 ) * BB( 2, 1 )
      !
      !write(*,*) "det:", det
      !
      allocate( complex( kind=prec ) :: self%I_BB( 2, 2 ) )
      !
      if( det /= 0 ) then
         self%I_BB( 1, 1 ) = BB( 2, 2 ) / det
         self%I_BB( 2, 2 ) = BB( 1, 1 ) / det
         self%I_BB( 1, 2 ) = -BB( 1, 2 ) / det
         self%I_BB( 2, 1 ) = -BB( 2, 1 ) / det
      else
         STOP "ReceiverFullImpedance.f90: Determinant is Zero!"
      endif
      !
      !write(*,*) "Inverse BB:"
      !write(*,*) self%I_BB( 1, 1 ), self%I_BB( 1, 2 )
      !write(*,*) self%I_BB( 2, 1 ), self%I_BB( 2, 2 )
      !
      if( .not. allocated( self%Z ) ) allocate( complex( kind=prec ) :: self%Z( 4 ) )
      !
      do j = 1,2
         do i = 1,2
            ij = 2 * ( i-1 ) + j
            self%Z( ij ) = self%EE( i, 1 ) * self%I_BB( 1, j ) + self%EE( i, 2 ) * self%I_BB( 2, j )
         enddo
      enddo
      !
      ! WRITE ON PredictedFile.dat
      call self%savePredictedData( transmitter )
      !
      deallocate( self%EE )
      deallocate( BB )
      deallocate( self%I_BB )
      deallocate( self%Z )
      !
   end subroutine predictedDataFullImpedance
   !
   !
   subroutine savePredictedDataFullImpedance( self, tx )
      implicit none
      !
      class( ReceiverFullImpedance_t ), intent( in ) :: self
      class( Transmitter_t ), intent( in )           :: tx
      !
      character(:), allocatable         :: type, code, component
      integer                           :: i, iDe
      real( kind=prec )                 :: period, real_part, imaginary, error
      real( kind=prec )                 :: latitude, longitude, xyz(3)
      !#Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
      !
      do i = 1, 4
          !
          type = "Output"
          iDe = self%predicted_data_entries%size() + i
          period = real( tx%period, kind=prec )
          code = trim( self%code )
          latitude = R_ZERO
          longitude = R_ZERO
          xyz = (/real( self%location( 1 ), kind=prec ), real( self%location( 2 ), kind=prec ), real( self%location( 3 ), kind=prec )/)
          component = trim( self%comp_names( i ) )
          real_part = real( self%Z( i ), kind=prec )
          imaginary = real( aimag( self%Z( i ) ), kind=prec )
          error = 1.0
          !
          call self%predicted_data_entries%add( DataEntryMT_t( iDe, type, period, code, latitude, longitude, xyz, component, real_part, imaginary, error ) )
          !
      enddo
	  !
   end subroutine savePredictedDataFullImpedance
   !
   subroutine writePredictedDataFullImpedance( self )
      implicit none
      !
      class( ReceiverFullImpedance_t ), intent( in ) :: self
      !
      class( DataEntry_t ), allocatable :: data_entry
      integer :: iDe, Nde
      !
      open( ioPredData, file = "predicted_data.dat", action="write", position="append" )
      !
      nDe = self%predicted_data_entries%size()
      !
      do iDe = 1, nDe
          !
          data_entry = self%predicted_data_entries%get( iDe )
          !
          !#Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
          write( ioPredData, "(1pe12.6, A8, f9.3, f9.3, f13.3, f13.3, f13.3, A4, 1pe16.6, 1pe16.6, 1pe16.6)" ) data_entry%period, data_entry%code, R_ZERO, R_ZERO, data_entry%xyz(1) / 10.0, data_entry%xyz(2) / 10.0, data_entry%xyz(3) / 10.0, data_entry%component, data_entry%real, data_entry%imaginary, 1.0
          !
      enddo
      !
      close( ioPredData )
      !
   end subroutine writePredictedDataFullImpedance
   !
   subroutine writeReceiverFullImpedance( self )
      implicit none
      !
      class( ReceiverFullImpedance_t ), intent( in ) :: self
      !
      integer                           :: iDg, nDg
      class( DataGroup_t ), allocatable :: data_group
      !
      nDg = self%getNDg()
      !
      write(*,*) "Write ReceiverFullImpedance_t: ", self%id,   &
      " N Data Groups: ", nDg
      !
      do iDg = 1, nDg
         data_group = self%get( iDg )
         call data_group%write()
      enddo
      !
   end subroutine writeReceiverFullImpedance
   !
end module ReceiverFullImpedance
