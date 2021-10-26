! *************
! 
! Derived class to define a Single E or B Field Receiver
! 
! Last modified at 30/06/2021 by Paulo Werdt
! 
! *************
! 
module ReceiverSingleField
   ! 
   use Receiver
   !
   type, extends( Receiver_t ), public :: ReceiverSingleField_t
      !
      real ( kind=prec )                 :: azimuth
      !
      contains
         !
         final :: ReceiverSingleField_dtor
         !
         procedure, public :: predictedData => predictedDataSingleField
         !
		 procedure, public :: writePredictedData => writePredictedDataSingleField
         procedure, public :: write => writeReceiverSingleField
         !
   end type ReceiverSingleField_t
   !
   interface ReceiverSingleField_t
      module procedure ReceiverSingleField_ctor
   end interface ReceiverSingleField_t
   !
contains
   !
   function ReceiverSingleField_ctor( id, location, azimuth ) result( self )
      !
      class( ReceiverSingleField_t ), pointer :: self
      integer, intent( in )                   :: id
      real( kind=prec ), intent( in )         :: location(3)
      real ( kind=prec ), intent( in )        :: azimuth
      !
      ! write(*,*) "Constructor ReceiverSingleField_t"
      !
      allocate( ReceiverSingleField_t :: self )
      !
      call self%init()
      !
      self%id = id
      self%location = location
      self%azimuth = azimuth
      !
      self%n_comp = 1
      self%is_complex = .TRUE.
      !
      !
      allocate( character(2) :: self%EHxy( 1 ) )
      !
      ! components required to get the full impdence tensor Z [Zxx, Zxy, Zyx, Zyy]
      !
      if( azimuth == 1.0 ) self%EHxy(1)="Ex"
      if( azimuth == 2.0 ) self%EHxy(1)="Ey"
      if( azimuth == 3.0 ) self%EHxy(1)="Bx"
      if( azimuth == 4.0 ) self%EHxy(1)="By"
      if( azimuth == 5.0 ) self%EHxy(1)="Bz"
      !
   end function ReceiverSingleField_ctor
   !
   subroutine ReceiverSingleField_dtor( self )
      implicit none
      !
      type( ReceiverSingleField_t ), intent( in out ) :: self
      !
      ! write(*,*) "Destructor ReceiverSingleField_t"
      !
      call self%dealloc()
      !
   end subroutine ReceiverSingleField_dtor
   !
   subroutine predictedDataSingleField( self, model_operator, transmitter )
      implicit none
	  class( ReceiverSingleField_t ), intent( inout )     :: self
      class( ModelOperator_t ), allocatable, intent( in ) :: model_operator
      class( Transmitter_t ), intent( in )                :: transmitter
      !
      complex(kind=prec) :: det, ctemp
      !
      write(*,*) "Implement predictedData ReceiverSingleField_t: ", self%id
      !
      allocate( complex(kind=prec) :: self%Z( 1 ) )
      !
      ! FIND BEST WAY TO IMPLEMENT
      !Z(1) = dotProd_noConj_scvector_f( Lex,ef%pol(1) )
      !
   end subroutine predictedDataSingleField
   !
   !
   subroutine writePredictedDataSingleField( self, tx )
      implicit none
      !
	  class( ReceiverSingleField_t ), intent( in ) :: self
      class( Transmitter_t ), intent( in )         :: tx
      !
      open( 666, file = 'predicted_data.dat', action='write', position='append' )
      !
      write( 666, '(1pe12.6, A8, f9.3, f9.3, f13.3, f13.3, f13.3, A4, 1pe16.6, 1pe16.6, 1pe16.6)' ) tx%period, self%code, R_ZERO, R_ZERO, self%location(1), self%location(2), self%location(3), self%comp_names( 1 ), aimag( self%Z( 1 ) ), dimag( self%Z( 1 )), 1.0
      write( 666, '(1pe12.6, A8, f9.3, f9.3, f13.3, f13.3, f13.3, A4, 1pe16.6, 1pe16.6, 1pe16.6)' ) tx%period, self%code, R_ZERO, R_ZERO, self%location(1), self%location(2), self%location(3), self%comp_names( 2 ), aimag( self%Z( 2 ) ), dimag( self%Z( 2 )), 1.0
      write( 666, '(1pe12.6, A8, f9.3, f9.3, f13.3, f13.3, f13.3, A4, 1pe16.6, 1pe16.6, 1pe16.6)' ) tx%period, self%code, R_ZERO, R_ZERO, self%location(1), self%location(2), self%location(3), self%comp_names( 3 ), aimag( self%Z( 3 ) ), dimag( self%Z( 3 )), 1.0
      write( 666, '(1pe12.6, A8, f9.3, f9.3, f13.3, f13.3, f13.3, A4, 1pe16.6, 1pe16.6, 1pe16.6)' ) tx%period, self%code, R_ZERO, R_ZERO, self%location(1), self%location(2), self%location(3), self%comp_names( 4 ), aimag( self%Z( 4 ) ), dimag( self%Z( 4 )), 1.0
      !
      close( 666 )
      !
   end subroutine writePredictedDataSingleField
   !
   subroutine writeReceiverSingleField( self )
      class( ReceiverSingleField_t ), intent( in ) :: self
      !
      integer                     :: iDg, nDg
      class( DataGroup_t ), pointer   :: data_group
      !
      nDg = self%getNDg()
      !
      write(*,*) "Write ReceiverSingleField_t: ", self%id,   &
      " N Data Groups: ", nDg
      !
      do iDg = 1, nDg
         data_group => self%get( iDg )
         call data_group%write()
      enddo
      !
   end subroutine writeReceiverSingleField
   !
end module ReceiverSingleField
