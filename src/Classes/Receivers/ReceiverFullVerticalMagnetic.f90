! *************
! 
! Derived class to define a Full_Impedance Receiver
! 
! Last modified at 08/06/2021 by Paulo Werdt
! 
! *************
! 
module ReceiverFullVerticalMagnetic
   !
   use Receiver
   !
   type, extends( Receiver_t ), public :: ReceiverFullVerticalMagnetic_t
      !
      ! PROPERTIES HERE
      !
      contains
         !
         final :: ReceiverFullVerticalMagnetic_dtor
         !
         procedure, public :: predictedData => predictedDataFullVerticalMagnetic
         !
         procedure, public :: write => writeReceiverFullVerticalMagnetic
         !
   end type ReceiverFullVerticalMagnetic_t
   !
   interface ReceiverFullVerticalMagnetic_t
      module procedure ReceiverFullVerticalMagnetic_ctor
   end interface ReceiverFullVerticalMagnetic_t
   !
contains
   !
   function ReceiverFullVerticalMagnetic_ctor( location ) result( self )
      implicit none
      !
      real( kind=prec ), intent( in )        :: location(3)
      type( ReceiverFullVerticalMagnetic_t ) :: self
      !
      ! write(*,*) "Constructor ReceiverFullVerticalMagnetic_t"
      !
      call self%init()
      !
      self%location = location
      !
      self%n_comp = 2
      self%is_complex = .TRUE.
      !
      allocate( character(2) :: self%EHxy( 3 ) )
      !
      ! components required to get the full impdence tensor Z [Zxx, Zxy, Zyx, Zyy]
      self%EHxy(1)="Bx"
      self%EHxy(2)="By"
      self%EHxy(3)="Bz"
      !
   end function ReceiverFullVerticalMagnetic_ctor
   !
   subroutine ReceiverFullVerticalMagnetic_dtor( self )
      implicit none
      !
      type( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
      !
      ! write(*,*) "Destructor ReceiverFullVerticalMagnetic_t"
      !
      call self%dealloc()
      !
   end subroutine ReceiverFullVerticalMagnetic_dtor
   !
   subroutine writeReceiverFullVerticalMagnetic( self )
      implicit none
      !
      class( ReceiverFullVerticalMagnetic_t ), intent( in ) :: self
      !
      integer                             :: iDg, nDg
      class( DataGroup_t ), allocatable   :: data_group
      !
      nDg = self%getNDg()
      !
      write(*,*) "Write ReceiverFullVerticalMagnetic_t: ", self%id,   &
      " N Data Groups: ", nDg
      !
      do iDg = 1, nDg
         data_group = self%get( iDg )
         call data_group%write()
      enddo
      !
   end subroutine writeReceiverFullVerticalMagnetic
   !
   subroutine predictedDataFullVerticalMagnetic( self, model_operator, transmitter )
      implicit none
      !
      class( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
      class( ModelOperator_t ), intent( in )                   :: model_operator
      class( Transmitter_t ), intent( in )                     :: transmitter
      !
      class( cVector_t ), allocatable :: e_tx_pol_1, e_tx_pol_2
      complex( kind=prec ), allocatable :: BB(:,:), det
      real( kind=prec )               :: omega
      integer                         :: i, j, ij
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
      allocate( complex( kind=prec ) :: BB( 3, 2 ) )
      !
      BB(1,1)= self%Lbx .dot. e_tx_pol_1
      BB(2,1)= self%Lby .dot. e_tx_pol_1
      BB(1,2)= self%Lbx .dot. e_tx_pol_2
      BB(2,2)= self%Lby .dot. e_tx_pol_2
      BB(3,1)= self%Lbz .dot. e_tx_pol_1
      BB(3,2)= self%Lbz .dot. e_tx_pol_2
      !
      deallocate( e_tx_pol_1 )
      deallocate( e_tx_pol_2 )
      !
      !invert horizontal B matrix using Kramer's rule.
      det = BB(1,1) * BB(2,2) - BB(1,2) * BB(2,1)
      !
      allocate( complex( kind=prec ) :: self%I_BB( 2, 2 ) )
      !
      if( det /= 0 ) then
         self%I_BB( 1, 1 ) = BB( 2, 2 ) / det
         self%I_BB( 2, 2 ) = BB( 1, 1 ) / det
         self%I_BB( 1, 2 ) = -BB( 1, 2 ) / det
         self%I_BB( 2, 1 ) = -BB( 2, 1 ) / det
      else
         STOP "ReceiverFullVerticalMagnetic.f90: Determinant is Zero!"
      endif
      !
      allocate( complex( kind=prec ) :: self%Z( 2 ) )
      !
      self%Z(1) = self%I_BB(3,1) * self%I_BB(1,1) + self%I_BB(3,2) * self%I_BB(2,1)
      self%Z(2) = self%I_BB(3,1) * self%I_BB(1,2) + self%I_BB(3,2) * self%I_BB(2,2)
      !
      ! WRITE ON PredictedFile.dat
      call self%savePredictedData( transmitter )
      !
      deallocate( BB )
      deallocate( self%I_BB )
      deallocate( self%Z )
      !
   end subroutine predictedDataFullVerticalMagnetic
   !
end module ReceiverFullVerticalMagnetic
