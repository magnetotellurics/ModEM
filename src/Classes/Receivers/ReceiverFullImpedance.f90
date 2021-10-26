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
      !
      class( ReceiverFullImpedance_t ), pointer :: self
      integer, intent( in )                     :: id
      real( kind=prec ), intent( in )           :: location(3)
      !
      allocate( ReceiverFullImpedance_t :: self )
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
      !write(*,*) "Destructor ReceiverFullImpedance_t"
      !
      call self%dealloc()
      !
   end subroutine ReceiverFullImpedance_dtor
   !
   subroutine predictedDataFullImpedance( self, model_operator, transmitter )
      !
      class( ReceiverFullImpedance_t ), intent( inout )   :: self
      class( ModelOperator_t ), allocatable, intent( in ) :: model_operator
      class( Transmitter_t ), pointer, intent( in )       :: transmitter
      !
      class( cVector_t ), allocatable :: e_tx_pol_1, e_tx_pol_2
      complex(kind=prec), allocatable :: BB(:,:), det
      real( kind=prec )               :: omega
      integer                         :: i, j, ij
      !
      omega = ( 2.0 * PI / transmitter%period )
      !
      ! Set Vectors Lex, Ley, Lbx, Lby
      call self%evaluationFunction( model_operator, omega )
      !
      ! get e_solution from the Tx 1st polarization
      !e_tx_pol_1 = transmitter%e_solution%get( 1 ) ! SHOULD BE LIKE THIS???
      allocate( e_tx_pol_1, source = transmitter%e_solution%get( 1 ) )
      !
      ! get e_solution from the Tx 2nd polarization
      !e_tx_pol_2 = transmitter%e_solution%get( 2 )
      allocate( e_tx_pol_2, source = transmitter%e_solution%get( 2 ) )
      !
      allocate( complex(kind=prec) :: self%EE( 2, 2 ) )
      !
      self%EE( 1, 1 ) = self%Lex .dot. e_tx_pol_1
      self%EE( 2, 1 ) = self%Ley .dot. e_tx_pol_1
      self%EE( 1, 2 ) = self%Lex .dot. e_tx_pol_2
      self%EE( 2, 2 ) = self%Ley .dot. e_tx_pol_2
      !
      allocate( complex(kind=prec) :: BB( 2, 2 ) )
      !
      BB( 1, 1 ) = self%Lbx .dot. e_tx_pol_1
      BB( 2, 1 ) = self%Lby .dot. e_tx_pol_1
      BB( 1, 2 ) = self%Lbx .dot. e_tx_pol_2
      BB( 2, 2 ) = self%Lby .dot. e_tx_pol_2
      !
     !write(*,*) "BB:"
     !write(*,*) BB( 1, 1 ), BB( 1, 2 )
     !write(*,*) BB( 2, 1 ), BB( 2, 2 )
      !
      !invert horizontal B matrix using Kramer's rule.
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
      if( .not. allocated( self%Z ) ) allocate( complex(kind=prec) :: self%Z( 4 ) )
      !
      do j = 1,2
         do i = 1,2
            ij = 2 * ( i-1 ) + j
            self%Z( ij ) = self%EE( i, 1 ) * self%I_BB( 1, j ) + self%EE( i, 2 ) * self%I_BB( 2, j )
         enddo
      enddo
      !
	  ! WRITE ON PredictedFile.dat
	  call self%writePredictedData( transmitter )
	  !
      deallocate( self%EE )
      deallocate( BB )
      deallocate( self%I_BB )
      deallocate( self%Z )
      !
   end subroutine predictedDataFullImpedance
   !
   !
   subroutine writePredictedDataFullImpedance( self, tx )
      implicit none
      !
	  class( ReceiverFullImpedance_t ), intent( in ) :: self
      class( Transmitter_t ), intent( in )           :: tx
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
   end subroutine writePredictedDataFullImpedance
   !
   subroutine writeReceiverFullImpedance( self )
      class( ReceiverFullImpedance_t ), intent( in ) :: self
      !
      integer                       :: iDg, nDg
      class( DataGroup_t ), pointer :: data_group
      !
      nDg = self%getNDg()
      !
      write(*,*) "Write ReceiverFullImpedance_t: ", self%id,   &
      " N Data Groups: ", nDg
      !
      do iDg = 1, nDg
         data_group => self%get( iDg )
         call data_group%write()
      enddo
      !
   end subroutine writeReceiverFullImpedance
   !
end module ReceiverFullImpedance
