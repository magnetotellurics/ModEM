! *************
! 
! Derived class to define a MT Transmitter
! 
! Last modified at 10/11/2021 by Paulo Werdt
! 
! *************
! 
module TransmitterMT
   ! 
   use Transmitter
   !
   type, extends( Transmitter_t ), public :: TransmitterMT_t
      !
      ! PROPERTIES HERE
      !
      contains
         !
         final :: TransmitterMT_dtor
         !
         procedure, public :: solveFWD  => solveFWDTransmitterMT
         !
         procedure, public :: getType => getTypeTransmitterMT
         procedure, public :: isEqual => isEqualTransmitterMT
         procedure, public :: write   => writeTransmitterMT
         !
   end type TransmitterMT_t
   !
   interface TransmitterMT_t
      module procedure TransmitterMT_ctor
   end interface TransmitterMT_t
   !
   contains
   !
   ! TransmitterMT constructor
   function TransmitterMT_ctor( id, period ) result ( self )
      implicit none
      !
      type( TransmitterMT_t ) :: self
      !
      integer, intent( in )           :: id
      real( kind=prec ), intent( in ) :: period
      !
      !write(*,*) "Constructor TransmitterMT_t"
      !
      call self%init()
      !
      self%id = id
      self%n_pol = 2
      self%period = period
      !
      self%DATA_TITLE = "Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
      !
   end function TransmitterMT_ctor
   !
   ! TransmitterMT destructor
   subroutine TransmitterMT_dtor( self )
      implicit none
      !
      type( TransmitterMT_t ), intent( inout ) :: self
      !
      !write(*,*) "Destructor TransmitterMT_t:", self%id
      !
      call self%dealloc()
      !
   end subroutine TransmitterMT_dtor
   !
   ! Set self%e_all from forward modelling solver
   subroutine solveFWDTransmitterMT( self )
      implicit none
      !
      class( TransmitterMT_t ), intent( inout ) :: self
      !
      integer           :: i_pol
      real( kind=prec ) :: omega
      !
      ! verbosis
      write( *, * ) "   Solving FWD for Tx", self%id
	  !
      omega = 2.0 * PI / self%period
      !
      ! Set ForwardSolver Frequency
      call self%forward_solver%setPeriod( self%period )
      !
      ! Loop over all polarizations (MT n_pol = 2)
      do i_pol = 1, self%n_pol
      !
         write(*,*) "MT Tx Solve for Polarization", i_pol
         !
         ! Set Source E
         call self%source%setE( omega, i_pol )
		 !
		 write(*,*) "2"
         ! Add polarization e_solution to self%e_all
         call self%e_all%add( self%forward_solver%getESolution( self%source, i_pol ) )
		 write(*,*) "3"
      !
      enddo
      !
   end subroutine solveFWDTransmitterMT
   !
   ! Get class string name
   function getTypeTransmitterMT( self ) result( type )
      implicit none
      !
      class( TransmitterMT_t ), intent( in ) :: self
      character(:), allocatable              :: type
      !
      type = "TransmitterMT_t"
      !
   end function getTypeTransmitterMT
   !
   ! Compare two transmitters
   function isEqualTransmitterMT( self, other ) result( equal )
      implicit none
      !
      class( TransmitterMT_t ), intent( in ) :: self
      class( Transmitter_t ), intent( in )   :: other
      logical                                :: equal
      !
      equal = .FALSE.
      !
      select type( other )
         !
         class is ( TransmitterMT_t )
            !
            if( ABS( self%period - other%period ) < TOL6 ) then
               equal = .TRUE.
            endif
            !
      end select
      !
   end function isEqualTransmitterMT
   !
   ! Print TransmitterMT info
   subroutine writeTransmitterMT( self )
      implicit none
      !
      class( TransmitterMT_t ), intent( in ) :: self
      !
      integer :: iRx, nRx
      !
      nRx = self%getNRx()
      !
      write(*,*) "Write TransmitterMT_t: ", self%id,   &
      " Period: ",   self%period,   &
      "fwd_key: ",   self%fwd_key(1), self%fwd_key(2), self%fwd_key(3), self%fwd_key(4),   &
                  self%fwd_key(5), self%fwd_key(6), self%fwd_key(7), self%fwd_key(8),   &
      " N Receivers: ", nRx
      !
      do iRx = 1, nRx
         write(*,*) "   ", self%get( iRx )
      enddo
      !
   end subroutine writeTransmitterMT
   !
end module TransmitterMT
