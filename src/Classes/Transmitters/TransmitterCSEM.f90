! *************
! 
! Derived class to define a CSEM Transmitter
! 
! Last modified at 08/06/2021 by Paulo Werdt
! 
! *************
! 
module TransmitterCSEM
   ! 
   use Transmitter 
   !
   type, extends( Transmitter_t ), public :: TransmitterCSEM_t
      !
      real( kind=prec ) :: location(3), azimuth
      !
      contains
         !
         final   :: TransmitterCSEM_dtor
         !
         procedure, public   :: solveFWD => solveFWDTransmitterCSEM
         procedure, public   :: getSource => getSourceTransmitterCSEM
         !
		 !procedure, public   :: sizeOf => sizeOfTransmitterCSEM
         procedure, public   :: isEqual => isEqualTransmitterCSEM
         procedure, public   :: write => writeTransmitterCSEM
         !
   end type TransmitterCSEM_t
   !
   interface TransmitterCSEM_t
      module procedure TransmitterCSEM_ctor
   end interface TransmitterCSEM_t
   !
contains
   !
   ! Parametrized constructor
   function TransmitterCSEM_ctor( id, period, location, type ) result ( self )
      !
      type( TransmitterCSEM_t ) :: self
      !
      integer, intent( in )                             :: id
      real( kind=prec ), intent( in )                   :: period
      real( kind=prec ), intent( in )                   :: location(3)
      character(:), allocatable, optional, intent( in ) :: type
      !
      ! write(*,*) "Constructor TransmitterCSEM_t"
      !
      call self%init()
      !
      self%id = id
      self%n_pol = 1
      self%period = period
      self%location = location
      !
      if( present( type ) ) then
         self%type = type
      else
         self%type = "TransmitterCSEM_t"
      endif
      !
   end function TransmitterCSEM_ctor
   !
   ! Destructor
   subroutine TransmitterCSEM_dtor( self )
      implicit none
      !
      type( TransmitterCSEM_t )   :: self
      !
      ! write(*,*) "Destructor TransmitterCSEM_t"
      !
      call self%dealloc()
      !
   end subroutine TransmitterCSEM_dtor
   !
   subroutine solveFWDTransmitterCSEM( self )
      !
      class( TransmitterCSEM_t ), intent( inout ) :: self
      !
      write(*,*) "implementing solveFWD TransmitterCSEM_t: ", self%id
      !
      !call self%e_solution%add( self%forward_solver%getESolution( 333, 1, 1, self%period  ) )
      !
   end subroutine solveFWDTransmitterCSEM
   !
   !
   subroutine getSourceTransmitterCSEM( self )
      !
      class( TransmitterCSEM_t ), intent(in)   :: self
      !
      write(*,*) "getSource TransmitterCSEM_t: ", self%location
      !
   end subroutine getSourceTransmitterCSEM
   !
   !
   function isEqualTransmitterCSEM( self, other ) result( equal )
      class( TransmitterCSEM_t ), intent( in )   :: self
      class( Transmitter_t ), intent( in )   :: other
      logical                           :: equal
      !
      equal = .FALSE.
      !
      select type( other )
         !
         class is ( TransmitterCSEM_t )
            !
            if( ABS( self%period - other%period ) < TOL6 .AND.   &
               self%location(1) == other%location(1) .AND.   &
               self%location(2) == other%location(2) .AND.   &
               self%location(3) == other%location(3) ) then
               !
               equal = .TRUE.
            endif
      !
      end select
      !
   end function isEqualTransmitterCSEM
   !
   subroutine writeTransmitterCSEM( self )
      !
      class( TransmitterCSEM_t ), intent(in)   :: self
      integer                           :: iRx, nRx
      !
      nRx = self%getNRx()
      !
      write(*,*) "Write TransmitterCSEM_t Id: ", self%id,   &
      ", Period: ",   self%period,   &
      ", Location: ",   self%location,   &
      "fwd_key: ",   self%fwd_key(1), self%fwd_key(2), self%fwd_key(3), self%fwd_key(4),   &
                  self%fwd_key(5), self%fwd_key(6), self%fwd_key(7), self%fwd_key(8),   &
      " N Receivers: ", nRx
      !
      do iRx = 1, nRx
         write(*,*) "   ", self%get( iRx )
      enddo
      !
   end subroutine writeTransmitterCSEM
   !
   !function sizeOfTransmitterCSEM( self ) result( size )
      !
      !class( TransmitterCSEM_t ), intent( in ) :: self
      !integer                                :: size
	  !
	  !size = sizeof( self%id ) + &
	         !sizeof( self%n_pol ) + &
			 !sizeof( self%fwd_key ) + &
			 !sizeof( self%type ) + &
	         !sizeof( self%period ) + &
			 !sizeof( self%forward_solver ) + &
			 !sizeof( self%e_all ) + &
			 !sizeof( self%receiver_indexes ) + &
			 !sizeof( self%DATA_TITLE )
			 !sizeof( self%location ) + &
			 !sizeof( self%azimuth )
			 !
	  !write( *, * ) "Size: ", size
	  !
   !end function sizeOfTransmitterCSEM
   !
end module TransmitterCSEM
