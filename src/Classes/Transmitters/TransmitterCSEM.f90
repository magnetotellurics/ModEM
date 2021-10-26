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
      real( kind=prec )   :: location(3)
      !
      real( kind=prec )   :: azimuth
      !
      contains
         !
         final   :: TransmitterCSEM_dtor
         !
         procedure, public   :: solveFWD => solveFWDTransmitterCSEM
         procedure, public   :: getSource => getSourceTransmitterCSEM
         !
		 procedure, public   :: getType => getTypeTransmitterCSEM
         procedure, public   :: isEqual => isEqualTransmitterCSEM
         procedure, public   :: write => writeTransmitterCSEM
         !
   end type TransmitterCSEM_t
   !
   interface Transmitter_t
      module procedure TransmitterCSEM_ctor
   end interface Transmitter_t
   !
contains
   !
   ! Parametrized constructor
   function TransmitterCSEM_ctor( id, period, location ) result ( self )
      !
      class( TransmitterCSEM_t ), pointer   :: self
      !
      integer, intent( in )            :: id
      real( kind=prec ), intent( in )      :: period
      real( kind=prec )               :: location(3)
      !
      ! write(*,*) "Constructor TransmitterCSEM_t"
      !
      allocate( TransmitterCSEM_t :: self )
      !
      call self%init()
      !
      self%id = id
      self%n_pol = 1
      self%period = period
      self%location = location
      !
      ! INSTANCIATE SOURCE CSEM ????
      !
      !source => Source_CSEM_t()
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
   subroutine solveFWDTransmitterCSEM( self, model_operator )
      !
      class( TransmitterCSEM_t ), intent(inout)  :: self
      class( ModelOperator_t ), allocatable, intent( in ) :: model_operator
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
   function getTypeTransmitterCSEM( self ) result( type )
      class( TransmitterCSEM_t ), intent( in ) :: self
      character(:), allocatable                :: type
      !
      type = "TransmitterCSEM_t"
      !
   end function getTypeTransmitterCSEM
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
end module TransmitterCSEM
