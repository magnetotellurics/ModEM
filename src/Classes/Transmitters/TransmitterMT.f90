! *************
! 
! Derived class to define a MT Transmitter
! 
! Last modified at 08/06/2021 by Paulo Werdt
! 
! *************
! 
module TransmitterMT
   ! 
   use Transmitter
   use Grid3D_SG
   use cVector3D_SG
   use Source_MT
   !
   type, extends( Transmitter_t ), public :: TransmitterMT_t
      !
      ! source polarizations for MT
      real( kind=prec ) :: pol! = {'X','Y'}
      ! 
      contains
         !
         final :: TransmitterMT_dtor
         !
         procedure, public :: solveFWD => solveFWDTransmitterMT
         procedure, public :: getSource => getSourceTransmitterMT
         !
         procedure, public :: getType => getTypeTransmitterMT
         procedure, public :: isEqual => isEqualTransmitterMT
         procedure, public :: write => writeTransmitterMT
         !
   end type TransmitterMT_t
   !
   interface Transmitter_t
      module procedure TransmitterMT_ctor
   end interface Transmitter_t
   !
contains
   !
   function TransmitterMT_ctor( id, period ) result ( self )
      !
      class( TransmitterMT_t ), pointer :: self
      integer, intent( in )             :: id
      real( kind=prec ), intent( in )   :: period
      !
      ! write(*,*) "Constructor TransmitterMT_t"
      !
      allocate( TransmitterMT_t :: self )
      !
      call self%init()
      !
      self%id = id
      self%n_pol = 2
      self%period = period
      !
      self%forward_solver => ForwardSolverFromFile_t( self%period, self%n_pol )
      !
     self%DATA_TITLE = "Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
     !
   end function TransmitterMT_ctor
   !
   ! Destructor
   subroutine TransmitterMT_dtor( self )
      implicit none
      !
      type( TransmitterMT_t ), intent( in out ) :: self
      !
      !write(*,*) "Destructor TransmitterMT_t"
      !
      call self%dealloc()
      !
   end subroutine TransmitterMT_dtor
   !
   subroutine solveFWDTransmitterMT( self, model_operator )
      !
      class( TransmitterMT_t ), intent(inout)  :: self
      class( ModelOperator_t ), allocatable, intent( in ) :: model_operator
      !
	  ! Define Source
      allocate( self%source, source=Source_MT_t( model_operator, 2.0 * PI / self%period, 'X' ) )
      !
	  ! Define Forward Solver ESolution
      select type( grid => model_operator%grid )
          class is( Grid3D_SG_t )
              allocate( self%forward_solver%e_solution, source=cVector3D_SG_t( grid, EDGE ) )
      end select
      !
      if ( .not. allocated( self%forward_solver%e_solution ) ) write(*, *) 'TX E SOLUTION NOT ALLOCATTED'
      if ( .not. self%forward_solver%e_solution%isAllocated ) write(*, *) 'TX E SOLUTION isAllocated = false'
      !
      ! Add first polarization to Tx e_solution
      call self%e_solution%add( self%forward_solver%getESolution( self%source ) )
      !
      ! Add second polarization to Tx e_solution
      call self%e_solution%add( self%forward_solver%getESolution( self%source ) )
      !
   end subroutine solveFWDTransmitterMT
   !
   !
   subroutine getSourceTransmitterMT( self )
      !
      class( TransmitterMT_t ), intent(in)   :: self
      !
      write(*,*) "To implement getSource TransmitterMT_t: ", self%id
      !
   end subroutine getSourceTransmitterMT
   !
   !
   function getTypeTransmitterMT( self ) result( type )
      class( TransmitterMT_t ), intent( in ) :: self
      character(:), allocatable              :: type
      !
      type = "TransmitterMT_t"
      !
   end function getTypeTransmitterMT
   !
   !
   function isEqualTransmitterMT( self, other ) result( equal )
      class( TransmitterMT_t ), intent( in )   :: self
      class( Transmitter_t ), intent( in )   :: other
      logical                           :: equal
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
   subroutine writeTransmitterMT( self )
      !
      class( TransmitterMT_t ), intent(in)   :: self
      integer                           :: iRx, nRx
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
