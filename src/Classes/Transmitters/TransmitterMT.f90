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
   use FileUnits
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
   function TransmitterMT_ctor( period, type ) result ( self )
      implicit none
      !
      type( TransmitterMT_t ) :: self
      !
      real( kind=prec ), intent( in ) :: period
      character(:), allocatable, optional, intent( in ) :: type
      !
      !write(*,*) "Constructor TransmitterMT_t"
      !
      call self%init()
      !
      self%n_pol = 2
      !
      allocate( cVector3D_SG_t :: self%e_all( self%n_pol ) )
      !
      self%period = period
      !
      if( present( type ) ) then
         self%type = type
      else
         self%type = "TransmitterMT_t"
      endif
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
      class( cVector_t ), allocatable :: e_solution
      !
      ! LATER SOUBROUTINE
      integer                 :: ioNum,iFreq,iMode
      character (len=20)     :: ModeName
      !
      ! verbosis
      write( *, * ) "   SolveFWD for Tx", self%id
      !
      omega = 2.0 * PI / self%period
      !
      open( ioESolution, file = 'e_solution', action='write', position='append', form ='unformatted' )
      !
      ! Loop over all polarizations (MT n_pol = 2)
      do i_pol = 1, self%n_pol
         !
         write(*,*) "MT Tx ", self%id, " Solve for Polarization", i_pol
         !
         ! Set Source E
         !   DO NOT WANT TO MODIFY SOURCE INSIDE FWDsolve (for inversion!)
         !    I guess this means we need to make some changes to source objects ...
         call self%source%setE( omega, i_pol )
         !
         if( allocated( e_solution ) ) deallocate( e_solution )
         allocate( e_solution, source = self%source%model_operator%createVector() )
         !
         call self%forward_solver%getESolution( self%source, e_solution )
         !
         if( i_pol == 1 ) then
            ModeName = "Ey"
         else
            ModeName = "Ex"
         endif
         !
         ! write the frequency header - 1 record
         write( ioESolution ) omega, self%id, i_pol, ModeName
         !
         call e_solution%write( ioESolution )
         !
         ! Add polarization e_solution to self%e_all
         self%e_all( i_pol ) = e_solution
         !
      enddo
      !
      close( ioESolution )
      !
   end subroutine solveFWDTransmitterMT
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
      integer :: iRx
      !
      write(*,*) "Write TransmitterMT_t: ", self%id,   &
      " Period: ",   self%period,   &
      "fwd_key: ",   self%fwd_key(1), self%fwd_key(2), self%fwd_key(3), self%fwd_key(4),   &
                  self%fwd_key(5), self%fwd_key(6), self%fwd_key(7), self%fwd_key(8),   &
      " N Receivers: ", size( self%receiver_indexes )
      !
      !do iRx = 1, size( self%receiver_indexes )
         !
         !write(*,*) self%receiver_indexes( iRx )
         !
      !enddo
      !
   end subroutine writeTransmitterMT
   !
end module TransmitterMT
