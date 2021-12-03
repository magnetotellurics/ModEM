! *************
! 
! Base class to define a Data Group
! 
! Last modified at 09/06/2021 by Paulo Werdt
! 
! *************
! 
module DataGroup
   !
   use Constants
   !
   type :: DataGroup_t
      !
      integer                                      :: id, n_data, id_rx, id_tx
      character(30), allocatable                   :: components(:)
      logical                                      :: is_complex
	  real( kind=prec ), dimension(:), allocatable :: reals, imaginaries, errors ! All Together
      !
      integer, private :: counter
	  !
   contains
      !
      final :: DataGroup_dtor
      !
      procedure, public :: add     => addDataDg
      procedure, public :: write   => writeDataGroup
      !
      procedure, public :: isEqual => isEqualDg
      !
   end type DataGroup_t
   !
   interface DataGroup_t
      module procedure DataGroup_ctor
   end interface DataGroup_t
   !
contains
   !
   ! Parametrized constructor
   function DataGroup_ctor( id, n_data ) result ( self )
      implicit none
      !
      type( DataGroup_t ) :: self
      !
      integer, intent( in )               :: id, n_data
      !
      self%id = id
      self%n_data = n_data
      !
      self%counter = 1
      !
      allocate( character(30) :: self%components( n_data ) )
      !
      allocate( self%reals( n_data ) )
      allocate( self%imaginaries( n_data ) )
      allocate( self%errors( n_data ) )
      !
   end function DataGroup_ctor
   !
   subroutine DataGroup_dtor( self )
      implicit none
      !
      type( DataGroup_t ), intent( in out ) :: self
      !
      !write(*,*) "Destructor DataGroup_t: ", self%id
      !
      deallocate( self%reals )
      deallocate( self%imaginaries )
      deallocate( self%errors )
      !
   end subroutine DataGroup_dtor
   !
   subroutine addDataDg( self, component, real, imaginary, error )
      implicit none
      class( DataGroup_t ), intent( inout )   :: self
      character(:), allocatable, intent( in )   :: component
      real( kind=prec ), intent( in )         :: real, imaginary, error
      !
      self%components( self%counter ) = component
      !
      self%reals( self%counter ) = real
      !
      self%imaginaries( self%counter ) = imaginary
      !
      self%errors( self%counter ) = error
      !
      self%counter = self%counter + 1
      !
   end subroutine addDataDg
   !
   subroutine writeDataGroup( self )
      class( DataGroup_t ), intent( in )   :: self
      integer                        :: i_data
      !
      write(*,*) "   Write DataGroup_t Id: ", self%id
      write(*,*) "      Receiver Id:   ", self%id_rx
      write(*,*) "      Transmitter Id:   ", self%id_tx
      !
      do i_data = 1, self%n_data
         !
         write(*,*) i_data, ":", self%components( i_data ), self%reals( i_data ), self%imaginaries( i_data ), self%errors( i_data )
         !
      enddo
      !
   end subroutine writeDataGroup
   !
   function isEqualDg( self, other ) result ( equal )
      class( DataGroup_t ), intent( in )   :: self
      class( DataGroup_t ), intent( in )   :: other
      logical                         :: equal
      !
      integer                        :: i_data
      !
      equal = .true.
      !
      if( self%id_tx /= other%id_tx .OR.   &
         self%id_rx /= other%id_rx ) then
         equal = .false.
         return
      end if
      !
      do i_data = 1, self%n_data
         !
         if( self%components( i_data ) /= other%components( i_data ) .OR.               &
            ABS( self%reals( i_data ) - other%reals( i_data ) ) >= TOL6 .OR.            &
            ABS( self%imaginaries( i_data ) - other%imaginaries( i_data ) ) >= TOL6 .OR.   &
            ABS( self%errors( i_data ) - other%errors( i_data ) ) >= TOL6 ) then
               equal = .false.
               exit
         end if
         !
      enddo
      !
   end function isEqualDg
   !
end module DataGroup
