!*************
!
! Base class to define a Receiver
!
! Last modified at 20/10/2021 by Paulo Werdt
!
!*************
!
module Receiver
   !
   use Transmitter
   use cVector3D_SG
   use ModelOperator
   use DataGroupArray
   use DataEntryArray
   use Grid3D_SG
   !
   type, abstract :: Receiver_t
      !
      integer                            :: id, n_comp
      !
      character(:), allocatable          :: code
      !
      real ( kind=prec )                 :: location(3)
      !
      logical                            :: is_complex = .FALSE., interpolation_set = .FALSE.
      !
      class( Grid_t ), pointer           :: grid
      !
      character(2), allocatable          :: EHxy(:)
      !
      character(3), allocatable          :: comp_names(:)
      !
      complex( kind=prec ), allocatable    :: I_BB(:,:), EE(:,:), Z(:)
      !
      class( cVector_t ), allocatable    :: Lex, Ley, Lez, Lbx, Lby, Lbz
      !
      class( DataGroupArray_t ), pointer :: data_groups
      !
	  class( DataEntryArray_t ), pointer :: predicted_data_entries
	  !
   contains
      !
      ! DEFERRED INTERFACES
      procedure( interface_predicted_data ), deferred, public :: predictedData
      !
	  procedure( interface_save_predicted_data_rx ), deferred, public  :: savePredictedData
      procedure( interface_write_predicted_data_rx ), deferred, public :: writePredictedData
      procedure( interface_write_rx ), deferred, public                :: write
      !
      ! CLASS PROCEDURES
      procedure, public :: evaluationFunction => evaluationFunctionRx
      !
      procedure, public :: init => initializeRx
      procedure, public :: dealloc => deallocateRx
      !
      procedure, public :: isEqual => isEqualRx
      !
      procedure, public :: has => hasDataGroupRx
      procedure, public :: add => addDataGroupRx
      procedure, public :: get => getDataGroupRx
      procedure, public :: getNdg => getNumberOfDataGroupRx
      !
   end type Receiver_t
   !
   abstract interface
      !
      subroutine interface_predicted_data( self, model_operator, transmitter )
         !
         import :: Receiver_t, ModelOperator_t, Transmitter_t
         !
         class( Receiver_t ), intent( inout )  :: self
         class( ModelOperator_t ),intent( in ) :: model_operator
         class( Transmitter_t ), intent( in )  :: transmitter
         !
      end subroutine interface_predicted_data
      !
      subroutine interface_save_predicted_data_rx( self, tx )
         !
         import :: Receiver_t, Transmitter_t
         !
         class( Receiver_t ), intent( in )    :: self
         class( Transmitter_t ), intent( in ) :: tx
         !
      end subroutine interface_save_predicted_data_rx
      !
      subroutine interface_write_predicted_data_rx( self )
         !
         import :: Receiver_t
         !
         class( Receiver_t ), intent( in ) :: self
         !
      end subroutine interface_write_predicted_data_rx
      !
      subroutine interface_write_rx( self )
         !
         import :: Receiver_t
         !
         class( Receiver_t ), intent(in) :: self
         !
      end subroutine interface_write_rx
      !
   end interface
   !
contains
   !
   subroutine initializeRx( self )
      implicit none
      !
      class( Receiver_t ), intent( inout ) :: self
      !
      self%id = 0
      !
      self%n_comp = 0
      !
      self%location = 0.0
      !
      self%grid => null()
      !
      allocate( self%data_groups, source = DataGroupArray_t() )
      !
	  allocate( self%predicted_data_entries, source = DataEntryArray_t() )
      !
   end subroutine initializeRx
   !
   subroutine deallocateRx( self )
      implicit none
      !
      class( Receiver_t ), intent( inout ) :: self
      !
      if( associated( self%grid ) ) deallocate( self%grid )
      !
      if( allocated( self%EHxy ) ) deallocate( self%EHxy )
      !
      if( allocated( self%I_BB ) ) deallocate( self%I_BB )
      !
      if( allocated( self%EE ) ) deallocate( self%EE )
      !
      if( allocated( self%Z ) ) deallocate( self%Z )
      !
      if( allocated( self%Lex ) ) deallocate( self%Lex )
      !
      if( allocated( self%Ley ) ) deallocate( self%Ley )
      !
      if( allocated( self%Lez ) ) deallocate( self%Lez )
      !
      if( allocated( self%Lbx ) ) deallocate( self%Lbx )
      !
      if( allocated( self%Lby ) ) deallocate( self%Lby )
      !
      if( allocated( self%Lbz ) ) deallocate( self%Lbz )
      !
      !deallocate( self%data_groups )
      !
   end subroutine deallocateRx
   !
   subroutine evaluationFunctionRx( self, model_operator, omega )
      implicit none
      !
      class( Receiver_t ), intent( inout )   :: self
      class( ModelOperator_t ), intent( in ) :: model_operator
      real( kind=prec ), intent( in )        :: omega
      !
      integer              :: k
      complex( kind=prec ) :: comega
      !
      ! THESE SHOULD BE rVECTORS ????
      class( cVector_t ), allocatable :: e, h
      class( cVector_t ), allocatable :: lh
      !
      comega = cmplx( 1./omega, 0.0, kind=prec )
      !
      do k = 1, self%n_comp
         !
         selectcase( self%EHxy(k) )
            !
            case( "Ex", "Ey" )
               !
                select type( grid => model_operator%grid )
                    class is( Grid3D_SG_t )
                        if( .not. allocated( e ) ) allocate( e, source = cVector3D_SG_t( grid, EDGE ) )
                    class default
                        stop "Receiver: Unclassified model_operator%grid for e"
                end select
               !
            case( "Bx", "By", "Bz" )
               !
                select type( grid => model_operator%grid )
                    class is( Grid3D_SG_t )
                        if( .not. allocated(h) ) allocate( h, source = cVector3D_SG_t( grid, FACE ) )
                    class default
                        stop "Receiver: Unclassified model_operator%grid for h"
                end select
               !
         end select
         !
      end do
      !
      do k = 1, self%n_comp
         !
         if( allocated( lh ) ) deallocate( lh )
         !
         selectcase( self%EHxy(k) )
            !
            case( "Ex" )
               call e%interpFunc( self%location, "x", self%Lex )
            !
            case( "Ey" )
               call e%interpFunc( self%location, "y", self%Ley )
            !
            case( "Ez" )
               call e%interpFunc( self%location, "z", self%Lez )
            !
            case( "Bx" )
			   !
               call h%interpFunc( self%location, "x", lh )
			   !
			   select type( lh )
                  class is(cVector3D_SG_t)
			         if( allocated( self%Lbx ) ) deallocate( self%Lbx )
			         allocate( self%Lbx, source = cVector3D_SG_t( lh%grid, EDGE ) )
					 !
                  class default
                     write(*, *) 'ERROR:Receiver::evaluationFunction:'
                     stop        '         Unkonow lh type'
               end select
               !
               call model_operator%multCurlT( lh, self%Lbx )
               call self%Lbx%mults( isign * ONE_I / comega )
               !
            case( "By" )
			   ! 
               call h%interpFunc( self%location, "y", lh )
			   !
			   select type( lh )
                  class is(cVector3D_SG_t)
			         if( allocated( self%Lby ) ) deallocate( self%Lby )
			         allocate( self%Lby, source = cVector3D_SG_t( lh%grid, EDGE ) )
					 !
                  class default
                     write(*, *) 'ERROR:Receiver::evaluationFunction:'
                     stop        '         Unkonow lh type'
               end select
               !
               call model_operator%multCurlT( lh, self%Lby )
               call self%Lby%mults( isign * ONE_I / comega )
            !
            case( "Bz" )
			   !
               call h%interpFunc( self%location, "z", lh )
			   !
			   select type( lh )
                  class is(cVector3D_SG_t)
			         if( allocated( self%Lbz ) ) deallocate( self%Lbz )
			         allocate( self%Lbz, source = cVector3D_SG_t( lh%grid, EDGE ) )
					 !
                  class default
                     write(*, *) 'ERROR:Receiver::evaluationFunction:'
                     stop        '         Unkonow lh type'
               end select
               !
               call model_operator%multCurlT( lh, self%Lbz )
               call self%Lbz%mults( isign * ONE_I / comega )
            !
         end select
         !
      end do
      !
      if( allocated( e ) ) deallocate( e )
      if( allocated( h ) ) deallocate( h )
      !
   end subroutine evaluationFunctionRx
   !
   function isEqualRx( self, other ) result( equal )
      implicit none
      !
      class( Receiver_t ), intent( in ) :: self
      class( Receiver_t ), intent( in ) :: other
      !
      logical :: equal
      !
      equal = .FALSE.
      !
      if( self%location(1) == other%location(1) .AND.   &
          self%location(2) == other%location(2) .AND.   &
          self%location(3) == other%location(3) ) then
         equal = .TRUE.
      endif
      !
   end function isEqualRx
   !
   function hasDataGroupRx( self, data_group ) result( found )
      implicit none
      !
      class( Receiver_t ), intent( in )  :: self
      class( DataGroup_t ), intent( in ) :: data_group
      !
      logical :: found
      integer :: iDg, nDg
      !
      found = .FALSE.
      !
      nDg = self%data_groups%size()
      !
      do iDg = 1, nDg
         !
         if( data_group%isEqual( self%data_groups%get( iDg ) ) ) then
            found = .TRUE.
         end if
      end do
      !
   end function hasDataGroupRx
   !
   subroutine addDataGroupRx( self, data_group )
      implicit none
      !
      class( Receiver_t ), intent( inout ) :: self
      class( DataGroup_t ), intent( in )   :: data_group
      !
      call self%data_groups%add( data_group )
      !
   end subroutine addDataGroupRx
   !
   function getDataGroupRx( self, index ) result( data_group )
      implicit none
      !
      class( Receiver_t ), intent( in ) :: self
      integer, intent( in )             :: index
      class( DataGroup_t ), allocatable :: data_group
      !
      data_group = self%data_groups%get( index )
      !
   end function getDataGroupRx
   !
   function getNumberOfDataGroupRx( self ) result( counter )
      implicit none
      !
      class( Receiver_t ), intent( in ) :: self
      integer                           :: counter
      !
      counter = self%data_groups%size()
      !
   end function getNumberOfDataGroupRx
   !
end module Receiver
