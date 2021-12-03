!*************
!
! Derived class to hold all the CSEM information from a valid data file line
!
! Last modified at 03/2021 by Paulo Werdt
!
!*************
!
module DataEntryCSEM
   !
   use DataEntry
   !
   type, extends( DataEntry_t ) :: DataEntryCSEM_t
      !
      character(:), allocatable :: dipole
      real( kind=prec )         :: moment, azimuth, dip, tx_xyz(3)
      !
   contains
      !
      procedure, public   :: write => writeDataEntryCSEM
	  procedure, public :: getCopy => getCopyDataEntryCSEM
      !
   end type DataEntryCSEM_t
   !
   interface DataEntryCSEM_t
      module procedure DataEntryCSEM_ctor
   end interface DataEntryCSEM_t
   !
contains
   !
   ! Parametrized constructor
   function DataEntryCSEM_ctor( id, type, dipole, period, moment, azimuth, &
      dip, tx_xyz, code, xyz, component, real, imaginary, error ) result( self )
      implicit none
      type( DataEntryCSEM_t ) :: self
      integer, intent( in )               :: id
      character(:), allocatable, intent( in )   :: type, dipole, code, component
      real( kind=prec ), intent( in )        :: period, moment, azimuth, dip, xyz(3), tx_xyz(3)
      real( kind=prec ), intent( in )         :: real, imaginary, error
      !
      !write(*,*) "Constructor DataEntryCSEM_t"
      !
      self%id = id
      self%type = type
      self%dipole = dipole
      self%period = period
      self%moment = moment
      self%azimuth = azimuth
      self%dip = dip
      self%tx_xyz = tx_xyz
      self%code = code
      self%xyz = xyz
      self%component = component
      self%real = real
      self%imaginary = imaginary
      self%error = error
      !
   end function DataEntryCSEM_ctor
   !
   function getCopyDataEntryCSEM( self ) result ( copy )
      implicit none
      !
      class( DataEntryCSEM_t ), intent( in ) :: self
	  class( DataEntry_t ), allocatable      :: copy
      !
      allocate( copy, source = DataEntryCSEM_t( self%id, self%type,   &
                self%dipole, self%period, self%moment, self%azimuth, self%dip, self%tx_xyz,   &
                self%code, self%xyz, self%component, self%real, self%imaginary, self%error ) )
      !
   end function getCopyDataEntryCSEM
   !
   subroutine writeDataEntryCSEM( self )
      class( DataEntryCSEM_t ), intent( in ) :: self
      !
      write(*,*) "Write DataEntryCSEM_t: ", self%id, self%type, self%dipole, self%period,      &
      self%moment, self%azimuth, self%dip, self%tx_xyz, self%code, self%xyz, self%component,   &
      self%real, self%imaginary, self%error
      !
   end subroutine writeDataEntryCSEM
   !
end module DataEntryCSEM
