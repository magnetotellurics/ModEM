! *************
! 
! Derived class to define a MT ForwardSolver that reads e solution from a file
! 
! Last modified at 19/10/2021 by Paulo Werdt
! 
! *************
! 
module ForwardSolverFromFile
   !
   use ForwardSolver
   use cVector3D_SG
   !
   type, extends( ForwardSolver_t ), public :: ForwardSolverFromFile_t
      !
      character(:), allocatable :: file_name
      integer                   :: IoE
      !
      contains
         !
         final :: ForwardSolverFromFile_dtor
         !
         procedure, public :: setPeriod => setPeriodForwardSolverFromFile
         procedure, public :: getESolution => getESolutionForwardSolverFromFile
         !
   end type ForwardSolverFromFile_t
   !
   interface ForwardSolverFromFile_t
      module procedure ForwardSolverFromFile_ctor
   end interface ForwardSolverFromFile_t
   !
contains
   !
   function ForwardSolverFromFile_ctor() result( self )
      !
      type( ForwardSolverFromFile_t ) :: self
      !
      !write(*,*) "Constructor ForwardSolverFromFile_t"
      !
      call self%init()
      !
      self%IoE = 901
      !
   end function ForwardSolverFromFile_ctor
   !
   ! Destructor
   subroutine ForwardSolverFromFile_dtor( self )
      implicit none
      !
      type( ForwardSolverFromFile_t ), intent( in out ) :: self
      !
      !write(*,*) "Destructor ForwardSolverFromFile_t"
      !
      call self%dealloc()
      !
   end subroutine ForwardSolverFromFile_dtor
   !
   subroutine setPeriodForwardSolverFromFile( self, period )
      implicit none
      !
      class( ForwardSolverFromFile_t ), intent( inout ) :: self
      real( kind=prec ), intent( in )                   :: period
      !
      self%period = period
   !
   end subroutine setPeriodForwardSolverFromFile
   !
   function getESolutionForwardSolverFromFile( self, source, polarization ) result( e_solution )
      implicit none
      !
      class( ForwardSolverFromFile_t ), intent( inout ) :: self
      class( Source_t ), intent( in )                   :: source
      integer, intent( in )                                :: polarization
      !
      class( cVector_t ), allocatable :: e_solution
      !
      character(80) :: grid_type, file_name
      complex          :: x, y, z
      integer          :: nx, ny, nz, io_stat
   !
      ! Construct the file name
      write ( file_name, '(a,I4.4,a,I1,a)' ) '../inputs/esol/E_solution_Per', int( self%period ), '_Pol', polarization, '.soln'
      !
      ! Open the File Unit for the file name
      open ( unit=self%ioE, file=file_name, status='unknown', form ='unformatted', iostat=io_stat )
      !
      ! If the file not exist
      if( io_stat /= 0 ) then
            write(*,*) "Unable to open [", trim(file_name), "], Stat: ", io_stat
      else
            !
            ! Save the file name
            self%file_name = trim( file_name )
            !
            ! Read and Save cVector e_solution
            call self%e_solution%Read( self%ioE )
            !
            ! Close the Unit
            close( self%ioE )
            !
            ! Increase the Unit
            self%ioE = self%ioE + 1
            !
            ! Set the e_solution result
            allocate( e_solution, source = self%e_solution )
            !
            ! Print the e_solution result
            write( *, * ) "    Polarization:", polarization
            select type( e_solution )
                  class is( cVector3D_SG_t )
                        write( *, * ) "         ", e_solution%nx, e_solution%ny, e_solution%nz, e_solution%gridType
                  class default
                        stop "Unclassified e_solution"
            end select
            !
      endif
      !
   end function getESolutionForwardSolverFromFile
   !
   !
   subroutine defineSource( self )
      !
      class( ForwardSolverFromFile_t ), intent(in)   :: self
      !
      write(*,*) "defineSource ForwardSolverFromFile: "
      !
   end subroutine defineSource
   !
end module ForwardSolverFromFile
