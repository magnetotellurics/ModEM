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
    use Grid3D_SG
    use cVector3D_SG
    use ModelOperator
    !
    type, extends( ForwardSolver_t ), public :: ForwardSolverFromFile_t
        !
		real( kind=prec ) :: period
		!
        class( ModelOperator_t ), pointer  :: model_operator
        !
        contains
            !
            final :: ForwardSolverFromFile_dtor
            !
            procedure, public :: setIterControl  => SetIterControlForwardSolverFromFile
            procedure, public :: initDiagnostics => initDiagnosticsForwardSolverFromFile
            procedure, public :: zeroDiagnostics => zeroDiagnosticsForwardSolverFromFile
            procedure, public :: setFrequency    => setFrequencyForwardSolverFromFile
            procedure, public :: setCond         => setCondForwardSolverFromFile
            procedure, public :: getESolution    => getESolutionForwardSolverFromFile
            !
    end type ForwardSolverFromFile_t
    !
    interface ForwardSolverFromFile_t
        module procedure ForwardSolverFromFile_ctor
    end interface ForwardSolverFromFile_t
    !
contains
    !
    function ForwardSolverFromFile_ctor( model_operator ) result( self )
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        type( ForwardSolverFromFile_t )                :: self
        !
        !write(*,*) "Constructor ForwardSolverFromFile_t"
        !
        self%model_operator => model_operator
        !
    end function ForwardSolverFromFile_ctor
    !
    ! Destructor
    subroutine ForwardSolverFromFile_dtor( self )
        implicit none
        !
        type( ForwardSolverFromFile_t ), intent( inout ) :: self
        !
        !write(*,*) "Destructor ForwardSolverFromFile_t"
        !
    end subroutine ForwardSolverFromFile_dtor
    !
    subroutine setFrequencyForwardSolverFromFile( self, period )
        implicit none
        !
        class( ForwardSolverFromFile_t ), intent( inout ) :: self
        real( kind=prec ), intent( in )                   :: period
        !
        self%period = period
    !
    end subroutine setFrequencyForwardSolverFromFile
    !
    !
    subroutine setCondForwardSolverFromFile( self, model_parameter )
        implicit none
          !
          class( ForwardSolverFromFile_t ), intent( inout ) :: self
          class( ModelParameter_t ), intent( in )           :: model_parameter
          !
    end subroutine setCondForwardSolverFromFile
    !
    !
    subroutine zeroDiagnosticsForwardSolverFromFile(self)
        implicit none
        !
        class( ForwardSolverFromFile_t ), intent( inout ) :: self
        !
        self%relResVec = R_ZERO
        call self%solver%zeroDiagnostics()
        !
    end subroutine zeroDiagnosticsForwardSolverFromFile
    !
    subroutine initDiagnosticsForwardSolverFromFile( self )
        implicit none
        !
        class( ForwardSolverFromFile_t ), intent( inout ) :: self
        !
        self%n_iter_actual = 0
        self%relResFinal = R_ZERO
        !
        if(allocated(self%relResVec)) deallocate(self%relResVec)
        allocate(self%relResVec(self%max_iter_total))
        !
    end subroutine initDiagnosticsForwardSolverFromFile
    !
    subroutine setIterControlForwardSolverFromFile( self, maxit, tolerance )
        implicit none
        !
        class( ForwardSolverFromFile_t ), intent( inout )  :: self
        integer, intent(in)                                :: maxit
        real(kind=prec), intent(in)                        :: tolerance
        !
        self%max_iter_total = maxit
        self%tolerance = tolerance
        !
        call self%solver%setParameters( maxit, tolerance )
        !
    end subroutine setIterControlForwardSolverFromFile
    !
    !**********
    !
    subroutine getESolutionForwardSolverFromFile( self, source, e_solution )
        implicit none
        !
        class( ForwardSolverFromFile_t ), intent( inout ) :: self
        class( Source_t ), intent( in )                   :: source
        class( cVector_t ), intent( inout )               :: e_solution
        !
        character(80) :: grid_type, file_name
        complex       :: x, y, z
        integer       :: nx, ny, nz, io_stat, ioE
        !
        ioE = 2022
        !
        ! Construct the file name
        write ( file_name, '(a,I4.4,a,I1,a)' ) '../inputs/esol/E_solution_Per', int( self%period ), '_Pol', source%polarization, '.soln'
        !
        ! Open the File Unit for the file name
        open ( unit=ioE, file=file_name, status='unknown', form ='unformatted', iostat=io_stat )
        !
        ! If the file not exist
        if( io_stat /= 0 ) then
            write(*,*) "Unable to open [", trim(file_name), "], Stat: ", io_stat
        else
            !
            ! Allocate e_solution based on the grid
            call e_solution%Read( ioE )
            !
            ! Close the Unit
            close( ioE )
            !
            ! Increase the Unit
            ioE = ioE + 1
            !
            ! Print the e_solution result
            write( *, * ) "     Polarization:", source%polarization
            !
            select type( e_solution )
                class is( cVector3D_SG_t )
                    write( *, * ) "            ", e_solution%nx, e_solution%ny, e_solution%nz, e_solution%gridType
                class default
                    stop "Unclassified ForwardSolverFromFile e_solution"
            end select
        endif
        !
    end subroutine getESolutionForwardSolverFromFile
    !
end module ForwardSolverFromFile
