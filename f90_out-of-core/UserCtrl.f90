! *****************************************************************************
module UserCtrl
  ! This module defines the derived data type structure with all filenames
  use utilities

  implicit none

  character*1, parameter  :: READ_WRITE = 'R'
  character*1, parameter	:: FORWARD = 'F'
  character*1, parameter	:: COMPUTE_J = 'J'
  character*1, parameter	:: MULT_BY_J = 'M'
  character*1, parameter	:: MULT_BY_J_T = 'T'
  character*1, parameter	:: MULT_BY_J_T_multi_Tx = 'x'
  character*1, parameter	:: INVERSE = 'I'
  character*1, parameter	:: APPLY_COV = 'C'
  character*1, parameter    :: EXTRACT_BC = 'b'
  character*1, parameter  :: TEST_GRAD = 'g'
  character*1, parameter  :: TEST_ADJ = 'A'
  character*1, parameter  :: TEST_SENS = 'S'

  ! ***************************************************************************
  ! * userdef_control contains the list of all essential input information currently
  ! * read in from fn_startup.
  type :: userdef_control

	! Options: FORWARD, COMPUTE_J, MULT_BY_J, MULT_BY_J_T,
	!          INVERSE, TEST_COV, READ_WRITE
	character(80)       :: job

	! File to set up inversion controls
	character(80)       :: rFile_invCtrl

	! File to set up forward solver controls
	character(80)       :: rFile_fwdCtrl

    ! Output file name for MPI nodes status info
    character(80)       :: wFile_MPI

	! Input files
	character(80)       :: rFile_Grid, rFile_Model, rFile_Data
	character(80)       :: rFile_dModel
  character(80)       :: rFile_EMsoln, rFile_EMrhs, rFile_Prior

	! Output files
	character(80)       :: wFile_Grid, wFile_Model, wFile_Data
	character(80)       :: wFile_dModel
	character(80)       :: wFile_EMsoln, wFile_EMrhs, wFile_Sens
	! Specify covariance configuration
	character(80)       :: rFile_Cov

	! Choose the inverse search algorithm
	character(80)       :: search

  ! Choose the sort of test / procedure variant you wish to perform
  character(80)       :: option
     !  Out-of-core solution file prefix
     character(80)       :: prefix



	! Specify damping parameter for the inversion
	real(8)             :: lambda

	! Misfit tolerance for the forward solver
	real(8)             :: eps
  ! Specify the magnitude for random perturbations
  real(8)             :: delta


	! Indicate how much output you want
	integer             :: output_level

     ! Out of core solution boolean
     logical             :: storeSolnsInFile
  end type userdef_control

  type (userdef_control),public:: cUserDef
Contains

  ! Initialize userCtrl. The defaults are not empty strings.
  ! This is because the inquire statements with empty file names
  ! are not understood on IBM machines.

  subroutine initUserCtrl(ctrl)

  	type(userdef_control), intent(out)   :: ctrl

    character(8) date
    character(10) time
  	integer(4) pid

  	ctrl%job = ''
  	ctrl%rFile_invCtrl = 'n'
  	ctrl%rFile_fwdCtrl = 'n'
  	ctrl%rFile_Grid = 'n'
  	ctrl%wFile_Grid = 'n'
  	ctrl%rFile_Model = 'n'
  	ctrl%wFile_Model = 'n'
  	ctrl%rFile_Data = 'n'
  	ctrl%wFile_Data = 'n'
  	ctrl%rFile_dModel = 'n'
  	ctrl%wFile_dModel = 'n'
  	ctrl%rFile_EMrhs = 'n'
    ctrl%wFile_EMrhs = 'n'
    ctrl%rFile_EMsoln = 'n'
  	ctrl%wFile_EMsoln = 'n'
  	ctrl%rFile_Prior = 'n'
  	ctrl%wFile_Sens = 'n'
  	ctrl%lambda = 10.
  	ctrl%eps = 1.0e-7
  	ctrl%rFile_Cov = 'n'
  	ctrl%search = 'NLCG'
  	ctrl%option = 'J'
  	ctrl%delta = 0.05
  	ctrl%output_level = 3
!  output core ctrl
        ctrl%prefix = 'n'
        ctrl%storeSolnsInFile = .false.
    ! Using process ID in MPI output file name has the advantage that
    ! the user may run several instances of the program in one directory
    ! simultaneously. Unfortunately, getpid is one of the portability
    ! routines that is not universally supported. Use start time instead.
    ! pid = getpid()
    ! write(ctrl%wFile_MPI,'(a13,i6.6,a5)') 'Nodes_Status_',pid,'.info'

    call date_and_time(date,time)
    write(ctrl%wFile_MPI,'(a13,a8,a1,a10,a5)') 'Nodes_Status_',date,'T',time,'.info'

  end subroutine initUserCtrl

  subroutine parseArgs(program,ctrl)

     character(1)     :: job
     integer (kind=4) :: iargc,narg,k
     integer          :: istat
     logical          :: exists
     character*80  gridType, header,arg, verbose, paramtype
     character*80, dimension(:), pointer :: temp
     character(*), intent(in)            :: program
     type(userdef_control), intent(out)  :: ctrl
     logical :: res
     integer  :: nfarg
     write(*,*) 'Copyright (c) 2004-2014 Oregon State University'
     write(*,*) 'AUTHORS  Gary Egbert, Anna Kelbert & Naser Meqbel'
     write(*,*) 'College of Earth, Ocean and Atmospheric Sciences'
     write(*,*)

     call initUserCtrl(ctrl)

     !  parse command line ...  for now just uses first argument
     !   to set job
     narg = iargc()
     !  quick fix against compilers which add additional parameters to mpirun;
     !   only read the first argument with a dash, or the argument -v followed
     !   by the verbose level; ignore the others
     verbose = 'regular'
     nfarg=0
     if(narg .gt. 1) then
       k=1
       search_arg: &
		   do
		       k=k+1
                      if (k .eq. narg) then
                         exit  search_arg
                         endif
               call getarg ( k, arg )
                if(arg(1:1).eq.'-') then
                  if (arg(2:2).eq.'v') then
                     call getarg ( k+1, verbose )
                     nfarg=nfarg+2
                  exit  search_arg
               else if ( arg(2:2).eq.'P') then
                  call getarg ( k+1, ctrl%prefix )
                  ctrl%storeSolnsInFile = .true.
                  nfarg=nfarg+2
                  end if
                end if
          end do search_arg
     end if
     narg=narg-nfarg
	 !  set the level of output based on the user input
	 select case (verbose)
	 case ('debug')
	   print *,'Output all information including debugging lines.'
	   ctrl%output_level = 5
	 case ('full')
	   print *,'Output full information to screen and to files.'
	   ctrl%output_level = 4
	 case ('regular')
	   print *,'Output information to files, and progress report to screen (default).'
	   ctrl%output_level = 3
	 case ('compact')
	   print *,'Output information to files, and compact summary to screen.'
	   ctrl%output_level = 2
	 case ('result')
	   print *,'Output minimal information to files, and result to screen.'
	   ctrl%output_level = 1
	 case ('none')
	   print *,'Output nothing at all except result to screen and to files.'
	   ctrl%output_level = 0
	 case default
	   ctrl%output_level = 3
	 end select
     if(narg .gt. 0) then
        call getarg(1,arg)
        if(arg(1:1).eq.'-') job = arg(2:2)
     else
        write(*,*) 'Usage: ',trim(program),' -[job] [args]'
        write(*,*)
        write(*,*) '[READ_WRITE]'
        write(*,*) ' -R  rFile_Model rFile_Data [wFile_Model wFile_Data]'
        write(*,*) '  Reads your input files and checks them for validity;'
        write(*,*) '  optionally also writes them out'
        write(*,*) '[FORWARD]'
        write(*,*) ' -F  rFile_Model rFile_Data wFile_Data [wFile_EMsoln rFile_fwdCtrl]'
        write(*,*) '  Calculates the predicted data and saves the EM solution'
        write(*,*) '[INVERSE]'
        write(*,*) ' -I NLCG rFile_Model rFile_Data [lambda eps]'
        write(*,*) '  Here, lambda = the initial damping parameter for inversion'
        write(*,*) '           eps = misfit tolerance for the forward solver'
        write(*,*) 'OR'
        write(*,*) ' -I NLCG rFile_Model rFile_Data [rFile_invCtrl rFile_fwdCtrl]'
        write(*,*) '  Optionally, may also supply'
        write(*,*) '      the model covariance configuration file   [rFile_Cov]'
        write(*,*) '      the starting model parameter perturbation [rFile_dModel]'
        write(*,*) '  Runs an inverse search to yield an inverse model at every iteration'
        write(*,*) '[COMPUTE_J]'
        write(*,*) ' -J  rFile_Model rFile_Data wFile_Sens [rFile_fwdCtrl]'
        write(*,*) '  Calculates and saves the full J(acobian)'
        write(*,*) '[MULT_BY_J]'
        write(*,*) ' -M  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl]'
        write(*,*) '  Multiplies a model by J to create a data vector'
        write(*,*) '[MULT_BY_J_T]'
        write(*,*) ' -T  rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]'
        write(*,*) '  Multiplies a data vector by J^T to create a model'
        write(*,*) '[MULT_BY_J_T_multi_Tx]'
        write(*,*) ' -x  rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]'
        write(*,*) '  Multiplies a data vector by J^T to output models for each transmitter'
        write(*,*) '[APPLY_COV]'
        write(*,*) ' -C FWD rFile_Model wFile_Model [rFile_Cov rFile_Prior]'
        write(*,*) '  Applies the model covariance to produce a smooth model output'
        write(*,*) '  Optionally, also specify the prior model to compute resistivities'
        write(*,*) '  from model perturbation: m = C_m^{1/2} \\tilde{m} + m_0'
        write(*,*) '[TEST_GRAD]'
        write(*,*) ' -g  rFile_Model rFile_Data rFile_dModel [rFile_fwdCtrl rFile_EMrhs]'
        write(*,*) '  Runs the ultimate test of the gradient computation based on'
        write(*,*) '  Taylor series approximation.'
        write(*,*) '[TEST_ADJ]'
        write(*,*) ' -A  J rFile_Model rFile_dModel rFile_Data [wFile_Model wFile_Data]'
        write(*,*) '  Tests the equality d^T J m = m^T J^T d for any model and data.'
        write(*,*) '  Optionally, outputs J m and J^T d.'
        write(*,*) '[TEST_SENS]'
        write(*,*) ' -S  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl]'
        write(*,*) '  Multiplies by the full Jacobian, row by row, to get d = J m.'
        write(*,*) '  Compare to the output of [MULT_BY_J] to test [COMPUTE_J]'
        write(*,*)
        write(*,*) 'Optional final argument -v [debug|full|regular|compact|result|none]'
        write(*,*) 'indicates the desired level of output to screen and to files.'
        write(*,*)
!  out-of-core prefix
        write(*,*) '-P [file_prefix]'
        write(*,*) 'indicates that the partial solutions should be stored as '
        write(*,*) 'temporary files on disk rather than in memory with file prefix'
        write(*,*) '[file_prefix]'
        write(*,*)
        stop
     endif

     ! extract all following command line arguments
     allocate(temp(1:narg-1))
     do k = 1,narg-1
       call getarg(k+1,temp(k))
     end do

     ! write feedback to standard output for records
     if (narg > 1) then
        write(6,*) 'Running job -',job,' with command line options:'
        do k = 1,narg-1
          write(6,'(x1a)',advance='no') trim(temp(k))
        end do
        write(6,*)
        write(6,*)
     end if

     ! set narg to the number of actual arguments (after the job is specified)
     narg = narg-1

     select case (job)

      case (READ_WRITE) !R
        if (narg < 2) then
           write(0,*) 'Usage: -R  rFile_Model rFile_Data [wFile_Model wFile_Data]'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	    end if
	    if (narg > 2) then
	       ctrl%wFile_Model = temp(3)
	    end if
	    if (narg > 3) then
	       ctrl%wFile_Data = temp(4)
	    end if

      case (FORWARD) !F
        if (narg < 3) then
           write(0,*) 'Usage: -F  rFile_Model rFile_Data wFile_Data [wFile_EMsoln rFile_fwdCtrl]'
           write(0,*)
           write(0,*) 'Here, rFile_fwdCtrl is the forward solver control file in the format'
           write(0,*)
           write(0,*) 'Number of QMR iters per divergence correction : 40'
           write(0,*) 'Maximum number of divergence correction calls : 20'
           write(0,*) 'Maximum number of divergence correction iters : 100'
           write(0,*) 'Misfit tolerance for EM forward solver        : 1.0e-7'
           write(0,*) 'Misfit tolerance for EM adjoint solver        : 1.0e-7'
           write(0,*) 'Misfit tolerance for divergence correction    : 1.0e-5'
           write(0,*) 'Optional EM solution file name for nested BC  : nested.esoln'
           write(0,*)
           write(0,*) 'To specify air layers, append one of these three options. Default ''mirror 10 3. 30.'' '
           write(0,*)
           write(0,*) 'Option 1:'
           write(0,*) 'Air layers mirror|fixed height|read from file : mirror'
           write(0,*) 'Number of air layers and min top dz in km     : 10 3. 30.'
           write(0,*)
           write(0,*) 'Option 2:'
           write(0,*) 'Air layers mirror|fixed height|read from file : fixed height'
           write(0,*) 'Number of air layers and max height in km     : 12 1000.'
           write(0,*)
           write(0,*) 'Option 3:'
           write(0,*) 'Air layers mirror|fixed height|read from file : read from file'
           write(0,*) 'Number of air layers and dz top to bottom km  : 10 500. 200. 100. 50. 20. 10. 5. 2. 1. 0.5'
           write(0,*)
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	       ctrl%wFile_Data = temp(3)
	    end if
	    if (narg > 3) then
	       ctrl%wFile_EMsoln = temp(4)
	    end if
	    if (narg > 4) then
	       ctrl%rFile_fwdCtrl = temp(5)
	    end if
        if (narg > 5) then
           ctrl%rFile_EMrhs = temp(6)
        end if

      case (COMPUTE_J) ! J
        if (narg < 3) then
           write(0,*) 'Usage: -J  rFile_Model rFile_Data wFile_Sens [rFile_fwdCtrl]'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	       ctrl%wFile_Sens = temp(3)
	    end if
	    if (narg > 3) then
	       ctrl%rFile_fwdCtrl = temp(4)
	    end if

      case (MULT_BY_J) ! M
        if (narg < 4) then
           write(0,*) 'Usage: -M  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl]'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_dModel = temp(2)
	       ctrl%rFile_Data = temp(3)
	       ctrl%wFile_Data = temp(4)
	    end if
	    if (narg > 4) then
	       ctrl%rFile_fwdCtrl = temp(5)
	    end if

      case (MULT_BY_J_T) ! T
        if (narg < 3) then
           write(0,*) 'Usage: -T  rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	       ctrl%wFile_dModel = temp(3)
	    end if
	    if (narg > 3) then
	       ctrl%rFile_fwdCtrl = temp(4)
	    end if

      case (MULT_BY_J_T_multi_Tx) ! x
        if (narg < 3) then
           write(0,*) 'Usage: -x  rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	       ctrl%wFile_dModel = temp(3)
	    end if
	    if (narg > 3) then
	       ctrl%rFile_fwdCtrl = temp(4)
	    end if

      case (INVERSE) ! I
        if (narg < 3) then
           write(0,*) 'Usage: -I NLCG rFile_Model rFile_Data [lambda eps]'
           write(0,*)
           write(0,*) 'Here, lambda = the initial damping parameter for inversion'
           write(0,*) '         eps = misfit tolerance for the forward solver'
           write(0,*) 'OR'
           write(0,*) 'Usage: -I NLCG rFile_Model rFile_Data [rFile_invCtrl rFile_fwdCtrl]'
           write(0,*)
           write(0,*) 'Here, rFile_invCtrl = the inversion control file in the format'
           write(0,*)
           write(0,*) 'Model and data output file name    : Example'
           write(0,*) 'Initial damping factor lambda      : 1.'
           write(0,*) 'To update lambda divide by         : 10.'
           write(0,*) 'Initial search step in model units : 100.'
           write(0,*) 'Restart when rms diff is less than : 2.0e-3'
           write(0,*) 'Exit search when rms is less than  : 1.05'
           write(0,*) 'Exit when lambda is less than      : 1.0e-4'
           write(0,*) 'Maximum number of iterations       : 120'
           write(0,*)
           write(0,*) '      rFile_fwdCtrl = the forward solver control file in the format'
           write(0,*)
           write(0,*) 'Number of QMR iters per divergence correction : 40'
           write(0,*) 'Maximum number of divergence correction calls : 20'
           write(0,*) 'Maximum number of divergence correction iters : 100'
           write(0,*) 'Misfit tolerance for EM forward solver        : 1.0e-7'
           write(0,*) 'Misfit tolerance for EM adjoint solver        : 1.0e-7'
           write(0,*) 'Misfit tolerance for divergence correction    : 1.0e-5'
           write(0,*) 'Optional EM solution file name for nested BC  : nested.esoln'
           write(0,*) 'Air layers mirror|fixed height|read from file : fixed height'
           write(0,*) 'Number of air layers and max height in km     : 12 1000'
           write(0,*)
           write(0,*) 'Optionally, may also supply'
           write(0,*)
           write(0,*) '      rFile_Cov     = the model covariance configuration file'
           write(0,*) '      rFile_dModel  = the starting model parameter perturbation'
           write(0,*)
           stop
        else
           ctrl%search = temp(1)
           select case (ctrl%search)
           case ('NLCG','DCG','Hybrid')
              	! write(0,*) 'Inverse search ',trim(ctrl%search),' selected.'
           case default
				write(0,*) 'Unknown inverse search. Usage: -I [NLCG | DCG | Hybrid]'
				stop
           end select
	       ctrl%rFile_Model = temp(2)
	       ctrl%rFile_Data = temp(3)
	    end if
	    if (narg > 3) then
          read(temp(4),*,iostat=istat) ctrl%lambda
          res=is_letter(temp(4))
          if (res) then
            ! check for the inverse solver configuration file
            ctrl%rFile_invCtrl = temp(4)
            inquire(FILE=ctrl%rFile_invCtrl,EXIST=exists)
            if (.not. exists) then
            	! problem - invalid argument
				write(0,*) 'Please specify a valid inverse control file or damping parameter'
				stop
            end if
          end if
        end if
        if (narg > 4) then
          read(temp(5),*,iostat=istat) ctrl%eps
          res=is_letter(temp(5))
          if (res) then
            ! check for the forward solver configuration file
            ctrl%rFile_fwdCtrl = temp(5)
            inquire(FILE=ctrl%rFile_fwdCtrl,EXIST=exists)
            if (.not. exists) then
            	! problem - invalid argument
				write(0,*) 'Please specify a valid forward solver control file or misfit tolerance'
				stop
            end if
          end if
        end if
	    if (narg > 5) then
	       ctrl%rFile_Cov = temp(6)
	    end if
        if (narg > 6) then
            ! check for the optional starting model file
            ctrl%rFile_dModel = temp(7)
            inquire(FILE=ctrl%rFile_dModel,EXIST=exists)
            if (.not. exists) then
            	! problem - invalid argument
				write(0,*) 'Please specify a valid starting model file'
				stop
            end if
        end if

      case (APPLY_COV) ! C
        if (narg < 3) then
           write(0,*) 'Usage: -C  [FWD|INV] rFile_Model wFile_Model [rFile_Cov rFile_Prior]'
           write(0,*)
           write(0,*) ' -C FWD rFile_Model wFile_Model [rFile_Cov rFile_Prior]'
           write(0,*) '  Applies the model covariance to produce a smooth model output'
           write(0,*) '  Optionally, also specify the prior model to compute resistivities'
           write(0,*) '  from model perturbation: m = C_m^{1/2} \\tilde{m} + m_0'
           write(0,*)
           write(0,*) ' -C INV rFile_Model wFile_Model [rFile_Cov rFile_Prior]'
           write(0,*) '  Applies the inverse of model covariance (if implemented)'
           write(0,*) '  Optionally, also specify the prior model to compute starting'
           write(0,*) '  perturbation for the inversion: \\tilde{m} = C_m^{-1/2} (m - m_0)'
           write(0,*) '  WARNING: may be poorly conditioned if the smoothing is too strong;'
           write(0,*) '  always check your model for white noise after using this option.'
           stop
        else
        ctrl%option = temp(1)
	    ctrl%rFile_Model = temp(2)
	    ctrl%wFile_Model = temp(3)
        end if
        if (narg > 3) then
            ctrl%rFile_Cov = temp(4)
        end if
        if (narg > 4) then
            ctrl%rFile_Prior = temp(5)
        end if

      case (EXTRACT_BC) ! b
        if (narg < 3) then
           write(0,*) 'Usage: -b  rFile_Model rFile_Data wFile_EMrhs [rFile_fwdCtrl]'
           write(0,*)
           write(0,*) '  Initializes the forward solver and extracts the boundary conditions,'
           write(0,*) '  writes to file.'
           stop
        else
        ctrl%rFile_Model = temp(1)
        ctrl%rFile_Data = temp(2)
        ctrl%wFile_EMrhs = temp(3)
        if (narg > 3) then
            ctrl%rFile_fwdCtrl = temp(4)
        end if
        end if

      case (TEST_GRAD) !g
        if (narg < 3) then
           write(0,*) 'Usage: -g  rFile_Model rFile_Data rFile_dModel [rFile_fwdCtrl rFile_EMrhs]'
           write(0,*)
           write(0,*) '  The ultimate test of the gradient computations. Based on the Taylor'
           write(0,*) '  series approximation:'
           write(0,*) '  Compute f(m0+dm) - f(m0)'
           write(0,*) '  Compute df/dm|_m0 x dm as a dot product'
           write(0,*) '  Compare the two resultant scalars'
           write(0,*)
           stop
        else
           ctrl%rFile_Model = temp(1)
           ctrl%rFile_Data = temp(2)
           ctrl%rFile_dModel = temp(3)
        end if
        if (narg > 3) then
           ctrl%rFile_fwdCtrl = temp(4)
        end if
        if (narg > 4) then
           ctrl%rFile_EMrhs = temp(5)
        end if

      case (TEST_ADJ) ! A
        if (narg < 3) then
           write(0,*) 'Usage: Test the adjoint implementation for each of the critical'
           write(0,*) '       operators in J = L S^{-1} P + Q'
           write(0,*)
           write(0,*) ' -A  J rFile_Model rFile_dModel rFile_Data [wFile_Model wFile_Data]'
           write(0,*) '  Tests the equality d^T J m = m^T J^T d for any model and data.'
           write(0,*) '  Optionally, outputs J m and J^T d.'
           write(0,*)
           write(0,*) ' -A  L rFile_Model rFile_EMsoln rFile_Data [wFile_EMrhs wFile_Data]'
           write(0,*) '  Tests the equality d^T L e = e^T L^T d for any EMsoln and data.'
           write(0,*) '  Optionally, outputs L e and L^T d.'
           write(0,*)
           write(0,*) ' -A  S rFile_Model rFile_EMrhs rFile_Data [wFile_EMsoln]'
           write(0,*) '  Tests the equality b^T S^{-1} b = b^T (S^{-1})^T b for any EMrhs.'
           write(0,*) '  For simplicity, use one EMrhs for forward and transpose solvers.'
           write(0,*) '  Data file only needed to set up dictionaries.'
           write(0,*) '  Optionally, outputs e = S^{-1} b.'
           write(0,*)
           write(0,*) ' -A  P rFile_Model rFile_dModel rFile_EMsoln rFile_Data [wFile_Model wFile_EMrhs]'
           write(0,*) '  Tests the equality e^T P m = m^T P^T e for any EMsoln and data.'
           write(0,*) '  The data template isn''t needed here except to set up the transmitters.'
           write(0,*) '  Optionally, outputs P m and P^T e.'
           write(0,*)
           write(0,*) ' -A  Q rFile_Model rFile_dModel rFile_Data [wFile_Model wFile_Data]'
           write(0,*) '  Tests the equality d^T Q m = m^T Q^T d for any model and data.'
           write(0,*) '  Optionally, outputs Q m and Q^T d.'
           write(0,*)
           write(0,*) 'Finally, generates random 5% perturbations, if implemented:'
           write(0,*) ' -A  m rFile_Model wFile_Model [delta]'
           write(0,*) ' -A  d rFile_Data wFile_Data [delta]'
           write(0,*) ' -A  e rFile_Model rFile_Data rFile_EMsoln wFile_EMsoln [delta]'
           write(0,*) ' -A  b rFile_Model rFile_Data rFile_EMrhs wFile_EMrhs [delta]'
           stop
        else
           ctrl%option = temp(1)
           ctrl%delta = 0.05
           select case (ctrl%option)
           ! tests of adjoint implementation ...
           case ('J')
                ctrl%rFile_Model = temp(2)
                ctrl%rFile_dModel = temp(3)
                ctrl%rFile_Data = temp(4)
                if (narg > 4) then
                    ctrl%wFile_Model = temp(5)
                endif
                if (narg > 5) then
                    ctrl%wFile_Data = temp(6)
                endif
           case ('L')
                ctrl%rFile_Model = temp(2)
                ctrl%rFile_EMsoln = temp(3)
                ctrl%rFile_Data = temp(4)
                if (narg > 4) then
                    ctrl%wFile_EMrhs = temp(5)
                endif
                if (narg > 5) then
                    ctrl%wFile_Data = temp(6)
                endif
           case ('S')
                ctrl%rFile_Model = temp(2)
                ctrl%rFile_EMrhs = temp(3)
                ctrl%rFile_Data = temp(4)
                if (narg > 4) then
                    ctrl%wFile_EMsoln = temp(5)
                endif
           case ('P')
                ctrl%rFile_Model = temp(2)
                ctrl%rFile_dModel = temp(3)
                ctrl%rFile_EMsoln = temp(4)
                if (narg < 5) then
                    write(0,*) 'Usage: -P rFile_Model rFile_mgCtrl rFile_dModel rFile_EMsoln rFile_Data [wFile_Model wFile_EMrhs]'
                    write(0,*) 'Please specify data template file to set up the transmitter dictionary'
                    stop
                endif
                ctrl%rFile_Data = temp(5)
                if (narg > 5) then
                    ctrl%wFile_Model = temp(6)
                endif
                if (narg > 6) then
                    ctrl%wFile_EMrhs = temp(7)
                endif
           case ('Q')
                ctrl%rFile_Model = temp(2)
                ctrl%rFile_dModel = temp(3)
                ctrl%rFile_Data = temp(4)
                if (narg > 4) then
                    ctrl%wFile_Model = temp(5)
                endif
                if (narg > 5) then
                    ctrl%wFile_Data = temp(6)
                endif
           ! random perturbations ... in principle, shouldn't need model and data
           ! to create random solution and RHS. But using these to create dictionaries.
           ! This is an artifact of reading routines and can later be fixed.
           case ('m')
                ctrl%rFile_Model = temp(2)
                ctrl%wFile_Model = temp(3)
                if (narg > 3) then
                    read(temp(4),*,iostat=istat) ctrl%delta
                endif
           case ('d')
                ctrl%rFile_Data = temp(2)
                ctrl%wFile_Data = temp(3)
                if (narg > 3) then
                    read(temp(4),*,iostat=istat) ctrl%delta
                endif
           case ('e')
                ctrl%rFile_Model = temp(2)
                ctrl%rFile_Data = temp(3)
                ctrl%rFile_EMsoln = temp(4)
                ctrl%wFile_EMsoln = temp(5)
                if (narg > 5) then
                    read(temp(6),*,iostat=istat) ctrl%delta
                endif
           case ('b')
                ctrl%rFile_Model = temp(2)
                ctrl%rFile_Data = temp(3)
                ctrl%rFile_EMrhs = temp(4)
                ctrl%wFile_EMrhs = temp(5)
                if (narg > 5) then
                    read(temp(6),*,iostat=istat) ctrl%delta
                endif
           case default
                write(0,*) 'Unknown operator. Usage: -A [J | L | S | P | Q] OR -A [m | d | e | b]'
                stop
           end select
        end if

      case (TEST_SENS) ! S
        if (narg < 4) then
           write(0,*) 'Usage: -S  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl]'
           stop
        else
           ctrl%rFile_Model = temp(1)
           ctrl%rFile_dModel = temp(2)
           ctrl%rFile_Data = temp(3)
           ctrl%wFile_Data = temp(4)
        end if
        if (narg > 4) then
           ctrl%rFile_fwdCtrl = temp(5)
        end if

      case default
         write(0,*) 'Unknown job. Please check your command line options'
         stop

     end select


     deallocate(temp)

     ! save this info for the main program
     ctrl%job = job
  end subroutine parseArgs
end module UserCtrl
