! *****************************************************************************
! mtInvSetup.f90 is derived from mtFwdSetup.f90
!   quick and dirty initialization, incorporating Kush's xml startup,
!   but getting rid of all of various reeadModelInput versions
!   ONLY supports Mackie format, non-xml startup
!          
module mtinvsetup

  use math_constants		! math/ physics constants
  use utilities
  use emsolve3d
  use ModelSpace
  use ioascii
  

  implicit none

  ! Initialization routines.
  public                :: Startup3D

	! main startup driver ...  this uses: 
	!   (some might need to be public)
  private               :: ReadStartFile,ReadSolverControl


  ! ***************************************************************************
     ! file names which may be read from startup ... 
     !                   only some are used in this version

     character (len=20),private	 :: version = ''
     character (len=80),private  :: fn_model = '', fn_air = '', fn_cper = ''
     character (len=80),private  :: fn_bc = '', fn_E0 = '', fn_diagn = ''
     character (len=80),private  :: fn_solverc = '', fn_soln = ''
     character (len=80),private  :: fn_z = '', fn_nested ='', fn_nbc = ''

     ! i/ o units
     integer,private			:: ioBC, ioEsol  
     integer,private			:: ioE0, ioDiag
     integer,private			:: ioZ, ioNest, ioNBC
     
     ! variables used for reading xml file
     logical, private  :: in_header = .false., in_earth = .false.
     logical, private  :: in_air = .false., in_freq = .false., in_bound = .false.
     logical, private  :: in_start = .false., in_diag = .false.
     logical, private  :: in_solver = .false., in_elect = .false.
     logical, private  :: in_imped = .false., in_nearth = .false.
     logical, private  :: in_nbc = .false.

Contains

  subroutine StartUp3D(StartFile,paramType,grid,solverParams,Sigma)

  !  this reads the non-xml startup file to get file name for Mackie
  !   style input (plus solver control file name; could also be used
  !   to provide other input or output file names, as it did for the
  !   forward mt code) . paramType is an input argument used to
  !   control type of returned modelParam object Sigma (log vs. linear
  !   conductivity)

     character(len=80),intent(in)	:: StartFile
     character(len=80),intent(in)	:: paramType
     type(grid_t), intent(out)	:: grid
     type(emsolve_control), intent(out)	:: solverParams		
     type(modelParam_t), intent(out)	:: Sigma
     !  reinstate use of this to pass output file names back to main
     ! type(output_files), intent(out)	:: outFiles
     
     !   local variables
     integer 				:: i1, i2, fidRM
     type(rscalar)			:: Cond

     ! I/O unit numbers ... most of these not now used
     ioBC = 11
     ioEsol = 21
     ioE0 = 31
     ioDiag = 41
     ioZ = 51
     ioNest = 61
     ioNBC = 71

     ! Get filenames from the .xml startup file
     i2 = len_trim(StartFile)
     i1 = i2 -3
     if (StartFile(i1:i2) == '.xml') then
        ! a .xml file is being used for StartFile
        call errstop('The XML startup format is not supported yet')
     else if (StartFile(i1:i2) /= '.xml') then
        Call ReadStartFile(StartFile)
     end if

     if (fn_model == '') then
        call errstop('Grid/conductivity model file not specified: stopping')
     end if
   
     ! Read input files and set up basic grid geometry, conductivities,
     ! and frequencies (stored in the transmitter dictionary, txDictMT)
     call ReadRMgridCond(fidRM,fn_model,grid,Cond)

     ! move cell conductivities read into rscalar object into a modelParam
     ! object ... this dance needed to keep modelParam attributes private
     !   First need to create model parameter
     call create_modelParam(grid,paramType,Sigma)
     call set_modelParam(Sigma,Cond,LINEAR)

     ! now done with Cond, so deallocate
     call deall_rscalar(Cond)

     if (fn_solverc /= '') then
        Call ReadSolverControl(solverParams)
     else
        solverParams%UseDefaults= .true.
     endif

  end subroutine StartUp3D
!********************************************************************************

  ! * ReadStartFile reads the file from the screen or hardwired that contains
 ! * information about the files that are need to run the 3D model for
  ! * electromagnetic induction. Note: This is for non-XML file in ascii format.
  ! * All the memroy allocation is done inside.
  subroutine ReadStartFile(cfileStartUp)

    implicit none
    character (len=80), intent(in)              :: cfileStartUp
    integer                                     :: ios

    OPEN (UNIT=10,FILE=cfileStartUp,STATUS='OLD',IOSTAT=ios)

    IF( ios/=0) THEN
       write(0,*) 'Error opening file:', cfileStartUp
    ENDIF

    ! This is the list of files that will be used for i/ o
    READ(10,'(a80)') fn_model
    READ(10,'(a80)') fn_air
    READ(10,'(a80)') fn_cper 
    READ(10,'(a80)') fn_E0
    READ(10,'(a80)') fn_diagn
    READ(10,'(a80)') fn_solverc
    READ(10,'(a80)') fn_soln
    READ(10,'(a80)') fn_z
    READ(10,'(a80)') fn_bc
!    READ(10,'(a80)') fn_nested
!    READ(10,'(a80)') fn_nbc

    CLOSE(10)

  end subroutine ReadStartFile  ! ReadStartFile

!*****************************************************************************
subroutine ReadSolverControl(solverParams)
  ! Purpose:  read solver control parameters from file fn_solverc

  implicit none
  type(emsolve_control), intent(out)	:: solverParams
  integer                               :: ios

  open (unit=10,file=fn_solverc, status='old',iostat=ios)

  if( ios/=0) then
      solverParams%UseDefaults = .true.
      write(0,*) 'Using default solver controls'
  else
      solverParams%UseDefaults = .false.
  endif

  ! number of QMR iterations before a call for divergence correction
  Read(10, *) solverParams%IterPerDivCor
  ! Maximum number of divergence correction calls
  Read(10,*) solverParams%MaxDivCor
  ! Maximum number of number of PCG iterations in each divergence correction
  Read(10,*) solverParams%MaxIterDivCor

  ! Error tolerance for divergence correction and EMsolve
  Read(10, *) solverParams%tolEM
  Read(10, *) solverParams%tolDivCor
  close(10)
  return
end subroutine ReadSolverControl

end module mtinvsetup
