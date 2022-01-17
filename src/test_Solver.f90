 program test_Solver
   !
   !use DataManager
   !use ModEMControlFile
   !
   use Constants
   use ModelReader
   use ModelReader_Weerachai
   use ModelOperator_MF
   use ModelOperator_File
   use ModelParameterCell_SG
   !
   use Grid3D_SG
   use CVector3D_SG
   use CScalar3D_SG
   !
   use Solver_QMR
   use Solver_PCG
   use SourceMT_1D
   use ForwardSolverIT
   use DivergenceCorrection
   use ForwardSolverIT_DC
   !
   !use ForwardSolverIT_DC
   !
   !use SourceMT_1D
   !use SourceMT_2D
   !
   !    THIS IS BASED ON test_Amult -- extended to also test source, solvers ...
   ! 
   class( Grid_t ), allocatable, target   :: main_grid
   class( ModelParameter_t ), allocatable :: model_parameter
   class( ModelOperator_t ), allocatable  :: model_operator
   class( Solver_t ), allocatable  :: slvrQMR,slvrPCG
   class( Source_t ), allocatable  :: src
   class( ForwardSolver_t ), allocatable  :: fwdIT, fwdIT_DC
   type( DivergenceCorrection_t)  :: divCor
   !   other things I make explicit types  -- seemed to work, but now not sure!
   class( CVector_t), allocatable   :: x, y
   class( CScalar_t), allocatable   :: phiIn, phiOut
   class( RVector_t), allocatable   :: xR, yR
   class( RScalar_t), allocatable   :: phiInR, phiOutR
   type( TAirLayers )   :: air_layer
   !
   character(:), allocatable :: control_file_name, model_file_name, data_file_name, modem_job
   character(:), allocatable :: xFile,yFile,gridType
   integer  :: printUnit, maxIter
   real(kind = prec) :: omega, T, tolerance
   !
   !   frequency is hard coded -- just test for a single frequency at a time
   T = .1
   omega = 2*pi/T
   !
   !   test job is also hard coded : options- Amult, QMR, RHS, MULT_DC, 
   !            LUsolve, PCG, FWD_IT, DC, FWD_IT_DC, CurlT
   modem_job = "CurlT"    
   fid = 1
   printUnit = 667   !   change this to get output y vector in a different ascii file
   !
   write ( *, * )
   write ( *, * ) "Start Solver test program"
   write ( *, * )
   !
   call handleModelFile()    !   this reads model file, sets up model_operator
   !    also creates x and y 

   call runTest()

   write ( *, * )
   write ( *, * ) "Finish test Solver for job ",modem_job
   write ( *, * )
   !
contains
   !
   subroutine handleModelFile()
      implicit none
      !
      ! It remains to standardize ????
      type( ModelReader_Weerachai_t ) :: model_reader
      !
      character(:), allocatable :: fnameA
      !    parameters for setting Air Layers for Tiny Model
      character(12) :: method = "fixed height"
      integer :: nzAir = 2
      real(kind=prec) :: maxHeight = 1.5  !   this should be in km, not meters
      !
      !fnameA = "/mnt/c/Users/protew/Desktop/ON/GITLAB_PROJECTS/modem-oo/inputs/Full_A_Matrix_TinyModel"
      fnameA = "/Users/garyegbert/Desktop/ModEM_ON/modem-oo/inputs/Full_A_Matrix_TinyModel"
      !model_file_name = "/mnt/c/Users/protew/Desktop/ON/GITLAB_PROJECTS/modem-oo/inputs/rFile_Model_Tiny"
      model_file_name = "/Users/garyegbert/Desktop/ModEM_ON/modem-oo/inputs/rFile_Model_Tiny"
      !
      write( *, * ) "   -> Model File: [", model_file_name, "]"
      !
      ! Read Grid and ModelParameter with ModelReader_Weerachai
      call model_reader%Read( model_file_name, main_grid, model_parameter ) 
      !
      select type( main_grid )
         !
         class is( Grid3D_SG_t )
            !   ADD Air Layers to grid -- this has to be done in general now!!!
            !    I am just setting this explicitly for the A from file case --
            !    More generally, could have a default ---  but I think we should
            !    set this explicitly outside of the gridReader procedure(s)

            call main_grid%SetupAirLayers(air_layer, method, nzAir,maxHeight)
            !   as coded have to use air_layer data structure to update grid
            call main_grid%UpdateAirLayers(air_layer%nz, air_layer%dz)

            !   create CVectors
            gridType = EDGE
            allocate(x, source = cVector3D_SG_t(main_grid,gridType))
            allocate(y, source = cVector3D_SG_t(main_grid,gridType))
            !   create RVectors (EDGE again)
            gridType = EDGE
            allocate(xR, source = rVector3D_SG_t(main_grid,gridType))
            allocate(yR, source = rVector3D_SG_t(main_grid,gridType))
            !   create CScalars
            gridType = NODE
            allocate(phiIn, source = cScalar3D_SG_t(main_grid,gridType))
            allocate(phiOut, source = cScalar3D_SG_t(main_grid,gridType))
            !   create RScalars
            allocate(phiInR, source = rScalar3D_SG_t(main_grid,gridType))
            allocate(phiOutR, source = rScalar3D_SG_t(main_grid,gridType))
            !
            ! Instantiate the ModelOperator object
            ! 
            !model_operator = ModelOperator_File_t( main_grid, fnameA )
            model_operator = ModelOperator_MF_t( main_grid )
            call model_parameter%setMetric( model_operator%metric )
            !
            !   complete model operator setup
            call model_operator%SetEquations()
            call model_operator%SetCond( model_parameter )
            !
            !   set pointer to metric elements in model parameter
            !    (note: this is needed because model_parameter does not depend
            !        on model_operator, which already has a pointer to metricElements)
            !
            !   Instantiate the Solver objects
            ! 
            slvrQMR = Solver_QMR_t(model_operator)
            slvrPCG = Solver_PCG_t(model_operator)
            !
            !   Instantiate the Source object
            ! 
            src = SourceMT_1D_t(model_operator,model_parameter)
            !
            !   Instantiate the forward modeling objects
            fwdIT = ForwardSolverIT_t(model_operator,QMR)
            call fwdIT%setCond(model_parameter)
            call fwdIT%setPeriod(T)
            fwdIT_DC = ForwardSolverIT_DC_t(model_operator,QMR)
            call fwdIT_DC%setCond(model_parameter)
            call fwdIT_DC%setPeriod(T)
            !
            !   Instnatiatee the DivergenceCorrection object
            !     note that this creates the PGC solver automatically
            !      creates preConditioner, and sets iteration controls to 
            !      default values
            divCor = DivergenceCorrection_t(model_operator)
            write(*,*) "divCor created"
            call divCor%setCond()   !   this has to be called AFTER setting
                                    !   conductivity in ModelOperator
         class default
             stop "Unclassified main_grid"
         !
      end select
      !
   end subroutine handleModelFile
   !
   subroutine runTest()
      implicit none
      !
      write(*,*) "JOB = ", modem_job
      select case(modem_job)
         case("Amult")
            !   TEST OF CURL-CURL OPERATOR
            !   read in cVector x used for test of A*x = y
            !   these are hard-coded at present
            xFile = "../inputs/Xvec_Tiny_1.dat"
            call readCVector()
            !    multiply by A
            call model_operator%Amult(omega,x,y)
            !
            !   set output file name for this test
            !      ... different for different model_operator types
            !           so results can be compared ...
            select type(model_operator)
               class is(ModelOperator_MF_t) 
               yFile = "../inputs/Yvec_Tiny_MF_1.dat"
            class is(ModelOperator_File_t)
               yFile = "../inputs/Yvec_Tiny_File_1.dat"
            end select
            call writeCVector()
         case("QMR")
            !   TEST OF QMR SOLVER
            !   read in cVector used for test -- rhs in A*y = x           
            xFile = "../inputs/RHS_Tiny.dat"
            call readCVector()
            !  create and setup Solver object ...
            call slvrQMR%SetDefaults()   !   set default convergence parameters
            !maxIter = 20
            !tolerance = 1d-7
            !call slvrQMR%setParameters(maxIter,tolerance)   !   set convergence parameters
            !   first test w/o preconditioner
            write(*,*) "before setting omega explicitly", slvrQMR%omega
            slvrQMR%omega = omega
            call slvrQMR%preconditioner%SetPreconditioner(omega)   !   set preconditioner
            !slvrQMR%preconditioner = PreConditioner_None_t()
            select type(slvrQMR)
               class is(Solver_QMR_t)
                  call slvrQMR%solve(x,y)
                  write(*,*) "n_iter",slvrQMR%n_iter
                  write(*,*) "relative residual",slvrQMR%relErr(slvrQMR%n_iter)
                  write(57,*) slvrQMR%relErr(1:slvrQMR%n_iter)

                  !   file name for output -- run with max_iter = 20 to make an
                  !     input for testing DC
                  if(slvrQMR%max_iter .eq. 20) then
                     yFile = "../inputs/QMR20.dat"
                  else
                     yFile = "../inputs/Soln_Tiny_QMR.dat"
                  endif
                  
               class default
                 stop "test program not coded for this solver type"
            end select
            call writeCVector()
         case("RHS")
           !   TEST OF RHS COMPUTATIONS -- does not test 1D modeling!
            !   read in array E from file, create and output rhs cVector
            xFile = "../inputs/E_Tiny.dat"
            call readCVector()
            !   cvector x into E in: source object
            allocate(src%E, source = x)
            write(*,*)  "src%E"
            !call src%E%print()
            !    compute RHS and output
            call src%SetRHS()
            y = src%bdry
            yFile = "../inputs/BDRYcompTiny.dat"
            call writeCVector()
            y = src%rhs
            yFile = "../inputs/RHScompTiny.dat"
            call writeCVector()
            !    compute E0 and output
            call src%setE0
            y = src%E0
            yFile = "../inputs/E0compTiny.dat"
            call writeCVector()
         case("MultDC")
            !    write out Sigma_E after setting ... 
            select type(model_operator)
               class is( ModelOperator_MF_t)
                  yR = model_operator%Sigma_E
                  yFile = "../inputs/Sigma_E.dat"
                  call writeRVector()
            end select
            !   TEST OF ModOp.divCgrad * phi_in
            !   read in cVector x used for test of A*x = y
            !   these are hard-coded at present
            xFile = "../inputs/PhiIn_Tiny.dat"
            call readCScalar()
            !    multiply by divCgrad
            call model_operator%divCgrad(phiIn,phiOut)
            !
            !   set output file name for this test
            !      ... different for different model_operator types
            !           so results can be compared ...
            yFile = "../inputs/PhiOut_Tiny.dat"
            call writeCScalar()
         case("LUsolve")
            xFile = "../inputs/PhiIn_Tiny.dat"
            call readCScalar()
            call slvrPCG%SetDefaults()   !   set default convergence parameters
            call slvrPCG%preconditioner%LUsolve(phiIn,phiOut)
            yFile = "../inputs/LU_Tiny_PCG.dat"
            call writeCScalar()
            
         case("PCG")
            !   TEST OF PCG SOLVER
            !   read in cScalar used for test -- rhs in A*y = x           
            xFile = "../inputs/PhiIn_Tiny.dat"
            call readCScalar()
            !  create and setup Solver object ...
            maxIter = 100
            tolerance = 1d-7
            call slvrPCG%setParameters(maxIter,tolerance)   !   set convergence parameters
                                          !   just testing -- defaults are 100, 1e-5
            call slvrPCG%preconditioner%SetPreconditioner(omega)   !   set preconditioner
                  !   NOTE: this will have to be done every time model parameter changes
            !   first test w/o preconditioner
            ! slvrPCG%omega = omega   !   don"t need to set omega in this case
            !slvrPCG%preconditioner = PreConditioner_None_t()
            select type(slvrPCG)
               class is(Solver_PCG_t)

                  call slvrPCG%solve(phiIn,phiOut)
                  write(*,*) "n_iter",slvrPCG%n_iter
                  write(*,*) "relative residual",slvrPCG%relErr(slvrPCG%n_iter+1)
                  write(57,*) slvrPCG%relErr(1:slvrPCG%n_iter+1)

                  !   file name for output
                  yFile = "../inputs/Soln_Tiny_PCG.dat"
               class default
                 stop "test program not coded for this solver type"
              end select
            call writeCScalar()
         case ("FWD_IT")
            !  finish setup of fwd object
            !call fwdIT%setPeriod( T )
            !   read in cVector used for test -- rhs in A*y = x           
            !    this file contains src%E for 1D problem
            xFile = "../inputs/E_Tiny.dat"
            call readCVector()
            !  copy cvector x into E in: source object
            allocate(src%E, source = x)
            !  set RHS
            call src%SetRHS()
            write(*,*) "RHS set up"
            y = src%rhs
            yFile = "../inputs/RHSfwdIT1.dat"
            call writeCVector()
            !
            call y%zeros   !   could try different initialization
                           !   E0 from input E_Tiny -- need to test
            write(*,*) "calling getEsolution"
            call fwdIT%getESolution( src, y )
            !
            yFile = "../inputs/Soln_Tiny_FWD_IT.dat"
            call writeCvector()
         case ("DC")
            !   read in cVector used for test -- in this case start
            !   with QMR solution after a small number of iterations (20)
            !   and then apply divergence correction
            xFile = "../inputs/QMR20.dat"
            call readCVector()
            !
            call divCor%DivCorr(x,y)
            !
            yFile = "../inputs/QMR20_DC.dat"
            call writeCVector()
         case ("FWD_IT_DC")
            !  finish setup of fwd object
            call fwdIT_DC%setPeriod( T )
            !   read in cVector used for test -- rhs in A*y = x           
            !    this file contains src%E for 1D problem
            xFile = "../inputs/E_Tiny.dat"
            call readCVector()
            !  copy cvector x into E in: source object
            allocate(src%E, source = x)
            !  set RHS
            call src%SetRHS()
            write(*,*) "RHS set up"
            y = src%rhs
            !
            maxIter = 20
            tolerance = 1d-7
            call fwdIT_DC%solver%setParameters(maxIter,tolerance)
            call y%zeros   !  try different initialization
            write(*,*) "calling getEsolution"
            call fwdIT_DC%getESolution( src, y )
            !
            ! OUTPUT THIS!
            write(*,*) "niter:  ",fwdIT_DC%n_iter_actual,   &
              "  Relative Residual",fwdIT_DC%relResFinal
            yFile = "../inputs/Soln_Tiny_FWD_IT_DC.dat"
            call writeCvector()
         case ("CurlT")
            !   TEST OF CURL transpose OPERATOR
            !   read in cVector x used for test of C'*x = y
            !     x has to be a "FACE" cVector
            xFile = "../inputs/Hvec_Tiny_1.dat"
            select type( main_grid )
               class is( Grid3D_SG_t )
                  deallocate(x)
                  gridType = FACE
                  allocate(x, source = cVector3D_SG_t(main_grid,gridType))
            end select
            call readCVector()
            !    multiply by C'
            call model_operator%MultCurlT(x,y)
            !
            !   set output file name for this test
            yFile = "../inputs/CurlTxH_Tiny1.dat"
            call writeCVector()
          end select

     end subroutine runTest
     !
     ! ********
     !
     subroutine readCVector()
        integer :: fid

        !   open input file, read in x
        fid = 55
        open(file = xFile,unit = fid, form="unformatted")
        select type(x)
           class is (cVector3D_SG_t)
              call x%read(fid)
           class default
              write(*,*)  "CVector of incorrect class"
              stop
        end select
        close(fid)
     end subroutine readCVector
     !
     ! ********
     !
     subroutine writeCVector()
        integer :: fid

        !  open output file, write out y 
        fid = 55
        open(file = yFile,unit = fid, form="unformatted")
        select type(y)
            class is (CVector3D_SG_t)
               call y%write(fid)
            class default
               write(*,*)  "CVector of incorrect class"
              stop
        end select
        close(fid)
     end subroutine writeCVector
     !
     ! ********
     !
     subroutine readCScalar()
        integer :: fid

        !   open input file, read in x
        fid = 55
        open(file = xFile,unit = fid, form="unformatted")
        select type(phiIn)
           class is (cScalar3D_SG_t)
              call phiIn%read(fid,"b")   !   read for CScalar requires
                                         !  specification of "b" for binary
              write(*,*) "phiIn read"
           class default
              write(*,*)  "CScalar of incorrect class"
              stop
        end select
        close(fid)
     end subroutine readCScalar
     !
     ! ********
     !
     subroutine writeCScalar()
        integer :: fid

        !  open output file, write out y 
        fid = 55
        open(file = yFile,unit = fid, form="unformatted")
        select type(phiOut)
            class is (CScalar3D_SG_t)
               call phiOut%write(fid,"b")
            class default
               write(*,*)  "CScalar of incorrect class"
              stop
        end select
        close(fid)
     end subroutine writeCScalar
     !
     ! ********
     !
     subroutine readRVector()
        integer :: fid

        !   open input file, read in x
        fid = 55
        open(file = xFile,unit = fid, form="unformatted")
        select type(xR)
           class is (rVector3D_SG_t)
              call xR%read(fid)
           class default
              write(*,*)  "CVector of incorrect class"
              stop
        end select
        close(fid)
     end subroutine readRVector
     !
     ! ********
     !
     subroutine writeRVector()
        integer :: fid

        !  open output file, write out y 
        fid = 55
        open(file = yFile,unit = fid, form="unformatted")
        select type(yR)
            class is (RVector3D_SG_t)
               call yR%write(fid)
            class default
               write(*,*)  "CVector of incorrect class"
              stop
        end select
        close(fid)
     end subroutine writeRVector
     !
     ! ********
     !
     subroutine readRScalar()
        integer :: fid

        !   open input file, read in x
        fid = 55
        open(file = xFile,unit = fid, form="unformatted")
        select type(phiInR)
           class is (rScalar3D_SG_t)
              call phiInR%read(fid,"b")   !   read for CScalar requires
                                         !  specification of "b" for binary
              write(*,*) "phiIn read"
           class default
              write(*,*)  "CScalar of incorrect class"
              stop
        end select
        close(fid)
     end subroutine readRScalar
     !
     ! ********
     !
     subroutine writeRScalar()
        integer :: fid

        !  open output file, write out y 
        fid = 55
        open(file = yFile,unit = fid, form="unformatted")
        select type(phiOutR)
            class is (RScalar3D_SG_t)
               call phiOutR%write(fid,"b")
            class default
               write(*,*)  "CScalar of incorrect class"
              stop
       end select
       close(fid)
    end subroutine writeRScalar

end program test_Solver
