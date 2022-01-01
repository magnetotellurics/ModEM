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
   !
   !use ForwardSolverFromFile
   !use ForwardSolverIT_DC
   !
   !use SourceMT_1D
   !use SourceMT_2D
   !
   !    THIS IS BASED ON test_Amult -- extended to also test source, solvers ...
   ! 
   class(Grid_t ), allocatable, target           :: main_grid
   class( ModelParameter_t ), allocatable :: model_parameter
   class( ModelOperator_t ), allocatable  :: model_operator
   class( Solver_t ), allocatable  :: slvrQMR,slvrPCG
   class( Source_t ), allocatable  :: src
   !   other things I make explicit types  -- seemed to work, but now not sure!
   class( CVector_t), allocatable   :: x, y
   class( CScalar_t), allocatable   :: phiIn, phiOut
   class( RVector_t), allocatable   :: xR, yR
   class( RScalar_t), allocatable   :: phiInR, phiOutR
   type( TAirLayers)   :: air_layer
   !
   character(:), allocatable :: control_file_name, model_file_name, data_file_name, modem_job
   character(:), allocatable :: xFile,yFile,gridType
   integer  :: printUnit
   real(kind = prec) :: omega
   !
   !   frequency is hard coded -- just test for a single frequency at a time
   omega = 2*pi/.1
   !
   !   test job is also hard coded : options- Amult, QMR, RHS
   modem_job = 'PCG'    
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
   write ( *, * ) "Finish test Solver"
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
      character(12) :: method = 'fixed height'
      integer :: nzAir = 2
      real(kind=prec) :: maxHeight = 1.5  !   this should be in km, not meters
      !
      !      fname = "/mnt/c/Users/protew/Desktop/ON/GITLAB_PROJECTS/modem-oo/inputs/Full_A_Matrix_TinyModel"
      fnameA = "/Users/garyegbert/Desktop/ModEM_ON/modem-oo/inputs/Full_A_Matrix_TinyModel"
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
            !   create RVectors
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
            !
            !   complete model operator setup
               call model_operator%SetEquations()
               call model_parameter%setMetric( model_operator%metric )
               call model_operator%SetCond( model_parameter )
            !
            !   Instantiate the Solver objects
            ! 
            slvrQMR = Solver_QMR_t(model_operator)
            slvrPCG = Solver_PCG_t(model_operator)
            !
            !   Instantiate the Source object
            ! 
            src = SourceMT_1D_t(model_operator,model_parameter)

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
      write(*,*) 'JOB = ', modem_job
      select case(modem_job)
         case("Amult")
            !   TEST OF CURL-CURL OPERATOR
            !   read in cVector x used for test of A*x = y
            !   these are hard-coded at present
            xFile = '../inputs/Xvec_Tiny_1.dat'
            call readCVector()
            !    multiply by A
            call model_operator%Amult(omega,x,y)
            !
            !   set output file name for this test
            !      ... different for different model_operator types
            !           so results can be compared ...
            select type(model_operator)
               class is(ModelOperator_MF_t) 
               yFile = '../inputs/Yvec_Tiny_MF_1.dat'
            class is(ModelOperator_File_t)
               yFile = '../inputs/Yvec_Tiny_File_1.dat'
            end select
            call writeCVector()
         case("QMR")
            !   TEST OF QMR SOLVER
            !   read in cVector used for test -- rhs in A*y = x           
            xFile = '../inputs/RHS_Tiny.dat'
            call readCVector()
            !  create and setup Solver object ...
            call slvrQMR%SetDefaults()   !   set default convergence parameters
            !   first test w/o preconditioner
            slvrQMR%omega = omega
            slvrQMR%preconditioner = PreConditioner_None_t()
            select type(slvrQMR)
               class is(Solver_QMR_t)
                  call slvrQMR%solve(x,y)
                  write(*,*) 'n_iter',slvrQMR%n_iter
                  write(*,*) 'relative residual',slvrQMR%relErr(slvrQMR%n_iter)
                  write(57,*) slvrQMR%relErr(1:slvrQMR%n_iter)

                  !   file name for output
                  yFile = '../inputs/Soln_Tiny_QMR.dat'
               class default
                 stop "test program not coded for this solver type"
            end select
            call writeCVector()
         case("RHS")
           !   TEST OF RHS COMPUTATIONS -- does not test 1D modeling!
            !   read in array E from file, create and output rhs cVector
            xFile = '../inputs/E_Tiny.dat'
            call readCVector()
            !   put cVector E into source object
            allocate(src%E, source = x)
            write(*,*)  'src%E'
            !call src%E%print()
            !    compute RHS and output
            call src%SetRHS()
            y = src%bdry
            yFile = '../inputs/BDRYcompTiny.dat'
            call writeCVector()
            y = src%rhs
            yFile = '../inputs/RHScompTiny.dat'
            call writeCVector()
            !    compute E0 and output
            call src%setE0
            y = src%E0
            yFile = '../inputs/E0compTiny.dat'
            call writeCVector()
          case("MultDC")
            !    write out Sigma_E after setting ... 
            select type(model_operator)
               class is( ModelOperator_MF_t)
                  yR = model_operator%Sigma_E
                  yFile = '../inputs/Sigma_E.dat'
                  call writeRVector()
            end select
            !   TEST OF ModOp.divCgrad * phi_in
            !   read in cVector x used for test of A*x = y
            !   these are hard-coded at present
            xFile = '../inputs/PhiIn_Tiny.dat'
            call readCScalar()
            !    multiply by divCgrad
            call model_operator%divCgrad(phiIn,phiOut)
            !
            !   set output file name for this test
            !      ... different for different model_operator types
            !           so results can be compared ...
            yFile = '../inputs/PhiOut_Tiny.dat'
            call writeCScalar()
         case("PCG")
            !   TEST OF PCG SOLVER
            !   read in cScalar used for test -- rhs in A*y = x           
            xFile = '../inputs/PhiIn_Tiny.dat'
            call readCScalar()
            !  create and setup Solver object ...
            call slvrPCG%SetDefaults()   !   set default convergence parameters
            !   first test w/o preconditioner
            ! slvrPCG%omega = omega   !   don't need to set omega in this case
            slvrPCG%preconditioner = PreConditioner_None_t()
            select type(slvrPCG)
               class is(Solver_PCG_t)
                  write(*,*) 'class is Solver_PCG_t'

                  call slvrPCG%solve(phiIn,phiOut)
                  write(*,*) 'n_iter',slvrPCG%n_iter
                  write(*,*) 'relative residual',slvrPCG%relErr(slvrPCG%n_iter)
                  write(57,*) slvrPCG%relErr(1:slvrPCG%n_iter)

                  !   file name for output
                  yFile = '../inputs/Soln_Tiny_PCG.dat'
               class default
                 stop "test program not coded for this solver type"
            end select
            call writeCScalar()
          end select

     end subroutine runTest
     !
     ! ********
     !
     subroutine readCVector()
        integer :: fid

        !   open input file, read in x
        fid = 55
        open(file = xFile,unit = fid, form='unformatted')
        select type(x)
           class is (cVector3D_SG_t)
              call x%read(fid)
           class default
              write(*,*)  'CVector of incorrect class'
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
        open(file = yFile,unit = fid, form='unformatted')
        select type(y)
            class is (CVector3D_SG_t)
               call y%write(fid)
            class default
               write(*,*)  'CVector of incorrect class'
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
        open(file = xFile,unit = fid, form='unformatted')
        select type(phiIn)
           class is (cScalar3D_SG_t)
              call phiIn%read(fid,'b')   !   read for CScalar requires
                                         !  specification of 'b' for binary
              write(*,*) 'phiIn read'
           class default
              write(*,*)  'CScalar of incorrect class'
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
        open(file = yFile,unit = fid, form='unformatted')
        select type(phiOut)
            class is (CScalar3D_SG_t)
               call phiOut%write(fid,'b')
            class default
               write(*,*)  'CScalar of incorrect class'
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
        open(file = xFile,unit = fid, form='unformatted')
        select type(xR)
           class is (rVector3D_SG_t)
              call xR%read(fid)
           class default
              write(*,*)  'CVector of incorrect class'
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
        open(file = yFile,unit = fid, form='unformatted')
        select type(yR)
            class is (RVector3D_SG_t)
               call yR%write(fid)
            class default
               write(*,*)  'CVector of incorrect class'
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
        open(file = xFile,unit = fid, form='unformatted')
        select type(phiInR)
           class is (rScalar3D_SG_t)
              call phiInR%read(fid,'b')   !   read for CScalar requires
                                         !  specification of 'b' for binary
              write(*,*) 'phiIn read'
           class default
              write(*,*)  'CScalar of incorrect class'
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
        open(file = yFile,unit = fid, form='unformatted')
        select type(phiOutR)
            class is (RScalar3D_SG_t)
               call phiOutR%write(fid,'b')
            class default
               write(*,*)  'CScalar of incorrect class'
              stop
       end select
       close(fid)
    end subroutine writeRScalar

end program test_Solver
