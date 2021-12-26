program test_MultA
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
   !
   !use ForwardSolverFromFile
   !use ForwardSolverIT_DC
   !
   !use SourceMT_1D
   !use SourceMT_2D
   !
   ! 
   class(Grid_t ), allocatable, target           :: main_grid
   class( ModelParameter_t ), allocatable :: model_parameter
   class( ModelOperator_t ), allocatable  :: model_operator
   !   other things I make explicit types
   class( CVector3D_SG_t), allocatable   :: x, y
   type( TAirLayers)   :: air_layer
   !
   character(:), allocatable :: control_file_name, model_file_name, data_file_name, modem_job
   character(:), allocatable :: xFile,yFile,gridType
   integer  :: fid
   real(kind = prec) :: omega
   !
   write ( *, * )
   write ( *, * ) "Start Amult test"
   write ( *, * )
   !
   call handleModelFile()    !   this reads model file, sets up model_operator
   !   read in cVector x (one for now)
   xFile = '../inputs/Xvec_Tiny_1.dat'
   select type(model_operator)
       class is(ModelOperator_MF_t) 
          yFile = '../inputs/Yvec_Tiny_MF_1.dat'
       class is(ModelOperator_File_t)
          yFile = '../inputs/Yvec_Tiny_File_1.dat'
   end select

   !   create CVectors in handleModelFile?
   !x = cVector3D_SG_t(main_grid,gridType)
   !y = cVector3D_SG_t(main_grid,gridType)
   !   open input file, read in x
   open(file = xFile,unit = fid, form='unformatted')
   select type(x)
      class is (cVector3D_SG_t)
         call x%read(fid)
      class default
         write(*,*)  'CVector of incorrect class'
         stop
   end select

   close(fid)

   !    multiply by A
   omega = R_ONE
   call model_operator%Amult(omega,x,y)
   !  open output file, write out y 
   open(file = yFile,unit = fid, form='unformatted')
   select type(y)
      class is (CVector3D_SG_t)
         call y%write(fid)
      class default
         write(*,*)  'CVector of incorrect class'
         stop
   end select
   close(fid)
   !   also print result to ascii file
   call y%print(667,'y = A*x')

   write ( *, * )
   write ( *, * ) "Finish ModEM-OO."
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
            ! Instantiate the ModelOperator object
            ! 
            !model_operator = ModelOperator_File_t( main_grid, fnameA )
            model_operator = ModelOperator_MF_t( main_grid )
            !
            call model_operator%SetEquations()
            !
            call model_parameter%setMetric( model_operator%metric )
            !
            call model_operator%SetCond( model_parameter )
            !
         class default
             stop "Unclassified main_grid"
         !
      end select
      !
   end subroutine handleModelFile
   !
   !
end program test_MultA
