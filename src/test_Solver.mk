# Makefile suited for building the test_Solver program
# Generated using: ./fmkmf.pl [OPTIONS] test_Solver.f90 > Makefile
# with command line options
# -p .:Classes:Classes/1D:Classes/2D:Classes/DataIO:Classes/DataSpace:Classes/ForwardSolver:Classes/Grid:Classes/ModelOperator:Classes/ModelParameter:Classes/Receivers:Classes/Solver:Classes/Source:Classes/Transmitters:Classes/Utils:Classes/Vectors
# -f90 gfortran (compiler)
# -opt -O3 -ffree-line-length-none (compiler optimisation)
# -lp /usr/lib64 (linking options: path to libraries)
# -l -llapack -lblas (linking options)
# -o objs/MT3D_OO (output directory for object files)

#  Uncomment these lines to make program for Solaris OS (legacy)
# F90 = f90
# FFLAGS = -dalign -g -C -w  -L/usr/local/lib
# LIBS = -xlic_lib=sunperf
#  Uncomment these lines to make program with g95
# include Makefile.local
# OBJDIR = ./objs/3D_MT/G95Debug
# F90 = g95
# FFLAGS = -O2
# FFLAGS = -g -ftrace=frame -fbounds-check
# MPIFLAGS = -cpp # for serial code
# MODULE = -fmod=$(OBJDIR)
# LIBS = -lblas -llapack
#  Uncomment these lines to make program with Intel compiler
# include Makefile.local
# OBJDIR = ./objs/3D_MT/IFortDebug
# F90 = ifort
# FFLAGS = -O3 -parallel -openmp #-heap-arrays
# FFLAGS = -debug all -check bounds -traceback -heap-arrays
# MPIFLAGS = -cpp # for serial code
# MODULE = -module $(OBJDIR)
# LIBS = -lblas -llapack
#  Uncomment these lines to make program with PGI compiler
# include Makefile.local
# OBJDIR = ./objs/3D_MT/PGIDebug
# F90 = pgf95  # mpif90
# FFLAGS = -O3
# FFLAGS = -g -Mprof=lines -Mbounds
# MPIFLAGS = -Mpreprocess # for serial code
# MPIFLAGS = -Bstatic  -Mipa=fast  -Mextend  -Kieee -Mpreprocess -DMPI
# MODULE = -module $(OBJDIR)
# LIBS = -llapack -lblas
# LIBS = -L/usr/lib64 -llapack -lblas -lpgftnrtl -Mprof=lines

# ------------------Macro-Defs---------------------
include Makefile.local
OBJDIR = ../../objs/MT3D_OO
F90 = gfortran 
FFLAGS = -O3 -ffree-line-length-none
MPIFLAGS =  
MODULE = -J $(OBJDIR)
LIBS_PATH = -L/usr/lib64
LIBS = -llapack -lblas

# -------------------End-macro-Defs---------------------------
OBJ = $(OBJDIR)/Constants.o $(OBJDIR)/Grid1D.o $(OBJDIR)/Grid2D.o $(OBJDIR)/Grid.o $(OBJDIR)/rScalar.o $(OBJDIR)/rVector.o $(OBJDIR)/ModelParameter1D.o $(OBJDIR)/Esoln2DTM.o $(OBJDIR)/ModelParameter2D.o $(OBJDIR)/MetricElements.o $(OBJDIR)/ModelParameter.o $(OBJDIR)/ModelReader.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/rScalar3D_SG.o $(OBJDIR)/MatUtils.o $(OBJDIR)/rVector3D_SG.o $(OBJDIR)/ModelParameterCell_SG.o $(OBJDIR)/ModelReader_Weerachai.o $(OBJDIR)/cScalar.o $(OBJDIR)/cScalar3D_SG.o $(OBJDIR)/cVector.o $(OBJDIR)/cVector3D_SG.o $(OBJDIR)/MetricElements_CSG.o $(OBJDIR)/ModelOperator.o $(OBJDIR)/ModelOperator_MF.o $(OBJDIR)/ModelOperator_File.o $(OBJDIR)/PreConditioner.o $(OBJDIR)/Solver.o $(OBJDIR)/PreConditioner_MF_CC.o $(OBJDIR)/PreConditioner_None.o $(OBJDIR)/Solver_QMR.o $(OBJDIR)/PreConditioner_MF_DC.o $(OBJDIR)/Solver_PCG.o $(OBJDIR)/Source.o $(OBJDIR)/Forward1D.o $(OBJDIR)/SourceMT_1D.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/DivergenceCorrection.o $(OBJDIR)/ForwardSolverIT.o $(OBJDIR)/ForwardSolverIT_DC.o $(OBJDIR)/test_Solver.o 


all: test_Solver 

# Here is the link step 
test_Solver: $(OBJDIR) $(OBJ) 
	 $(F90) -o $(OUTDIR)/test_Solver $(OBJ) $(LIBS_PATH) $(LIBS)

# Here are the compile steps 

$(OBJDIR): 
	mkdir -p $(OBJDIR)

$(OBJDIR)/Constants.o:Classes/Utils/Constants.f90  
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Utils/Constants.f90 -o $(OBJDIR)/Constants.o

$(OBJDIR)/Grid1D.o:Classes/1D/Grid1D.f90 $(OBJDIR)/Constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/1D/Grid1D.f90 -o $(OBJDIR)/Grid1D.o

$(OBJDIR)/Grid2D.o:Classes/2D/Grid2D.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid1D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/2D/Grid2D.f90 -o $(OBJDIR)/Grid2D.o

$(OBJDIR)/Grid.o:Classes/Grid/Grid.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid1D.o $(OBJDIR)/Grid2D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Grid/Grid.f90 -o $(OBJDIR)/Grid.o

$(OBJDIR)/rScalar.o:Classes/Vectors/rScalar.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/Constants.o $(OBJDIR)/Constants.o $(OBJDIR)/Constants.o $(OBJDIR)/Constants.o $(OBJDIR)/Constants.o $(OBJDIR)/Constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Vectors/rScalar.f90 -o $(OBJDIR)/rScalar.o

$(OBJDIR)/rVector.o:Classes/Vectors/rVector.f90 $(OBJDIR)/Constants.o $(OBJDIR)/rScalar.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Vectors/rVector.f90 -o $(OBJDIR)/rVector.o

$(OBJDIR)/ModelParameter1D.o:Classes/1D/ModelParameter1D.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid1D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/1D/ModelParameter1D.f90 -o $(OBJDIR)/ModelParameter1D.o

$(OBJDIR)/Esoln2DTM.o:Classes/2D/Esoln2DTM.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid2D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/2D/Esoln2DTM.f90 -o $(OBJDIR)/Esoln2DTM.o

$(OBJDIR)/ModelParameter2D.o:Classes/2D/ModelParameter2D.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid2D.o $(OBJDIR)/Esoln2DTM.o $(OBJDIR)/Grid1D.o $(OBJDIR)/ModelParameter1D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/2D/ModelParameter2D.f90 -o $(OBJDIR)/ModelParameter2D.o

$(OBJDIR)/MetricElements.o:Classes/ModelOperator/MetricElements.f90 $(OBJDIR)/Grid.o $(OBJDIR)/rVector.o $(OBJDIR)/rScalar.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelOperator/MetricElements.f90 -o $(OBJDIR)/MetricElements.o

$(OBJDIR)/ModelParameter.o:Classes/ModelParameter/ModelParameter.f90 $(OBJDIR)/Constants.o $(OBJDIR)/rScalar.o $(OBJDIR)/rVector.o $(OBJDIR)/Grid2D.o $(OBJDIR)/ModelParameter1D.o $(OBJDIR)/ModelParameter2D.o $(OBJDIR)/Grid.o $(OBJDIR)/MetricElements.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelParameter/ModelParameter.f90 -o $(OBJDIR)/ModelParameter.o

$(OBJDIR)/ModelReader.o:Classes/ModelParameter/ModelReader.f90 $(OBJDIR)/Grid.o $(OBJDIR)/ModelParameter.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelParameter/ModelReader.f90 -o $(OBJDIR)/ModelReader.o

$(OBJDIR)/Grid3D_SG.o:Classes/Grid/Grid3D_SG.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/Grid1D.o $(OBJDIR)/Grid2D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Grid/Grid3D_SG.f90 -o $(OBJDIR)/Grid3D_SG.o

$(OBJDIR)/rScalar3D_SG.o:Classes/Vectors/rScalar3D_SG.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/rScalar.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Vectors/rScalar3D_SG.f90 -o $(OBJDIR)/rScalar3D_SG.o

$(OBJDIR)/MatUtils.o:Classes/Utils/MatUtils.f90 $(OBJDIR)/Constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Utils/MatUtils.f90 -o $(OBJDIR)/MatUtils.o

$(OBJDIR)/rVector3D_SG.o:Classes/Vectors/rVector3D_SG.f90 $(OBJDIR)/Constants.o $(OBJDIR)/MatUtils.o $(OBJDIR)/Grid.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/rScalar.o $(OBJDIR)/rScalar3D_SG.o $(OBJDIR)/rVector.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Vectors/rVector3D_SG.f90 -o $(OBJDIR)/rVector3D_SG.o

$(OBJDIR)/ModelParameterCell_SG.o:Classes/ModelParameter/ModelParameterCell_SG.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/rScalar.o $(OBJDIR)/rVector.o $(OBJDIR)/rScalar3D_SG.o $(OBJDIR)/rVector3D_SG.o $(OBJDIR)/ModelParameter.o $(OBJDIR)/Grid1D.o $(OBJDIR)/Grid2D.o $(OBJDIR)/ModelParameter1D.o $(OBJDIR)/ModelParameter2D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelParameter/ModelParameterCell_SG.f90 -o $(OBJDIR)/ModelParameterCell_SG.o

$(OBJDIR)/ModelReader_Weerachai.o:Classes/ModelParameter/ModelReader_Weerachai.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/rScalar3D_SG.o $(OBJDIR)/ModelParameter.o $(OBJDIR)/ModelReader.o $(OBJDIR)/ModelParameter.o $(OBJDIR)/ModelParameterCell_SG.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelParameter/ModelReader_Weerachai.f90 -o $(OBJDIR)/ModelReader_Weerachai.o

$(OBJDIR)/cScalar.o:Classes/Vectors/cScalar.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/rScalar.o $(OBJDIR)/Constants.o $(OBJDIR)/Constants.o $(OBJDIR)/Constants.o $(OBJDIR)/Constants.o $(OBJDIR)/Constants.o $(OBJDIR)/Constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Vectors/cScalar.f90 -o $(OBJDIR)/cScalar.o

$(OBJDIR)/cScalar3D_SG.o:Classes/Vectors/cScalar3D_SG.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/cScalar.o $(OBJDIR)/rScalar3D_SG.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Vectors/cScalar3D_SG.f90 -o $(OBJDIR)/cScalar3D_SG.o

$(OBJDIR)/cVector.o:Classes/Vectors/cVector.f90 $(OBJDIR)/Constants.o $(OBJDIR)/rScalar.o $(OBJDIR)/rVector.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Vectors/cVector.f90 -o $(OBJDIR)/cVector.o

$(OBJDIR)/cVector3D_SG.o:Classes/Vectors/cVector3D_SG.f90 $(OBJDIR)/Constants.o $(OBJDIR)/MatUtils.o $(OBJDIR)/Grid.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/rScalar.o $(OBJDIR)/rScalar3D_SG.o $(OBJDIR)/rVector.o $(OBJDIR)/rVector3D_SG.o $(OBJDIR)/cVector.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Vectors/cVector3D_SG.f90 -o $(OBJDIR)/cVector3D_SG.o

$(OBJDIR)/MetricElements_CSG.o:Classes/ModelOperator/MetricElements_CSG.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/MetricElements.o $(OBJDIR)/rVector3D_SG.o $(OBJDIR)/rScalar3D_SG.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelOperator/MetricElements_CSG.f90 -o $(OBJDIR)/MetricElements_CSG.o

$(OBJDIR)/ModelOperator.o:Classes/ModelOperator/ModelOperator.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/cVector.o $(OBJDIR)/cScalar.o $(OBJDIR)/MetricElements.o $(OBJDIR)/ModelParameter.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelOperator/ModelOperator.f90 -o $(OBJDIR)/ModelOperator.o

$(OBJDIR)/ModelOperator_MF.o:Classes/ModelOperator/ModelOperator_MF.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/rScalar3D_SG.o $(OBJDIR)/cScalar3D_SG.o $(OBJDIR)/cVector3D_SG.o $(OBJDIR)/rVector.o $(OBJDIR)/rVector3D_SG.o $(OBJDIR)/MetricElements_CSG.o $(OBJDIR)/ModelParameter.o $(OBJDIR)/ModelOperator.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelOperator/ModelOperator_MF.f90 -o $(OBJDIR)/ModelOperator_MF.o

$(OBJDIR)/ModelOperator_File.o:Classes/ModelOperator/ModelOperator_File.f90 $(OBJDIR)/ModelOperator_MF.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelOperator/ModelOperator_File.f90 -o $(OBJDIR)/ModelOperator_File.o

$(OBJDIR)/PreConditioner.o:Classes/ModelOperator/PreConditioner.f90 $(OBJDIR)/cVector.o $(OBJDIR)/cScalar.o $(OBJDIR)/ModelOperator_MF.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelOperator/PreConditioner.f90 -o $(OBJDIR)/PreConditioner.o

$(OBJDIR)/Solver.o:Classes/Solver/Solver.f90 $(OBJDIR)/Constants.o $(OBJDIR)/cVector.o $(OBJDIR)/ModelOperator.o $(OBJDIR)/PreConditioner.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Solver/Solver.f90 -o $(OBJDIR)/Solver.o

$(OBJDIR)/PreConditioner_MF_CC.o:Classes/ModelOperator/PreConditioner_MF_CC.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/cVector.o $(OBJDIR)/cVector3D_SG.o $(OBJDIR)/cScalar3D_SG.o $(OBJDIR)/ModelOperator_MF.o $(OBJDIR)/PreConditioner.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelOperator/PreConditioner_MF_CC.f90 -o $(OBJDIR)/PreConditioner_MF_CC.o

$(OBJDIR)/PreConditioner_None.o:Classes/ModelOperator/PreConditioner_None.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/cVector.o $(OBJDIR)/cVector3D_SG.o $(OBJDIR)/cScalar3D_SG.o $(OBJDIR)/ModelOperator_MF.o $(OBJDIR)/PreConditioner.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelOperator/PreConditioner_None.f90 -o $(OBJDIR)/PreConditioner_None.o

$(OBJDIR)/Solver_QMR.o:Classes/Solver/Solver_QMR.f90 $(OBJDIR)/Solver.o $(OBJDIR)/cVector.o $(OBJDIR)/ModelOperator_MF.o $(OBJDIR)/PreConditioner_MF_CC.o $(OBJDIR)/PreConditioner_None.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Solver/Solver_QMR.f90 -o $(OBJDIR)/Solver_QMR.o

$(OBJDIR)/PreConditioner_MF_DC.o:Classes/ModelOperator/PreConditioner_MF_DC.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid.o $(OBJDIR)/cScalar.o $(OBJDIR)/cScalar3D_SG.o $(OBJDIR)/ModelOperator_MF.o $(OBJDIR)/PreConditioner.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ModelOperator/PreConditioner_MF_DC.f90 -o $(OBJDIR)/PreConditioner_MF_DC.o

$(OBJDIR)/Solver_PCG.o:Classes/Solver/Solver_PCG.f90 $(OBJDIR)/Solver.o $(OBJDIR)/cScalar.o $(OBJDIR)/ModelOperator_MF.o $(OBJDIR)/PreConditioner_MF_DC.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Solver/Solver_PCG.f90 -o $(OBJDIR)/Solver_PCG.o

$(OBJDIR)/Source.o:Classes/Source/Source.f90 $(OBJDIR)/cVector.o $(OBJDIR)/ModelOperator.o $(OBJDIR)/ModelParameter.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Source/Source.f90 -o $(OBJDIR)/Source.o

$(OBJDIR)/Forward1D.o:Classes/1D/Forward1D.f90 $(OBJDIR)/Constants.o $(OBJDIR)/Grid1D.o $(OBJDIR)/ModelParameter1D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/1D/Forward1D.f90 -o $(OBJDIR)/Forward1D.o

$(OBJDIR)/SourceMT_1D.o:Classes/Source/SourceMT_1D.f90 $(OBJDIR)/Constants.o $(OBJDIR)/cVector3D_SG.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/Source.o $(OBJDIR)/ModelOperator.o $(OBJDIR)/ModelParameter1D.o $(OBJDIR)/Forward1D.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/Source/SourceMT_1D.f90 -o $(OBJDIR)/SourceMT_1D.o

$(OBJDIR)/ForwardSolver.o:Classes/ForwardSolver/ForwardSolver.f90 $(OBJDIR)/Constants.o $(OBJDIR)/cVector.o $(OBJDIR)/Source.o $(OBJDIR)/Solver.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ForwardSolver/ForwardSolver.f90 -o $(OBJDIR)/ForwardSolver.o

$(OBJDIR)/DivergenceCorrection.o:Classes/ForwardSolver/DivergenceCorrection.f90 $(OBJDIR)/Solver_PCG.o $(OBJDIR)/ModelOperator.o $(OBJDIR)/cScalar3D_SG.o $(OBJDIR)/cVector3D_SG.o $(OBJDIR)/Source.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ForwardSolver/DivergenceCorrection.f90 -o $(OBJDIR)/DivergenceCorrection.o

$(OBJDIR)/ForwardSolverIT.o:Classes/ForwardSolver/ForwardSolverIT.f90 $(OBJDIR)/ForwardSolver.o $(OBJDIR)/DivergenceCorrection.o $(OBJDIR)/ModelOperator.o $(OBJDIR)/Solver_QMR.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ForwardSolver/ForwardSolverIT.f90 -o $(OBJDIR)/ForwardSolverIT.o

$(OBJDIR)/ForwardSolverIT_DC.o:Classes/ForwardSolver/ForwardSolverIT_DC.f90 $(OBJDIR)/ForwardSolverIT.o $(OBJDIR)/DivergenceCorrection.o $(OBJDIR)/ModelOperator.o $(OBJDIR)/Solver_QMR.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Classes/ForwardSolver/ForwardSolverIT_DC.f90 -o $(OBJDIR)/ForwardSolverIT_DC.o

$(OBJDIR)/test_Solver.o:test_Solver.f90 $(OBJDIR)/Constants.o $(OBJDIR)/ModelReader.o $(OBJDIR)/ModelReader_Weerachai.o $(OBJDIR)/ModelOperator_MF.o $(OBJDIR)/ModelOperator_File.o $(OBJDIR)/ModelParameterCell_SG.o $(OBJDIR)/Grid3D_SG.o $(OBJDIR)/cVector3D_SG.o $(OBJDIR)/cScalar3D_SG.o $(OBJDIR)/Solver_QMR.o $(OBJDIR)/Solver_PCG.o $(OBJDIR)/SourceMT_1D.o $(OBJDIR)/ForwardSolverIT.o $(OBJDIR)/DivergenceCorrection.o $(OBJDIR)/ForwardSolverIT_DC.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) test_Solver.f90 -o $(OBJDIR)/test_Solver.o

# Type " make clean " to get rid of all object and module files 
clean: 
	cd $(OBJDIR); \
	rm -f *~ *.o *.obj *.mod *.d *.s00 *.dbg *.stackdump \
	`find . -mindepth 1 -name "*~"` 

cleanall: clean 
	rm -f $(OUTDIR)/test_Solver 

src: clean 
	tar cvfz $(ARCHIVE).tgz * 

  
