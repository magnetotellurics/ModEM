# Makefile suited for building the Mod3DMT program
# Generated using: ./fmkmf.pl [OPTIONS] Mod3DMT.f90 > Makefile
# with command line options
# -p .:INV:LAPACK:SENS:UTILS:FIELDS/FiniteDiff3D:3D_MT:3D_MT/DICT:3D_MT/modelParam:3D_MT/FWD:3D_MT/FWD/Mod2d:3D_MT/ioMod
# -f90 gfortran (compiler)
# -opt -O3 -ffree-line-length-none (compiler optimisation)
# -lp /usr/lib64 (linking options: path to libraries)
# -l -llapack -lblas (linking options)
# -o ./objs/3D_MT/GFortRelease (output directory for object files)

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
include Makefile.local
OBJDIR = ./objs/3D_MT/IFortReleaseMPI
F90 = mpiifort #mpif90
FFLAGS = -O3 -w -xSSE4.2 -std03 -heap-arrays
MPIFLAGS = -cpp -DMPI
MODULE = -module $(OBJDIR)
LIBS_PATH = -L/opt/intel/compilers_and_libraries_2019.0.117/linux/mkl/lib/intel64_lin
LIBS = -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -qopenmp -lpthread
OF = NLCG_OF2.f90
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
# include Makefile.local
# OBJDIR = ./objs/3D_MT/GFortRelease
# F90 = gfortran 
# FFLAGS = -O3 -ffree-line-length-none
# MPIFLAGS = -x f95-cpp-input
# MODULE = --sysroot=$(OBJDIR)
# LIBS_PATH = -L/usr/lib64
# LIBS = -llapack -lblas

# -------------------End-macro-Defs---------------------------
OBJ = $(OBJDIR)/math_constants.o $(OBJDIR)/elements.o $(OBJDIR)/utilities.o $(OBJDIR)/file_units.o $(OBJDIR)/anisConstraint.o $(OBJDIR)/DataSpace.o $(OBJDIR)/GridDef.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_scalar.o $(OBJDIR)/sg_spherical.o $(OBJDIR)/GridCalc.o $(OBJDIR)/sg_sparse_vector.o $(OBJDIR)/MPI_declaration.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/SensMatrix.o $(OBJDIR)/sg_boundary.o $(OBJDIR)/transmitters.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/receivers.o $(OBJDIR)/dataTypes.o $(OBJDIR)/SolverSens.o $(OBJDIR)/sg_diff_oper.o $(OBJDIR)/Fwd1DANI.o $(OBJDIR)/Fwd2DANI.o $(OBJDIR)/boundary_Pek.o $(OBJDIR)/modelOperator3D.o $(OBJDIR)/EMfieldInterp.o $(OBJDIR)/DataFunc.o $(OBJDIR)/DataSens.o $(OBJDIR)/EMsolve3D.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/SensComp.o $(OBJDIR)/SymmetryTest.o $(OBJDIR)/UserCtrl.o $(OBJDIR)/ioAscii.o $(OBJDIR)/DataIO.o $(OBJDIR)/Main.o $(OBJDIR)/MPI_sub.o $(OBJDIR)/MPI_main.o $(OBJDIR)/NLCG.o $(OBJDIR)/LBFGS.o $(OBJDIR)/DCG.o $(OBJDIR)/Mod3DMT.o 


all: Mod3DMT 

# Here is the link step 
Mod3DMT: $(OBJDIR) $(OBJ) 
	 $(F90) -o $(OUTDIR)/MT3DANI_P6 $(OBJ) $(LIBS_PATH) $(LIBS)

# Here are the compile steps 

$(OBJDIR): 
	mkdir -p $(OBJDIR)

$(OBJDIR)/math_constants.o:math_constants.f90  
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) math_constants.f90 -o $(OBJDIR)/math_constants.o

$(OBJDIR)/elements.o:elements.f90  
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) elements.f90 -o $(OBJDIR)/elements.o

$(OBJDIR)/utilities.o:utilities.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) utilities.f90 -o $(OBJDIR)/utilities.o

$(OBJDIR)/file_units.o:file_units.f90  
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) file_units.f90 -o $(OBJDIR)/file_units.o

$(OBJDIR)/anisConstraint.o:anisConstraint.f90 $(OBJDIR)/math_constants.o
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) anisConstraint.f90 -o $(OBJDIR)/anisConstraint.o

$(OBJDIR)/DataSpace.o:DataSpace.f90 $(OBJDIR)/utilities.o $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) DataSpace.f90 -o $(OBJDIR)/DataSpace.o

$(OBJDIR)/GridDef.o:GridDef.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) GridDef.f90 -o $(OBJDIR)/GridDef.o

$(OBJDIR)/sg_vector.o:sg_vector.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/GridDef.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) sg_vector.f90 -o $(OBJDIR)/sg_vector.o

$(OBJDIR)/sg_scalar.o:sg_scalar.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/GridDef.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) sg_scalar.f90 -o $(OBJDIR)/sg_scalar.o

$(OBJDIR)/sg_spherical.o:sg_spherical.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/sg_scalar.o $(OBJDIR)/sg_vector.o
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) sg_spherical.f90 -o $(OBJDIR)/sg_spherical.o

$(OBJDIR)/GridCalc.o:GridCalc.f90 $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_scalar.o $(OBJDIR)/elements.o $(OBJDIR)/sg_spherical.o
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) GridCalc.f90 -o $(OBJDIR)/GridCalc.o

$(OBJDIR)/sg_sparse_vector.o:sg_sparse_vector.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/GridDef.o $(OBJDIR)/sg_vector.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) sg_sparse_vector.f90 -o $(OBJDIR)/sg_sparse_vector.o

$(OBJDIR)/MPI_declaration.o:MPI_declaration.f90  
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) MPI_declaration.f90 -o $(OBJDIR)/MPI_declaration.o

$(OBJDIR)/ModelSpace.o:ModelSpace.f90 $(OBJDIR)/GridCalc.o $(OBJDIR)/file_units.o $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/sg_scalar.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_sparse_vector.o $(OBJDIR)/anisConstraint.o modelCov/RecursiveAR.hd ModelMap.inc modelCov/RecursiveAR.inc modelParamIO/WS.inc
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) ModelSpace.f90 -o $(OBJDIR)/ModelSpace.o

$(OBJDIR)/SensMatrix.o:SensMatrix.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/file_units.o $(OBJDIR)/utilities.o $(OBJDIR)/DataSpace.o $(OBJDIR)/ModelSpace.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) SensMatrix.f90 -o $(OBJDIR)/SensMatrix.o

$(OBJDIR)/sg_boundary.o:sg_boundary.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/GridDef.o $(OBJDIR)/sg_vector.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) sg_boundary.f90 -o $(OBJDIR)/sg_boundary.o

$(OBJDIR)/transmitters.o:transmitters.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) transmitters.f90 -o $(OBJDIR)/transmitters.o

$(OBJDIR)/SolnSpace.o:SolnSpace.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_boundary.o $(OBJDIR)/sg_sparse_vector.o $(OBJDIR)/transmitters.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) SolnSpace.f90 -o $(OBJDIR)/SolnSpace.o

$(OBJDIR)/receivers.o:receivers.f90 $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) receivers.f90 -o $(OBJDIR)/receivers.o

$(OBJDIR)/dataTypes.o:dataTypes.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) dataTypes.f90 -o $(OBJDIR)/dataTypes.o

$(OBJDIR)/SolverSens.o:SolverSens.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/transmitters.o $(OBJDIR)/dataTypes.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) SolverSens.f90 -o $(OBJDIR)/SolverSens.o

$(OBJDIR)/sg_diff_oper.o:sg_diff_oper.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_scalar.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) sg_diff_oper.f90 -o $(OBJDIR)/sg_diff_oper.o

$(OBJDIR)/Fwd1DANI.o:Fwd1DANI.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Fwd1DANI.f90 -o $(OBJDIR)/Fwd1DANI.o

$(OBJDIR)/Fwd2DANI.o:Fwd2DANI.f90 $(OBJDIR)/Fwd1DANI.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Fwd2DANI.f90 -o $(OBJDIR)/Fwd2DANI.o

$(OBJDIR)/boundary_Pek.o:boundary_Pek.f90 $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_scalar.o $(OBJDIR)/sg_boundary.o $(OBJDIR)/Fwd2DANI.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) boundary_Pek.f90 -o $(OBJDIR)/boundary_Pek.o

$(OBJDIR)/modelOperator3D.o:modelOperator3D.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/GridCalc.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_boundary.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/boundary_Pek.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) modelOperator3D.f90 -o $(OBJDIR)/modelOperator3D.o

$(OBJDIR)/EMfieldInterp.o:EMfieldInterp.f90 $(OBJDIR)/utilities.o $(OBJDIR)/sg_sparse_vector.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/modelOperator3D.o
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) EMfieldInterp.f90 -o $(OBJDIR)/EMfieldInterp.o

$(OBJDIR)/DataFunc.o:DataFunc.f90 $(OBJDIR)/EMfieldInterp.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/receivers.o $(OBJDIR)/transmitters.o $(OBJDIR)/dataTypes.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) DataFunc.f90 -o $(OBJDIR)/DataFunc.o

$(OBJDIR)/DataSens.o:DataSens.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/DataSpace.o $(OBJDIR)/DataFunc.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) DataSens.f90 -o $(OBJDIR)/DataSens.o

$(OBJDIR)/EMsolve3D.o:EMsolve3D.f90 $(OBJDIR)/sg_boundary.o $(OBJDIR)/sg_diff_oper.o $(OBJDIR)/sg_sparse_vector.o $(OBJDIR)/modelOperator3D.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/sg_vector.o $(OBJDIR)/sg_scalar.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) EMsolve3D.f90 -o $(OBJDIR)/EMsolve3D.o

$(OBJDIR)/ForwardSolver.o:ForwardSolver.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/DataFunc.o $(OBJDIR)/DataSpace.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/EMsolve3D.o $(OBJDIR)/transmitters.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) ForwardSolver.f90 -o $(OBJDIR)/ForwardSolver.o

$(OBJDIR)/DataIO.o:DataIO.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/file_units.o $(OBJDIR)/utilities.o $(OBJDIR)/DataSpace.o $(OBJDIR)/transmitters.o $(OBJDIR)/receivers.o $(OBJDIR)/dataTypes.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) DataIO.f90 -o $(OBJDIR)/DataIO.o

$(OBJDIR)/SensComp.o:SensComp.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/SensMatrix.o $(OBJDIR)/DataSens.o $(OBJDIR)/SolverSens.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/DataIO.o
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) SensComp.f90 -o $(OBJDIR)/SensComp.o

$(OBJDIR)/SymmetryTest.o:SymmetryTest.f90 $(OBJDIR)/SensComp.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) SymmetryTest.f90 -o $(OBJDIR)/SymmetryTest.o

$(OBJDIR)/UserCtrl.o:./UserCtrl.f90 $(OBJDIR)/utilities.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) ./UserCtrl.f90 -o $(OBJDIR)/UserCtrl.o

$(OBJDIR)/ioAscii.o:ioAscii.f90 $(OBJDIR)/GridDef.o $(OBJDIR)/math_constants.o $(OBJDIR)/EMsolve3D.o $(OBJDIR)/DataSpace.o $(OBJDIR)/DataFunc.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/transmitters.o $(OBJDIR)/receivers.o $(OBJDIR)/dataTypes.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) ioAscii.f90 -o $(OBJDIR)/ioAscii.o

$(OBJDIR)/Main.o:Main.f90 $(OBJDIR)/ModelSpace.o $(OBJDIR)/DataSpace.o $(OBJDIR)/DataFunc.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/SensMatrix.o $(OBJDIR)/UserCtrl.o $(OBJDIR)/ioAscii.o $(OBJDIR)/DataIO.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Main.f90 -o $(OBJDIR)/Main.o

$(OBJDIR)/MPI_sub.o:MPI_sub.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/utilities.o $(OBJDIR)/SolnSpace.o $(OBJDIR)/UserCtrl.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/MPI_declaration.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) MPI_sub.f90 -o $(OBJDIR)/MPI_sub.o

$(OBJDIR)/MPI_main.o:MPI_main.f90 $(OBJDIR)/math_constants.o $(OBJDIR)/file_units.o $(OBJDIR)/utilities.o $(OBJDIR)/DataSens.o $(OBJDIR)/SolverSens.o $(OBJDIR)/ForwardSolver.o $(OBJDIR)/SensComp.o $(OBJDIR)/MPI_declaration.o $(OBJDIR)/MPI_sub.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) MPI_main.f90 -o $(OBJDIR)/MPI_main.o

$(OBJDIR)/NLCG.o:$(OF) $(OBJDIR)/SensComp.o $(OBJDIR)/DataIO.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) $(OF) -o $(OBJDIR)/NLCG.o

$(OBJDIR)/LBFGS.o:LBFGS.f90 $(OBJDIR)/SensComp.o $(OBJDIR)/DataIO.o $(OBJDIR)/MPI_main.o $(OBJDIR)/MPI_sub.o
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) LBFGS.f90 -o $(OBJDIR)/LBFGS.o
$(OBJDIR)/DCG.o:DCG.f90 $(OBJDIR)/SensComp.o $(OBJDIR)/Main.o $(OBJDIR)/utilities.o $(OBJDIR)/math_constants.o 
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) DCG.f90 -o $(OBJDIR)/DCG.o

$(OBJDIR)/Mod3DMT.o:Mod3DMT.f90 $(OBJDIR)/SensComp.o $(OBJDIR)/SymmetryTest.o $(OBJDIR)/Main.o $(OBJDIR)/NLCG.o $(OBJDIR)/LBFGS.o $(OBJDIR)/DCG.o
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) Mod3DMT.f90 -o $(OBJDIR)/Mod3DMT.o

# Type " make clean " to get rid of all object and module files 
clean: 
	cd $(OBJDIR); \
	rm -f *~ *.o *.obj *.mod *.d *.s00 *.dbg *.stackdump \
	`find . -mindepth 1 -name "*~"` 

cleanall: clean 
	rm -f $(OUTDIR)/MT3DANI_P6 

src: clean 
	tar cvfz $(ARCHIVE).tgz * 

  
