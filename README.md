
ModEM_OO - User Guide | Package Instructions - 04/09/2024

    ModEM (Modular Electromagnetic) is a software for inversion of 3D EM datasets,
    originally developed at Oregon State University using Fortran 90,
    with numerous extensions of the code by various collaborators over the past decade.
    In order to put together the best approaches developed since its conception,
    into a leaner and more easily extended/maintained code, a new ModEM version 
    has been implemented, using Fortran 2003 Object Oriented features.
    Result of a partnership project between Observatório Nacional and Shell-Brasil,
    this version (referred to as ModEM_OO) implements new capabilities including:

        - Controlled Source Electromagnetic (CSEM) datasets.
        - Anisotropic models with Vertical Transverse Isotropy (VTI).
        - The Multi-Resolution (MR) Grids described in Cherevotova et al. (2018).

    The code is of course not nearly as throughly tested as the stable Mod3DMT code previously
    distributed to academic users. Additional work will be required to improved efficiency
    and parallelization. Furthermore, due to copyright issues we do not provide the 1D solver
    required for CSEM (primary field computations). Further details on making use of CSEM
    capabilities are provided in a separate document.

    Here we present a quick User Guide for the first release of ModEM_OO.
    Note that in the following, model and data files are in the standard ASCII formats,
    same as used for the stable Mod3DMT code.


    1. Getting Started

        1.1. Extract ModEM_OO-1.0.0/

            1.1.1. For Unix systems, follow the steps from 1.2. using a terminal.

            1.1.2. For Microsoft systems (Windows OS),
                   we recommend installing a latest Windows Subsystem for Linux (WSL).
                   Preferably above Ubuntu 20.04.5 LTS OS. Freely obtained from Microsoft Store.
                   Then follow the next steps using the WSL command prompt.

        1.2. Installing ModEM dependencies (make; gfortran; MPI; blas, lapack and fftw3 libs).

            - Run (or follow) the Shell script provided in the scripts/ folder:

                bash ../scripts/install_dependencies.sh

        1.3. Compiling|Installing ModEM_OO.

            - Must be done within the src/ folder.

            - To ensure permissions, run: 

                chmod 777 Build/*

            1.3.1. Creating Makefiles and Compiling executables:

                - To create both ModEM_SERIAL and ModEM_MPI executables at once, run:

                    bash ../scripts/compile_serial_mpi.sh

                Notes:
                    - Users who just want to use ModEM_OO can go to step 2.
                    - Executables and Makefiles can be renamed at will.
                    - It is possible that some users may have to create Makefiles,
                      and then edit path names for libraries and include files.

                - To create only ModEM_MPI executable, run:

                    bash ../scripts/compile_mpi.sh

                - To create only ModEM_SERIAL executable, run:

                    bash ../scripts/compile_serial.sh

                - To generate one Makefile (Serial or MPI), run:

                    ./Build/Configure.modem.ubuntu.gfortran MakefileSerial SERIAL ModEM.f90

                    ./Build/Configure.modem.ubuntu.gfortran MakefileMPI MPI ModEM.f90

                - To create a ModEM executable according to a Makefile, run:

                    make -f <MAKEFILE>

                    - For Mac OS only (Serial or MPI), run:

                        make -f Build/MakefileSerial_MAC

                        make -f Build/MakefileMPI_MAC ????


    2. Running ModEM.

        Code execution is controlled (as in Mod3DMT) through command line arguments, and control
        files setting options and parameters. In contrast to Mod3DMT, command line arguments can
        be given in any order. Note also that a range of options for grids, model operators,
        solvers, and forward modeling approaches are supported in a single executable and these are
        set in the forward control file, which will thus almost always be necessary.

        - To use the Serial Version (from Makefile or Script respectively), just run:

            ./ModEM or ./ModEM_SERIAL

        - For MPI version, use mpirun, providing the number of processes (NPROCS >= 2):

            mpirun -np <NPROCS> ./ModEM_MPI

        Note: If everything was done correctly, minimal usage instructions will then appear on the screen.
                Exactly as shown in 3.1.

        2.1. ModEM command line arguments.

            - A full list of supported arguments (3.2.) can be obtained by running:

                ./ModEM --help

        2.2. ModEM Main Jobs.

            - Different jobs are called using specific lone flags (shown in 3.2.),
              which will tell the program to perform one of the following tasks:

            2.2.1. [-f]    Forward Modeling (FWD).

                - To perform a FWD job using default parameters,
                    just Model and Data Files are needed, preceded by -m and -d flags respectively, as follows:

                    ./ModEM -f -m <MODEL_PATH> -d <DATA_PATH>

                - This generates an output Data File, named 'all_predicted_data.dat',
                    in the same format as the inputed one, containing the data predicted by the FWD.

                2.2.1.1. FWD supported options.

                    - Other file paths can be specified, preceded by flags, for a FWD command line:

                        [-pd], [--predicted]   Sets the path to output the Predicted Data File.

                        [-es], [--esolution]   Enables and sets a path to output a binary file,
                                               containing solved electric fields for all transmitters.

                        [-cf], [--ctrl_fwd]    Sets the path to a input Control File with FWD parameters.
                                               (Shown in 2.2.5.1.)

            2.2.2. [-i]    Inversion (INV).

                - To perform a INV job using default parameters,
                  just Model and Data Files paths are required after -m and -d flags respectively, as follows:

                    ./ModEM -i -m <MODEL_PATH> -d <DATA_PATH>

                - This generates an output Folder, named by the INV type plus execution date and time.
                  Containing INV log files and the output files resultant from each iteration.

                    - PredictedData_it#.dat:        ????
                    - ResidualData_it#.res:         ????
                    - SigmaModel_it#.rho:           ????
                    - PerturbationModel_it#.prm:    ????

                2.2.2.1. INV supported options.

                    - Other file paths can be specified, preceded by flags, for an INV command line:

                        [-pm], [--pmodel]     Sets an input perturbation model file path ????

                        [-c],  [--cov]        Sets the path to input a Covariance File, in ModEM format ????

                        [-o],  [--outdir]     Sets the path to the output Folder.

                        [-ci], [--ctrl_inv]   Sets the path to a input Control File,
                                              with specific INV parameters, as detailed in 2.2.5.2.

            2.2.3. [-j]    Jacobian Multiplication (jMult).

                - To run a default jMult job, two Model Files are required: Starting and Perturbation;
                  plus one Data File. Preceded by the flags as bellow:

                    ./ModEM -j -m <START_MODEL_PATH> -pm <PERT_MODEL_PATH> -d <DATA_PATH>

                - This job results one Data File, named 'jmhat.dat',
                  formatted as the inputed one, containing jMult calculations ????

                2.2.3.1. jMult supported options.

                    - More flagged file paths can be used in a jMult command line:

                        [-jm], [--jmhat]      Sets the path to the resultant JmHat Data File.

                        [-cf], [--ctrl_fwd]   Sets the path to an input Control File with FWD parameters.
                                              (Shown in 2.2.5.1.)

            2.2.4. [-jt]    Transposed J Multiplication (jMult_T).

                - To perform jMult_T job using default parameters, Model and Data Files are needed,
                  preceded by -m and -d flags respectively, as follows:

                    ./ModEM -jt -m <MODEL_PATH> -d <DATA_PATH>

                - Outputting a single Model File, named 'dsigma.rho',
                  in the same dimensions as the inputed one, containing ???? calculated by jMult_T.

                2.2.4.1. jMult_T supported options.

                    - Other paths preceded by flags fit on a jMult_T execution line:

                        [-dm], [--dmodel]     Sets the path to output dSigma ???? Model File.

                        [-cf], [--ctrl_fwd]   Sets the path to a input Control File with FWD parameters.
                                              (Shown in 2.2.5.1.)

            2.2.5. Control File Templates.

                - FWD and INV can have their control parameters modified through input text files.
                    Each featuring specific variables, that can be created by running:

                    ./ModEM --template

                - This create two distinct text files:

                2.2.5.1. Forward Modeling Control File - 'control.fwd'

                    - Making it possible to edit core features that control FWD, such as:

                        model_operator_type [MF|SP|SP2] : MF 
                        grid_format : 0,a,1,b,2,c
                        model_method [mirror|fixed height] : fixed height
                        model_n_air_layer [10]             : 10
                        model_max_height [200.]            : 200.
                        source_type_mt [1D|2D]                                              : 1D
                        source_type_csem [EM1D|Dipole1D]                                    : EM1D
                        get_1d_from [Fixed_Value|Geometric_Mean|Mean_around_Tx|Tx_Position] : Geometric_Mean
                        solver_type [QMR|BICG]         : BICG
                        forward_solver_type [IT|IT_DC] : IT_DC
                        max_solver_iters [80]          : 80
                        max_solver_calls [20]          : 20
                        max_divcor_iters [100]         : 100
                        tolerance_solver [1E-7]        : 1E-7
                        tolerance_divcor [1E-5]        : 1E-5

                2.2.5.2. Inversion Control File - 'control.inv'

                    - Allowing modification of values of the variables that govern INV iterations, such as:

                        inversion_type [DCG|NLCG]       : NLCG
                        joint_type [Unweighted|TxBased] : Unweighted
                        max_inv_iters [100]             : 100
                        max_grad_iters [20]             : 20
                        error_tol [1E-3]                : 1E-3
                        rms_tol [1.05]                  : 1.05
                        lambda [1.]                     : 1.
                        lambda_tol [1.0e-4]             : 1.0e-4
                        lambda_div [10.]                : 10.
                        startdm [20.]                   : 20.

                Notes:
                    - Templates are generated containing the exact same default values used by ModEM_OO.
                    - Allowed values and ranges are described in brackets for each parameter,
                      setting them outside of these throws runtime errors.
                    - Any parameter can be removed or commented out, causing the use of its default value.


    3. ModEM_OO Standards.

        3.1. ModEM_OO_1.0.0 Minimal Usage:

                    Forward Modeling (FWD):
                        <ModEM_OO> -f -m <rFile_Model> -d <rFile_Data>
                        Output:
                        - 'all_predicted_data.dat' or the path specified by      [-pd]

                    Inversion (INV):
                        <ModEM_OO> -i -m <rFile_Model> -d <rFile_Data>
                        Output:
                        - directory named 'Output_<date>_<time>' or specified by [-o]

                    Jacobian Multiplication (JMult):
                        <ModEM_OO> -j -m <rFile_Model> -pm <rFile_pModel> -d <rFile_Data>
                        Output:
                        - 'jmhat.dat' or the path specified by                   [-jm]

                    Transposed J Multiplication (JMult_T):
                        <ModEM_OO> -jt -m <rFile_Model> -d <rFile_Data>
                        Output:
                        - 'dsigma.rho' or the path specified by                  [-dm]

                    Other options:
                        <ModEM_OO> -h or <ModEM_OO> --help

            Note: See Mod3DMT documentation (at docs/papers/ModEM_UserGuide.pdf)
                  for further details about Jmult, JMult_T options.

        3.2. ModEM_OO_1.0.0 Options:

                    Flags to define a job:
                        [-f],  [--forward]   :  Forward Modeling.
                        [-j],  [--jmult]     :  Jacobian Multiplication.
                        [-jt], [--jmult_t]   :  Transposed Jacobian Multiplication.
                        [-i],  [--inversion] :  Inversion.
                        [-s],  [--split]     :  Split one input model into single axes.

                    Other arguments:
                        [-d],  [--data]      :  Flag to precede data file path.
                        [-m],  [--model]     :  Flag to precede model file path.
                        [-pm], [--pmodel]    :  Flag to precede perturbation model file path.
                        [-c],  [--cov]       :  Flag to precede covariance file path.
                        [-cf], [--ctrl_fwd]  :  Flag to precede forward control file path.
                        [-ci], [--ctrl_inv]  :  Flag to precede inversion control file path.
                        [-o],  [--outdir]    :  Flag to precede output directory path.
                        [-dm], [--dmodel]    :  Flag to precede output dsigma model file path.
                        [-pd], [--predicted] :  Flag to precede output predicted data file path.
                        [-jm], [--jmhat]     :  Flag to precede output JmHat data file path.
                        [-es], [--esolution] :  Flag to precede binary output e-solution file path.
                        [-v],  [--version]   :  Print version.
                        [-h],  [--help]      :  Print this information.
                        [-tmp],[--template]  :  Create control file templates.

        3.3. Using VTI anisotropic Models.

            - ModEM standard format presents its Model Files with a four line header,
                describing grid dimensions (nx, ny, nz) and cell sizes for each dimension.
                Followed by a block of nx*ny*nz resistivity values, ended by 3D origin and 1D rotation lines.
                VTI Models take into account different horizontal and vertical resistivities,
                It is enabled in ModEM_OO by the following changes to its Model Files:

                    ---------------------------------------------
                    1st header line:  nx    ny    nz   ? type
                    ---------------------------------------------
                          Isotropic: <nx>  <ny>  <nz>  0 LOGE
                    Anisotropic VTI: <nx>  <ny>  <nz>  0 LOGE VTI

            - Then instead of one, two blocks of values in sequence with the same dimension (nx, ny, nz).
              The first referring to horizontal resistivities and the second to the vertical ones.

        3.4. Using Multi-Resolution grids.

            - ModEM's first approach to work with Multi-Resolution grids, aims to gain performance
              by reducing the grid spacing at depth. By coarsening some layers of the original grid,
              enlarging its cells by a factor of 2, keeping best resolution for interesting layers.
              This follows the development given in Cherevatova et al., (2018).

              Resolution of an input grid can be modified through the 'grid_format' parameter, 
              enabled if present and uncommented in the FWD Control File.
              This variable expects an array with G [factor, depth] integer pairs, where:

                  - G      = number of layer groups, into which the original grid will be subdivided.
                  - factor = coarse intensity, groups are coarsened exponentially according to 2** factor,
                             where 0 represents the highest resolution. 
                  - depth  = amount of layers in a group.

                  For instance:

                      grid_format : 0,a,1,b,2,c

                  - Sets 3 layer groups, with a, b, c layers of depth
                    and coarseness factors of 0, 1 and 2, respectively.

                      grid_format : 0,10,1,10,0,10,2,10

                  - Subdivides the grid in four groups with ten layers of depth.
                    Leaving the first and third groups with the original (best) resolution.
                    Reducing by half the resolution of the second group,
                    and making the last ten layers four times coarser.

            - Notes:
                   - At present, it is not possible to coarsen the grid air layers,
                     which are always added to the first group of layers at the top of the model.
                   - Coarseness can only change by one level (factor of 2) between adjacent grids.

    SOMEWHERE NEED TO SAY SOMETHING ABOUT 1D background model for CSEM -- I.e., how to get Dipole1D and install.
    And how do we deal with makefile, etc. for a version w/o CSEM ????


    4. Package extra features.

        4.1. Test Sets containing validated examples of ModEM_OO input files:
                models, data, covariances and controls, can be found in the inputs/ folder.

        4.2. Published articles and additional information about ModEM can be checked at the docs/ folder.

        4.3. ModEM add-on programs are provided in the tools/ folder, such as:

            4.3.1. 3DGrid ????

                - NEED TO ASK NASER.   From my perspective make  NO mention of 3DGrid.   It is not my program!

            4.3.2. A JavaScript tool for electrical solution visualization ????
                Certainly not now!   

        4.4. Extensive Doxygen Source Code documentation can be found at:

            https://on.multiphysics.gitlab.io/modem-oo/

        4.5. GitLab releases of packages at:

            https://gitlab.com/on.multiphysics/modem-oo/-/releases/

        4.6. Contacts.

            sergio@on.br                   - Sergio Fontes    - Project's Coordinator.
            
            egbert@coas.oregonstate.edu    - Gary Egbert      - Creator of ModEM.
            
            paulowerdt@on.br               - Paulo Werdt      - ModEM_OO Developer.
