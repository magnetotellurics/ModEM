!
!> Prototype serial version for testing the ModEM-OO program
!
!> with three main jobs implemented: ForwardModelling, Adjoint and JMult_T
!
!> Obs.: To merge with the MPI version in the future
!
program TestSerial
    !
    use ModEMControlFile
    !
    use ModelReader_Weerachai
    !
    use DataFileStandard
    !
    use Sensitivity
    !
    use ModelCovarianceRec
    !
    real( kind=prec ) :: t_start, t_finish
    !
    !> Start the program and runtime count
    call cpu_time( t_start )
    !
    modem_job = "unknown"
    !
    call setupDefaultParameters()
    !
    !> Validate arguments, set model_file_name, data_file_name, control_file_name, etc...
    call handleArguments()
    !
    write( *, * )
    write( *, * ) "Start ModEM-OO."
    write( *, * )
    !
    !> Check parameters at the control file
    if( has_control_file ) call handleControlFile()
    !
    !> Execute the job specified in the arguments
    call handleJob()
    !
    !> Deallocate remaining main program memory
    call garbageCollector()
    !
    !> End runtime countdown
    call cpu_time( t_finish )
    !
    write( *, * )
    write( *, "(A18, F10.2, A3)" ) "Finish ModEM-OO (", t_finish - t_start, "s)"
    write( *, * )
    !
contains
    !
    !> Run ForwardModelling: Calculate ESolution for all transmitters and...
    !>     Receive a flag indicating whether it will be used for Adjoint:
    !>         - False or None: Calculate only the predicted data for each transmitter-receiver pair.
    !>         - True: Calculate LRows in the receivers after calculating predicted data.
    !> Obs.: Require the previous definition of the global ForwardSolver (createForwardSolver())
    !
    subroutine ForwardModeling( sigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        !
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        integer :: i_tx, n_tx, i_rx
        !
        ! Verbose
        write( *, * ) "          - Start Forward Modeling"
        !
        !>
        n_tx = size( transmitters )
        !
        ! Verbose
        write( *, * ) "          > Write ESolution to file: [", e_solution_file_name, "]"
        !
        !> Write the first header of the ESolution binary file, according to the first transmitter
        call writeEsolutionHeader( n_tx, transmitters(1)%Tx%n_pol )
        !
        !> Loop over all Transmitters
        do i_tx = 1, n_tx
            !
            !> Pointer to the Transmitter
            Tx => getTransmitter( i_tx )
            !
            !> Set Transmitter's ForwardSolver Period) and Conductivity
            call Tx%forward_solver%setFrequency( sigma, Tx%period )
            !
            !> Instantiate Transmitter's Source - According to transmitter type or chosen via control file
            select type( Tx )
                !
                class is( TransmitterMT_t )
                    !
                    call Tx%setSource( SourceMT_1D_t( model_operator, sigma, Tx%period ) )
                    !
                class is( TransmitterCSEM_t )
                    !
                    call Tx%setSource( SourceCSEM_Dipole1D_t( model_operator, sigma, Tx%period, Tx%location, Tx%dip, Tx%azimuth, Tx%moment ) )
                    !
                class default
                    stop "Error: ForwardModeling > Unclassified Transmitter"
                !
            end select
            !
            !> Build Source E according to source type
            call Tx%source%createE()
            !
            !> Solve ESolution for this Transmitter
            call Tx%solve()
            !
            ! Verbose
            write( *, * ) "                   - Receivers calculation"
            !
            !> Loop for each Receiver related to this Transmitter
            do i_rx = 1, size( Tx%receiver_indexes )
                !
                !> Pointer to the Tx Receiver
                Rx => getReceiver( Tx%receiver_indexes( i_rx ) )
                !
                call Rx%predictedData( Tx )
                !
                call all_predicted_data( i_tx )%set( Rx%data_group )
                !
            enddo
            !
        enddo
        !
        ! Verbose
        write( *, * ) "          - Finish Forward Modeling"
        !
    end subroutine ForwardModeling
    !
    !> Routine to run a full ForwardModeling job and deliver the result (PredictedData) in a text file
    subroutine jobForwardModeling()
        implicit none
        !
        ! Verbose
        write( *, * ) "     - Start jobForwardModeling"
        !
        !> Reads Model File: instantiates Grid, ModelOperator and ModelParameter
        if( .NOT. has_model_file ) then 
            stop "Error: jobForwardModeling > Missing Model file!"
        else
            call handleModelFile()
            !
        endif
        !
        !> Reads Data File: instantiates and builds the Data relation between Txs and Rxs
        if( .NOT. has_data_file ) then 
            stop "Error: jobForwardModeling > Missing Data file!"
        else
            call handleDataFile()
        endif
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        call createForwardSolver()
        !
        !> Run ForwardModelling only to calculate predicted data
        call ForwardModeling( sigma0 )
        !
        !> Write all_predicted_data, with its proper Rx headers, to the file <predicted_data_file_name>
        call writeDataGroupArray( all_predicted_data, predicted_data_file_name )
        !
        ! Verbose
        write( *, * ) "     - Finish jobForwardModeling"
        !
    end subroutine jobForwardModeling
    !
    !> Routine to run a full JMult job and deliver the result (JmHat_data Data) in a text file
    subroutine JobJMult()
        implicit none
        !
        !> ????
        type( DataGroupTx_t ), allocatable, dimension(:) :: JmHat_data
        !
        ! Verbose
        write( *, * ) "     - Start JobJMult"
        !
        !> Read Prior Model File: instantiates PModelParameter
        if( .NOT. has_pmodel_file ) then 
            stop "Error: JobJMult > Missing Prior Model file!"
        else
            !
            call handlePModelFile()
            !
            !> Read Model File: instantiates Grid, ModelOperator and ModelParameter
            if( .NOT. has_model_file ) then 
                stop "Error: JobJMult > Missing Model file!"
            else
                call handleModelFile()
                !
                call pmodel%setMetric( model_operator%metric )
                !
            endif
            !
        endif
        !
        !> Reads Data File: instantiates and builds the Data relation between Txs and Rxs
        if( .NOT. has_data_file ) then 
            stop "Error: JobJMult > Missing Data file!"
        else
            call handleDataFile()
        endif
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        call createForwardSolver()
        !
        !> Run ForwardModelling for calculate Predicted Data and LRows
        call ForwardModeling( sigma0 )
        !
        !> Write all_predicted_data, with its proper Rx headers, to the file <predicted_data_file_name>
        call writeDataGroupArray( all_predicted_data, predicted_data_file_name )
        !
        !> Initialize JmHat_data data array in the same format as the predicted data array
        allocate( JmHat_data, source = all_predicted_data )
        !
        !> Fill JmHat_data array from JMult routine
        call JMult( sigma0, pmodel, JmHat_data )
        !
        !> Write JmHat_data, with its proper Rx headers, to the file <JmHat_data_file_name>
        call writeDataGroupArray( JmHat_data, JmHat_data_file_name )
        !
        !> Deallocates specific variables from the adjoint job
        deallocate( JmHat_data )
        !
        ! Verbose
        write( *, * ) "     - Finish JobJMult"
        !
    end subroutine JobJMult
    !
    !> Routine to run a full JMult_T job and deliver the result (DSigma) in a text file
    subroutine JobJMult_T()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: dsigma
        !
        !> ????
        type( DataGroupTx_t ), allocatable, dimension(:) :: all_measured_data
        !
        ! Verbose
        write( *, * ) "     - Start JobJMult_T"
        !
        !> Read Model File: instantiates Grid, ModelOperator and ModelParameter
        if( has_model_file ) then 
            !
            call handleModelFile()
            !
        else
            stop "Error: JobJMult_T > Missing Model file!"
        endif
        !
        !> Read Data File: instantiates and builds the Data relation between Txs and Txs
        !> Create DataGroupTxArray for measured data
        if( has_data_file ) then 
            !
            call handleDataFile()
            !
            all_measured_data = all_predicted_data
            !
        else
            stop "Error: JobJMult_T > Missing Data file!"
        endif
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        call createForwardSolver()
        !
        !> Run ForwardModelling for calculate Predicted Data and LRows
        call ForwardModeling( sigma0 )
        !
        !> Get sum of all Tx's dsigma from sensitivity JMult_T
        call JMult_T( sigma0, all_measured_data, dsigma )
        !
        !> Write dsigma to <dsigma_file_name> file path
        call dsigma%write()
        !
        !> Flush local memory
        deallocate( dsigma )
        !
        ! Verbose
        write( *, * ) "     - Finish JobJMult_T"
        !
    end subroutine JobJMult_T
    !
    !> Routine to run a full Inversion Job - Minimize Residual
    !> Where:
    !>    SIGMA (M) = Production model (for predicted data and final inversion model)
    !>    PMODEL = Perturbation model  (if exist -dm readed input model)
    !>    SIGMA0 = Readed input model  (-m)
    !>    DSIGMA = From JMult_T        (????)
    subroutine jobInversionDCG()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: sigma, dsigma, JTP_model
        !
        class( ModelCovarianceRec_t ), allocatable :: model_cov
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: res, all_measured_data, JmHat_data, resJmHat_data, JJTp_data
        type( DataGroupTx_t ), allocatable, dimension(:) :: p_data, r_data, x_data, aux_data, lambdaP
        !
        integer :: iter_dcg, iter_cg, Ndata
        !
        real( kind=prec ) :: alpha, beta, error_cg, residual_rmsd
        !
        real( kind=prec ) :: b_norm, r_norm, r_norm_pre, lambda, SS
        !
        ! Verbose
        write( *, * ) "     - Start jobInversionDCG"
        !
        !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
        if( has_model_file ) then 
            !
            call handleModelFile()
            !
            !> Initialize sigma as the input model
            allocate( sigma, source = sigma0 )
            !
            !> Instantiate ModelCovariance
            allocate( model_cov, source = ModelCovarianceRec_t( sigma ) )
            !
            !> Initialize pmodel with Zeros
            pmodel = sigma0
            !
            call pmodel%zeros()
            !
        else
            stop "Error: jobInversionDCG > Missing Model file!"
        endif
        !
        !> Read Perturbation Model File: instantiate pmodel (NOT USING RIGHT NOW ????)
        if( has_pmodel_file ) then 
            !
            call handlePModelFile()
            !
            call pmodel%setMetric( model_operator%metric )
            !
            dsigma = model_cov%multBy_Cm( pmodel )
            !
            call sigma%add( dsigma )
            !
        endif
        !
        !> Reads Data File: instantiates and builds the Data relation between Txs and Rxs
        if( has_data_file ) then 
            !
            call handleDataFile()
            !
            !> Initialize array with measure data
            all_measured_data = all_predicted_data
            !
        else
            stop "Error: jobInversionDCG > Missing Data file!"
        endif
        !
        !> Instantiate ForwardSolver - Specific type via control file
        call createForwardSolver()
        !
        !> Run ForwardModelling for calculate Predicted Data
        call ForwardModeling( sigma )
        !
        !> Write all_predicted_data, with its proper Rx headers, to the file <predicted_data_file_name>
        call writeDataGroupArray( all_predicted_data, predicted_data_file_name )
        !
        !> Create the array with residual data (all_residual_data)
        !call setResidualData()
        !
        !> Calculate the first rmsd from the residual data
        !residual_rmsd = rmsdDataGroupTxArray( all_residual_data )
        !
            ! initialize res
            res = all_measured_data
            
            call subDataGroupTxArray( res, all_predicted_data )

            ! normalize residuals, compute sum of squares
            all_residual_data=res
            call normalizeDataGroupTxArray( all_residual_data, 2 )
            SS = dotProdDataGroupTxArray( res, all_residual_data )
            Ndata = countDataGroupTxArray( res )

            ! if required, compute the Root Mean Squared misfit
            residual_rmsd = sqrt( SS / Ndata )
        !
        !
        !> Initialize the JmHat data array in the same format as the predicted data array
        JmHat_data = all_predicted_data
        !
        iter_dcg = 1
        !
        lambda = 10.
        !
        ! Verbose
        write( *, * ) "### START DCG LOOP RMSD, LAMBDA: ", residual_rmsd, lambda
        !stop
        !
        !> #### LOOP DCG: Iterations for Inversion Convergence
        do while( iter_dcg < 4 )
            !
            !> Calculate JmHat_data
            call JMult( sigma0, pmodel, JmHat_data )
            !
            !> resJmHat_data (b) = res + JmHat_data
            resJmHat_data = all_residual_data
            !
            call addDataGroupTxArray( resJmHat_data, JmHat_data )
            !
            ! Verbose
            write( *, * ) "### START CG LOOP"
            !
            iter_cg = 1
            error_cg = 1.0
            !
            p_data = resJmHat_data
            r_data = all_residual_data
            x_data = all_measured_data
            !
            !>
            call resetDataGroupTxArray( x_data )
            !
            !> Calculate b_norm
            b_norm = dotProdDataGroupTxArray( resJmHat_data, resJmHat_data )
            !
            !> Calculate r_norm
            r_norm = dotProdDataGroupTxArray( r_data, r_data )
            !
            !> #### LOOP CG_DS_standard ( ARBITRARY VALUES FOR NOW ???? )
            do while( error_cg .GT. 10E-4 .AND. iter_cg .LT. 20 )
                !
                !> Initialize JTP_model, and calculate it from JMult_T
                JTP_model = sigma
                !
                !> Calculate JTP_model
                call JMult_T( sigma0, resJmHat_data, JTP_model )
                !
                !> Smooth JTP_model
                JTP_model = model_cov%multBy_Cm( JTP_model )
                !
                !> Initialize JJTp_data, and calculate it from JMult
                JJTp_data = resJmHat_data
                !
                call JMult( sigma0, JTP_model, JJTp_data )
                !
                !> THIS OPERATIONS DO NOTHING IF ERRORS ARE 1.0
                call normalizeDataGroupTxArray( JJTp_data, 1 )
                !
                !call scMult_dataVectorMTX( lambda, p, lambdaP ) ????
                lambdaP = p_data
                call multDataGroupTxArray( lambdaP, lambda )
                !
                !call linComb_dataVectorMTX( ONE, Ap, ONE, lambdaP, Ap ) ????
                call addDataGroupTxArray( JJTp_data, lambdaP )
                !
                !> Compute alpha: alpha= (r^T r) / (p^T Ap)    
                alpha = r_norm / dotProdDataGroupTxArray( p_data, JJTp_data )
                !
                !> Compute new x: x = x + alpha * p
                call multAddDataGroupTxArray( x_data, p_data, alpha )
                !
                !> Compute new r: r = r - alpha * all_residual_data
                call multAddDataGroupTxArray( r_data, JJTp_data, -alpha )
                !
                !> Calculate r_norm and beta
                r_norm_pre = r_norm
                !
                r_norm = dotProdDataGroupTxArray( r_data, r_data )
                !
                !> Compute beta: beta= r_norm /r_norm_previous
                beta = r_norm / r_norm_pre
                !
                ! Compute new p: p = r + beta * p
                allocate( aux_data, source = r_data ) !Check this
                call multAddDataGroupTxArray( aux_data, p_data, beta )
                p_data = aux_data
                deallocate( aux_data )
                !
                error_cg = r_norm / b_norm 
                !
                ! Verbose
                write( *, * ) "ITER_CG, r_norm, b_norm, error_cg: ", iter_cg, r_norm, b_norm, error_cg
                !
                iter_cg = iter_cg + 1
                !
            enddo
            !
            !> Calculate dsigma from JMult_T
            call JMult_T( sigma0, x_data, dsigma )
            !
            !> Save dsigma in pmodel, to use in the next DCG step
            pmodel = dsigma
            !
            !> Smooth dsigma and add it to sigma
            dsigma = model_cov%multBy_Cm( dsigma )
            !
            call sigma%add( dsigma )
            !
            !> Run ForwardModelling for calculate Predicted Data with new sigma
            call ForwardModeling( sigma )
            !
            !> Create the array with residual data (all_residual_data)
            !call setResidualData()
            !
            !> Residual RMSD 
            !residual_rmsd = rmsdDataGroupTxArray( all_residual_data )
        !
            ! initialize res
            res = all_measured_data
            
            call subDataGroupTxArray( res, all_predicted_data )

            ! normalize residuals, compute sum of squares
            all_residual_data=res
            call normalizeDataGroupTxArray( all_residual_data, 2 )
            SS = dotProdDataGroupTxArray( res, all_residual_data )
            Ndata = countDataGroupTxArray( res )

            ! if required, compute the Root Mean Squared misfit
            residual_rmsd = sqrt( SS / Ndata )
        !
            ! Verbose
            write( *, * ) "          - iter_dcg: ", iter_dcg, residual_rmsd
            !
            iter_dcg = iter_dcg + 1
            !
            !> WRITE DSIGMA FOREACH STEP ????
        enddo
        !
        !> Write dsigma to <dsigma_file_name> file path
        call sigma%write()
        !
        !> Write the data gradient, with its proper Rx headers, to the file <JmHat_data_file_name>
        call writeDataGroupArray( JmHat_data, JmHat_data_file_name )
        !
        !> Flush local memory
        deallocate( sigma, JmHat_data )
        !
        ! Verbose
        write( *, * ) "     - Finish jobInversionDCG"
        !
    end subroutine jobInversionDCG
    !
    !> No subroutine briefing
    subroutine createForwardSolver()
        implicit none
        !
        integer :: i_tx
        !
        !>
        if( allocated( forward_solver ) ) deallocate( forward_solver )
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        select case ( forward_solver_type )
            !
            case( FWD_IT_DC )
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
            case default
                !
                write( *, * ) "Warning: createForwardSolver > Undefined forward_solver, using IT_DC"
                !
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
        end select
        !
        !> Loop over all Transmitters
        do i_tx = 1, size( transmitters )
            !
            transmitters( i_tx )%Tx%forward_solver => forward_solver
            !
        enddo
        !
    end subroutine createForwardSolver
    !
    !> Allocate and Initialize the predicted data array, 
    !> according to the arrangement of the Transmitter-Receiver pairs of the input
    subroutine initPredictedDataArray()
        implicit none
        !
        !> Auxiliary variable to group data under a single transmitter index
        type( DataGroupTx_t ) :: tx_data
        !
        !> Local indexes
        integer :: i_tx, i_data
        !
        !> If is the first time, create the predicted data array, 
        !> according to the arrangement of the Transmitter-Receiver pairs of the input
        if( .NOT. allocated( all_predicted_data ) ) then
            !
            write( *, * ) "     - Create predicted data array"
            !
            !> Create an array of DataGroupTx to store the predicted data in the same format as the measured data
            !> Enabling the use of grouped predicted data in future jobs (all_predicted_data)
            do i_tx = 1, size( transmitters )
                !
                tx_data = DataGroupTx_t( i_tx )
                !
                do i_data = 1, size( measured_data )
                    !
                    if( measured_data( i_data )%i_tx == i_tx ) then
                        !
                        call tx_data%put( measured_data( i_data ) )
                        !
                    endif
                enddo
                !
                call updateDataGroupTxArray( all_predicted_data, tx_data )
                !
            enddo
            !
        else
            stop "Error: initPredictedDataArray > Unnecessary recreation of predicted data array"
        endif
        !
    end subroutine initPredictedDataArray
    !
    !> No subroutine briefing
    subroutine handleJob()
        implicit none
        !
        select case ( modem_job )
            !
            case ( "forward" )
                !
                call jobForwardModeling()
                !
            case ( "JMult" )
                !
                call JobJMult()
                !
            case ( "JMult_T" )
                !
                call JobJMult_T()
                !
            case ( "inversion" )
                !
                call jobInversionDCG()
                !
            case default
                !
                write( *, * ) "Error: unknown job: [", modem_job, "]"
                call printHelp()
                stop
            !
        end select
        !
    end subroutine handleJob
    !
    !> No subroutine briefing
    subroutine handleControlFile()
        implicit none
        !
        type( ModEMControlFile_t ) :: control_file
        !
        write( *, * ) "     < Control File: [", control_file_name, "]"
        !
        !> Instantiate the ControlFile object
        !> Reads control file and sets the options in the Constants module
        control_file = ModEMControlFile_t( ioStartup, control_file_name )
        !
    end subroutine handleControlFile
    !
    !> No subroutine briefing
    subroutine handleModelFile()
        implicit none
        !
        type( ModelReader_Weerachai_t ) :: model_reader
        type( TAirLayers ) :: air_layer
        !
        integer :: i
        !
        ! Verbose
        write( *, * ) "     < Model File: [", model_file_name, "]"
        !
        !> Read Grid and ModelParameter with ModelReader_Weerachai
        call model_reader%Read( model_file_name, main_grid, sigma0 ) 
        !
        !> Instantiate the ModelOperator object
        select type( main_grid )
            !
            class is( Grid3D_SG_t )
                !
                call main_grid%SetupAirLayers( air_layer, model_method, model_n_air_layer, model_max_height )
                !
                call main_grid%UpdateAirLayers( air_layer%nz, air_layer%dz )
                !
                write( *, * ) "          Air Layers [i, dz(i)]:"
                !
                do i = air_layer%nz, 1, -(1)
                    write( *, * ) "               ", i, air_layer%dz(i)
                enddo
                !
                write( *, * ) "          Air layers from the method: ", air_layer%method, "."
                !
                write( *, * ) "          Top of the air layers is at ", sum( air_layer%Dz ) / 1000, " km."
                !
                write( *, "(A, F10.2, A2, F10.2, A2, F10.2, A4, F10.2)" ) "           o(x_data,y,z) * rotDeg: (", main_grid%ox, ", ", main_grid%oy, ", ", main_grid%oz, ") * ", main_grid%rotDeg
                !
                allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
                !
                call model_operator%SetEquations()
                !
                call sigma0%setMetric( model_operator%metric )
                !
                call model_operator%SetCond( sigma0 )
                !
            class default
                stop "Error: handleModelFile > Unclassified main_grid"
            !
        end select
        !
    end subroutine handleModelFile
    !
    !> No subroutine briefing
    subroutine handlePModelFile()
        implicit none
        !
        type( ModelReader_Weerachai_t ) :: model_reader
        !
        class( Grid_t ), allocatable :: prior_grid
        !
        ! Verbose
        write( *, * ) "     < Prior Model File: [", pmodel_file_name, "]"
        !
        !> Read prior_grid and pmodel with ModelReader_Weerachai
        call model_reader%Read( pmodel_file_name, prior_grid, pmodel ) 
        !
        deallocate( prior_grid )
        !
    end subroutine handlePModelFile
    !
    !> No subroutine briefing
    subroutine handleDataFile()
        implicit none
        !
        integer :: irx, nrx
        class( Receiver_t ), pointer :: Rx
        !
        !> Local object to dealt data, self-destructs at the end of the subroutine
        type( DataFileStandard_t ) :: data_file_standard
        !
        write( *, * ) "     < Data File: [", data_file_name, "]"
        !
        data_file_standard = DataFileStandard_t( ioStartup, data_file_name )
        !
        nrx = size( receivers )
        !
        ! Verbose
        if( nrx == data_file_standard%nRx ) then
            !
            write( *, * ) "          Checked ", nrx, " Receivers."
            !
        else
            !
            write( *, * ) "Error: Number of Rx mismatched from Header :[", nrx, " and ", data_file_standard%nRx, "]"
            stop
            !
        endif
        !
        if( size( transmitters ) == data_file_standard%nTx ) then
             !
             call printTransmitterArray()
             !
        else
             !
             write( *, * ) "Number of Tx mismatched from Header :[", size( transmitters ), " and ", data_file_standard%nTx, "]"
             stop
             !
        endif
        !
        write( *, * ) "          Checked ", size( measured_data ), " DataGroups."
        !
        !> Initialize the predicted data array
        call initPredictedDataArray()
        !
        write( *, * ) "     - Create Rx evaluation vectors"
        !
        !> Calculate evaluation vectors for the Receivers
        do irx = 1, nrx
            !
            Rx => getReceiver( irx )
            !
            call Rx%evaluationFunction( model_operator )
            !
        enddo
        !
    end subroutine handleDataFile
    !
    !> No subroutine briefing
    subroutine handleArguments()
        implicit none
        !
        character( len=200 ) :: argument
        integer :: argument_index
        !
        if( command_argument_count() == 0 ) then
            !
            call printHelp()
            call printUsage()
            stop
            !
        else
            !
            argument_index = 1
            !
            do while( argument_index <= command_argument_count() ) 
                 !
                 call get_command_argument( argument_index, argument )
                 !
                 select case ( argument )
                      !
                      case ( "-c", "--control" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         control_file_name = trim( argument )
                         !
                         if( len( control_file_name ) > 0 ) has_control_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-d", "--data" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         data_file_name = trim( argument )
                         !
                         if( len( data_file_name ) > 0 ) has_data_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-f", "--forward" )
                         !
                         modem_job = "forward"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-i", "--inversion" )
                         !
                         modem_job = "inversion"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-j", "--jmult" )
                         !
                         modem_job = "JMult"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-jt", "--jmult_t" )
                         !
                         modem_job = "JMult_T"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-m", "--model" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         model_file_name = trim( argument )
                         !
                         if( len( model_file_name ) > 0 ) has_model_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-pm", "--pmodel" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         pmodel_file_name = trim( argument )
                         !
                         if( len( pmodel_file_name ) > 0 ) has_pmodel_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-dm", "--dsigma" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         dsigma_file_name = trim( argument )
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-pd", "--predicted" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         predicted_data_file_name = trim( argument )
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-gd", "--gradient" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         JmHat_data_file_name = trim( argument )
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-es", "--esolution" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         e_solution_file_name = trim( argument )
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-v", "--version" )
                         !
                         write( *, * ) "    + ModEM-OO version 1.0.0"
                         stop
                         !
                      case ( "-h", "--help" )
                         !
                         call printHelp()
                         call printUsage()
                         stop
                         !
                      case ( "-vb", "--verbose" )
                         !
                         stop "Error: handleArguments > Verbose level not implemented yet!"
                         !
                      case default
                         !
                         write( *, * ) "Error: Unknown Argument: [", trim( argument ), "]"
                         call printHelp()
                         stop
                         !
                 end select
                 !
                 argument = ""
                 !
            enddo
            !
        endif
        !
    end subroutine handleArguments
    !
    !> No subroutine briefing
    subroutine setupDefaultParameters()
        implicit none
        !
        ! I|O
        predicted_data_file_name = "predicted_data.dat"
        JmHat_data_file_name     = "JmHat.dat"
        e_solution_file_name     = "esolution.bin"
        dsigma_file_name         = "dsigma.mod"
        has_control_file         = .FALSE.
        has_model_file           = .FALSE.
        has_pmodel_file          = .FALSE.
        has_data_file            = .FALSE.
        verbosis                 = .FALSE.
        !
        ! Solvers
        QMR_iters = 40
        BCG_iters = 80
        max_divcor = 20
        max_divcor_iters = 100
        tolerance_divcor = 1E-5
        tolerance_qmr = 1E-7
        forward_solver_type = FWD_IT_DC
        !
        ! Source
        source_type = SRC_MT_1D
        get_1D_from = "Geometric_mean"
        !
        ! Model
        model_method      = MM_METHOD_FIXED_H
        model_n_air_layer = 10
        model_max_height  = 200.0
        !
    end subroutine setupDefaultParameters
    !
    !> No subroutine briefing
    subroutine garbageCollector()
        implicit none
        !
        !> Flush memory used by job data arrays
        deallocate( measured_data )
        !
        !> Flush memory used by job data arrays
        deallocate( all_predicted_data )
        !
        !> Deallocate global array of Receivers
        call deallocateReceiverArray()
        !
        !> Deallocate global array of Transmitters
        call deallocateTransmitterArray()
        !
        !> Deallocate global components
        deallocate( forward_solver, model_operator, sigma0, main_grid )
        !
        !> Deallocate global pmodel, if its the case
        if( allocated( pmodel ) ) deallocate( pmodel )
        !
        !> Flush memory used by main program control variables and flags
        if( allocated( forward_solver_type ) ) deallocate( forward_solver_type )
        if( allocated( source_type ) ) deallocate( source_type )
        if( allocated( model_method ) ) deallocate( model_method )
        if( allocated( get_1D_from ) ) deallocate( get_1D_from )
        if( allocated( predicted_data_file_name ) ) deallocate( predicted_data_file_name )
        if( allocated( JmHat_data_file_name ) ) deallocate( JmHat_data_file_name )
        if( allocated( e_solution_file_name ) ) deallocate( e_solution_file_name )
        !
        if( allocated( control_file_name ) ) deallocate( control_file_name )
        if( allocated( model_file_name ) ) deallocate( model_file_name )
        if( allocated( pmodel_file_name ) ) deallocate( pmodel_file_name )
        if( allocated( data_file_name ) ) deallocate( data_file_name )
        if( allocated( dsigma_file_name ) ) deallocate( dsigma_file_name )
        if( allocated( modem_job ) ) deallocate( modem_job )
        !
    end subroutine garbageCollector
    !
    !> No subroutine briefing
    subroutine writeDataGroupArray( target_tx_data_array, file_name )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: target_tx_data_array
        character(*), intent( in ) :: file_name
        !
        type( DataGroup_t ), allocatable, dimension(:) :: to_write_data
        !
        class( Transmitter_t ), pointer :: transmitter
        class( Receiver_t ), pointer :: receiver
        type( DataGroup_t ), pointer :: data_group
        !
        integer :: receiver_type, i, j, ios
        !
        ! Verbose
        write( *, * ) "     > Write Data to file: [", file_name, "]"
        !
        !> Put the data in the same input format
        allocate( to_write_data, source = measured_data )
        !
        do i = 1, size( target_tx_data_array )
            do j = 1, size( target_tx_data_array(i)%data )
                call setDataGroup( to_write_data, target_tx_data_array(i)%data(j) )
            enddo
        enddo
        !
        receiver_type = 0
        !
        open( unit = ioPredData, file = file_name, action = "write", form = "formatted", iostat = ios )
        !
        if( ios == 0 ) then
            !
            do i = 1, size( to_write_data )
                !
                data_group => getDataGroupByIndex( to_write_data, i )
                !
                receiver => getReceiver( data_group%i_rx )
                !
                call writePredictedDataHeader( receiver, receiver_type )
                !
                transmitter => getTransmitter( data_group%i_tx )
                !
                do j = 1, data_group%n_comp
                    !
                    select type( transmitter )
                        !
                        class is( TransmitterMT_t )
                            !
                            write( ioPredData, "(es12.6, 1X, A, 1X, f15.3, f15.3, f15.3, f15.3, f15.3, 1X, A, 1X, es16.6, es16.6, es16.6)" ) transmitter%period, trim(receiver%code), R_ZERO, R_ZERO, receiver%location(1), receiver%location(2), receiver%location(3), trim(data_group%components(j)%str), data_group%reals(j), data_group%imaginaries(j), data_group%errors(j)
                            !
                        class is( TransmitterCSEM_t )
                            !
                            write( ioPredData, "(A, 1X, es12.6, f15.3, f15.3, f15.3, f15.3, f15.3, f15.3, 1X, A, 1X, f15.3, f15.3, f15.3, 1X, A, 1X, es16.6, es16.6, es16.6)" ) trim(transmitter%dipole), transmitter%period, transmitter%moment, transmitter%azimuth, transmitter%dip, transmitter%location(1), transmitter%location(2), transmitter%location(3), trim(receiver%code), receiver%location(1), receiver%location(2), receiver%location(3), trim(data_group%components(j)%str), data_group%reals(j), data_group%imaginaries(j), data_group%errors(j)
                            !
                        class default
                            stop "Error: Unclassified data_group!"
                        !
                    end select
                    !
                enddo
                !
            enddo
            !
            deallocate( to_write_data )
            !
            close( ioPredData )
            !
        else
            write( *, * ) "Error opening [", file_name, "] in writeDataGroupArray!"
            stop
        endif
        !
    end subroutine writeDataGroupArray
    !
    !> No subroutine briefing
    subroutine writePredictedDataHeader( receiver, receiver_type )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: receiver
        !
        integer, intent( inout ) :: receiver_type
        !
        if( receiver_type /= receiver%rx_type ) then
            !
            select case( receiver%rx_type )
                !
                case( 1, 11, 12 )
                    !
                    write( ioPredData, "(4A, 40A)" ) "#    ", DATA_FILE_TITLE_MT
                    !
                    write( ioPredData, "(74A)" ) "#    Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
                    !
                case( 2 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Full_Interstation_TF"
                case( 3 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Off_Diagonal_Rho_Phase"
                case( 4 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Phase_Tensor"
                case( 5 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Off_Diagonal_Impedance"
                case( 6, 7, 8, 9, 10 )
                    !
                    write( ioPredData, "(4A, 40A)" ) "#    ", DATA_FILE_TITLE_CSEM
                    !
                    write( ioPredData, "(125A)" ) "#    Tx_Dipole Tx_Period(s) Tx_Moment(Am) Tx_Azi Tx_Dip Tx_X(m) Tx_Y(m) Tx_Z(m) Code X(m) Y(m) Z(m) Component Real Imag Error"
                    !
                case default
                    write( *, * ) "Unknown receiver type :[", receiver%rx_type, "]"
                    stop "Error: test_FWD.f90: writePredictedDataHeader()"
                !
            end select
            !
            write( ioPredData, "(4A, 100A)" ) ">    ", getStringReceiverType( receiver%rx_type )
            write( ioPredData, "(4A, 100A)" ) ">    ", "exp(-i\omega t)"
            write( ioPredData, "(4A, 100A)" ) ">    ", "[V/m]/[T]"
            write( ioPredData, "(7A, 100A)" ) ">        ", "0.00"
            write( ioPredData, "(7A, 100A)" ) ">        ", "0.000    0.000"
            write( ioPredData, "(A3, i8, i8)" ) ">        ", size( transmitters ), size( receivers )
            !
            receiver_type = receiver%rx_type
            !
        endif
        !
    end subroutine writePredictedDataHeader
    !
    !> No subroutine briefing
    subroutine writeEsolutionHeader( nTx, nMode )
        implicit none
        !
        integer, intent( in ) :: nTx, nMode
        integer :: ios
        character(len=20) :: version
        !
        version = "Modem-OO"
        !
        open( ioESolution, file = e_solution_file_name, action = "write", form = "unformatted", iostat = ios)
        !
        if( ios == 0 ) then
            !
            write( ioESolution ) version, nTx, nMode, &
            main_grid%nx, main_grid%ny, main_grid%nz, main_grid%nzAir, &
            main_grid%ox, main_grid%oy, main_grid%oz, main_grid%rotdeg
            !
            write( ioESolution ) main_grid%dx
            write( ioESolution ) main_grid%dy
            write( ioESolution ) main_grid%dz
            !
            close( ioESolution )
            !
        else
            !
            write( *, * ) "Error opening file in writeEsolutionHeader [", e_solution_file_name, "]!"
            stop
            !
        endif
        !
        !
    end subroutine writeEsolutionHeader
    !
    !> No subroutine briefing
    subroutine printUsage()
        implicit none
        !
        write( *, * ) "ModEM Minimal Usage:"
        write( *, * ) ""
        write( *, * ) "    Forward Modeling:"
        write( *, * ) "        <ModEM> -f -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Outputs:"
        write( *, * ) "        - predicted_data.dat or the path specified by [-pd]"
        write( *, * ) "        - esolution.bin or the path specified by      [-es]"
        write( *, * ) ""
        write( *, * ) "    JMult:"
        write( *, * ) "        <ModEM> -j -m <rFile_Model> -pm <rFile_pModel> -d <rFile_Data>"
        write( *, * ) "    Outputs:"
        write( *, * ) "        - JmHat.dat or the path specified by          [-jm]"
        write( *, * ) "        - predicted_data.dat or the path specified by [-pd]"
        write( *, * ) "        - esolution.bin or the path specified by      [-es]"
        write( *, * ) ""
        write( *, * ) "    JMult_T:"
        write( *, * ) "        <ModEM> -jt -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Output:"
        write( *, * ) "        - dsigma.mod or the path specified by         [-dm]"
        write( *, * ) "        - predicted_data.dat or the path specified by [-pd]"
        write( *, * ) "        - esolution.bin or the path specified by      [-es]"
        write( *, * ) ""
        write( *, * ) "    Inversion:"
        write( *, * ) "        <ModEM> -i -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Output:"
        write( *, * ) "        - dsigma.mod or the path specified by         [-dm]"
        !
    end subroutine printUsage
    !
    !> No subroutine briefing
    subroutine printHelp()
        implicit none
        !
        write( *, * ) "ModEM Options:"
        write( *, * ) ""
        write( *, * ) "    Flags to define a job:"
        write( *, * ) "        [-f],  [--forward]   :  Forward Modeling."
        write( *, * ) "        [-j],  [--jmult]     :  JMult."
        write( *, * ) "        [-jt], [--jmult_t]   :  Transposed JMult."
        write( *, * ) "        [-i],  [--inversion] :  Inversion."
        write( *, * )
        write( *, * ) "    Other arguments:"
        write( *, * ) "        [-d],  [--data]      :  Flag for the input data file path."
        write( *, * ) "        [-m],  [--model]     :  Flag for the input model file path."
        write( *, * ) "        [-pm], [--pmodel]    :  Flag for the input perturbation model file path."
        write( *, * ) "        [-c],  [--control]   :  Flag for the user control file path."
        write( *, * ) "        [-dm], [--dmodel]    :  Flag for the output dsigma model file path."
        write( *, * ) "        [-pd], [--predicted] :  Flag for the output predicted data file path."
        write( *, * ) "        [-jm], [--jmhat]     :  Flag for the output JmHat data file path."
        write( *, * ) "        [-es], [--esolution] :  Flag for the binary output e-solution file path."
        write( *, * ) "        [-v],  [--version]   :  Print version."
        write( *, * ) "        [-vb], [--verbose]   :  Print runtime information."
        write( *, * ) "        [-h],  [--help]      :  Print usage information."
        !
        write( *, * ) ""
        write( *, * ) "Version 1.0.0"
        !
    end subroutine printHelp
    !
end program TestSerial
!