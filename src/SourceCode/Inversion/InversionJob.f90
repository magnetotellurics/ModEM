!
!> Module with a single job ????
!
module InversionJob
    !
    use InversionDCG
    use InversionNLCG
    !
    !> Public module routines
    public :: jobInversion
    !
contains
    !
    !> Routine to run a full Inversion Job - Minimize Residual data error
    !> Where:
    !>     sigma  = Input|Start model (for predicted data and final inversion model)
    !>     pmodel = Perturbation model (if exist -dm read input model)
    !>     sigma0 = Read input model (-m)
    !>     dsigma = Production model (data gradient From serialJMult_T)
    !
    subroutine jobInversion()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: sigma, dsigma
        !
        class( Inversion_t ), allocatable :: inversion
        !
        ! Verbose
        write( *, * ) "     - Start jobInversion"
        !
        !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
        if( has_model_file ) then 
            !
            call handleModelFile( sigma )
            !
            !> Instantiate ModelCovariance
            allocate( model_cov, source = ModelCovarianceRec_t( sigma ) )
            !
            if( has_cov_file ) then 
                !
                call model_cov%read_CmSqrt( cov_file_name )
                !
                write( *, * ) "     < Cov File: [", cov_file_name, "]"
                !
            else
				call warning( "jobInversion > Missing Covariance file!" )
                !write( *, * ) achar(27)//"[91m# Warning:"//achar(27)//"[0m jobInversion > Missing Covariance file!"
            endif
            !
            !> Initialize pmodel with Zeros
            allocate( dsigma, source = sigma )
            !
            call dsigma%zeros
            !
        else
            stop "Error: jobInversion > Missing Model file!"
        endif
        !
        !> Read Perturbation Model File: instantiate pmodel (NOT USING RIGHT NOW ????)
        if( has_pmodel_file ) then 
            !
            deallocate( dsigma )
            !
            call handlePModelFile( dsigma )
            !
            call dsigma%setMetric( model_operator%metric )
            !
            call model_cov%multBy_Cm( dsigma )
            !
            call sigma%linComb( ONE, ONE, dsigma )
            !
        endif
        !
        !> Read Data File: instantiate and build the Data relation between Txs and Rxs
        if( has_data_file ) then 
            !
            call handleDataFile
            !
        else
            stop "Error: jobInversion > Missing Data file!"
        endif
        !
#ifdef MPI
        !
        call broadcastBasicComponents
        !
#else
        !
        call createDistributeForwardSolver
        !
#endif
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        select case( inversion_type )
            !
            case( DCG )
                !
                allocate( inversion, source = InversionDCG_t() )
                !
            case( NLCG )
                !
                allocate( inversion, source = InversionNLCG_t() )
                !
            case default
                !
                stop "Error: jobInversion > Undefined inversion_type"
                !
        end select
        !
        call inversion%solve( all_measured_data, sigma, dsigma )
        !
#ifdef MPI
        call broadcastFinish
#endif
        !
        deallocate( sigma, dsigma, inversion )
        !
        ! Verbose
        write( *, * ) "     - Finish jobInversion"
        !
    end subroutine jobInversion
    !
end module InversionJob
!