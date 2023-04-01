!
!> Derived class to define a InversionNLCG
!
module InversionNLCG
    !
    use Inversion
    !
    type, extends( Inversion_t ) :: InversionNLCG_t
        !
        ! the condition to identify when the inversion stalls
        real( kind=prec ) :: fdiffTol
        ! exit if lambda < lambda_tol approx. 1e-4
        real( kind=prec ) :: lambda_tol
        ! set lambda_i = lambda_{i-1}/k when the inversion stalls
        real( kind=prec ) :: k
        ! the factor that ensures sufficient decrease in the line search
        real( kind=prec ) :: c
        ! the factor that ensures curvature condition in the line search
        real( kind=prec ) :: c2
        ! restart CG every nCGmax iterations to ensure conjugacy
        integer :: nCGmax
        ! the starting step for the line search
        real( kind=prec ) :: alpha_1
        ! maximum initial delta mHat (overrides alpha_1)
        real( kind=prec ) :: startdm
        ! optional relaxation parameter (Renormalized Steepest Descent algorithm)
        real( kind=prec ) :: gamma
        !
        contains
            !
            final :: InversionNLCG_dtor
            !
            procedure, public :: solve => solveInversionNLCG
            !
            procedure, private :: gradient, func, updateDampingParameter, lineSearchCubic
            !
    end type InversionNLCG_t
    !
    private :: weightGradrients, writeHeaders, cdInvMult, outputFilesInversionNLCG
    !
    interface InversionNLCG_t
        module procedure InversionNLCG_ctor
    end interface InversionNLCG_t
    !
contains
    !
    !> No function briefing
    !
    function InversionNLCG_ctor() result( self )
        implicit none
        !
        type( InversionNLCG_t ) :: self
        !
        !write( *, * ) "Constructor InversionNLCG_t"
        !
        call self%init
        !
        !> NLCG Default Parameters
        ! the factor that ensures sufficient decrease in the line search >=1e-4
        self%c = 1.0e-4
        ! restart CG every nCGmax iterations to ensure conjugacy
        self%nCGmax = 8
        ! the starting step for the line search
        self%alpha_1 = 20.
        ! optional relaxation parameter (Renormalized Steepest Descent algorithm)
        self%gamma = 0.99
        !
        !> NLCG Parameters that can be defined via control file
        ! exit if lambda < lambdaTol approx. 1e-8
        self%lambda_tol = 1.0e-4
        ! set lambda_i = lambda_{i-1}/k when the inversion stalls
        self%k = 10.    !lambda_div
        ! maximum initial delta mHat (overrides alpha_1)
        self%startdm = 20.
        ! inversion stalls when abs(rms - rmsPrev) < fdiffTol (2e-3 works well)
        self%fdiffTol = 2.0e-3
        !
        if( has_inv_control_file ) then
            !
            if( allocated( inv_control_file%lambda_tol ) ) &
                read( inv_control_file%lambda_tol, * ) self%lambda_tol
            !
            if( allocated( inv_control_file%lambda_div ) ) &
                read( inv_control_file%lambda_div, * ) self%k
            !
            if( allocated( inv_control_file%startdm ) ) &
                read( inv_control_file%startdm, * ) self%startdm
            !
            if( allocated( inv_control_file%fdiffTol ) ) &
                read( inv_control_file%fdiffTol, * ) self%fdiffTol
            !
        endif
        !
        write( *, "( A45, es20.2 )" ) "lambda_tol = ", self%lambda_tol
        !
        write( *, "( A45, es20.2 )" ) "lambda_div = ", self%k
        !
        write( *, "( A45, es20.2 )" ) "startdm = ", self%startdm
        !
        write( *, "( A45, es20.2 )" ) "fdiffTol = ", self%fdiffTol
        !
        !> Free the memory used by the global control file, which is no longer useful
        if( allocated( inv_control_file ) ) deallocate( inv_control_file )
        !
    end function InversionNLCG_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    !
    subroutine InversionNLCG_dtor( self )
        implicit none
        !
        type( InversionNLCG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor InversionNLCG_t"
        !
    end subroutine InversionNLCG_dtor
    !
    !> No subroutine briefing
    !
    subroutine solveInversionNLCG( self, all_data, sigma, dsigma )
        implicit none
        !
        class( InversionNLCG_t ), intent( inout ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_data
        class( ModelParameter_t ), allocatable, intent( in ) :: sigma
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        !> initial step size in the line search direction in model units
        real( kind=prec ) :: startdm
        !> flavor is a string that specifies the algorithm to use
        character(80) :: flavor = "Cubic"
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: dHat, res
        class( ModelParameter_t ), allocatable :: mHat, grad, g, h, gPrev
        real( kind=prec ) :: r_value, valuePrev, rms
        real( kind=prec ) :: rmsPrev, alpha, beta
        real( kind=prec ) :: gnorm, mNorm, Nmodel
        real( kind=prec ) :: grad_dot_h, g_dot_g
        real( kind=prec ) :: g_dot_gPrev, g_dot_h
        real( kind=prec ) :: gPrev_dot_gPrev, h_dot_g, h_dot_gPrev
        integer :: iter, nCG, nLS, nfunc, ios, i_sol
        !
        !>
        call createOutputDirectory()
        !
        call writeHeaders
        !
        ! Verbose
        write( *, * ) "     - Start Inversion NLCG, output files in [", trim( outdir_name ), "]"
        !
        ! initialize the line search
        alpha = self%alpha_1
        !
        startdm = self%startdm
        !
        !>
        open( unit = ioInvLog, file = trim( outdir_name )//"/NLCG.log", status="unknown", position="append", iostat=ios )
        !
        if( ios == 0 ) then
            !
            write( *, "( a50, es12.5 )" ) "The initial damping parameter lambda is ", self%lambda
            write( *, "( a64, f12.6 )" ) "The initial line search step size (in model units) is ", startdm
            !
            write( ioInvLog, "( a41, es8.1 )" ) "The initial damping parameter lambda is ", self%lambda
            write( ioInvLog, "( a55, f12.6 )" ) "The initial line search step size (in model units) is ", startdm
            !
            ! starting model contains the rough deviations from the prior
            allocate( mHat, source = dsigma )
            !
            !  compute the penalty functional and predicted data
            i_sol = 0
            !
            call func( self, all_data, sigma, mHat, r_value, mNorm, dHat, i_sol, rms )
            !
            nfunc = 1
            !
            ! output (smoothed) initial model and responses for later reference
            call model_cov%multBy_CmSqrt( mHat, dsigma )
            !
            call dsigma%linComb( ONE, ONE, sigma )
            !
            !> compute gradient of the full penalty functional
            call gradient( self, all_data, sigma, mHat, grad, dHat, i_sol )
            !
            gnorm = sqrt( grad%dotProd( grad ) )
            !
            write( *, "( a42, es12.5 )" ) "    GRAD: initial norm of the gradient is", gnorm
            write( ioInvLog, "( a42, es12.5 )" ) "     GRAD: initial norm of the gradient is", gnorm
            !
            if( gnorm < TOL6 ) then
                stop "Error: NLCGsolver: Problem with your gradient computations: first gradient is zero"
            else
                !
                alpha = startdm / gnorm
                !
                write( *, "( a39, es12.5 )" ) "The initial value of alpha updated to ", alpha
                write( ioInvLog, "( a39, es12.5 )" ) "The initial value of alpha updated to ", alpha
                !
            endif
            !
            !> initialize CG: g = - grad; h = g
            nCG = 0
            iter = 0
            !
            allocate( g, source = grad )
            !
            call g%linComb( MinusONE, R_ZERO, grad )
            !
            allocate( h, source = g )
            !
            do! while( rms .GE. self%rms_tol .AND. iter .LT. self%max_inv_iters )
                !
                !  test for convergence ...
                if( rms .LT. self%rms_tol .OR. iter .GE. self%max_inv_iters ) then
                    exit
                endif
                !
                iter = iter + 1
                !
                ! save the values of the functional and the directional derivative
                rmsPrev = rms
                !
                valuePrev = r_value
                !
                grad_dot_h = grad%dotProd( h )
                !
                ! at the end of line search, set mHat to the new r_value
                ! mHat = mHat + alpha*h  and evaluate gradient at new mHat
                ! data and solnVector only needed for output
                !
                write( *, "(a23)" ) "Starting line search..."
                write( ioInvLog, "(a23)" ) "Starting line search..."
                !
                select case( flavor )
                    !
                    case( "Cubic" )
                        !
                        call self%lineSearchCubic( all_data, sigma, h, alpha, mHat, r_value, grad, rms, nLS, dHat )
                        !
                    case( "Quadratic" )
                        !
                        !call self%lineSearchQuadratic( all_data,sigma,h,alpha,mHat,r_value,grad,rms,nLS,dHat )
                        !
                    case( "Wolfe" )
                        !
                        !call self%lineSearchWolfe( all_data,sigma,h,alpha,mHat,r_value,grad,rms,nLS,dHat )
                        !
                    case default
                        stop "Error: NLCGsolver: Unknown line search requested in NLCG"
                    !
                end select
                !
                nfunc = nfunc + nLS
                !
                if( allocated( gPrev ) ) deallocate( gPrev )
                allocate( gPrev, source = g )
                !
                if( allocated( g ) ) deallocate( g )
                allocate( g, source = grad )
                !
                call g%linComb( MinusONE, R_ZERO, grad )
                !
                ! compute the starting step for the next line search
                alpha = 2 * ( r_value - valuePrev ) / grad_dot_h
                !
                ! adjust the starting step to ensure super linear convergence properties
                alpha = ( ONE + 0.01 ) * alpha
                !
                write( *, "( a25, i5 )" ) "Completed NLCG iteration ", iter
                write( ioInvLog, "( a25, i5 )" ) "Completed NLCG iteration ", iter
                ! 
                Nmodel = mHat%countModel()
                !
                mNorm = mHat%dotProd( mHat ) / Nmodel
                !
                write( *, * ) "     lambda, alpha, r_value, mNorm, rms: ", self%lambda, alpha, r_value, mNorm, rms
                write( ioInvLog, * ) "     lambda, alpha, r_value, mNorm, rms: ", self%lambda, alpha, r_value, mNorm, rms
                !
                ! write out the intermediate model solution and responses
                call model_cov%multBy_CmSqrt( mHat, dsigma )
                !
                call dsigma%linComb( ONE, ONE, sigma )
                !
                res = all_data
                !
                call linCombData( ONE, all_data, MinusONE, dHat, res )
                !
                call outputFilesInversionNLCG( iter, dHat, res, dsigma, mHat )
                !
                !> if alpha is too small, we are not making progress: update lambda
                if( abs( rmsPrev - rms ) < self%fdiffTol ) then
                    !
                    ! update lambda, penalty functional and gradient
                    call self%updateDampingParameter( mHat, r_value, grad )
                    !
                    ! check that lambda is still at a reasonable r_value
                    if( self%lambda < self%lambda_tol ) then
                        write( *, "(a55)" ) "Unable to get out of a local minimum. Exiting..."
                        write( ioInvLog, "(a55)" ) "Unable to get out of a local minimum. Exiting..."
                        exit
                    endif
                    !
                    !> update alpha
                    gnorm = sqrt( grad%dotProd( grad ) )
                    !
                    write( *, "(a34,es12.5)" ) "The norm of the last gradient is ", gnorm
                    write( ioInvLog, "(a34,es12.5)" ) "The norm of the last gradient is ", gnorm
                    !
                    !> alpha = min(self%alpha_1,startdm/gnorm)
                    alpha = min( ONE, startdm ) / gnorm
                    !
                    write( *, "( a48, es12.5 )" ) "The value of line search step alpha updated to ", alpha
                    write( ioInvLog, "( a48, es12.5 )" ) "The value of line search step alpha updated to ", alpha
                    !
                    !> g = - grad
                    g = grad
                    !
                    call g%linComb( MinusONE, R_ZERO, grad )
                    !
                    !> restart
                    write( *, * ) "Restarting NLCG with the damping parameter updated"
                    write( *, * ) "lambda, alpha, r_value, mNorm, rms: ", self%lambda, alpha, r_value, mNorm, rms
                    !
                    write( ioInvLog, * ) "Restarting NLCG with the damping parameter updated"
                    write( ioInvLog, * ) "lambda, alpha, r_value, mNorm, rms: ", self%lambda, alpha, r_value, mNorm, rms
                    !
                    if( allocated( h ) ) deallocate( h )
                    allocate( h, source = g )
                    !
                    nCG = 0
                    !
                    cycle
                    !
                endif
                !
                g_dot_g = g%dotProd( g )
                !
                g_dot_gPrev = g%dotProd( gPrev )
                !
                gPrev_dot_gPrev = gPrev%dotProd( gPrev )
                !
                g_dot_h = g%dotProd( h )
                !
                !> Polak-Ribiere variant
                beta = ( g_dot_g - g_dot_gPrev ) / gPrev_dot_gPrev
                !
                !> restart CG if the orthogonality conditions fail. Using the fact that
                !> h_{i+1} = g_{i+1} + beta * h_i. In order for the next directional
                !> derivative = -g_{i+1}.dot.h_{i+1} to be negative, the condition
                !> g_{i+1}.dot.(g_{i+1}+beta*h_i) > 0 must hold. Alternatively, books
                !> say we can take beta > 0 (didn"t work as well)
                !> if((beta.lt.R_ZERO).or.(g_dot_g + beta*g_dot_h .le. R_ZERO)&
                !>    .AND.(nCG .ge. self%nCGmax)) then  !PR+
                if( g_dot_g + beta * g_dot_h .LE. R_ZERO .AND. nCG .GE. self%nCGmax ) then  !PR
                    !
                    ! restart
                    write( *, "(a45)" ) "Restarting NLCG to restore orthogonality"
                    write( ioInvLog, "(a45)" ) "Restarting NLCG to restore orthogonality"
                    !
                    nCG = 0
                    !
                    beta = R_ZERO
                    !
                else
                    !
                    nCG = nCG + 1
                    !
                endif
                !
                if( allocated( h ) ) deallocate( h )
                allocate( h, source = g )
                !
                call h%linComb( ONE, beta, h )
                !
            enddo
            !
            !> multiply by C^{1/2} and add m_0
            call model_cov%multBy_CmSqrt( mHat, dsigma )
            !
            call dsigma%linComb( ONE, ONE, sigma )
            !
            all_data = dHat
            !
            ! cleaning up
            !call deallocateDataGroupTxArray( dHat )
            !call deallocateDataGroupTxArray( res )
            !
            write( *, "( a25, i5, a25, i5 )" ) "NLCG iterations:", iter," function evaluations:", nfunc
            write( ioInvLog, "( a25, i5, a25, i5 )" ) "NLCG iterations:", iter," function evaluations:", nfunc
            !
            close( ioInvLog )
            !
            ! Verbose
            write( *, * ) "     - Finish Inversion NLCG, output files in [", trim( outdir_name ), "]"
            !
            if( allocated( gPrev ) ) deallocate( gPrev )
            deallocate( mHat, grad, g, h )
            !
        else
            !
            write( *, * ) "Error opening [", trim( outdir_name )//"/NLCG.log", "] in writeData!"
            stop
            !
        endif
        !
    end subroutine solveInversionNLCG
    !
    !> Computes the gradient of the penalty functional,
    !> using EM solution (e_all) and the predicted data (dHat)
    !> Here, mHat denotes the non-regularized model parameter that
    !> is normally referred to as \tilde{dsigma} = C_m^{-1/2}(dsigma - m_0),
    !> and the gradient is computed with respect to \tilde{dsigma}.
    !> Before calling this routine, the forward solver must be run:
    !> call CmSqrtMult(mHat,dsigma)
    !> call linComb(ONE,dsigma,ONE,sigma,dsigma)
    !> call fwdPred(dsigma,dHat,e_all)
    !
    subroutine gradient( self, all_data, sigma, mHat, grad, dHat, i_sol )
        implicit none
        !
        class( InversionNLCG_t ), intent( in ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( in ) :: sigma, mHat
        class( ModelParameter_t ), allocatable, intent( inout ) :: grad
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: dHat
        integer, intent( in ) :: i_sol
        !
        real( kind=prec ) :: Ndata, Nmodel
        type( DataGroupTx_t ), allocatable, dimension(:) :: res
        class( ModelParameter_t ), allocatable :: dsigma, JTd, CmJTd
        class( Scalar_t ), allocatable, dimension(:) :: s_hat
        !
        ! compute the smoothed model parameter vector
        call model_cov%multBy_CmSqrt( mHat, dsigma )
        !
        ! overwriting the input with output
        call dsigma%linComb( ONE, ONE, sigma )
        !
        ! initialize res
        res = all_data
        !
        ! compute residual: res = (all_data-dHat)/Ndata
        !call linCombData( ONE, all_data, MinusONE, dHat, res )
        call subData( res, dHat )
        !
        Ndata = countValues( dHat )
        !
        call cdInvMult( res )
        !
        !> ????
        allocate( rScalar3D_SG_t :: s_hat( size( all_data ) ) )
        !
#ifdef MPI
        call masterJMult_T( dsigma, res, JTd, i_sol, s_hat )
#else
        call serialJMult_T( dsigma, res, JTd, i_sol, s_hat )
#endif
        !
        !> ????
        call weightGradrients( s_hat, all_data, dHat, JTd )
        !
        call model_cov%multBy_CmSqrt( JTd, CmJTd )
        !
        if( allocated( grad ) ) deallocate( grad )
        allocate( grad, source = CmJTd )
        !
        ! compute the number of data and model parameters for scaling
        Nmodel = mHat%countModel()
        !
        ! multiply by 2 (to be consistent with the formula)
        ! and add the gradient of the model norm
        !
        !call linComb(MinusTWO/Ndata,CmJTd,TWO*lambda/Nmodel,mHat,grad)
        call grad%linComb( MinusTWO / Ndata, TWO * self%lambda / Nmodel, mHat )
        !
        !call deallocateDataGroupTxArray( res )
        !
        deallocate( dsigma, JTd, CmJTd )
        !
    end subroutine gradient
    !
    !> ????
    !
    subroutine writeHeaders()
        implicit none
        !
        integer :: i_tx, n_tx, ios
        !
        n_tx = size( transmitters )
        !
        open( unit = ioGradNorm, file = trim( outdir_name )//"/GradNorm.log", status="unknown", position="append", iostat=ios )
        !
        open( unit = ioGradRMS, file = trim( outdir_name )//"/GradRMS.log", status="unknown", position="append", iostat=ios )
        !
        do i_tx = 1, n_tx + 1
            !
            if( i_tx == n_tx + 1 ) then
                !
                write( ioGradNorm, "( A12 )" ) "Grad"
                write( ioGradRMS, "( A12 )" ) "RMS"
                !
            else
                !
                !> Instantiate Transmitter's Source - According to transmitter type and chosen via control file
                select type( Tx => getTransmitter( i_tx ) )
                    !
                    class is( TransmitterMT_t )
                        !
                        write( ioGradNorm, "( A12, i3, A2 )", advance = "no" ) "MT", i_tx, ", "
                        write( ioGradRMS, "( A12, i3, A2 )", advance = "no" ) "MT", i_tx, ", "
                        !
                    class is( TransmitterCSEM_t )
                        !
                        write( ioGradNorm, "( A12, i3, A2 )", advance = "no" ) "CSEM", i_tx, ", "
                        write( ioGradRMS, "( A12, i3, A2 )", advance = "no" ) "CSEM", i_tx, ", "
                        !
                    class default
                        stop "Error: writeHeaders > Unclassified Transmitter"
                    !
                end select
                !
            endif
        enddo
        !
        close( ioGradNorm )
        !
        close( ioGradRMS )
        !
    end subroutine writeHeaders
    !
    subroutine weightGradrients( s_hat, d, dHat, JTd )
        implicit none
        !
        class( Scalar_t ), allocatable, dimension(:), intent( in ) :: s_hat
        type( DataGroupTx_t ), dimension(:), intent( in ) :: d, dHat
        class( ModelParameter_t ), intent( inout ) :: JTd
        !
        integer :: Ndata, i_tx, n_tx, ios
        class( Scalar_t ), allocatable :: temp_scalar
        real( kind=prec ) :: grad_norm, rms, SS
        real( kind=prec ) :: sum_tx_grad_norm, mt_grad_norm, csem_grad_norm
        real( kind=prec ) :: sum_tx_rms, mt_rms, csem_rms
        type( DataGroupTx_t ), allocatable, dimension(:) :: res, Nres
        real( kind=prec ), allocatable, dimension(:) :: tx_weights, tx_grad_norms, tx_rms
        !
        !> Open a file for output the operations realized on the gradient
        open( unit = ioGradLog, file = trim( outdir_name )//"/GradNLCG.log", status="unknown", position="append", iostat=ios )
        !
        open( unit = ioGradNorm, file = trim( outdir_name )//"/GradNorm.log", status="unknown", position="append", iostat=ios )
        !
        open( unit = ioGradRMS, file = trim( outdir_name )//"/GradRMS.log", status="unknown", position="append", iostat=ios )
        !
        if( ios == 0 ) then
            !
            write( ioGradLog, * ) "##########################"
            !
            res = d
            !
            !> res = res - dHat
            call subData( res, dHat )
            !
            call cdInvMult( res, Nres )
            !
            !> Zero the variables for summation
            sum_tx_grad_norm = 0.0
            sum_tx_rms = 0.0
            mt_grad_norm = 0.0
            csem_grad_norm = 0.0
            mt_rms = 0.0
            csem_rms = 0.0
            !
            !> Allocate arrays to store values for each of the n_tx transmitters
            n_tx = size( transmitters )
            !
            allocate( tx_weights( n_tx ) )
            allocate( tx_grad_norms( n_tx ) )
            allocate( tx_rms( n_tx ) )
            !
            do i_tx = 1, n_tx
                !
                tx_grad_norms( i_tx ) = sqrt( s_hat( i_tx )%dotProd( s_hat( i_tx ) ) )
                !
                write( ioGradNorm, "( f15.3, A1 )", advance = "no" ) tx_grad_norms( i_tx ), ","
                !
                tx_rms( i_tx ) = sqrt( res( i_tx )%dotProd( Nres( i_tx ) ) / ( size( Nres( i_tx )%data ) * Nres( i_tx )%data(1)%n_comp * 2 ) )
                !
                write( ioGradRMS, "( f15.3, A1 )", advance = "no" ) tx_rms( i_tx ), ","
                !
                !> Instantiate Transmitter's Source - According to transmitter type and chosen via control file
                select type( Tx => getTransmitter( i_tx ) )
                    !
                    class is( TransmitterMT_t )
                        !
                        mt_grad_norm = mt_grad_norm + tx_grad_norms( i_tx )
                        mt_rms = mt_rms + tx_rms( i_tx )
                        !
                        write( *, * ) "MT > GRAD_NORM, RMS: ", tx_grad_norms( i_tx ), tx_rms( i_tx )
                        write( ioGradLog, * ) "MT > GRAD_NORM, RMS: ", tx_grad_norms( i_tx ), tx_rms( i_tx )
                        !
                    class is( TransmitterCSEM_t )
                        !
                        csem_grad_norm = csem_grad_norm + tx_grad_norms( i_tx )
                        csem_rms = csem_rms + tx_rms( i_tx )
                        !
                        write( *, * ) "CSEM > GRAD_NORM, RMS: ", tx_grad_norms( i_tx ), tx_rms( i_tx )
                        write( ioGradLog, * ) "CSEM > GRAD_NORM, RMS: ", tx_grad_norms( i_tx ), tx_rms( i_tx )
                        !
                    class default
                        stop "Error: handleGradVectorsNLCG > Unclassified Transmitter"
                    !
                end select
                !
                sum_tx_grad_norm = sum_tx_grad_norm + tx_grad_norms( i_tx )
                sum_tx_rms = sum_tx_rms + tx_rms( i_tx )
                !
            enddo
            !
            !> Get gradient conductivity
            call JTd%getCond( temp_scalar )
            !
            grad_norm = real( sqrt( temp_scalar%dotProd( temp_scalar ) ), kind=prec )
            !
            write( *, * ) "MODEL GRAD NORM: ", grad_norm
            write( ioGradLog, * ) "MODEL GRAD NORM: ", grad_norm
            !
            write( *, * ) "SUM TX GRAD NORM = MT + CSEM", sum_tx_grad_norm, mt_grad_norm, csem_grad_norm
            write( ioGradLog, * ) "SUM TX GRAD NORM = MT + CSEM", sum_tx_grad_norm, mt_grad_norm, csem_grad_norm
            !
            SS = dotProdData( res, Nres )
            !
            Ndata = countValues( res )
            !
            rms = sqrt( SS / Ndata )
            !
            write( *, * ) "RMS: ", rms
            !
            write( ioGradRMS, "( f15.3 )" ) rms
            !
            close( ioGradRMS )
            !
            write( ioGradLog, * ) "RMS: ", rms
            !
            write( *, * ) "AVG/SUM TX RMS = MT + CSEM", sum_tx_rms/n_tx, sum_tx_rms, mt_rms, csem_rms
            write( ioGradLog, * ) "AVG/SUM TX RMS = MT + CSEM", sum_tx_rms/n_tx, sum_tx_rms, mt_rms, csem_rms
            !
            do i_tx = 1, n_tx
                !
                tx_weights( i_tx ) = sum_tx_grad_norm / ( tx_grad_norms( i_tx ) * n_tx )
                !
            enddo
            !
            write( *, * ) "TX WEIGHTS: ", tx_weights
            write( ioGradLog, * ) "TX WEIGHTS: ", tx_weights
            !
            if( joint_type /= INV_UNWEIGHTED ) then
                !
                !> Set back the weighted gradient conductivity
                call JTd%zeros()
                do i_tx = 1, n_tx
                    !
                    temp_scalar = s_hat( i_tx )
                    !
                    call temp_scalar%mult( tx_weights( i_tx ) )
                    !
                    call JTd%addCond( temp_scalar )
                    !
                enddo
                !
            endif
            !
            deallocate( temp_scalar )
            !
            !> Get gradient conductivity
            call JTd%getCond( temp_scalar )
            !
            grad_norm = real( sqrt( temp_scalar%dotProd( temp_scalar ) ), kind=prec )
            !
            
            write( *, * ) "FINAL GRAD NORM: ", grad_norm
            !
            write( ioGradNorm, "( f15.3 )" ) grad_norm
            !
            close( ioGradNorm )
            !
            write( ioGradLog, * ) "FINAL GRAD NORM: ", grad_norm
            !
            deallocate( tx_weights, tx_grad_norms, tx_rms, temp_scalar )
            !
            close( ioGradLog )
            !
        else
            !
            write( *, * ) "Error opening [", trim( outdir_name )//"/GradNLCG.log", "] in handleGradVectorsNLCG!"
            stop
            !
        endif
        !
    end subroutine weightGradrients
    !
    !> Compute the full penalty functional F
    !> Also output the predicted data and the EM solution
    !> that can be used for evaluating the gradient
    !
    subroutine func( self, all_data, sigma, mHat, F, mNorm, dHat, i_sol, rms )
        implicit none
        !
        class( InversionNLCG_t ), intent( inout ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( in ) :: sigma, mHat
        real( kind=prec ), intent( out ) :: F, mNorm
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: dHat
        integer, intent( inout ) :: i_sol
        real( kind=prec ), optional, intent( out ) :: rms
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: res, Nres
        class( ModelParameter_t ), allocatable :: dsigma, JTd
        real( kind=prec ) :: SS, angle2, angle1, diff, diff1
        integer :: Ndata, Nmodel, j, i, isite
        !
        ! compute the smoothed model parameter vector
        call model_cov%multBy_CmSqrt( mHat, dsigma )
        !
        ! overwriting input with output
        call dsigma%linComb( ONE, ONE, sigma )
        !
        ! initialize dHat
        dHat = all_data
        !
#ifdef MPI
        call masterForwardModelling( dsigma, dHat, i_sol )
#else
        call serialForwardModeling( dsigma, dHat, i_sol )
#endif
        !
        !> initialize res
        res = all_data
        !
        !call linCombData( ONE, all_data, MinusONE, dHat, res )
        call subData( res, dHat )
        !
        !> normalize residuals, compute sum of squares
        call cdInvMult( res, Nres )
        !
        SS = dotProdData( res, Nres )
        !
        Ndata = countValues( res )
        !
        mNorm = mHat%dotProd( mHat )
        !
        Nmodel = mHat%countModel()
        !
        !> penalty functional = sum of squares + scaled model norm
        F = SS / Ndata + ( self%lambda * mNorm / Nmodel )
        !
        !> scale mNorm for output
        mNorm = mNorm / Nmodel
        !
        ! if required, compute the Root Mean Squared misfit
        if( present( rms ) ) then
            rms = sqrt( SS / Ndata )
        endif
        !
        !call deallocateDataGroupTxArray( res )
        !
        !call deallocateDataGroupTxArray( Nres )
        !
        deallocate( dsigma )
        !
    end subroutine func
    !
    !> Divides by the data covariance C_d, which is a diagonal
    !> operator. Divides by the variances (squared error bars)
    !> and scales by the number of data (degrees of freedom).
    !
    subroutine cdInvMult( d_in, d_out )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: d_in
        type( DataGroupTx_t ), allocatable, dimension(:), optional, intent( out ) :: d_out
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: all_data
        !
        all_data = d_in
        !
        call normalizeData( all_data, 2 )
        !
        if( present( d_out ) ) then
            d_out = all_data
        else
            d_in = all_data
        endif
        !
        !call deallocateDataGroupTxArray( all_data )
        !
    end subroutine cdInvMult
    !
    !>
    subroutine updateDampingParameter( self, mHat, F, grad )
        implicit none
        !
        class( InversionNLCG_t ), intent( inout ) :: self
        class( ModelParameter_t ), allocatable, intent( in ) :: mHat
        real( kind=prec ), intent( inout ) :: F
        class( ModelParameter_t ), allocatable, intent( inout ) :: grad
        !
        real( kind=prec ) :: SS, mNorm, Nmodel
        class( ModelParameter_t ), allocatable :: dSS
        !
        ! compute the model norm
        mNorm = mHat%dotProd( mHat )
        !
        Nmodel = mHat%countModel()
        !
        ! (scaled) sum of squares = penalty functional - scaled model norm
        SS = F - ( self%lambda * mNorm / Nmodel )
        !
        ! initialize
        allocate( dSS, source = mHat )
        !
        !> subtract the model norm derivative from the gradient of the penalty functional
        dSS = grad
        !
        call dSS%linComb( ONE, MinusTWO * self%lambda / Nmodel, mHat )
        !
        ! update the damping parameter lambda
        self%lambda = self%lambda / self%k
        !
        ! penalty functional = (scaled) sum of squares + scaled model norm
        F = SS + ( self%lambda * mNorm / Nmodel )
        !
        ! add the model norm derivative to the gradient of the penalty functional
        grad = dSS
        !
        call grad%linComb( ONE, TWO * self%lambda / Nmodel, mHat )
        !
        deallocate( dSS )
        !
    end subroutine updateDampingParameter
    !
    !> Line search that is based on the Numerical Recipes and on the
    !> text by Michael Ferris, Chapter 3, p 59. We only test the sufficient
    !> decrease (Armijo) condition (ignoring the curvature condition).
    !> We first interpolate using a quadratic approximation; if the
    !> solution does not satisfy the condition, we backtrack using
    !> cubic interpolation. This strategy only requires one gradient
    !> evaluation and is very efficient when computing gradients is
    !> expensive.
    !
    !> The initial step size is set outside of this routine (in the NLCG)
    !> but these are the good choices (ref. Michael Ferris, Chapter 3, p 59):
    !> alpha_1 = alpha_{k-1} dotProd(grad_{k-1},h_{k-1})/dotProd(grad_k,h_k})
    !> or interpolate the quadratic to f(m_{k-1}), f(m_k) and
    !> dotProd(grad_{k-1},h_{k-1}) and find the minimizer
    !> alpha_1 = 2(f_k - f_{k-1})/dotProd(grad_{k-1},h_{k-1}),
    !> the update alpha_1 <- min(1.00,1.01 * alpha_1).
    !
    !> Set f(alpha) = func(mHat + alpha*h). Fit the quadratic
    !>     f_q(alpha) = a alpha^2 + b alpha + f(0)
    !> using the information f(0), f"(0) and f(alpha_1) to obtain
    !> a = (f(alpha_1) - f(0) - f"(0) alpha_1)/(alpha_1 * alpha_1),
    !> b = f"(0).
    !> Then, the minimum point of the quadratic is alpha_q = -b/(2a),
    !> assuming that a > 0. If this try is not successful, fit a cubic
    !>     f_c(alpha) = a alpha^3 + b alpha^2 + f"(0) alpha + f(0)
    !> using f(0), f"(0), f(alpha_1) and f(alpha_q). Repeat as necessary.
    !> Here, a and b are as described in the code.
    !> A new cubic is not identical to a previous curve since f_c is only
    !> an approximation to f: in general, f(alpha_c) /= f_c(alpha_c),
    !> hence the new point does not lie on the approximating curve.
    !
    !> Our solution has to satisfy the sufficient decrease condition
    !>     f(alpha) < f(0) + c alpha f"(0).
    !
    !> The optional relaxation parameter gamma is needed for algorithms
    !> like the Renormalised Steepest Descent (RSD). See the dynamical
    !> systems in optimisation research (Pronzato et al [2000, 2001]).
    !> To the best of my knowledge, it is not useful for NLCG.
    !
    subroutine lineSearchCubic( self, all_data, sigma, h, alpha, mHat, f, grad, rms, niter, dHat, gamma )
        implicit none
        !
        class( InversionNLCG_t ), intent( inout ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( in ) :: sigma, h    ! search direction
        real( kind=prec ), intent( inout ) :: alpha ! step size
        class( ModelParameter_t ), allocatable, intent( inout ) :: mHat
        real( kind=prec ), intent( inout ) :: f
        class( ModelParameter_t ), allocatable, intent( inout ) :: grad
        real( kind=prec ), intent( out ) :: rms
        integer, intent( out ) :: niter
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: dHat
        !
        ! optionally add relaxation (e.g. for Renormalized Steepest Descent)
        real( kind=prec ), intent( in ), optional :: gamma
        !
        real( kind=prec ) :: alpha_1, alpha_i, alpha_j, mNorm
        logical :: starting_guess, relaxation
        real( kind=prec ) :: eps, k, c, a, b, q1, q2, q3
        real( kind=prec ) :: g_0, f_0, f_1, f_i, f_j, rms_1, mNorm_1
        class( ModelParameter_t ), allocatable :: mHat_0, mHat_1
        type( DataGroupTx_t ), allocatable, dimension(:) :: dHat_1
        integer :: i_sol
        !
        ! parameters
        c = self%c
        !
        ! initialize the line search
        niter = 0
        !
        allocate( mHat_0, source = mHat )
        !
        f_0 = f
        !
        starting_guess = .FALSE.
        !
        ! g_0 is the directional derivative f"(0) = (df/dm).dot.h
        g_0 = grad%dotProd( h )
        !
        ! alpha_1 is the initial step size, which is set in NLCG
        alpha_1 = alpha
        !
        ! with relaxation, we specify gamma = 1 - eps, eps > 0 small; then the final
        ! solution is f(gamma*alpha) = func(mHat + gamma*alpha*h)
        if( present( gamma ) ) then
            relaxation = .TRUE.
        else
            relaxation = .FALSE.
        endif
        !
        ! compute the trial mHat, f, dHat, e_all, rms
        allocate( mHat_1, source = mHat_0 )
        !
        call mHat_1%linComb( ONE, alpha_1, h )
        !
        i_sol = 1
        !
        call self%func( all_data, sigma, mHat_1, f_1, mNorm_1, dHat_1, i_sol, rms_1 )
        !
        write( *, * ) "STARTLS > lambda, alpha, f_1, mNorm_1, rms_1:", self%lambda, alpha, f_1, mNorm_1, rms_1
        write( ioInvLog, * ) "STARTLS > lambda, alpha, f_1, mNorm_1, rms_1:", self%lambda, alpha, f_1, mNorm_1, rms_1
        !
        niter = niter + 1
        !
        if( f_1 - f_0 >= R_LARGE ) then
            write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m Try a smaller starting r_value of alpha."
            write( *, * ) "Exiting..."
            stop
        endif
        !
        ! try fitting a quadratic
        a = ( f_1 - f_0 - g_0 * alpha_1 ) / ( alpha_1 ** 2 )
        !
        b = g_0
        !
        ! if the curvature is -ve, there is no minimum; take the initial guess
        if( a < 0 ) then
            !
            starting_guess = .TRUE.
            alpha = alpha_1
            dHat = dHat_1
            i_sol = 1
            !
            if( allocated( mHat ) ) deallocate( mHat )
            allocate( mHat, source = mHat_1 )
            !
            rms = rms_1
            f = f_1
            !
            ! compute the gradient and exit
            if( relaxation ) then
                !
                if( allocated( mHat ) ) deallocate( mHat )
                allocate( mHat, source = mHat_0 )
                !
                call mHat%linComb( ONE, gamma * alpha, h )
                !
                i_sol = 0
                !
                call self%func( all_data, sigma, mHat, f, mNorm, dHat, i_sol, rms )
                !
                write( *, * ) "lambda, gamma * alpha, f, mNorm, rms:", self%lambda, gamma * alpha, f, mNorm, rms
                write( ioInvLog, * ) "lambda, gamma * alpha, f, mNorm, rms:", self%lambda, gamma * alpha, f, mNorm, rms
                !
            endif
            !
            call self%gradient( all_data, sigma, mHat, grad, dHat, i_sol )
            !
            write( *, * ) "Quadratic has no minimum, exiting line search"
            write( ioInvLog, * ) "Quadratic has no minimum, exiting line search"
            !
            !call deallocateDataGroupTxArray( dHat_1 )
            !
            deallocate( mHat_0, mHat_1 )
            !
            return
            !
        endif
        !
        ! otherwise compute the functional at the minimizer of the quadratic
        alpha = - b / ( TWO * a )
        !
        if( allocated( mHat ) ) deallocate( mHat )
        allocate( mHat, source = mHat_0 )
        !
        call mHat%linComb( ONE, alpha, h )
        !
        i_sol = 0
        !
        call func( self, all_data, sigma, mHat, f, mNorm, dHat, i_sol, rms )
        !
        write( *, * ) "QUADLS: lambda, alpha, f, mNorm, rms:", self%lambda, alpha, f, mNorm, rms
        write( ioInvLog, * ) "QUADLS: lambda, alpha, f, mNorm, rms:", self%lambda, alpha, f, mNorm, rms
        !
        niter = niter + 1
        !
        ! check whether the solution satisfies the sufficient decrease condition
        if( f < f_0 + c * alpha * g_0 ) then
            !
            ! if the initial guess was better than what we found, take it
            if( f_1 < f ) then
                !
                starting_guess = .TRUE.
                alpha = alpha_1
                dHat = dHat_1
                i_sol = 1
                !
                if( allocated( mHat ) ) deallocate( mHat )
                allocate( mHat, source = mHat_1 )
                !
                rms = rms_1
                f = f_1
                !
            endif
            !
            ! compute the gradient and exit
            if( relaxation ) then
                !
                if( allocated( mHat ) ) deallocate( mHat )
                allocate( mHat, source = mHat_0 )
                !
                call mHat%linComb( ONE, gamma * alpha, h )
                !
                i_sol = 0
                !
                call self%func( all_data, sigma, mHat, f, mNorm, dHat, i_sol, rms )
                !
                write( *, * ) "QUADLS: lambda, gamma*alpha, f, mNorm, rms:", self%lambda, gamma*alpha, f, mNorm, rms
                write( ioInvLog, * ) "QUADLS: lambda, gamma*alpha, f, mNorm, rms:", self%lambda, gamma*alpha, f, mNorm, rms
                !
            endif
            !
            call self%gradient( all_data, sigma, mHat, grad, dHat, i_sol )
            !
            write( *, * ) "Sufficient decrease condition satisfied, exiting line search"
            write( ioInvLog, * ) "Sufficient decrease condition satisfied, exiting line search"
            !
            !call deallocateDataGroupTxArray( dHat_1 )
            !
            deallocate( mHat_0, mHat_1 )
            !
            return
            !
        endif
        !
        ! this should not happen, but in practice it is possible to end up with
        ! a function increase at this point (e.g. in the current global code).
        ! Most likely, this is due to an inaccuracy in the gradient computations.
        ! In this case, we avoid an infinite loop by exiting line search.
        ! It is also possible that both f_1 and f are worse than the starting r_value!
        ! Then, take whichever is smaller. Ideally, want to decrease the tolerance
        ! for gradient computations if this happens.
        if( f > f_0 ) then
            !
            write( *, * ) "Unable to fit a quadratic due to bad gradient estimate, exiting line search"
            write( ioInvLog, * ) "Unable to fit a quadratic due to bad gradient estimate, exiting line search"
            !
        else
            !
            ! fit a cubic and backtrack (initialize)
            alpha_i = alpha_1
            f_i = f_1
            !
            alpha_j = alpha
            f_j = f
            !
            fit_cubic: do
                !
                ! compute the minimizer
                q1 = f_i - f_0 - g_0 * alpha_i
                !
                q2 = f_j - f_0 - g_0 * alpha_j
                !
                q3 = alpha_i**2 * alpha_j**2 * (alpha_j - alpha_i)
                !
                a = (alpha_i**2 * q2 - alpha_j**2 * q1)/q3
                !
                b = (alpha_j**3 * q1 - alpha_i**3 * q2)/q3
                !
                alpha = (- b + sqrt(b*b - 3*a*g_0))/(3*a)
                !
                ! compute the penalty functional
                !
                if( allocated( mHat ) ) deallocate( mHat )
                allocate( mHat, source = mHat_0 )
                !
                call mHat%linComb( ONE, alpha, h )
                !
                i_sol = 0
                !
                call self%func( all_data, sigma, mHat, f, mNorm, dHat, i_sol, rms )
                !
                write( *, * ) "CUBICLS: lambda, alpha, f, mNorm, rms:", self%lambda, alpha, f, mNorm, rms
                write( ioInvLog, * ) "CUBICLS: lambda, alpha, f, mNorm, rms:", self%lambda, alpha, f, mNorm, rms
                !
                niter = niter + 1
                !
                ! check whether the solution satisfies the sufficient decrease condition
                if( f < f_0 + c * alpha * g_0 ) then
                    exit
                endif
                !
                ! if not, iterate, using the two most recent values of f & alpha
                alpha_i = alpha_j
                f_i = f_j
                !
                alpha_j = alpha
                f_j = f
                !
                ! check that the function still decreases to avoid infinite loops in case of a bug
                if( abs( f_j - f_i ) < TOL8 ) then
                    !
                    write( *, * ) "Warning: exiting cubic search since the function no longer decreases!"
                    write( ioInvLog, * ) "Warning: exiting cubic search since the function no longer decreases!"
                    !
                    exit
                    !
                endif
                !
            enddo fit_cubic
            !
        endif
        !
        if( f_1 < f ) then
            !
            starting_guess = .TRUE.
            !
        endif
        !
        ! if the initial guess was better than what we found, take it
        if( starting_guess ) then
            !
            alpha = alpha_1
            dHat = dHat_1
            i_sol = 1
            !
            if( allocated( mHat ) ) deallocate( mHat )
            allocate( mHat, source = mHat_1 )
            !
            rms = rms_1
            f = f_1
            !
        endif
        !
        ! compute gradient of the full penalty functional and exit
        if( relaxation ) then
            !
            if( allocated( mHat ) ) deallocate( mHat )
            allocate( mHat, source = mHat_0 )
            !
            call mHat%linComb( ONE, gamma*alpha, h )
            !
            i_sol = 0
            !
            call self%func( all_data, sigma, mHat, f,mNorm, dHat, i_sol, rms )
            !
            write( *, * ) "RELAX: lambda, gamma*alpha, f, mNorm, rms:", self%lambda, gamma*alpha, f, mNorm, rms
            write( ioInvLog, * ) "RELAX: lambda, gamma*alpha, f, mNorm, rms:", self%lambda, gamma*alpha, f, mNorm, rms
            !
        endif
        !
        call self%gradient( all_data, sigma, mHat, grad, dHat, i_sol )
        !
        write( *, * ) "Gradient computed, line search finished"
        write( ioInvLog, * ) "Gradient computed, line search finished"
        !
        !call deallocateDataGroupTxArray(dHat_1)
        !
        deallocate( mHat_0, mHat_1 )
        !
    end subroutine lineSearchCubic
    !
    !> ????
    !
    subroutine outputFilesInversionNLCG( iter, dHat, res, dsigma, mHat )
        implicit none
        !
        integer, intent( in ) :: iter
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: dHat, res
        class( ModelParameter_t ), intent( in ) :: dsigma, mHat
        !
        character(100) :: out_file_name
        character(3) :: char3
        !
        write( char3, "(i3.3)" ) iter
        !
        !> Write predicted data for this NLCG iteration
        out_file_name = trim( outdir_name )//"/PredictedData_NLCG_"//char3//".dat"
        !
        call writeData( dHat, trim( out_file_name ) )
        !
        !> Write residual data for this NLCG iteration
        out_file_name = trim( outdir_name )//"/ResidualData_NLCG_"//char3//".dat"
        !
        call writeData( res, trim( out_file_name ) )
        !
        !> Write model for this NLCG iteration
        out_file_name = trim( outdir_name )//"/SigmaModel_NLCG_"//char3//".rho"
        !
        call dsigma%write( trim( out_file_name ) )
        !
        !> Write perturbation model for this NLCG iteration
        out_file_name = trim( outdir_name )//"/PerturbationModel_NLCG_"//char3//".rho"
        !
        call mHat%write( trim( out_file_name ) )
        !
    end subroutine outputFilesInversionNLCG
    !
end module InversionNLCG
!