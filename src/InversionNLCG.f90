!
!> Module with the InversionNLCG routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module InversionNLCG
    !
    use Sensitivity
    !
    type :: NLCGiterControl_t
        !
        ! maximum number of iterations in one call to iterative solver
        integer :: maxIter
        ! convergence criteria: return from solver if rms < rmsTol
        real( kind=prec ) :: rmsTol
        ! the condition to identify when the inversion stalls
        real( kind=prec ) :: fdiffTol
        ! initial r_value of lambda (will not override the NLCG input argument)
        real( kind=prec ) :: lambda
        ! exit if lambda < lambdaTol approx. 1e-4
        real( kind=prec ) :: lambdaTol
        ! set lambda_i = lambda_{i-1}/k when the inversion stalls
        real( kind=prec ) :: k
        ! the factor that ensures sufficient decrease in the line search
        real( kind=prec ) :: c
        ! the factor that ensures curvature condition in the line search
        real( kind=prec ) :: c2
        ! restart CG every nCGmax iterations to ensure conjugacy
        integer :: nCGmax
        ! restart CG if orthogonality is lost (not necessarily needed)
        ! real( kind=prec ) :: delta ! 0.5
        ! the starting step for the line search
        real( kind=prec ) :: alpha_1
        ! if alpha_{i+1} < alpha_i * k_{alpha}, set alpha_{i+1} = alpha_i/2
        ! real( kind=prec ) :: alpha_k ! 0.1
        ! if alpha_{i+1} - alpha_i < tol_{alpha}, set alpha_{i+1} = alpha_i/2
        ! real( kind=prec ) :: alpha_tol ! 1.0e-2
        ! maximum initial delta mHat (overrides alpha_1)
        real( kind=prec ) :: startdm
        ! optional relaxation parameter (Renormalized Steepest Descent algorithm)
        real( kind=prec ) :: gamma
        ! model and data output file name
        character(80) :: fname
        !
    end type NLCGiterControl_t
    !
    type( NLCGiterControl_t ), private, save :: NLCGiterControl
    !
contains
    !
    !>
    subroutine set_NLCGiterControl( NLCGiterControl )
        implicit none
        !
        type( NLCGiterControl_t ), intent( inout ) :: NLCGiterControl
        !
        ! maximum number of iterations in one call to iterative solver
        NLCGiterControl%maxIter = 600
        ! convergence criteria: return from solver if rms < rmsTol
        NLCGiterControl%rmsTol = 1.05
        ! inversion stalls when abs(rms - rmsPrev) < fdiffTol (2e-3 works well)
        NLCGiterControl%fdiffTol = 2.0e-3
        ! initial r_value of lambda (will not override the NLCG input argument)
        NLCGiterControl%lambda = 1.
        ! exit if lambda < lambdaTol approx. 1e-4
        NLCGiterControl%lambdaTol = 1.0e-8
        ! set lambda_i = lambda_{i-1}/k when the inversion stalls
        NLCGiterControl%k = 10.
        ! the factor that ensures sufficient decrease in the line search >=1e-4
        NLCGiterControl%c = 1.0e-4
        ! restart CG every nCGmax iterations to ensure conjugacy
        NLCGiterControl%nCGmax = 8
        ! the starting step for the line search
        NLCGiterControl%alpha_1 = 20.
        ! maximum initial delta mHat (overrides alpha_1)
        NLCGiterControl%startdm = 20.
        ! optional relaxation parameter (Renormalized Steepest Descent algorithm)
        NLCGiterControl%gamma = 0.99
        ! model and data output file name
        NLCGiterControl%fname = 'Modular'
        !
    end subroutine set_NLCGiterControl
    !
    !>  computes inverse solution minimizing penalty functional
    !>  for fixed r_value of regularization parameter, using
    !>  a variant of non-linear conjugate gradient search.
    !>  Various flavours of the algorithm and of the line search
    !>  can be called from this routine
    !>
    !> Note about the starting model:
    !> The starting model has to be in the smoothed model space,
    !> i.e. of the form m = C_m^{1/2} \tilde{m} + m_0.
    !> In order to compute \tilde{m} from the starting model,
    !> C_m^{-1/2} has to be implemented. To avoid this issue
    !> altogether, we are always starting from the prior,
    !> with \tilde{m} = 0. However, in general we could also
    !> start with the result of a previous search.
    !
    subroutine NLCGsolver( d, m0, m )
        implicit none
        !
        !> d is data; on output it contains the responses for the inverse model
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: d
        !
        !> m0 is prior model parameter
        class( ModelParameter_t ), allocatable, intent( in ) :: m0
        !> m is solution parameter ... on input m contains starting guess
        class( ModelParameter_t ), allocatable, intent( inout ) :: m
        !
        !> initial step size in the line search direction in model units
        real( kind=prec ) :: startdm
        !> flavor is a string that specifies the algorithm to use
        character(80) :: flavor = 'Cubic'
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: dHat, res
        class( ModelParameter_t ), allocatable :: mHat, grad, g, h, gPrev
        real( kind=prec ) :: r_value, valuePrev, rms
        real( kind=prec ) :: rmsPrev, alpha, beta
        real( kind=prec ) :: gnorm, mNorm, Nmodel
        real( kind=prec ) :: grad_dot_h, g_dot_g
        real( kind=prec ) :: g_dot_gPrev,g_dot_h
        real( kind=prec ) :: gPrev_dot_gPrev 
        real( kind=prec ) :: h_dot_g, h_dot_gPrev
        integer :: iter, nCG, nLS, nfunc
        type( ESolMTx ) :: e_all
        !
        !>
        call createOutputDirectory()
        !
        ! Verbose
        write( *, * ) "     - Start Inversion NLCG, output files in [", trim( outdir_name ), "]"
        !
        call set_NLCGiterControl( NLCGiterControl )
        !
        ! initialize the line search
        alpha = NLCGiterControl%alpha_1
        !
        startdm = NLCGiterControl%startdm
        !
        write( *, * ) "lambda, startdm: ", NLCGiterControl%lambda, startdm
        !
        ! starting model contains the rough deviations from the prior
        allocate( mHat, source = m )
        !
        !  compute the penalty functional and predicted data
        e_all%SolnIndex = 0
        !
        call func( NLCGiterControl%lambda, d, m0, mHat, r_value, mNorm, dHat, e_all, rms )
        !
        nfunc = 1
        !
        ! output (smoothed) initial model and responses for later reference
        m = model_cov%multBy_CmSqrt( mHat )
        !
        call m%linComb( ONE, ONE, m0 )
        !
        !> compute gradient of the full penalty functional
        call gradient( NLCGiterControl%lambda, d, m0, mHat, grad, dHat, e_all )
        !
        gnorm = sqrt( grad%dotProd( grad ) )
        !
        !write( *, * ) "gnorm: ", gnorm
        !stop
        !
        if ( gnorm < TOL6 ) then
            stop "Error: NLCGsolver: Problem with your gradient computations: first gradient is zero"
        else
            !
            alpha = startdm / gnorm
            !
            !write( *, * ) "alpha: ", alpha
            !stop
            !
        endif
        !
        !> initialize CG: g = - grad; h = g
        nCG = 0
        !
        iter = 0
        !
        allocate( g, source = grad )
        !
        call g%linComb( MinusONE, R_ZERO, grad )
        !
        allocate( h, source = g )
        !
        do
            !  test for convergence ...
            if( rms .LT. NLCGiterControl%rmsTol .OR. iter .GE. NLCGiterControl%maxIter ) then
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

            !write( *, * ) "grad_dot_h: ", grad_dot_h
            !stop

            ! at the end of line search, set mHat to the new r_value
            ! mHat = mHat + alpha*h  and evaluate gradient at new mHat
            ! data and solnVector only needed for output
            write( *, * ) "Starting line search..."
            !
            select case ( flavor )
                !
                case ( 'Cubic' )
                    call lineSearchCubic( NLCGiterControl%lambda, d, m0, h, alpha, mHat, r_value, grad, rms, nLS, dHat, e_all )
                    !call deall(e_all)
                case ('Quadratic')
                    !call lineSearchQuadratic(NLCGiterControl%lambda,d,m0,h,alpha,mHat,r_value,grad,rms,nLS,dHat,e_all)
                    !call deall(e_all)
                case ('Wolfe')
                    !call lineSearchWolfe(NLCGiterControl%lambda,d,m0,h,alpha,mHat,r_value,grad,rms,nLS,dHat,e_all)
                    !call deall(e_all)
                case default
                    stop "Error: NLCGsolver: Unknown line search requested in NLCG"
                !
            end select
            !
            nfunc = nfunc + nLS
            !
            if( allocated( gPrev ) ) then
                gPrev = g
            else
                allocate( gPrev, source = g )
            endif
            !
            g = grad
            !
            call g%linComb( MinusONE, R_ZERO, grad )
            !
            ! compute the starting step for the next line search
            alpha = 2 * ( r_value - valuePrev ) / grad_dot_h
            !
            ! adjust the starting step to ensure super linear convergence properties
            alpha = ( ONE + 0.01 ) * alpha
            !
            write( *, * ) "Completed NLCG iteration ", iter
            !
            Nmodel = mHat%countModel()
            !
            mNorm = mHat%dotProd( mHat ) / Nmodel
            !
            write( *, * ) "     lambda, alpha, r_value, mNorm, rms: ", NLCGiterControl%lambda, alpha, r_value, mNorm, rms
            !
            ! write out the intermediate model solution and responses
            m = model_cov%multBy_CmSqrt( mHat )
            !
            call m%linComb( ONE, ONE, m0 )
            !
            res = d
            !
            call linCombDataGroupTxArray( ONE, d, MinusONE, dHat, res )
            !
            !> if alpha is too small, we are not making progress: update lambda
            if( abs( rmsPrev - rms ) < NLCGiterControl%fdiffTol ) then
                !
                ! update lambda, penalty functional and gradient
                call update_damping_parameter( NLCGiterControl%lambda, mHat, r_value, grad )
                !
                ! check that lambda is still at a reasonable r_value
                if( NLCGiterControl%lambda < NLCGiterControl%lambdaTol ) then
                    stop "Error: NLCGsolver: Unable to get out of a local minimum."
                    exit
                endif
                !
                !> update alpha
                gnorm = sqrt( grad%dotProd( grad ) )
                !
                write( *, * ) "gnorm: ", gnorm
                !
                !> alpha = min(NLCGiterControl%alpha_1,startdm/gnorm)
                alpha = min( ONE, startdm ) / gnorm
                !
                write( *, * ) "alpha: ", alpha
                !
                !> g = - grad
                g = grad
                !
                call g%linComb( MinusONE, R_ZERO, grad )
                !
                !> restart
                write( *, * ) "Restarting NLCG with the damping parameter updated"
                write( *, * ) "lambda, alpha, r_value, mNorm, rms: ", NLCGiterControl%lambda, alpha, r_value, mNorm, rms
                !
                h = g
                !
                nCG = 0
                !
                cycle !????
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
            !> say we can take beta > 0 (didn't work as well)
            !> if ((beta.lt.R_ZERO).or.(g_dot_g + beta*g_dot_h .le. R_ZERO)&
            !>    .and.(nCG .ge. NLCGiterControl%nCGmax)) then  !PR+
            if( g_dot_g + beta * g_dot_h .LE. R_ZERO .AND. nCG .GE. NLCGiterControl%nCGmax ) then  !PR
                !
                ! restart
                write( *, * ) "Restarting NLCG to restore orthogonality"
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
            h = g
            !
            call h%linComb( ONE, beta, h )
            !
        end do
        !
        !> multiply by C^{1/2} and add m_0
        m = model_cov%multBy_CmSqrt( mHat )
        !
        call m%linComb( ONE, ONE, m0 )
        !
        d = dHat
        !
        ! cleaning up
        call deallocateDataGroupTxArray( dHat )
        call deallocateDataGroupTxArray( res )
        !
        ! Verbose
        write( *, * ) "     - Finish Inversion NLCG, output files in [", trim( outdir_name ), "]"
        !
        deallocate( mHat, grad, g, h, gPrev )
        !
    end subroutine NLCGsolver
    !
    !> Computes the gradient of the penalty functional,
    !> using EM solution (e_all) and the predicted data (dHat)
    !> Here, mHat denotes the non-regularized model parameter that
    !> is normally referred to as \tilde{m} = C_m^{-1/2}(m - m_0),
    !> and the gradient is computed with respect to \tilde{m}.
    !> Before calling this routine, the forward solver must be run:
    !> call CmSqrtMult(mHat,m)
    !> call linComb(ONE,m,ONE,m0,m)
    !> call fwdPred(m,dHat,e_all)
    !
    subroutine gradient( lambda, d, m0, mHat, grad, dHat, e_all )
        implicit none
        !
        real( kind=prec ), intent( in ) :: lambda
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: d
        class( ModelParameter_t ), allocatable, intent( in ) :: m0, mHat
        class( ModelParameter_t ), allocatable, intent( inout ) :: grad
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: dHat
        type( ESolMTx ), intent( inout ) :: e_all
        !
        real( kind=prec ) :: Ndata, Nmodel, angle1, angle2, diff, diff1
        type( DataGroupTx_t ), allocatable, dimension(:) :: res
        class( ModelParameter_t ), allocatable :: m, JTd, CmJTd
        integer :: j, i, icomp, isite
        !
        ! compute the smoothed model parameter vector
        allocate( m, source = model_cov%multBy_CmSqrt( mHat ) )
        !
        ! overwriting the input with output
        call m%linComb( ONE, ONE, m0 )
        !
        ! initialize res
        res = d
        !
        ! compute residual: res = (d-dHat)/Ndata
        !call linCombDataGroupTxArray( ONE, d, MinusONE, dHat, res )
        call subDataGroupTxArray( res, dHat )
        !
        Ndata = countValuesGroupTxArray( dHat )
        !
        !write( *, * ) "Ndata: ", Ndata
        !stop
        !
        call CdInvMult( res )
        !
        call JMult_T( m, res, JTd )
        !
        allocate( CmJTd, source = model_cov%multBy_CmSqrt( JTd ) )
        !
        ! compute the number of data and model parameters for scaling
        Nmodel = mHat%countModel()
        !
        ! multiply by 2 (to be consistent with the formula)
        ! and add the gradient of the model norm
        !
        if( allocated( grad ) ) then
            grad = CmJTd
        else
            allocate( grad, source = CmJTd )
        endif
        !
        call grad%linComb( MinusTWO / Ndata, TWO * lambda / Nmodel, mHat )
        !
        call deallocateDataGroupTxArray( res )
        !
        deallocate( m, JTd, CmJTd )
        !
    end subroutine gradient
    !
    !> Compute the full penalty functional F
    !> Also output the predicted data and the EM solution
    !> that can be used for evaluating the gradient
    subroutine func( lambda, d, m0, mHat, F, mNorm, dHat, e_all, rms )
        implicit none
        !
        real( kind=prec ), intent( in ) :: lambda
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: d
        class( ModelParameter_t ), allocatable, intent( in ) :: m0, mHat
        real( kind=prec ), intent( out ) :: F, mNorm
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: dHat
        type( ESolMTx ), optional, intent( inout ) :: e_all
        real( kind=prec ), optional, intent( out ) :: rms
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: res, Nres
        class( ModelParameter_t ), allocatable :: m, JTd
        real( kind=prec ) :: SS, angle2, angle1, diff, diff1
        integer :: Ndata, Nmodel, j, i, isite
        !
        ! compute the smoothed model parameter vector
        allocate( m, source = model_cov%multBy_CmSqrt( mHat ) )
        !
        ! overwriting input with output
        call m%linComb( ONE, ONE, m0 )
        !
        ! initialize dHat
        dHat = d
        !
        call runForwardModeling( m, dHat, e_all )
        !
        !> initialize res
        res = d
        !
        !call linCombDataGroupTxArray( ONE, d, MinusONE, dHat, res )
        call subDataGroupTxArray( res, dHat )
        !
        !> normalize residuals, compute sum of squares
        call CdInvMult( res, Nres )
        !
        SS = dotProdDataGroupTxArray( res, Nres )
        Ndata = countValuesGroupTxArray( res )
        !
        !> compute the model norm
        mNorm = mHat%dotProd( mHat )
        Nmodel = mHat%countModel()
        !
        !> penalty functional = sum of squares + scaled model norm
        F = SS / Ndata + ( lambda * mNorm / Nmodel )
        !
        !> scale mNorm for output
        mNorm = mNorm / Nmodel
        !
        ! if required, compute the Root Mean Squared misfit
        if( present( rms ) ) then
            rms = sqrt( SS / Ndata )
        endif
   
        call deallocateDataGroupTxArray( res )
        !
        call deallocateDataGroupTxArray( Nres )
        !
        deallocate( m )
        !
    end subroutine func
    !
    !> Divides by the data covariance C_d, which is a diagonal
    !> operator. Divides by the variances (squared error bars)
    !> and scales by the number of data (degrees of freedom).
    subroutine CdInvMult( d_in, d_out )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: d_in
        type( DataGroupTx_t ), allocatable, dimension(:), optional, intent( out ) :: d_out
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: d
        !
        d = d_in
        !
        call normalizeDataGroupTxArray( d, 2 )
        !
        if( present( d_out ) ) then
            d_out = d
        else
            d_in = d
        endif
        !
        call deallocateDataGroupTxArray( d )
        !
    end subroutine CdInvMult
    !
    !>
    subroutine update_damping_parameter( lambda, mHat, F, grad )
        implicit none
        !
        real( kind=prec ), intent( inout ) :: lambda
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
        SS = F - ( lambda * mNorm / Nmodel )
        !
        ! initialize
        allocate( dSS, source = mHat )
        !
        !> subtract the model norm derivative from the gradient of the penalty functional
        dSS = grad
        !
        call dSS%linComb( ONE, MinusTWO * lambda / Nmodel, mHat )
        !
        ! update the damping parameter lambda
        lambda = lambda / NLCGiterControl%k
        !
        ! penalty functional = (scaled) sum of squares + scaled model norm
        F = SS + ( lambda * mNorm / Nmodel )
        !
        ! add the model norm derivative to the gradient of the penalty functional
        grad = dSS
        !
        call grad%linComb( ONE, TWO * lambda / Nmodel, mHat )
        !
        deallocate( dSS )
        !
    end subroutine update_damping_parameter
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
    !> using the information f(0), f'(0) and f(alpha_1) to obtain
    !> a = (f(alpha_1) - f(0) - f'(0) alpha_1)/(alpha_1 * alpha_1),
    !> b = f'(0).
    !> Then, the minimum point of the quadratic is alpha_q = -b/(2a),
    !> assuming that a > 0. If this try is not successful, fit a cubic
    !>     f_c(alpha) = a alpha^3 + b alpha^2 + f'(0) alpha + f(0)
    !> using f(0), f'(0), f(alpha_1) and f(alpha_q). Repeat as necessary.
    !> Here, a and b are as described in the code.
    !> A new cubic is not identical to a previous curve since f_c is only
    !> an approximation to f: in general, f(alpha_c) /= f_c(alpha_c),
    !> hence the new point does not lie on the approximating curve.
    !
    !> Our solution has to satisfy the sufficient decrease condition
    !>     f(alpha) < f(0) + c alpha f'(0).
    !
    !> The optional relaxation parameter gamma is needed for algorithms
    !> like the Renormalised Steepest Descent (RSD). See the dynamical
    !> systems in optimisation research (Pronzato et al [2000, 2001]).
    !> To the best of my knowledge, it is not useful for NLCG.
    subroutine lineSearchCubic( lambda, d, m0, h, alpha, mHat, f, grad, &
    rms, niter, dHat, e_all, gamma )
        implicit none
        !
        real( kind=prec ), intent( in ) :: lambda
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: d
        class( ModelParameter_t ), allocatable, intent( in ) :: m0, h  ! search direction
        real( kind=prec ), intent(inout) :: alpha ! step size
        class( ModelParameter_t ), allocatable, intent( inout ) :: mHat
        real( kind=prec ), intent(inout) :: f
        class( ModelParameter_t ), allocatable, intent( inout ) :: grad
        real( kind=prec ), intent( out ) :: rms
        integer, intent( out ) :: niter
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: dHat
        type( ESolMTx ), intent( inout ) :: e_all
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
        type( ESolMTx ) :: eAll_1
        !
        ! parameters
        c = NLCGiterControl%c
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
        ! g_0 is the directional derivative f'(0) = (df/dm).dot.h
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
        eAll_1%SolnIndex = 1
        !
        call func( lambda, d, m0, mHat_1, f_1, mNorm_1, dHat_1, eAll_1, rms_1 )
        !
        write( *, * ) "lambda, alpha, f_1, mNorm_1, rms_1:", lambda, alpha, f_1, mNorm_1, rms_1
        !
        niter = niter + 1
        !
        if ( f_1 - f_0 >= LARGE_REAL ) then
            write( *, * ) "Error: Try a smaller starting r_value of alpha."
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
            !
            alpha = alpha_1
            !
            dHat = dHat_1
            !
            e_all = eAll_1
            !
            e_all%SolnIndex = 1
            !
            mHat = mHat_1
            !
            rms = rms_1
            !
            f = f_1
            !
            ! compute the gradient and exit
            if( relaxation ) then
                !
                mHat = mHat_0
                !
                call mHat%linComb( ONE, gamma * alpha, h )
                !
                e_all%SolnIndex = 0
                !
                call func( lambda, d, m0, mHat, f, mNorm, dHat, e_all, rms )
                !
                write( *, * ) "lambda, gamma*alpha, f, mNorm, rms:", lambda, gamma * alpha, f, mNorm, rms
                !
            endif
            !
            call gradient( lambda, d, m0, mHat, grad, dHat, e_all )
            !
            write( *, * ) "Quadratic has no minimum, exiting line search"
            !
            call deallocateDataGroupTxArray( dHat_1 )
            !
            deallocate( mHat_0, mHat_1 )
            !
            !call deall_solnVectorMTX( eAll_1 )
            !
            return
            !
        endif
        !
        ! otherwise compute the functional at the minimizer of the quadratic
        alpha = - b / ( TWO * a )
        !
        mHat = mHat_0
        !
        call mHat%linComb( ONE, alpha, h )
        !
        e_all%SolnIndex = 0
        !
        call func( lambda, d, m0, mHat, f, mNorm, dHat, e_all, rms )
        !
        write( *, * ) "QUADLS: lambda, alpha, f, mNorm, rms:", lambda, alpha, f, mNorm, rms
        !
        niter = niter + 1
        !
        ! check whether the solution satisfies the sufficient decrease condition
        if( f < f_0 + c * alpha * g_0 ) then
            !
            ! if the initial guess was better than what we found, take it
            if ( f_1 < f ) then
                starting_guess = .TRUE.
                alpha = alpha_1
                dHat = dHat_1
                e_all = eAll_1
                e_all%SolnIndex = 1
                mHat = mHat_1
                rms = rms_1
                f = f_1
            endif
            !
            ! compute the gradient and exit
            if( relaxation ) then
                !
                mHat = mHat_0
                !
                call mHat%linComb( ONE, gamma * alpha, h )
                !
                e_all%SolnIndex = 0
                !
                call func( lambda, d, m0, mHat, f, mNorm, dHat, e_all, rms )
                !
                write( *, * ) "QUADLS: lambda, gamma*alpha, f, mNorm, rms:", lambda, gamma*alpha, f, mNorm, rms
                !
            endif
            !
            call gradient( lambda, d, m0, mHat, grad, dHat, e_all )
            !
            write( *, * ) "Sufficient decrease condition satisfied, exiting line search"
            !
            call deallocateDataGroupTxArray( dHat_1 )
            !
            deallocate( mHat_0, mHat_1 )
            !
            !call deall_solnVectorMTX(eAll_1)
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
            !
        else
            !
            ! fit a cubic and backtrack (initialize)
            alpha_i = alpha_1
            f_i = f_1
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
                mHat = mHat_0
                !
                call mHat%linComb( ONE, alpha, h )
                !
                e_all%SolnIndex = 0
                !
                call func( lambda, d, m0, mHat, f, mNorm, dHat, e_all, rms )
                !
                write( *, * ) "CUBICLS: lambda, alpha, f, mNorm, rms:", lambda, alpha, f, mNorm, rms
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
                !
                f_i = f_j
                !
                alpha_j = alpha
                !
                f_j = f
                !
                ! check that the function still decreases to avoid infinite loops in case of a bug
                if( abs( f_j - f_i ) < TOL8 ) then
                    !
                    write( *, * ) "Warning: exiting cubic search since the function no longer decreases!"
                    !
                    exit
                    !
                endif
                !
            end do fit_cubic
            !
        endif
        !
        if( f_1 < f ) then
            starting_guess = .TRUE.
        endif
        !
        ! if the initial guess was better than what we found, take it
        if( starting_guess ) then
            alpha = alpha_1
            !
            dHat = dHat_1
            !
            e_all = eAll_1
            !
            e_all%SolnIndex = 1
            !
            mHat = mHat_1
            !
            rms = rms_1
            !
            f = f_1
        endif
        !
        ! compute gradient of the full penalty functional and exit
        if( relaxation ) then
            !
            mHat = mHat_0
            !
            call mHat%linComb( ONE, gamma*alpha, h )
            !
            e_all%SolnIndex = 0
            !
            call func( lambda, d, m0, mHat, f,mNorm,dHat,e_all,rms)
            !
            write( *, * ) "RELAX: lambda, gamma*alpha, f, mNorm, rms:", lambda, gamma*alpha, f, mNorm, rms
            !
        endif
        !
        call gradient( lambda, d, m0, mHat, grad, dHat, e_all )
        !
        write( *, * ) "Gradient computed, line search finished"
        !
        call deallocateDataGroupTxArray(dHat_1)
        !
        deallocate( mHat_0, mHat_1 )
        !
        !call deall_solnVectorMTX(eAll_1)
        !
    end subroutine lineSearchCubic
    !
    !> ????
    subroutine outputFiles_NLCG( NLCG_iter, all_predicted_data, res, m, mHat )
        implicit none
        !
        integer, intent( in ) :: NLCG_iter
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_predicted_data, res
        class( ModelParameter_t ), intent( in ) :: m, mHat
        !
        character(100) :: out_file_name
        character(8) str_date
        character(6) str_time
        character(3) :: char3
        !
        write( char3, "(i3.3)" ) NLCG_iter
        !
        !> Write predicted data for this NLCG iteration
        out_file_name = trim( outdir_name )//"/PredictedData_NLCG_"//char3//".dat"
        !
        call writeDataGroupTxArray( all_predicted_data, trim( out_file_name ) )
        !
        !> Write residual data for this NLCG iteration
        out_file_name = trim( outdir_name )//"/ResidualData_NLCG_"//char3//".dat"
        !
        call writeDataGroupTxArray( res, trim( out_file_name ) )
        !
        !> Write model for this NLCG iteration
        out_file_name = trim( outdir_name )//"/SigmaModel_NLCG_"//char3//".rho"
        !
        call m%write( trim( out_file_name ) )
        !
        !> Write perturbation model for this NLCG iteration
        out_file_name = trim( outdir_name )//"/PerturbationModel_NLCG_"//char3//".rho"
        !
        call mHat%write( trim( out_file_name ) )
        !
    end subroutine outputFiles_NLCG
    !
end module InversionNLCG
!