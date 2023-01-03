!
!> Module with the InversionDCG routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module InversionDCG
    !
    use ForwardModeling
    use Sensitivity
    !
    type :: IterControl_t
        !
        ! maximum number of iterations in one call to iterative solver
        integer :: maxIt
        ! convergence criteria: return from solver if relative error < tol
        real( kind=prec ) :: tol
        ! actual number of iterations before return
        integer :: niter
        ! relative error for each iteration
        real( kind=prec ), pointer, dimension(:) :: rerr
        ! logical variable indicating if algorithm "failed"
        logical :: failed = .FALSE.
        !
    end type IterControl_t
    !
contains
    !
    !>
    subroutine setIterControl( CGiter )
        implicit none
        !
        type( IterControl_t), intent( inout ) :: CGiter
        !
        CGiter%maxit = 20
        CGiter%tol = 10E-4
        CGiter%niter = 0
        allocate( CGiter%rerr( 0 : CGiter%maxit ) )
        !
    end subroutine setIterControl
    !
    !>
    subroutine DCGsolver( d, m0, m, lambda )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: d
        class( ModelParameter_t ), allocatable, intent( in ) :: m0
        class( ModelParameter_t ), allocatable, intent( inout ) :: m
        real( kind=prec ), intent( inout ) :: lambda
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: b, dx, d_pred, res, JmHat
        class( ModelParameter_t ), allocatable :: mHat, Cm_mHat
        real( kind=prec ) :: rmsd
        type( IterControl_t ) :: CGiter
        integer :: DCG_iter
        !
        lambda = 10.0
        !
        call setIterControl( CGiter )
        !
        allocate( mHat, source = m )
        !
        allocate( Cm_mHat, source = m )
        !
        m = model_cov%multBy_Cm( mHat ) 
        !
        !call linComb(ONE,m,ONE,m0,m)
        call m%linComb( ONE, ONE, m0 )
        !
        JmHat = d
        dx = d
        b = d
        res = d
        !
        call zerosDataGroupTxArray( JmHat )
        !
        call zerosDataGroupTxArray( b )
        !
        d_pred = d
        !
        call Calc_FWD( lambda, d, m, d_pred, res, rmsd )
        !
        DCG_iter = 1
        !
        write( *, * ) "#START DCG LOOP > rmsd: ", rmsd
        !
        do
            !
            call JMult( m, mHat, JmHat )
            !
            b = d
            !
            call linCombDataGroupTxArray( ONE, res, ONE, JmHat, b )
            !
            call normalizeDataGroupTxArray( b, 1 )
            !
            call CG_DS_standard( b, dx, m, d, lambda, CGiter )
            !
            call normalizeWithDataGroupTxArray( 1, d, dx )
            !
            call JMult_T( m, dx, mHat )
            !
            Cm_mHat = model_cov%multBy_Cm( mHat )
            !
            mHat = Cm_mHat
            !
            m = m0
            !
            call m%linComb( ONE, ONE, mHat )
            !
            call Calc_FWD( lambda, d, m, d_pred, res, rmsd )
            !
            write( *, * ) "#DCG_iter > rmsd: ", DCG_iter, rmsd
            !
            if( rmsd .LT. 1.05 .OR. DCG_iter .GE. 3 ) then
                exit
            endif
            !
            DCG_iter = DCG_iter + 1
            !
        end do
        !
        call deallocateDataGroupTxArray( JmHat )
        call deallocateDataGroupTxArray( b )
        call deallocateDataGroupTxArray( res )
        !
        deallocate( mHat, Cm_mHat )
        !
    end subroutine DCGsolver
    !
    !>
    subroutine DCGsolverLanczos( d, m0, m, lambda )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: d
        class( ModelParameter_t ), allocatable, intent( in ) :: m0
        class( ModelParameter_t ), allocatable, intent( inout ) :: m
        real( kind=prec ), intent( inout ) :: lambda
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: JmHat, b, dx, d_pred, res
        class( ModelParameter_t ), allocatable :: mHat, Cm_mHat
        real( kind=prec ) :: rmsd
        integer :: DS_iter
        type( IterControl_t ) :: CGiter
        !
        allocate( mHat, source = m )
        allocate( Cm_mHat, source = m )
        !
        m = model_cov%multBy_Cm( mHat )
        !
        call m%linComb( ONE, ONE, m0 )
        !
        JmHat = d
        dx = d
        b = d
        res = d
        !
        call zerosDataGroupTxArray( JmHat )
        !
        call zerosDataGroupTxArray( b )
        !
        call setIterControl( CGiter )
        !
        d_pred = d
        !
        call Calc_FWD( lambda, d, m, d_pred, res, rmsd )
        !
        write( *, * ) "lambda, rmsd: ", lambda, rmsd
        !
        do DS_iter = 1, 5
            !
            ! Compute the right hand side vector (b) for the CG solver.
            ! b= (d-dPred)+ J(m-m0)
            !
            if ( DS_iter .GT. 1 ) then
                call Jmult( m, mHat, JmHat )
            endif
            !
            b = d
            write( *, * ) "Norm JmHat: ", dotProdDataGroupTxArray( JmHat, JmHat )
            !
            call linCombDataGroupTxArray( ONE, res, ONE, JmHat, b )
            !
            call normalizeDataGroupTxArray( b, 1 )
            !
            !  call CG_DS(b,dx,m,m0,d,lambda,CGiter,d_Pred_m0,rms,mhat)
            !call Lanczos_DS (b,m,m0,d,d_Pred_m0,lambda,mhat,res,CGiter,DS_iter,rms,Jm0)
            !  call Multi_Trans_DS (b,dx,m,m0,d,lambda,mhat,res,CGiter,DS_iter,rms)
            ! call Lanczos_CG_DS(b,dx,m,m0,d,lambda,CGiter,DS_iter,rms)
            call CG_DS_standard( b, dx, m, d, lambda, CGiter )
            !
            call normalizeWithDataGroupTxArray( 1, d, dx )
            !
            call JMult_T( m, dx, mHat )
            !
            Cm_mHat = model_cov%multBy_Cm( mHat )
            !
            mHat = Cm_mHat
            !
            !call linComb_modelParam( ONE, m0, ONE, mHat, m )
            m = m0
            call m%linComb( ONE, ONE, mHat )
            !
            call Calc_FWD( lambda, d, m, d_pred, res, rmsd )
            !
            write( *, * ) "DS_iter, lambda, rmsd, CGiter%niter: ", DS_iter, lambda, rmsd, CGiter%niter
            !
        end do
        !
        deallocate( mHat )
        deallocate( Cm_mHat )
        !
    end subroutine DCGsolverLanczos
    !
    !>
    subroutine Calc_FWD( lambda, d, m, d_pred, res, rmsd )
        implicit none
        !
        real( kind=prec ), intent( in ) :: lambda
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: d
        class( ModelParameter_t ), allocatable, intent( in ) :: m
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: d_pred
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: res
        real( kind=prec ), intent( inout ) :: rmsd
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: Nres
        real( kind=prec ) :: SS
        integer :: Ndata
        !
        d_pred = d
        !
        call runForwardModeling( m, d_pred )
        !
        res = d
        !
        call linCombDataGroupTxArray( ONE, d, MinusONE, d_pred, res )
        !
        Nres = res
        !
        call normalizeDataGroupTxArray( Nres, 2 )
        !
        SS = dotProdDataGroupTxArray( res, Nres )
        !
        Ndata = countDataGroupTxArray( res )
        !
        rmsd = sqrt( SS / Ndata )
        !
        call deallocateDataGroupTxArray( Nres )
        !
    end subroutine Calc_FWD
    !
    !>
    subroutine CG_DS_standard( b, x, m, d, lambda, CGiter )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: b
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: x
        class( ModelParameter_t ), allocatable, intent( in ) :: m
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: d
        real( kind=prec ), intent( in ) :: lambda
        type( IterControl_t ), intent( inout ) :: CGiter
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: r, p, Ap
        real( kind=prec ) :: alpha, beta, r_norm_pre, r_norm, b_norm
        integer :: cg_iter
        !
        r = b
        !
        p = r
        !
        Ap = d
        !
        b_norm = dotProdDataGroupTxArray( b, b )
        !
        call zerosDataGroupTxArray( x )
        !
        r_norm = dotProdDataGroupTxArray( r, r )
        !
        cg_iter = 1
        !
        CGiter%rerr(cg_iter) = r_norm / b_norm
        !
        loop: do while ( CGiter%rerr(cg_iter) .GT. CGiter%tol .AND. cg_iter .LT.CGiter%maxIt )
            !
            write( *, * ) "#CG-Iter, Error, Lambda: ", cg_iter, CGiter%rerr(cg_iter), lambda
            !
            ! Compute matrix-vector product A*p and save the result in Ap  
            call MultA_DS( p, m, d, lambda, Ap )
            !
            call setErrorBarDataGroupTxArray( r, .FALSE. )
            call setErrorBarDataGroupTxArray( p, .FALSE. )
            call setErrorBarDataGroupTxArray( x, .FALSE. )
            call setErrorBarDataGroupTxArray( Ap, .FALSE. )
            !
            ! Compute alpha: alpha= (r^T r) / (p^T Ap)    
            alpha = r_norm / dotProdDataGroupTxArray( p, Ap )
            !
            ! Compute new x: x = x + alpha*p
            call scMultAddDataGroupTxArray( alpha, p, x )
            !
            ! Compute new r: r = r - alpha*Ap
            call scMultAddDataGroupTxArray( -alpha, Ap, r ) 
            !
            r_norm_pre = r_norm
            r_norm = dotProdDataGroupTxArray( r, r )
            !
            ! Compute beta: beta= r_norm /r_norm_previous
            beta = r_norm / r_norm_pre
            
            ! Compute new p: p = r + beta*p    
            call linCombDataGroupTxArray( ONE, r, beta, p, p )
            !
            cg_iter = cg_iter + 1
            !
            CGiter%rerr(cg_iter) = r_norm / b_norm 
            !
        enddo loop
        !
        CGiter%niter = cg_iter
        !
    end subroutine CG_DS_standard
    !
    !> ????
    subroutine MultA_DS( p, m, d, lambda, Ap )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: p
        class( ModelParameter_t ), allocatable, intent( in ) :: m
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: d
        real( kind=prec ), intent( in ) :: lambda
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: Ap
        !
        class( ModelParameter_t ), allocatable :: JTp, CmJTp_temp
        type( DataGroupTx_t ), allocatable, dimension(:) :: lambdaP, p_temp
        integer :: i, j, k, iDt
        !
        allocate( JTp, source = m )
        !
        allocate( CmJTp_temp, source = m )
        !
        p_temp = p
        !
        lambdaP = p
        !
        call normalizeWithDataGroupTxArray( 1, d, p_temp )
        !
        call JMult_T( m, p_temp, JTp )
        !
        CmJTp_temp = model_cov%multBy_Cm( JTp )
        !
        deallocate( JTp )
        !
        Ap = d
        !
        call JMult( m, CmJTp_temp, Ap )
        !
        deallocate( CmJTp_temp )
        !
        call normalizeWithDataGroupTxArray( 1, d, Ap )
        !
        p_temp = p
        !
        call setErrorBarDataGroupTxArray( p_temp, .FALSE. )
        !
        call scMultAddDataGroupTxArray( lambda, p_temp, lambdaP )
        !
        call setErrorBarDataGroupTxArray( lambdaP, .FALSE. )
        !
        call linCombDataGroupTxArray( ONE, Ap, ONE, lambdaP, Ap )
        !
    end subroutine MultA_DS
    !
end module InversionDCG
!