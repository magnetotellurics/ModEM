!
!> Module with the InversionDCG routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module InversionDCG
    !
    use ForwardModeling
    use Sensitivity
    !
    !>
    type  :: DCGiterControl_t
        !
        integer :: maxIter
        !
        real( kind=prec ) :: rmsTol, lambda
        !
    end type DCGiterControl_t
    !
    type( DCGiterControl_t ), private, save :: DCGiterControl
    !
    !>
    type :: IterControl_t
        !
        integer :: max_it, n_iter
        !
        real( kind=prec ) :: tol
        !
        real( kind=prec ), pointer, dimension(:) :: r_err
        !
    end type IterControl_t
    !
contains
    !
    !>
    subroutine set_DCGiterControl( DCGiterControl )
        implicit none
        !
        type( DCGiterControl_t ), intent( inout ) :: DCGiterControl
        !
        DCGiterControl%maxIter = 3
        !
        DCGiterControl%rmsTol = 1.05
        !
        DCGiterControl%lambda = 10.
        !
    end subroutine set_DCGiterControl
    !
    !>
    subroutine setIterControl( cg_iter )
        implicit none
        !
        type( IterControl_t), intent( inout ) :: cg_iter
        !
        cg_iter%max_it = 3
        !
        cg_iter%tol = 10E-4
        !
        cg_iter%n_iter = 0
        !
        allocate( cg_iter%r_err( 0 : cg_iter%max_it ) )
        !
    end subroutine setIterControl
    !
    !>
    subroutine DCGsolver( d, m0, m )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: d
        class( ModelParameter_t ), allocatable, intent( in ) :: m0
        class( ModelParameter_t ), allocatable, intent( inout ) :: m
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: b, dx, d_pred, res, JmHat
        class( ModelParameter_t ), allocatable :: mHat, Cm_mHat
        real( kind=prec ) :: rmsd, lambda
        type( IterControl_t ) :: CGiter
        integer :: DCG_iter
        !
        !>
        call set_DCGiterControl( DCGiterControl )
        !
        call setIterControl( CGiter )
        !
        lambda = 10.
        !
        allocate( mHat, source = m )
        !
        allocate( Cm_mHat, source = m )
        !
        m = model_cov%multBy_Cm( mHat ) 
        !
        call m%linComb( ONE, ONE, m0 )
        !
        JmHat = d
        !
        dx = d
        !
        b = d
        !
        call zerosDataGroupTxArray( JmHat )
        !
        call zerosDataGroupTxArray( b )
        !
        call Calc_FWD( d, m, d_pred, res, rmsd )
        !
        DCG_iter = 1
        !
        write( *, * ) "#START DCG LOOP > residual rmsd: ", rmsd
        !
        write( 8888, * ) "DCG_iter,    rmsd,    rmsTol,    maxIter"
        !
        write( 8888, * ) DCG_iter, rmsd, DCGiterControl%rmsTol, DCGiterControl%maxIter
        !
        !>
        dcg_loop : do
            !
            JmHat = d
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
            call Calc_FWD( d, m, d_pred, res, rmsd )
            !
            write( *, * ) "#DCG_iter > residual rmsd: ", DCG_iter, rmsd
			!
			write( 8888, * ) DCG_iter, rmsd
			!
            !>
            if( rmsd .LT. DCGiterControl%rmsTol .OR. DCG_iter .GE. DCGiterControl%maxIter ) then
                exit
            end if
            !
            DCG_iter = DCG_iter + 1
            !
        end do dcg_loop
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
    subroutine Calc_FWD( d, m, d_pred, res, rmsd )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: d
        class( ModelParameter_t ), allocatable, intent( in ) :: m
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: d_pred
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
        write( *, * ) "     #Calc_FWD > SS, Ndata, rmsd: ", SS, Ndata, rmsd
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
        CGiter%r_err(cg_iter) = r_norm / b_norm
        !
        write( *, * ) "     #START CG > Error, Lambda: ", cg_iter, CGiter%r_err(cg_iter), lambda
        !
        write( 8888, * ) "ii,    alpha,    beta,    b_norm,    r_norm,    CGiter%rerr(ii),    lambda"
        !
        write( 8888, * ) cg_iter, alpha, beta, b_norm, r_norm, CGiter%r_err(cg_iter), lambda
        !
        cg_loop : do while ( CGiter%r_err(cg_iter) .GT. CGiter%tol .AND. cg_iter .LT. CGiter%max_it )
            ! 
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
            !
            r_norm = dotProdDataGroupTxArray( r, r )
            !
            beta = r_norm / r_norm_pre
            !
            ! Compute new p: p = r + beta*p    
            call linCombDataGroupTxArray( ONE, r, beta, p, p )
            !
            cg_iter = cg_iter + 1
            !
            CGiter%r_err(cg_iter) = r_norm / b_norm 
            !
            write( *, * ) "       #CG-Iter, Error, Lambda: ", cg_iter, CGiter%r_err(cg_iter), lambda
            !
            write( 8888, * ) cg_iter, alpha, beta, b_norm, r_norm, CGiter%r_err(cg_iter), lambda
            !
        enddo cg_loop
        !
        CGiter%n_iter = cg_iter
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
        call linCombDataGroupTxArray( R_ZERO, d, ONE, p_temp, p_temp ) 
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
        call scMultDataGroupTxArray( lambda, p, lambdaP )
        !
        call setErrorBarDataGroupTxArray( lambdaP, .FALSE. )
        !
        call linCombDataGroupTxArray( ONE, Ap, ONE, lambdaP, Ap )
        !
    end subroutine MultA_DS
    !
end module InversionDCG
!