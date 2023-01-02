!
!> Module with the InversionDCG routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module InversionDCG
    !
    use Sensitivity
    !
    type  :: IterControl_t
        !
        ! maximum number of iterations in one call to iterative solver
        integer :: maxIt
        ! convergence criteria: return from solver if relative error < tol
        real( kind=prec ) :: tol
        ! actual number of iterations before return
        integer :: niter
        ! relative error for each iteration
        real( kind=prec ) , pointer, dimension(:) :: rerr
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
        real(kind=prec) :: alpha, beta, r_norm_pre, r_norm, b_norm, error
        integer :: i, j, k, ii, iDt
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
        ii = 1
        CGiter%rerr(ii) = r_norm / b_norm
        !
        loop: do while ( CGiter%rerr(ii) .GT. CGiter%tol .AND. ii .LT.CGiter%maxIt )
            !
            ! Compute matrix-vector product A*p and save the result in Ap  
            call MultA_DS( p, m, d, lambda, Ap )
            !
            do i = 1, size( x )
                do iDt = 1, size( x(i)%data )
                    r(i)%data(iDt)%error_bar= .FALSE.
                    p(i)%data(iDt)%error_bar= .FALSE.
                    x(i)%data(iDt)%error_bar= .FALSE.
                    Ap(i)%data(iDt)%error_bar= .FALSE.
                enddo
            enddo
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
            ii = ii + 1
            !
            CGiter%rerr(ii) = r_norm / b_norm 
            !
            write( *, * ) "CG-Iter, error, Lambda: ", ii, CGiter%rerr(ii), lambda
            !
        enddo loop
        !
        CGiter%niter = ii
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
        !call normalize_with_dataVecMTX( p_temp, d, 1)
        call normalizeWithDataGroupTxArray( 1, d, p_temp )
        ! Compute   J^T  Cd^(-1/2) p
        !
        call JMult_T( m, p_temp, JTp )
        !
        ! Compute  Cm  J^T  Cd^(-1/2) p 
        CmJTp_temp = model_cov%multBy_Cm( JTp )
        !
        deallocate( JTp )
        !
        ! Compute J Cm  J^T  Cd^(-1/2) p = Ap 
        Ap = d
        !
        call JMult( m, CmJTp_temp, Ap )
        !
        deallocate( CmJTp_temp )
        !
        !> Normalize: Cd^(-1/2)*Ap
        !call normalize_with_dataVecMTX(Ap,d,1)
        call normalizeWithDataGroupTxArray( 1, d, Ap )
        !
        call scMultAddDataGroupTxArray( lambda, p, lambdaP )
        !
        !> Add Cd^(-1/2)*Ap*Cd^(-1/2) to lambda*p
        do i = 1, size( lambdaP )
            do iDt = 1, size( lambdaP(i)%data )
                lambdaP(i)%data(iDt)%error_bar= .FALSE.
            enddo
        enddo
        !
        call linCombDataGroupTxArray( ONE, Ap, ONE, lambdaP, Ap )
        !
    end subroutine MultA_DS
    !
end module InversionDCG
!