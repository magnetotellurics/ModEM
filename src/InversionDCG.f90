!
!> Module InversionDCG 
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
        real( kind=prec ), allocatable, dimension(:) :: r_err
        !
    end type IterControl_t
    !
contains
    !
    !>
    subroutine setDCGiterControl( DCGiterControl )
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
    end subroutine setDCGiterControl
    !
    !>
    subroutine setIterControl( cg_iter )
        implicit none
        !
        type( IterControl_t), intent( inout ) :: cg_iter
        !
        cg_iter%max_it = 20
        !
        cg_iter%tol = 10E-4
        !
        cg_iter%n_iter = 0
        !
        allocate( cg_iter%r_err( cg_iter%max_it ) )
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
        real( kind=prec ) :: F, mNorm
        type( DataGroupTx_t ), allocatable, dimension(:) :: b, dx, d_pred, res, JmHat
        class( ModelParameter_t ), allocatable :: mHat
        real( kind=prec ) :: rms
        type( IterControl_t ) :: CGiter
        integer :: DCG_iter, ios
        !
        !>
        call createOutputDirectory()
        !
        !>
        call setDCGiterControl( DCGiterControl )
        !
        !>
        call setIterControl( CGiter )
        !
        !>
        open( unit = ioInvLog, file = trim( outdir_name )//"/DCG.log", status="unknown", position="append", iostat=ios )
        !
        if( ios == 0 ) then
            !
            allocate( mHat, source = m )
            !
            call model_cov%multBy_Cm( m ) 
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
            !> Write in DCG.log
            write( ioInvLog, "( a41, es12.5 )" ) "The initial damping parameter lambda is ", DCGiterControl%lambda
            !
            call Calc_FWD( DCGiterControl%lambda, d, m, mHat, d_pred, res, F, mNorm, rms )
            !
            !> Write in DCG.log
            write( ioInvLog, "( a10, a3, es12.5, a4, es12.5, a5, f12.5, a8, es12.5 )" ) "START:", " f=", f, " m2=", mNorm, " rms=", rms, " lambda=", DCGiterControl%lambda
            !
            DCG_iter = 1
            !
            !> Print
            write( *, "( a38, f12.5)" ) "            Start_DCG : Residual rms=", rms
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
                call CG_DS_standard( b, dx, m, d, DCGiterControl%lambda, CGiter )
                !
                call normalizeWithDataGroupTxArray( 1, d, dx )
                !
                call JMult_T( m, dx, mHat )
                !
                call model_cov%multBy_Cm( mHat )
                !
                m = m0
                !
                call m%linComb( ONE, ONE, mHat )
                !
                call Calc_FWD( DCGiterControl%lambda, d, m, mHat, d_pred, res, F, mNorm, rms )
                !
                call outputFiles_DCG( DCG_iter, d_pred, res, m, mHat )
                !
                !> Write / Print DCG.log
                write( *, "( a20, i5, a16, f12.5)" ) "            DCG_iter", DCG_iter, ": Residual rms=", rms
                !
                write( ioInvLog, "( a25, i5 )" ) "Completed DCG iteration ", DCG_iter
                write( ioInvLog, "( a10, a3, es12.5, a4, es12.5, a5, f12.5, a8, es12.5 )" ) "with:", " f=", f, " m2=", mNorm, " rms=", rms, " lambda=", DCGiterControl%lambda
                !
                !>
                if( rms .LT. DCGiterControl%rmsTol .OR. DCG_iter .GE. DCGiterControl%maxIter ) then
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
            deallocate( mHat )
            !
        else
            !
            write( *, * ) "Error opening [", trim( outdir_name )//"/DCG.log", "] in writeDataGroupTxArray!"
            stop
            !
        endif
        !
    end subroutine DCGsolver
    !
    !>
    subroutine Calc_FWD( lambda, d, m, mHat, d_pred, res, F, mNorm, rms )
        implicit none
        !
        real( kind=prec ), intent( in ) :: lambda
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: d
        class( ModelParameter_t ), allocatable, intent( in ) :: m, mHat
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: d_pred
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: res
        real( kind=prec ), intent( out ) :: F, mNorm
        real( kind=prec ), intent( inout ) :: rms
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: Nres
        real( kind=prec ) :: SS
        integer :: Ndata, Nmodel
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
        Ndata = countValuesGroupTxArray( res )
        !
        mNorm = mHat%dotProd( mHat )
        !
        Nmodel = mHat%countModel()
        !
        F = SS / Ndata + ( lambda * mNorm / Nmodel )
        !
        rms = sqrt( SS / Ndata )
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
        !> Write / Print DCG.log
        write( ioInvLog, "(a18)" ) "Relative CG-error:"
        write( ioInvLog, "( a9, i5, a10, es12.5, a10, es12.5 )" ) "CG-Iter= ", cg_iter, ", error = ", CGiter%r_err(cg_iter), " Lambda= ", lambda
        !
        write( *, "( a22, i5, a8, es12.5 )" ) "               CG_iter", cg_iter, ": Error=", CGiter%r_err(cg_iter)
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
            !> Write / Print DCG.log
            write( ioInvLog, "( a9, i5, a10, es12.5, a10, es12.5 )" ) "CG-Iter= ", cg_iter, ", error = ", CGiter%r_err(cg_iter), " Lambda= ", lambda
            !
            write( *, "( a22, i5, a8, es12.5, a7, es12.5, a8, es12.5 )" ) "               CG_iter", cg_iter, ": Alpha=", alpha, ", Beta=", beta, ", Error=", CGiter%r_err( cg_iter )
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
        class( ModelParameter_t ), allocatable :: JTp
        type( DataGroupTx_t ), allocatable, dimension(:) :: lambdaP, p_temp
        !
        p_temp = p
        !
        lambdaP = p
        !
        call normalizeWithDataGroupTxArray( 1, d, p_temp )
        !
        call JMult_T( m, p_temp, JTp )
        !
        call model_cov%multBy_Cm( JTp )
        !
        Ap = d
        !
        call JMult( m, JTp, Ap )
        !
        deallocate( JTp )
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
    !> ????
    subroutine outputFiles_DCG( DCG_iter, d_pred, res, m, mHat )
        implicit none
        !
        integer, intent( in ) :: DCG_iter
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: d_pred, res
        class( ModelParameter_t ), intent( in ) :: m, mHat
        !
        character(100) :: out_file_name
        character(8) str_date
        character(6) str_time
        character(3) :: char3
        !
        write( char3, "(i3.3)" ) DCG_iter
        !
        !> Write predicted data for this DCG iteration
        out_file_name = trim( outdir_name )//"/PredictedData_DCG_"//char3//".dat"
        !
        call writeDataGroupTxArray( d_pred, trim( out_file_name ) )
        !
        !> Write residual data for this DCG iteration
        out_file_name = trim( outdir_name )//"/ResidualData_DCG_"//char3//".dat"
        !
        call writeDataGroupTxArray( res, trim( out_file_name ) )
        !
        !> Write model for this DCG iteration
        out_file_name = trim( outdir_name )//"/SigmaModel_DCG_"//char3//".rho"
        !
        call m%write( trim( out_file_name ) )
        !
        !> Write perturbation model for this DCG iteration
        out_file_name = trim( outdir_name )//"/PerturbationModel_DCG_"//char3//".rho"
        !
        call mHat%write( trim( out_file_name ) )
        !
    end subroutine outputFiles_DCG
    !
    !> ????
    subroutine createOutputDirectory()
        implicit none
        !
        character(8) str_date
        character(6) str_time
        !
        !>
        if( .NOT. has_outdir_name ) then
            !
            !>
            call date_and_time( str_date, str_time )
            !
            write( outdir_name, "(a11, a8, a1, a6)" ) "DCG_Output_", str_date, "_", str_time
            !
            !write( *, * ) "outdir_name: [", trim( outdir_name ), "]"
            !
        endif
        !
        call system( "mkdir -p "//outdir_name )
        !
    end subroutine createOutputDirectory
    !
end module InversionDCG
!