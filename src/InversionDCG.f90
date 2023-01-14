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
    subroutine DCGsolver( all_data, sigma, dsigma )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_data
        class( ModelParameter_t ), allocatable, intent( in ) :: sigma
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        real( kind=prec ) :: F, mNorm
        type( DataGroupTx_t ), allocatable, dimension(:) :: b, dx, all_predicted_data, res, JmHat
        class( ModelParameter_t ), allocatable :: mHat
        real( kind=prec ) :: rms
        type( IterControl_t ) :: CGiter
        integer :: DCG_iter, ios
        !
        !>
        call createOutputDirectory()
        !
        ! Verbose
        write( *, * ) "     - Start Inversion DCG, output files in [", trim( outdir_name ), "]"
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
            allocate( mHat, source = dsigma )
            !
            call model_cov%multBy_Cm( dsigma ) 
            !
            call dsigma%linComb( ONE, ONE, sigma )
            !
            JmHat = all_data
            !
            dx = all_data
            !
            b = all_data
            !
            call zerosDataGroupTxArray( JmHat )
            !
            call zerosDataGroupTxArray( b )
            !
            !> Write in DCG.log
            write( ioInvLog, "( a41, es12.5 )" ) "The initial damping parameter lambda is ", DCGiterControl%lambda
            !
            call Calc_FWD( DCGiterControl%lambda, all_data, dsigma, mHat, all_predicted_data, res, F, mNorm, rms )
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
                JmHat = all_data
                !
                call JMult( dsigma, mHat, JmHat )
                !
                b = all_data
                !
                call linCombDataGroupTxArray( ONE, res, ONE, JmHat, b )
                !
                call normalizeDataGroupTxArray( b, 1 )
                !
                call CG_DS_standard( b, dx, dsigma, all_data, DCGiterControl%lambda, CGiter )
                !
                call normalizeWithDataGroupTxArray( 1, all_data, dx )
                !
                call JMult_T( dsigma, dx, mHat )
                !
                call model_cov%multBy_Cm( mHat )
                !
                dsigma = sigma
                !
                call dsigma%linComb( ONE, ONE, mHat )
                !
                call Calc_FWD( DCGiterControl%lambda, all_data, dsigma, mHat, all_predicted_data, res, F, mNorm, rms )
                !
                call outputFiles_DCG( DCG_iter, all_predicted_data, res, dsigma, mHat )
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
            ! Verbose
            write( *, * ) "     - Finish Inversion DCG, output files in [", trim( outdir_name ), "]"
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
    subroutine Calc_FWD( lambda, all_data, dsigma, mHat, all_predicted_data, res, F, mNorm, rms )
        implicit none
        !
        real( kind=prec ), intent( in ) :: lambda
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( in ) :: dsigma, mHat
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: all_predicted_data
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: res
        real( kind=prec ), intent( out ) :: F, mNorm
        real( kind=prec ), intent( inout ) :: rms
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: Nres
        real( kind=prec ) :: SS
        integer :: Ndata, Nmodel
        !
        all_predicted_data = all_data
        !
        call runForwardModeling( dsigma, all_predicted_data )
        !
        res = all_data
        !
        call linCombDataGroupTxArray( ONE, all_data, MinusONE, all_predicted_data, res )
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
    subroutine CG_DS_standard( b, x, dsigma, all_data, lambda, CGiter )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: b
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: x
        class( ModelParameter_t ), allocatable, intent( in ) :: dsigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_data
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
        Ap = all_data
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
        cg_loop : do while ( CGiter%r_err(cg_iter) .GT. CGiter%tol .AND. cg_iter .LE. CGiter%max_it )
            ! 
            call MultA_DS( p, dsigma, all_data, lambda, Ap )
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
    subroutine MultA_DS( p, dsigma, all_data, lambda, Ap )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: p
        class( ModelParameter_t ), allocatable, intent( in ) :: dsigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
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
        call normalizeWithDataGroupTxArray( 1, all_data, p_temp )
        !
        call JMult_T( dsigma, p_temp, JTp )
        !
        call model_cov%multBy_Cm( JTp )
        !
        Ap = all_data
        !
        call JMult( dsigma, JTp, Ap )
        !
        deallocate( JTp )
        !
        call normalizeWithDataGroupTxArray( 1, all_data, Ap )
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
    subroutine outputFiles_DCG( DCG_iter, all_predicted_data, res, dsigma, mHat )
        implicit none
        !
        integer, intent( in ) :: DCG_iter
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_predicted_data, res
        class( ModelParameter_t ), intent( in ) :: dsigma, mHat
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
        call writeDataGroupTxArray( all_predicted_data, trim( out_file_name ) )
        !
        !> Write residual data for this DCG iteration
        out_file_name = trim( outdir_name )//"/ResidualData_DCG_"//char3//".dat"
        !
        call writeDataGroupTxArray( res, trim( out_file_name ) )
        !
        !> Write model for this DCG iteration
        out_file_name = trim( outdir_name )//"/SigmaModel_DCG_"//char3//".rho"
        !
        call dsigma%write( trim( out_file_name ) )
        !
        !> Write perturbation model for this DCG iteration
        out_file_name = trim( outdir_name )//"/PerturbationModel_DCG_"//char3//".rho"
        !
        call mHat%write( trim( out_file_name ) )
        !
    end subroutine outputFiles_DCG
    !
end module InversionDCG
!