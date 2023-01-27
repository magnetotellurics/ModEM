!
!> Derived class to define a InversionDCG
!
module InversionDCG
    !
    use Inversion
    !
    type, extends( Inversion_t ) :: InversionDCG_t
        !
        !> No derived properties
        !
        contains
            !
            final :: InversionDCG_dtor
            !
            procedure, public :: solve => solveInversionDCG
            !
            procedure, public :: outputFiles => outputFilesInversionDCG
            !
            procedure, private :: Calc_FWD, CG_DS_standard, MultA_DS
            !
    end type InversionDCG_t
    !
    interface InversionDCG_t
        module procedure InversionDCG_ctor
    end interface InversionDCG_t
    !
contains
    !
    !> No function briefing
    function InversionDCG_ctor() result( self )
        implicit none
        !
        type( InversionDCG_t ) :: self
        !
        write( *, * ) "Constructor InversionDCG_t"
        !
        call self%init()
        !
        self%max_iter = 3
        !
        self%rms_tol = 1.05
        !
        self%lambda = 10.
        !
        allocate( self%r_err( self%max_cg_iter ) )
        !
    end function InversionDCG_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    subroutine InversionDCG_dtor( self )
        implicit none
        !
        type( InversionDCG_t ), intent( inout ) :: self
        !
        write( *, * ) "Destructor InversionDCG_t"
        !
        call self%dealloc()
        !
    end subroutine InversionDCG_dtor
    !
    !>
    subroutine solveInversionDCG( self, all_data, sigma, dsigma )
        implicit none
        !
        class( InversionDCG_t ), intent( inout ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( in ) :: sigma
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        real( kind=prec ) :: F, mNorm
        type( DataGroupTx_t ), allocatable, dimension(:) :: b, dx, all_predicted_data, res, JmHat
        class( ModelParameter_t ), allocatable :: mHat
        real( kind=prec ) :: rms
        integer :: DCG_iter, ios
        !
        !>
        call createOutputDirectory()
        !
        ! Verbose
        write( *, * ) "     - Start Inversion DCG, output files in [", trim( outdir_name ), "]"
        !
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
            dx = all_data
            !
            b = all_data
            !
            call zerosDataGroupTxArray( b )
            !
            !> Write in DCG.log
            write( ioInvLog, "( a41, es12.5 )" ) "The initial damping parameter lambda is ", self%lambda
            !
            call self%Calc_FWD( all_data, dsigma, mHat, all_predicted_data, res, F, mNorm, rms )
            !
            !> Write in DCG.log
            write( ioInvLog, "( a10, a3, es12.5, a4, es12.5, a5, f18.5, a8, es12.5 )" ) "START:", " f=", f, " m2=", mNorm, " rms=", rms, " lambda=", self%lambda
            !
            DCG_iter = 1
            !
            !> Print
            write( *, "( a38, f18.5)" ) "            Start_DCG : Residual rms=", rms
            !
            !>
            dcg_loop : do
                !
#ifdef MPI
                call masterJMult( dsigma, mHat, JmHat )
#else
                call serialJMult( dsigma, mHat, JmHat, self%new_sigma )
#endif
                !
                b = all_data
                !
                call linCombDataGroupTxArray( ONE, res, ONE, JmHat, b )
                !
                call normalizeDataGroupTxArray( b, 1 )
                !
                call CG_DS_standard( self, b, dx, dsigma, all_data )
                !
                call normalizeWithDataGroupTxArray( 1, all_data, dx )
                !
#ifdef MPI
                call masterJMult_T( dsigma, dx, mHat )
#else
                call serialJMult_T( dsigma, dx, mHat, self%new_sigma )
#endif
                !
                call model_cov%multBy_Cm( mHat )
                !
                dsigma = sigma
                !
                call dsigma%linComb( ONE, ONE, mHat )
                !
                call self%Calc_FWD( all_data, dsigma, mHat, all_predicted_data, res, F, mNorm, rms )
                !
                call self%outputFiles( DCG_iter, all_predicted_data, res, dsigma, mHat )
                !
                !> Write / Print DCG.log
                write( *, "( a20, i5, a16, f18.5)" ) "            DCG_iter", DCG_iter, ": Residual rms=", rms
                !
                write( ioInvLog, "( a25, i5 )" ) "Completed DCG iteration ", DCG_iter
                write( ioInvLog, "( a10, a3, es12.5, a4, es12.5, a5, f18.5, a8, es12.5 )" ) "with:", " f=", f, " m2=", mNorm, " rms=", rms, " lambda=", self%lambda
                !
                !>
                if( rms .LT. self%rms_tol .OR. DCG_iter .GE. self%max_iter ) then
                    exit
                end if
                !
                DCG_iter = DCG_iter + 1
                !
            end do dcg_loop
            !
            !call deallocateDataGroupTxArray( JmHat )
            !call deallocateDataGroupTxArray( b )
            !call deallocateDataGroupTxArray( res )
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
    end subroutine solveInversionDCG
    !
    !>
    subroutine Calc_FWD( self, all_data, dsigma, mHat, all_predicted_data, res, F, mNorm, rms )
        implicit none
        !
        class( InversionDCG_t ), intent( inout ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( in ) :: dsigma, mHat
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_predicted_data
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
#ifdef MPI
        !
        call masterForwardModelling( dsigma, all_predicted_data )
        !
#else
        !
        call serialForwardModeling( dsigma, all_predicted_data, self%new_sigma )
        !
#endif
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
        F = SS / Ndata + ( self%lambda * mNorm / Nmodel )
        !
        rms = sqrt( SS / Ndata )
        !
        !call deallocateDataGroupTxArray( Nres )
        !
    end subroutine Calc_FWD
    !
    !>
    subroutine CG_DS_standard( self, b, x, dsigma, all_data )
        implicit none
        !
        class( InversionDCG_t ), intent( inout ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: b
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: x
        class( ModelParameter_t ), allocatable, intent( in ) :: dsigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: r, p, Ap
        real( kind=prec ) :: alpha, beta, r_norm_pre, r_norm, b_norm
        integer :: iter
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
        iter = 1
        !
        self%r_err(iter) = r_norm / b_norm
        !
        !> Write / Print DCG.log
        write( ioInvLog, "(a18)" ) "Relative CG-error:"
        write( ioInvLog, "( a9, i5, a10, es12.5, a10, es12.5 )" ) "CG-Iter= ", iter, ", error = ", self%r_err(iter), " Lambda= ", self%lambda
        !
        write( *, "( a22, i5, a8, es12.5 )" ) "               CG_iter", iter, ": Error=", self%r_err(iter)
        !
        cg_loop : do while ( self%r_err(iter) .GT. self%tol .AND. iter .LT. self%max_cg_iter )
            ! 
            call self%MultA_DS( p, dsigma, all_data, Ap )
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
            iter = iter + 1
            !
            self%r_err(iter) = r_norm / b_norm 
            !
            !> Write / Print DCG.log
            write( ioInvLog, "( a9, i5, a10, es12.5, a10, es12.5 )" ) "CG-Iter= ", iter, ", error = ", self%r_err(iter), " Lambda= ", self%lambda
            !
            write( *, "( a22, i5, a8, es12.5, a7, es12.5, a8, es12.5 )" ) "               CG_iter", iter, ": Alpha=", alpha, ", Beta=", beta, ", Error=", self%r_err( iter )
            !
        enddo cg_loop
        !
        self%n_iter = iter
        !
    end subroutine CG_DS_standard
    !
    !> ????
    subroutine MultA_DS( self, p, dsigma, all_data, Ap )
        implicit none
        !
        class( InversionDCG_t ), intent( inout ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: p
        class( ModelParameter_t ), allocatable, intent( in ) :: dsigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
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
#ifdef MPI
        call masterJMult_T( dsigma, p_temp, JTp )
#else
        call serialJMult_T( dsigma, p_temp, JTp, self%new_sigma )
#endif
        !
        call model_cov%multBy_Cm( JTp )
        !
        Ap = all_data
        !
#ifdef MPI
        call masterJMult( dsigma, JTp, Ap )
#else
        call serialJMult( dsigma, JTp, Ap, self%new_sigma )
#endif
        !
        deallocate( JTp )
        !
        call normalizeWithDataGroupTxArray( 1, all_data, Ap )
        !
        call scMultDataGroupTxArray( self%lambda, p, lambdaP )
        !
        call setErrorBarDataGroupTxArray( lambdaP, .FALSE. )
        !
        call linCombDataGroupTxArray( ONE, Ap, ONE, lambdaP, Ap )
        !
    end subroutine MultA_DS
    !
    !> ????
    subroutine outputFilesInversionDCG( self, iter, all_predicted_data, res, dsigma, mHat )
        implicit none
        !
        class( InversionDCG_t ), intent( inout ) :: self
        integer, intent( in ) :: iter
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_predicted_data, res
        class( ModelParameter_t ), intent( in ) :: dsigma, mHat
        !
        character(100) :: out_file_name
        character(8) str_date
        character(6) str_time
        character(3) :: char3
        !
        write( char3, "(i3.3)" ) iter
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
    end subroutine outputFilesInversionDCG
    !
end module InversionDCG
!