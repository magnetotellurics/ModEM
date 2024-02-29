!
!> Derived class to define a InversionDCG
!
module InversionDCG
    !
    use Inversion
    !
    type, extends( Inversion_t ) :: InversionDCG_t
        !
        integer :: max_grad_iters
        !
        real( kind=prec ) :: error_tol
        !
        real( kind=prec ), allocatable, dimension(:) :: r_err
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
    !
    function InversionDCG_ctor() result( self )
        implicit none
        !
        type( InversionDCG_t ) :: self
        !
        !write( *, * ) "Constructor InversionDCG_t"
        !
        call self%init
        !
        self%max_grad_iters = 20
        self%error_tol = 10E-4
        !
        if( has_inv_control_file ) then
            !
            if( allocated( inv_control_file%max_grad_iters ) ) &
                read( inv_control_file%max_grad_iters, * ) self%max_grad_iters
            !
            if( allocated( inv_control_file%error_tol ) ) &
                read( inv_control_file%error_tol, * ) self%error_tol
            !
        endif
        !
        write( *, "( A45, I20 )" ) "max_grad_iters = ", self%max_grad_iters
        !
        write( *, "( A45, es20.2 )" ) "error_tol = ", self%error_tol
        !
        !> Free the memory used by the global control file, which is no longer useful
        if( allocated( inv_control_file ) ) deallocate( inv_control_file )
        !
        allocate( self%r_err( self%max_grad_iters ) )
        !
    end function InversionDCG_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    !
    subroutine InversionDCG_dtor( self )
        implicit none
        !
        type( InversionDCG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor InversionDCG_t"
        !
        deallocate( self%r_err )
        !
    end subroutine InversionDCG_dtor
    !
    !> No subroutine briefing
    !
    subroutine solveInversionDCG( self, all_data, sigma, dsigma )
        implicit none
        !
        class( InversionDCG_t ), intent( inout ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_data
        class( ModelParameter_t ), allocatable, intent( in ) :: sigma
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        real( kind=prec ) :: F, mNorm
        type( DataGroupTx_t ), allocatable, dimension(:) :: b, dx, all_predicted_data, res, JmHat
        class( ModelParameter_t ), allocatable :: mHat
        integer :: ios
        !
        !>
        call createOutputDirectory()
        !
        ! Verbose
        write( *, * ) "           > DCG output files in [", trim( outdir_name ), "]"
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
            call zerosData( b )
            !
            !> Write in DCG.log
            write( ioInvLog, "( a41, es12.5 )" ) "The initial damping parameter lambda is ", self%lambda
            !
            call self%Calc_FWD( all_data, dsigma, mHat, all_predicted_data, res, F, mNorm )
            !
            !> Write in DCG.log
            write( ioInvLog, "( a10, a3, es12.5, a4, es12.5, a5, f18.5, a8, es12.5 )" ) "START:", " f=", f, " m2=", mNorm, " rms=", self%rms, " lambda=", self%lambda
            !
            self%iter = 1
            !
            !> Print
            write( *, "( a38, f18.5)" ) "            Start_DCG : Residual rms=", self%rms
            !
            !>
            dcg_loop : do
                !
#ifdef MPI
                call masterJMult( dsigma, mHat, JmHat )
#else
                call serialJMult( dsigma, mHat, JmHat )
#endif
                !
                b = all_data
                !
                call linCombData( ONE, res, ONE, JmHat, b )
                !
                call normalizeData( b, 1 )
                !
                call CG_DS_standard( self, b, dx, dsigma, all_data )
                !
                call normalizeDataWith( 1, all_data, dx )
                !
#ifdef MPI
                call masterJMult_T( dsigma, dx, mHat )
#else
                call serialJMult_T( dsigma, dx, mHat )
#endif
                !
                call model_cov%multBy_Cm( mHat )
                !
                dsigma = sigma
                !
                call dsigma%linComb( ONE, ONE, mHat )
                !
                call self%Calc_FWD( all_data, dsigma, mHat, all_predicted_data, res, F, mNorm )
                !
                call self%outputFiles( all_predicted_data, res, dsigma, mHat )
                !
                !> Write / Print DCG.log
                write( *, "( a20, i5, a16, f18.5)" ) "            self%iter", self%iter, ": Residual rms=", self%rms
                !
                write( ioInvLog, "( a25, i5 )" ) "Completed DCG iteration ", self%iter
                write( ioInvLog, "( a10, a3, es12.5, a4, es12.5, a5, f18.5, a8, es12.5 )" ) "with:", " f=", f, " m2=", mNorm, " rms=", self%rms, " lambda=", self%lambda
                !
                !>
                if( self%rms .LT. self%rms_tol .OR. self%iter .GE. self%max_iters ) then
                    exit
                endif
                !
                self%iter = self%iter + 1
                !
            enddo dcg_loop
            !
            close( ioInvLog )
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
            write( *, * ) "Error opening [", trim( outdir_name )//"/DCG.log", "] in writeData!"
            stop
            !
        endif
        !
    end subroutine solveInversionDCG
    !
    !> No subroutine briefing
    !
    subroutine Calc_FWD( self, all_data, dsigma, mHat, all_predicted_data, res, F, mNorm )
        implicit none
        !
        class( InversionDCG_t ), intent( inout ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( in ) :: dsigma, mHat
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_predicted_data
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: res
        real( kind=prec ), intent( out ) :: F, mNorm
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
        call serialForwardModeling( dsigma, all_predicted_data )
        !
#endif
        !
        res = all_data
        !
        call linCombData( ONE, all_data, MinusONE, all_predicted_data, res )
        !
        Nres = res
        !
        call normalizeData( Nres, 2 )
        !
        SS = dotProdData( res, Nres )
        !
        Ndata = countValues( res )
        !
        mNorm = mHat%dotProd( mHat )
        !
        Nmodel = mHat%countModel()
        !
        F = SS / Ndata + ( self%lambda * mNorm / Nmodel )
        !
        self%rms = sqrt( SS / Ndata )
        !
    end subroutine Calc_FWD
    !
    !> No subroutine briefing
    !
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
        real( kind=prec ) :: r_norm_pre, r_norm, b_norm
        !
        r = b
        !
        p = r
        !
        Ap = all_data
        !
        b_norm = dotProdData( b, b )
        !
        call zerosData( x )
        !
        r_norm = dotProdData( r, r )
        !
        self%iter = 1
        !
        self%r_err( self%iter ) = r_norm / b_norm
        !
        !> Write / Print DCG.log
        write( ioInvLog, "(a18)" ) "Relative CG-error:"
        write( ioInvLog, "( a9, i5, a10, es12.5, a10, es12.5 )" ) "CG-Iter= ", self%iter, ", error = ", self%r_err( self%iter ), " Lambda= ", self%lambda
        !
        cg_loop : do while( self%r_err( self%iter ) .GT. self%error_tol .AND. self%iter .LT. self%max_grad_iters )
            ! 
            call self%MultA_DS( p, dsigma, all_data, Ap )
            !
            call setErrorBar( r, .FALSE. )
            call setErrorBar( p, .FALSE. )
            call setErrorBar( x, .FALSE. )
            call setErrorBar( Ap, .FALSE. )
            !
            ! Compute self%alpha: alpha= (r^T r) / (p^T Ap)
            self%alpha = r_norm / dotProdData( p, Ap )
            !
            ! Compute new x: x = x + alpha*p
            call scMultAddData( self%alpha, p, x )
            !
            ! Compute new r: r = r - alpha*Ap
            call scMultAddData( -self%alpha, Ap, r ) 
            !
            r_norm_pre = r_norm
            !
            r_norm = dotProdData( r, r )
            !
            self%beta = r_norm / r_norm_pre
            !
            ! Compute new p: p = r + beta*p    
            call linCombData( ONE, r, self%beta, p, p )
            !
            self%iter = self%iter + 1
            !
            self%r_err( self%iter ) = r_norm / b_norm 
            !
            !> Write / Print DCG.log
            write( ioInvLog, "( a9, i5, a10, es12.5, a10, es12.5 )" ) "CG-Iter= ", self%iter, ", error = ", self%r_err( self%iter ), " Lambda= ", self%lambda
            !
        enddo cg_loop
        !
        self%n_iter = self%iter
        !
    end subroutine CG_DS_standard
    !
    !> No subroutine briefing
    !
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
        call normalizeDataWith( 1, all_data, p_temp )
        !
#ifdef MPI
        call masterJMult_T( dsigma, p_temp, JTp )
#else
        call serialJMult_T( dsigma, p_temp, JTp )
#endif
        !
        call model_cov%multBy_Cm( JTp )
        !
        Ap = all_data
        !
#ifdef MPI
        call masterJMult( dsigma, JTp, Ap )
#else
        call serialJMult( dsigma, JTp, Ap )
#endif
        !
        deallocate( JTp )
        !
        call normalizeDataWith( 1, all_data, Ap )
        !
        call scMultData( self%lambda, p, lambdaP )
        !
        call setErrorBar( lambdaP, .FALSE. )
        !
        call linCombData( ONE, Ap, ONE, lambdaP, Ap )
        !
    end subroutine MultA_DS
    !
    !> No subroutine briefing
    !
    subroutine outputFilesInversionDCG( self, all_predicted_data, res, dsigma, mHat )
        implicit none
        !
        class( InversionDCG_t ), intent( in ) :: self
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_predicted_data, res
        class( ModelParameter_t ), intent( in ) :: dsigma, mHat
        !
        character(100) :: out_file_name
        character(3) :: char3
        !
        write( char3, "(i3.3)" ) self%iter
        !
        !> Write predicted data for this DCG iteration
        out_file_name = trim( outdir_name )//"/PredictedData_DCG_"//char3//".dat"
        !
        call writeData( all_predicted_data, trim( out_file_name ) )
        !
        !> Write residual data for this DCG iteration
        out_file_name = trim( outdir_name )//"/ResidualData_DCG_"//char3//".res"
        !
        call writeData( res, trim( out_file_name ) )
        !
        !> Write model for this DCG iteration
        out_file_name = trim( outdir_name )//"/SigmaModel_DCG_"//char3//".rho"
        !
        call dsigma%write( trim( out_file_name ) )
        !
        !> Write perturbation model for this DCG iteration
        out_file_name = trim( outdir_name )//"/PerturbationModel_DCG_"//char3//".prm"
        !
        call mHat%write( trim( out_file_name ) )
        !
    end subroutine outputFilesInversionDCG
    !
end module InversionDCG
!