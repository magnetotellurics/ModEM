!
!> Abstract Base class to encapsulate the basic components needed to perform ModEM Inversion
!
module Inversion
    !
    use ModelParameter
    use DataGroupTxArray
    use Sensitivity
    use TransmitterMT
    use TransmitterCSEM
    !
    type, abstract :: Inversion_t
        !
        integer :: iter, n_iter, max_iters
        !
        real( kind=prec ) :: rms_tol, lambda
        real( kind=prec ) :: alpha, beta, rms
        !
        contains
            !
            procedure( interface_solve_inversion ), deferred, public :: solve
            !
            procedure( interface_output_files_inversion ), deferred, public :: outputFiles
            !
            procedure, public :: baseInit => initializeInversion
            !
            procedure, public :: printInvPlot
            !
    end type Inversion_t
    !
    public :: printLog, printFuncPlot
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_solve_inversion( self, all_data, sigma, dsigma )
            import :: Inversion_t, DataGroupTx_t, ModelParameter_t
            !
            class( Inversion_t ), intent( inout ) :: self
            type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_data
            class( ModelParameter_t ), allocatable, intent( in ) :: sigma
            class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
            !
        end subroutine interface_solve_inversion
        !
        !> No interface subroutine briefing
        subroutine interface_output_files_inversion( self, all_predicted_data, res, dsigma, mHat )
            import :: Inversion_t, DataGroupTx_t, ModelParameter_t
            !
            class( Inversion_t ), intent( in ) :: self
            type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_predicted_data, res
            class( ModelParameter_t ), intent( inout ) :: dsigma, mHat
            !
        end subroutine interface_output_files_inversion
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine initializeInversion( self )
        implicit none
        !
        class( Inversion_t ), intent( inout ) :: self
        !
        self%max_iters = 5
        self%rms_tol = 1.05
        self%lambda = 10.
        !
        if( has_inv_control_file ) then
            !
            if( allocated( inv_control_file%max_inv_iters ) ) &
                read( inv_control_file%max_inv_iters, * ) self%max_iters
            !
            if( allocated( inv_control_file%rms_tol ) ) &
                read( inv_control_file%rms_tol, * ) self%rms_tol
            !
            if( allocated( inv_control_file%lambda ) ) &
                read( inv_control_file%lambda, * ) self%lambda
            !
        endif
        !
        write( *, "( A45 )" ) "Using Inversion parameters:"
        !
        write( *, "( A45, I20 )" ) "max_inv_iters = ", self%max_iters
        !
        write( *, "( A45, es20.2 )" ) "rms_tol = ", self%rms_tol
        !
        write( *, "( A45, es20.2 )" ) "lambda = ", self%lambda
        !
        self%n_iter = 0
        !
    end subroutine initializeInversion
    !
    !> Compute the full penalty functional F
    !> Also output the predicted data and the EM solution
    !> that can be used for evaluating the gradient
    !> Assuming that the model norm is already scaled by Nmodel
    !
    subroutine printInvPlot( self, g_norm )
        implicit none
        !
        class( Inversion_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: g_norm
        !
        integer :: ios
        !
        open( unit = ioInvPlot, file = trim( outdir_name )//"/"//trim( outdir_name )//".inv_plot", &
        status="unknown", position="append", iostat=ios )
        !
        if( ios == 0 ) then
            !
            write( ioInvPlot, * ) self%iter, ", ", self%alpha, ", ", &
                                  self%beta, ", ", g_norm, ", ",  self%rms
            !
            close( ioInvPlot )
        else
            call errStop( "printFuncPlot > cant open ioInvPlot" )
        endif
        !
    end subroutine printInvPlot
    !
    subroutine printLog( comment, lambda, alpha, f, mNorm, rms, on_screen )
        implicit none
        !
        character(*), intent( in ) :: comment
        real( kind=prec ), intent( in ) :: lambda, alpha, f, mNorm, rms
        logical, intent( in ) :: on_screen
        !
        integer :: io_unit, ios
        logical :: is_open
        !
        if( on_screen ) then
            io_unit = 0
        else
            !
            io_unit = ioInvLog
            !
            inquire( file = trim( outdir_name )//"/NLCG.log", opened = is_open )
            !
            if( .NOT. is_open ) then
                !
                open( unit = io_unit, file = trim( outdir_name )//"/NLCG.log", status="unknown", position="append", iostat=ios )
                !
            endif
            !
        endif
        !
        write( io_unit, "(a10)", advance="no" ) trim(comment)//":"
        write( io_unit, "(a3,es15.3)", advance="no" ) " f=", f
        write( io_unit, "(a4,es15.3)", advance="no" ) " m2=", mNorm
        write( io_unit, "(a5,f15.3)", advance="no" ) " rms=", rms
        write( io_unit, "(a8,es15.3)", advance="no" ) " lambda=", lambda
        write( io_unit, "(a7,es15.3)") " alpha=", alpha
        !
        if( .NOT. on_screen ) then
            !
            close( io_unit )
            !
            open( unit = io_unit, file = trim( outdir_name )//"/NLCG.log", status="unknown", position="append", iostat=ios )
            !
        endif
        !
    end subroutine printLog
    !
    subroutine printFuncPlot( SS, Ndata, m_norm, n_model, F, rms )
        implicit none
        !
        real( kind=prec ), intent( in ) :: SS
        integer, intent( in ) :: Ndata
        real( kind=prec ), intent( in ) :: m_norm
        integer, intent( in ) :: n_model
        real( kind=prec ), intent( in ) :: F, rms
        !
        integer :: ios
        !
        open( unit = ioFuncPlot, &
        file = trim( outdir_name )//"/"//trim( outdir_name )//".func_plot", &
        status="unknown", position="append", iostat=ios )
        !
        if( ios == 0 ) then
            !
            write( ioFuncPlot, * ) SS, ", ", Ndata, ", ", m_norm, ", ", n_model, ", ", F, ", ", rms
            !
            close( ioFuncPlot )
        else
            call errStop( "printFuncPlot > cant open ioFuncPlot" )
        endif
        !
    end subroutine printFuncPlot
    !
end module Inversion
!