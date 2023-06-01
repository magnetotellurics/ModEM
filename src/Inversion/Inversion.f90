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
            procedure, public :: init => initializeInversion
            !
            procedure( interface_solve_inversion ), deferred, public :: solve
            !
            procedure( interface_output_files_inversion ), deferred, public :: outputFiles
            !
    end type Inversion_t
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
            class( ModelParameter_t ), intent( in ) :: dsigma, mHat
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
    ! !
    ! !> Compute the full penalty functional F
    ! !> Also output the predicted data and the EM solution
    ! !> that can be used for evaluating the gradient
    ! !> Assuming that the model norm is already scaled by Nmodel
    ! !
    ! subroutine printf( comment, lambda, alpha, f, mNorm, rms, logfile )
        ! implicit none
        ! !
        ! character(*), intent( in ) :: comment
        ! real( kind=prec ), intent( in ) :: lambda, alpha, f, mNorm, rms
        ! character(*), intent( in ), optional :: logfile
        ! integer :: io_unit, ios
        ! logical :: is_opened
        ! !
        ! if( present( logfile ) ) then
            ! !
            ! io_unit = ioLog
            ! !
            ! inquire( file = logfile, opened=is_opened )
            ! !
            ! if(.not. opened) then
                ! open (unit=ioLog,file=logfile,status='unknown',position='append',iostat=ios)
            ! endif
        ! else
            ! io_unit = 6
        ! endif
        ! !
        ! write(io_unit,'(a10)',advance='no') trim(comment)//':'
        ! write(io_unit,'(a3,es12.6)',advance='no') ' f=',f
        ! write(io_unit,'(a4,es12.6)',advance='no') ' m2=',mNorm
        ! write(io_unit,'(a5,f11.6)',advance='no') ' rms=',rms
        ! write(io_unit,'(a8,es12.6)',advance='no') ' lambda=',lambda
        ! write(io_unit,'(a7,es12.6)') ' alpha=',alpha
        ! !
        ! ! flush(io_unit): this has the effect of flushing the buffer
        ! if (present(logfile)) then
            ! close(io_unit)
            ! open( unit=ioLog,file=logfile,status='old',position='append',iostat=ios )
        ! end if
        ! !
    ! end subroutine printf
    ! !
end module Inversion
!