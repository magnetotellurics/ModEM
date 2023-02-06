!
!> Abstract Base class to define a Inversion
!
module Inversion
    !
    use ModelParameter
    use DataGroupTxArray
    use Sensitivity
    !
    type, abstract :: Inversion_t
        !
        integer :: max_grad_iters, max_inv_iters, n_inv_iter
        !
        real( kind=prec ) :: tolerance_error, tolerance_rms, lambda
        !
        real( kind=prec ), allocatable, dimension(:) :: r_err
        !
        logical :: new_sigma
        !
        contains
            !
            procedure, public :: init => initializeInversion
            !
            procedure, public :: dealloc => deallocateInversion
            !
            procedure( interface_solve_inversion ), deferred, public :: solve
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
        subroutine interface_output_files_inversion( self, iter, all_predicted_data, res, dsigma, mHat )
            import :: Inversion_t, DataGroupTx_t, ModelParameter_t
            !
            class( Inversion_t ), intent( in ) :: self
            integer, intent( in ) :: iter
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
    subroutine initializeInversion( self )
        implicit none
        !
        class( Inversion_t ), intent( inout ) :: self
        !
        self%max_grad_iters = 20
        self%tolerance_error = 10E-4
        !
        if( has_inv_control_file ) then
            !
            if( allocated( inv_control_file%max_grad_iters ) ) &
                read( inv_control_file%max_grad_iters, * ) self%max_grad_iters
            !
            if( allocated( inv_control_file%tolerance_error ) ) &
                read( inv_control_file%tolerance_error, * ) self%tolerance_error
            !
        endif
        !
        write( *, "( A45 )" ) "Using inversion parameters:"
        !
        write( *, "( A45, I20 )" ) "max_grad_iters = ", self%max_grad_iters
        !
        write( *, "( A45, es20.2 )" ) "tolerance_error = ", self%tolerance_error
        !
        self%max_inv_iters = 0
        !
        self%n_inv_iter = 0
        !
        self%tolerance_rms = R_ZERO
        !
        self%lambda = R_ZERO
        !
        self%new_sigma = .TRUE.
        !
    end subroutine initializeInversion
    !
    !> No subroutine briefing
    subroutine deallocateInversion( self )
        implicit none
        !
        class( Inversion_t ), intent( inout ) :: self
        !
        deallocate( self%r_err )
        !
    end subroutine deallocateInversion
    !
end module Inversion
!