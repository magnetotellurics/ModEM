!
!> Abstract Base class to encapsulate the basic components needed to perform ModEM Inversion
!
module Inversion
    !
    use ModelParameter
    use DataGroupTxArray
    use Sensitivity
    !
    type, abstract :: Inversion_t
        !
        integer :: max_inv_iters, n_inv_iter
        !
        real( kind=prec ) :: tolerance_rms, lambda
        !
        contains
            !
            procedure, public :: init => initializeInversion
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
    !
    subroutine initializeInversion( self )
        implicit none
        !
        class( Inversion_t ), intent( inout ) :: self
        !
        self%max_inv_iters = 5
        self%tolerance_rms = 1.05
        self%lambda = 10.
        !
        if( has_inv_control_file ) then
            !
            if( allocated( inv_control_file%max_inv_iters ) ) &
                read( inv_control_file%max_inv_iters, * ) self%max_inv_iters
            !
            if( allocated( inv_control_file%tolerance_rms ) ) &
                read( inv_control_file%tolerance_rms, * ) self%tolerance_rms
            !
            if( allocated( inv_control_file%lambda ) ) &
                read( inv_control_file%lambda, * ) self%lambda
            !
        endif
        !
        write( *, "( A45, I20 )" ) "max_inv_iters = ", self%max_inv_iters
        !
        write( *, "( A45, es20.2 )" ) "tolerance_rms = ", self%tolerance_rms
        !
        write( *, "( A45, es20.2 )" ) "lambda = ", self%lambda
        !
        self%max_inv_iters = 0
        !
        self%n_inv_iter = 0
        !
        self%tolerance_rms = R_ZERO
        !
        self%lambda = R_ZERO
        !
    end subroutine initializeInversion
    !
end module Inversion
!