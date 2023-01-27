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
        integer :: max_cg_iter, max_iter, n_iter
        !
        real( kind=prec ) :: tol, rms_tol, lambda
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
            type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: all_data
            class( ModelParameter_t ), allocatable, intent( in ) :: sigma
            class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
            !
        end subroutine interface_solve_inversion
        !
        !> No interface subroutine briefing
        subroutine interface_output_files_inversion( self, iter, all_predicted_data, res, dsigma, mHat )
            import :: Inversion_t, DataGroupTx_t, ModelParameter_t
            !
            class( Inversion_t ), intent( inout ) :: self
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
        self%max_iter = 0
        !
        self%max_cg_iter = 20
        !
        self%n_iter = 0
        !
        self%tol = 10E-4
        !
        self%rms_tol = R_ZERO
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