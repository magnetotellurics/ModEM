!
!> Abstract base class to define a ModelOperator
!
module ModelOperator
    !
    use MetricElements
    use ModelParameter
    !
    character(:), allocatable :: model_operator_type
    character( len=15 ), parameter :: MODELOP_MF = "MatrixFreeField"
    character( len=17 ), parameter :: MODELOP_SP = "SparseMatrixField"
    character( len=19 ), parameter :: MODELOP_SP2 = "SparseMatrixFieldV2"
    !
    type, abstract :: ModelOperator_t
        !
        class( MetricElements_t ), allocatable :: metric
        !
        class( Vector_t ), allocatable :: Adiag
        !
        integer :: mKey(8)
        !
        logical :: eqset, is_allocated
        !
        contains 
            !
            !> Abstract Interfaces
            !
            !> Setup
            procedure( interface_set_equations_model_operator ), deferred, public :: setEquations
            procedure( interface_set_cond_model_operator ), deferred, public :: setCond
            !
            procedure( interface_divcor_setup_model_operator ), deferred, public :: divCorSetUp
            !
            !> Operations
            procedure( interface_amult_model_operator ), deferred, public :: amult
            procedure( interface_multaib_model_operator ), deferred, public :: multAib
            !
            procedure( interface_div_model_operator ), deferred, public :: div
            procedure( interface_divc_model_operator ), deferred, public :: divC
            procedure( interface_divc_grad_model_operator ), deferred, public :: divCGrad
            !
            procedure( interface_grad_model_operator ), deferred, public :: grad
            !
            !> Miscellaneous
            procedure( interface_print_model_operator ), deferred, public :: print
            !
            !> Base procedures
            procedure, public :: baseInit => baseInit_ModelOperator
            procedure, public :: baseDealloc => baseDealloc_ModelOperator
            !
            procedure, public :: multCurlT => multCurlT_ModelOperator
            !
    end type ModelOperator_t
    !
    !> Public Global Generic ModelOperator object
    class( ModelOperator_t ), allocatable :: model_operator
    !
    public :: Maxwell
    !
    abstract interface
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_equations_model_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( inout ) :: self
            !
        end subroutine interface_set_equations_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_cond_model_operator( self, sigma, omega_in ) 
            import :: ModelOperator_t, ModelParameter_t, prec
            !
            class( ModelOperator_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in ) :: sigma
            real( kind=prec ), intent( in ) :: omega_in
            !
        end subroutine interface_set_cond_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_divcor_setup_model_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( inout ) :: self
            !
        end subroutine interface_divcor_setup_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_amult_model_operator( self, in_e, out_e, omega, adjoint )
            import :: ModelOperator_t, prec, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in ) :: in_e
            class( Vector_t ), intent( inout ) :: out_e
            real( kind=prec ), intent( in ) :: omega
            logical, intent( in ) :: adjoint
            !
        end subroutine interface_amult_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_multaib_model_operator( self, in_e, out_e )
            import :: ModelOperator_t, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in ) :: in_e
            class( Vector_t ), intent( inout ) :: out_e
            !
        end subroutine interface_multaib_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_div_model_operator( self, in_e, out_phi )
            import :: ModelOperator_t, Vector_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in ) :: in_e
            class( Scalar_t ), intent( inout ) :: out_phi
            !
        end subroutine interface_div_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_divc_model_operator( self, in_e, out_phi )
            import :: ModelOperator_t, Vector_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in) :: in_e
            class( Scalar_t ), intent( inout ) :: out_phi
            !
        end subroutine interface_divc_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_divc_grad_model_operator( self, in_phi, out_phi )
            import :: ModelOperator_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Scalar_t ), intent( in ) :: in_phi
            class( Scalar_t ), intent( inout ) :: out_phi
            !
        end subroutine interface_divc_grad_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_grad_model_operator( self, in_phi, out_e )
            import :: ModelOperator_t, Scalar_t, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Scalar_t ), intent( in ) :: in_phi
            class( Vector_t ), intent( inout ) :: out_e
            !
        end subroutine interface_grad_model_operator
        !!
        !> No interface subroutine briefing
        !
        subroutine interface_print_model_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            !
        end subroutine interface_print_model_operator
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine baseInit_ModelOperator( self )
        implicit none
        !
        class( ModelOperator_t ), intent( inout ) :: self
        !
        self%eqset = .FALSE.
        !
        self%is_allocated = .FALSE.
        !
        call date_and_time( values=self%mKey )
        !
    end subroutine baseInit_ModelOperator
    !
    !> No subroutine briefing
    subroutine baseDealloc_ModelOperator( self )
        implicit none
        !
        class( ModelOperator_t ), intent( inout ) :: self
        !
        if( allocated( self%metric ) ) deallocate( self%metric )
        !
        if( allocated( self%Adiag ) ) deallocate( self%Adiag )
        !
    end subroutine baseDealloc_ModelOperator
    !
    !> No subroutine briefing
    !
    subroutine multCurlT_ModelOperator( self, in_e, out_e )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: in_e
        class( Vector_t ), allocatable, intent( out ) :: out_e
        !
        integer :: ix, iy, iz
        complex( kind=prec ), allocatable, dimension(:,:,:) :: in_e_x, in_e_y, in_e_z
        complex( kind=prec ), allocatable, dimension(:,:,:) :: out_e_v
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "multCurlT_ModelOperator > in_e not allocated" )
        endif
        !
        call self%metric%createVector( complex_t, EDGE, out_e )
        call out_e%zeros
        !
        call in_e%div( self%metric%face_area )
        !
        in_e_x = in_e%getX()
        in_e_y = in_e%getY()
        in_e_z = in_e%getZ()
        !
        !> Ex
        out_e_v = out_e%getX()
        !
        do iy = 2, in_e%Ny
            do iz = 2, in_e%Nz
                out_e_v(:, iy, iz) = (in_e_z(:, iy, iz) - &
                in_e_z(:, iy - 1, iz)) - &
                (in_e_y(:, iy, iz) - in_e_y(:, iy, iz - 1))
            enddo
        enddo
        !
        call out_e%setX( out_e_v )
        !
        !> Ey
        out_e_v = out_e%getY()
        !
        do iz = 2, in_e%Nz
            do ix = 2, in_e%Nx
                out_e_v(ix, :, iz) = (in_e_x(ix, :, iz) - &
                in_e_x(ix, :, iz - 1)) - &
                (in_e_z(ix, :, iz) - in_e_z(ix - 1, :, iz))
            enddo
        enddo
        !
        call out_e%setY( out_e_v )
        !
        !> Ez
        out_e_v = out_e%getZ()
        !
        do ix = 2, in_e%Nx
            do iy = 2, in_e%Ny
                out_e_v(ix,iy,:) = (in_e_y(ix, iy, :) - &
                in_e_y(ix - 1, iy, :)) - &
                (in_e_x(ix, iy, :) - in_e_x(ix, iy - 1, :))
            enddo
        enddo
        !
        call out_e%setZ( out_e_v )
        !
        call out_e%mult( self%metric%edge_length )
        !
    end subroutine multCurlT_ModelOperator
    !
end module ModelOperator
