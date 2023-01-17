!
!> Abstract Base class to define a ModEM Field
!
module Field
    !
    use Constants
    use Grid
    !
    !> STORAGE STATES: 0 - Full, 1 - Column
    type, abstract :: Field_t
        !
        class( Grid_t ), pointer :: grid
        !
        character( len=4 ) :: grid_type
        !
        integer :: nx, ny, nz
        !
        logical :: is_allocated
        !
    contains
        !
        procedure( interface_read_field ), deferred, public :: read
        procedure( interface_write_field ), deferred, public :: write
        !
        procedure( interface_set_all_boundary_field ), deferred, public :: setAllBoundary
        procedure( interface_set_one_boundary_field ), deferred, public :: setOneBoundary
        procedure( interface_set_all_interior_field ), deferred, public :: setAllInterior
        procedure( interface_int_bdry_indices_field ), deferred, public :: intBdryIndices
        !
        procedure( interface_length_field ), deferred, public :: length
        !
        procedure( interface_get_array_field ), deferred, public :: getArray
        procedure( interface_set_array_field ), deferred, public :: setArray
        !
        procedure( interface_zeros_field ), deferred, public :: zeros
        procedure( interface_add_field ), deferred, public :: add
        !
        procedure( interface_sub_value_field ), deferred, public :: subValue
        procedure( interface_sub_field_field ), deferred, public :: subField
        generic :: sub => subValue, subField
        !
        procedure( interface_field_mult_by_field ), deferred, public :: multByField
        procedure( interface_field_mult_by_value ), deferred, public :: multByValue
        generic :: mult => multByField, multByValue
        !
        procedure( interface_field_div_by_field ), deferred, public :: divByField
        procedure( interface_field_div_by_value ), deferred, public :: divByValue
        generic :: div => divByField, divByValue
        !
        procedure( interface_conjugate_field ), deferred, public :: conjugate
        !
        procedure( interface_lin_comb_field ), deferred, public :: linComb
        !
        procedure( interface_mult_add_field ), deferred, public :: multAdd
        !
        procedure( interface_dot_product_field ), deferred, public :: dotProd
        generic :: operator(.dot.) => dotProd
        !
        procedure( interface_sum_edges_field ), deferred, public :: sumEdges
        !
        procedure( interface_get_real_field ), deferred, public :: getReal
        !
        procedure( interface_copy_from_field ), deferred, public :: copyFrom
        generic :: assignment(=) => copyFrom
        !
        procedure( interface_print_field ), deferred, public :: print
        !
        !> Field procedures
        procedure, public :: init => initializeField
        !
        procedure, public :: isCompatible => isCompatibleField
        !
        procedure, public :: boundary => boundaryField
        !
        procedure, public :: interior => interiorField
        !
    end type Field_t
    !
    abstract interface
        !
        ! I/O operations
        !
        !> No interface subroutine briefing
        subroutine interface_read_field( self, funit, ftype )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            integer, intent( in ) :: funit
            character(:), allocatable, intent( in ), optional :: ftype
        end subroutine interface_read_field
        !
        !> No interface subroutine briefing
        subroutine interface_write_field( self, funit, ftype )
            import :: Field_t
            class( Field_t ), intent( in ) :: self
            integer, intent( in ) :: funit
            character(:), allocatable, intent( in ), optional :: ftype
        end subroutine interface_write_field
        !
        ! Boundary operations
        !
        !> No interface subroutine briefing
        subroutine interface_set_all_boundary_field( self, cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_set_all_boundary_field
        !
        !> No interface subroutine briefing
        subroutine interface_set_one_boundary_field( self, bdry, cvalue, int_only )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            character(:), allocatable, intent( in ) :: bdry
            complex( kind=prec ), intent( in ) :: cvalue
            logical, intent( in ), optional :: int_only
        end subroutine interface_set_one_boundary_field
        !
        !> No interface subroutine briefing
        subroutine interface_set_all_interior_field( self, cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_set_all_interior_field
        !
        !> No interface subroutine briefing
        subroutine interface_int_bdry_indices_field( self, ind_i, ind_b )
            import :: Field_t, prec
            class( Field_t ), intent( in ) :: self
            integer, allocatable, intent( out ) :: ind_i(:), ind_b(:)
        end subroutine interface_int_bdry_indices_field
        !
        ! Dimensioning operations
        !
        !> No interface function briefing
        function interface_length_field( self ) result( field_length )
            import :: Field_t
            class( Field_t ), intent( in ) :: self
            integer :: field_length
        end function interface_length_field
        !
        !> No interface subroutine briefing
        subroutine interface_get_array_field( self, array )
            import :: Field_t, prec
            class( Field_t ), intent( in ) :: self
            complex( kind=prec ), allocatable, dimension(:), intent( out ) :: array
        end subroutine interface_get_array_field
        !
        !> No interface subroutine briefing
        subroutine interface_set_array_field( self, array )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), allocatable, dimension(:), intent( inout ) :: array
        end subroutine interface_set_array_field
        !
        ! Arithmetic/algebraic operations
        !
        !> No interface subroutine briefing
        subroutine interface_zeros_field( self )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
        end subroutine interface_zeros_field
        !
        !> No interface subroutine briefing
        subroutine interface_add_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_add_field
        !
        !> No interface subroutine briefing
        subroutine interface_sub_value_field( self, cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_sub_value_field
        !
        !> No interface subroutine briefing
        subroutine interface_sub_field_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_sub_field_field
        !
        !> No interface subroutine briefing
        subroutine interface_field_mult_by_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_field_mult_by_field
        !
        !> No interface subroutine briefing
        subroutine interface_field_mult_by_value( self, cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_field_mult_by_value
        !
        !> No interface subroutine briefing
        subroutine interface_field_div_by_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_field_div_by_field
        !
        !> No interface subroutine briefing
        subroutine interface_field_div_by_value( self, cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_field_div_by_value
        !
        !> No interface subroutine briefing
        subroutine interface_conjugate_field( self )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
        end subroutine interface_conjugate_field
        !
        !> No interface subroutine briefing
        subroutine interface_lin_comb_field( self, rhs, c1, c2 )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
            complex( kind=prec ), intent( in ) :: c1, c2
        end subroutine interface_lin_comb_field
        !
        !> No interface subroutine briefing
        subroutine interface_mult_add_field( self, cvalue, rhs )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_mult_add_field
        !
        !> No interface subroutine briefing
        subroutine interface_sum_edges_field( self, cell_obj, interior_only )
            import :: Field_t
            class( Field_t ), intent( in ) :: self
            class( Field_t ), allocatable, intent( inout ) :: cell_obj
            logical, optional, intent( in ) :: interior_only
        end subroutine interface_sum_edges_field
        !
        !> No interface function briefing
        function interface_dot_product_field( self, rhs ) result( cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( in ) :: self, rhs
            complex( kind=prec ) :: cvalue
        end function interface_dot_product_field
        !
        !> No interface subroutine briefing
        subroutine interface_get_real_field( self, r_field )
            import :: Field_t
            class( Field_t ), intent( in ) :: self
            class( Field_t ), allocatable, intent( out ) :: r_field
        end subroutine interface_get_real_field
        !
        !> No interface subroutine briefing
        subroutine interface_copy_from_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_copy_from_field
        !
        !> No interface function briefing
        function interface_is_compatible_field( self, rhs ) result( is_compatible )
            import :: Field_t
            class( Field_t ), intent( in ) :: self, rhs
            logical :: is_compatible
        end function interface_is_compatible_field
        !
        !> No interface subroutine briefing
        subroutine interface_print_field( self, io_unit, title, append )
            import :: Field_t
            class( Field_t ), intent( in ) :: self
            integer, intent( in ), optional :: io_unit
            character(*), intent( in ), optional :: title
            logical, intent( in ), optional :: append
        end subroutine interface_print_field
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine initializeField( self )
        implicit none
        !
        class( Field_t ), intent( inout ) :: self
        !
        self%grid => null()
        !
        self%grid_type = ""
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        !
        self%is_allocated = .FALSE.
        !
    end subroutine initializeField
    !
    !> No function briefing
    function isCompatibleField( self, rhs ) result( is_compatible )
        implicit none
        !
        class( Field_t ), intent( in ) :: self, rhs
        !
        logical :: is_compatible
        !
        is_compatible = .FALSE.
        !
        !write( *, * ) "SELF: ", self%nx, self%ny, self%nz, self%grid_type
        !write( *, * ) "RHS : ", rhs%nx, rhs%ny, rhs%nz, rhs%grid_type
        !
        if( self%nx == rhs%nx .AND. self%ny == rhs%ny .AND. self%nz == rhs%nz .AND. &
            self%grid_type == rhs%grid_type ) then
            is_compatible = .TRUE.
        endif
        !
    end function isCompatibleField
    !
    !> No subroutine briefing
    subroutine boundaryField( self, boundary )
        implicit none
        !
        class( Field_t ), intent( in ) :: self
        class( Field_t ), allocatable, intent( inout ) :: boundary
        !
        allocate( boundary, source = self )
        !
        call boundary%setAllInterior( C_ZERO )
       !
    end subroutine boundaryField
    !
    !> No subroutine briefing
    subroutine interiorField( self, interior )
        implicit none
        !
        class( Field_t ), intent( in ) :: self
        class( Field_t ), allocatable, intent( inout ) :: interior
        !
        allocate( interior, source = self )
        !
        call interior%setAllboundary( C_ZERO )
        !
    end subroutine interiorField
    !
end module Field
