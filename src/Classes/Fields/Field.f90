!
module Field
    !
    use Constants
    use Grid
    !
    ! STARORE STAGES: 0 - Full, 1 - Column
    !
    type, abstract :: Field_t
        !
        class( Grid_t ), pointer :: grid
        !
        character( len=4 ) :: grid_type
        !
        integer :: nx, ny, nz, store_state
        !
        logical :: is_allocated
        !
    contains
        !
        procedure( interface_read_field ), deferred, public  :: read
        procedure( interface_write_field ), deferred, public :: write
        !
        procedure( interface_set_all_boundary_field ), deferred, public :: setAllBoundary
        procedure( interface_set_one_boundary_field ), deferred, public :: setOneBoundary
        procedure( interface_set_all_interior_field ), deferred, public :: setAllInterior
        procedure( interface_int_bdry_indices_field ), deferred, public :: intBdryIndices
        !
        procedure( interface_length_field ), deferred, public :: length
        !
        procedure( interface_get_real_array_field ), deferred, public    :: getRealArray
        procedure( interface_get_complex_array_field ), deferred, public :: getComplexArray
        generic :: getArray => getRealArray, getComplexArray
        !
        procedure( interface_set_real_array_field ), deferred, public    :: setRealArray
        procedure( interface_set_complex_array_field ), deferred, public :: setComplexArray
        generic :: setArray => setRealArray, setComplexArray
        !
        procedure( interface_zeros_field ), deferred, public :: zeros
        procedure( interface_add_field ), deferred, public   :: add
        procedure( interface_sub_field ), deferred, public   :: sub
        !
        procedure( interface_copy_from_field ), deferred, public :: copyFrom
        generic :: assignment(=) => copyFrom
        !
        procedure( interface_print_field ), deferred, public :: print
        !
        procedure, public :: init => initializeField
        !
        procedure, public :: isCompatible => isCompatibleField
        !
    end type Field_t
    !
    abstract interface
        !
        ! I/O operation
        subroutine interface_read_field( self, funit, ftype )
            import :: Field_t
            class( Field_t ), intent( inout )    :: self
            integer, intent( in )                :: funit
            character(:), allocatable, intent( in ), optional :: ftype
        end subroutine interface_read_field
        !
        subroutine interface_write_field( self, funit, ftype )
            import :: Field_t
            class( Field_t ), intent( in ) :: self
            integer, intent( in )          :: funit
            character(:), allocatable, intent( in ), optional :: ftype
        end subroutine interface_write_field
        !
        ! Boundary operations
        subroutine interface_set_all_boundary_field( self, cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout )  :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_set_all_boundary_field
        !
        subroutine interface_set_one_boundary_field( self, bdry, cvalue, int_only )
            import :: Field_t, prec
            class( Field_t ), intent( inout )       :: self
            character(:), allocatable, intent( in ) :: bdry
            complex( kind=prec ), intent( in )      :: cvalue
            logical, intent( in ), optional         :: int_only
        end subroutine interface_set_one_boundary_field
        !
        subroutine interface_set_all_interior_field( self, cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout )  :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_set_all_interior_field
        !
        subroutine interface_int_bdry_indices_field( self, ind_i, ind_b )
            import :: Field_t, prec
            class( Field_t ), intent( in )      :: self
            integer, allocatable, intent( out ) :: ind_i(:), ind_b(:)
        end subroutine interface_int_bdry_indices_field
        !
        ! Dimensioning operations
        function interface_length_field( self ) result( field_length )
            import :: Field_t
            class( Field_t ), intent( in ) :: self
            integer                        :: field_length
        end function interface_length_field
        !
        subroutine interface_get_real_array_field( self, array )
            import :: Field_t, prec
            class( Field_t ), intent( in )                :: self
            real( kind=prec ), allocatable, intent( out ) :: array(:)
        end subroutine interface_get_real_array_field
        !
        subroutine interface_get_complex_array_field( self, array )
            import :: Field_t, prec
            class( Field_t ), intent( in )                   :: self
            complex( kind=prec ), allocatable, intent( out ) :: array(:)
        end subroutine interface_get_complex_array_field
        !
        subroutine interface_set_real_array_field( self, array )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            real( kind=prec ), intent( in )   :: array(:)
        end subroutine interface_set_real_array_field
        !
        subroutine interface_set_complex_array_field( self, array )
            import :: Field_t, prec
            class( Field_t ), intent( inout )  :: self
            complex( kind=prec ), intent( in ) :: array(:)
        end subroutine interface_set_complex_array_field
        !
        ! Arithmetic/algebraic operations
        subroutine interface_zeros_field( self )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
        end subroutine interface_zeros_field
        !
        subroutine interface_add_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in )    :: rhs
        end subroutine interface_add_field
        !
        subroutine interface_sub_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in )    :: rhs
        end subroutine interface_sub_field
        !
        ! Miscellaneous
        subroutine interface_copy_from_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in )    :: rhs
        end subroutine interface_copy_from_field
        !
        function interface_is_compatible_field( self, rhs ) result( is_compatible )
            import :: Field_t
            class( Field_t ), intent( in ) :: self, rhs
            logical :: is_compatible
        end function interface_is_compatible_field
        !
        subroutine interface_print_field( self, io_unit, title, append )
            import :: Field_t
            class( Field_t ), intent( in )                    :: self
            integer, intent( in ), optional                   :: io_unit
            character(:), allocatable, intent( in ), optional :: title
            logical, intent( in ), optional                   :: append
        end subroutine interface_print_field
        !
    end interface
    !
contains
    !
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
        self%store_state = 0
        !
        self%is_allocated = .FALSE.
        !
    end subroutine initializeField
    !
    function isCompatibleField( self, rhs ) result( is_compatible )
        implicit none
        !
        class( Field_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        logical :: is_compatible
        !
        is_compatible = .FALSE.
        !
        if( self%nx == rhs%nx .AND. self%ny == rhs%ny .AND. self%nz == rhs%nz .AND. &
            self%grid_type == rhs%grid_type ) then
            is_compatible = .TRUE.
        end if
        !
    end function isCompatibleField
    !
end module Field
