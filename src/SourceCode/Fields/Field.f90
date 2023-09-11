!
!> Abstract Base class to define a ModEM Field
!>     store_state: 1 - compound, 2 - singleton
!
module Field
    !
    use Utilities
    use Grid
    !
    !> Field Store States
    integer, parameter :: compound = 1
    integer, parameter :: singleton = 2
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
            !> Field Interfaces
            !
            !> Boundary operations
            procedure( interface_set_all_boundary_field ), deferred, public :: setAllBoundary
            procedure( interface_set_one_boundary_field ), deferred, public :: setOneBoundary
            !
            !> Dimensioning operations
            procedure( interface_length_field ), deferred, public :: length
            !
            !> Arithmetic/algebraic unary operations
            procedure( interface_zeros_field ), deferred, public :: zeros
            !
            procedure( interface_conjugate_field ), deferred, public :: conjugate
            !
            !> Arithmetic/algebraic binary operations
            procedure( interface_add_field ), deferred, public :: add
            !
            procedure( interface_lin_comb_field ), deferred, public :: linComb
            !
            procedure( interface_sub_value_field ), deferred, public :: subValue
            procedure( interface_sub_field_field ), deferred, public :: subField
            generic :: sub => subValue, subField
            !
            procedure( interface_field_mult_by_field ), deferred, public :: multByField
            procedure( interface_field_mult_by_complex ), deferred, public :: multByComplex
            procedure( interface_field_mult_by_real ), deferred, public :: multByReal
            generic :: mult => multByField, multByComplex, multByReal
            !
            procedure( interface_mult_add_field ), deferred, public :: multAdd
            !
            procedure( interface_dot_product_field ), deferred, public :: dotProd
            !
            procedure( interface_field_div_by_field ), deferred, public :: divByField
            procedure( interface_field_div_by_value ), deferred, public :: divByValue
            generic :: div => divByField, divByValue
            !
            !> Getters & Setters
            procedure( interface_get_array_field ), deferred, public :: getArray
            procedure( interface_set_array_field ), deferred, public :: setArray
            !
            !> Miscellaneous
            procedure( interface_copy_from_field ), deferred, public :: copyFrom
            generic :: assignment(=) => copyFrom
            !
            procedure( interface_deallocate_other_state_field ), deferred, public :: deallOtherState
            !
            !> Field procedures
            procedure, public :: baseInit => initialize_Field
            procedure, public :: baseDealloc => deallocate_Field
            !
            procedure, public :: switchStoreState => switchStoreState_Field
            !
            procedure, public :: isCompatible => isCompatible_Field
            !
            procedure, public :: setIndexArrays => setIndexArrays_Field
            !
            procedure, public :: indInterior => indInterior_Field
            procedure, public :: indBoundary => indBoundary_Field
            procedure, public :: indActive => indActive_Field
            !
            !> I/O operations
            procedure( interface_read_field ), deferred, public :: read
            procedure( interface_write_field ), deferred, public :: write
            procedure( interface_print_field ), deferred, public :: print
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
            class( Field_t ), intent( inout ) :: self
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
            character(*), intent( in ) :: bdry
            complex( kind=prec ), intent( in ) :: cvalue
            logical, intent( in ), optional :: int_only
        end subroutine interface_set_one_boundary_field
        !
        ! Dimensioning operations
        !
        !> No interface function briefing
        !
        function interface_length_field( self ) result( field_length )
            import :: Field_t
            class( Field_t ), intent( in ) :: self
            integer :: field_length
        end function interface_length_field
        !
        !> No interface function briefing
        !
        function interface_get_array_field( self ) result( array )
            import :: Field_t, prec
            class( Field_t ), intent( in ) :: self
            complex( kind=prec ), allocatable, dimension(:) :: array
        end function interface_get_array_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_array_field( self, array )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), dimension(:), intent( in ) :: array
        end subroutine interface_set_array_field
        !
        ! Arithmetic/algebraic operations
        !
        !> No interface subroutine briefing
        !
        subroutine interface_zeros_field( self )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
        end subroutine interface_zeros_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_add_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_add_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_sub_value_field( self, cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_sub_value_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_sub_field_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_sub_field_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_field_mult_by_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_field_mult_by_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_field_mult_by_complex( self, cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_field_mult_by_complex
        !
        !> No interface subroutine briefing
        !
        subroutine interface_field_mult_by_real( self, rvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            real( kind=prec ), intent( in ) :: rvalue
        end subroutine interface_field_mult_by_real
        !
        !> No interface subroutine briefing
        !
        subroutine interface_field_div_by_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_field_div_by_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_field_div_by_value( self, cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_field_div_by_value
        !
        !> No interface subroutine briefing
        !
        subroutine interface_conjugate_field( self )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
        end subroutine interface_conjugate_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_lin_comb_field( self, rhs, c1, c2 )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
            complex( kind=prec ), intent( in ) :: c1, c2
        end subroutine interface_lin_comb_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_mult_add_field( self, cvalue, rhs )
            import :: Field_t, prec
            class( Field_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_mult_add_field
        !
        !> No interface function briefing
        !
        function interface_dot_product_field( self, rhs ) result( cvalue )
            import :: Field_t, prec
            class( Field_t ), intent( in ) :: self
            class( Field_t ), intent( in ) :: rhs
            complex( kind=prec ) :: cvalue
        end function interface_dot_product_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_copy_from_field( self, rhs )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
        end subroutine interface_copy_from_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_deallocate_other_state_field( self )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
        end subroutine interface_deallocate_other_state_field
        !
        !> No interface function briefing
        !
        function interface_is_compatible_field( self, rhs ) result( is_compatible )
            import :: Field_t
            class( Field_t ), intent( in ) :: self, rhs
            logical :: is_compatible
        end function interface_is_compatible_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_print_field( self, io_unit, title, append )
            import :: Field_t
            class( Field_t ), intent( in ) :: self
            integer, intent( in ), optional :: io_unit
            character(*), intent( in ), optional :: title
            logical, intent( in ), optional :: append
        end subroutine interface_print_field
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_boundary_interior_field( self )
            import :: Field_t
            class( Field_t ), intent( inout ) :: self
        end subroutine interface_set_boundary_interior_field
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine initialize_Field( self )
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
        self%store_state = compound
        !
        self%is_allocated = .FALSE.
        !
    end subroutine initialize_Field
    !
    !> No subroutine briefing
    subroutine deallocate_Field( self )
        implicit none
        !
        class( Field_t ), intent( inout ) :: self
        !
    end subroutine deallocate_Field
    !
    !> No subroutine briefing
    !
    subroutine switchStoreState_Field( self, store_state )
        implicit none
        !
        class( Field_t ), intent( inout ) :: self
        integer, intent( in ) :: store_state
        !
        complex( kind=prec ), allocatable, dimension(:) :: field_array
        !
        if( self%store_state /= store_state ) then
            !
            field_array = self%getArray()
            !
            select case( self%store_state )
                !
                case( compound )
                    !
                    self%store_state = singleton
                    !
                case( singleton )
                    !
                    self%store_state = compound
                    !
                case default
                    call errStop( "switchStoreState_Field > store_state should be 'singleton' or 'compound'" )
                !
            end select
            !
            call self%setArray( field_array )
            !
        else
            !
            aux_counter = aux_counter + 1
            !
        endif
        !
    end subroutine switchStoreState_Field
    !
    !> No subroutine briefing
    !
    function isCompatible_Field( self, rhs ) result( is_compatible )
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
    end function isCompatible_Field
    !
    !> Defines the index arrays: ind_interior and ind_boundary.
    !>     Create copy with zeros and value boundaries with C_ONE.
    !>     Take two sizes and allocate the two arrays.
    !>     Fills the two arrays with their proper indices.
    !
    subroutine setIndexArrays_Field( self, n_full, ind_boundary, ind_interior, ind_active, xy_in )
        implicit none
        !
        class( Field_t ), intent( in ) :: self
        integer, intent( inout ) :: n_full
        integer, dimension(:), allocatable, intent( out ) :: ind_boundary, ind_interior
        integer, dimension(:), allocatable, intent( out ), optional :: ind_active
        logical, intent( in ), optional :: xy_in
        !
        integer :: i, j, k, int_size, bdry_size
        class( Field_t ), allocatable :: temp_field
        complex( kind=prec ), dimension(:), allocatable :: c_array
        !
        allocate( temp_field, source = self )
        !
        n_full = temp_field%length()
        !
        call temp_field%zeros
        !
        call temp_field%setAllBoundary( C_ONE )
        !
        c_array = temp_field%getArray()
        !
        deallocate( temp_field )
        !
        int_size = 0
        bdry_size = 0
        do i = 1, size( c_array )
            if( c_array(i) == C_ONE ) then
                bdry_size = bdry_size + 1
            else
                int_size = int_size + 1
            endif
        enddo
        !
        allocate( ind_boundary( bdry_size ) )
        !
        allocate( ind_interior( int_size ) )
        !
        j = 1
        k = 1
        do i = 1, size( c_array )
            if( c_array(i) == C_ONE ) then
                ind_boundary(j) = i
                j = j + 1
            else
                ind_interior(k) = i
                k = k + 1
            endif
        enddo
        !
        !write( *, * ) self%grid_type, size( ind_boundary ), size( ind_interior )
        !        
    end subroutine setIndexArrays_Field
    !
    ! No function briefing
    !
    function indBoundary_Field( self ) result( ind_boundary )
        implicit none
        !
        class( Field_t ), intent( in ) :: self
        !
        integer, dimension(:), allocatable :: ind_boundary
        !
        select case( self%grid_type )
            !
            case( EDGE )
                !
                ind_boundary = self%grid%EDGEb
                !
            case( FACE )
                !
                ind_boundary = self%grid%FACEb
                !
            case( NODE )
                !
                ind_boundary = self%grid%NODEb
                !
            case( CELL, CELL_EARTH )
                !
                call errStop( "CELL/CELL_EARTH indBoundary need to be implement" )
                !
            case default
                call errStop( "indBoundary > Invalid grid type ["//self%grid_type//"]" )
        end select 
        !
    end function indBoundary_Field
    !
    ! No function briefing
    !
    function indInterior_Field( self ) result( ind_interior )
        implicit none
        !
        class( Field_t ), intent( in ) :: self
        !
        integer, dimension(:), allocatable :: ind_interior
        !
        select case( self%grid_type )
            !
            case( EDGE )
                !
                ind_interior = self%grid%EDGEi
                !
            case( FACE )
                !
                ind_interior = self%grid%FACEi
                !
            case( NODE )
                !
                ind_interior = self%grid%NODEi
                !
            case( CELL, CELL_EARTH )
                !
                call errStop( "CELL/CELL_EARTH indInterior need to be implement" )
                !
            case default
                call errStop( "indInterior > Invalid grid type ["//self%grid_type//"]" )
        end select 
        !
    end function indInterior_Field
    !
    ! No function briefing
    !
    function indActive_Field( self ) result( ind_active )
        implicit none
        !
        class( Field_t ), intent( in ) :: self
        !
        integer, dimension(:), allocatable :: ind_active
        !
        select case( self%grid_type )
            !
            case( EDGE )
                !
                ind_active = self%grid%EDGEa
                !
            case( FACE )
                !
                ind_active = self%grid%FACEa
                !
            case( NODE )
                !
                ind_active = self%grid%NODEa
                !
            case( CELL, CELL_EARTH )
                !
                call errStop( "CELL/CELL_EARTH indActive need to be implement" )
                !
            case default
                call errStop( "indActive > Invalid grid type ["//self%grid_type//"]" )
        end select 
        !
    end function indActive_Field
    !
end module Field
