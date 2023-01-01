!
!> Class to provide a dynamic and polymorphic field_array of Field_t objects
!
module FieldArray
    !
    use Constants
    use Field
    !
    !> Allocatable Field element of the array, for Old Fortran polymorphism !!!
    type, public :: Fd_t
        !
        class( Field_t ), allocatable :: Fd
        !
    end type Fd_t
    !
    public :: getField, printFieldArray
    public :: updateFieldArray, deallocateFieldArray
    !
contains
    !
    !> Add a new Field_t and initialize it if necessary
    subroutine updateFieldArray( field_array, new_field )
        implicit none
        !
        type( Fd_t ), allocatable, dimension(:) :: field_array
        class( Field_t ), intent( in ) :: new_field
        !
        integer :: iFd, nFd
        type( Fd_t ), allocatable, dimension(:) :: temp_array
        type( Fd_t ), allocatable :: temp_field
        !
        if( .NOT. allocated( field_array ) ) then
            allocate( field_array(1) )
            allocate( Fd_t :: temp_field )
            temp_field%Fd = new_field
            field_array(1) = temp_field
            deallocate( temp_field )
        else
            !
            nFd = size( field_array )
            !
            allocate( temp_array( nFd + 1 ) )
            temp_array( 1 : nFd ) = field_array
            allocate( Fd_t :: temp_field )
            temp_field%Fd = new_field
            !
            temp_array( nFd + 1 ) = temp_field
            !
            if( allocated( field_array ) ) deallocate( field_array )
            allocate( field_array, source = temp_array )
            !
            deallocate( temp_field )
            deallocate( temp_array )
            !
        endif
        !
    end subroutine updateFieldArray
    !
    !> No function briefing
    function getField( field_array, iFd ) result( field )
        implicit none
        !
        type( Fd_t ), target, allocatable, dimension(:) :: field_array
        integer :: iFd
        !
        class( Field_t ), pointer :: field
        !
        field => field_array( iFd )%Fd
        !
    end function getField
    !
    !> No subroutine briefing
    subroutine deallocateFieldArray( field_array )
        implicit none
        !
        type( Fd_t ), allocatable, dimension(:) :: field_array
        !
        integer :: nfield, ifield
        !
        !write( *, * ) "deallocateFieldArray:", size( field_array )
        !
        nfield = size( field_array )
        !
        if( nfield == 1 ) then
            deallocate( field_array(1)%Fd )
        else
            do ifield = nfield, 1, -(1)
                deallocate( field_array(ifield)%Fd )
            enddo
        endif
        !
        deallocate( field_array )
        !
    end subroutine deallocateFieldArray
    !
    !> Prints the content of the field_array on screen
    subroutine printFieldArray( field_array )
        implicit none
        !
        type( Fd_t ), allocatable, dimension(:) :: field_array
        !
        integer :: ifield
        !
        write( *, * ) size( field_array ), " FieldArray_t:"
        !
        do ifield = 1, size( field_array )
            call field_array(ifield)%Fd%print()
        enddo
        !
    end subroutine printFieldArray
    !
end module FieldArray
