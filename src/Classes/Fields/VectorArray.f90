!
!> Class to provide a dynamic and polymorphic vector_array of Vector_t objects
!
module VectorArray
    !
    use Constants
    use Vector
    !
    !> Allocatable Vector element of the array, for Old Fortran polymorphism !!!
    type, public :: Vt_t
        !
        class( Vector_t ), allocatable :: Vt
        !
    end type Vt_t
    !
    public :: getVector, printVectorArray
    public :: updateVectorArray, deallocateVectorArray
    !
contains
    !
    !> Add a new Vector_t and initialize it if necessary
    subroutine updateVectorArray( vector_array, new_vector )
        implicit none
        !
        type( Vt_t ), allocatable, dimension(:) :: vector_array
        class( Vector_t ), intent( in ) :: new_vector
        !
        integer :: iVt, nVt
        type( Vt_t ), allocatable, dimension(:) :: temp_array
        type( Vt_t ), allocatable :: temp_vector
        !
        if( .NOT. allocated( vector_array ) ) then
            allocate( vector_array(1) )
            allocate( Vt_t :: temp_vector )
            allocate( temp_vector%Vt, source = new_vector )
            vector_array(1) = temp_vector
            deallocate( temp_vector )
        else
            !
            nVt = size( vector_array )
            !
            allocate( temp_array( nVt + 1 ) )
            temp_array( 1 : nVt ) = vector_array
            allocate( Vt_t :: temp_vector )
            allocate( temp_vector%Vt, source = new_vector )
            !
            temp_array( nVt + 1 ) = temp_vector
            !
            if( allocated( vector_array ) ) deallocate( vector_array )
            allocate( vector_array, source = temp_array )
            !
        endif
            !
            if( allocated( temp_vector ) ) deallocate( temp_vector )
            if( allocated( temp_array ) ) deallocate( temp_array )
            !
    end subroutine updateVectorArray
    !
    !> No function briefing
    function getVector( vector_array, iVt ) result( vector )
        implicit none
        !
        type( Vt_t ), target, allocatable, dimension(:) :: vector_array
        integer :: iVt
        !
        class( Vector_t ), pointer :: vector
        !
        vector => vector_array( iVt )%Vt
        !
    end function getVector
    !
    !> No subroutine briefing
    subroutine deallocateVectorArray( vector_array )
        implicit none
        !
        type( Vt_t ), allocatable, dimension(:) :: vector_array
        !
        integer :: nvector, ivector
        !
        !write( *, * ) "deallocateVectorArray:", size( vector_array )
        !
        nvector = size( vector_array )
        !
        if( nvector == 1 ) then
            deallocate( vector_array(1)%Vt )
        else
            do ivector = nvector, 1, -(1)
                deallocate( vector_array(ivector)%Vt )
            enddo
        endif
        !
        deallocate( vector_array )
        !
    end subroutine deallocateVectorArray
    !
    !> Prints the content of the vector_array on screen
    subroutine printVectorArray( vector_array )
        implicit none
        !
        type( Vt_t ), allocatable, dimension(:) :: vector_array
        !
        integer :: ivector
        !
        write( *, * ) size( vector_array ), " VectorArray_t:"
        !
        do ivector = 1, size( vector_array )
            call vector_array(ivector)%Vt%print()
        enddo
        !
    end subroutine printVectorArray
    !
end module VectorArray
