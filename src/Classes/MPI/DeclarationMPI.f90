!*************
!
! MPI declarations
!
! Last modified at 03/2022 by Paulo Werdt
!
!*************
!
module DeclarationMPI
    !
    use, intrinsic :: iso_c_binding, only: c_ptr, c_sizeof, c_f_pointer
    !
    use mpi_f08
    !
    use Grid3D_SG
    use rVector3D_SG
    use rScalar3D_SG
    use cVector3D_SG
    use cScalar3D_SG
    use MetricElements_CSG
    use ModelOperator_MF
    use ModelParameterCell_SG
    !
    type( MPI_Comm ) :: main_comm, child_comm
    !
    integer :: mpi_rank, mpi_size, node_rank, node_size, &
               nodestringlen, ierr
    !
    type( MPI_Win )  :: shared_window
    !
    integer( MPI_ADDRESS_KIND ) :: shared_window_size
    !
    character, dimension(:), pointer :: shared_buffer
    character, dimension(:), pointer :: job_buffer
    !
    integer                          :: disp_unit = 1
    integer                          :: shared_buffer_size = 1
    integer                          :: job_size_bytes = 1
    !
    type( c_ptr ) :: shared_c_ptr
    !
    integer    :: field_derived_type
    integer, parameter :: real_3d_sg = 1
    integer, parameter :: complex_3d_sg = 2
    !
    integer    :: grid_derived_type
    integer, parameter :: grid_3d_sg = 1
    !
    integer    :: metric_derived_type
    integer, parameter :: metric_csg = 1
    !
    integer    :: model_operator_derived_type
    integer, parameter :: model_operator_mf = 1
    !
    integer    :: model_parameter_derived_type
    integer, parameter :: model_parameter_cell_sg = 1
    !
    character*( MPI_MAX_PROCESSOR_NAME ) :: node_name
    !
    ! PROGRAM GLOBAL VARIABLES
    integer :: tag = 666, master_id = 0, max_procs = 4
    !
    character*15    :: job_master = "MASTER_JOB", job_finish = "STOP_JOBS", job_done = "FINISH_JOB"
    !
    ! STRUCT job_info
    type :: struct_job_info
        SEQUENCE
        character*50 :: name
        character*20 :: transmitter_derived_type
        integer      :: tx_index
        integer      :: worker_rank
    end type struct_job_info
    !
    type( struct_job_info ), save :: job_info
    !
    contains
    !
    ! ALLOCATE shared_buffer
    subroutine allocateSharedBuffer( grid, model_operator, model_parameter )
        implicit none
        !
        class( Grid_t ), intent( in )           :: grid
        class( ModelOperator_t ), intent( in )  :: model_operator
		class( ModelParameter_t ), intent( in ) :: model_parameter
        !
        !
        shared_buffer_size = shared_buffer_size + allocateGridBuffer( grid )
        !
        write( *, * ) "$$$$$$$$$ allocateSharedBuffer GRID OK: ", shared_buffer_size
        !
        shared_buffer_size = shared_buffer_size + allocateModelOperatorBuffer( model_operator )
        !
        write( *, * ) "$$$$$$$$$ allocateSharedBuffer MODEL OPERATOR OK: ", shared_buffer_size
        !
        shared_buffer_size = shared_buffer_size + allocateModelParameterBuffer( model_parameter )
        !
        write( *, * ) "$$$$$$$$$ allocateSharedBuffer MODEL PARAMETER OK: ", shared_buffer_size
        !
        if( .NOT. associated( shared_buffer ) ) then
            allocate( shared_buffer( shared_buffer_size ) )
        endif
        !
    end subroutine allocateSharedBuffer
    !
    ! PACK job_info STRUCT TO Grid Object
    subroutine packSharedBuffer( grid, model_operator, model_parameter )
        implicit none
        !
        class( Grid_t ), intent( in )          :: grid
        class( ModelOperator_t ), intent( in )   :: model_operator
		class( ModelParameter_t ), intent( in ) :: model_parameter
        !
        !
        integer :: index = 1
        !
        call packGridBuffer( grid, index )
        !
        write( *, * ) "$$$$$$$$$ packSharedBuffer GRID OK: ", index
        !
        call packModelOperatorBuffer( grid, model_operator, index )
        !
        write( *, * ) "$$$$$$$$$ packSharedBuffer MODEL OPERATOR OK: ", index
        !
        call packModelParameterBuffer( model_parameter, index )
        !
        write( *, * ) "$$$$$$$$$ packSharedBuffer MODEL PARAMETER OK: ", index
        !
    end subroutine packSharedBuffer
    !
    ! PACK job_info STRUCT TO Grid Object
    subroutine unpackSharedBuffer( buffer_size, grid, model_operator, model_parameter )
        implicit none
        !
        class( Grid_t ), allocatable, intent( inout )           :: grid
        class( ModelOperator_t ), allocatable, intent( inout )  :: model_operator
		class( ModelParameter_t ), allocatable, intent( inout ) :: model_parameter
        !
        !
        integer, intent( in ) :: buffer_size
        !
        integer :: index = 1
        !
        shared_buffer_size = buffer_size
        !
        grid = unpackGridBuffer( index )
        !
        write( *, * ) "$$$$$$$$$ unpackSharedBuffer GRID OK: ", index
        !
        model_operator = unpackModelOperatorBuffer( grid, index )
        !
        write( *, * ) "$$$$$$$$$ unpackSharedBuffer MODEL OPERATOR OK: ", index
        !
        model_parameter = unpackModelParameterBuffer( grid, model_operator%metric, index )
        !
        write( *, * ) "$$$$$$$$$ unpackSharedBuffer MODEL PARAMETER OK: ", index
        !
    end subroutine unpackSharedBuffer
    !
    !
    function allocateRScalarBuffer( scalar ) result( scalar_size_bytes )
        implicit none
        !
        class( rScalar_t ), intent( in ):: scalar
        !
        integer :: i, nbytes(5), scalar_size_bytes
        !
        scalar_size_bytes = 0
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, child_comm, nbytes(1), ierr )
        !
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, child_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 7, MPI_INTEGER, child_comm, nbytes(4), ierr )
        !
        select type( scalar )
           !
           class is( rScalar3D_SG_t )
              !
              call MPI_PACK_SIZE( scalar%Nxyz, MPI_REAL, child_comm, nbytes(5), ierr )
              !
           class default
              stop "allocateRScalarBuffer: Unclassified scalar"
           !
        end select
        !
        do i = 1, size( nbytes )
           scalar_size_bytes = scalar_size_bytes + nbytes(i)
        end do
        !
    end function allocateRScalarBuffer
    !
    !
    subroutine packRScalarBuffer( scalar, index )
        implicit none
        !
        class( rScalar_t ), intent( in ) :: scalar
        integer, intent( inout )        :: index
        !
        real( kind=prec ), allocatable    :: r_array(:)
        !
        select type( scalar )
           !
           class is( rScalar3D_SG_t )
              !
              call MPI_PACK( real_3d_sg, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call MPI_PACK( scalar%isAllocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%gridType, 1, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%NdV(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%Nxyz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call scalar%getArray( r_array )
              call MPI_PACK( r_array(1), scalar%Nxyz, MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
           class default
              stop "packRScalarBuffer: Unclassified scalar"
           !
        end select
        !
    end subroutine packRScalarBuffer
    !
    function unpackRScalarBuffer( grid, index ) result( scalar )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        integer, intent( inout )              :: index
        !
        class( rScalar_t ), allocatable :: scalar
        !
        real( kind=prec ), allocatable       :: r_array(:)
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, field_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        !
        select case( field_derived_type )
            !
            case( real_3d_sg )
                !
                allocate( rScalar3D_SG_t :: scalar )
                !
                select type( scalar )
                    !
                    class is ( rScalar3D_SG_t )
                        !
                        select type( grid )
                           !
                           class is( Grid3D_SG_t )
                              !
                              scalar%grid => grid
                              !
                           class default
                              stop "unpackRScalarBuffer: Unclassified grid"
                           !
                        end select
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%isAllocated, 1, MPI_LOGICAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%gridType, 1, MPI_CHARACTER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%nx, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%ny, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%nz, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%NdV(1), 3, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%Nxyz, 1, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( r_array( scalar%Nxyz ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), scalar%Nxyz, MPI_REAL, child_comm, ierr )
                        call scalar%setArray( r_array )
                        !
                    class default
                        stop "unpackRScalarBuffer: Unclassified scalar"
                    !
                end select
                !
            case default
                stop "unpackRScalarBuffer: Unclassified scalar"
            !
        end select
        !
    end function unpackRScalarBuffer
    !
    !
    function allocateCScalarBuffer( scalar ) result( scalar_size_bytes )
        implicit none
        !
        class( cScalar_t ), intent( in ):: scalar
        !
        integer :: i, nbytes(5), scalar_size_bytes
        !
        scalar_size_bytes = 0
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, child_comm, nbytes(1), ierr )
        !
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, child_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 7, MPI_INTEGER, child_comm, nbytes(4), ierr )
        !
        select type( scalar )
           !
           class is( cScalar3D_SG_t )
              !
              call MPI_PACK_SIZE( scalar%Nxyz, MPI_COMPLEX, child_comm, nbytes(5), ierr )
              !
           class default
              stop "allocateCScalarBuffer: Unclassified scalar"
           !
        end select
        !
        do i = 1, size( nbytes )
           scalar_size_bytes = scalar_size_bytes + nbytes(i)
        end do
        !
    end function allocateCScalarBuffer
    !
    !
    subroutine packCScalarBuffer( scalar, index )
        implicit none
        !
        class( cScalar_t ), intent( in ) :: scalar
        integer, intent( inout )        :: index
        !
        complex( kind=prec ), allocatable :: c_array(:)
        !
        select type( scalar )
           !
           class is( cScalar3D_SG_t )
              !
              call MPI_PACK( complex_3d_sg, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call MPI_PACK( scalar%isAllocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%gridType, 1, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%NdV(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%Nxyz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call scalar%getArray( c_array )
              call MPI_PACK( c_array(1), scalar%Nxyz, MPI_COMPLEX, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
           class default
              stop "packCScalarBuffer: Unclassified scalar"
           !
        end select
        !
    end subroutine packCScalarBuffer
    !
    function unpackCScalarBuffer( grid, index ) result( scalar )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        integer, intent( inout )              :: index
        !
        class( cScalar_t ), allocatable :: scalar
        !
        complex( kind=prec ), allocatable    :: c_array(:)
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, field_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        !
        select case( field_derived_type )
            !
            case( complex_3d_sg )
                !
                allocate( cScalar3D_SG_t :: scalar )
                !
                select type( scalar )
                    !
                    class is ( cScalar3D_SG_t )
                        !
                        select type( grid )
                           !
                           class is( Grid3D_SG_t )
                              !
                              scalar%grid => grid
                              !
                           class default
                              stop "unpackCScalarBuffer: Unclassified grid"
                           !
                        end select
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%isAllocated, 1, MPI_LOGICAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%gridType, 1, MPI_CHARACTER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%nx, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%ny, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%nz, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%NdV(1), 3, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%Nxyz, 1, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( c_array( scalar%Nxyz ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, c_array(1), scalar%Nxyz, MPI_COMPLEX, child_comm, ierr )
                        call scalar%setArray( c_array )
                        !
                    class default
                        stop "unpackCScalarBuffer: Unclassified scalar"
                    !
                end select
                !
            case default
                stop "unpackCScalarBuffer: Unclassified scalar"
            !
        end select
        !
    end function unpackCScalarBuffer
    !
    function allocateRVectorBuffer( vector ) result( vector_size_bytes )
        implicit none
        !
        class( rVector_t ), intent( in ):: vector
        !
        integer :: i, nbytes(5), vector_size_bytes
        !
        vector_size_bytes = 0
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, child_comm, nbytes(1), ierr )
        !
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, child_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 15, MPI_INTEGER, child_comm, nbytes(4), ierr )
        !
        select type( vector )
           !
           class is( rVector3D_SG_t )
              !
              call MPI_PACK_SIZE( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_REAL, child_comm, nbytes(5), ierr )
              !
           class default
              stop "allocateRVectorBuffer: Unclassified vector"
           !
        end select
        !
        do i = 1, size( nbytes )
           vector_size_bytes = vector_size_bytes + nbytes(i)
        end do
        !
    end function allocateRVectorBuffer
    !
    !
    subroutine packRVectorBuffer( vector, index )
        implicit none
        !
        class( rVector_t ), intent( in ) :: vector
        integer, intent( inout )        :: index
        !
        real( kind=prec ), allocatable    :: r_array(:)
        !
        select type( vector )
           !
           class is( rVector3D_SG_t )
              !
              call MPI_PACK( real_3d_sg, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call MPI_PACK( vector%isAllocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%gridType, 1, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdX(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdY(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdZ(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%Nxyz(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call vector%getArray( r_array )
              call MPI_PACK( r_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
           class default
              stop "packRVectorBuffer: Unclassified vector"
           !
        end select
        !
    end subroutine packRVectorBuffer
    !
    function unpackRVectorBuffer( grid, index ) result( vector )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        integer, intent( inout )              :: index
        !
        class( rVector_t ), allocatable :: vector
        !
        real( kind=prec ), allocatable  :: r_array(:)
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, field_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        !
        select case( field_derived_type )
            !
            case( real_3d_sg )
                !
                allocate( rVector3D_SG_t :: vector )
                !
                select type( vector )
                    !
                    class is ( rVector3D_SG_t )
                        !
                        select type( grid )
                           !
                           class is( Grid3D_SG_t )
                              !
                              vector%grid => grid
                              !
                           class default
                              stop "unpackRVectorBuffer: Unclassified grid"
                           !
                        end select
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%isAllocated, 1, MPI_LOGICAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%gridType, 1, MPI_CHARACTER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%nx, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%ny, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%nz, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdX(1), 3, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdY(1), 3, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdZ(1), 3, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%Nxyz(1), 3, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( r_array( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3) ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_REAL, child_comm, ierr )
                        call vector%setArray( r_array )
                        !
                    class default
                        stop "unpackRVectorBuffer: Unclassified vector"
                    !
                end select
                !
            case default
                stop "unpackRVectorBuffer: Unclassified vector"
            !
        end select
        !
    end function unpackRVectorBuffer
    !
    function allocateCVectorBuffer( vector ) result( vector_size_bytes )
        implicit none
        !
        class( cVector_t ), intent( in ):: vector
        !
        integer :: i, nbytes(5), vector_size_bytes
        !
        vector_size_bytes = 0
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, child_comm, nbytes(1), ierr )
        !
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, child_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 15, MPI_INTEGER, child_comm, nbytes(4), ierr )
        !
        select type( vector )
           !
           class is( cVector3D_SG_t )
              !
              call MPI_PACK_SIZE( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_COMPLEX, child_comm, nbytes(5), ierr )
              !
           class default
              stop "allocateCVectorBuffer: Unclassified vector"
           !
        end select
        !
        do i = 1, size( nbytes )
           vector_size_bytes = vector_size_bytes + nbytes(i)
        end do
        !
    end function allocateCVectorBuffer
    !
    !
    subroutine packCVectorBuffer( vector, index )
        implicit none
        !
        class( cVector_t ), intent( in ) :: vector
        integer, intent( inout )         :: index
        !
        complex( kind=prec ), allocatable :: c_array(:)
        !
        select type( vector )
           !
           class is( cVector3D_SG_t )
              !
              call MPI_PACK( complex_3d_sg, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call MPI_PACK( vector%isAllocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%gridType, 1, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdX(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdY(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdZ(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%Nxyz(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call vector%getArray( c_array )
              call MPI_PACK( c_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_COMPLEX, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
           class default
              stop "packCVectorBuffer: Unclassified vector"
           !
        end select
        !
    end subroutine packCVectorBuffer
    !
    function unpackCVectorBuffer( grid, index ) result( vector )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        integer, intent( inout )              :: index
        !
        class( cVector_t ), allocatable :: vector
        !
        complex( kind=prec ), allocatable :: c_array(:)
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, field_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        !
        select case( field_derived_type )
            !
            case( complex_3d_sg )
                !
                allocate( cVector3D_SG_t :: vector )
                !
                select type( vector )
                    !
                    class is ( cVector3D_SG_t )
                        !
                        select type( grid )
                           !
                           class is( Grid3D_SG_t )
                              !
                              vector%grid => grid
                              !
                           class default
                              stop "unpackCVectorBuffer: Unclassified grid"
                           !
                        end select
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%isAllocated, 1, MPI_LOGICAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%gridType, 1, MPI_CHARACTER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%nx, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%ny, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%nz, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdX(1), 3, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdY(1), 3, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdZ(1), 3, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%Nxyz(1), 3, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( c_array( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3) ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, c_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_COMPLEX, child_comm, ierr )
                        call vector%setArray( c_array )
                        !
                    class default
                        stop "unpackCVectorBuffer: Unclassified vector"
                    !
                end select
                !
            case default
                stop "unpackCVectorBuffer: Unclassified vector"
            !
        end select
        !
    end function unpackCVectorBuffer
    !
    ! ALLOCATE shared_buffer
    function allocateGridBuffer( grid ) result( grid_size_bytes )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        integer :: i, nbytes(24), grid_size_bytes
        !
        grid_size_bytes = 0
        !
        ! SIZES FOR THE FUTURE
        call MPI_PACK_SIZE( 19, MPI_INTEGER, child_comm, nbytes(1), ierr )
        !
        call MPI_PACK_SIZE( 80, MPI_CHARACTER, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 4, MPI_REAL, child_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 5, MPI_INTEGER, child_comm, nbytes(4), ierr )
        call MPI_PACK_SIZE( size( grid%dx ), MPI_REAL, child_comm, nbytes(5), ierr )
        call MPI_PACK_SIZE( size( grid%dy ), MPI_REAL, child_comm, nbytes(6), ierr )
        call MPI_PACK_SIZE( size( grid%dz ), MPI_REAL, child_comm, nbytes(7), ierr )
        call MPI_PACK_SIZE( size( grid%dxInv ), MPI_REAL, child_comm, nbytes(8), ierr )
        call MPI_PACK_SIZE( size( grid%dyInv ), MPI_REAL, child_comm, nbytes(9), ierr )
        call MPI_PACK_SIZE( size( grid%dzInv ), MPI_REAL, child_comm, nbytes(10), ierr )
        call MPI_PACK_SIZE( size( grid%delX ), MPI_REAL, child_comm, nbytes(11), ierr )
        call MPI_PACK_SIZE( size( grid%delY ), MPI_REAL, child_comm, nbytes(12), ierr )
        call MPI_PACK_SIZE( size( grid%delZ ), MPI_REAL, child_comm, nbytes(13), ierr )
        call MPI_PACK_SIZE( size( grid%delXInv ), MPI_REAL, child_comm, nbytes(14), ierr )
        call MPI_PACK_SIZE( size( grid%delYInv ), MPI_REAL, child_comm, nbytes(15), ierr )
        call MPI_PACK_SIZE( size( grid%delZInv ), MPI_REAL, child_comm, nbytes(16), ierr )
        call MPI_PACK_SIZE( size( grid%xEdge ), MPI_REAL, child_comm, nbytes(17), ierr )
        call MPI_PACK_SIZE( size( grid%yEdge ), MPI_REAL, child_comm, nbytes(18), ierr )
        call MPI_PACK_SIZE( size( grid%zEdge ), MPI_REAL, child_comm, nbytes(19), ierr )
        call MPI_PACK_SIZE( size( grid%xCenter ), MPI_REAL, child_comm, nbytes(20), ierr )
        call MPI_PACK_SIZE( size( grid%yCenter ), MPI_REAL, child_comm, nbytes(21), ierr )
        call MPI_PACK_SIZE( size( grid%zCenter ), MPI_REAL, child_comm, nbytes(22), ierr )
        call MPI_PACK_SIZE( 1, MPI_REAL, child_comm, nbytes(23), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(24), ierr )
        !
        do i = 1, size( nbytes )
           grid_size_bytes = grid_size_bytes + nbytes(i)
        end do
        !
    end function allocateGridBuffer
    !
    ! PACK job_info STRUCT TO Grid Object
    subroutine packGridBuffer( grid, index )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        integer, intent( inout )      :: index
        !
        select type( grid )
           !
           class is( Grid3D_SG_t )
                !
                call MPI_PACK( grid_3d_sg, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                ! SIZES FOR THE FUTURE
                call MPI_PACK( size( grid%dx ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%dy ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%dz ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%dxInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%dyInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%dzInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%delX ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%delY ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%delZ ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%delXInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%delYInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%delZInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%xEdge ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%yEdge ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%zEdge ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%xCenter ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%yCenter ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( grid%zCenter ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( grid%geometry, 80, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%ox, 1, MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%oy, 1, MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%oz, 1, MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%rotDeg, 1, MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%nzAir, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%nzEarth, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dx(1), size( grid%dx ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dy(1), size( grid%dy ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dz(1), size( grid%dz ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dxInv(1), size( grid%dxInv ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dyInv(1), size( grid%dyInv ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dzInv(1), size( grid%dzInv ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delX(1), size( grid%delX ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delY(1), size( grid%delY ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delZ(1), size( grid%delZ ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delXInv(1), size( grid%delXInv ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delYInv(1), size( grid%delYInv ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delZInv(1), size( grid%delZInv ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%xEdge(1), size( grid%xEdge ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%yEdge(1), size( grid%yEdge ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%zEdge(1), size( grid%zEdge ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%xCenter(1), size( grid%xCenter ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%yCenter(1), size( grid%yCenter ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%zCenter(1), size( grid%zCenter ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%zAirThick, 1, MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%allocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
           class default
              stop "packGridBuffer: Unclassified grid"
           !
        end select
        !
    end subroutine packGridBuffer
    !
    ! UNPACK job_buffer TO job_info STRUCT
    function unpackGridBuffer( index ) result ( grid )
        implicit none
        !
        integer, intent( inout )     :: index
        !
        class( Grid_t ), allocatable :: grid
        !
        integer :: grid_dx , grid_dy, grid_dz, grid_dxInv, grid_dyInv, grid_dzInv, &
                   grid_delX, grid_delY, grid_delZ, grid_delXInv, grid_delYInv, grid_delZInv, &
                   grid_xEdge, grid_yEdge, grid_zEdge, grid_xCenter, grid_yCenter, grid_zCenter
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        !
        select case( grid_derived_type )
           !
           case( grid_3d_sg )
                !
                allocate( Grid3D_SG_t :: grid )
                !
                select type( grid )
                   !
                   class is( Grid3D_SG_t )
                        !
                        ! SIZES NOW
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_dx, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_dy, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_dz, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_dxInv, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_dyInv, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_dzInv, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_delX, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_delY, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_delZ, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_delXInv, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_delYInv, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_delZInv, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_xEdge, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_yEdge, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_zEdge, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_xCenter, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_yCenter, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid_zCenter, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%geometry, 80, MPI_CHARACTER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%ox, 1, MPI_REAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%oy, 1, MPI_REAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%oz, 1, MPI_REAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%rotDeg, 1, MPI_REAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%nx, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%ny, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%nz, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%nzAir, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%nzEarth, 1, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( grid%dx( grid_dx ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dx(1), grid_dx, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%dy( grid_dy ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dy(1), grid_dy, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%dz( grid_dz ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dz(1), grid_dz, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%dxInv( grid_dxInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dxInv(1), grid_dxInv, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%dyInv( grid_dyInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dyInv(1), grid_dyInv, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%dzInv( grid_dzInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dzInv(1), grid_dzInv, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%delX( grid_delX ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delX(1), grid_delX, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%delY( grid_delY ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delY(1), grid_delY, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%delZ( grid_delZ ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delZ(1), grid_delZ, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%delXInv( grid_delXInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delXInv(1), grid_delXInv, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%delYInv( grid_delYInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delYInv(1), grid_delYInv, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%delZInv( grid_delZInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delZInv(1), grid_delZInv, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%xEdge( grid_xEdge ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%xEdge(1), grid_xEdge, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%yEdge( grid_yEdge ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%yEdge(1), grid_yEdge, MPI_REAL,child_comm, ierr )
                        !
                        allocate( grid%zEdge( grid_zEdge ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%zEdge(1), grid_zEdge, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%xCenter( grid_xCenter ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%xCenter(1), grid_xCenter, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%yCenter( grid_yCenter ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%yCenter(1), grid_yCenter, MPI_REAL, child_comm, ierr )
                        !
                        allocate( grid%zCenter( grid_zCenter ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%zCenter(1), grid_zCenter, MPI_REAL, child_comm, ierr )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%zAirThick, 1, MPI_REAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%allocated, 1, MPI_LOGICAL, child_comm, ierr )
                        !
                   class default
                      stop "unpackGridBuffer: Unclassified grid"
                   !
                end select
                !
           case default
              stop "unpackGridBuffer: Unclassified grid"
           !
        end select
        !
    end function unpackGridBuffer
    !
    !
    function allocateMetricElementsBuffer( metric ) result( metric_size_bytes )
        implicit none
        !
        class( MetricElements_t ), allocatable :: metric
        integer :: i, nbytes(2), metric_size_bytes
        !
        metric_size_bytes = 0
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, child_comm, nbytes(1), ierr )
        !
        metric_size_bytes = metric_size_bytes + allocateRVectorBuffer( metric%EdgeLength )
        !
        metric_size_bytes = metric_size_bytes + allocateRVectorBuffer( metric%FaceArea )
        !
        metric_size_bytes = metric_size_bytes + allocateRVectorBuffer( metric%DualFaceArea )
        !
        metric_size_bytes = metric_size_bytes + allocateRVectorBuffer( metric%DualEdgeLength )
        !
        metric_size_bytes = metric_size_bytes + allocateRScalarBuffer( metric%Vnode )
        !
        metric_size_bytes = metric_size_bytes + allocateRScalarBuffer( metric%Vcell )
        !
        metric_size_bytes = metric_size_bytes + allocateRVectorBuffer( metric%Vedge )
        !
        call MPI_PACK_SIZE( 3, MPI_INTEGER, child_comm, nbytes(2), ierr )
        !
        do i = 1, size( nbytes )
           metric_size_bytes = metric_size_bytes + nbytes(i)
        end do
        !
    end function allocateMetricElementsBuffer
    !
    !
    subroutine packMetricBuffer( metric, index )
        implicit none
        !
        class( MetricElements_t ), allocatable :: metric
        integer, intent( inout )               :: index
        !
        !
        select type( metric )
           !
           class is( MetricElements_CSG_t )
                !
                call MPI_PACK( metric_csg, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call packRVectorBuffer( metric%EdgeLength, index )
                !
                call packRVectorBuffer( metric%FaceArea, index )
                !
                call packRVectorBuffer( metric%DualFaceArea, index )
                !
                call packRVectorBuffer( metric%DualEdgeLength, index )
                !
                call packRScalarBuffer( metric%Vnode, index )
                !
                call packRScalarBuffer( metric%Vcell, index )
                !
                call packRVectorBuffer( metric%Vedge, index )
                !
                call MPI_PACK( metric%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( metric%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( metric%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
           class default
              stop "packMetricBuffer: Unclassified metric"
           !
        end select
        !
    end subroutine packMetricBuffer
    !
    !
    function unpackMetricBuffer( grid, index ) result ( metric )
        implicit none
        !
        class( Grid_t ), target, intent( in )  :: grid
        integer, intent( inout )               :: index
        !
        class( MetricElements_t ), allocatable :: metric
        !
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, metric_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        !
        select case( metric_derived_type )
           !
           case( metric_csg )
                !
                if( allocated( metric ) ) deallocate( metric )
                allocate( MetricElements_CSG_t :: metric )
                !
                select type( metric )
                   !
                   class is( MetricElements_CSG_t )
                       !
                       metric%grid => grid
                       !
                       allocate( metric%EdgeLength, source = unpackRVectorBuffer( grid, index ) )
                       !
                       allocate( metric%FaceArea, source = unpackRVectorBuffer( grid, index ) )
                       !
                       allocate( metric%DualFaceArea, source = unpackRVectorBuffer( grid, index ) )
                       !
                       allocate( metric%DualEdgeLength, source = unpackRVectorBuffer( grid, index ) )
                       !
                       allocate( metric%Vnode, source = unpackRScalarBuffer( grid, index ) )
                       !
                       allocate( metric%Vcell, source = unpackRScalarBuffer( grid, index ) )
                       !
                       allocate( metric%Vedge, source = unpackRVectorBuffer( grid, index ) )
                       !
                       call MPI_UNPACK( shared_buffer, shared_buffer_size, index, metric%nx, 1, MPI_INTEGER, child_comm, ierr )
                       call MPI_UNPACK( shared_buffer, shared_buffer_size, index, metric%ny, 1, MPI_INTEGER, child_comm, ierr )
                       call MPI_UNPACK( shared_buffer, shared_buffer_size, index, metric%nz, 1, MPI_INTEGER, child_comm, ierr )
                       !
                   class default
                      stop "unpackMetricBuffer: Unclassified metric"
                   !
                end select
                !
           case default
              stop "unpackMetricBuffer: Unclassified metric"
           !
        end select
        !
    end function unpackMetricBuffer
    !
    !
    function allocateModelOperatorBuffer( model_operator ) result( model_operator_size_bytes )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        integer :: i, nbytes(19), model_operator_size_bytes
        !
        model_operator_size_bytes = 0
        !
        select type( model_operator )
           !
           class is( ModelOperator_MF_t )
                !
                ! SIZES FOR THE FUTURE
                call MPI_PACK_SIZE( 31, MPI_INTEGER, child_comm, nbytes(1), ierr )
                !
                model_operator_size_bytes = model_operator_size_bytes + allocateMetricElementsBuffer( model_operator%metric )
                !
                call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 2, MPI_LOGICAL, child_comm, nbytes(3), ierr )
                call MPI_PACK_SIZE( 11, MPI_INTEGER, child_comm, nbytes(4), ierr )
                call MPI_PACK_SIZE( size( model_operator%xXY ), MPI_REAL, child_comm, nbytes(5), ierr )
                call MPI_PACK_SIZE( size( model_operator%xXZ ), MPI_REAL, child_comm, nbytes(6), ierr )
                call MPI_PACK_SIZE( size( model_operator%xY ), MPI_REAL, child_comm, nbytes(7), ierr )
                call MPI_PACK_SIZE( size( model_operator%xZ ), MPI_REAL, child_comm, nbytes(8), ierr )
                call MPI_PACK_SIZE( size( model_operator%xXO ), MPI_REAL, child_comm, nbytes(9), ierr )
                call MPI_PACK_SIZE( size( model_operator%yYX ), MPI_REAL, child_comm, nbytes(10), ierr )
                call MPI_PACK_SIZE( size( model_operator%yYZ ), MPI_REAL, child_comm, nbytes(11), ierr )
                call MPI_PACK_SIZE( size( model_operator%yX ), MPI_REAL, child_comm, nbytes(12), ierr )
                call MPI_PACK_SIZE( size( model_operator%yZ ), MPI_REAL, child_comm, nbytes(13), ierr )
                call MPI_PACK_SIZE( size( model_operator%yYO ), MPI_REAL, child_comm, nbytes(14), ierr )
                call MPI_PACK_SIZE( size( model_operator%zZX ), MPI_REAL, child_comm, nbytes(15), ierr )
                call MPI_PACK_SIZE( size( model_operator%zZY ), MPI_REAL, child_comm, nbytes(16), ierr )
                call MPI_PACK_SIZE( size( model_operator%zX ), MPI_REAL, child_comm, nbytes(17), ierr )
                call MPI_PACK_SIZE( size( model_operator%zY ), MPI_REAL, child_comm, nbytes(18), ierr )
                call MPI_PACK_SIZE( size( model_operator%zZO ), MPI_REAL, child_comm, nbytes(19), ierr )
                !
                model_operator_size_bytes = model_operator_size_bytes + allocateRVectorBuffer( model_operator%Sigma_E )
                !
                model_operator_size_bytes = model_operator_size_bytes + allocateRVectorBuffer( model_operator%db1 )
                !
                model_operator_size_bytes = model_operator_size_bytes + allocateRVectorBuffer( model_operator%db2 )
                !
                model_operator_size_bytes = model_operator_size_bytes + allocateRScalarBuffer( model_operator%c )
                !
                do i = 1, size( nbytes )
                   model_operator_size_bytes = model_operator_size_bytes + nbytes(i)
                end do
                !
           class default
              stop "allocateModelOperatorBuffer: Unclassified model_operator"
           !
        end select
        !
    end function allocateModelOperatorBuffer
    !
    function BiArrayToArray( d2_array ) result ( d1_array )
        !
        real( kind=prec ), intent( in ) :: d2_array(:,:)
        real( kind=prec ), allocatable  :: d1_array(:)
        !
        integer :: size_array, dim_d2_array(2)
        !
        dim_d2_array = shape( d2_array )
        size_array = product( dim_d2_array )
        !
        allocate( d1_array( size_array ) )
        d1_array = (/reshape( d1_array, (/size_array, 1/))/)
        !
   end function BiArrayToArray
   !
   !
   function arrayToBiArray( d1_array, x, y ) result ( d2_array )
        !
        real( kind=prec ), intent( in ) :: d1_array(:)
        integer, intent( in ) :: x, y
        !
        real( kind=prec ), allocatable  :: d2_array(:,:)
        !
        allocate( d2_array( x, y ) )
        !
        d2_array = reshape( d1_array, (/x, y/) )
        !
   end function arrayToBiArray
   !
   !
    ! PACK job_info STRUCT TO Grid Object
    subroutine packModelOperatorBuffer( grid, model_operator, index )
        implicit none
        !
        class( Grid_t ), intent( in )          :: grid
        class( ModelOperator_t ), intent( in ) :: model_operator
        integer, intent( inout )      :: index
        !
        real( kind=prec ), allocatable    :: r_array(:)
        !
        select type( model_operator )
           !
           class is( ModelOperator_MF_t )
                !
                call MPI_PACK( model_operator_mf, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                ! SIZES FOR THE FUTURE
                call MPI_PACK( size( model_operator%xXY, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%xXY, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%xXZ, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%xXZ, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%xY, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%xY, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%xZ, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%xZ, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%xXO, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%xXO, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%yYX, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%yYX, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%yYZ, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%yYZ, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%yX, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%yX, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%yZ, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%yZ, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%yYO, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%yYO, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%zZX, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%zZX, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%zZY, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%zZY, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%zX, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%zX, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%zY, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%zY, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%zZO, dim=1 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( model_operator%zZO, dim=2 ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call packMetricBuffer( model_operator%metric, index )
                !
                call MPI_PACK( model_operator%is_allocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_operator%eqset, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_operator%mKey(1), 8, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_operator%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_operator%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_operator%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%xXY )
                call MPI_PACK( r_array(1), product( shape( model_operator%xXY ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%xXZ )
                call MPI_PACK( r_array(1), product( shape( model_operator%xXZ ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%xY )
                call MPI_PACK( r_array(1), product( shape( model_operator%xY ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%xZ )
                call MPI_PACK( r_array(1), product( shape( model_operator%xZ ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%xXO )
                call MPI_PACK( r_array(1), product( shape( model_operator%xXO ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%yYX )
                call MPI_PACK( r_array(1), product( shape( model_operator%yYX ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%yYZ )
                call MPI_PACK( r_array(1), product( shape( model_operator%yYZ ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%yX )
                call MPI_PACK( r_array(1), product( shape( model_operator%yX ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%yZ )
                call MPI_PACK( r_array(1), product( shape( model_operator%yZ ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%yYO )
                call MPI_PACK( r_array(1), product( shape( model_operator%yYO ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%zZX )
                call MPI_PACK( r_array(1), product( shape( model_operator%zZX ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%zZY )
                call MPI_PACK( r_array(1), product( shape( model_operator%zZY ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%zX )
                call MPI_PACK( r_array(1), product( shape( model_operator%zX ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%zY )
                call MPI_PACK( r_array(1), product( shape( model_operator%zY ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%zZO )
                call MPI_PACK( r_array(1), product( shape( model_operator%zZO ) ), MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call packRVectorBuffer( model_operator%Sigma_E, index )
                !
                call packRVectorBuffer( model_operator%db1, index )
                !
                call packRVectorBuffer( model_operator%db2, index )
                !
                call packRScalarBuffer( model_operator%c, index )
                !
           class default
              stop "allocateModelOperatorBuffer: Unclassified model_operator"
           !
        end select

    end subroutine packModelOperatorBuffer
    !
    ! UNPACK job_buffer TO job_info STRUCT
    function unpackModelOperatorBuffer( grid, index ) result ( model_operator )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        integer, intent( inout )     :: index
        !
        class( ModelOperator_t ), allocatable :: model_operator
        !
        real( kind=prec ), allocatable    :: r_array(:)
        !
        integer :: model_operator_xXY_dim1, model_operator_xXY_dim2, model_operator_xXZ_dim1
        integer :: model_operator_xXZ_dim2, model_operator_xY_dim1, model_operator_xY_dim2
        integer :: model_operator_xZ_dim1, model_operator_xZ_dim2, model_operator_xXO_dim1
        integer :: model_operator_xXO_dim2, model_operator_yYX_dim1, model_operator_yYX_dim2
        integer :: model_operator_yYZ_dim1, model_operator_yYZ_dim2, model_operator_yX_dim1
        integer :: model_operator_yX_dim2, model_operator_yZ_dim1, model_operator_yZ_dim2
        integer :: model_operator_yYO_dim1, model_operator_yYO_dim2, model_operator_zZX_dim1
        integer :: model_operator_zZX_dim2, model_operator_zZY_dim1, model_operator_zZY_dim2
        integer :: model_operator_zX_dim1, model_operator_zX_dim2, model_operator_zY_dim1
        integer :: model_operator_zY_dim2, model_operator_zZO_dim1, model_operator_zZO_dim2
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        !
        ! SIZES NOW
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_xXY_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_xXY_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_xXZ_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_xXZ_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_xY_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_xY_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_xZ_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_xZ_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_xXO_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_xXO_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_yYX_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_yYX_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_yYZ_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_yYZ_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_yX_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_yX_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_yZ_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_yZ_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_yYO_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_yYO_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_zZX_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_zZX_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_zZY_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_zZY_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_zX_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_zX_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_zY_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_zY_dim2, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_zZO_dim1, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator_zZO_dim2, 1, MPI_INTEGER, child_comm, ierr )
        !
        select case( model_operator_derived_type )
           !
           case ( model_operator_mf )
                !
                allocate( ModelOperator_MF_t :: model_operator )
                !
                select type( model_operator )
                !
                   class is( ModelOperator_MF_t )
                        !
                        model_operator%grid => grid
                        !
                        model_operator%metric = unpackMetricBuffer( grid, index )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator%is_allocated, 1, MPI_LOGICAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator%eqset, 1, MPI_LOGICAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator%mKey(1), 8, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator%nx, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator%ny, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_operator%nz, 1, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( r_array( model_operator_xXY_dim1 * model_operator_xXY_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_xXY_dim1 * model_operator_xXY_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%xXY = arrayToBiArray( r_array, model_operator_xXY_dim1, model_operator_xXY_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_xXZ_dim1 * model_operator_xXZ_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_xXZ_dim1 * model_operator_xXZ_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%xXZ = arrayToBiArray( r_array, model_operator_xXZ_dim1, model_operator_xXZ_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_xY_dim1 * model_operator_xY_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_xY_dim1 * model_operator_xY_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%xY = arrayToBiArray( r_array, model_operator_xY_dim1, model_operator_xY_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_xZ_dim1 * model_operator_xZ_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_xZ_dim1 * model_operator_xZ_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%xZ = arrayToBiArray( r_array, model_operator_xZ_dim1, model_operator_xZ_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_xXO_dim1 * model_operator_xXO_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_xXO_dim1 * model_operator_xXO_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%xXO = arrayToBiArray( r_array, model_operator_xXO_dim1, model_operator_xXO_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_yYX_dim1 * model_operator_yYX_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_yYX_dim1 * model_operator_yYX_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%yYX = arrayToBiArray( r_array, model_operator_yYX_dim1, model_operator_yYX_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_yYZ_dim1 * model_operator_yYZ_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_yYZ_dim1 * model_operator_yYZ_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%yYZ = arrayToBiArray( r_array, model_operator_yYZ_dim1, model_operator_yYZ_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_yX_dim1 * model_operator_yX_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_yX_dim1 * model_operator_yX_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%yX = arrayToBiArray( r_array, model_operator_yX_dim1, model_operator_yX_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_yZ_dim1 * model_operator_yZ_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_yZ_dim1 * model_operator_yZ_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%yZ = arrayToBiArray( r_array, model_operator_yZ_dim1, model_operator_yZ_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_yYO_dim1 * model_operator_yYO_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_yYO_dim1 * model_operator_yYO_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%yYO = arrayToBiArray( r_array, model_operator_yYO_dim1, model_operator_yYO_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_zZX_dim1 * model_operator_zZX_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_zZX_dim1 * model_operator_zZX_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%zZX = arrayToBiArray( r_array, model_operator_zZX_dim1, model_operator_zZX_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_zZY_dim1 * model_operator_zZY_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_zZY_dim1 * model_operator_zZY_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%zZY = arrayToBiArray( r_array, model_operator_zZY_dim1, model_operator_zZY_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_zX_dim1 * model_operator_zX_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_zX_dim1 * model_operator_zX_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%zX = arrayToBiArray( r_array, model_operator_zX_dim1, model_operator_zX_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_zY_dim1 * model_operator_zY_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_zY_dim1 * model_operator_zY_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%zY = arrayToBiArray( r_array, model_operator_zY_dim1, model_operator_zY_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_zZO_dim1 * model_operator_zZO_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_zZO_dim1 * model_operator_zZO_dim2, MPI_REAL, child_comm, ierr )
                        model_operator%zZO = arrayToBiArray( r_array, model_operator_zZO_dim1, model_operator_zZO_dim2 )
                        deallocate( r_array )
                        !
                        model_operator%Sigma_E = unpackRVectorBuffer( grid, index )
                        !
                        model_operator%db1 = unpackRVectorBuffer( grid, index )
                        !
                        model_operator%db2 = unpackRVectorBuffer( grid, index )
                        !
                        model_operator%c = unpackRScalarBuffer( grid, index )
                        !
                    class default
                        stop "unpackModelOperatorBuffer: Unclassified model_operator"
                    !
                end select
                !
           case default
              stop "unpackModelOperatorBuffer: Unclassified model_operator"
           !
        end select
        !
    end function unpackModelOperatorBuffer
    !
    !
    function allocateModelParameterBuffer( model_parameter ) result( model_parameter_size_bytes )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: model_parameter
        integer :: i, nbytes(4), model_parameter_size_bytes
        !
        model_parameter_size_bytes = 0
        !
        select type( model_parameter )
           !
           class is( ModelParameterCell_SG_t )
                !
                call MPI_PACK_SIZE( 9, MPI_INTEGER, child_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( 80, MPI_CHARACTER, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 1, MPI_REAL, child_comm, nbytes(3), ierr )
                call MPI_PACK_SIZE( 2, MPI_LOGICAL, child_comm, nbytes(4), ierr )
                !
                model_parameter_size_bytes = model_parameter_size_bytes + allocateGridBuffer( model_parameter%paramGrid )
                !
                model_parameter_size_bytes = model_parameter_size_bytes + allocateRScalarBuffer( model_parameter%cellCond )
                !
                do i = 1, size( nbytes )
                   model_parameter_size_bytes = model_parameter_size_bytes + nbytes(i)
                end do
                !
           class default
              stop "allocateModelParameterBuffer: Unclassified model_parameter"
           !
        end select
        !
    end function allocateModelParameterBuffer
    !
    !
    subroutine packModelParameterBuffer( model_parameter, index )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: model_parameter
        integer, intent( inout )                :: index
        !
        select type( model_parameter )
           !
           class is( ModelParameterCell_SG_t )
                !
                call MPI_PACK( model_parameter_cell_sg, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( model_parameter%mKey(1), 8, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_parameter%paramType, 80, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_parameter%airCond, 1, MPI_REAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_parameter%zeroValued, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_parameter%isAllocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call packGridBuffer( model_parameter%paramGrid, index )
                !
                call packRScalarBuffer( model_parameter%cellCond, index )
                !
           class default
              stop "allocateModelParameterBuffer: Unclassified model_parameter"
           !
        end select

    end subroutine packModelParameterBuffer
    !
    ! UNPACK job_buffer TO job_info STRUCT
    function unpackModelParameterBuffer( grid, metric, index ) result ( model_parameter )
        implicit none
        !
        class( Grid_t ), target, intent( in )           :: grid
        class( MetricElements_t ), target, intent( in ) :: metric
        integer, intent( inout )                        :: index
        !
        class( ModelParameter_t ), allocatable :: model_parameter
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        !
        select case( model_parameter_derived_type )
           !
           case ( model_parameter_cell_sg )
                !
                allocate( ModelParameterCell_SG_t :: model_parameter )
                !
                select type( model_parameter )
                !
                   class is( ModelParameterCell_SG_t )
                        !
                        model_parameter%grid => grid
                        !
                        model_parameter%metric => metric
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%mKey(1), 8, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%paramType, 80, MPI_CHARACTER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%airCond, 1, MPI_REAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%zeroValued, 1, MPI_LOGICAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%isAllocated, 1, MPI_LOGICAL, child_comm, ierr )
                        !
                        select type( param_grid => unpackGridBuffer( index ) )
                        !
                           class is( Grid3D_SG_t )
                                !
                                model_parameter%paramGrid = param_grid
                                !
                                model_parameter%cellCond = unpackRScalarBuffer( param_grid, index )
                                !
                            class default
                                stop "unpackModelParameterBuffer: Unclassified grid"
                            !
                        end select
                        !
                    class default
                        stop "unpackModelParameterBuffer: Unclassified model_parameter"
                    !
                end select
                !
           case default
              stop "unpackModelParameterBuffer: Unclassified model_parameter"
           !
        end select
        !
    end function unpackModelParameterBuffer
    !
    ! ALLOCATE job_buffer
    subroutine allocateJobBuffer
        !
        integer nbytes1, nbytes2
        !
        call MPI_PACK_SIZE( 70, MPI_CHARACTER, child_comm, nbytes1, ierr )
        call MPI_PACK_SIZE( 2, MPI_INTEGER, child_comm, nbytes2, ierr )
        !
        job_size_bytes = ( nbytes1 + nbytes2 ) + 1
        !
        if( .NOT. associated( job_buffer ) ) then
            allocate( job_buffer( job_size_bytes ) )
        endif
        !
    end subroutine allocateJobBuffer
    !
    ! PACK job_info STRUCT TO job_buffer
    subroutine packJobBuffer
        !
        integer :: index
        !
        index = 1
        !
        call MPI_PACK( job_info%name, 50, MPI_CHARACTER, job_buffer, job_size_bytes, index, child_comm, ierr )
        call MPI_PACK( job_info%transmitter_derived_type, 20, MPI_CHARACTER, job_buffer, job_size_bytes, index, child_comm, ierr )
        call MPI_PACK( job_info%tx_index, 1, MPI_INTEGER, job_buffer, job_size_bytes, index, child_comm, ierr )
        call MPI_PACK( job_info%worker_rank, 1, MPI_INTEGER, job_buffer, job_size_bytes, index, child_comm, ierr )
        !
    end subroutine packJobBuffer
    !
    ! UNPACK job_buffer TO job_info STRUCT
    subroutine unpackJobBuffer
        !
        integer :: index
        !
        index = 1
        !
        call MPI_UNPACK( job_buffer, job_size_bytes, index, job_info%name, 50, MPI_CHARACTER, child_comm, ierr )
        call MPI_UNPACK( job_buffer, job_size_bytes, index, job_info%transmitter_derived_type, 20, MPI_CHARACTER, child_comm, ierr )
        call MPI_UNPACK( job_buffer, job_size_bytes, index, job_info%tx_index , 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( job_buffer, job_size_bytes, index, job_info%worker_rank , 1, MPI_INTEGER, child_comm, ierr )
        !
    end subroutine unpackJobBuffer
    !
    ! RECEIVE job_info FROM ANY TARGET
    subroutine receiveFromAny()
        !
        write( *, * ) "<<<< ", mpi_rank, " JOB: ", job_info%name, " FROM: ", job_info%worker_rank
        !
        call allocateJobBuffer
        call MPI_RECV( job_buffer, job_size_bytes, MPI_PACKED, MPI_ANY_SOURCE, tag, child_comm, MPI_STATUS_IGNORE, ierr )
        call unpackJobBuffer
        !
    end subroutine receiveFromAny
    !
    ! RECEIVE job_info FROM target_id
    subroutine receiveFrom( target_id )
        !
        integer, intent( in )    :: target_id
        !
        write( *, * ) "<<<< ", mpi_rank, " JOB: ", job_info%name, " FROM: ", target_id
        !
        call allocateJobBuffer
        call MPI_RECV( job_buffer, job_size_bytes, MPI_PACKED, target_id, tag, child_comm, MPI_STATUS_IGNORE, ierr )
        call unpackJobBuffer
        !
    end subroutine receiveFrom
    !
    ! SEND job_info FROM target_id
    subroutine sendTo( target_id )
        !
        integer, intent( in )    :: target_id
        !
        write( *, * ) ">>>> ", mpi_rank, " JOB: ", job_info%name, " TO: ", target_id
        !
        call allocateJobBuffer
        call packJobBuffer
        call MPI_SEND( job_buffer, job_size_bytes, MPI_PACKED, target_id, tag, child_comm, ierr )
        !
    end subroutine sendTo
    !
end module DeclarationMPI