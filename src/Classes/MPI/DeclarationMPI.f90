!*************
!
! MPI declarations
!
!*************
!
module DeclarationMPI
    !
    use, intrinsic :: iso_c_binding, only: c_ptr, c_sizeof, c_f_pointer
    !
    use Constants
    use FileUnits
    !
    use ModEMControlFile
    !
    use Grid3D_SG
    !
    use ModelReader
    use ModelReader_Weerachai
    use ModelOperator_MF
    use ModelOperator_File
    use ModelParameterCell_SG
    !
    use DataFileStandard
    !
    use ForwardSolverIT_DC
    !
    use SourceMT_1D
    use SourceMT_2D
    use SourceCSEM_Dipole1D
    !
    use DataHandleFArray
    use DataHandleMT
    use DataHandleCSEM
    !
    !
    include 'mpif.h'
    !
    integer :: main_comm, child_comm
    !
    integer :: mpi_rank, mpi_size, node_rank, node_size, &
               nodestringlen, ierr
    !
    character*( MPI_MAX_PROCESSOR_NAME ) :: node_name
    !
    integer                          :: shared_window
    integer( KIND=MPI_ADDRESS_KIND ) :: shared_window_size
    integer                          :: shared_disp_unit = 1
    !
    character, dimension(:), pointer     :: shared_buffer
    !
    character, dimension(:), allocatable :: fwd_info_buffer, predicted_data_buffer
    !
    integer :: shared_buffer_size = 1
    integer :: fwd_info_buffer_size = 1
    integer :: predicted_data_buffer_size = 1
    !
    type( c_ptr ) :: shared_c_ptr
    !
    integer :: field_derived_type
    integer, parameter :: real_3d_sg = 1
    integer, parameter :: complex_3d_sg = 2
    !
    integer :: grid_derived_type
    integer, parameter :: grid_3d_sg = 1
    !
    integer :: model_parameter_derived_type
    integer, parameter :: model_parameter_cell_sg = 1
    !
    integer :: transmitter_derived_type, transmitters_size
    integer, parameter :: transmitter_mt = 1
    integer, parameter :: transmitter_csem = 2
    !
    integer :: receiver_derived_type, receivers_size
    integer, parameter :: receiver_full_impedance = 1
    integer, parameter :: receiver_full_vertical_magnetic = 2
    integer, parameter :: receiver_off_diagonal_impedance = 3
    integer, parameter :: receiver_single_field = 4
    !
    integer :: data_derived_type
    integer, parameter :: data_mt = 1
    integer, parameter :: data_csem = 2
    !
    class( Grid_t ), allocatable, target   :: main_grid
    class( ModelParameter_t ), allocatable :: model_parameter
    !
    ! PROGRAM GLOBAL VARIABLES
    integer :: tag = 2022, master_id = 0
    !
    character*15    :: job_master = "MASTER_JOB", job_finish = "STOP_JOBS", job_fwd_done = "FINISH_FWD_JOB"
    character*15    :: job_share_memory = "SHARE_MEMORY", job_forward = "JOB_FORWARD"
    !
    ! STRUCT job_info
    type :: FWDInfo_t
        !
        SEQUENCE
        !
        character*15 :: job_name
        integer      :: worker_rank
        integer      :: tx_index
        integer      :: n_data
        integer      :: data_size
        !
    end type FWDInfo_t
    !
    type( FWDInfo_t ), save :: fwd_info
    !
    !
    contains
    !
    ! ALLOCATE shared_buffer
    subroutine allocateSharedBuffer()
        implicit none
        !
        integer :: i, last_size, nbytes(8)
        !
        call MPI_PACK_SIZE( 4, MPI_INTEGER, child_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 2, MPI_DOUBLE_PRECISION, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( len( e_solution_file_name ), MPI_CHARACTER, child_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( len( model_method ), MPI_CHARACTER, child_comm, nbytes(4), ierr )
        call MPI_PACK_SIZE( len( forward_solver_type ), MPI_CHARACTER, child_comm, nbytes(5), ierr )
        call MPI_PACK_SIZE( len( source_type ), MPI_CHARACTER, child_comm, nbytes(6), ierr )
        !
        shared_buffer_size = shared_buffer_size + allocateGridBuffer( main_grid )
        !
        write( *, "(A60, i8)" ) "MPI Allocated main_grid size:", shared_buffer_size
        last_size = shared_buffer_size
        !
        shared_buffer_size = shared_buffer_size + allocateModelParameterBuffer()
        !
        write( *, "(A60, i8)" ) "MPI Allocated model_parameter size:", shared_buffer_size - last_size
        last_size = shared_buffer_size
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, child_comm, nbytes(7), ierr )
        !
        do i = 1, size( transmitters )
            shared_buffer_size = shared_buffer_size + allocateTransmitterBuffer( getTransmitter( i ) )
            !
            !write( *, "(A50, i8)" ) "MPI Allocated tx size: ", shared_buffer_size - last_size
            !last_size = shared_buffer_size
            !
        end do
        !
        write( *, "(A60, i8)" ) "MPI Allocated transmitters size:", shared_buffer_size - last_size
        last_size = shared_buffer_size
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, child_comm, nbytes(8), ierr )
        !
        do i = 1, size( receivers )
            shared_buffer_size = shared_buffer_size + allocateReceiverBuffer( getReceiver(i) )
            !
            !write( *, "(A60, i8)" ) "MPI Allocated rx size:", shared_buffer_size - last_size
            !last_size = shared_buffer_size
            !
        end do
        !
        write( *, "(A60, i8)" ) "MPI Allocated receivers size:", shared_buffer_size - last_size
        !
        do i = 1, size( nbytes )
            shared_buffer_size = shared_buffer_size + nbytes(i)
        end do
        !
        write( *, "(A60, i8)" ) "MPI Allocated total size =", shared_buffer_size
        !
        allocate( shared_buffer( shared_buffer_size ) )
        !
    end subroutine allocateSharedBuffer
    !
    !
    subroutine packSharedBuffer()
        implicit none
        !
        integer :: i, index
        !
        index = 1
        !
        call MPI_PACK( len( e_solution_file_name ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( len( model_method ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( len( forward_solver_type ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( len( source_type ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( model_n_air_layer, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( model_max_height, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( e_solution_file_name, len( e_solution_file_name ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( model_method, len( model_method ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( forward_solver_type, len( forward_solver_type ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( source_type, len( source_type ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        !
        call packGridBuffer( main_grid, index )
        !
        call packModelParameterBuffer( index )
        !
        call MPI_PACK( size( transmitters ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        !
        do i = 1, size( transmitters )
            !
            call packTransmitterBuffer( getTransmitter( i ), index )
            !
        end do
        !
        call MPI_PACK( size( receivers ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        !
        do i = 1, size( receivers )
            !
            call packReceiverBuffer( getReceiver( i ), index )
            !
        end do
        !
    end subroutine packSharedBuffer
    !
    subroutine unpackSharedBuffer( buffer_size )
        implicit none
        !
        integer, intent( in ) :: buffer_size
        !
        integer :: i, aux_size, n_e_solution_file_name, n_model_method, n_forward_solver_type, n_source_type, index
        !
        index = 1
        !
        shared_buffer_size = buffer_size
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, n_e_solution_file_name, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, n_model_method, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, n_forward_solver_type, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, n_source_type, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_n_air_layer, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_max_height, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
        !
        allocate( character( n_e_solution_file_name ) :: e_solution_file_name )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, e_solution_file_name, n_e_solution_file_name, MPI_CHARACTER, child_comm, ierr )
        !
        allocate( character( n_model_method ) :: model_method )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_method, n_model_method, MPI_CHARACTER, child_comm, ierr )
        !
        allocate( character( n_forward_solver_type ) :: forward_solver_type )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, forward_solver_type, n_forward_solver_type, MPI_CHARACTER, child_comm, ierr )
        !
        allocate( character( n_source_type ) :: source_type )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, source_type, n_source_type, MPI_CHARACTER, child_comm, ierr )
        !
        allocate( main_grid, source = unpackGridBuffer( index ) )
        !
        allocate( model_parameter, source = unpackModelParameterBuffer( index ) )
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, aux_size, 1, MPI_INTEGER, child_comm, ierr )
        !
        do i = 1, aux_size
            !
            call updateTransmitterArray( unpackTransmitterBuffer( index ) )
            !
        end do
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, aux_size, 1, MPI_INTEGER, child_comm, ierr )
        !
        do i = 1, aux_size
            !
            ierr = updateReceiverArray( unpackReceiverBuffer( index ) )
            !
        end do
        !
    end subroutine unpackSharedBuffer
    !
    !
    function allocateRScalarBuffer( scalar ) result( scalar_size_bytes )
        implicit none
        !
        class( rScalar_t ), intent( in ):: scalar
        !
        integer :: i, nbytes(4), scalar_size_bytes
        !
        scalar_size_bytes = 0
        !
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, child_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 7, MPI_INTEGER, child_comm, nbytes(3), ierr )
        !
        select type( scalar )
           !
           class is( rScalar3D_SG_t )
              !
              call MPI_PACK_SIZE( scalar%Nxyz, MPI_DOUBLE_PRECISION, child_comm, nbytes(4), ierr )
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
        integer, intent( inout )         :: index
        !
        real( kind=prec ), allocatable   :: r_array(:)
        !
        select type( scalar )
           !
           class is( rScalar3D_SG_t )
              !
              call MPI_PACK( scalar%gridType, 4, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%is_allocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%NdV(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%Nxyz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call scalar%getArray( r_array )
              call MPI_PACK( r_array(1), scalar%Nxyz, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
        class( Grid_t ), intent( in ) :: grid
        integer, intent( inout )      :: index
        !
        class( rScalar_t ), allocatable :: scalar
        !
        real( kind=prec ), allocatable  :: r_array(:)
        !
        character( len=4 ) :: gridType
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, gridType, 4, MPI_CHARACTER, child_comm, ierr )
        !
        select type( grid )
            !
            class is( Grid3D_SG_t )
                !
                allocate( scalar, source = rScalar3D_SG_t( grid, gridType ) )
                !
            class default
                stop "unpackRScalarBuffer: Unclassified grid"
                !
        end select
        !
        select type( scalar )
            !
            class is ( rScalar3D_SG_t )
                !
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%is_allocated, 1, MPI_LOGICAL, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%nx, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%ny, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%nz, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%NdV(1), 3, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%Nxyz, 1, MPI_INTEGER, child_comm, ierr )
                !
                allocate( r_array( scalar%Nxyz ) )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), scalar%Nxyz, MPI_DOUBLE_PRECISION, child_comm, ierr )
                call scalar%setArray( r_array )
                !
            class default
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
        integer :: i, nbytes(4), scalar_size_bytes
        !
        scalar_size_bytes = 0
        !
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, child_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 7, MPI_INTEGER, child_comm, nbytes(3), ierr )
        !
        select type( scalar )
           !
           class is( cScalar3D_SG_t )
              !
              call MPI_PACK_SIZE( scalar%Nxyz, MPI_DOUBLE_COMPLEX, child_comm, nbytes(4), ierr )
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
              call MPI_PACK( scalar%gridType, 4, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%is_allocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%NdV(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( scalar%Nxyz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call scalar%getArray( c_array )
              call MPI_PACK( c_array(1), scalar%Nxyz, MPI_DOUBLE_COMPLEX, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
        class( Grid_t ), intent( in ) :: grid
        integer, intent( inout )      :: index
        !
        class( cScalar_t ), allocatable :: scalar
        !
        complex( kind=prec ), allocatable    :: c_array(:)
        !
        character( len=4 ) :: gridType
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, gridType, 4, MPI_CHARACTER, child_comm, ierr )
        !
        select type( grid )
            !
            class is( Grid3D_SG_t )
                !
                allocate( scalar, source = cScalar3D_SG_t( grid, gridType ) )
                !
            class default
                stop "unpackCScalarBuffer: Unclassified grid"
                !
        end select
        !
        select type( scalar )
            !
            class is ( cScalar3D_SG_t )
                !
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%is_allocated, 1, MPI_LOGICAL, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%nx, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%ny, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%nz, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%NdV(1), 3, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%Nxyz, 1, MPI_INTEGER, child_comm, ierr )
                !
                allocate( c_array( scalar%Nxyz ) )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, c_array(1), scalar%Nxyz, MPI_DOUBLE_COMPLEX, child_comm, ierr )
                call scalar%setArray( c_array )
                !
            class default
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
        integer :: i, nbytes(4), vector_size_bytes
        !
        vector_size_bytes = 0
        !
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, child_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 15, MPI_INTEGER, child_comm, nbytes(3), ierr )
        !
        select type( vector )
           !
           class is( rVector3D_SG_t )
              !
              call MPI_PACK_SIZE( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_PRECISION, child_comm, nbytes(4), ierr )
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
        integer, intent( inout )         :: index
        !
        real( kind=prec ), allocatable   :: r_array(:)
        !
        select type( vector )
           !
           class is( rVector3D_SG_t )
              !
              call MPI_PACK( vector%gridType, 4, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%is_allocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdX(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdY(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdZ(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%Nxyz(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call vector%getArray( r_array )
              call MPI_PACK( r_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
        class( Grid_t ), intent( in ) :: grid
        integer, intent( inout )      :: index
        !
        class( rVector_t ), allocatable :: vector
        !
        real( kind=prec ), allocatable  :: r_array(:)
        !
        character( len=4 ) :: gridType
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, gridType, 4, MPI_CHARACTER, child_comm, ierr )
        !
        select type( grid )
            !
            class is( Grid3D_SG_t )
                !
                allocate( vector, source = rVector3D_SG_t( grid, gridType ) )
                !
            class default
                stop "unpackRVectorBuffer: Unclassified grid"
                !
        end select
        !
        select type( vector )
            !
            class is ( rVector3D_SG_t )
                !
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%is_allocated, 1, MPI_LOGICAL, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%nx, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%ny, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%nz, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdX(1), 3, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdY(1), 3, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdZ(1), 3, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%Nxyz(1), 3, MPI_INTEGER, child_comm, ierr )
                !
                allocate( r_array( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3) ) )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_PRECISION, child_comm, ierr )
                call vector%setArray( r_array )
                !
            class default
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
        integer :: i, nbytes(4), vector_size_bytes
        !
        vector_size_bytes = 0
        !
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, child_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 15, MPI_INTEGER, child_comm, nbytes(3), ierr )
        !
        select type( vector )
           !
           class is( cVector3D_SG_t )
              !
              call MPI_PACK_SIZE( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_COMPLEX, child_comm, nbytes(4), ierr )
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
              call MPI_PACK( vector%gridType, 4, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%is_allocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdX(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdY(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%NdZ(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              call MPI_PACK( vector%Nxyz(1), 3, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
              !
              call vector%getArray( c_array )
              call MPI_PACK( c_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_COMPLEX, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
        class( Grid_t ), intent( in ) :: grid
        integer, intent( inout )      :: index
        !
        class( cVector_t ), allocatable :: vector
        !
        complex( kind=prec ), allocatable :: c_array(:)
        !
        character( len=4 ) :: gridType
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, gridType, 4, MPI_CHARACTER, child_comm, ierr )
        !
        select type( grid )
            !
            class is( Grid3D_SG_t )
                !
                allocate( vector, source = cVector3D_SG_t( grid, gridType ) )
                !
            class default
                stop "unpackCVectorBuffer: Unclassified grid"
                !
        end select
        !
        select type( vector )
            !
            class is ( cVector3D_SG_t )
                !
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%is_allocated, 1, MPI_LOGICAL, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%nx, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%ny, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%nz, 1, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdX(1), 3, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdY(1), 3, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%NdZ(1), 3, MPI_INTEGER, child_comm, ierr )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%Nxyz(1), 3, MPI_INTEGER, child_comm, ierr )
                !
                allocate( c_array( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3) ) )
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, c_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_COMPLEX, child_comm, ierr )
                call vector%setArray( c_array )
                !
            class default
                stop "unpackCVectorBuffer: Unclassified vector"
            !
        end select
        !
    end function unpackCVectorBuffer
    !
    function allocateCSparseVectorBuffer( vector ) result( vector_size_bytes )
        implicit none
        !
        type( cSparsevector3D_SG_t ), intent( in ):: vector
        !
        integer :: i, nbytes(7), vector_size_bytes
        !
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, vector_size_bytes, ierr )
        !
        if( vector%is_allocated ) then
            !
            call MPI_PACK_SIZE( 4, MPI_CHARACTER, child_comm, nbytes(1), ierr )
            call MPI_PACK_SIZE( 6, MPI_INTEGER, child_comm, nbytes(2), ierr )
            call MPI_PACK_SIZE( size( vector%i ), MPI_INTEGER, child_comm, nbytes(3), ierr )
            call MPI_PACK_SIZE( size( vector%j ), MPI_INTEGER, child_comm, nbytes(4), ierr )
            call MPI_PACK_SIZE( size( vector%k ), MPI_INTEGER, child_comm, nbytes(5), ierr )
            call MPI_PACK_SIZE( size( vector%xyz ), MPI_INTEGER, child_comm, nbytes(6), ierr )
            call MPI_PACK_SIZE( size( vector%c ), MPI_DOUBLE_COMPLEX, child_comm, nbytes(7), ierr )
            !
            do i = 1, size( nbytes )
               vector_size_bytes = vector_size_bytes + nbytes(i)
            end do
            !
        endif
        !
    end function allocateCSparseVectorBuffer
    !
    !
    subroutine packCSparseVectorBuffer( vector, index )
        implicit none
        !
        !
        type( cSparsevector3D_SG_t ), intent( in ) :: vector
        !
        integer, intent( inout ) :: index
        !
        call MPI_PACK( vector%is_allocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        !
        if( vector%is_allocated ) then
            !
            call MPI_PACK( vector%gridType, 4, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( vector%nCoeff, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( size( vector%i ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( size( vector%j ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( size( vector%k ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( size( vector%xyz ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( size( vector%c ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( vector%i(1), size( vector%i ), MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( vector%j(1), size( vector%j ), MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( vector%k(1), size( vector%k ), MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( vector%xyz(1), size( vector%xyz ), MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            call MPI_PACK( vector%c(1), size( vector%c ), MPI_DOUBLE_COMPLEX, shared_buffer, shared_buffer_size, index, child_comm, ierr )
            !
        endif
        !
    end subroutine packCSparseVectorBuffer
    !
    function unpackCSparseVectorBuffer( index ) result( vector )
        implicit none
        !
        integer, intent( inout )     :: index
        type( cSparsevector3D_SG_t ) :: vector
        !
        !
        integer :: vector_n_i, vector_n_j, vector_n_k, vector_n_xyz, vector_n_c
        !
        !
        vector = cSparsevector3D_SG_t()
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%is_allocated, 1, MPI_LOGICAL, child_comm, ierr )
        !
        if( vector%is_allocated ) then
            !
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%gridType, 4, MPI_CHARACTER, child_comm, ierr )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%nCoeff, 1, MPI_INTEGER, child_comm, ierr )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector_n_i, 1, MPI_INTEGER, child_comm, ierr )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector_n_j, 1, MPI_INTEGER, child_comm, ierr )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector_n_k, 1, MPI_INTEGER, child_comm, ierr )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector_n_xyz, 1, MPI_INTEGER, child_comm, ierr )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector_n_c, 1, MPI_INTEGER, child_comm, ierr )
            !
            allocate( vector%i( vector_n_i ) )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%i(1), vector_n_i, MPI_INTEGER, child_comm, ierr )
            !
            allocate( vector%j( vector_n_j ) )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%j(1), vector_n_j, MPI_INTEGER, child_comm, ierr )
            !
            allocate( vector%k( vector_n_k ) )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%k(1), vector_n_k, MPI_INTEGER, child_comm, ierr )
            !
            allocate( vector%xyz( vector_n_xyz ) )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%xyz(1), vector_n_xyz, MPI_INTEGER, child_comm, ierr )
            !
            allocate( vector%c( vector_n_c ) )
            call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%c(1), vector_n_c, MPI_DOUBLE_COMPLEX, child_comm, ierr )
            !
        endif
        !
    end function unpackCSparseVectorBuffer
    !
    ! ALLOCATE shared_buffer
    function allocateGridBuffer( main_grid ) result( grid_size_bytes )
        implicit none
        !
        class( Grid_t ), intent( in ) :: main_grid
        integer :: i, nbytes(24), grid_size_bytes
        !
        grid_size_bytes = 0
        !
        ! SIZES FOR THE FUTURE
        call MPI_PACK_SIZE( 19, MPI_INTEGER, child_comm, nbytes(1), ierr )
        !
        call MPI_PACK_SIZE( 80, MPI_CHARACTER, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 4, MPI_DOUBLE_PRECISION, child_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 5, MPI_INTEGER, child_comm, nbytes(4), ierr )
        call MPI_PACK_SIZE( size( main_grid%dx ), MPI_DOUBLE_PRECISION, child_comm, nbytes(5), ierr )
        call MPI_PACK_SIZE( size( main_grid%dy ), MPI_DOUBLE_PRECISION, child_comm, nbytes(6), ierr )
        call MPI_PACK_SIZE( size( main_grid%dz ), MPI_DOUBLE_PRECISION, child_comm, nbytes(7), ierr )
        call MPI_PACK_SIZE( size( main_grid%dxInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(8), ierr )
        call MPI_PACK_SIZE( size( main_grid%dyInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(9), ierr )
        call MPI_PACK_SIZE( size( main_grid%dzInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(10), ierr )
        call MPI_PACK_SIZE( size( main_grid%delX ), MPI_DOUBLE_PRECISION, child_comm, nbytes(11), ierr )
        call MPI_PACK_SIZE( size( main_grid%delY ), MPI_DOUBLE_PRECISION, child_comm, nbytes(12), ierr )
        call MPI_PACK_SIZE( size( main_grid%delZ ), MPI_DOUBLE_PRECISION, child_comm, nbytes(13), ierr )
        call MPI_PACK_SIZE( size( main_grid%delXInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(14), ierr )
        call MPI_PACK_SIZE( size( main_grid%delYInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(15), ierr )
        call MPI_PACK_SIZE( size( main_grid%delZInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(16), ierr )
        call MPI_PACK_SIZE( size( main_grid%xEdge ), MPI_DOUBLE_PRECISION, child_comm, nbytes(17), ierr )
        call MPI_PACK_SIZE( size( main_grid%yEdge ), MPI_DOUBLE_PRECISION, child_comm, nbytes(18), ierr )
        call MPI_PACK_SIZE( size( main_grid%zEdge ), MPI_DOUBLE_PRECISION, child_comm, nbytes(19), ierr )
        call MPI_PACK_SIZE( size( main_grid%xCenter ), MPI_DOUBLE_PRECISION, child_comm, nbytes(20), ierr )
        call MPI_PACK_SIZE( size( main_grid%yCenter ), MPI_DOUBLE_PRECISION, child_comm, nbytes(21), ierr )
        call MPI_PACK_SIZE( size( main_grid%zCenter ), MPI_DOUBLE_PRECISION, child_comm, nbytes(22), ierr )
        call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, child_comm, nbytes(23), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes(24), ierr )
        !
        do i = 1, size( nbytes )
           grid_size_bytes = grid_size_bytes + nbytes(i)
        end do
        !
    end function allocateGridBuffer
    !
    !
    subroutine packGridBuffer( main_grid, index )
        implicit none
        !
        class( Grid_t ), intent( in ) :: main_grid
        integer, intent( inout )      :: index
        !
        select type( main_grid )
           !
           class is( Grid3D_SG_t )
                !
                call MPI_PACK( grid_3d_sg, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                ! SIZES FOR THE FUTURE
                call MPI_PACK( size( main_grid%dx ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%dy ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%dz ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%dxInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%dyInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%dzInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%delX ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%delY ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%delZ ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%delXInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%delYInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%delZInv ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%xEdge ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%yEdge ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%zEdge ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%xCenter ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%yCenter ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( main_grid%zCenter ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( main_grid%geometry, 80, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%ox, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%oy, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%oz, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%rotDeg, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%nzAir, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%nzEarth, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( main_grid%dx(1), size( main_grid%dx ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%dy(1), size( main_grid%dy ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%dz(1), size( main_grid%dz ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%dxInv(1), size( main_grid%dxInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%dyInv(1), size( main_grid%dyInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%dzInv(1), size( main_grid%dzInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%delX(1), size( main_grid%delX ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%delY(1), size( main_grid%delY ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%delZ(1), size( main_grid%delZ ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%delXInv(1), size( main_grid%delXInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%delYInv(1), size( main_grid%delYInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%delZInv(1), size( main_grid%delZInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%xEdge(1), size( main_grid%xEdge ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%yEdge(1), size( main_grid%yEdge ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%zEdge(1), size( main_grid%zEdge ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%xCenter(1), size( main_grid%xCenter ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%yCenter(1), size( main_grid%yCenter ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%zCenter(1), size( main_grid%zCenter ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( main_grid%zAirThick, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( main_grid%is_allocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
           class default
              stop "packGridBuffer: Unclassified main_grid"
           !
        end select
        !
    end subroutine packGridBuffer
    !
    !
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
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%geometry, 80, MPI_CHARACTER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%ox, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%oy, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%oz, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%rotDeg, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%nx, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%ny, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%nz, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%nzAir, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%nzEarth, 1, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( grid%dx( grid_dx ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dx(1), grid_dx, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%dy( grid_dy ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dy(1), grid_dy, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%dz( grid_dz ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dz(1), grid_dz, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%dxInv( grid_dxInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dxInv(1), grid_dxInv, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%dyInv( grid_dyInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dyInv(1), grid_dyInv, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%dzInv( grid_dzInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%dzInv(1), grid_dzInv, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%delX( grid_delX ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delX(1), grid_delX, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%delY( grid_delY ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delY(1), grid_delY, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%delZ( grid_delZ ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delZ(1), grid_delZ, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%delXInv( grid_delXInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delXInv(1), grid_delXInv, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%delYInv( grid_delYInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delYInv(1), grid_delYInv, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%delZInv( grid_delZInv ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%delZInv(1), grid_delZInv, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%xEdge( grid_xEdge ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%xEdge(1), grid_xEdge, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%yEdge( grid_yEdge ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%yEdge(1), grid_yEdge, MPI_DOUBLE_PRECISION,child_comm, ierr )
                        !
                        allocate( grid%zEdge( grid_zEdge ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%zEdge(1), grid_zEdge, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%xCenter( grid_xCenter ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%xCenter(1), grid_xCenter, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%yCenter( grid_yCenter ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%yCenter(1), grid_yCenter, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( grid%zCenter( grid_zCenter ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%zCenter(1), grid_zCenter, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%zAirThick, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, grid%is_allocated, 1, MPI_LOGICAL, child_comm, ierr )
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
    function allocateModelParameterBuffer() result( model_parameter_size_bytes )
        implicit none
        !
        integer :: i, nbytes(4), model_parameter_size_bytes
        !
        model_parameter_size_bytes = 0
        !
        select type( model_parameter )
           !
           class is( ModelParameterCell_SG_t )
                !
                call MPI_PACK_SIZE( 10, MPI_INTEGER, child_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( len( model_parameter%paramType ), MPI_CHARACTER, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, child_comm, nbytes(3), ierr )
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
    subroutine packModelParameterBuffer( index )
        implicit none
        !
        integer, intent( inout ) :: index
        !
        select type( model_parameter )
           !
           class is( ModelParameterCell_SG_t )
                !
                call MPI_PACK( model_parameter_cell_sg, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( len( model_parameter%paramType ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_parameter%mKey(1), 8, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_parameter%paramType, len( model_parameter%paramType ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_parameter%airCond, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_parameter%zeroValued, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( model_parameter%is_allocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
    !
    function unpackModelParameterBuffer( index ) result ( model_parameter )
        implicit none
        !
        integer, intent( inout ) :: index
        !
        class( ModelParameter_t ), allocatable :: model_parameter
        !
        integer :: param_type_size
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
                        model_parameter%grid => main_grid
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, param_type_size, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%mKey(1), 8, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( character( param_type_size ) :: model_parameter%paramType )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%paramType, param_type_size, MPI_CHARACTER, child_comm, ierr )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%airCond, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%zeroValued, 1, MPI_LOGICAL, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%is_allocated, 1, MPI_LOGICAL, child_comm, ierr )
                        !
                        allocate( model_parameter%paramGrid, source = unpackGridBuffer( index ) )
                        !
                        allocate( rScalar3D_SG_t :: model_parameter%cellCond )
                        model_parameter%cellCond = unpackRScalarBuffer( main_grid, index )
                        !
                        call model_parameter%SetSigMap( model_parameter%paramType )
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
    !
    function allocateTransmitterBuffer( transmitter ) result( transmitter_size_bytes )
        implicit none
        !
        class( Transmitter_t ), intent( in ) :: transmitter
        integer :: i, transmitter_size_bytes
        integer, allocatable, dimension(:) :: nbytes
        !
        transmitter_size_bytes = 0
        !
        select type( transmitter )
           !
           class is( TransmitterMT_t )
                !
                allocate( nbytes(3) )
                !
                call MPI_PACK_SIZE( 12, MPI_INTEGER, child_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( size( transmitter%receiver_indexes ), MPI_INTEGER, child_comm, nbytes(3), ierr )
                !
            class is( TransmitterCSEM_t )
                !
                allocate( nbytes(4) )
                !
                call MPI_PACK_SIZE( 13, MPI_INTEGER, child_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( 7, MPI_DOUBLE_PRECISION, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( size( transmitter%receiver_indexes ), MPI_INTEGER, child_comm, nbytes(3), ierr )
                call MPI_PACK_SIZE( len( transmitter%dipole ), MPI_CHARACTER, child_comm, nbytes(4), ierr )
                !
           class default
              stop "allocateTransmitterBuffer: Unclassified transmitter"
           !
        end select
        !
        do i = 1, size( nbytes )
            transmitter_size_bytes = transmitter_size_bytes + nbytes(i)
        end do
        !
    end function allocateTransmitterBuffer
    !
    !
    subroutine packTransmitterBuffer( transmitter, index )
        implicit none
        !
        class( Transmitter_t ), intent( in ) :: transmitter
        integer, intent( inout )             :: index
        !
        integer :: i
        !
        select type( transmitter )
           !
           class is( TransmitterMT_t )
                !
                ! TYPE
                call MPI_PACK( transmitter_mt, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( transmitter%id, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%n_pol, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%fwd_key(1), 8, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( transmitter%receiver_indexes ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%period, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%receiver_indexes(1), size( transmitter%receiver_indexes ), MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
           class is( TransmitterCSEM_t )
                !
                ! TYPE
                call MPI_PACK( transmitter_csem, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( transmitter%id, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%n_pol, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%fwd_key(1), 8, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( transmitter%receiver_indexes ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( len( transmitter%dipole ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%period, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%location(1), 3, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%azimuth, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%dip, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%moment, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%receiver_indexes(1), size( transmitter%receiver_indexes ), MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%dipole, len( transmitter%dipole ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
           class default
              stop "allocateTransmitterBuffer: Unclassified transmitter"
           !
        end select

    end subroutine packTransmitterBuffer
    !
    !
    function unpackTransmitterBuffer( index ) result ( transmitter )
        implicit none
        !
        integer, intent( inout )            :: index
        !
        class( Transmitter_t ), allocatable :: transmitter
        !
        integer :: transmitter_receiver_indexes, transmitter_dipole
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        !
        select case( transmitter_derived_type )
           !
           case ( transmitter_mt )
                !
                allocate( TransmitterMT_t :: transmitter )
                !
                select type( transmitter )
                   !
                   class is( TransmitterMT_t )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%id, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%n_pol, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%fwd_key(1), 8, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_receiver_indexes, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%period, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( transmitter%receiver_indexes( transmitter_receiver_indexes ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%receiver_indexes(1), transmitter_receiver_indexes, MPI_INTEGER, child_comm, ierr )
                        !
                    class default
                        stop "unpackTransmitterBuffer: Unclassified transmitter 1!"
                    !
                end select
                !
            case ( transmitter_csem )
                !
                allocate( TransmitterCSEM_t :: transmitter )
                !
                select type( transmitter )
                   !
                   class is( TransmitterCSEM_t )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%id, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%n_pol, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%fwd_key(1), 8, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_receiver_indexes, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_dipole, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%period, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%location(1), 3, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%azimuth, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%dip, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%moment, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( transmitter%receiver_indexes( transmitter_receiver_indexes ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%receiver_indexes(1), transmitter_receiver_indexes, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( character( transmitter_dipole ) :: transmitter%dipole )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%dipole, transmitter_dipole, MPI_CHARACTER, child_comm, ierr )
                        !
                    class default
                        stop "unpackTransmitterBuffer: Unclassified transmitter 2"
                    !
                end select
                !
           case default
              write( *, * ) "unpackTransmitterBuffer: Unclassified transmitter 3", transmitter_derived_type
              stop
           !
        end select
        !
    end function unpackTransmitterBuffer
    !
    !
    function allocateReceiverBuffer( receiver ) result( receiver_size_bytes )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: receiver
        integer :: receiver_size_bytes
        !
        integer :: i, nbytes(4)
        !
        !
        select type( receiver )
            !
            class is( ReceiverFullImpedance_t )
                !
                call MPI_PACK_SIZE( 4, MPI_INTEGER, child_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 2, MPI_LOGICAL, child_comm, nbytes(3), ierr )
                call MPI_PACK_SIZE( len( receiver%code ), MPI_CHARACTER, child_comm, nbytes(4), ierr )
                !
            class is( ReceiverFullVerticalMagnetic_t )
                !
                call MPI_PACK_SIZE( 4, MPI_INTEGER, child_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 2, MPI_LOGICAL, child_comm, nbytes(3), ierr )
                call MPI_PACK_SIZE( len( receiver%code ), MPI_CHARACTER, child_comm, nbytes(4), ierr )
                !
            class is( ReceiverOffDiagonalImpedance_t )
                !
                call MPI_PACK_SIZE( 4, MPI_INTEGER, child_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 2, MPI_LOGICAL, child_comm, nbytes(3), ierr )
                call MPI_PACK_SIZE( len( receiver%code ), MPI_CHARACTER, child_comm, nbytes(4), ierr )
                !
            class is( ReceiverSingleField_t )
                !
                call MPI_PACK_SIZE( 4, MPI_INTEGER, child_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( 4, MPI_DOUBLE_PRECISION, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 2, MPI_LOGICAL, child_comm, nbytes(3), ierr )
                call MPI_PACK_SIZE( len( receiver%code ), MPI_CHARACTER, child_comm, nbytes(4), ierr )
                !
           class default
              stop "allocateReceiverBuffer: Unclassified receiver"
           !
        end select
        !
        receiver_size_bytes = 0
        receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lex )
        receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Ley )
        receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lez )
        receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lbx )
        receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lby )
        receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lbz )
        !
        do i = 1, size( nbytes )
            receiver_size_bytes = receiver_size_bytes + nbytes(i)
        end do
        !
    end function allocateReceiverBuffer
    !
    !
    subroutine packReceiverBuffer( receiver, index )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: receiver
        integer, intent( inout )             :: index
        !
        select type( receiver )
            !
            class is( ReceiverFullImpedance_t )
                !
                call MPI_PACK( receiver_full_impedance, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( receiver%id, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%is_complex, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%interpolation_set, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%code, len( receiver%code ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
            class is( ReceiverFullVerticalMagnetic_t )
                !
                call MPI_PACK( receiver_full_vertical_magnetic, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( receiver%id, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%is_complex, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%interpolation_set, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%code, len( receiver%code ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
            class is( ReceiverOffDiagonalImpedance_t )
                !
                call MPI_PACK( receiver_off_diagonal_impedance, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( receiver%id, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%is_complex, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%interpolation_set, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%code, len( receiver%code ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
            class is( ReceiverSingleField_t )
                !
                call MPI_PACK( receiver_single_field, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( receiver%id, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%azimuth, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%is_complex, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%interpolation_set, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( receiver%code, len( receiver%code ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
           class default
              stop "packReceiverBuffer: Unclassified receiver"
           !
        end select
        !
        call packCSparseVectorBuffer( receiver%Lex, index )
        call packCSparseVectorBuffer( receiver%Ley, index )
        call packCSparseVectorBuffer( receiver%Lez, index )
        call packCSparseVectorBuffer( receiver%Lbx, index )
        call packCSparseVectorBuffer( receiver%Lby, index )
        call packCSparseVectorBuffer( receiver%Lbz, index )
        !
    end subroutine packReceiverBuffer
    !
    !
    function unpackReceiverBuffer( index ) result ( receiver )
        implicit none
        !
        integer, intent( inout )         :: index
        !
        class( Receiver_t ), allocatable :: receiver
        !
        integer :: receiver_id, receiver_type, code_size
        !
        character(:), allocatable :: code
        real( kind=prec )         :: receiver_location(3), receiver_azymuth
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver_id, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver_type, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, code_size, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver_location(1), 3, MPI_DOUBLE_PRECISION, child_comm, ierr )
        !
        select case( receiver_derived_type )
           !
           case( receiver_full_impedance )
                !
                allocate( receiver, source = ReceiverFullImpedance_t( receiver_location, receiver_type ) )
                !
           case( receiver_full_vertical_magnetic )
                !
                allocate( receiver, source = ReceiverFullVerticalMagnetic_t( receiver_location, receiver_type ) )
                !
           case( receiver_off_diagonal_impedance )
                !
                allocate( receiver, source = ReceiverOffDiagonalImpedance_t( receiver_location, receiver_type ) )
                !
           case( receiver_single_field )
                !
                call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver_azymuth, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                !
                allocate( receiver, source = ReceiverSingleField_t( receiver_location, receiver_azymuth, receiver_type ) )
                !
           case default
              stop "unpackReceiverBuffer: Unclassified receiver"
           !
        end select
        !
        receiver%id = receiver_id
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver%is_complex, 1, MPI_LOGICAL, child_comm, ierr )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver%interpolation_set, 1, MPI_LOGICAL, child_comm, ierr )
        !
        if( allocated( receiver%code ) ) deallocate( receiver%code )
        allocate( character( code_size ) :: receiver%code )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver%code, code_size, MPI_CHARACTER, child_comm, ierr )
        !
        if( allocated( receiver%Lex ) ) deallocate( receiver%Lex )
        allocate( receiver%Lex, source = unpackCSparseVectorBuffer( index ) )
        !
        if( allocated( receiver%Ley ) ) deallocate( receiver%Ley )
        allocate( receiver%Ley, source = unpackCSparseVectorBuffer( index ) )
        !
        if( allocated( receiver%Lez ) ) deallocate( receiver%Lez )
        allocate( receiver%Lez, source = unpackCSparseVectorBuffer( index ) )
        !
        if( allocated( receiver%Lbx ) ) deallocate( receiver%Lbx )
        allocate( receiver%Lbx, source = unpackCSparseVectorBuffer( index ) )
        !
        if( allocated( receiver%Lby ) ) deallocate( receiver%Lby )
        allocate( receiver%Lby, source = unpackCSparseVectorBuffer( index ) )
        !
        if( allocated( receiver%Lbz ) ) deallocate( receiver%Lbz )
        allocate( receiver%Lbz, source = unpackCSparseVectorBuffer( index ) )
        !
    end function unpackReceiverBuffer
    !
    subroutine createDataBuffer()
        implicit none
        !
        predicted_data_buffer_size = fwd_info%data_size
        !
        if( .NOT. allocated( predicted_data_buffer ) ) then
            allocate( predicted_data_buffer( predicted_data_buffer_size ) )
        endif
        !
    end subroutine createDataBuffer
    !
    subroutine allocateDataBuffer( data_handles )
        implicit none
        !
        type( Dh_t ), dimension(:), intent( in ) :: data_handles
        !
        class( DataHandle_t ), allocatable :: data_handle
        !
        !
        integer :: i, j
        integer, allocatable, dimension(:) :: nbytes
        !
        predicted_data_buffer_size = 1
        !
        do i = 1, size( data_handles )
            !
            data_handle = getDataHandle( data_handles, i )
            !
            select type( data_handle )
                !
                class is( DataHandleMT_t )
                    !
                    if( .NOT. allocated( nbytes ) ) allocate( nbytes(4) )
                    !
                    call MPI_PACK_SIZE( 4, MPI_INTEGER, child_comm, nbytes(1), ierr )
                    call MPI_PACK_SIZE( len( data_handle%code ), MPI_CHARACTER, child_comm, nbytes(2), ierr )
                    call MPI_PACK_SIZE( len( data_handle%component ), MPI_CHARACTER, child_comm, nbytes(3), ierr )
                    call MPI_PACK_SIZE( 6, MPI_DOUBLE_PRECISION, child_comm, nbytes(4), ierr )
                    !
                class is( DataHandleCSEM_t )
                    !
                    if( .NOT. allocated( nbytes ) ) allocate( nbytes(5) )
                    !
                    call MPI_PACK_SIZE( 5, MPI_INTEGER, child_comm, nbytes(1), ierr )
                    call MPI_PACK_SIZE( len( data_handle%code ), MPI_CHARACTER, child_comm, nbytes(2), ierr )
                    call MPI_PACK_SIZE( len( data_handle%component ), MPI_CHARACTER, child_comm, nbytes(3), ierr )
                    call MPI_PACK_SIZE( len( data_handle%dipole ), MPI_CHARACTER, child_comm, nbytes(4), ierr )
                    call MPI_PACK_SIZE( 12, MPI_DOUBLE_PRECISION, child_comm, nbytes(5), ierr )
                    !
                class default
                    stop "allocateDataBuffer: Unclassified data handle"
                    !
            end select
            !
            do j = 1, size( nbytes )
                predicted_data_buffer_size = predicted_data_buffer_size + nbytes(j)
            end do
            !
        end do
        !
        if( allocated( predicted_data_buffer ) ) deallocate( predicted_data_buffer )
        allocate( predicted_data_buffer( predicted_data_buffer_size ) )
        !
    end subroutine allocateDataBuffer
    !
    !
    subroutine packDataBuffer( data_handles )
        implicit none
        !
        type( Dh_t ), dimension(:), intent( in ) :: data_handles
        !
        class( DataHandle_t ), allocatable :: data_handle
        !
        integer :: i, index
        !
        index = 1
        !
        do i = 1, size( data_handles )
            !
            data_handle = getDataHandle( data_handles, i )
            !
            select type( data_handle )
                !
                class is( DataHandleMT_t )
                    !
                    call MPI_PACK( data_mt, 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%rx_type, 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( len( data_handle%code ), 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( len( data_handle%component ), 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%code, len( data_handle%code ), MPI_CHARACTER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%component, len( data_handle%component ), MPI_CHARACTER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    !
                    call MPI_PACK( data_handle%period, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%rx_location(1), 3, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%rvalue, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%imaginary, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    !
                class is( DataHandleCSEM_t )
                    !
                    call MPI_PACK( data_csem, 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%rx_type, 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( len( data_handle%code ), 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( len( data_handle%component ), 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( len( data_handle%dipole ), 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%code, len( data_handle%code ), MPI_CHARACTER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%component, len( data_handle%component ), MPI_CHARACTER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%dipole, len( data_handle%dipole ), MPI_CHARACTER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    !!
                    call MPI_PACK( data_handle%period, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%rx_location(1), 3, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%rvalue, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%imaginary, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%tx_location(1), 3, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%azimuth, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%dip, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    call MPI_PACK( data_handle%moment, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
                    !
                class default
                    stop "allocateDataBuffer: Unclassified data handle"
                    !
            end select
            !
        end do
        !
    end subroutine packDataBuffer
    !
    ! UNPACK predicted_data_buffer TO predicted_data STRUCT
    function unpackDataBuffer() result( data_handles )
        implicit none
        !
        type( Dh_t ), allocatable, dimension(:) :: data_handles
        !
        class( DataHandle_t ), allocatable :: data_handle
        !
        integer :: i, data_handles_code, data_handles_component, data_handles_dipole, index
        !
        index = 1
        !
        do i = 1, fwd_info%n_data
            !
            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_derived_type, 1, MPI_INTEGER, child_comm, ierr )
            !
            if( allocated( data_handle ) ) deallocate( data_handle )
            !
            select case( data_derived_type )
               !
               case( data_mt )
                    !
                    allocate( DataHandleMT_t :: data_handle )
                    !
                    select type( data_handle )
                        !
                        class is( DataHandleMT_t )
                            !
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%rx_type, 1, MPI_INTEGER, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handles_code, 1, MPI_INTEGER, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handles_component, 1, MPI_INTEGER, child_comm, ierr )
                            !
                            allocate( character( data_handles_code ) :: data_handle%code )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%code, data_handles_code, MPI_CHARACTER, child_comm, ierr )
                            !
                            allocate( character( data_handles_component ) :: data_handle%component )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%component, data_handles_component, MPI_CHARACTER, child_comm, ierr )
                            !
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%period, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%rx_location(1), 3, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%rvalue, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%imaginary, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            !
                    end select
                    !
               case( data_csem )
                    !
                    allocate( DataHandleCSEM_t :: data_handle )
                    !
                    select type( data_handle )
                        !
                        class is( DataHandleCSEM_t )
                            !
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%rx_type, 1, MPI_INTEGER, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handles_code, 1, MPI_INTEGER, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handles_component, 1, MPI_INTEGER, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handles_dipole, 1, MPI_INTEGER, child_comm, ierr )
                            !
                            allocate( character( data_handles_code ) :: data_handle%code )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%code, data_handles_code, MPI_CHARACTER, child_comm, ierr )
                            !
                            allocate( character( data_handles_component ) :: data_handle%component )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%component, data_handles_component, MPI_CHARACTER, child_comm, ierr )
                            !
                            allocate( character( data_handles_dipole ) :: data_handle%dipole )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%dipole, data_handles_dipole, MPI_CHARACTER, child_comm, ierr )
                            !
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%period, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%rx_location(1), 3, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%rvalue, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%imaginary, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%tx_location(1), 3, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%azimuth, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%dip, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_handle%moment, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                            !
                    end select
                    !
               case default
                  stop "unpackDataBuffer: Unknow data case"
               !
            end select
            !
            call updateDataHandleArray( data_handles, data_handle )
            !
        end do
        !
    end function unpackDataBuffer
    !
    ! RECEIVE predicted_data FROM ANY TARGET
    function receiveData() result( data_handles )
        implicit none
        !
        type( Dh_t ), allocatable, dimension(:) :: data_handles
        !
        call createDataBuffer
        call MPI_RECV( predicted_data_buffer, predicted_data_buffer_size, MPI_PACKED, fwd_info%worker_rank, MPI_ANY_TAG, child_comm, MPI_STATUS_IGNORE, ierr )
        !
        data_handles = unpackDataBuffer()
        !
    end function receiveData
    !
    ! SEND fwd_info FROM target_id
    subroutine sendData( data_handles )
        !
        type( Dh_t ), dimension(:), intent( in ) :: data_handles
        !
        call packDataBuffer( data_handles )
        !
        call MPI_SEND( predicted_data_buffer, predicted_data_buffer_size, MPI_PACKED, master_id, tag, child_comm, ierr )
        !
    end subroutine sendData
    !
    ! ALLOCATE fwd_info_buffer
    subroutine allocateFWDInfoBuffer
        !
        integer nbytes1, nbytes2, nbytes3
        !
        call MPI_PACK_SIZE( 15, MPI_CHARACTER, child_comm, nbytes1, ierr )
        call MPI_PACK_SIZE( 4, MPI_INTEGER, child_comm, nbytes2, ierr )
        !
        fwd_info_buffer_size = ( nbytes1 + nbytes2 ) + 1
        !
        if( allocated( fwd_info_buffer ) ) deallocate( fwd_info_buffer )
        allocate( fwd_info_buffer( fwd_info_buffer_size ) )
        !
    end subroutine allocateFWDInfoBuffer
    !
    ! PACK fwd_info STRUCT TO fwd_info_buffer
    subroutine packFWDInfoBuffer
        !
        integer :: index
        !
        index = 1
        !
        call MPI_PACK( fwd_info%job_name, 15, MPI_CHARACTER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%worker_rank, 1, MPI_INTEGER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%tx_index, 1, MPI_INTEGER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%n_data, 1, MPI_INTEGER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%data_size, 1, MPI_INTEGER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        !
    end subroutine packFWDInfoBuffer
    !
    ! UNPACK fwd_info_buffer TO fwd_info STRUCT
    subroutine unpackFWDInfoBuffer
        !
        integer :: index
        !
        index = 1
        !
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%job_name, 15, MPI_CHARACTER, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%worker_rank, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%tx_index, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%n_data, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%data_size, 1, MPI_INTEGER, child_comm, ierr )
        !
    end subroutine unpackFWDInfoBuffer
    !
    ! RECEIVE fwd_info FROM ANY TARGET
    subroutine receiveFromAny()
        !
        call allocateFWDInfoBuffer
        call MPI_RECV( fwd_info_buffer, fwd_info_buffer_size, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, child_comm, MPI_STATUS_IGNORE, ierr )
        call unpackFWDInfoBuffer
        !
        write( *, * ) mpi_rank, " RECV ", fwd_info%job_name, " FROM ", fwd_info%worker_rank
        !
    end subroutine receiveFromAny
    !
    ! RECEIVE fwd_info FROM target_id
    subroutine receiveFrom( target_id )
        !
        integer, intent( in )    :: target_id
        !
        call allocateFWDInfoBuffer
        call MPI_RECV( fwd_info_buffer, fwd_info_buffer_size, MPI_PACKED, target_id, MPI_ANY_TAG, child_comm, MPI_STATUS_IGNORE, ierr )
        call unpackFWDInfoBuffer
        !
        write( *, * ) mpi_rank, " RECV ", fwd_info%job_name, " FROM ", target_id
        !
    end subroutine receiveFrom
    !
    ! SEND fwd_info FROM target_id
    subroutine sendTo( target_id )
        !
        integer, intent( in )    :: target_id
        !
        call allocateFWDInfoBuffer
        call packFWDInfoBuffer
        call MPI_SEND( fwd_info_buffer, fwd_info_buffer_size, MPI_PACKED, target_id, tag, child_comm, ierr )
        !
        write( *, * ) mpi_rank, " SEND ", fwd_info%job_name, " TO ", target_id
        !
    end subroutine sendTo
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
end module DeclarationMPI
