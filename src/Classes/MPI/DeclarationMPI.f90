!
!> No module briefing
!
module DeclarationMPI
    !
    use ModEMControlFile
    !
    use InversionDCG
    use InversionNLCG
    !
    include 'mpif.h'
    !
    !> MPI variables
    integer :: main_comm, mpi_rank, mpi_size, ierr
    !
    integer :: tag = 2022, master_id = 0
    !
    !> MPI communication buffers and sizes
    character, dimension(:), allocatable :: fwd_buffer, job_info_buffer, data_buffer, model_buffer
    integer :: fwd_buffer_size, job_info_buffer_size, data_buffer_size, model_buffer_size
    !
    !> Flags for Polymorphic Objects
    integer :: scalar_derived_type
    integer, parameter :: real_scalar = 1
    integer, parameter :: complex_scalar = 2
    !
    integer :: vector_derived_type
    integer, parameter :: real_vector = 1
    integer, parameter :: complex_vector = 2
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
    !> Labels for ModEM jobs
    character( len=15 ) :: job_master = "MASTER_JOB", job_done = "FINISH_JOB", job_finish = "STOP_JOBS"
    character( len=15 ) :: job_share_memory = "SHARE_MEMORY", job_em_solve = "JOB_EM_SOLVE", job_forward = "JOB_FORWARD"
    character( len=15 ) :: job_jmult = "JOB_JMULT", job_jmult_t = "JOB_JMULT_T", job_inversion = "JOB_INVERSION"
    !
    !> Struct JobInfo_t:
    !> Gather MPI information necessary for the execution of the different ModEM jobs.
    type :: JobInfo_t
        !
        SEQUENCE
        !
        character( len=15 ) :: job_name
        integer :: worker_rank, i_tx, buffer_size
        !
    end type JobInfo_t
    !
    type( JobInfo_t ) :: job_info
    !
    !> Time counters
    real( kind=prec ) :: t_start, t_finish
    !
contains
    !
    !> Allocates initial memory buffer for ForwardModelling
    !> With a preset size (for workers)
    subroutine constructorMPI()
        implicit none
        !
        main_comm = MPI_COMM_WORLD
        !
        call MPI_Init( ierr )
        !
        !> Set mpi_size with np from mpirun for mpi_comm_world
        call MPI_Comm_size( main_comm, mpi_size, ierr )
        !
        if( mpi_size < 2 ) then
        write( *, * ) "Error: Minimum of two processes required!"
        call MPI_Finalize( ierr )
        stop
        endif 
        !
        !> Set mpi_rank with process id for mpi_comm_world
        call MPI_Comm_rank( main_comm, mpi_rank, ierr )
        !
    end subroutine constructorMPI
    !
    !> Allocates initial memory buffer for ForwardModelling
    !> With a preset size (for workers)
    subroutine createFwdBuffer()
        implicit none
        !
        fwd_buffer_size = job_info%buffer_size
        !
        if( allocated( fwd_buffer ) ) deallocate( fwd_buffer )
        !
        allocate( fwd_buffer( fwd_buffer_size ) )
        !
        fwd_buffer = ""
        !
    end subroutine createFwdBuffer
    !
    !> Allocates initial memory buffer for ForwardModelling
    !> Making room for necessary information about Grid, Model and arrays of Transmitters and Receivers.
    subroutine allocateFwdBuffer()
        implicit none
        !
        integer :: i, last_size, nbytes(10)
        !
        fwd_buffer_size = 1
        !
        write( *, "(A45)" ) "MPI FWD Allocation Sizes:"
        !
        call MPI_PACK_SIZE( 10, MPI_INTEGER, main_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, main_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( len( e_solution_file_name ), MPI_CHARACTER, main_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( len( model_method ), MPI_CHARACTER, main_comm, nbytes(4), ierr )
        call MPI_PACK_SIZE( len( forward_solver_type ), MPI_CHARACTER, main_comm, nbytes(5), ierr )
        call MPI_PACK_SIZE( len( source_type ), MPI_CHARACTER, main_comm, nbytes(6), ierr )
        call MPI_PACK_SIZE( len( get_1D_from ), MPI_CHARACTER, main_comm, nbytes(7), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, nbytes(8), ierr )
        !
        fwd_buffer_size = fwd_buffer_size + allocateGridBuffer( main_grid )
        !
        write( *, "(A45, i8)" ) "Main Grid = ", fwd_buffer_size
        last_size = fwd_buffer_size
        !
        fwd_buffer_size = fwd_buffer_size + allocateModelParameterBuffer( sigma0 )
        !
        write( *, "(A45, i8)" ) "Model Parameter = ", fwd_buffer_size - last_size
        last_size = fwd_buffer_size
        !
        if( has_pmodel_file ) then
            !
            fwd_buffer_size = fwd_buffer_size + allocateModelParameterBuffer( pmodel )
            !
            write( *, "(A45, i8)" ) "Prior Model Parameter = ", fwd_buffer_size - last_size
            last_size = fwd_buffer_size
            !
        endif
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, main_comm, nbytes(9), ierr )
        !
        do i = 1, size( transmitters )
            !
            fwd_buffer_size = fwd_buffer_size + allocateTransmitterBuffer( getTransmitter(i) )
            !
        enddo
        !
        write( *, "(A45, i8)" ) "Transmitters Array = ", fwd_buffer_size - last_size
        last_size = fwd_buffer_size
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, main_comm, nbytes(10), ierr )
        !
        do i = 1, size( receivers )
            !
            fwd_buffer_size = fwd_buffer_size + allocateReceiverBuffer( getReceiver(i) )
            !
        enddo
        !
        write( *, "(A45, i8)" ) "Receivers Array = ", fwd_buffer_size - last_size
        !
        do i = 1, size( nbytes )
             fwd_buffer_size = fwd_buffer_size + nbytes(i)
        enddo
        !
        write( *, "(A45, i8)" ) "Total = ", fwd_buffer_size
        !
        if( allocated( fwd_buffer ) ) deallocate( fwd_buffer )
        !
        allocate( fwd_buffer( fwd_buffer_size ) )
        !
        fwd_buffer = ""
        !
    end subroutine allocateFwdBuffer
    !
    !> Pack initial memory buffer for ForwardModelling
    !> Gathering in the same place the necessary information about Grid, Model and arrays of Transmitters and Receivers.
    subroutine packFwdBuffer()
        implicit none
        !
        integer :: i, index
        !
        index = 1
        !
        call MPI_PACK( QMR_iters, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( BCG_iters, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( max_divcor_calls, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( max_divcor_iters, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( len( e_solution_file_name ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( len( model_method ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( len( forward_solver_type ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( len( source_type ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( len( get_1D_from ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( model_n_air_layer, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( model_max_height, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( tolerance_divcor, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( tolerance_qmr, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( e_solution_file_name, len( e_solution_file_name ), MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( model_method, len( model_method ), MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( forward_solver_type, len( forward_solver_type ), MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( source_type, len( source_type ), MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( get_1D_from, len( get_1D_from ), MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( has_pmodel_file, 1, MPI_LOGICAL, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        !
        call packGridBuffer( main_grid, index )
        !
        call packModelParameterBuffer( sigma0, index )
        !
        if( has_pmodel_file ) then
            !
            call packModelParameterBuffer( pmodel, index )
            !
        endif
        !
        call MPI_PACK( size( transmitters ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        !
        do i = 1, size( transmitters )
             !
             call packTransmitterBuffer( getTransmitter(i), index )
             !
        enddo
        !
        call MPI_PACK( size( receivers ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        !
        do i = 1, size( receivers )
             !
             call packReceiverBuffer( getReceiver(i), index )
             !
        enddo
        !
    end subroutine packFwdBuffer
    !
    !> Unpack initial memory buffer for ForwardModelling how were they packaged
    !> Instantiating Grid, Model and arrays of Transmitters and Receivers.
    subroutine unpackFwdBuffer()
        implicit none
        !
        integer :: i, tx_id, aux_size, n_e_solution_file_name, n_model_method, n_forward_solver_type, n_source_type, n_get_1d_from, index
        !
        index = 1
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, QMR_iters, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, BCG_iters, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, max_divcor_calls, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, max_divcor_iters, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, n_e_solution_file_name, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, n_model_method, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, n_forward_solver_type, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, n_source_type, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, n_get_1d_from, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, model_n_air_layer, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, model_max_height, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, tolerance_divcor, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, tolerance_qmr, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
        !
        allocate( character( n_e_solution_file_name ) :: e_solution_file_name )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, e_solution_file_name, n_e_solution_file_name, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_model_method ) :: model_method )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, model_method, n_model_method, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_forward_solver_type ) :: forward_solver_type )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, forward_solver_type, n_forward_solver_type, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_source_type ) :: source_type )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, source_type, n_source_type, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_get_1d_from ) :: get_1D_from )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, get_1D_from, n_get_1d_from, MPI_CHARACTER, main_comm, ierr )
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, has_pmodel_file, 1, MPI_LOGICAL, main_comm, ierr )
        !
        allocate( main_grid, source = unpackGridBuffer( index ) )
        !
        allocate( sigma0, source = unpackModelParameterBuffer( index ) )
        !
        if( has_pmodel_file ) then
            !
            allocate( pmodel, source = unpackModelParameterBuffer( index ) )
            !
        endif
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, aux_size, 1, MPI_INTEGER, main_comm, ierr )
        !
        do i = 1, aux_size
            !
            tx_id = updateTransmitterArray( unpackTransmitterBuffer( index ) )
            !
        enddo
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, aux_size, 1, MPI_INTEGER, main_comm, ierr )
        !
        do i = 1, aux_size
            !
            ierr = updateReceiverArray( unpackReceiverBuffer( index ) )
            !
        enddo
        !
    end subroutine unpackFwdBuffer
    !
    !> Allocate the buffer for a single Scalar Field
    function allocateScalarBuffer( scalar ) result( scalar_size_bytes )
        implicit none
        !
        class( Scalar_t ), intent( in ):: scalar
        !
        integer :: i, nbytes(5), scalar_size_bytes
        !
        scalar_size_bytes = 0
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, main_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, main_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 7, MPI_INTEGER, main_comm, nbytes(4), ierr )
        call MPI_PACK_SIZE( scalar%Nxyz, MPI_DOUBLE_COMPLEX, main_comm, nbytes(5), ierr )
        !
        do i = 1, size( nbytes )
            scalar_size_bytes = scalar_size_bytes + nbytes(i)
        enddo
        !
    end function allocateScalarBuffer
    !
    !> Pack the info for a single Scalar Field
    subroutine packScalarBuffer( scalar, index )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: scalar
        integer, intent( inout ) :: index
        !
        complex( kind=prec ), allocatable :: aux_array(:)
        !
        select type( scalar )
            !
            class is( cScalar3D_SG_t )
                !
                call MPI_PACK( complex_scalar, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
            class is( rScalar3D_SG_t )
                !
                call MPI_PACK( real_scalar, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
            class default
               stop "packScalarBuffer: Unclassified scalar"
            !
        end select
        !
        call MPI_PACK( scalar%grid_type, 4, MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%is_allocated, 1, MPI_LOGICAL, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%nx, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%ny, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%nz, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%NdV(1), 3, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%Nxyz, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        !
        call scalar%getArray( aux_array )
        !
        call MPI_PACK( aux_array(1), scalar%Nxyz, MPI_DOUBLE_COMPLEX, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        !
        deallocate( aux_array )
        !
    end subroutine packScalarBuffer
    !
    !> No subroutine briefing
    subroutine unpackScalarBuffer( scalar, grid, index )
        implicit none
        !
        class( Scalar_t ), intent( inout ) :: scalar
        class( Grid_t ), intent( in ) :: grid
        integer, intent( inout ) :: index
        !
        complex( kind=prec ), allocatable :: aux_array(:)
        !
        character( len=4 ) :: grid_type
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, scalar_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_type, 4, MPI_CHARACTER, main_comm, ierr )
        !
        select type( grid )
             !
             class is( Grid3D_SG_t )
                !
                select case( scalar_derived_type )
                    !
                    case( real_scalar )
                        !
                        scalar = rScalar3D_SG_t( grid, grid_type )
                    !
                    case( complex_scalar )
                        !
                        scalar = cScalar3D_SG_t( grid, grid_type )
                        !
                    case default
                        stop "unpackScalarBuffer: Unknown case"
                        !
                end select
                !
             class default
                stop "unpackScalarBuffer: Unclassified grid"
                !
        end select
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, scalar%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, scalar%nx, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, scalar%ny, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, scalar%nz, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, scalar%NdV(1), 3, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, scalar%Nxyz, 1, MPI_INTEGER, main_comm, ierr )
        !
        allocate( aux_array( scalar%Nxyz ) )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, aux_array(1), scalar%Nxyz, MPI_DOUBLE_COMPLEX, main_comm, ierr )
        call scalar%setArray( aux_array )
        !
        deallocate( aux_array )
        !
    end subroutine unpackScalarBuffer
    !
    !> Allocate the buffer for a single Vector Field
    function allocateVectorBuffer( vector ) result( vector_size_bytes )
        implicit none
        !
        class( Vector_t ), intent( in ):: vector
        !
        integer :: i, nbytes(5), vector_size_bytes
        !
        vector_size_bytes = 0
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, main_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, main_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 15, MPI_INTEGER, main_comm, nbytes(4), ierr )
        !
        call MPI_PACK_SIZE( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_COMPLEX, main_comm, nbytes(5), ierr )
        !
        do i = 1, size( nbytes )
            vector_size_bytes = vector_size_bytes + nbytes(i)
        enddo
        !
    end function allocateVectorBuffer
    !
    !> No subroutine briefing
    subroutine packVectorBuffer( vector, index )
        implicit none
        !
        class( Vector_t ), intent( in ) :: vector
        integer, intent( inout ) :: index
        !
        complex( kind=prec ), allocatable :: aux_array(:)
        !
        select type( vector )
            !
            class is( cVector3D_SG_t )
                !
                call MPI_PACK( complex_vector, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
            class is( rVector3D_SG_t )
                !
                call MPI_PACK( real_vector, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
            class default
               stop "packVectorBuffer: Unclassified vector"
            !
        end select
        !
        call MPI_PACK( vector%grid_type, 4, MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%is_allocated, 1, MPI_LOGICAL, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%nx, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%ny, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%nz, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%NdX(1), 3, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%NdY(1), 3, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%NdZ(1), 3, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%Nxyz(1), 3, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        !
        call vector%getArray( aux_array )
        call MPI_PACK( aux_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_COMPLEX, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        !
        deallocate( aux_array )
        !
    end subroutine packVectorBuffer
    !
    !> No subroutine briefing
    subroutine unpackVectorBuffer( vector, grid, index )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        integer, intent( inout ) :: index
        !
        class( Vector_t ), allocatable :: vector
        !
        complex( kind=prec ), allocatable :: aux_array(:)
        !
        character( len=4 ) :: grid_type
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_type, 4, MPI_CHARACTER, main_comm, ierr )
        !
        select type( grid )
             !
             class is( Grid3D_SG_t )
                !
                select case( vector_derived_type )
                    !
                    case( real_vector )
                        !
                        vector = rVector3D_SG_t( grid, grid_type )
                    !
                    case( complex_vector )
                        !
                        vector = cVector3D_SG_t( grid, grid_type )
                        !
                    case default
                        stop "unpackVectorBuffer: Unknown case"
                        !
                end select
                !
             class default
                stop "unpackVectorBuffer: Unclassified grid"
                !
        end select
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%nx, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%ny, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%nz, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%NdX(1), 3, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%NdY(1), 3, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%NdZ(1), 3, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%Nxyz(1), 3, MPI_INTEGER, main_comm, ierr )
        !
        allocate( aux_array( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3) ) )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, aux_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_COMPLEX, main_comm, ierr )
        call vector%setArray( aux_array )
        !
        deallocate( aux_array )
        !
    end subroutine unpackVectorBuffer
    !
    !> Allocate the buffer for a single cVectorSparse3D_SG
    function allocateCSparseVectorBuffer( vector ) result( vector_size_bytes )
        implicit none
        !
        type( cVectorSparse3D_SG_t ), intent( in ):: vector
        !
        integer :: i, nbytes(7), vector_size_bytes
        !
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, vector_size_bytes, ierr )
        !
        if( vector%is_allocated ) then
             !
             call MPI_PACK_SIZE( 4, MPI_CHARACTER, main_comm, nbytes(1), ierr )
             call MPI_PACK_SIZE( 6, MPI_INTEGER, main_comm, nbytes(2), ierr )
             call MPI_PACK_SIZE( size( vector%i ), MPI_INTEGER, main_comm, nbytes(3), ierr )
             call MPI_PACK_SIZE( size( vector%j ), MPI_INTEGER, main_comm, nbytes(4), ierr )
             call MPI_PACK_SIZE( size( vector%k ), MPI_INTEGER, main_comm, nbytes(5), ierr )
             call MPI_PACK_SIZE( size( vector%xyz ), MPI_INTEGER, main_comm, nbytes(6), ierr )
             call MPI_PACK_SIZE( size( vector%c ), MPI_DOUBLE_COMPLEX, main_comm, nbytes(7), ierr )
             !
             do i = 1, size( nbytes )
                 vector_size_bytes = vector_size_bytes + nbytes(i)
             enddo
             !
        endif
        !
    end function allocateCSparseVectorBuffer
    !
    !> No subroutine briefing
    subroutine packCSparseVectorBuffer( vector, index )
        implicit none
        !
        !
        type( cVectorSparse3D_SG_t ), intent( in ) :: vector
        !
        integer, intent( inout ) :: index
        !
        call MPI_PACK( vector%is_allocated, 1, MPI_LOGICAL, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        !
        if( vector%is_allocated ) then
             !
             call MPI_PACK( vector%grid_type, 4, MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( vector%nCoeff, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( size( vector%i ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( size( vector%j ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( size( vector%k ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( size( vector%xyz ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( size( vector%c ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( vector%i(1), size( vector%i ), MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( vector%j(1), size( vector%j ), MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( vector%k(1), size( vector%k ), MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( vector%xyz(1), size( vector%xyz ), MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             call MPI_PACK( vector%c(1), size( vector%c ), MPI_DOUBLE_COMPLEX, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
             !
        endif
        !
    end subroutine packCSparseVectorBuffer
    !
    !> No function briefing
    function unpackCSparseVectorBuffer( index ) result( vector )
        implicit none
        !
        integer, intent( inout ) :: index
        type( cVectorSparse3D_SG_t ) :: vector
        !
        !
        integer :: vector_n_i, vector_n_j, vector_n_k, vector_n_xyz, vector_n_c
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
        !
        if( vector%is_allocated ) then
            !
            vector%grid => main_grid
            !
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%grid_type, 4, MPI_CHARACTER, main_comm, ierr )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%nCoeff, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector_n_i, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector_n_j, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector_n_k, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector_n_xyz, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector_n_c, 1, MPI_INTEGER, main_comm, ierr )
            !
            allocate( vector%i( vector_n_i ) )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%i(1), vector_n_i, MPI_INTEGER, main_comm, ierr )
            !
            allocate( vector%j( vector_n_j ) )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%j(1), vector_n_j, MPI_INTEGER, main_comm, ierr )
            !
            allocate( vector%k( vector_n_k ) )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%k(1), vector_n_k, MPI_INTEGER, main_comm, ierr )
            !
            allocate( vector%xyz( vector_n_xyz ) )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%xyz(1), vector_n_xyz, MPI_INTEGER, main_comm, ierr )
            !
            allocate( vector%c( vector_n_c ) )
            call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, vector%c(1), vector_n_c, MPI_DOUBLE_COMPLEX, main_comm, ierr )
            !
        endif
        !
    end function unpackCSparseVectorBuffer
    !
    !> No function briefing
    function allocateGridBuffer( grid ) result( grid_size_bytes )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        integer :: i, nbytes(24), grid_size_bytes
        !
        grid_size_bytes = 0
        !
        ! SIZES FOR THE FUTURE
        call MPI_PACK_SIZE( 19, MPI_INTEGER, main_comm, nbytes(1), ierr )
        !
        call MPI_PACK_SIZE( 80, MPI_CHARACTER, main_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 4, MPI_DOUBLE_PRECISION, main_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 5, MPI_INTEGER, main_comm, nbytes(4), ierr )
        call MPI_PACK_SIZE( size( grid%dx ), MPI_DOUBLE_PRECISION, main_comm, nbytes(5), ierr )
        call MPI_PACK_SIZE( size( grid%dy ), MPI_DOUBLE_PRECISION, main_comm, nbytes(6), ierr )
        call MPI_PACK_SIZE( size( grid%dz ), MPI_DOUBLE_PRECISION, main_comm, nbytes(7), ierr )
        call MPI_PACK_SIZE( size( grid%dxInv ), MPI_DOUBLE_PRECISION, main_comm, nbytes(8), ierr )
        call MPI_PACK_SIZE( size( grid%dyInv ), MPI_DOUBLE_PRECISION, main_comm, nbytes(9), ierr )
        call MPI_PACK_SIZE( size( grid%dzInv ), MPI_DOUBLE_PRECISION, main_comm, nbytes(10), ierr )
        call MPI_PACK_SIZE( size( grid%delX ), MPI_DOUBLE_PRECISION, main_comm, nbytes(11), ierr )
        call MPI_PACK_SIZE( size( grid%delY ), MPI_DOUBLE_PRECISION, main_comm, nbytes(12), ierr )
        call MPI_PACK_SIZE( size( grid%delZ ), MPI_DOUBLE_PRECISION, main_comm, nbytes(13), ierr )
        call MPI_PACK_SIZE( size( grid%delXInv ), MPI_DOUBLE_PRECISION, main_comm, nbytes(14), ierr )
        call MPI_PACK_SIZE( size( grid%delYInv ), MPI_DOUBLE_PRECISION, main_comm, nbytes(15), ierr )
        call MPI_PACK_SIZE( size( grid%delZInv ), MPI_DOUBLE_PRECISION, main_comm, nbytes(16), ierr )
        call MPI_PACK_SIZE( size( grid%xEdge ), MPI_DOUBLE_PRECISION, main_comm, nbytes(17), ierr )
        call MPI_PACK_SIZE( size( grid%yEdge ), MPI_DOUBLE_PRECISION, main_comm, nbytes(18), ierr )
        call MPI_PACK_SIZE( size( grid%zEdge ), MPI_DOUBLE_PRECISION, main_comm, nbytes(19), ierr )
        call MPI_PACK_SIZE( size( grid%xCenter ), MPI_DOUBLE_PRECISION, main_comm, nbytes(20), ierr )
        call MPI_PACK_SIZE( size( grid%yCenter ), MPI_DOUBLE_PRECISION, main_comm, nbytes(21), ierr )
        call MPI_PACK_SIZE( size( grid%zCenter ), MPI_DOUBLE_PRECISION, main_comm, nbytes(22), ierr )
        call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, main_comm, nbytes(23), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, nbytes(24), ierr )
        !
        do i = 1, size( nbytes )
            grid_size_bytes = grid_size_bytes + nbytes(i)
        enddo
        !
    end function allocateGridBuffer
    !
    !> No subroutine briefing
    subroutine packGridBuffer( grid, index )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        integer, intent( inout ) :: index
        !
        select type( grid )
            !
            class is( Grid3D_SG_t )
                !
                call MPI_PACK( grid_3d_sg, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                ! SIZES FOR THE FUTURE
                call MPI_PACK( size( grid%dx ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%dy ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%dz ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%dxInv ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%dyInv ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%dzInv ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%delX ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%delY ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%delZ ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%delXInv ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%delYInv ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%delZInv ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%xEdge ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%yEdge ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%zEdge ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%xCenter ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%yCenter ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%zCenter ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( grid%geometry, 80, MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%ox, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%oy, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%oz, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%rotDeg, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%nx, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%ny, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%nz, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%nzAir, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%nzEarth, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( grid%dx(1), size( grid%dx ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%dy(1), size( grid%dy ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%dz(1), size( grid%dz ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%dxInv(1), size( grid%dxInv ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%dyInv(1), size( grid%dyInv ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%dzInv(1), size( grid%dzInv ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%delX(1), size( grid%delX ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%delY(1), size( grid%delY ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%delZ(1), size( grid%delZ ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%delXInv(1), size( grid%delXInv ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%delYInv(1), size( grid%delYInv ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%delZInv(1), size( grid%delZInv ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%xEdge(1), size( grid%xEdge ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%yEdge(1), size( grid%yEdge ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%zEdge(1), size( grid%zEdge ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%xCenter(1), size( grid%xCenter ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%yCenter(1), size( grid%yCenter ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%zCenter(1), size( grid%zCenter ), MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( grid%zAirThick, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%is_allocated, 1, MPI_LOGICAL, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
            class default
               stop "packGridBuffer: Unclassified grid"
            !
        end select
        !
    end subroutine packGridBuffer
    !
    !> No function briefing
    function unpackGridBuffer( index ) result( grid )
        implicit none
        !
        integer, intent( inout ) :: index
        !
        class( Grid_t ), allocatable :: grid
        !
        integer :: grid_dx, grid_dy, grid_dz, grid_dxInv, grid_dyInv, grid_dzInv, &
                   grid_delX, grid_delY, grid_delZ, grid_delXInv, grid_delYInv, grid_delZInv, &
                   grid_xEdge, grid_yEdge, grid_zEdge, grid_xCenter, grid_yCenter, grid_zCenter
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_derived_type, 1, MPI_INTEGER, main_comm, ierr )
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
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_dx, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_dy, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_dz, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_dxInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_dyInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_dzInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_delX, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_delY, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_delZ, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_delXInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_delYInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_delZInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_xEdge, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_yEdge, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_zEdge, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_xCenter, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_yCenter, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid_zCenter, 1, MPI_INTEGER, main_comm, ierr )
                        !
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%geometry, 80, MPI_CHARACTER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%ox, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%oy, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%oz, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%rotDeg, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%nx, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%ny, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%nz, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%nzAir, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%nzEarth, 1, MPI_INTEGER, main_comm, ierr )
                        !
                        allocate( grid%dx( grid_dx ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%dx(1), grid_dx, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%dy( grid_dy ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%dy(1), grid_dy, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%dz( grid_dz ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%dz(1), grid_dz, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%dxInv( grid_dxInv ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%dxInv(1), grid_dxInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%dyInv( grid_dyInv ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%dyInv(1), grid_dyInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%dzInv( grid_dzInv ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%dzInv(1), grid_dzInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%delX( grid_delX ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%delX(1), grid_delX, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%delY( grid_delY ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%delY(1), grid_delY, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%delZ( grid_delZ ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%delZ(1), grid_delZ, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%delXInv( grid_delXInv ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%delXInv(1), grid_delXInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%delYInv( grid_delYInv ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%delYInv(1), grid_delYInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%delZInv( grid_delZInv ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%delZInv(1), grid_delZInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%xEdge( grid_xEdge ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%xEdge(1), grid_xEdge, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%yEdge( grid_yEdge ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%yEdge(1), grid_yEdge, MPI_DOUBLE_PRECISION,main_comm, ierr )
                        !
                        allocate( grid%zEdge( grid_zEdge ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%zEdge(1), grid_zEdge, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%xCenter( grid_xCenter ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%xCenter(1), grid_xCenter, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%yCenter( grid_yCenter ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%yCenter(1), grid_yCenter, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%zCenter( grid_zCenter ) )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%zCenter(1), grid_zCenter, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%zAirThick, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, grid%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
                        !
                   class default
                      stop "unpackGridBuffer: Unclassified grid"
                   !
                end select
                !
            case default
               stop "unpackGridBuffer: Grid incorrectly typed "
            !
        end select
        !
    end function unpackGridBuffer
    !
    !
    function allocateModelParameterBuffer( target_model_param ) result( model_parameter_size_bytes )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: target_model_param
        !
        integer :: i, nbytes(4), model_parameter_size_bytes
        !
        model_parameter_size_bytes = 0
        !
        select type( target_model_param )
            !
            class is( ModelParameterCell_SG_t )
                !
                call MPI_PACK_SIZE( 10, MPI_INTEGER, main_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( len( target_model_param%param_type ), MPI_CHARACTER, main_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, main_comm, nbytes(3), ierr )
                call MPI_PACK_SIZE( 2, MPI_LOGICAL, main_comm, nbytes(4), ierr )
                !
                model_parameter_size_bytes = model_parameter_size_bytes + allocateGridBuffer( target_model_param%param_grid )
                !
                model_parameter_size_bytes = model_parameter_size_bytes + allocateScalarBuffer( target_model_param%cell_cond )
                !
                do i = 1, size( nbytes )
                    model_parameter_size_bytes = model_parameter_size_bytes + nbytes(i)
                enddo
                !
            class default
               stop "allocateModelParameterBuffer: Unclassified target_model_param"
            !
        end select
        !
    end function allocateModelParameterBuffer
    !
    !
    subroutine packModelParameterBuffer( target_model_param, index )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: target_model_param
        !
        integer, intent( inout ) :: index
        !
        select type( target_model_param )
            !
            class is( ModelParameterCell_SG_t )
                !
                call MPI_PACK( model_parameter_cell_sg, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( len( target_model_param%param_type ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( target_model_param%mKey(1), 8, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( target_model_param%param_type, len( target_model_param%param_type ), MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( target_model_param%air_cond, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( target_model_param%zero_valued, 1, MPI_LOGICAL, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( target_model_param%is_allocated, 1, MPI_LOGICAL, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call packGridBuffer( target_model_param%param_grid, index )
                !
                call packScalarBuffer( target_model_param%cell_cond, index )
                !
            class default
               stop "allocateModelParameterBuffer: Unclassified target_model_param"
            !
        end select
        !
    end subroutine packModelParameterBuffer
    !
    !
    function unpackModelParameterBuffer( index ) result( target_model_param )
        implicit none
        !
        integer, intent( inout ) :: index
        !
        class( ModelParameter_t ), allocatable :: target_model_param
        !
        integer :: param_type_size
        !
        param_type_size = 0
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, model_parameter_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        !
        select case( model_parameter_derived_type )
            !
            case ( model_parameter_cell_sg )
                !
                allocate( ModelParameterCell_SG_t :: target_model_param )
                !
                select type( target_model_param )
                !
                    class is( ModelParameterCell_SG_t )
                        !
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, param_type_size, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, target_model_param%mKey(1), 8, MPI_INTEGER, main_comm, ierr )
                        !
                        allocate( character( param_type_size ) :: target_model_param%param_type )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, target_model_param%param_type, param_type_size, MPI_CHARACTER, main_comm, ierr )
                        !
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, target_model_param%air_cond, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, target_model_param%zero_valued, 1, MPI_LOGICAL, main_comm, ierr )
                        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, target_model_param%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
                        !
                        allocate( target_model_param%param_grid, source = unpackGridBuffer( index ) )
                        !
                        call unpackScalarBuffer( target_model_param%cell_cond, main_grid, index )
                        !
                        call target_model_param%SetSigMap( target_model_param%param_type )
                        !
                    class default
                        stop "unpackModelParameterBuffer: Unclassified target_model_param"
                    !
                end select
                !
            case default
               stop "unpackModelParameterBuffer: Unclassified target_model_param"
            !
        end select
        !
    end function unpackModelParameterBuffer
    !
    !> No function briefing
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
                call MPI_PACK_SIZE( 12, MPI_INTEGER, main_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, main_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( size( transmitter%receiver_indexes ), MPI_INTEGER, main_comm, nbytes(3), ierr )
                !
             class is( TransmitterCSEM_t )
                !
                allocate( nbytes(4) )
                !
                call MPI_PACK_SIZE( 13, MPI_INTEGER, main_comm, nbytes(1), ierr )
                call MPI_PACK_SIZE( 7, MPI_DOUBLE_PRECISION, main_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( size( transmitter%receiver_indexes ), MPI_INTEGER, main_comm, nbytes(3), ierr )
                call MPI_PACK_SIZE( len( transmitter%dipole ), MPI_CHARACTER, main_comm, nbytes(4), ierr )
                !
            class default
               stop "allocateTransmitterBuffer: Unclassified transmitter"
            !
        end select
        !
        do i = 1, size( nbytes )
            transmitter_size_bytes = transmitter_size_bytes + nbytes(i)
        enddo
        !
    end function allocateTransmitterBuffer
    !
    !> No subroutine briefing
    subroutine packTransmitterBuffer( transmitter, index )
        implicit none
        !
        class( Transmitter_t ), intent( in ) :: transmitter
        integer, intent( inout ) :: index
        !
        integer :: i
        !
        select type( transmitter )
            !
            class is( TransmitterMT_t )
                !
                !> TYPE
                call MPI_PACK( transmitter_mt, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( transmitter%i_tx, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%n_pol, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%fwd_key(1), 8, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( transmitter%receiver_indexes ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%period, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%receiver_indexes(1), size( transmitter%receiver_indexes ), MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
            class is( TransmitterCSEM_t )
                !
                !> TYPE
                call MPI_PACK( transmitter_csem, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( transmitter%i_tx, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%n_pol, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%fwd_key(1), 8, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( transmitter%receiver_indexes ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( len( transmitter%dipole ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%period, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%location(1), 3, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%azimuth, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%dip, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%moment, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%receiver_indexes(1), size( transmitter%receiver_indexes ), MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%dipole, len( transmitter%dipole ), MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
            class default
               stop "allocateTransmitterBuffer: Unclassified transmitter"
            !
        end select
        !
    end subroutine packTransmitterBuffer
    !
    !> No function briefing
    function unpackTransmitterBuffer( index ) result ( transmitter )
        implicit none
        !
        integer, intent( inout ) :: index
        !
        class( Transmitter_t ), allocatable :: transmitter
        !
        integer :: transmitter_receiver_indexes, transmitter_dipole
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        !
        select case( transmitter_derived_type )
            !
            case ( transmitter_mt )
                !
                allocate( TransmitterMT_t :: transmitter )
                !
            case ( transmitter_csem )
                !
                allocate( TransmitterCSEM_t :: transmitter )
                !
            case default
               write( *, * ) "unpackTransmitterBuffer: Unknown transmitter case: ", transmitter_derived_type
               stop
            !
        end select
        !
        select type( transmitter )
            !
            class is( TransmitterMT_t )
                !
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%i_tx, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%n_pol, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%fwd_key(1), 8, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter_receiver_indexes, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%period, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                !
                allocate( transmitter%receiver_indexes( transmitter_receiver_indexes ) )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%receiver_indexes(1), transmitter_receiver_indexes, MPI_INTEGER, main_comm, ierr )
            !
            class is( TransmitterCSEM_t )
                !
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%i_tx, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%n_pol, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%fwd_key(1), 8, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter_receiver_indexes, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter_dipole, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%period, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%location(1), 3, MPI_DOUBLE_PRECISION, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%azimuth, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%dip, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%moment, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                !
                allocate( transmitter%receiver_indexes( transmitter_receiver_indexes ) )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%receiver_indexes(1), transmitter_receiver_indexes, MPI_INTEGER, main_comm, ierr )
                !
                if( allocated( transmitter%dipole ) ) deallocate( transmitter%dipole )
                allocate( character( transmitter_dipole ) :: transmitter%dipole )
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, transmitter%dipole, transmitter_dipole, MPI_CHARACTER, main_comm, ierr )
                !
            class default
                stop "unpackTransmitterBuffer: Unclassified transmitter!"
            !
        end select
        !
    end function unpackTransmitterBuffer
    !
    !> No function briefing
    function allocateReceiverBuffer( receiver ) result( receiver_size_bytes )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: receiver
        integer :: receiver_size_bytes
        !
        integer :: i, nbytes(4)
        !
        receiver_size_bytes = 0
        !
        call MPI_PACK_SIZE( 4, MPI_INTEGER, main_comm, nbytes(1), ierr )
        !
        select type( receiver )
             !
             class is( ReceiverFullImpedance_t )
                !
                call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, main_comm, nbytes(2), ierr )
                !
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lex )
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Ley )
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lbx )
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lby )
                !
             class is( ReceiverFullVerticalMagnetic_t )
                !
                call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, main_comm, nbytes(2), ierr )
                !
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lbx )
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lby )
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lbz )
                !
             class is( ReceiverOffDiagonalImpedance_t )
                !
                call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, main_comm, nbytes(2), ierr )
                !
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lex )
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Ley )
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lbx )
                receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lby )
                !
             class is( ReceiverSingleField_t )
                !
                call MPI_PACK_SIZE( 4, MPI_DOUBLE_PRECISION, main_comm, nbytes(2), ierr )
                !
                if( receiver%azimuth == 1.0 ) receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lex )
                if( receiver%azimuth == 2.0 ) receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Ley )
                if( receiver%azimuth == 3.0 ) receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lbx )
                if( receiver%azimuth == 4.0 ) receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lby )
                if( receiver%azimuth == 5.0 ) receiver_size_bytes = receiver_size_bytes + allocateCSparseVectorBuffer( receiver%Lbz )
                !
            class default
               stop "allocateReceiverBuffer: Unclassified receiver"
            !
        end select
        !
        call MPI_PACK_SIZE( 2, MPI_LOGICAL, main_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( len( receiver%code ), MPI_CHARACTER, main_comm, nbytes(4), ierr )
        !
        do i = 1, size( nbytes )
             receiver_size_bytes = receiver_size_bytes + nbytes(i)
        enddo
        !
    end function allocateReceiverBuffer
    !
    !> No subroutine briefing
    subroutine packReceiverBuffer( receiver, index )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: receiver
        integer, intent( inout ) :: index
        !
        select type( receiver )
             !
             class is( ReceiverFullImpedance_t )
                !
                call MPI_PACK( receiver_full_impedance, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( receiver%i_rx, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call packCSparseVectorBuffer( receiver%Lex, index )
                call packCSparseVectorBuffer( receiver%Ley, index )
                call packCSparseVectorBuffer( receiver%Lbx, index )
                call packCSparseVectorBuffer( receiver%Lby, index )
                !
             class is( ReceiverFullVerticalMagnetic_t )
                !
                call MPI_PACK( receiver_full_vertical_magnetic, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( receiver%i_rx, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call packCSparseVectorBuffer( receiver%Lbx, index )
                call packCSparseVectorBuffer( receiver%Lby, index )
                call packCSparseVectorBuffer( receiver%Lbz, index )
                !
             class is( ReceiverOffDiagonalImpedance_t )
                !
                call MPI_PACK( receiver_off_diagonal_impedance, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( receiver%i_rx, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call packCSparseVectorBuffer( receiver%Lex, index )
                call packCSparseVectorBuffer( receiver%Ley, index )
                call packCSparseVectorBuffer( receiver%Lbx, index )
                call packCSparseVectorBuffer( receiver%Lby, index )
                !
             class is( ReceiverSingleField_t )
                !
                call MPI_PACK( receiver_single_field, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( receiver%i_rx, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%azimuth, 1, MPI_DOUBLE_PRECISION, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
                !
                if( receiver%azimuth == 1.0 ) call packCSparseVectorBuffer( receiver%Lex, index )
                if( receiver%azimuth == 2.0 ) call packCSparseVectorBuffer( receiver%Ley, index )
                if( receiver%azimuth == 3.0 ) call packCSparseVectorBuffer( receiver%Lbx, index )
                if( receiver%azimuth == 4.0 ) call packCSparseVectorBuffer( receiver%Lby, index )
                if( receiver%azimuth == 5.0 ) call packCSparseVectorBuffer( receiver%Lbz, index )
                !
            class default
               stop "packReceiverBuffer: Unclassified receiver"
            !
        end select
        !
        call MPI_PACK( receiver%is_complex, 1, MPI_LOGICAL, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( receiver%interpolation_set, 1, MPI_LOGICAL, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        call MPI_PACK( receiver%code, len( receiver%code ), MPI_CHARACTER, fwd_buffer, fwd_buffer_size, index, main_comm, ierr )
        !
    end subroutine packReceiverBuffer
    !
    !> No function briefing
    function unpackReceiverBuffer( index ) result ( receiver )
        implicit none
        !
        integer, intent( inout ) :: index
        !
        class( Receiver_t ), allocatable :: receiver
        !
        integer :: receiver_id, receiver_type, code_size
        !
        character(:), allocatable :: code
        real( kind=prec ) :: receiver_location(3), receiver_azymuth
        !
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, receiver_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, receiver_id, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, receiver_type, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, code_size, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, receiver_location(1), 3, MPI_DOUBLE_PRECISION, main_comm, ierr )
        !
        select case( receiver_derived_type )
            !
            case( receiver_full_impedance )
                !
                allocate( receiver, source = ReceiverFullImpedance_t( receiver_location, receiver_type ) )
                !
                receiver%Lex = unpackCSparseVectorBuffer( index )
                receiver%Ley = unpackCSparseVectorBuffer( index )
                receiver%Lbx = unpackCSparseVectorBuffer( index )
                receiver%Lby = unpackCSparseVectorBuffer( index )
                !
            case( receiver_full_vertical_magnetic )
                !
                allocate( receiver, source = ReceiverFullVerticalMagnetic_t( receiver_location, receiver_type ) )
                !
                receiver%Lbx = unpackCSparseVectorBuffer( index )
                receiver%Lby = unpackCSparseVectorBuffer( index )
                receiver%Lbz = unpackCSparseVectorBuffer( index )
                !
            case( receiver_off_diagonal_impedance )
                !
                allocate( receiver, source = ReceiverOffDiagonalImpedance_t( receiver_location, receiver_type ) )
                !
                receiver%Lex = unpackCSparseVectorBuffer( index )
                receiver%Ley = unpackCSparseVectorBuffer( index )
                receiver%Lbx = unpackCSparseVectorBuffer( index )
                receiver%Lby = unpackCSparseVectorBuffer( index )
                !
            case( receiver_single_field )
                !
                call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, receiver_azymuth, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                !
                allocate( receiver, source = ReceiverSingleField_t( receiver_location, receiver_azymuth, receiver_type ) )
                !
                if( receiver_azymuth == 1.0 ) receiver%Lex = unpackCSparseVectorBuffer( index )
                if( receiver_azymuth == 2.0 ) receiver%Ley = unpackCSparseVectorBuffer( index )
                if( receiver_azymuth == 3.0 ) receiver%Lbx = unpackCSparseVectorBuffer( index )
                if( receiver_azymuth == 4.0 ) receiver%Lby = unpackCSparseVectorBuffer( index )
                if( receiver_azymuth == 5.0 ) receiver%Lbz = unpackCSparseVectorBuffer( index )
                !
            case default
               stop "unpackReceiverBuffer: Unclassified receiver"
            !
        end select
        !
        receiver%i_rx = receiver_id
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, receiver%is_complex, 1, MPI_LOGICAL, main_comm, ierr )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, receiver%interpolation_set, 1, MPI_LOGICAL, main_comm, ierr )
        !
        if( allocated( receiver%code ) ) deallocate( receiver%code )
        allocate( character( code_size ) :: receiver%code )
        call MPI_UNPACK( fwd_buffer, fwd_buffer_size, index, receiver%code, code_size, MPI_CHARACTER, main_comm, ierr )
        !
    end function unpackReceiverBuffer
    !
    !> RECEIVE predicted_data FROM ANY TARGET
    subroutine receiveFwdBuffer( target_id )
        implicit none
        !
        integer, intent( in ) :: target_id
        !
        call createFwdBuffer()
        !
        call MPI_RECV( fwd_buffer, fwd_buffer_size, MPI_PACKED, target_id, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        !
        call unpackFwdBuffer()
        !
        deallocate( fwd_buffer )
        !
    end subroutine receiveFwdBuffer
    !
    !> SEND job_info FROM target_id
    subroutine sendFwdBuffer( target_id )
        implicit none
        !
        integer, intent( in ) :: target_id
        !
        call MPI_SEND( fwd_buffer, fwd_buffer_size, MPI_PACKED, target_id, tag, main_comm, ierr )
        !
    end subroutine sendFwdBuffer
    !
    !>
    subroutine allocateDataBuffer( tx_data )
        implicit none
        !
        type( DataGroupTx_t ), intent( in ) :: tx_data
        !
        integer :: i, j, int_byte, nbytes(4)
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, main_comm, int_byte, ierr )
        !
        data_buffer_size = int_byte + 1
        !
        do i = 1, size( tx_data%data )
             !
             call MPI_PACK_SIZE( 2, MPI_INTEGER, main_comm, nbytes(1), ierr )
             call MPI_PACK_SIZE( tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, main_comm, nbytes(2), ierr )
             call MPI_PACK_SIZE( tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, main_comm, nbytes(3), ierr )
             call MPI_PACK_SIZE( tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, main_comm, nbytes(4), ierr )
             !
             do j = 1, size( nbytes )
                data_buffer_size = data_buffer_size + nbytes(j)
             enddo
             !
        enddo
        !
        if( allocated( data_buffer ) ) deallocate( data_buffer )
        allocate( data_buffer( data_buffer_size ) )
        data_buffer = ""
        !
    end subroutine allocateDataBuffer
    !
    !> No subroutine briefing
    subroutine packDataBuffer( tx_data )
        implicit none
        !
        type( DataGroupTx_t ), intent( in ) :: tx_data
        !
        integer :: i, index
        !
        index = 1
        !
        call MPI_PACK( tx_data%i_tx, 1, MPI_INTEGER, data_buffer, data_buffer_size, index, main_comm, ierr )
        !
        do i = 1, size( tx_data%data )
             !
             call MPI_PACK( tx_data%data(i)%i_rx, 1, MPI_INTEGER, data_buffer, data_buffer_size, index, main_comm, ierr )
             call MPI_PACK( tx_data%data(i)%i_tx, 1, MPI_INTEGER, data_buffer, data_buffer_size, index, main_comm, ierr )
             !
             call MPI_PACK( tx_data%data(i)%reals(1), tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, data_buffer, data_buffer_size, index, main_comm, ierr )
             call MPI_PACK( tx_data%data(i)%imaginaries(1), tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, data_buffer, data_buffer_size, index, main_comm, ierr )
             call MPI_PACK( tx_data%data(i)%errors(1), tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, data_buffer, data_buffer_size, index, main_comm, ierr )
             !
        enddo
        !
    end subroutine packDataBuffer
    !
    !> UNPACK data_buffer TO predicted_data STRUCT
    subroutine unpackDataBuffer( tx_data )
        implicit none
        !
        type( DataGroupTx_t ), intent( inout ) :: tx_data
        !
        integer :: i, index
        !
        index = 1
        !
        call MPI_UNPACK( data_buffer, data_buffer_size, index, tx_data%i_tx, 1, MPI_INTEGER, main_comm, ierr )
        !
        do i = 1, size( tx_data%data )
             !
             call MPI_UNPACK( data_buffer, data_buffer_size, index, tx_data%data(i)%i_rx, 1, MPI_INTEGER, main_comm, ierr )
             call MPI_UNPACK( data_buffer, data_buffer_size, index, tx_data%data(i)%i_tx, 1, MPI_INTEGER, main_comm, ierr )
             call MPI_UNPACK( data_buffer, data_buffer_size, index, tx_data%data(i)%reals(1), tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, main_comm, ierr )
             call MPI_UNPACK( data_buffer, data_buffer_size, index, tx_data%data(i)%imaginaries(1), tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, main_comm, ierr )
             call MPI_UNPACK( data_buffer, data_buffer_size, index, tx_data%data(i)%errors(1), tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, main_comm, ierr )
             !
        enddo
        !
    end subroutine unpackDataBuffer
    !
    !> RECEIVE predicted_data FROM ANY TARGET
    subroutine receiveData( tx_data, target_id )
        implicit none
        !
        type( DataGroupTx_t ), intent( inout ) :: tx_data
        !
        integer, intent( in ) :: target_id
        !
        call allocateDataBuffer( tx_data )
        !
        call MPI_RECV( data_buffer, data_buffer_size, MPI_PACKED, target_id, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        !
        call unpackDataBuffer( tx_data )
        !
        deallocate( data_buffer )
        !
    end subroutine receiveData
    !
    !> SEND job_info FROM target_id
    subroutine sendData( tx_data, target_id )
        !
        type( DataGroupTx_t ), intent( in ) :: tx_data
        !
        integer, intent( in ) :: target_id
        !
        call allocateDataBuffer( tx_data )
        !
        call packDataBuffer( tx_data )
        !
        call MPI_SEND( data_buffer, data_buffer_size, MPI_PACKED, target_id, tag, main_comm, ierr )
        !
        deallocate( data_buffer )
        !
    end subroutine sendData
    !
    subroutine allocateModelBuffer( ccond )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: ccond
        !
        integer :: i, nbytes(4)
        !
        model_buffer_size = 1
        !
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, main_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 7, MPI_INTEGER, main_comm, nbytes(3), ierr )
        !
        select type( ccond )
            !
            class is( rScalar3D_SG_t )
               !
               call MPI_PACK_SIZE( ccond%Nxyz, MPI_DOUBLE_COMPLEX, main_comm, nbytes(4), ierr )
               !
            class default
               stop "allocateModelBuffer: Unclassified ccond"
            !
        end select
        !
        do i = 1, size( nbytes )
            model_buffer_size = model_buffer_size + nbytes(i)
        enddo
        !
        if( allocated( model_buffer ) ) deallocate( model_buffer )
        allocate( model_buffer( model_buffer_size ) )
        model_buffer = ""
        !
    end subroutine allocateModelBuffer
    !
    !> No subroutine briefing
    subroutine packModelBuffer( ccond )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: ccond
        !
        integer :: index
        !
        complex( kind=prec ), allocatable :: aux_array(:)
        !
        index = 1
        !
        select type( ccond )
            !
            class is( rScalar3D_SG_t )
                !
                call MPI_PACK( ccond%grid_type, 4, MPI_CHARACTER, model_buffer, model_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%is_allocated, 1, MPI_LOGICAL, model_buffer, model_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%nx, 1, MPI_INTEGER, model_buffer, model_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%ny, 1, MPI_INTEGER, model_buffer, model_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%nz, 1, MPI_INTEGER, model_buffer, model_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%NdV(1), 3, MPI_INTEGER, model_buffer, model_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%Nxyz, 1, MPI_INTEGER, model_buffer, model_buffer_size, index, main_comm, ierr )
                !
                call ccond%getArray( aux_array )
                call MPI_PACK( aux_array(1), ccond%Nxyz, MPI_DOUBLE_COMPLEX, model_buffer, model_buffer_size, index, main_comm, ierr )
                !
                deallocate( aux_array )
                !
            class default
               stop "packModelBuffer: Unclassified ccond"
            !
        end select
        !
    end subroutine packModelBuffer
    !
    !> UNPACK model_buffer TO predicted_data STRUCT
    subroutine unpackModelBuffer( ccond )
        implicit none
        !
        class( Scalar_t ), allocatable, intent( inout ) :: ccond
        !
        complex( kind=prec ), allocatable :: aux_array(:)
        !
        character( len=4 ) :: grid_type
        !
        integer :: index
        !
        index = 1
        !
        select type( ccond )
            !
            class is ( rScalar3D_SG_t )
                !
                call MPI_UNPACK( model_buffer, model_buffer_size, index, grid_type, 4, MPI_CHARACTER, main_comm, ierr )
                call MPI_UNPACK( model_buffer, model_buffer_size, index, ccond%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
                call MPI_UNPACK( model_buffer, model_buffer_size, index, ccond%nx, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( model_buffer, model_buffer_size, index, ccond%ny, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( model_buffer, model_buffer_size, index, ccond%nz, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( model_buffer, model_buffer_size, index, ccond%NdV(1), 3, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( model_buffer, model_buffer_size, index, ccond%Nxyz, 1, MPI_INTEGER, main_comm, ierr )
                !
                allocate( aux_array( ccond%Nxyz ) )
                call MPI_UNPACK( model_buffer, model_buffer_size, index, aux_array(1), ccond%Nxyz, MPI_DOUBLE_COMPLEX, main_comm, ierr )
                call ccond%setArray( aux_array )
                !
                deallocate( aux_array )
                !
            class default
                stop "unpackModelBuffer: Unclassified ccond"
            !
        end select
        !
    end subroutine unpackModelBuffer
    !
    !> RECEIVE predicted_data FROM ANY TARGET
    subroutine receiveModel( ccond, target_id )
        implicit none
        !
        class( Scalar_t ), allocatable, intent( inout ) :: ccond
        !
        integer, intent( in ) :: target_id
        !
        call allocateModelBuffer( ccond )
        !
        call MPI_RECV( model_buffer, model_buffer_size, MPI_PACKED, target_id, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        !
        call unpackModelBuffer( ccond )
        !
        deallocate( model_buffer )
        !
    end subroutine receiveModel
    !
    !> SEND job_info FROM target_id
    subroutine sendModel( ccond, target_id )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: ccond
        !
        integer, intent( in ) :: target_id
        !
        call allocateModelBuffer( ccond )
        !
        call packModelBuffer( ccond )
        !
        call MPI_SEND( model_buffer, model_buffer_size, MPI_PACKED, target_id, tag, main_comm, ierr )
        !
        deallocate( model_buffer )
        !
    end subroutine sendModel
    !
    !> ALLOCATE job_info_buffer
    subroutine allocateFWDInfoBuffer
        !
        integer nbytes1, nbytes2, nbytes3
        !
        call MPI_PACK_SIZE( 15, MPI_CHARACTER, main_comm, nbytes1, ierr )
        call MPI_PACK_SIZE( 3, MPI_INTEGER, main_comm, nbytes2, ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, nbytes3, ierr )
        !
        job_info_buffer_size = ( nbytes1 + nbytes2 + nbytes3 ) + 1
        !
        if( allocated( job_info_buffer ) ) deallocate( job_info_buffer )
        allocate( job_info_buffer( job_info_buffer_size ) )
        job_info_buffer = ""
        !
    end subroutine allocateFWDInfoBuffer
    !
    !> PACK job_info STRUCT TO job_info_buffer
    subroutine packFWDInfoBuffer
        !
        integer :: index
        !
        index = 1
        !
        call MPI_PACK( job_info%job_name, 15, MPI_CHARACTER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%worker_rank, 1, MPI_INTEGER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%i_tx, 1, MPI_INTEGER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%buffer_size, 1, MPI_INTEGER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        !
    end subroutine packFWDInfoBuffer
    !
    !> UNPACK job_info_buffer TO job_info STRUCT
    subroutine unpackFWDInfoBuffer
        !
        integer :: index
        !
        index = 1
        !
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%job_name, 15, MPI_CHARACTER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%worker_rank, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%i_tx, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%buffer_size, 1, MPI_INTEGER, main_comm, ierr )
        !
    end subroutine unpackFWDInfoBuffer
    !
    !> RECEIVE job_info FROM ANY TARGET
    subroutine receiveFromAny()
        !
        call allocateFWDInfoBuffer
        call MPI_RECV( job_info_buffer, job_info_buffer_size, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        call unpackFWDInfoBuffer
        !
        !write( *, * ) mpi_rank, " RECV ", job_info%job_name, " FROM ", job_info%worker_rank
        !
        deallocate( job_info_buffer )
        !
    end subroutine receiveFromAny
    !
    !> RECEIVE job_info FROM target_id
    subroutine receiveFrom( target_id )
        !
        integer, intent( in ) :: target_id
        !
        call allocateFWDInfoBuffer
        call MPI_RECV( job_info_buffer, job_info_buffer_size, MPI_PACKED, target_id, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        call unpackFWDInfoBuffer
        !
        !write( *, * ) mpi_rank, " RECV ", job_info%job_name, " FROM ", target_id
        !
        deallocate( job_info_buffer )
        !
    end subroutine receiveFrom
    !
    !> SEND job_info FROM target_id
    subroutine sendTo( target_id )
        !
        integer, intent( in ) :: target_id
        !
        call allocateFWDInfoBuffer
        call packFWDInfoBuffer
        call MPI_SEND( job_info_buffer, job_info_buffer_size, MPI_PACKED, target_id, tag, main_comm, ierr )
        !
        !write( *, * ) mpi_rank, " SEND ", job_info%job_name, " TO ", target_id
        !
        deallocate( job_info_buffer )
        !
    end subroutine sendTo
    !
    !> Return a one-dimensional array from a two-dimensional one
    function BiArrayToArray( d2_array ) result( d1_array )
        !
        real( kind=prec ), intent( in ) :: d2_array(:,:)
        real( kind=prec ), allocatable :: d1_array(:)
        !
        integer :: size_array, dim_d2_array(2)
        !
        dim_d2_array = shape( d2_array )
        size_array = product( dim_d2_array )
        !
        allocate( d1_array( size_array ) )
        d1_array = (/reshape( d1_array, (/size_array, 1/) )/)
        !
    end function BiArrayToArray
    !
    !> Return a two-dimensional array from a one-dimensional and its sizes
    function arrayToBiArray( d1_array, x, y ) result ( d2_array )
        !
        real( kind=prec ), intent( in ) :: d1_array(:)
        integer, intent( in ) :: x, y
        !
        real( kind=prec ), allocatable :: d2_array(:,:)
        !
        allocate( d2_array( x, y ) )
        !
        d2_array = reshape( d1_array, (/x, y/) )
        !
    end function arrayToBiArray
    !
end module DeclarationMPI
!