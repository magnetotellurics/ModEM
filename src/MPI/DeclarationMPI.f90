!
!> Module with MPI send, receive, pack, unpack routines
!
module DeclarationMPI
    !
    use CoreComponents
    !
    !include 'mpif.h'
    use mpi
    !use mpi_f08
    !
    !> MPI variables
    integer :: main_comm, mpi_rank, mpi_size, ierr
    !
    integer :: tag = 2023, master_id = 0
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
    integer :: model_derived_type
    integer, parameter :: model_cell_sg = 1
    integer, parameter :: model_cell_sg_vti = 2
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
    character( len=15 ) :: job_master = "MASTER_JOB", job_finish = "STOP_JOBS", job_done = "FINISH_JOB"
    character( len=15 ) :: job_em_solve = "JOB_EM_SOLVE", job_forward = "JOB_FORWARD", job_inversion = "JOB_INVERSION"
    character( len=15 ) :: job_jmult = "JOB_JMULT", job_jmult_t = "JOB_JMULT_T"
    character( len=15 ) :: job_basic_components = "HANDLE_FWD_COMP"
    character( len=15 ) :: job_sigma_model = "HANDLE_SIGMA"
    character( len=15 ) :: job_dsigma_model = "HANDLE_DSIGMA"
    !
    !> Structure to gather necessary MPI information for the execution of the different ModEM jobs.
    type :: JobInfo_t
        !
        SEQUENCE
        !
        character( len=15 ) :: job_name
        integer :: worker_rank, i_tx
        integer :: data_size, model_size, basic_comp_size
        integer :: inv_iter, sol_index
        real( kind=prec ) :: rms_tol, lambda
        !
    end type JobInfo_t
    !
    type( JobInfo_t ) :: job_info
    !
    !> Global models
    class( ModelParameter_t ), allocatable :: sigma, dsigma
    !
    !> MPI communication buffers and sizes
    character, dimension(:), allocatable :: job_info_buffer
    character, dimension(:), allocatable :: basic_comp_buffer
    character, dimension(:), allocatable :: grid_buffer
    character, dimension(:), allocatable :: conductivity_buffer
    character, dimension(:), allocatable :: data_buffer
    character, dimension(:), allocatable :: model_buffer
    !
    integer :: job_info_buffer_size, conductivity_buffer_size
    !
contains
    !
    !> Initialize the Message Passing Interface (MPI).
    !> Define the variables:
    !> main_comm: Main communicator.
    !> mpi_size: Number of processes used (passed by mpirun -np)
    !> mpi_rank: Process identifier index.
    !
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
            !
            write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m Minimum of two processes required!"
            write( *, * ) "        1 Master and 1 Worker."
            !
            call MPI_Finalize( ierr )
            !
            stop
            !
        endif 
        !
        !> Set mpi_rank with process id for mpi_comm_world
        call MPI_Comm_rank( main_comm, mpi_rank, ierr )
        !
    end subroutine constructorMPI
    !
    !> Allocates initial memory buffer for ForwardModelling
    !> With a preset size (for workers)
    !
    subroutine createBasicComponentsBuffer()
        implicit none
        !
        if( allocated( basic_comp_buffer ) ) deallocate( basic_comp_buffer )
        !
        allocate( basic_comp_buffer( job_info%basic_comp_size ) )
        !
        basic_comp_buffer = ""
        !
    end subroutine createBasicComponentsBuffer
    !
    !> Allocates initial memory buffer for ForwardModelling
    !> Making room for necessary information about Grid, Model and arrays of Transmitters and Receivers.
    !
    function allocateBasicComponentsBuffer() result( basic_comp_size )
        implicit none
        !
        integer :: basic_comp_size
        !
        integer :: i, last_size, nbytes(11)
        !
        basic_comp_size = 1
        !
        write( *, "(A45)" ) "Component's memory in bytes:"
        !
        call MPI_PACK_SIZE( 15, MPI_INTEGER, main_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, main_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( len( model_operator_type ), MPI_CHARACTER, main_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( len( model_method ), MPI_CHARACTER, main_comm, nbytes(4), ierr )
        call MPI_PACK_SIZE( len( source_type_mt ), MPI_CHARACTER, main_comm, nbytes(5), ierr )
        call MPI_PACK_SIZE( len( source_type_csem ), MPI_CHARACTER, main_comm, nbytes(6), ierr )
        call MPI_PACK_SIZE( len( get_1d_from ), MPI_CHARACTER, main_comm, nbytes(7), ierr )
        call MPI_PACK_SIZE( len( forward_solver_type ), MPI_CHARACTER, main_comm, nbytes(8), ierr )
        call MPI_PACK_SIZE( len( inversion_type ), MPI_CHARACTER, main_comm, nbytes(9), ierr )
        call MPI_PACK_SIZE( len( joint_type ), MPI_CHARACTER, main_comm, nbytes(10), ierr )
        call MPI_PACK_SIZE( len( e_solution_file_name ), MPI_CHARACTER, main_comm, nbytes(11), ierr )
        !
        basic_comp_size = basic_comp_size + allocateGridBuffer( main_grid, .TRUE. )
        !
        write( *, "(A45, i8)" ) "Main Grid = ", basic_comp_size
        !
        last_size = basic_comp_size
        !
        do i = 1, size( transmitters )
            !
            basic_comp_size = basic_comp_size + allocateTransmitterBuffer( getTransmitter(i) )
            !
        enddo
        !
        write( *, "(A45, i8)" ) "Transmitters Array = ", basic_comp_size - last_size
        !
        last_size = basic_comp_size
        !
        do i = 1, size( receivers )
            !
            basic_comp_size = basic_comp_size + allocateReceiverBuffer( getReceiver(i) )
            !
        enddo
        !
        write( *, "(A45, i8)" ) "Receivers Array = ", basic_comp_size - last_size
        !
        do i = 1, size( nbytes )
             basic_comp_size = basic_comp_size + nbytes(i)
        enddo
        !
        write( *, "(A45, i8)" ) "Total = ", basic_comp_size
        !
        if( allocated( basic_comp_buffer ) ) deallocate( basic_comp_buffer )
        !
        allocate( basic_comp_buffer( basic_comp_size ) )
        !
        basic_comp_buffer = ""
        !
    end function allocateBasicComponentsBuffer
    !
    !> Pack initial memory buffer for ForwardModelling
    !> Gathering in the same place the necessary information about Grid, Model and arrays of Transmitters and Receivers.
    subroutine packBasicComponentsBuffer
        implicit none
        !
        integer :: i, index
        !
        index = 1
        !
        ! 15 Integers
        call MPI_PACK( model_n_air_layer, 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( max_solver_iters, 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( max_divcor_calls, 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( max_divcor_iters, 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( len( model_operator_type ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( len( model_method ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( len( source_type_mt ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( len( source_type_csem ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( len( get_1d_from ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( len( forward_solver_type ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( len( inversion_type ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( len( joint_type ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( len( e_solution_file_name ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        ! Arrays size
        call MPI_PACK( size( transmitters ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( size( receivers ), 1, MPI_INTEGER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        !
        ! 3 Reals
        call MPI_PACK( model_max_height, 1, MPI_DOUBLE_PRECISION, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( tolerance_solver, 1, MPI_DOUBLE_PRECISION, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( tolerance_divcor, 1, MPI_DOUBLE_PRECISION, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        !
        ! 9 Strings
        call MPI_PACK( model_operator_type, len( model_operator_type ), MPI_CHARACTER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( model_method, len( model_method ), MPI_CHARACTER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( source_type_mt, len( source_type_mt ), MPI_CHARACTER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( source_type_csem, len( source_type_csem ), MPI_CHARACTER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( get_1d_from, len( get_1d_from ), MPI_CHARACTER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( forward_solver_type, len( forward_solver_type ), MPI_CHARACTER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( inversion_type, len( inversion_type ), MPI_CHARACTER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( joint_type, len( joint_type ), MPI_CHARACTER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        call MPI_PACK( e_solution_file_name, len( e_solution_file_name ), MPI_CHARACTER, basic_comp_buffer, job_info%basic_comp_size, index, main_comm, ierr )
        !
        call packGridBuffer( main_grid, basic_comp_buffer, job_info%basic_comp_size, index )
        !
        do i = 1, size( transmitters )
             !
             call packTransmitterBuffer( getTransmitter(i), basic_comp_buffer, job_info%basic_comp_size, index )
             !
        enddo
        !
        do i = 1, size( receivers )
             !
             call packReceiverBuffer( getReceiver(i), basic_comp_buffer, job_info%basic_comp_size, index )
             !
        enddo
        !
    end subroutine packBasicComponentsBuffer
    !
    !> Unpack initial memory buffer for ForwardModelling how were they packaged
    !> Instantiating Grid, Model and arrays of Transmitters and Receivers.
    subroutine unpackBasicComponentsBuffer()
        implicit none
        !
        integer :: i, index, tx_id, n_transmitters, n_receivers
        integer :: n_model_operator_type, n_model_method, n_source_type_mt, n_source_type_csem
        integer :: n_get_1d_from, n_forward_solver_type, n_inversion_type
        integer :: n_joint_type, n_e_solution_file_name
        !
        index = 1
        !
        ! 15 Integers
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, model_n_air_layer, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, max_solver_iters, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, max_divcor_calls, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, max_divcor_iters, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_model_operator_type, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_model_method, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_source_type_mt, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_source_type_csem, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_get_1d_from, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_forward_solver_type, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_inversion_type, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_joint_type, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_e_solution_file_name, 1, MPI_INTEGER, main_comm, ierr )
        !
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_transmitters, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, n_receivers, 1, MPI_INTEGER, main_comm, ierr )
        !
        ! 3 Reals
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, model_max_height, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, tolerance_solver, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, tolerance_divcor, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
        !
        ! 9 Strings
        allocate( character( n_model_operator_type ) :: model_operator_type )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, model_operator_type, n_model_operator_type, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_model_method ) :: model_method )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, model_method, n_model_method, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_source_type_mt ) :: source_type_mt )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, source_type_mt, n_source_type_mt, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_source_type_csem ) :: source_type_csem )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, source_type_csem, n_source_type_csem, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_get_1d_from ) :: get_1d_from )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, get_1d_from, n_get_1d_from, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_forward_solver_type ) :: forward_solver_type )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, forward_solver_type, n_forward_solver_type, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_inversion_type ) :: inversion_type )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, inversion_type, n_inversion_type, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_joint_type ) :: joint_type )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, joint_type, n_joint_type, MPI_CHARACTER, main_comm, ierr )
        !
        allocate( character( n_e_solution_file_name ) :: e_solution_file_name )
        call MPI_UNPACK( basic_comp_buffer, job_info%basic_comp_size, index, e_solution_file_name, n_e_solution_file_name, MPI_CHARACTER, main_comm, ierr )
        !
        ! Grid
        call unpackGridBuffer( main_grid, basic_comp_buffer, job_info%basic_comp_size, index )
        !
        do i = 1, n_transmitters
            !
            ierr = updateTransmitterArray( unpackTransmitterBuffer( basic_comp_buffer, job_info%basic_comp_size, index ) )
            !
        enddo
        !
        do i = 1, n_receivers
            !
            ierr = updateReceiverArray( unpackReceiverBuffer( basic_comp_buffer, job_info%basic_comp_size, index ) )
            !
        enddo
        !
    end subroutine unpackBasicComponentsBuffer
    !
    !> Receive grid from any target
    subroutine receiveBasicComponents( target_id )
        implicit none
        !
        integer, intent( in ) :: target_id
        !
        call createBasicComponentsBuffer
        !
        call MPI_RECV( basic_comp_buffer, job_info%basic_comp_size, MPI_PACKED, target_id, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        !
        call unpackBasicComponentsBuffer
        !
        deallocate( basic_comp_buffer )
        !
    end subroutine receiveBasicComponents
    !
    !> Send basic components to target_id
    subroutine sendBasicComponents( target_id )
        implicit none
        !
        integer, intent( in ) :: target_id
        !
        call packBasicComponentsBuffer
        !
        call MPI_SEND( basic_comp_buffer, job_info%basic_comp_size, MPI_PACKED, target_id, tag, main_comm, ierr )
        !
    end subroutine sendBasicComponents
    !
    !> No subroutine briefing
    !
    function allocateGridBuffer( grid, is_embedded ) result( grid_buffer_size )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        logical, intent( in ) :: is_embedded
        !
        integer :: grid_buffer_size
        !
        integer :: i, nbytes(24)
        !
        if( is_embedded ) then
            !
            grid_buffer_size = 0
            !
        else
            !
            grid_buffer_size = 1
            !
        endif
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
        call MPI_PACK_SIZE( size( grid%del_x ), MPI_DOUBLE_PRECISION, main_comm, nbytes(11), ierr )
        call MPI_PACK_SIZE( size( grid%del_y ), MPI_DOUBLE_PRECISION, main_comm, nbytes(12), ierr )
        call MPI_PACK_SIZE( size( grid%del_z ), MPI_DOUBLE_PRECISION, main_comm, nbytes(13), ierr )
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
            grid_buffer_size = grid_buffer_size + nbytes(i)
        enddo
        !
        if( .NOT. is_embedded ) then
            !
            if( allocated( grid_buffer ) ) deallocate( grid_buffer )
            !
            allocate( grid_buffer( grid_buffer_size ) )
            !
            grid_buffer = ""
            !
        endif
        !
    end function allocateGridBuffer
    !
    !> No subroutine briefing
    subroutine packGridBuffer( grid, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        character, dimension(:), allocatable, intent( inout ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        select type( grid )
            !
            class is( Grid3D_SG_t )
                !
                call MPI_PACK( grid_3d_sg, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                ! SIZES FOR THE FUTURE
                call MPI_PACK( size( grid%dx ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%dy ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%dz ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%dxInv ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%dyInv ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%dzInv ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%del_x ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%del_y ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%del_z ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%delXInv ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%delYInv ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%delZInv ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%xEdge ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%yEdge ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%zEdge ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%xCenter ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%yCenter ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( grid%zCenter ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( grid%geometry, 80, MPI_CHARACTER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%ox, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%oy, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%oz, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%rotDeg, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%nx, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%ny, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%nz, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%nzAir, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%nzEarth, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( grid%dx(1), size( grid%dx ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%dy(1), size( grid%dy ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%dz(1), size( grid%dz ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%dxInv(1), size( grid%dxInv ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%dyInv(1), size( grid%dyInv ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%dzInv(1), size( grid%dzInv ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%del_x(1), size( grid%del_x ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%del_y(1), size( grid%del_y ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%del_z(1), size( grid%del_z ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%delXInv(1), size( grid%delXInv ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%delYInv(1), size( grid%delYInv ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%delZInv(1), size( grid%delZInv ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%xEdge(1), size( grid%xEdge ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%yEdge(1), size( grid%yEdge ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%zEdge(1), size( grid%zEdge ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%xCenter(1), size( grid%xCenter ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%yCenter(1), size( grid%yCenter ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%zCenter(1), size( grid%zCenter ), MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( grid%zAirThick, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( grid%is_allocated, 1, MPI_LOGICAL, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
            class default
               stop "packGridBuffer: Unclassified grid"
            !
        end select
        !
    end subroutine packGridBuffer
    !
    !> No subroutine briefing
    !
    subroutine unpackGridBuffer( grid, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        class( Grid_t ), allocatable, intent( inout ) :: grid
        character, dimension(:), allocatable, intent( in ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        integer :: grid_dx, grid_dy, grid_dz, grid_dxInv, grid_dyInv, grid_dzInv, &
                   grid_delX, grid_delY, grid_delZ, grid_delXInv, grid_delYInv, grid_delZInv, &
                   grid_xEdge, grid_yEdge, grid_zEdge, grid_xCenter, grid_yCenter, grid_zCenter
        !
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        !
        select case( grid_derived_type )
            !
            case( grid_3d_sg )
                !
                if( allocated( grid ) ) deallocate( grid )
                allocate( Grid3D_SG_t :: grid )
                !
                select type( grid )
                   !
                   class is( Grid3D_SG_t )
                        !
                        ! SIZES NOW
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_dx, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_dy, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_dz, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_dxInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_dyInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_dzInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_delX, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_delY, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_delZ, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_delXInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_delYInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_delZInv, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_xEdge, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_yEdge, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_zEdge, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_xCenter, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_yCenter, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_zCenter, 1, MPI_INTEGER, main_comm, ierr )
                        !
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%geometry, 80, MPI_CHARACTER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%ox, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%oy, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%oz, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%rotDeg, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%nx, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%ny, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%nz, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%nzAir, 1, MPI_INTEGER, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%nzEarth, 1, MPI_INTEGER, main_comm, ierr )
                        !
                        allocate( grid%dx( grid_dx ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%dx(1), grid_dx, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%dy( grid_dy ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%dy(1), grid_dy, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%dz( grid_dz ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%dz(1), grid_dz, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%dxInv( grid_dxInv ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%dxInv(1), grid_dxInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%dyInv( grid_dyInv ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%dyInv(1), grid_dyInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%dzInv( grid_dzInv ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%dzInv(1), grid_dzInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%del_x( grid_delX ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%del_x(1), grid_delX, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%del_y( grid_delY ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%del_y(1), grid_delY, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%del_z( grid_delZ ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%del_z(1), grid_delZ, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%delXInv( grid_delXInv ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%delXInv(1), grid_delXInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%delYInv( grid_delYInv ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%delYInv(1), grid_delYInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%delZInv( grid_delZInv ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%delZInv(1), grid_delZInv, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%xEdge( grid_xEdge ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%xEdge(1), grid_xEdge, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%yEdge( grid_yEdge ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%yEdge(1), grid_yEdge, MPI_DOUBLE_PRECISION,main_comm, ierr )
                        !
                        allocate( grid%zEdge( grid_zEdge ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%zEdge(1), grid_zEdge, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%xCenter( grid_xCenter ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%xCenter(1), grid_xCenter, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%yCenter( grid_yCenter ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%yCenter(1), grid_yCenter, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        allocate( grid%zCenter( grid_zCenter ) )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%zCenter(1), grid_zCenter, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        !
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%zAirThick, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
                        !
                   class default
                      stop "unpackGridBuffer: Unclassified grid"
                   !
                end select
                !
            case default
               stop "unpackGridBuffer: Grid unknown case"
            !
        end select
        !
    end subroutine unpackGridBuffer
    !
    !> Allocate the model buffer with the predefined size in job_info
    subroutine createModelBuffer()
        implicit none
        !
        if( allocated( model_buffer ) ) deallocate( model_buffer )
        !
        allocate( model_buffer( job_info%model_size ) )
        !
        model_buffer = ""
        !
    end subroutine createModelBuffer
    !
    !> ????
    function allocateModelBuffer( model, is_embedded ) result( model_buffer_size )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: model
        logical, intent( in ) :: is_embedded
        !
        integer :: model_buffer_size
        !
        integer :: i, nbytes(4)
        !
        if( is_embedded ) then
            !
            model_buffer_size = 0
            !
        else
            !
            model_buffer_size = 1
            !
        endif
        !
        call MPI_PACK_SIZE( 10, MPI_INTEGER, main_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( len( model%param_type ), MPI_CHARACTER, main_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, main_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, nbytes(4), ierr )
        !
        select type( model )
            !
            class is( ModelParameterCell_SG_t )
                !
                model_buffer_size = model_buffer_size + allocateGridBuffer( model%param_grid, .TRUE. )
                !
                model_buffer_size = model_buffer_size + allocateScalarBuffer( model%cell_cond )
                !
            class is( ModelParameterCell_SG_VTI_t )
                !
                model_buffer_size = model_buffer_size + allocateGridBuffer( model%param_grid, .TRUE. )
                !
                model_buffer_size = model_buffer_size + allocateScalarBuffer( model%cell_cond_h )
                !
                model_buffer_size = model_buffer_size + allocateScalarBuffer( model%cell_cond_v )
                !
            class default
               stop "allocateModelBuffer: Unclassified model"
            !
        end select
        !
        do i = 1, size( nbytes )
            model_buffer_size = model_buffer_size + nbytes(i)
        enddo
        !
        if( .NOT. is_embedded ) then
            !
            if( allocated( model_buffer ) ) deallocate( model_buffer )
            !
            allocate( model_buffer( model_buffer_size ) )
            !
            model_buffer = ""
            !
        endif
        !
    end function allocateModelBuffer
    !
    !
    subroutine packModelBuffer( model, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: model
        character, dimension(:), allocatable, intent( inout ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        select type( model )
            !
            class is( ModelParameterCell_SG_t )
                !
                call MPI_PACK( model_cell_sg, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				!
				call MPI_PACK( len( model%param_type ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				call MPI_PACK( model%mKey(1), 8, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				call MPI_PACK( model%param_type, len( model%param_type ), MPI_CHARACTER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				call MPI_PACK( model%air_cond, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				call MPI_PACK( model%is_allocated, 1, MPI_LOGICAL, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				!
                call packGridBuffer( model%param_grid, parent_buffer, parent_buffer_size, index )
                !
                call packScalarBuffer( model%cell_cond, parent_buffer, parent_buffer_size, index )
                !
            class is( ModelParameterCell_SG_VTI_t )
                !
                call MPI_PACK( model_cell_sg_vti, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				!
				call MPI_PACK( len( model%param_type ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				call MPI_PACK( model%mKey(1), 8, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				call MPI_PACK( model%param_type, len( model%param_type ), MPI_CHARACTER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				call MPI_PACK( model%air_cond, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				call MPI_PACK( model%is_allocated, 1, MPI_LOGICAL, parent_buffer, parent_buffer_size, index, main_comm, ierr )
				!
                call packGridBuffer( model%param_grid, parent_buffer, parent_buffer_size, index )
                !
                call packScalarBuffer( model%cell_cond_h, parent_buffer, parent_buffer_size, index )
                !
                call packScalarBuffer( model%cell_cond_v, parent_buffer, parent_buffer_size, index )
                !
            class default
               stop "packModelBuffer: Unclassified model"
            !
        end select
        !
    end subroutine packModelBuffer
    !
    !
    subroutine unpackModelBuffer( model, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        class( ModelParameter_t ), allocatable, intent( inout ) :: model
        character, dimension(:), allocatable, intent( in ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        integer :: param_type_size
        !
        param_type_size = 0
        !
        if( allocated( model ) ) deallocate( model )
        !
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, model_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        !
        select case( model_derived_type )
            !
            case( model_cell_sg )
                !
                if( allocated( model ) ) deallocate( model )
                allocate( ModelParameterCell_SG_t :: model )
                !
                select type( model )
                    !
                    class is( ModelParameterCell_SG_t )
						!
						call MPI_UNPACK( parent_buffer, parent_buffer_size, index, param_type_size, 1, MPI_INTEGER, main_comm, ierr )
						call MPI_UNPACK( parent_buffer, parent_buffer_size, index, model%mKey(1), 8, MPI_INTEGER, main_comm, ierr )
						!
						allocate( character( param_type_size ) :: model%param_type )
						call MPI_UNPACK( parent_buffer, parent_buffer_size, index, model%param_type, param_type_size, MPI_CHARACTER, main_comm, ierr )
						!
						call MPI_UNPACK( parent_buffer, parent_buffer_size, index, model%air_cond, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
						call MPI_UNPACK( parent_buffer, parent_buffer_size, index, model%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
						!
                        call unpackGridBuffer( model%param_grid, parent_buffer, parent_buffer_size, index )
                        !
                        call unpackScalarBuffer( model%cell_cond, main_grid, parent_buffer, parent_buffer_size, index )
                        !
                        call model%setsigMap( model%param_type )
                        !
                end select
                !
            case( model_cell_sg_vti )
                !
                if( allocated( model ) ) deallocate( model )
                allocate( ModelParameterCell_SG_VTI_t :: model )
                !
                select type( model )
                    !
                    class is( ModelParameterCell_SG_VTI_t )
						!
						call MPI_UNPACK( parent_buffer, parent_buffer_size, index, param_type_size, 1, MPI_INTEGER, main_comm, ierr )
						call MPI_UNPACK( parent_buffer, parent_buffer_size, index, model%mKey(1), 8, MPI_INTEGER, main_comm, ierr )
						!
						allocate( character( param_type_size ) :: model%param_type )
						call MPI_UNPACK( parent_buffer, parent_buffer_size, index, model%param_type, param_type_size, MPI_CHARACTER, main_comm, ierr )
						!
						call MPI_UNPACK( parent_buffer, parent_buffer_size, index, model%air_cond, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
						call MPI_UNPACK( parent_buffer, parent_buffer_size, index, model%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
						!
                        call unpackGridBuffer( model%param_grid, parent_buffer, parent_buffer_size, index )
                        !
                        call unpackScalarBuffer( model%cell_cond_h, main_grid, parent_buffer, parent_buffer_size, index )
                        !
                        call unpackScalarBuffer( model%cell_cond_v, main_grid, parent_buffer, parent_buffer_size, index )
                        !
                        call model%setsigMap( model%param_type )
                        !
                end select
                !
            case default
               stop "unpackModelBuffer: Unclassified model"
            !
        end select
        !
    end subroutine unpackModelBuffer
    !
    !> Receive model from any target
    subroutine receiveModel( model, target_id )
        implicit none
        !
        class( ModelParameter_t ), allocatable, intent( inout ) :: model
        integer, intent( in ) :: target_id
        !
        integer :: index
        !
        index = 1
        !
        call createModelBuffer
        !
        call MPI_RECV( model_buffer, job_info%model_size, MPI_PACKED, target_id, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        !
        call unpackModelBuffer( model, model_buffer, job_info%model_size, index )
        !
        deallocate( model_buffer )
        !
    end subroutine receiveModel
    !
    !> Send model to target_id
    subroutine sendModel( model, target_id )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: model
        integer, intent( in ) :: target_id
        !
        integer :: index
        !
        index = 1
        !
        call packModelBuffer( model, model_buffer, job_info%model_size, index )
        !
        call MPI_SEND( model_buffer, job_info%model_size, MPI_PACKED, target_id, tag, main_comm, ierr )
        !
    end subroutine sendModel
    !
    !> No subroutine briefing
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
    subroutine packTransmitterBuffer( transmitter, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        class( Transmitter_t ), intent( in ) :: transmitter
        character, dimension(:), allocatable, intent( inout ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        integer :: i
        !
        select type( transmitter )
            !
            class is( TransmitterMT_t )
                !
                !> TYPE
                call MPI_PACK( transmitter_mt, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( transmitter%i_tx, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%n_pol, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%fwd_key(1), 8, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( transmitter%receiver_indexes ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%period, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%receiver_indexes(1), size( transmitter%receiver_indexes ), MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
            class is( TransmitterCSEM_t )
                !
                !> TYPE
                call MPI_PACK( transmitter_csem, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( transmitter%i_tx, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%n_pol, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%fwd_key(1), 8, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( size( transmitter%receiver_indexes ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( len( transmitter%dipole ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%period, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%location(1), 3, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%azimuth, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%dip, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%moment, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%receiver_indexes(1), size( transmitter%receiver_indexes ), MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( transmitter%dipole, len( transmitter%dipole ), MPI_CHARACTER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
            class default
               stop "allocateTransmitterBuffer: Unclassified transmitter"
            !
        end select
        !
    end subroutine packTransmitterBuffer
    !
    !> No subroutine briefing
    !
    function unpackTransmitterBuffer( parent_buffer, parent_buffer_size, index ) result( transmitter )
        implicit none
        !
        character, dimension(:), allocatable, intent( in ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        class( Transmitter_t ), allocatable :: transmitter
        !
        integer :: transmitter_receiver_indexes, transmitter_dipole
        !
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        !
        select case( transmitter_derived_type )
            !
            case( transmitter_mt )
                !
                allocate( TransmitterMT_t :: transmitter )
                !
            case( transmitter_csem )
                !
                allocate( TransmitterCSEM_t :: transmitter )
                !
            case default
               write( *, * ) "unpackTransmitterBuffer: Unknown transmitter case: ", transmitter_derived_type
               stop
            !
        end select
        !
        transmitter%i_sol = 0
        !
        select type( transmitter )
            !
            class is( TransmitterMT_t )
                !
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%i_tx, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%n_pol, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%fwd_key(1), 8, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter_receiver_indexes, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%period, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                !
                allocate( transmitter%receiver_indexes( transmitter_receiver_indexes ) )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%receiver_indexes(1), transmitter_receiver_indexes, MPI_INTEGER, main_comm, ierr )
            !
            class is( TransmitterCSEM_t )
                !
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%i_tx, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%n_pol, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%fwd_key(1), 8, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter_receiver_indexes, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter_dipole, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%period, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%location(1), 3, MPI_DOUBLE_PRECISION, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%azimuth, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%dip, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%moment, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                !
                allocate( transmitter%receiver_indexes( transmitter_receiver_indexes ) )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%receiver_indexes(1), transmitter_receiver_indexes, MPI_INTEGER, main_comm, ierr )
                !
                if( allocated( transmitter%dipole ) ) deallocate( transmitter%dipole )
                allocate( character( transmitter_dipole ) :: transmitter%dipole )
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, transmitter%dipole, transmitter_dipole, MPI_CHARACTER, main_comm, ierr )
                !
            class default
                stop "unpackTransmitterBuffer: Unclassified transmitter!"
            !
        end select
        !
    end function unpackTransmitterBuffer
    !
    !> No subroutine briefing
    !
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
    subroutine packReceiverBuffer( receiver, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: receiver
        character, dimension(:), allocatable, intent( inout ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        select type( receiver )
             !
             class is( ReceiverFullImpedance_t )
                !
                call MPI_PACK( receiver_full_impedance, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( receiver%i_rx, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call packCSparseVectorBuffer( receiver%Lex, parent_buffer, parent_buffer_size, index )
                call packCSparseVectorBuffer( receiver%Ley, parent_buffer, parent_buffer_size, index )
                call packCSparseVectorBuffer( receiver%Lbx, parent_buffer, parent_buffer_size, index )
                call packCSparseVectorBuffer( receiver%Lby, parent_buffer, parent_buffer_size, index )
                !
             class is( ReceiverFullVerticalMagnetic_t )
                !
                call MPI_PACK( receiver_full_vertical_magnetic, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( receiver%i_rx, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call packCSparseVectorBuffer( receiver%Lbx, parent_buffer, parent_buffer_size, index )
                call packCSparseVectorBuffer( receiver%Lby, parent_buffer, parent_buffer_size, index )
                call packCSparseVectorBuffer( receiver%Lbz, parent_buffer, parent_buffer_size, index )
                !
             class is( ReceiverOffDiagonalImpedance_t )
                !
                call MPI_PACK( receiver_off_diagonal_impedance, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( receiver%i_rx, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call packCSparseVectorBuffer( receiver%Lex, parent_buffer, parent_buffer_size, index )
                call packCSparseVectorBuffer( receiver%Ley, parent_buffer, parent_buffer_size, index )
                call packCSparseVectorBuffer( receiver%Lbx, parent_buffer, parent_buffer_size, index )
                call packCSparseVectorBuffer( receiver%Lby, parent_buffer, parent_buffer_size, index )
                !
             class is( ReceiverSingleField_t )
                !
                call MPI_PACK( receiver_single_field, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                call MPI_PACK( receiver%i_rx, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%rx_type, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                call MPI_PACK( receiver%azimuth, 1, MPI_DOUBLE_PRECISION, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
                if( receiver%azimuth == 1.0 ) &
                    call packCSparseVectorBuffer( receiver%Lex, parent_buffer, parent_buffer_size, index )
                if( receiver%azimuth == 2.0 ) &
                    call packCSparseVectorBuffer( receiver%Ley, parent_buffer, parent_buffer_size, index )
                if( receiver%azimuth == 3.0 ) &
                    call packCSparseVectorBuffer( receiver%Lbx, parent_buffer, parent_buffer_size, index )
                if( receiver%azimuth == 4.0 ) &
                    call packCSparseVectorBuffer( receiver%Lby, parent_buffer, parent_buffer_size, index )
                if( receiver%azimuth == 5.0 ) &
                    call packCSparseVectorBuffer( receiver%Lbz, parent_buffer, parent_buffer_size, index )
                !
            class default
               stop "packReceiverBuffer: Unclassified receiver"
            !
        end select
        !
        call MPI_PACK( receiver%is_complex, 1, MPI_LOGICAL, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( receiver%interpolation_set, 1, MPI_LOGICAL, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( receiver%code, len( receiver%code ), MPI_CHARACTER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        !
    end subroutine packReceiverBuffer
    !
    !> No subroutine briefing
    !
    function unpackReceiverBuffer( parent_buffer, parent_buffer_size, index ) result( receiver )
        implicit none
        !
        character, dimension(:), allocatable, intent( in ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        class( Receiver_t ), allocatable :: receiver
        !
        integer :: receiver_id, receiver_type, code_size
        !
        character(:), allocatable :: code
        real( kind=prec ) :: receiver_location(3), receiver_azymuth
        !
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, receiver_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, receiver_id, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, receiver_type, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, code_size, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, receiver_location(1), 3, MPI_DOUBLE_PRECISION, main_comm, ierr )
        !
        select case( receiver_derived_type )
            !
            case( receiver_full_impedance )
                !
                allocate( receiver, source = ReceiverFullImpedance_t( receiver_location, receiver_type ) )
                !
                call unpackCSparseVectorBuffer( receiver%Lex, main_grid, parent_buffer, parent_buffer_size, index )
                call unpackCSparseVectorBuffer( receiver%Ley, main_grid, parent_buffer, parent_buffer_size, index )
                call unpackCSparseVectorBuffer( receiver%Lbx, main_grid, parent_buffer, parent_buffer_size, index )
                call unpackCSparseVectorBuffer( receiver%Lby, main_grid, parent_buffer, parent_buffer_size, index )
                !
            case( receiver_full_vertical_magnetic )
                !
                allocate( receiver, source = ReceiverFullVerticalMagnetic_t( receiver_location, receiver_type ) )
                !
                call unpackCSparseVectorBuffer( receiver%Lbx, main_grid, parent_buffer, parent_buffer_size, index )
                call unpackCSparseVectorBuffer( receiver%Lby, main_grid, parent_buffer, parent_buffer_size, index )
                call unpackCSparseVectorBuffer( receiver%Lbz, main_grid, parent_buffer, parent_buffer_size, index )
                !
            case( receiver_off_diagonal_impedance )
                !
                allocate( receiver, source = ReceiverOffDiagonalImpedance_t( receiver_location, receiver_type ) )
                !
                call unpackCSparseVectorBuffer( receiver%Lex, main_grid, parent_buffer, parent_buffer_size, index )
                call unpackCSparseVectorBuffer( receiver%Ley, main_grid, parent_buffer, parent_buffer_size, index )
                call unpackCSparseVectorBuffer( receiver%Lbx, main_grid, parent_buffer, parent_buffer_size, index )
                call unpackCSparseVectorBuffer( receiver%Lby, main_grid, parent_buffer, parent_buffer_size, index )
                !
            case( receiver_single_field )
                !
                call MPI_UNPACK( parent_buffer, parent_buffer_size, index, receiver_azymuth, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
                !
                allocate( receiver, source = ReceiverSingleField_t( receiver_location, receiver_azymuth, receiver_type ) )
                !
                if( receiver_azymuth == 1.0 ) &
                    call unpackCSparseVectorBuffer( receiver%Lex, main_grid, parent_buffer, parent_buffer_size, index )
                if( receiver_azymuth == 2.0 ) &
                    call unpackCSparseVectorBuffer( receiver%Ley, main_grid, parent_buffer, parent_buffer_size, index )
                if( receiver_azymuth == 3.0 ) &
                    call unpackCSparseVectorBuffer( receiver%Lbx, main_grid, parent_buffer, parent_buffer_size, index )
                if( receiver_azymuth == 4.0 ) &
                    call unpackCSparseVectorBuffer( receiver%Lby, main_grid, parent_buffer, parent_buffer_size, index )
                if( receiver_azymuth == 5.0 ) &
                    call unpackCSparseVectorBuffer( receiver%Lbz, main_grid, parent_buffer, parent_buffer_size, index )
                !
            case default
               stop "unpackReceiverBuffer: Unclassified receiver"
            !
        end select
        !
        receiver%i_rx = receiver_id
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, receiver%is_complex, 1, MPI_LOGICAL, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, receiver%interpolation_set, 1, MPI_LOGICAL, main_comm, ierr )
        !
        if( allocated( receiver%code ) ) deallocate( receiver%code )
        allocate( character( code_size ) :: receiver%code )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, receiver%code, code_size, MPI_CHARACTER, main_comm, ierr )
        !
    end function unpackReceiverBuffer
    !
    !> ????
    subroutine createDataBuffer()
        implicit none
        !
        if( allocated( data_buffer ) ) deallocate( data_buffer )
        !
        allocate( data_buffer( job_info%data_size ) )
        !
        data_buffer = ""
        !
    end subroutine createDataBuffer
    !
    !>
    function allocateDataBuffer( tx_data ) result( data_buffer_size )
        implicit none
        !
        type( DataGroupTx_t ), intent( in ) :: tx_data
        !
        integer :: data_buffer_size
        !
        integer :: i, j, int_byte, nbytes(5)
        !
        call MPI_PACK_SIZE( 2, MPI_INTEGER, main_comm, int_byte, ierr )
        !
        data_buffer_size = int_byte + 1
        !
        do i = 1, size( tx_data%data )
             !
             call MPI_PACK_SIZE( 4, MPI_INTEGER, main_comm, nbytes(1), ierr )
             call MPI_PACK_SIZE( 3, MPI_LOGICAL, main_comm, nbytes(2), ierr )
             call MPI_PACK_SIZE( tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, main_comm, nbytes(3), ierr )
             call MPI_PACK_SIZE( tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, main_comm, nbytes(4), ierr )
             call MPI_PACK_SIZE( tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, main_comm, nbytes(5), ierr )
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
    end function allocateDataBuffer
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
        call MPI_PACK( tx_data%i_tx, 1, MPI_INTEGER, data_buffer, job_info%data_size, index, main_comm, ierr )
        call MPI_PACK( size( tx_data%data ), 1, MPI_INTEGER, data_buffer, job_info%data_size, index, main_comm, ierr )
        !
        do i = 1, size( tx_data%data )
            !
            call MPI_PACK( tx_data%data(i)%n_comp, 1, MPI_INTEGER, data_buffer, job_info%data_size, index, main_comm, ierr )
            call MPI_PACK( tx_data%data(i)%i_dg, 1, MPI_INTEGER, data_buffer, job_info%data_size, index, main_comm, ierr )
            call MPI_PACK( tx_data%data(i)%i_rx, 1, MPI_INTEGER, data_buffer, job_info%data_size, index, main_comm, ierr )
            call MPI_PACK( tx_data%data(i)%i_tx, 1, MPI_INTEGER, data_buffer, job_info%data_size, index, main_comm, ierr )
            call MPI_PACK( tx_data%data(i)%is_allocated, 1, MPI_LOGICAL, data_buffer, job_info%data_size, index, main_comm, ierr )
            call MPI_PACK( tx_data%data(i)%is_complex, 1, MPI_LOGICAL, data_buffer, job_info%data_size, index, main_comm, ierr )
            call MPI_PACK( tx_data%data(i)%error_bar, 1, MPI_LOGICAL, data_buffer, job_info%data_size, index, main_comm, ierr )
            !
            call MPI_PACK( tx_data%data(i)%reals(1), tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, data_buffer, job_info%data_size, index, main_comm, ierr )
            call MPI_PACK( tx_data%data(i)%imaginaries(1), tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, data_buffer, job_info%data_size, index, main_comm, ierr )
            call MPI_PACK( tx_data%data(i)%errors(1), tx_data%data(i)%n_comp, MPI_DOUBLE_PRECISION, data_buffer, job_info%data_size, index, main_comm, ierr )
            !
        enddo
        !
    end subroutine packDataBuffer
    !
    !> UNPACK data_buffer TO all_predicted_data STRUCT
    subroutine unpackDataBuffer( tx_data )
        implicit none
        !
        type( DataGroupTx_t ), intent( inout ) :: tx_data
        !
        type( DataGroup_t ) :: data_group
        !
        integer :: i, i_tx, n_data, n_comp, index
        !
        index = 1
        !
        call MPI_UNPACK( data_buffer, job_info%data_size, index, i_tx, 1, MPI_INTEGER, main_comm, ierr )
        !
        tx_data = DataGroupTx_t( i_tx )
        !
        call MPI_UNPACK( data_buffer, job_info%data_size, index, n_data, 1, MPI_INTEGER, main_comm, ierr )
        !
        do i = 1, n_data
            !
            call MPI_UNPACK( data_buffer, job_info%data_size, index, n_comp, 1, MPI_INTEGER, main_comm, ierr )
            !
            call MPI_UNPACK( data_buffer, job_info%data_size, index, data_group%i_dg, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( data_buffer, job_info%data_size, index, data_group%i_rx, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( data_buffer, job_info%data_size, index, data_group%i_tx, 1, MPI_INTEGER, main_comm, ierr )
            !
            call MPI_UNPACK( data_buffer, job_info%data_size, index, data_group%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
            call MPI_UNPACK( data_buffer, job_info%data_size, index, data_group%is_complex, 1, MPI_LOGICAL, main_comm, ierr )
            call MPI_UNPACK( data_buffer, job_info%data_size, index, data_group%error_bar, 1, MPI_LOGICAL, main_comm, ierr )
            !
            data_group%n_comp = n_comp
            !
            if( allocated( data_group%reals ) ) deallocate( data_group%reals )
            allocate( data_group%reals( n_comp ) )
            !
            call MPI_UNPACK( data_buffer, job_info%data_size, index, data_group%reals(1), n_comp, MPI_DOUBLE_PRECISION, main_comm, ierr )
            !
            if( allocated( data_group%imaginaries ) ) deallocate( data_group%imaginaries )
            allocate( data_group%imaginaries( n_comp ) )
            !
            call MPI_UNPACK( data_buffer, job_info%data_size, index, data_group%imaginaries(1), n_comp, MPI_DOUBLE_PRECISION, main_comm, ierr )
            !
            if( allocated( data_group%errors ) ) deallocate( data_group%errors )
            allocate( data_group%errors( n_comp ) )
            !
            call MPI_UNPACK( data_buffer, job_info%data_size, index, data_group%errors(1), n_comp, MPI_DOUBLE_PRECISION, main_comm, ierr )
            !
            call tx_data%put( data_group )
            !
        enddo
        !
    end subroutine unpackDataBuffer
    !
    !> Receive a DataGroupTx from any target
    subroutine receiveData( tx_data, target_id )
        implicit none
        !
        type( DataGroupTx_t ), intent( inout ) :: tx_data
        !
        integer, intent( in ) :: target_id
        !
        call createDataBuffer
        !
        call MPI_RECV( data_buffer, job_info%data_size, MPI_PACKED, target_id, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        !
        call unpackDataBuffer( tx_data )
        !
        deallocate( data_buffer )
        !
    end subroutine receiveData
    !
    !> Send a DataGroupTx to any target
    subroutine sendData( tx_data, target_id )
        !
        type( DataGroupTx_t ), intent( in ) :: tx_data
        !
        integer, intent( in ) :: target_id
        !
        call packDataBuffer( tx_data )
        !
        call MPI_SEND( data_buffer, job_info%data_size, MPI_PACKED, target_id, tag, main_comm, ierr )
        !
        deallocate( data_buffer )
        !
    end subroutine sendData
    !
    subroutine allocateConductivityBuffer( ccond )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: ccond
        !
        integer :: i, nbytes(4)
        !
        conductivity_buffer_size = 1
        !
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, main_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 8, MPI_INTEGER, main_comm, nbytes(3), ierr )
        !
        select type( ccond )
            !
            class is( rScalar3D_SG_t )
               !
               call MPI_PACK_SIZE( ccond%Nxyz, MPI_DOUBLE_COMPLEX, main_comm, nbytes(4), ierr )
               !
            class default
               stop "allocateConductivityBuffer: Unclassified ccond"
            !
        end select
        !
        do i = 1, size( nbytes )
            conductivity_buffer_size = conductivity_buffer_size + nbytes(i)
        enddo
        !
        if( allocated( conductivity_buffer ) ) deallocate( conductivity_buffer )
        allocate( conductivity_buffer( conductivity_buffer_size ) )
        !
        conductivity_buffer = ""
        !
    end subroutine allocateConductivityBuffer
    !
    !> No subroutine briefing
    subroutine packConductivityBuffer( ccond )
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
                call MPI_PACK( ccond%grid_type, 4, MPI_CHARACTER, conductivity_buffer, conductivity_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%is_allocated, 1, MPI_LOGICAL, conductivity_buffer, conductivity_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%nx, 1, MPI_INTEGER, conductivity_buffer, conductivity_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%ny, 1, MPI_INTEGER, conductivity_buffer, conductivity_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%nz, 1, MPI_INTEGER, conductivity_buffer, conductivity_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%NdV(1), 3, MPI_INTEGER, conductivity_buffer, conductivity_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%Nxyz, 1, MPI_INTEGER, conductivity_buffer, conductivity_buffer_size, index, main_comm, ierr )
                call MPI_PACK( ccond%store_state, 1, MPI_INTEGER, conductivity_buffer, conductivity_buffer_size, index, main_comm, ierr )
                !
                aux_array = ccond%getArray()
                !
                call MPI_PACK( aux_array(1), ccond%Nxyz, MPI_DOUBLE_COMPLEX, conductivity_buffer, conductivity_buffer_size, index, main_comm, ierr )
                !
                deallocate( aux_array )
                !
            class default
               stop "packConductivityBuffer: Unclassified ccond"
            !
        end select
        !
    end subroutine packConductivityBuffer
    !
    !> UNPACK conductivity_buffer TO all_predicted_data STRUCT
    subroutine unpackConductivityBuffer( ccond )
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
                call MPI_UNPACK( conductivity_buffer, conductivity_buffer_size, index, grid_type, 4, MPI_CHARACTER, main_comm, ierr )
                call MPI_UNPACK( conductivity_buffer, conductivity_buffer_size, index, ccond%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
                call MPI_UNPACK( conductivity_buffer, conductivity_buffer_size, index, ccond%nx, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( conductivity_buffer, conductivity_buffer_size, index, ccond%ny, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( conductivity_buffer, conductivity_buffer_size, index, ccond%nz, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( conductivity_buffer, conductivity_buffer_size, index, ccond%NdV(1), 3, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( conductivity_buffer, conductivity_buffer_size, index, ccond%Nxyz, 1, MPI_INTEGER, main_comm, ierr )
                call MPI_UNPACK( conductivity_buffer, conductivity_buffer_size, index, ccond%store_state, 1, MPI_INTEGER, main_comm, ierr )
                !
                allocate( aux_array( ccond%Nxyz ) )
                !
                call MPI_UNPACK( conductivity_buffer, conductivity_buffer_size, index, aux_array(1), ccond%Nxyz, MPI_DOUBLE_COMPLEX, main_comm, ierr )
                !
                call ccond%setArray( aux_array )
                !
                deallocate( aux_array )
                !
            class default
                stop "unpackConductivityBuffer: Unclassified ccond"
            !
        end select
        !
    end subroutine unpackConductivityBuffer
    !
    !> Receive conductivity from any target
    subroutine receiveConductivity( ccond, target_id )
        implicit none
        !
        class( Scalar_t ), allocatable, intent( inout ) :: ccond
        !
        integer, intent( in ) :: target_id
        !
        call allocateConductivityBuffer( ccond )
        !
        call MPI_RECV( conductivity_buffer, conductivity_buffer_size, MPI_PACKED, target_id, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        !
        call unpackConductivityBuffer( ccond )
        !
        deallocate( conductivity_buffer )
        !
    end subroutine receiveConductivity
    !
    !> SEND job_info FROM target_id
    subroutine sendConductivity( ccond, target_id )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: ccond
        !
        integer, intent( in ) :: target_id
        !
        call allocateConductivityBuffer( ccond )
        !
        call packConductivityBuffer( ccond )
        !
        call MPI_SEND( conductivity_buffer, conductivity_buffer_size, MPI_PACKED, target_id, tag, main_comm, ierr )
        !
        deallocate( conductivity_buffer )
        !
    end subroutine sendConductivity
    !
    !> ALLOCATE job_info_buffer
    subroutine allocateJobInfoBuffer
        !
        integer nbytes1, nbytes2, nbytes3
        !
        call MPI_PACK_SIZE( 15, MPI_CHARACTER, main_comm, nbytes1, ierr )
        call MPI_PACK_SIZE( 7, MPI_INTEGER, main_comm, nbytes2, ierr )
        call MPI_PACK_SIZE( 2, MPI_DOUBLE_PRECISION, main_comm, nbytes3, ierr )
        !
        job_info_buffer_size = ( nbytes1 + nbytes2 + nbytes3 ) + 1
        !
        if( allocated( job_info_buffer ) ) deallocate( job_info_buffer )
        allocate( job_info_buffer( job_info_buffer_size ) )
        !
        job_info_buffer = ""
        !
    end subroutine allocateJobInfoBuffer
    !
    !> PACK job_info STRUCT TO job_info_buffer
    subroutine packJobInfoBuffer
        !
        integer :: index
        !
        index = 1
        !
        call MPI_PACK( job_info%job_name, 15, MPI_CHARACTER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%worker_rank, 1, MPI_INTEGER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%i_tx, 1, MPI_INTEGER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%data_size, 1, MPI_INTEGER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%model_size, 1, MPI_INTEGER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%basic_comp_size, 1, MPI_INTEGER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%inv_iter, 1, MPI_INTEGER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%sol_index, 1, MPI_INTEGER, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%rms_tol, 1, MPI_DOUBLE_PRECISION, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        call MPI_PACK( job_info%lambda, 1, MPI_DOUBLE_PRECISION, job_info_buffer, job_info_buffer_size, index, main_comm, ierr )
        !
    end subroutine packJobInfoBuffer
    !
    !> UNPACK job_info_buffer TO job_info STRUCT
    subroutine unpackJobInfoBuffer
        !
        integer :: index
        !
        index = 1
        !
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%job_name, 15, MPI_CHARACTER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%worker_rank, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%i_tx, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%data_size, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%model_size, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%basic_comp_size, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%inv_iter, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%sol_index, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%rms_tol, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
        call MPI_UNPACK( job_info_buffer, job_info_buffer_size, index, job_info%lambda, 1, MPI_DOUBLE_PRECISION, main_comm, ierr )
        !
    end subroutine unpackJobInfoBuffer
    !
    !> RECEIVE job_info FROM ANY TARGET
    subroutine receiveFromAny()
        !
        call allocateJobInfoBuffer
        !
        call MPI_RECV( job_info_buffer, job_info_buffer_size, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        !
        call unpackJobInfoBuffer
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
        call allocateJobInfoBuffer
        call MPI_RECV( job_info_buffer, job_info_buffer_size, MPI_PACKED, target_id, MPI_ANY_TAG, main_comm, MPI_STATUS_IGNORE, ierr )
        call unpackJobInfoBuffer
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
        call allocateJobInfoBuffer
        call packJobInfoBuffer
        call MPI_SEND( job_info_buffer, job_info_buffer_size, MPI_PACKED, target_id, tag, main_comm, ierr )
        !
        !write( *, * ) mpi_rank, " SEND ", job_info%job_name, " TO ", target_id
        !
        deallocate( job_info_buffer )
        !
    end subroutine sendTo
    !
    !> Allocate the buffer for a single Scalar Field
    function allocateScalarBuffer( scalar ) result( scalar_buffer_size )
        implicit none
        !
        class( Scalar_t ), intent( in ):: scalar
        !
        integer :: scalar_buffer_size
        !
        integer :: i, nbytes(5)
        !
        scalar_buffer_size = 1
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, main_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 4, MPI_CHARACTER, main_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 9, MPI_INTEGER, main_comm, nbytes(4), ierr )
        call MPI_PACK_SIZE( scalar%Nxyz, MPI_DOUBLE_COMPLEX, main_comm, nbytes(5), ierr )
        !
        do i = 1, size( nbytes )
            scalar_buffer_size = scalar_buffer_size + nbytes(i)
        enddo
        !
    end function allocateScalarBuffer
    !
    !> Pack the info for a single Scalar Field
    subroutine packScalarBuffer( scalar, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: scalar
        character, dimension(:), allocatable, intent( inout ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        complex( kind=prec ), allocatable :: aux_array(:)
        !
        select type( scalar )
            !
            class is( cScalar3D_SG_t )
                !
                call MPI_PACK( complex_scalar, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
            class is( rScalar3D_SG_t )
                !
                call MPI_PACK( real_scalar, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
            class default
               stop "packScalarBuffer: Unclassified scalar"
            !
        end select
        !
        call MPI_PACK( scalar%grid_type, 4, MPI_CHARACTER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%is_allocated, 1, MPI_LOGICAL, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%nx, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%ny, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%nz, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%NdV(1), 3, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%Nxyz, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( scalar%store_state, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        !
        aux_array = scalar%getArray()
        !
        call MPI_PACK( aux_array(1), scalar%Nxyz, MPI_DOUBLE_COMPLEX, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        !
        deallocate( aux_array )
        !
    end subroutine packScalarBuffer
    !
    !> No subroutine briefing
    subroutine unpackScalarBuffer( scalar, grid, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        class( Scalar_t ), allocatable, intent( inout ) :: scalar
        class( Grid_t ), intent( in ) :: grid
        character, dimension(:), allocatable, intent( in ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        complex( kind=prec ), allocatable, dimension(:) :: aux_array
        !
        character( len=4 ) :: grid_type
        !
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, scalar_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        !
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_type, 4, MPI_CHARACTER, main_comm, ierr )
        !
        select type( grid )
             !
             class is( Grid3D_SG_t )
                !
                select case( scalar_derived_type )
                    !
                    case( real_scalar )
                        !
                        allocate( scalar, source = rScalar3D_SG_t( grid, grid_type ) )
                        !
                    case( complex_scalar )
                        !
                        allocate( scalar, source = cScalar3D_SG_t( grid, grid_type ) )
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
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, scalar%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, scalar%nx, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, scalar%ny, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, scalar%nz, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, scalar%NdV(1), 3, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, scalar%Nxyz, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, scalar%store_state, 1, MPI_INTEGER, main_comm, ierr )
        !
        allocate( aux_array( scalar%Nxyz ) )
        !
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, aux_array(1), scalar%Nxyz, MPI_DOUBLE_COMPLEX, main_comm, ierr )
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
        call MPI_PACK_SIZE( 17, MPI_INTEGER, main_comm, nbytes(4), ierr )
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
    subroutine packVectorBuffer( vector, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        class( Vector_t ), intent( in ) :: vector
        character, dimension(:), allocatable, intent( inout ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        complex( kind=prec ), allocatable, dimension(:) :: aux_array
        !
        select type( vector )
            !
            class is( cVector3D_SG_t )
                !
                call MPI_PACK( complex_vector, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
            class is( rVector3D_SG_t )
                !
                call MPI_PACK( real_vector, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
                !
            class default
               stop "packVectorBuffer: Unclassified vector"
            !
        end select
        !
        call MPI_PACK( vector%grid_type, 4, MPI_CHARACTER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%is_allocated, 1, MPI_LOGICAL, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%nx, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%ny, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%nz, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%NdX(1), 3, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%NdY(1), 3, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%NdZ(1), 3, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%Nxyz(1), 3, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        call MPI_PACK( vector%store_state, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        !
        aux_array = vector%getArray()
        !
        call MPI_PACK( aux_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_COMPLEX, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        !
        deallocate( aux_array )
        !
    end subroutine packVectorBuffer
    !
    !> No subroutine briefing
    subroutine unpackVectorBuffer( vector, grid, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        class( Vector_t ), intent( inout ) :: vector
        class( Grid_t ), intent( in ) :: grid
        character, dimension(:), allocatable, intent( in ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        complex( kind=prec ), allocatable, dimension(:) :: aux_array
        !
        character( len=4 ) :: grid_type
        !
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector_derived_type, 1, MPI_INTEGER, main_comm, ierr )
        !
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, grid_type, 4, MPI_CHARACTER, main_comm, ierr )
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
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector%nx, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector%ny, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector%nz, 1, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector%NdX(1), 3, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector%NdY(1), 3, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector%NdZ(1), 3, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector%Nxyz(1), 3, MPI_INTEGER, main_comm, ierr )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector%store_state, 1, MPI_INTEGER, main_comm, ierr )
        !
        allocate( aux_array( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3) ) )
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, aux_array(1), vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_COMPLEX, main_comm, ierr )
        call vector%setArray( aux_array )
        !
        deallocate( aux_array )
        !
    end subroutine unpackVectorBuffer
    !
    !> Allocate the buffer for a single cVectorSparse3D_SG
    function allocateCSparseVectorBuffer( sp_vector ) result( vector_size_bytes )
        implicit none
        !
        type( cVectorSparse3D_SG_t ), intent( in ):: sp_vector
        !
        integer :: i, nbytes(7), vector_size_bytes
        !
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, main_comm, vector_size_bytes, ierr )
        !
        if( sp_vector%is_allocated ) then
             !
             call MPI_PACK_SIZE( 4, MPI_CHARACTER, main_comm, nbytes(1), ierr )
             call MPI_PACK_SIZE( 6, MPI_INTEGER, main_comm, nbytes(2), ierr )
             call MPI_PACK_SIZE( size( sp_vector%i ), MPI_INTEGER, main_comm, nbytes(3), ierr )
             call MPI_PACK_SIZE( size( sp_vector%j ), MPI_INTEGER, main_comm, nbytes(4), ierr )
             call MPI_PACK_SIZE( size( sp_vector%k ), MPI_INTEGER, main_comm, nbytes(5), ierr )
             call MPI_PACK_SIZE( size( sp_vector%xyz ), MPI_INTEGER, main_comm, nbytes(6), ierr )
             call MPI_PACK_SIZE( size( sp_vector%c ), MPI_DOUBLE_COMPLEX, main_comm, nbytes(7), ierr )
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
    subroutine packCSparseVectorBuffer( sp_vector, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        type( cVectorSparse3D_SG_t ), intent( in ) :: sp_vector
        character, dimension(:), allocatable, intent( inout ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        call MPI_PACK( sp_vector%is_allocated, 1, MPI_LOGICAL, parent_buffer, parent_buffer_size, index, main_comm, ierr )
        !
        if( sp_vector%is_allocated ) then
             !
             call MPI_PACK( sp_vector%grid_type, 4, MPI_CHARACTER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( sp_vector%nCoeff, 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( size( sp_vector%i ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( size( sp_vector%j ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( size( sp_vector%k ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( size( sp_vector%xyz ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( size( sp_vector%c ), 1, MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( sp_vector%i(1), size( sp_vector%i ), MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( sp_vector%j(1), size( sp_vector%j ), MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( sp_vector%k(1), size( sp_vector%k ), MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( sp_vector%xyz(1), size( sp_vector%xyz ), MPI_INTEGER, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             call MPI_PACK( sp_vector%c(1), size( sp_vector%c ), MPI_DOUBLE_COMPLEX, parent_buffer, parent_buffer_size, index, main_comm, ierr )
             !
        endif
        !
    end subroutine packCSparseVectorBuffer
    !
    !> No subroutine briefing
    !
    subroutine unpackCSparseVectorBuffer( sp_vector, grid, parent_buffer, parent_buffer_size, index )
        implicit none
        !
        type( cVectorSparse3D_SG_t ), intent( inout ) :: sp_vector
        class( Grid_t ), target, intent( in ) :: grid
        character, dimension(:), allocatable, intent( in ) :: parent_buffer
        integer, intent( in ) :: parent_buffer_size
        integer, intent( inout ) :: index
        !
        integer :: vector_n_i, vector_n_j, vector_n_k, vector_n_xyz, vector_n_c
        !
        call MPI_UNPACK( parent_buffer, parent_buffer_size, index, sp_vector%is_allocated, 1, MPI_LOGICAL, main_comm, ierr )
        !
        if( sp_vector%is_allocated ) then
            !
            sp_vector%grid => grid
            !
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, sp_vector%grid_type, 4, MPI_CHARACTER, main_comm, ierr )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, sp_vector%nCoeff, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector_n_i, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector_n_j, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector_n_k, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector_n_xyz, 1, MPI_INTEGER, main_comm, ierr )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, vector_n_c, 1, MPI_INTEGER, main_comm, ierr )
            !
            allocate( sp_vector%i( vector_n_i ) )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, sp_vector%i(1), vector_n_i, MPI_INTEGER, main_comm, ierr )
            !
            allocate( sp_vector%j( vector_n_j ) )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, sp_vector%j(1), vector_n_j, MPI_INTEGER, main_comm, ierr )
            !
            allocate( sp_vector%k( vector_n_k ) )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, sp_vector%k(1), vector_n_k, MPI_INTEGER, main_comm, ierr )
            !
            allocate( sp_vector%xyz( vector_n_xyz ) )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, sp_vector%xyz(1), vector_n_xyz, MPI_INTEGER, main_comm, ierr )
            !
            allocate( sp_vector%c( vector_n_c ) )
            call MPI_UNPACK( parent_buffer, parent_buffer_size, index, sp_vector%c(1), vector_n_c, MPI_DOUBLE_COMPLEX, main_comm, ierr )
            !
        endif
        !
    end subroutine unpackCSparseVectorBuffer
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