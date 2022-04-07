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
    !
    use rVector3D_SG
    use rScalar3D_SG
    use cVector3D_SG
    use cScalar3D_SG
    !
    use MetricElements_CSG
    !
    use ModelOperator_MF
    !
    use ModelParameterCell_SG
    !
    use TransmitterFArray
    use TransmitterMT
    use TransmitterCSEM
    !
    use ReceiverFArray
    use ReceiverFullImpedance
    use ReceiverFullVerticalMagnetic
    use ReceiverOffDiagonalImpedance
    use ReceiverSingleField
    !
    use PredictedDataHandle
    !
    type( MPI_Comm ) :: main_comm, child_comm
    !
    integer :: mpi_rank, mpi_size, node_rank, node_size, &
               nodestringlen, ierr
    !
    character*( MPI_MAX_PROCESSOR_NAME ) :: node_name
    !
    type( MPI_Win )             :: shared_window
    integer( MPI_ADDRESS_KIND ) :: shared_window_size
    integer                     :: shared_disp_unit = 1
    !
    character, dimension(:), pointer :: shared_buffer
    character, dimension(:), pointer :: fwd_info_buffer
    character, dimension(:), pointer :: predicted_data_buffer
    !
    integer :: shared_buffer_size = 1
    integer :: fwd_info_buffer_size = 1
    integer :: predicted_data_buffer_size = 1
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
    integer    :: transmitter_derived_type, transmitters_size
    integer, parameter :: transmitter_mt = 1
    integer, parameter :: transmitter_csem = 2
    !
    integer    :: receiver_derived_type, receivers_size
    integer, parameter :: receiver_full_impedance = 1
    integer, parameter :: receiver_full_vertical_magnetic = 2
    integer, parameter :: receiver_off_diagonal_impedance = 3
    integer, parameter :: receiver_single_field = 4
    !
    ! PROGRAM GLOBAL VARIABLES
    integer :: tag = 2022, master_id = 0
    !
    character*20    :: job_master = "MASTER_JOB", job_finish = "STOP_JOBS", job_fwd_done = "FINISH_FWD_JOB"
    character*20    :: job_share_memory = "SHARE_MEMORY", job_forward = "JOB_FORWARD"
    !
    ! STRUCT job_info
    type :: FWDInfo_t
        !
        SEQUENCE
        !
        character*20 :: job_name
        integer      :: worker_rank
        integer      :: tx_index
        integer      :: n_data
        integer      :: data_size
        logical      :: tx_changed
        character*25 :: forward_solver_type
        character*15 :: source_type
        !
    end type FWDInfo_t
    !
    type( FWDInfo_t ), save :: fwd_info
    !
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
        integer :: i, last_size, nbytes(2)
        !
        shared_buffer_size = shared_buffer_size + allocateGridBuffer( grid )
        !
        write( *, "(A50, i8)" ) "MPI Allocated grid size:", shared_buffer_size
        last_size = shared_buffer_size
        !
        shared_buffer_size = shared_buffer_size + allocateModelOperatorBuffer( model_operator )
        !
        write( *, "(A50, i8)" ) "MPI Allocated model_operator size:", shared_buffer_size - last_size
        last_size = shared_buffer_size
        !
        shared_buffer_size = shared_buffer_size + allocateModelParameterBuffer( model_parameter )
        !
        write( *, "(A50, i8)" ) "MPI Allocated model_parameter size:", shared_buffer_size - last_size
        last_size = shared_buffer_size
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, child_comm, nbytes(1), ierr )
        !
        do i = 1, size( transmitters )
            shared_buffer_size = shared_buffer_size + allocateTransmitterBuffer( getTransmitter( i ) )
            !
            !write( *, "(A50, i8)" ) "MPI Allocated tx size: ", shared_buffer_size - last_size
            !last_size = shared_buffer_size
            !
        end do
        !
        write( *, "(A50, i8)" ) "MPI Allocated transmitters size:", shared_buffer_size - last_size
        last_size = shared_buffer_size
        !
        call MPI_PACK_SIZE( 1, MPI_INTEGER, child_comm, nbytes(2), ierr )
        !
        do i = 1, size( receivers )
            shared_buffer_size = shared_buffer_size + allocateReceiverBuffer( getReceiver(i) )
            !
            !write( *, "(A50, i8)" ) "MPI Allocated rx size:", shared_buffer_size - last_size
            !last_size = shared_buffer_size
            !
        end do
        !
        write( *, "(A50, i8)" ) "MPI Allocated receivers size:", shared_buffer_size - last_size
        !
        do i = 1, size( nbytes )
            shared_buffer_size = shared_buffer_size + nbytes(i)
        end do
        !
        write( *, "(A50, i8)" ) "MPI Allocated total size =", shared_buffer_size
        !
        if( .NOT. associated( shared_buffer ) ) then
            allocate( shared_buffer( shared_buffer_size ) )
        endif
        !
    end subroutine allocateSharedBuffer
    !
    !
    subroutine packSharedBuffer( grid, model_operator, model_parameter )
        implicit none
        !
        class( Grid_t ), intent( in )           :: grid
        class( ModelOperator_t ), intent( in )  :: model_operator
        class( ModelParameter_t ), intent( in ) :: model_parameter
        !
        integer :: i, index
        !
        index = 1
        !
        call packGridBuffer( grid, index )
        !
        call packModelOperatorBuffer( model_operator, index )
        !
        call packModelParameterBuffer( model_parameter, index )
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
    subroutine unpackSharedBuffer( buffer_size, grid, model_operator, model_parameter )
        implicit none
        !
        class( Grid_t ), allocatable, intent( inout )           :: grid
        class( ModelOperator_t ), allocatable, intent( inout )  :: model_operator
        class( ModelParameter_t ), allocatable, intent( inout ) :: model_parameter
        !
        class( Transmitter_t ), allocatable :: transmitter
        class( Receiver_t ), allocatable    :: receiver
        !
        !
        integer, intent( in ) :: buffer_size
        !
        integer :: i, aux_size, index
        !
        index = 1
        !
        shared_buffer_size = buffer_size
        !
        grid = unpackGridBuffer( index )
        !
        model_operator = unpackModelOperatorBuffer( grid, index )
        !
        model_parameter = unpackModelParameterBuffer( grid, model_operator%metric, index )
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, aux_size, 1, MPI_INTEGER, child_comm, ierr )
        !
        do i = 1, aux_size
            !
            call unpackTransmitterBuffer( grid, transmitter, index )
            !
            call updateTransmitterArray( transmitter )
            !
            deallocate( transmitter )
            !
        end do
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, aux_size, 1, MPI_INTEGER, child_comm, ierr )
        !
        do i = 1, aux_size
            !
            call unpackReceiverBuffer( grid, receiver, index )
            !
            call updateReceiverArray( receiver )
            !
            deallocate( receiver )
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
              call MPI_PACK_SIZE( scalar%Nxyz, MPI_DOUBLE_PRECISION, child_comm, nbytes(5), ierr )
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
              call MPI_PACK( scalar%gridType, 4, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%gridType, 4, MPI_CHARACTER, child_comm, ierr )
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
              call MPI_PACK_SIZE( scalar%Nxyz, MPI_DOUBLE_COMPLEX, child_comm, nbytes(5), ierr )
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
              call MPI_PACK( scalar%gridType, 4, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, scalar%gridType, 4, MPI_CHARACTER, child_comm, ierr )
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
              call MPI_PACK_SIZE( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_PRECISION, child_comm, nbytes(5), ierr )
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
              call MPI_PACK( vector%gridType, 4, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%gridType, 4, MPI_CHARACTER, child_comm, ierr )
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
              call MPI_PACK_SIZE( vector%Nxyz(1) + vector%Nxyz(2) + vector%Nxyz(3), MPI_DOUBLE_COMPLEX, child_comm, nbytes(5), ierr )
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
              call MPI_PACK( vector%gridType, 4, MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, vector%gridType, 4, MPI_CHARACTER, child_comm, ierr )
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
        call MPI_PACK_SIZE( 4, MPI_DOUBLE_PRECISION, child_comm, nbytes(3), ierr )
        call MPI_PACK_SIZE( 5, MPI_INTEGER, child_comm, nbytes(4), ierr )
        call MPI_PACK_SIZE( size( grid%dx ), MPI_DOUBLE_PRECISION, child_comm, nbytes(5), ierr )
        call MPI_PACK_SIZE( size( grid%dy ), MPI_DOUBLE_PRECISION, child_comm, nbytes(6), ierr )
        call MPI_PACK_SIZE( size( grid%dz ), MPI_DOUBLE_PRECISION, child_comm, nbytes(7), ierr )
        call MPI_PACK_SIZE( size( grid%dxInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(8), ierr )
        call MPI_PACK_SIZE( size( grid%dyInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(9), ierr )
        call MPI_PACK_SIZE( size( grid%dzInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(10), ierr )
        call MPI_PACK_SIZE( size( grid%delX ), MPI_DOUBLE_PRECISION, child_comm, nbytes(11), ierr )
        call MPI_PACK_SIZE( size( grid%delY ), MPI_DOUBLE_PRECISION, child_comm, nbytes(12), ierr )
        call MPI_PACK_SIZE( size( grid%delZ ), MPI_DOUBLE_PRECISION, child_comm, nbytes(13), ierr )
        call MPI_PACK_SIZE( size( grid%delXInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(14), ierr )
        call MPI_PACK_SIZE( size( grid%delYInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(15), ierr )
        call MPI_PACK_SIZE( size( grid%delZInv ), MPI_DOUBLE_PRECISION, child_comm, nbytes(16), ierr )
        call MPI_PACK_SIZE( size( grid%xEdge ), MPI_DOUBLE_PRECISION, child_comm, nbytes(17), ierr )
        call MPI_PACK_SIZE( size( grid%yEdge ), MPI_DOUBLE_PRECISION, child_comm, nbytes(18), ierr )
        call MPI_PACK_SIZE( size( grid%zEdge ), MPI_DOUBLE_PRECISION, child_comm, nbytes(19), ierr )
        call MPI_PACK_SIZE( size( grid%xCenter ), MPI_DOUBLE_PRECISION, child_comm, nbytes(20), ierr )
        call MPI_PACK_SIZE( size( grid%yCenter ), MPI_DOUBLE_PRECISION, child_comm, nbytes(21), ierr )
        call MPI_PACK_SIZE( size( grid%zCenter ), MPI_DOUBLE_PRECISION, child_comm, nbytes(22), ierr )
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
                call MPI_PACK( grid%ox, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%oy, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%oz, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%rotDeg, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%nx, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%ny, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%nz, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%nzAir, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%nzEarth, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( grid%dx(1), size( grid%dx ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dy(1), size( grid%dy ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dz(1), size( grid%dz ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dxInv(1), size( grid%dxInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dyInv(1), size( grid%dyInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%dzInv(1), size( grid%dzInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delX(1), size( grid%delX ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delY(1), size( grid%delY ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delZ(1), size( grid%delZ ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delXInv(1), size( grid%delXInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delYInv(1), size( grid%delYInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%delZInv(1), size( grid%delZInv ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%xEdge(1), size( grid%xEdge ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%yEdge(1), size( grid%yEdge ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%zEdge(1), size( grid%zEdge ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%xCenter(1), size( grid%xCenter ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%yCenter(1), size( grid%yCenter ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%zCenter(1), size( grid%zCenter ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( grid%zAirThick, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( grid%allocated, 1, MPI_LOGICAL, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
           class default
              stop "packGridBuffer: Unclassified grid"
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
        integer :: i, nbytes(18), model_operator_size_bytes
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
                call MPI_PACK_SIZE( 2, MPI_LOGICAL, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 11, MPI_INTEGER, child_comm, nbytes(3), ierr )
                call MPI_PACK_SIZE( size( model_operator%xXY ), MPI_DOUBLE_PRECISION, child_comm, nbytes(4), ierr )
                call MPI_PACK_SIZE( size( model_operator%xXZ ), MPI_DOUBLE_PRECISION, child_comm, nbytes(5), ierr )
                call MPI_PACK_SIZE( size( model_operator%xY ), MPI_DOUBLE_PRECISION, child_comm, nbytes(6), ierr )
                call MPI_PACK_SIZE( size( model_operator%xZ ), MPI_DOUBLE_PRECISION, child_comm, nbytes(7), ierr )
                call MPI_PACK_SIZE( size( model_operator%xXO ), MPI_DOUBLE_PRECISION, child_comm, nbytes(8), ierr )
                call MPI_PACK_SIZE( size( model_operator%yYX ), MPI_DOUBLE_PRECISION, child_comm, nbytes(9), ierr )
                call MPI_PACK_SIZE( size( model_operator%yYZ ), MPI_DOUBLE_PRECISION, child_comm, nbytes(10), ierr )
                call MPI_PACK_SIZE( size( model_operator%yX ), MPI_DOUBLE_PRECISION, child_comm, nbytes(11), ierr )
                call MPI_PACK_SIZE( size( model_operator%yZ ), MPI_DOUBLE_PRECISION, child_comm, nbytes(12), ierr )
                call MPI_PACK_SIZE( size( model_operator%yYO ), MPI_DOUBLE_PRECISION, child_comm, nbytes(13), ierr )
                call MPI_PACK_SIZE( size( model_operator%zZX ), MPI_DOUBLE_PRECISION, child_comm, nbytes(14), ierr )
                call MPI_PACK_SIZE( size( model_operator%zZY ), MPI_DOUBLE_PRECISION, child_comm, nbytes(15), ierr )
                call MPI_PACK_SIZE( size( model_operator%zX ), MPI_DOUBLE_PRECISION, child_comm, nbytes(16), ierr )
                call MPI_PACK_SIZE( size( model_operator%zY ), MPI_DOUBLE_PRECISION, child_comm, nbytes(17), ierr )
                call MPI_PACK_SIZE( size( model_operator%zZO ), MPI_DOUBLE_PRECISION, child_comm, nbytes(18), ierr )
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
    !
    subroutine packModelOperatorBuffer( model_operator, index )
        implicit none
        !
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
                call MPI_PACK( r_array(1), product( shape( model_operator%xXY ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%xXZ )
                call MPI_PACK( r_array(1), product( shape( model_operator%xXZ ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%xY )
                call MPI_PACK( r_array(1), product( shape( model_operator%xY ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%xZ )
                call MPI_PACK( r_array(1), product( shape( model_operator%xZ ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%xXO )
                call MPI_PACK( r_array(1), product( shape( model_operator%xXO ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%yYX )
                call MPI_PACK( r_array(1), product( shape( model_operator%yYX ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%yYZ )
                call MPI_PACK( r_array(1), product( shape( model_operator%yYZ ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%yX )
                call MPI_PACK( r_array(1), product( shape( model_operator%yX ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%yZ )
                call MPI_PACK( r_array(1), product( shape( model_operator%yZ ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%yYO )
                call MPI_PACK( r_array(1), product( shape( model_operator%yYO ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%zZX )
                call MPI_PACK( r_array(1), product( shape( model_operator%zZX ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%zZY )
                call MPI_PACK( r_array(1), product( shape( model_operator%zZY ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%zX )
                call MPI_PACK( r_array(1), product( shape( model_operator%zX ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%zY )
                call MPI_PACK( r_array(1), product( shape( model_operator%zY ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                r_array = BiArrayToArray( model_operator%zZO )
                call MPI_PACK( r_array(1), product( shape( model_operator%zZO ) ), MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
    !
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
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_xXY_dim1 * model_operator_xXY_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%xXY = arrayToBiArray( r_array, model_operator_xXY_dim1, model_operator_xXY_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_xXZ_dim1 * model_operator_xXZ_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_xXZ_dim1 * model_operator_xXZ_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%xXZ = arrayToBiArray( r_array, model_operator_xXZ_dim1, model_operator_xXZ_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_xY_dim1 * model_operator_xY_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_xY_dim1 * model_operator_xY_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%xY = arrayToBiArray( r_array, model_operator_xY_dim1, model_operator_xY_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_xZ_dim1 * model_operator_xZ_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_xZ_dim1 * model_operator_xZ_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%xZ = arrayToBiArray( r_array, model_operator_xZ_dim1, model_operator_xZ_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_xXO_dim1 * model_operator_xXO_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_xXO_dim1 * model_operator_xXO_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%xXO = arrayToBiArray( r_array, model_operator_xXO_dim1, model_operator_xXO_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_yYX_dim1 * model_operator_yYX_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_yYX_dim1 * model_operator_yYX_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%yYX = arrayToBiArray( r_array, model_operator_yYX_dim1, model_operator_yYX_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_yYZ_dim1 * model_operator_yYZ_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_yYZ_dim1 * model_operator_yYZ_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%yYZ = arrayToBiArray( r_array, model_operator_yYZ_dim1, model_operator_yYZ_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_yX_dim1 * model_operator_yX_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_yX_dim1 * model_operator_yX_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%yX = arrayToBiArray( r_array, model_operator_yX_dim1, model_operator_yX_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_yZ_dim1 * model_operator_yZ_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_yZ_dim1 * model_operator_yZ_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%yZ = arrayToBiArray( r_array, model_operator_yZ_dim1, model_operator_yZ_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_yYO_dim1 * model_operator_yYO_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_yYO_dim1 * model_operator_yYO_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%yYO = arrayToBiArray( r_array, model_operator_yYO_dim1, model_operator_yYO_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_zZX_dim1 * model_operator_zZX_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_zZX_dim1 * model_operator_zZX_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%zZX = arrayToBiArray( r_array, model_operator_zZX_dim1, model_operator_zZX_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_zZY_dim1 * model_operator_zZY_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_zZY_dim1 * model_operator_zZY_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%zZY = arrayToBiArray( r_array, model_operator_zZY_dim1, model_operator_zZY_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_zX_dim1 * model_operator_zX_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_zX_dim1 * model_operator_zX_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%zX = arrayToBiArray( r_array, model_operator_zX_dim1, model_operator_zX_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_zY_dim1 * model_operator_zY_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_zY_dim1 * model_operator_zY_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        model_operator%zY = arrayToBiArray( r_array, model_operator_zY_dim1, model_operator_zY_dim2 )
                        deallocate( r_array )
                        !
                        allocate( r_array( model_operator_zZO_dim1 * model_operator_zZO_dim2 ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, r_array(1), model_operator_zZO_dim1 * model_operator_zZO_dim2, MPI_DOUBLE_PRECISION, child_comm, ierr )
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
                call MPI_PACK( model_parameter%airCond, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
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
    !
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
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, model_parameter%airCond, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
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
                allocate( nbytes(5) )
                !
                call MPI_PACK_SIZE( 15, MPI_INTEGER, child_comm, nbytes(1), ierr )
                !
                call MPI_PACK_SIZE( len( transmitter%type_name ), MPI_CHARACTER, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, child_comm, nbytes(3), ierr )
                !
                !do i = 1, size( transmitter%e_all )
                   !transmitter_size_bytes = transmitter_size_bytes + allocateCVectorBuffer( transmitter%e_all(i) )
                !end do
                !
                call MPI_PACK_SIZE( size( transmitter%receiver_indexes ), MPI_INTEGER, child_comm, nbytes(4), ierr )
                call MPI_PACK_SIZE( len( transmitter%DATA_TITLE ), MPI_CHARACTER, child_comm, nbytes(5), ierr )
                !
            class is( TransmitterCSEM_t )
                !
                allocate( nbytes(7) )
                !
                call MPI_PACK_SIZE( 15, MPI_INTEGER, child_comm, nbytes(1), ierr )
                !
                call MPI_PACK_SIZE( len( transmitter%type_name ), MPI_CHARACTER, child_comm, nbytes(2), ierr )
                call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, child_comm, nbytes(3), ierr )
                !
                !do i = 1, size( transmitter%e_all )
                   !transmitter_size_bytes = transmitter_size_bytes + allocateCVectorBuffer( transmitter%e_all(i) )
                !end do
                !
                call MPI_PACK_SIZE( size( transmitter%receiver_indexes ), MPI_INTEGER, child_comm, nbytes(4), ierr )
                call MPI_PACK_SIZE( len( transmitter%DATA_TITLE ), MPI_CHARACTER, child_comm, nbytes(5), ierr )
                call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, child_comm, nbytes(6), ierr )
                call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, child_comm, nbytes(7), ierr )
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
                ! SYZES
                call MPI_PACK( transmitter%id, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%n_pol, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%fwd_key(1), 8, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( len( transmitter%type_name ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( transmitter%e_all ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( size( transmitter%receiver_indexes ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( len( transmitter%DATA_TITLE ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%type_name, len( transmitter%type_name ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                call MPI_PACK( transmitter%period, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                !do i = 1, size( transmitter%e_all )
                   !call packCVectorBuffer( transmitter%e_all(i), index )
                !end do
                !
                call MPI_PACK( transmitter%receiver_indexes(1), size( transmitter%receiver_indexes ), MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%DATA_TITLE, len( transmitter%DATA_TITLE ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
           class is( TransmitterCSEM_t )
                !
                ! TYPE
                call MPI_PACK( transmitter_csem, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                ! SYZES
                call MPI_PACK( transmitter%id, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%n_pol, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%fwd_key(1), 8, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( len( transmitter%type_name ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( transmitter%e_all ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( size( transmitter%receiver_indexes ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( len( transmitter%DATA_TITLE ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%type_name, len( transmitter%type_name ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%period, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
                !do i = 1, size( transmitter%e_all )
                   !call packCVectorBuffer( transmitter%e_all(i), index )
                !end do
                !
                call MPI_PACK( transmitter%receiver_indexes(1), size( transmitter%receiver_indexes ), MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%DATA_TITLE, len( transmitter%DATA_TITLE ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%location(1), 3, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                call MPI_PACK( transmitter%azimuth, 1, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
           class default
              stop "allocateTransmitterBuffer: Unclassified transmitter"
           !
        end select

    end subroutine packTransmitterBuffer
    !
    !
    subroutine unpackTransmitterBuffer( grid, transmitter, index )
        implicit none
        !
        class( Grid_t ), target, intent( in )                :: grid
        class( Transmitter_t ), allocatable, intent( inout ) :: transmitter
        integer, intent( inout )                             :: index
        !
        integer :: i, transmitter_type, transmitter_e_all, transmitter_receiver_indexes, transmitter_DATA_TITLE
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
                        ! SYZES
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%id, 1, MPI_INTEGER, child_comm, ierr )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%n_pol, 1, MPI_INTEGER, child_comm, ierr )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%fwd_key(1), 8, MPI_INTEGER, child_comm, ierr )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_type, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_e_all, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_receiver_indexes, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_DATA_TITLE, 1, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( character( transmitter_type ) :: transmitter%type_name )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%type_name, transmitter_type, MPI_CHARACTER, child_comm, ierr )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%period, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( cVector3D_SG_t :: transmitter%e_all( transmitter_e_all ) )
                        !do i = 1, transmitter_e_all
                            !write(*,*) unpackCVectorBuffer( grid, index )
                        !end do
                        !
                        allocate( transmitter%receiver_indexes( transmitter_receiver_indexes ) )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%receiver_indexes(1), transmitter_receiver_indexes, MPI_INTEGER, child_comm, ierr )
                        !
                        allocate( character( transmitter_DATA_TITLE ) :: transmitter%DATA_TITLE )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%DATA_TITLE, transmitter_DATA_TITLE, MPI_CHARACTER, child_comm, ierr )
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
                        ! SYZES
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%id, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%n_pol, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%fwd_key(1), 8, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_type, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_e_all, 1, MPI_INTEGER, child_comm, ierr )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_receiver_indexes, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter_DATA_TITLE, 1, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%type_name, transmitter_type, MPI_CHARACTER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%period, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        !
                        allocate( cVector3D_SG_t :: transmitter%e_all( transmitter_e_all ) )
                        !do i = 1, transmitter_e_all
                            !allocate( transmitter%e_all(i), source = unpackCVectorBuffer( grid, index ) )
                        !end do
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%receiver_indexes(1), transmitter_receiver_indexes, MPI_INTEGER, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%DATA_TITLE, transmitter_DATA_TITLE, MPI_CHARACTER, child_comm, ierr )
                        !
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%location(1), 3, MPI_DOUBLE_PRECISION, child_comm, ierr )
                        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, transmitter%azimuth, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
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
    end subroutine unpackTransmitterBuffer
    !
    !
    function allocateReceiverBuffer( receiver ) result( receiver_size_bytes )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: receiver
        integer :: i, receiver_size_bytes
        integer :: nbytes(3)
        !
        receiver_size_bytes = 0
        !
        call MPI_PACK_SIZE( 3, MPI_INTEGER, child_comm, nbytes(1), ierr )
        call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, child_comm, nbytes(2), ierr )
        call MPI_PACK_SIZE( len( receiver%code ), MPI_CHARACTER, child_comm, nbytes(3), ierr )
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
           class is( ReceiverFullVerticalMagnetic_t )
                !
                call MPI_PACK( receiver_full_vertical_magnetic, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
            class is( ReceiverOffDiagonalImpedance_t )
                !
                call MPI_PACK( receiver_off_diagonal_impedance, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
            class is( ReceiverSingleField_t )
                !
                call MPI_PACK( receiver_single_field, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
                !
           class default
              stop "allocateReceiverBuffer: Unclassified receiver"
           !
        end select
        !
        call MPI_PACK( receiver%id, 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( len( receiver%code ), 1, MPI_INTEGER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        call MPI_PACK( receiver%code, len( receiver%code ), MPI_CHARACTER, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        !
        call MPI_PACK( receiver%location(1), 3, MPI_DOUBLE_PRECISION, shared_buffer, shared_buffer_size, index, child_comm, ierr )
        !
    end subroutine packReceiverBuffer
    !
    !
    subroutine unpackReceiverBuffer( grid, receiver, index )
        implicit none
        !
        class( Grid_t ), target, intent( in )             :: grid
        class( Receiver_t ), allocatable, intent( inout ) :: receiver
        integer, intent( inout )                          :: index
        !
        integer :: receiver_id, code_size
        !
        character(:), allocatable :: code
        real( kind=prec )         :: receiver_location(3)
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver_derived_type, 1, MPI_INTEGER, child_comm, ierr )
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver_id, 1, MPI_INTEGER, child_comm, ierr )
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, code_size, 1, MPI_INTEGER, child_comm, ierr )
        !
        allocate( character( code_size ) :: code )
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, code, code_size, MPI_CHARACTER, child_comm, ierr )
        !
        call MPI_UNPACK( shared_buffer, shared_buffer_size, index, receiver_location(1), 3, MPI_DOUBLE_PRECISION, child_comm, ierr )
        !
        select case( receiver_derived_type )
           !
           case( receiver_full_impedance )
                !
                allocate( receiver, source = ReceiverFullImpedance_t( receiver_location ) )
                !
           case( receiver_full_vertical_magnetic )
                !
                allocate( receiver, source = ReceiverFullVerticalMagnetic_t( receiver_location ) )
                !
           case( receiver_off_diagonal_impedance )
                !
                allocate( receiver, source = ReceiverOffDiagonalImpedance_t( receiver_location ) )
                !
           case( receiver_single_field )
                !
                allocate( receiver, source = ReceiverFullImpedance_t( receiver_location ) )
                !
           case default
              stop "allocateReceiverBuffer: Unclassified receiver"
           !
        end select
        !
        receiver%id = receiver_id
        receiver%code = code
        !
    end subroutine unpackReceiverBuffer
    !
    subroutine createDataBuffer()
        implicit none
        !
        predicted_data_buffer_size = fwd_info%data_size
        !
        if( .NOT. associated( predicted_data_buffer ) ) then
            allocate( predicted_data_buffer( predicted_data_buffer_size ) )
        endif
        !
    end subroutine createDataBuffer
    !
    subroutine allocateDataBuffer( data_entries )
        implicit none
        !
        type( PredictedDataHandle_t ), dimension(:), intent( in ) :: data_entries
        !
        integer i, j, nbytes(7)
        !
        predicted_data_buffer_size = 1
        !
        do i = 1, size( data_entries )
            !
            call MPI_PACK_SIZE( 3, MPI_INTEGER, child_comm, nbytes(1), ierr )
            call MPI_PACK_SIZE( len( data_entries(i)%code ), MPI_CHARACTER, child_comm, nbytes(2), ierr )
            call MPI_PACK_SIZE( len( data_entries(i)%component ), MPI_CHARACTER, child_comm, nbytes(3), ierr )
            call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, child_comm, nbytes(4), ierr )
            call MPI_PACK_SIZE( 3, MPI_DOUBLE_PRECISION, child_comm, nbytes(5), ierr )
            call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, child_comm, nbytes(6), ierr )
            call MPI_PACK_SIZE( 1, MPI_DOUBLE_PRECISION, child_comm, nbytes(7), ierr )
            !
            do j = 1, size( nbytes )
                predicted_data_buffer_size = predicted_data_buffer_size + nbytes(j)
            end do
            !
        end do
        !
        if( .NOT. associated( predicted_data_buffer ) ) then
            allocate( predicted_data_buffer( predicted_data_buffer_size ) )
        endif
        !
    end subroutine allocateDataBuffer
    !
    !
    subroutine packDataBuffer( data_entries )
        implicit none
        !
        type( PredictedDataHandle_t ), dimension(:), intent( in ) :: data_entries
        !
        integer :: i, index
        !
        index = 1
        !
        do i = 1, size( data_entries )
            !
            call MPI_PACK( data_entries(i)%rx_id, 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
            call MPI_PACK( len( data_entries(i)%code ), 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
            call MPI_PACK( len( data_entries(i)%component ), 1, MPI_INTEGER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
            call MPI_PACK( data_entries(i)%code, len( data_entries(i)%code ), MPI_CHARACTER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
            call MPI_PACK( data_entries(i)%component, len( data_entries(i)%component ), MPI_CHARACTER, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
            !
            call MPI_PACK( data_entries(i)%period, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
            call MPI_PACK( data_entries(i)%xyz(1), 3, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
            call MPI_PACK( data_entries(i)%real, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
            call MPI_PACK( data_entries(i)%imaginary, 1, MPI_DOUBLE_PRECISION, predicted_data_buffer, predicted_data_buffer_size, index, child_comm, ierr )
            !
        end do
        !
    end subroutine packDataBuffer
    !
    ! UNPACK predicted_data_buffer TO predicted_data STRUCT
    function unpackDataBuffer() result( data_entries )
        implicit none
        !
        type( PredictedDataHandle_t ), allocatable, dimension(:) :: data_entries
        !
        integer :: i, data_entries_code, data_entries_component, index
        !
        index = 1
        !
        if( allocated( data_entries ) ) deallocate( data_entries )
        allocate( data_entries( fwd_info%n_data ) )
        !
        do i = 1, size( data_entries )
            !
            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_entries(i)%rx_id, 1, MPI_INTEGER, child_comm, ierr )
            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_entries_code, 1, MPI_INTEGER, child_comm, ierr )
            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_entries_component, 1, MPI_INTEGER, child_comm, ierr )
            !
            allocate( character( data_entries_code ) :: data_entries(i)%code )
            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_entries(i)%code, data_entries_code, MPI_CHARACTER, child_comm, ierr )
            !
            allocate( character( data_entries_component ) :: data_entries(i)%component )
            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_entries(i)%component, data_entries_component, MPI_CHARACTER, child_comm, ierr )
            !
            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_entries(i)%period, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_entries(i)%xyz(1), 3, MPI_DOUBLE_PRECISION, child_comm, ierr )
            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_entries(i)%real, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
            call MPI_UNPACK( predicted_data_buffer, predicted_data_buffer_size, index, data_entries(i)%imaginary, 1, MPI_DOUBLE_PRECISION, child_comm, ierr )
            !
        end do
        !
    end function unpackDataBuffer
    !
    ! RECEIVE predicted_data FROM ANY TARGET
    function receiveData() result( data_entries )
        implicit none
        !
        type( PredictedDataHandle_t ), allocatable, dimension(:) :: data_entries
        !
        call createDataBuffer
        call MPI_RECV( predicted_data_buffer, predicted_data_buffer_size, MPI_PACKED, fwd_info%worker_rank, MPI_ANY_TAG, child_comm, MPI_STATUS_IGNORE, ierr )
        !
        data_entries = unpackDataBuffer()
        !
    end function receiveData
    !
    ! SEND fwd_info FROM target_id
    subroutine sendData( data_entries )
        !
        type( PredictedDataHandle_t ), allocatable, dimension(:), intent( in ) :: data_entries
        !
        call packDataBuffer( data_entries )
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
        call MPI_PACK_SIZE( 20, MPI_CHARACTER, child_comm, nbytes1, ierr )
        call MPI_PACK_SIZE( 4, MPI_INTEGER, child_comm, nbytes2, ierr )
        call MPI_PACK_SIZE( 1, MPI_LOGICAL, child_comm, nbytes3, ierr )
        call MPI_PACK_SIZE( 40, MPI_CHARACTER, child_comm, nbytes4, ierr )
        !
        fwd_info_buffer_size = ( nbytes1 + nbytes2 + nbytes3 + nbytes4 ) + 1
        !
        if( .NOT. associated( fwd_info_buffer ) ) then
            allocate( fwd_info_buffer( fwd_info_buffer_size ) )
        endif
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
        call MPI_PACK( fwd_info%job_name, 20, MPI_CHARACTER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%worker_rank, 1, MPI_INTEGER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%tx_index, 1, MPI_INTEGER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%n_data, 1, MPI_INTEGER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%data_size, 1, MPI_INTEGER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%tx_changed, 1, MPI_LOGICAL, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%forward_solver_type, 25, MPI_CHARACTER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
        call MPI_PACK( fwd_info%source_type, 15, MPI_CHARACTER, fwd_info_buffer, fwd_info_buffer_size, index, child_comm, ierr )
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
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%job_name, 20, MPI_CHARACTER, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%worker_rank, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%tx_index, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%n_data, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%data_size, 1, MPI_INTEGER, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%tx_changed, 1, MPI_LOGICAL, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%forward_solver_type, 25, MPI_CHARACTER, child_comm, ierr )
        call MPI_UNPACK( fwd_info_buffer, fwd_info_buffer_size, index, fwd_info%source_type, 15, MPI_CHARACTER, child_comm, ierr )
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