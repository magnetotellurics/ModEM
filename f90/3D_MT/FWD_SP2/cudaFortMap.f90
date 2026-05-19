module cudaFortMap
   ! A BAND-AID MODULE TO MAP CUDA C INTERFACES TO FORTRAN
   ! WITH ISO_C_BINDING
   use math_constants   ! math/ physics constants
   use iso_c_binding

   implicit none
   save
   ! ======================= CUDA enumerators ======================== !
   ! for C these are already setup in cudart/cublas/cusparse headers
   ! however, our fortran subroutines know nothing about those
   ! as such they need to be explicitly setup here
   ! memcpy to host or to device
   integer*8         :: cudaMemcpyDeviceToHost 
   integer*8         :: cudaMemcpyHostToDevice
   integer*8         :: cudaMemcpyDeviceToDevice
   parameter (cudaMemcpyHostToDevice=1)
   parameter (cudaMemcpyDeviceToHost=2)
   parameter (cudaMemcpyDeviceToDevice=3)
   integer*4         :: cudaSuccess
   parameter (cudaSuccess=0)
   ! matrix operation (whether do transpose)
   ! note this is important as the tests confirm the non-transposed
   ! operation is much faster than transposed ones on GPU
   ! need to avoid the transposed operations (that's why I didn't use QMR)
   integer*4         :: CUSPARSE_OPERATION_NON_TRANSPOSE
   integer*4         :: CUSPARSE_OPERATION_TRANSPOSE
   integer*4         :: CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE
   parameter (CUSPARSE_OPERATION_NON_TRANSPOSE=0)
   parameter (CUSPARSE_OPERATION_TRANSPOSE=1)
   parameter (CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE=2)
   ! matrix type (symmetric may save some spaces)
   integer*4         :: CUSPARSE_MATRIX_TYPE_GENERAL
   integer*4         :: CUSPARSE_MATRIX_TYPE_SYMMETRIC
   integer*4         :: CUSPARSE_MATRIX_TYPE_HERMITIAN
   integer*4         :: CUSPARSE_MATRIX_TYPE_TRIANGULAR
   parameter (CUSPARSE_MATRIX_TYPE_GENERAL=0)
   parameter (CUSPARSE_MATRIX_TYPE_SYMMETRIC=1)
   parameter (CUSPARSE_MATRIX_TYPE_HERMITIAN=2)
   parameter (CUSPARSE_MATRIX_TYPE_TRIANGULAR=3)
   ! index base - we all use 1 in fortran
   integer*4         :: CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_INDEX_BASE_ONE
   parameter (CUSPARSE_INDEX_BASE_ZERO=0)
   parameter (CUSPARSE_INDEX_BASE_ONE=1)
   ! fill mode - upper or lower triangular 
   integer*4         :: CUSPARSE_FILL_MODE_LOWER
   integer*4         :: CUSPARSE_FILL_MODE_UPPER
   parameter (CUSPARSE_FILL_MODE_LOWER=0) 
   parameter (CUSPARSE_FILL_MODE_UPPER=1)
   ! mat index type
   integer*4         :: CUSPARSE_INDEX_16U, CUSPARSE_INDEX_32I
   integer*4         :: CUSPARSE_INDEX_64I
   parameter (CUSPARSE_INDEX_16U=1)
   parameter (CUSPARSE_INDEX_32I=2)
   parameter (CUSPARSE_INDEX_64I=3)
   ! mat value type
   integer*4         :: CUDA_R_32F, CUDA_C_32F, CUDA_R_64F, CUDA_C_64F
   parameter (CUDA_R_32F=0)
   parameter (CUDA_C_32F=4)
   parameter (CUDA_R_64F=1)
   parameter (CUDA_C_64F=5)
   ! diag type  unit or not
   integer*4         :: CUSPARSE_DIAG_TYPE_UNIT 
   integer*4         :: CUSPARSE_DIAG_TYPE_NON_UNIT 
   parameter (CUSPARSE_DIAG_TYPE_UNIT=1)
   parameter (CUSPARSE_DIAG_TYPE_NON_UNIT=0)
   ! Algorithm type for SpMV
   integer*4         :: CUSPARSE_SPMV_ALG_DEFAULT
   integer*4         :: CUSPARSE_SPMV_COO_ALG1
   integer*4         :: CUSPARSE_SPMV_CSR_ALG1
   integer*4         :: CUSPARSE_SPMV_CSR_ALG2
   integer*4         :: CUSPARSE_SPMV_COO_ALG2
   parameter (CUSPARSE_SPMV_ALG_DEFAULT=0)
   parameter (CUSPARSE_SPMV_COO_ALG1=1)
   parameter (CUSPARSE_SPMV_CSR_ALG1=2) 
   parameter (CUSPARSE_SPMV_CSR_ALG2=3)
   parameter (CUSPARSE_SPMV_COO_ALG2=4)
   ! Algorithm type for SpSV
   integer*4         :: CUSPARSE_SPSV_ALG_DEFAULT
   parameter (CUSPARSE_SPSV_ALG_DEFAULT=0)
   ! policy for cusparse solve
   integer*4         :: CUSPARSE_SOLVE_POLICY_NO_LEVEL
   integer*4         :: CUSPARSE_SOLVE_POLICY_USE_LEVEL
   parameter (CUSPARSE_SOLVE_POLICY_NO_LEVEL=0)
   parameter (CUSPARSE_SOLVE_POLICY_USE_LEVEL=1)
   ! matrix attributes for SpMat
   integer*4         :: CUSPARSE_SPMAT_FILL_MODE
   integer*4         :: CUSPARSE_SPMAT_DIAG_TYPE
   parameter (CUSPARSE_SPMAT_FILL_MODE=0)
   parameter (CUSPARSE_SPMAT_DIAG_TYPE=1)
   ! cudaStream flag    
   integer*4         :: cudaStreamNonBlocking
   integer*4         :: cudaStreamDefault
   parameter (cudaStreamNonBlocking=1)
   parameter (cudaStreamDefault=0)
   ! nccl enumerators
   integer*4         :: ncclInt8
   integer*4         :: ncclChar
   integer*4         :: ncclInt
   integer*4         :: ncclFloat
   integer*4         :: ncclFloat64
   integer*4         :: ncclDouble
   integer*4         :: NCCL_UNIQUE_ID_BYTES
   parameter (ncclInt8=0)
   parameter (ncclChar=0)
   parameter (ncclInt=2)
   parameter (ncclFloat=7)
   parameter (ncclFloat64=8)
   parameter (ncclDouble=8)
   parameter (NCCL_UNIQUE_ID_BYTES=128)
   ! index base - we all use 1 in fortran
   ! additional C_ONE in FP32 - useful in mixed precision calculations
   complex(kind=SP)  :: C_ONEs
   parameter (C_ONEs=(1.0, 0.0))
   ! ======================= NCCL C derived types ===================== !
   ! ncclUniqueId
   type, bind(c) :: ncclUniqueId
       character(c_char) :: internal(NCCL_UNIQUE_ID_BYTES)
   end type
   ! ncclComm
   type, bind(c) :: ncclComm
       type(c_ptr) :: member
   end type
   ! ncclResult
   type, bind(c) :: ncclResult
       integer(c_int) :: member
   end type
   ! ======================= CUDA C interfaces ======================== !
   interface
   ! I only included some "might-be-useful" interfaces here...
   ! which ranges from basic cuda memory manipulation to cublas and cusparse
   ! feel free to add more if you need other interfaces, but keep the 
   ! original naming convention.
   ! Hao
   ! 
   ! ==================================================================== !
   ! ====================== BASIC CUDA INTERFACES ======================= !
   ! ==================================================================== !
   ! cudaSetDevice
   integer (c_int) function cudaSetDevice( device_idx ) &
    &              bind (C, name="cudaSetDevice" ) 
     ! set the current working device according to device_idx
     use iso_c_binding
     implicit none
     integer (c_int), value  :: device_idx
   end function cudaSetDevice

   ! cudaGetDeviceCount
   integer (c_int) function cudaGetDeviceCount( dev_count ) &
    &              bind (C, name="cudaGetDeviceCount" ) 
     ! get the number of available devices 
     use iso_c_binding
     implicit none
     type(c_ptr), value  :: dev_count
   end function cudaGetDeviceCount

   ! cudaMemset
   integer (c_int) function cudaMemset( devPtr,value, count ) &
    &              bind (C, name="cudaMemset" ) 
     ! fill the count bytes of device memory pointed to by devPr with value
     ! doesn't seems quite useful unless you are filling everything with 
     ! zeros
     use iso_c_binding
     implicit none
     ! devPtr is the *device* pointer to memory 
     type (c_ptr),value  :: devPtr
     integer(c_int), value :: value
     integer(c_size_t), value :: count
   end function cudaMemset

   ! cudaMalloc
   integer (c_int) function cudaMalloc ( devPtr, count ) &
    &              bind (C, name="cudaMalloc" ) 
     ! allocate count bytes of memory pointed by devPtr 
     use iso_c_binding
     implicit none
     ! devPtr is the *device* pointer to memory 
     type (c_ptr)  :: devPtr
     integer (c_size_t), value :: count
   end function cudaMalloc

   ! cudaFree
   integer (c_int) function cudaFree(devPtr) & 
    &              bind(C, name="cudaFree")
     ! free the chunk of memory used by buffer
     use iso_c_binding
     implicit none
     ! devPtr is the *device* pointer to memory that you want to free
     type (c_ptr),value :: devPtr
   end function cudaFree

   ! cudaMallocHost
   integer (c_int) function cudaMallocHost(hostPtr, count ) &
    &              bind (C, name="cudaMallocHost" ) 
     ! allocate count bytes of memory pointed by devPtr 
     use iso_c_binding
     implicit none
     ! devPtr is the *host* pointer to memory 
     type (c_ptr)  :: hostPtr
     integer (c_size_t), value :: count
   end function cudaMallocHost

   ! cudaFreeHost
   integer (c_int) function cudaFreeHost(hostPtr) & 
    &              bind(C, name="cudaFreeHost")
     ! free the chunk of memory used by buffer
     use iso_c_binding
     implicit none
     ! devPtr is the *host* pointer to memory that you want to free
     type (c_ptr),value :: hostPtr
   end function cudaFreeHost

   ! cudaMemcpy
   integer (c_int) function cudaMemcpy ( dst, src, count, kind ) & 
    &              bind (C, name="cudaMemcpy" )
     use iso_c_binding
     implicit none
     ! copy a chunk of memory from *src* to *dst*
     type (c_ptr), value :: dst, src
     ! count is the size of memory (with unit of size_t)
     integer (c_size_t), value :: count, kind 
     ! kind specifies the direction 
     ! cudaMemcpyHostToDevice = 1
     ! cudaMemcpyDeviceToHost = 2
     ! cudaMemcpyDeviceToDevice = 3
   end function cudaMemcpy

   ! cudaMemcpyAsync
   integer (c_int) function cudaMemcpyAsync( dst, src, count, kind, stream) &
    &              bind (C, name="cudaMemcpyAsync" )
     use iso_c_binding
     implicit none
     ! copy a chunk of memory from *src* to *dst*
     type (c_ptr), value       :: dst, src
     ! count is the size of memory (with unit of size_t)
     integer (c_size_t), value :: count, kind 
     ! kind specifies the direction 
     ! cudaMemcpyHostToDevice = 1
     ! cudaMemcpyDeviceToHost = 2
     ! cudaMemcpyDeviceToDevice = 3
     type(c_ptr),value         :: stream ! streamid
   end function cudaMemcpyAsync

   ! cudaHostRegister
   integer (c_int) function cudaHostRegister(hostPtr, size, flag)  &
    &     bind(C, name="cudaHostRegister")
     use iso_c_binding
     implicit none
     type(c_ptr),value       :: hostPtr  ! host pointer
     integer(c_size_t),value :: size  ! memory in bytes
     integer, value          :: flag ! mem type 
   end function cudaHostRegister

   ! cudaHostUnregister
   integer (c_int) function cudaHostUnregister(hostPtr)  &
    &     bind(C, name="cudaHostUnregister")
     use iso_c_binding
     implicit none
     type(c_ptr),value       :: hostPtr  ! host pointer
   end function cudaHostUnregister


   ! cudaMemGetInfo
   integer (c_int) function cudaMemGetInfo(free, total)  &
    &     bind(C, name="cudaMemGetInfo")
     use iso_c_binding
     implicit none
     type(c_ptr),value :: free  ! free memory in bytes
     type(c_ptr),value :: total ! total memory in bytes
   end function cudaMemGetInfo

   ! cudaDeviceSynchronize
   integer(c_int) function cudaDeviceSynchronize() &
    &            bind(C,name="cudaDeviceSynchronize")
     use iso_c_binding
     implicit none
   end function cudaDeviceSynchronize

   ! cudaStreamCreate
   integer(c_int) function cudaStreamCreate(stream) & 
    &            bind(C,name="cudaStreamCreate")
     use iso_c_binding
     implicit none
     type(c_ptr)::stream ! streamid
   end function cudaStreamCreate

   ! cudaStreamCreateWithFlags
   integer(c_int) function cudaStreamCreateWithFlags(stream,flag ) & 
    &            bind(C,name="cudaStreamCreate")
     ! this creates cuda stream (with non-blocking flag)
     use iso_c_binding
     implicit none
     type(c_ptr)::stream ! streamid, out
     integer(c_int),value :: flag ! in
     ! can be 
     ! cudaStreamDefault     
     ! cudaStreamNonBlocking
   end function cudaStreamCreateWithFlags

   ! cudaStreamWaitEvent
   integer(c_int) function cudaStreamWaitEvent(cstream,event,flag) &
    &            bind(C,name="cudaStreamWaitEvent")
    ! this let the stream to wait for a event to finish  
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: cstream ! streamid, in
    type(c_ptr), value    :: event   ! event ID, in
    integer(c_int),value  :: flag    ! in
   end function cudaStreamWaitEvent

   ! cudaStreamDestroy
   integer(c_int) function cudaStreamDestroy(stream) & 
    &            bind(C,name="cudaStreamDestroy")
     ! this eliminates a cuda stream
     use iso_c_binding
     implicit none
     type(c_ptr), value ::stream ! streamid
   end function cudaStreamDestroy

   ! ==================================================================== !
   ! ======================= NCCL API INTERFACES ======================== !
   ! ==================================================================== !

   ! ncclGetErrorString
   character(c_char) function ncclGetErrorString(ierr) &
    &              bind (C, name="ncclGetErrorString" )
     ! getting error string
     use iso_c_binding
     implicit none
     integer(c_int), value :: ierr
   end function ncclGetErrorString

   ! ncclGetUniqueId
   integer (c_int) function ncclGetUniqueId(uid) &
    &              bind (C, name="ncclGetUniqueId" )
     ! get the unique id for a device
     use iso_c_binding
     import ncclUniqueId
     implicit none
     type(ncclUniqueID)          :: uid ! address of
   end function ncclGetUniqueId

   ! ncclCommInitRank
   integer (c_int) function ncclCommInitRank(comm, nDevice, uid, rank) &
    &              bind (C, name="ncclCommInitRank" )
     ! initialize the Communicator for a certain rank
     use iso_c_binding
     import ncclUniqueId
     import ncclComm
     implicit none
     type(ncclComm)         :: comm ! address of
     ! type(c_ptr)                :: commPtr ! address of
     type(ncclUniqueId),value   :: uid
     integer(c_int), value  :: nDevice
     integer(c_int), value  :: rank
   end function ncclCommInitRank
   
   ! ncclCommInitAll
   integer (c_int) function ncclCommInitAll(comm, nDevice, devs) &
    &              bind (C, name="ncclCommInitAll" )
     ! initialize the Communicator for all devices
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm)         :: comm
     integer(c_int), value  :: nDevice
     integer(c_int)         :: devs(nDevice)
   end function ncclCommInitAll

   ! ncclCommUserRank
   integer (c_int) function ncclCommUserRank(comm, rank) &
    &              bind (C, name="ncclCommUserRank" )
     ! initialize the Communicator for a certain rank
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm), value  :: comm
     type(c_ptr)            :: rank
   end function ncclCommUserRank

   ! ncclCommCount
   integer (c_int) function ncclCommCount(comm, count) &
    &              bind (C, name="ncclCommCount" )
     ! initialize the Communicator for a certain rank
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm), value  :: comm
     type(c_ptr)            :: count
   end function ncclCommCount

   integer(c_int) function ncclSend(sendbuff, count, datatype, &
    &              peer, comm, stream)  bind (C, name="ncclSend" )
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm),value      :: comm
     ! type(c_ptr),value   :: comm
     ! return the device assoicated with the come
     type (c_ptr), value       :: sendbuff,  stream! pointer
     integer (c_size_t), value :: count
     integer (c_int), value :: datatype,  peer
   end function ncclSend
   
   ! ncclRecv
   integer(c_int) function ncclRecv(recvbuff, count, datatype, &
    &              peer, comm, stream)  bind (C, name="ncclRecv" )
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm),value      :: comm
     ! type(c_ptr),value   :: comm
     ! return the device assoicated with the come
     type (c_ptr), value       :: recvbuff,  stream! pointer
     integer (c_size_t), value :: count
     integer (c_int), value    :: datatype,  peer
   end function ncclRecv


   ! ncclAllReduce
   integer(c_int) function ncclAllReduce(sendbuff, recvbuff, count, datatype, &
    &              comm, stream)  bind (C, name="ncclAllReduce" )
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm), value     :: comm
     ! type(c_ptr),value   :: comm
     ! return the device assoicated with the comm
     type (c_ptr), value       :: sendbuff, recvbuff, stream! pointer
     integer (c_size_t), value :: count
     integer (c_int), value    ::  datatype
   end function ncclAllReduce

   ! ncclAllGather
   integer(c_int) function ncclAllGather(sendbuff, recvbuff, count, datatype, &
    &              comm, stream)  bind (C, name="ncclAllGather" )
     use iso_c_binding
     import ncclComm
     implicit none
     ! return the device assoicated with the comm
     type (c_ptr), value       :: sendbuff, recvbuff
     integer (c_size_t), value :: count
     integer (c_int), value    :: datatype
     type(ncclComm),value      :: comm
     type (c_ptr), value       :: stream! pointer
   end function ncclAllGather

   ! ncclAllGatherV
   integer(c_int) function ncclAllGatherV(sendbuff, sendcount, datatype, &
    &              recvbuff, recvcounts, displs, root, n, comm, stream)  &
    &              bind (C, name="ncclAllGatherV" )
     use iso_c_binding
     import ncclComm
     implicit none
     ! return the device assoicated with the comm
     type (c_ptr), value       :: sendbuff, recvbuff
     integer (c_size_t), value :: sendcount
     integer (c_int), value    :: datatype, root, n
     integer (c_size_t)        :: recvcounts(n), displs(n)
     type(ncclComm),value      :: comm
     type (c_ptr), value       :: stream! pointer
   end function ncclAllGatherV

   ! ncclCommFinalize
   integer (c_int) function ncclCommFinalize(comm) &
    &              bind (C, name="ncclCommFinalize" )
     ! destroy the Communicator
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm),value      :: comm
     ! type(c_ptr),value          :: comm
   end function ncclCommFinalize

   ! ncclCommDestroy
   integer (c_int) function ncclCommDestroy(comm) &
    &              bind (C, name="ncclCommDestroy" )
     ! destroy the Communicator
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm),value      :: comm
   end function ncclCommDestroy

   ! ==================================================================== !
   ! ===================== CUDA EVENT INTERFACES ======================== !
   ! ==================================================================== !
   ! not really useful for computation, needed for performance assessment

   ! cudaEventCreate
   integer(c_int) function cudaEventCreate(event) &
    &            bind(C,name="cudaEventCreate")
    ! this creates a cuda event 
    use iso_c_binding
    implicit none
    type(c_ptr)           :: event ! event ID
   end function cudaEventCreate

   ! cudaEventCreateWithFlags
   integer(c_int) function cudaEventCreateWithFlags(event, flag) &
    &            bind(C,name="cudaEventCreateWithFlags")
    ! this creates a cuda event with given flags
    use iso_c_binding
    implicit none
    type(c_ptr)           :: event ! event ID
    integer(c_int), value :: flag  ! flag 
    ! could be cudaEventDefault, cudaEventBlockingSync
    ! cudaEventDisableTiming
   end function cudaEventCreateWithFlags

   ! cudaEventRecord
   integer(c_int) function cudaEventRecord(event, stream) &
    &            bind(C,name="cudaEventRecord")
    ! this records a cuda event (in a given stream)
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: event ! event ID
    type(c_ptr), value    :: stream  ! cuda stream ID
   end function cudaEventRecord

   ! cudaEventElapsedTime
   integer(c_int) function cudaEventElapsedTime(time, starts, ends) &
    &            bind(C,name="cudaEventElapsedTime")
    ! this calculates the time betweent two recoreded cuda events
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: time  ! (pointer) in milliseconds
    type(c_ptr), value    :: starts! starting event ID
    type(c_ptr), value    :: ends  ! ending event ID
   end function cudaEventElapsedTime

   ! cudaEventQuery
   integer(c_int) function cudaEventQuery(event) &
    &            bind(C,name="cudaEventQuery")
    ! this queries the status of the device preceding a cuda event 
    ! that is recorded by cudaEventRecord
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: event ! event ID
   end function cudaEventQuery

   ! cudaEventDestroy
   integer(c_int) function cudaEventDestroy(event) &
    &            bind(C,name="cudaEventDestroy")
    ! this distroys the specified cuda event 
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: event ! event ID
   end function cudaEventDestroy

   ! ==================================================================== !
   ! ======================= CUBLAS INTERFACES ========================== !
   ! ==================================================================== !
   ! cublasCreate - need to be called to initialize cublas
   integer(c_int) function cublasCreate(handle) &
    &            bind(C,name="cublasCreate_v2")
     use iso_c_binding
     implicit none
     type(c_ptr):: handle
   end function cublasCreate

   ! cublasDestroy - need to be called after cublas ends
   integer(c_int) function cublasDestroy(handle) &
    &            bind(C,name="cublasDestroy_v2")
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
   end function cublasDestroy

   ! cublasSetStream
   integer(c_int) function cublasSetStream(handle,stream) &
    &            bind(C,name="cublasSetStream_v2")
     ! this sets the stream to be used by the cublas lib
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     type(c_ptr),value::stream
   end function cublasSetStream

   ! cublasGetStream
   integer(c_int) function cublasGetStream(handle) &
    &            bind(C,name="cublasGetStream")
     ! this gets the stream to be used by the cublas lib
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
   end function cublasGetStream

   ! cublasGetVector
   integer(c_int) function cublasGetVector(n, elemSize, x, incx, y, incy) &
    &            bind(C,name="cublasGetVector")
     ! this copies n elements from a vector x in GPU memory space to a
     ! vector y in the host memory
     ! each element costs elemSize in the memory
     ! note this is column based if the x/y are matrices
     use iso_c_binding
     implicit none
     integer(c_int), value :: n, elemSize
     type(c_ptr),    value :: x, y
     integer(c_int), value :: incx, incy !strides between elements of x/y
   end function cublasGetVector

   ! cublasSetVector
   integer(c_int) function cublasSetVector(n, elemSize, x, incx, y, incy) &
    &            bind(C,name="cublasSetVector")
     ! this copies n elements from a vector x in host memory space to a
     ! vector y in the GPU memory
     ! each element costs elemSize in the memory
     ! note this is column based if the x/y are matrices
     use iso_c_binding
     implicit none
     integer(c_int), value :: n, elemSize
     type(c_ptr),    value :: x, y
     integer(c_int), value :: incx, incy !strides between elements of x/y
   end function cublasSetVector

   ! cublasSetMatrix
   integer(c_int) function cublasSetMatrix(nrow, ncol, elemSize, &
    &            A, lda, B, ldb) &
    &            bind(C,name="cublasSetMatrix")
     ! this copies nrow by ncol  elements from a matrix A in host memory 
     ! to a matrix B in the GPU memory
     ! each element costs elemSize in the memory
     use iso_c_binding
     implicit none
     integer(c_int), value :: nrow, ncol, elemSize 
     type(c_ptr),    value :: A, B ! A-> src B-> dest
     ! the lda and ldb are the size of the leading dimension (rows) of A/B 
     integer(c_int)       :: lda, ldb
   end function cublasSetMatrix

   ! cublasGetMatrix
   integer(c_int) function cublasGetMatrix(nrow, ncol, elemSize, &
    &            A, lda, B, ldb) &
    &            bind(C,name="cublasGetMatrix")
     ! this copies nrow by ncol  elements from a matrix A in GPU memory 
     ! to a matrix B in the host memory
     ! each element costs elemSize in the memory
     use iso_c_binding
     implicit none
     integer(c_int), value :: nrow, ncol, elemSize 
     type(c_ptr),    value :: A, B ! A-> src B-> dest
     ! the lda and ldb are the size of the leading dimension (rows) of A/B 
     integer(c_int)        :: lda, ldb
   end function cublasGetMatrix

   ! cublasDaxpy
   integer(c_int) function cublasDaxpy(handle,n,alpha,x,incx,y,incy) &
    &            bind(C,name="cublasDaxpy_v2")
     ! compute y = y + a*x with double precision
     ! note that x and y should be located in GPU memory
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! cublas context
     integer(c_int),value  :: n   ! n --> size of x/y
     real(c_double)        ::alpha        ! scaler a
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y in/out
     integer(c_int),value  :: incy ! strides between elements of x
   end function cublasDaxpy

   ! cublasZaxpy
   integer(c_int) function cublasZaxpy(handle,n,alpha,x,incx,y,incy) &
    &             bind(C,name="cublasZaxpy_v2")
     ! compute y = y + a*x with complex double precision
     ! note that x and y should be located in GPU memory
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! cublas context
     integer(c_int),value  :: n   ! n --> size of x/y
     complex(c_double)::alpha        ! complex scaler a
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y in/out
     integer(c_int),value  :: incy ! strides between elements of x
   end function cublasZaxpy

   ! cublasDcopy
   integer(c_int) function cublasDcopy(handle,n,x,incx,y,incy) &
    &             bind(C,name="cublasDcopy_v2")
     ! compute y = x with double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! cublas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr), value    :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
   end function cublasDcopy

   ! cublasZaxpy
   integer(c_int) function cublasZcopy(handle,n,x,incx,y,incy) &
    &             bind(C,name="cublasZcopy_v2")
     ! compute y = x with complex double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! cublas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
   end function cublasZcopy

   ! cublasDdot
   integer(c_int) function cublasDdot(handle,n,x,incx,y,incy,result) &
    &            bind(C,name="cublasDdot_v2")
     ! compute result = x dot y with double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! cublas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
     type(c_ptr),value     :: result  ! output result (can be in host mem)
   end function cublasDdot

   ! cublasZdot
   integer(c_int) function cublasZdot(handle,n,x,incx,y,incy,result) &
    &            bind(C,name="cublasZdotc_v2")
     ! compute result = x dot y with complex double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! cublas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
     type(c_ptr),value     :: result  ! output result (can be in host mem)
   end function cublasZdot

   ! cublasDnrm2
   integer(c_int) function cublasDnrm2(handle,n,x,incx,norm) &
    &            bind(C,name="cublasDnrm2_v2")
     ! compute result = norm(x) in double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! cublas context
     integer(c_int),value  :: n   ! n --> size of x
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: norm    ! output result (can be in host mem)
   end function cublasDnrm2

   ! cublasZnrm2 there is no such thing like znrm2!
   integer(c_int) function cublasZnrm2(handle,n,x,incx,norm) &
    &            bind(C,name="cublasDznrm2_v2")
     ! compute result = norm(x) in complex double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! cublas context
     integer(c_int),value  :: n   ! n --> size of x
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: norm    ! output result (can be in host mem)
   end function cublasZnrm2


   ! cublasDscal
   integer(c_int) function cublasDscal(handle,n,alpha,x,incx) &
    &            bind(C,name="cublasDscal_v2")
     ! this scales the vector x by the scalar alpha 
     ! x = x/alpha in double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! cublas context
     integer(c_int), value :: n   ! n --> size of x
     real(c_double)        :: alpha   ! double scaler alpha
     type(c_ptr),value     :: x       ! vector x (in/out)
     integer(c_int),value  :: incx ! strides between elements of x
   end function cublasDscal

   ! cublasZscal
   integer(c_int) function cublasZscal(handle,n,alpha,x,incx) &
    &            bind(C,name="cublasZscal_v2")
     ! this scales the vector x by the scalar alpha 
     ! x = x/alpha in double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! cublas context
     integer(c_int), value :: n   ! n --> size of x
     complex(c_double):: alpha   ! complex double scaler alpha
     type(c_ptr),value     :: x       ! vector x (in/out)
     integer(c_int),value  :: incx ! strides between elements of x
   end function cublasZscal

   ! ==================================================================== !
   ! ====================== CUSPARSE INTERFACES ========================= !
   ! ==================================================================== !

   ! cusparseCreate - need to be called to initialize cusparse
   integer(c_int) function cusparseCreate(cusparseHandle) & 
    &            bind(C,name="cusparseCreate")
     use iso_c_binding
     implicit none
     type(c_ptr)::cusparseHandle  ! cusparse context
   end function cusparseCreate

   ! cusparseDestroy - need to be called after the cusparse operation ends
   integer(c_int) function cusparseDestroy(cusparseHandle) &
    &            bind(C,name="cusparseDestroy")
     use iso_c_binding
     implicit none
     type(c_ptr),value::cusparseHandle   ! cusparse context
   end function cusparseDestroy

   ! cusparseGetStream
   integer(c_int) function cusparseGetStream(cusparseHandle,stream) &
    &            bind(C,name="cusparseGetStream")
     use iso_c_binding
     implicit none
     type(c_ptr),value :: cusparseHandle !handle to cuSPARSE context
     type(c_ptr),value :: stream !out, NULL if the stream is not set
   end function cusparseGetStream

   ! cusparseSetStream
   integer(c_int) function cusparseSetStream(cusparseHandle,stream) &
    &            bind(C,name="cusparseSetStream")
     ! this sets the stream to be used by the cusparse lib
     use iso_c_binding
     implicit none
     type(c_ptr),value :: cusparseHandle !handle to cuSPARSE context
     type(c_ptr),value :: stream !in, stream id 
   end function cusparseSetStream

   ! cusparseCreateMatDescr
   integer(c_int) function cusparseCreateMatDescr(descrA) &
    &            bind(C,name="cusparseCreateMatDescr")
     ! this creates an empty sparse mat with descrA
     use iso_c_binding
     implicit none
     type(c_ptr):: descrA ! the matrix descr
   end function cusparseCreateMatDescr

   ! cusparseDestroyMatDescr
   integer(c_int) function cusparseDestroyMatDescr(descrA) &
    &            bind(C,name="cusparseDestroyMatDescr")
     ! this distroys the mat descrA 
     use iso_c_binding
     implicit none
     type(c_ptr), value :: descrA ! the matrix descr
   end function cusparseDestroyMatDescr

   ! cusparseDestroySpMat
   integer(c_int) function cusparseDestroySpMat(SpMatA) &
    &            bind(C,name="cusparseDestroySpMat")
     ! this distroys the mat descrA and releases the memory
     ! note this is different from the cusparseDestroyMatDescr...
     use iso_c_binding
     implicit none
     type(c_ptr), value :: SpMatA ! the matrix descr
   end function cusparseDestroySpMat

   ! cusparseGetMatType
   integer(c_int) function cusparseGetMatType(descrA) &
    &            bind(C,name="cusparseGetMatType")
     ! this gets the MatType from the mat descrA
     ! the type can be:
     ! CUSPARSE_MATRIX_TYPE_GENERAL
     ! CUSPARSE_MATRIX_TYPE_SYMMETRIC
     ! CUSPARSE_MATRIX_TYPE_HERMITIAN
     ! CUSPARSE_MATRIX_TYPE_TRIANGULAR
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function cusparseGetMatType

   ! cusparseSetMatType
   integer(c_int) function cusparseSetMatType(descrA, type) &
    &            bind(C,name="cusparseGetMatType")
     ! this sets the MatType for the mat descrA
     ! the type can be:
     ! CUSPARSE_MATRIX_TYPE_GENERAL
     ! CUSPARSE_MATRIX_TYPE_SYMMETRIC
     ! CUSPARSE_MATRIX_TYPE_HERMITIAN
     ! CUSPARSE_MATRIX_TYPE_TRIANGULAR
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int), value:: type !in
   end function cusparseSetMatType

   ! cusparseGetMatIndexBase
   integer(c_int) function cusparseGetMatIndexBase(descrA)             &
    &            bind(C,name="cusparseGetMatIndexBase")
     ! this gets the index base for the mat descrA
     ! could be zero (c convention) or one (fortran convention)
     ! CUSPARSE_INDEX_BASE_ZERO
     ! CUSPARSE_INDEX_BASE_ONE
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function cusparseGetMatIndexBase

   ! cusparseSetMatIndexBase
   integer(c_int) function cusparseSetMatIndexBase(descrA,idxbase)      &
    &            bind(C,name="cusparseSetMatIndexBase")
     ! this sets the index base for the mat descrA
     ! could be zero (c convention) or one (fortran convention)
     ! CUSPARSE_INDEX_BASE_ZERO
     ! CUSPARSE_INDEX_BASE_ONE
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int), value:: idxbase! the matrix index base
   end function cusparseSetMatIndexBase


   ! cusparseGetMatFillMode
   integer(c_int) function cusparseGetMatFillMode(descrA) &
    &            bind(C,name="cusparseGetMatFillMode")
     ! this gets the fill mode for the mat descrA
     ! could be 
     ! CUSPARSE_FILL_MODE_LOWER
     ! CUSPARSE_FILL_MODE_UPPER
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function cusparseGetMatFillMode

   ! cusparseSetMatFillMode
   integer(c_int) function cusparseSetMatFillMode(descrA,fillmode) &
    &            bind(C,name="cusparseSetMatFillMode")
     ! this sets the fill mode for the mat descrA
     ! could be 
     ! CUSPARSE_FILL_MODE_LOWER
     ! CUSPARSE_FILL_MODE_UPPER
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int),value :: fillmode !in
   end function cusparseSetMatFillMode

   ! cusparseGetMatDiagType
   integer(c_int) function cusparseGetMatDiagType(descrA) &
    &            bind(C,name="cusparseGetMatDiagType")
     ! this gets the Mat diag type for the mat descrA
     ! could be 
     ! CUSPARSE_DIAG_TYPE_NON_UNIT
     ! CUSPARSE_DIAG_TYPE_UNIT
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function cusparseGetMatDiagType

   ! cusparseGetMatDiagType
   integer(c_int) function cusparseSetMatDiagType(descrA,diagtype) &
    &            bind(C,name="cusparseSetMatDiagType")
     ! this sets the Mat diag type for the mat descrA
     ! could be 
     ! CUSPARSE_DIAG_TYPE_NON_UNIT
     ! CUSPARSE_DIAG_TYPE_UNIT
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int),value :: diagtype ! in
   end function cusparseSetMatDiagType

   ! cusparseCreateCsr
   integer(c_int) function cusparseCreateCsr(SpMatA, &
    &           nRow,nCol,nnz,row,col,val,&
    &           rowIndType, ColIndType, idxbase, valueType) &
    &           bind(C,name="cusparseCreateCsr")
     ! this creates a CSR sparse matrix with given parameters
     use iso_c_binding
     implicit none
     type(c_ptr)         ::SpMatA ! the SP matrix descr (output)
     integer(c_int),value::nRow  ! number of rows 
     integer(c_int),value::nCol  ! number of cols
     integer(c_int),value::nnz   ! number of none zero elements
     type(c_ptr),value::row      ! row indices, size of rows+1 
     type(c_ptr),value::col      ! col indices, size of nnz 
     type(c_ptr),value::val      ! values, size of nnz
     ! could be zero (c convention) or one (fortran convention)
     ! CUSPARSE_INDEX_BASE_ZERO
     ! CUSPARSE_INDEX_BASE_ONE
     integer(c_int),value::rowIndType  ! data type of row
     integer(c_int),value::colIndType  ! data type of col
     integer(c_int),value::idxbase     ! base index of col and row 
     integer(c_int),value::valueType   ! datatype of val
   end function cusparseCreateCsr

   ! cusparseCsrGet
   integer(c_int) function cusparseCsrGet(SpMatA, &
    &           nRow,nCol,nnz,row,col,val,&
    &           rowIndType, ColIndType, idxbase, valueType) &
    &           bind(C,name="cusparseCsrGet")
     ! this gets all the parameters of a CSR sparse matrix 
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     integer(c_int),value::nRow  ! number of rows 
     integer(c_int),value::nCol  ! number of cols
     integer(c_int),value::nnz   ! number of none zero elements
     type(c_ptr),value::row      ! row indices, size of rows+1 
     type(c_ptr),value::col      ! col indices, size of nnz 
     type(c_ptr),value::val      ! values, size of nnz
     ! could be zero (c convention) or one (fortran convention)
     ! CUSPARSE_INDEX_BASE_ZERO
     ! CUSPARSE_INDEX_BASE_ONE
     integer(c_int),value::rowIndType  ! data type of row
     integer(c_int),value::colIndType  ! data type of col
     integer(c_int),value::idxbase     ! base index of col and row 
     integer(c_int),value::valueType   ! datatype of val
   end function cusparseCsrGet
   
   ! cusparseCreateDnVec
   integer(c_int) function cusparseCreateDnVec(vec, n, &
    &           val, valueType) bind(C,name="cusparseCreateDnVec")
     ! this creates a Vec datatype used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr)         ::vec    ! the vec descr (out)
     integer(c_int),value::n      ! size of the vec (in)
     type(c_ptr),value   ::val       ! values, size of n (in) on device
     integer(c_int),value::valueType   ! datatype of val (in)
   end function cusparseCreateDnVec

   ! cusparseDestroyDnVec
   integer(c_int) function cusparseDestroyDnVec(vec) &
    &             bind(C,name="cusparseDestroyDnVec")
     ! this destroys a Vec datatype used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr), value  ::  vec    ! the vec descr (in)
   end function cusparseDestroyDnVec

   ! cusparseDnVecGet
   integer(c_int) function cusparseDnVecGet(vec, n, &
    &           val, valueType) bind(C,name="cusparseDnVecGet")
     ! this gets all the field of the Vec datatype used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr)   ,value::vec    ! the vec descr (in)
     integer(c_int),value::n      ! size of the vec (out)
     type(c_ptr)         ::val       ! values, size of n (out)
     integer(c_int),value::valueType   ! datatype of val (out)
   end function cusparseDnVecGet

   ! cusparseDnVecGetValues
   integer(c_int) function cusparseDnVecGetValues(vec,  &
    &           val) bind(C,name="cusparseDnVecGetValues")
     ! this gets the Vec datatype values used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr), value::vec    ! the vec descr (in)
     type(c_ptr)       ::val    ! values, size of n (out)
   end function cusparseDnVecGetValues

   ! cusparseDnVecSetValues
   integer(c_int) function cusparseDnVecSetValues(vec,  &
    &           val) bind(C,name="cusparseDnVecSetValues")
     ! this sets the Vec datatype values used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr), value :: vec    ! the vec descr (in)
     type(c_ptr), value :: val    ! values, size of n (in)
   end function cusparseDnVecSetValues


   ! cusparseSpMatGetSize
   integer(c_int) function cusparseSpMatGetSize(SpMatA, &
    &           nRow,nCol,nnz)  bind(C,name="cusparseSpMatGetSize")
     ! this gets the size of a sparse matrix 
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     integer(c_int),value::nRow  ! number of rows 
     integer(c_int),value::nCol  ! number of cols
     integer(c_int),value::nnz   ! number of none zero elements
   end function cusparseSpMatGetSize

   ! cusparseSpMatGetValues
   integer(c_int) function cusparseSpMatGetValues(SpMatA, &
    &           val) bind(C,name="cusparseSpMatGetValues")
     ! this gets the val of a sparse matrix
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     type(c_ptr),value::val      ! values, size of nnz (output)
   end function cusparseSpMatGetValues

   ! cusparseSpMatSetValues
   integer(c_int) function cusparseSpMatSetValues(SpMatA, &
    &           val) bind(C,name="cusparseSpMatSetValues")
     ! this gets the val of a sparse matrix
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     type(c_ptr),value::val      ! values, size of nnz (input)
   end function cusparseSpMatSetValues

   ! cusparseSpMatSetAttribute
   integer(c_int) function cusparseSpMatSetAttribute(SpMatA, attribute, &
    &           data,dataSize) bind(C,name="cusparseSpMatSetAttribute")
     ! this gets the val of a sparse matrix
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     integer(c_int), value:: attribute! parameter type
     ! can be
     ! CUSPARSE_SPMAT_FILL_MODE or CUSPARSE_SPMAT_DIAG_TYPE
     integer (c_int)        :: data     ! parameter value
     ! can be 
     ! CUSPARSE_FILL_MODE_LOWER   CUSPARSE_FILL_MODE_UPPER
     ! CUSPARSE_DIAG_TYPE_NON_UNIT   CUSPARSE_DIAG_TYPE_UNIT
     integer (c_int), value :: dataSize ! parameter value size
   end function cusparseSpMatSetAttribute

   ! cusparseSpMV_bufferSize
   integer(c_int) function cusparseSpMV_bufferSize(handle,opA,    &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, bufferSize) bind(C,name="cusparseSpMV_bufferSize")
     ! this returns the buffersize needed for matrix-vector multiplition
     ! note that we need to manually allocate the pBuffer 
     ! according to the size
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     real(c_double)    ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     real(c_double)    ::beta ! scaler b
     type(c_ptr), value:: vecY ! vector Y (input)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function cusparseSpMV_bufferSize

   ! cusparseSpMV_bufferSize_c
   integer(c_int) function cusparseSpMV_bufferSize_cmplx(handle,opA, &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, bufferSize) bind(C,name="cusparseSpMV_bufferSize")
     ! this returns the buffersize needed for matrix-vector multiplition
     ! note that we need to manually allocate the pBuffer 
     ! according to the size
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_double)::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     complex(c_double)::beta ! scaler b
     type(c_ptr), value:: vecY ! vector Y (input)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function cusparseSpMV_bufferSize_cmplx
   
   ! cusparseSpMV
   integer(c_int) function cusparseSpMV(handle,opA, &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, pBuffer) bind(C,name="cusparseSpMV")
     ! this calculates the matrix-vector multiplition
     ! y = alpha*A*x + beta*y (this is called Axpy in most libs)
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by cusparseSpMV_bufferSize
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     real(c_double)    ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     real(c_double)    ::beta ! scaler b
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value   :: pBuffer
   end function cusparseSpMV

   ! cusparseSpMV_cmplx
   integer(c_int) function cusparseSpMV_cmplx(handle,opA, &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, pBuffer) bind(C,name="cusparseSpMV")
     ! this calculates the matrix-vector multiplition
     ! y = alpha*A*x + beta*y (this is called Axpy in most libs)
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by cusparseSpMV_bufferSize
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_double)::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     complex(c_double)::beta ! scaler b
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value   :: pBuffer
   end function cusparseSpMV_cmplx

   ! cusparseSpSV_createDescr
   integer(c_int) function cusparseSpSV_createDescr(spsvDescr  &
    &            ) bind(C,name="cusparseSpSV_createDescr")
     ! this creates a handle for the triangle solver cusparseSpSV
     use iso_c_binding
     implicit none
     type(c_ptr):: spsvDescr ! spsv context handle
   end function cusparseSpSV_createDescr

   ! cusparseSpSV_destroyDescr
   integer(c_int) function cusparseSpSV_destroyDescr(spsvDescr &
    &            ) bind(C,name="cusparseSpSV_destroyDescr")
     ! this destroys a handle for the triangle solver cusparseSpSV
     ! and releases the memory associated
     use iso_c_binding
     implicit none
     type(c_ptr), value :: spsvDescr ! spsv context handle
   end function cusparseSpSV_destroyDescr

   ! cusparseSpSV_buffersize
   integer(c_int) function cusparseSpSV_bufferSize(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            bufferSize) bind(C,name="cusparseSpSV_bufferSize")
     ! this computes the buffersize needed by the triangle solver 
     ! cusparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     real(c_double)    ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function cusparseSpSV_bufferSize

   ! cusparseSpSV_buffersize_cmplx
   integer(c_int) function cusparseSpSV_bufferSize_cmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            bufferSize) bind(C,name="cusparseSpSV_bufferSize")
     ! this computes the buffersize needed by the triangle solver cusparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size
     use iso_c_binding
     ! complex double version
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_float) ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function cusparseSpSV_bufferSize_cmplx

   ! cusparseSpSV_buffersize_dcmplx
   integer(c_int) function cusparseSpSV_bufferSize_dcmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            bufferSize) bind(C,name="cusparseSpSV_bufferSize")
     ! this computes the buffersize needed by the triangle solver cusparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size
     use iso_c_binding
     ! complex double version
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_double) ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function cusparseSpSV_bufferSize_dcmplx

   ! cusparseSpSV_analysis
   integer(c_int) function cusparseSpSV_analysis(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            pBuffer ) bind(C,name="cusparseSpSV_analysis")
     ! this does the analysis phase needed by the triangle solver cusparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by cusparseSpMV_bufferSize
     ! double version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     real(c_double)      ::alpha ! scaler a
     type(c_ptr),value   :: matA ! the matrix descr (input)
     type(c_ptr),value   :: vecX ! vector X (input)
     type(c_ptr),value   :: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
     type(c_ptr),value:: pBuffer ! external buffer
   end function cusparseSpSV_analysis

   ! cusparseSpSV_analysis_cmplx
   integer(c_int) function cusparseSpSV_analysis_cmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            pBuffer ) bind(C,name="cusparseSpSV_analysis")
     ! this does the analysis phase needed by the triangle solver cusparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by cusparseSpMV_bufferSize
     ! complex double version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_float) ::alpha ! scaler a
     type(c_ptr),value   :: matA ! the matrix descr (input)
     type(c_ptr),value   :: vecX ! vector X (input)
     type(c_ptr),value   :: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
     type(c_ptr),value:: pBuffer ! external buffer
   end function cusparseSpSV_analysis_cmplx

   ! cusparseSpSV_analysis_dcmplx
   integer(c_int) function cusparseSpSV_analysis_dcmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            pBuffer ) bind(C,name="cusparseSpSV_analysis")
     ! this does the analysis phase needed by the triangle solver cusparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by cusparseSpMV_bufferSize
     ! complex float version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_double)   ::alpha ! scaler a
     type(c_ptr),value   :: matA ! the matrix descr (input)
     type(c_ptr),value   :: vecX ! vector X (input)
     type(c_ptr),value   :: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
     type(c_ptr),value:: pBuffer ! external buffer
   end function cusparseSpSV_analysis_dcmplx

   ! cusparseSpSV_solve
   integer(c_int) function cusparseSpSV_solve(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, &
    &            spsvDescr) bind(C,name="cusparseSpSV_solve")
     ! this solves the triangle system 
     ! y = Lx or y = Ux
     ! double version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     real(c_double)      ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
   end function cusparseSpSV_solve

   ! cusparseSpSV_solve_cmplx
   integer(c_int) function cusparseSpSV_solve_cmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, &
    &            spsvDescr) bind(C,name="cusparseSpSV_solve")
     ! this solves the triangle system 
     ! y = Lx or y = Ux
     ! complex float version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_float) ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
   end function cusparseSpSV_solve_cmplx

   ! cusparseSpSV_solve_dcmplx
   integer(c_int) function cusparseSpSV_solve_dcmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, &
    &            spsvDescr) bind(C,name="cusparseSpSV_solve")
     ! this solves the triangle system 
     ! y = Lx or y = Ux 
     ! complex double version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_double) ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
   end function cusparseSpSV_solve_dcmplx
   
   
   ! cusparseCreateCsrilu0Info
   integer(c_int) function cusparseCreateCsrilu0Info(info) &
    &           bind(C,name="cusparseCreateCsrilu02Info")
     ! this create the info pointer for the ilu0 context
     use iso_c_binding
     implicit none
     type(c_ptr)::info
   end function cusparseCreateCsrilu0Info

   ! cusparseDestroyCsrilu0Info
   integer(c_int) function cusparseDestroyCsrilu0Info(info) &
    &           bind(C,name="cusparseDestroyCsrilu02Info")
     ! this destroys the info pointer for the ilu0 context
     ! and frees its memory
     use iso_c_binding
     implicit none
     type(c_ptr), value ::info
   end function cusparseDestroyCsrilu0Info

   ! cusparseDcsrilu0_bufferSize
   integer(c_int) function cusparseDcsrilu0_bufferSize(handle,m,nnz,descrA,&
    &           valA,rowA,colA,info,bufferSize) &
    &           bind(C,name="cusparseDcsrilu02_bufferSize")
     ! this calculate the buffersize needed for ilu0 operation 
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::m
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr),value::info
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function cusparseDcsrilu0_bufferSize

   ! cusparseZcsrilu0_bufferSize
   integer(c_int) function cusparseZcsrilu0_bufferSize(handle,m,nnz,descrA,&
    &           valA,rowA,colA,info,bufferSize) &
    &           bind(C,name="cusparseZcsrilu02_bufferSize")
     ! this calculate the buffersize needed for ilu0 operation 
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::m
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr),value::info
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function cusparseZcsrilu0_bufferSize

   ! cusparseDcsrilu0_analysis
   integer(c_int) function cusparseDcsrilu0_analysis(handle,n,nnz,descrA,&
    &           valA,rowA,colA,info,policy,pBuffer) &
    &           bind(C,name="cusparseDcsrilu02_analysis")
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::n
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr), value::info
     integer(c_int),value::policy
     type(c_ptr),value   ::pBuffer
   end function cusparseDcsrilu0_analysis

   ! cusparseZcsrilu0_analysis
   integer(c_int) function cusparseZcsrilu0_analysis(handle,n,nnz,descrA,&
    &           valA,rowA,colA,info,policy,pBuffer) &
    &           bind(C,name="cusparseZcsrilu02_analysis")
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::n
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr), value::info
     integer(c_int),value::policy
     type(c_ptr),value   ::pBuffer
   end function cusparseZcsrilu0_analysis

   ! cusparseDcsrilu0
   integer(c_int) function cusparseDcsrilu0(handle,m,nnz, &
    &           descrA,valA,rowA,colA,info, policy, pBuffer) &
    &           bind(C,name="cusparseDcsrilu02")
     ! this performs the incomplete LU factorization A = LU 
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::m
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA  ! output
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr), value::info
     integer(c_int),value::policy
     ! in - note this is the address (ptr) to the buffer
     type(c_ptr),value   ::pBuffer 
   end function cusparseDcsrilu0

   ! cusparseZcsrilu0
   integer(c_int) function cusparseZcsrilu0(handle,m,nnz, &
    &           descrA,valA,rowA,colA,info, policy, pBuffer) &
    &           bind(C,name="cusparseZcsrilu02")
     ! this performs the incomplete LU factorization A = LU 
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::m
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA  ! output
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr), value::info
     integer(c_int),value::policy
     ! in - note this is the address (ptr) to the buffer
     type(c_ptr),value   ::pBuffer 
   end function cusparseZcsrilu0

   ! cusparseDcsrilu0_zeroPivot
   integer(c_int) function cusparseXcsrilu0_zeroPivot(handle,info,loc) &
    &           bind(C,name="cusparseXcsrilu02_zeroPivot")
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     type(c_ptr),value::info
     type(c_ptr),value::loc !output, in hostmem
   end function cusparseXcsrilu0_zeroPivot

   ! =========================================================================!
   ! ====================== CUSTOM KERNEL INTERFACES =========================!
   ! =========================================================================!

   ! kernelc_s2d 
   subroutine kernelc_s2d( single, double, count, stream ) & 
    &              bind (C, name="kernelc_s2d" )
     use iso_c_binding
     implicit none
     ! force convert a chunk of memory from *single* to *double*
     type (c_ptr), value :: single, double ! pointers 
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
     type (c_ptr),value     :: stream ! cudaStream
   end subroutine kernelc_s2d

   ! kernelc_d2s 
   subroutine kernelc_d2s( double, single, count, stream ) & 
    &              bind (C, name="kernelc_d2s" )
     use iso_c_binding
     implicit none
     ! force convert a chunk of memory from *single* to *double*
     type (c_ptr),value   :: double, single ! pointers
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
     type (c_ptr),value     :: stream ! cudaStream
   end subroutine kernelc_d2s

   ! kernelc_hadar
   subroutine kernelc_hadar( a, b, c, count ) & 
    &              bind (C, name="kernelc_hadar" )
     use iso_c_binding
     implicit none
     ! perform hadamard multiply c = a*b, (real)
     type (c_ptr),value   :: a, b, c ! pointers
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
   end subroutine kernelc_hadar

   ! kernelc_hadac
   subroutine kernelc_hadac( a, b, c, count, cstream ) & 
    &              bind (C, name="kernelc_hadac" )
     use iso_c_binding
     implicit none
     ! perform hadamard multiply c = a*b, (complex)
     ! force convert two doubles into one complex double
     type (c_ptr),value     :: a, b, c ! pointers
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
     type (c_ptr),value     :: cstream ! cudaStream
   end subroutine kernelc_hadac

   ! kernelc_xpbyc
   subroutine kernelc_xpbyc( a, b, c, count , cstream) & 
    &              bind (C, name="kernelc_xpbyc" )
     use iso_c_binding
     implicit none
     ! perform y = x + by, (complex)
     ! force convert a chunk of memory from *single* to *double*
     type (c_ptr),value          :: a, c ! pointers
     complex(c_double), value    :: b ! scaler b
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value      :: count 
     type (c_ptr),value          :: cstream ! cudaStream
   end subroutine kernelc_xpbyc

   ! kernelc_update_p
   subroutine kernelc_update_pc( r, v, beta, omega, p, count, cstream ) & 
    &              bind (C, name="kernelc_update_pc" )
     use iso_c_binding
     implicit none
     ! perform y = x + by, (complex)
     ! force convert a chunk of memory from *single* to *double*
     type (c_ptr),value          :: p, r, v ! pointers
     complex(c_double), value    :: beta ! scaler beta
     complex(c_double), value    :: omega ! scaler omega
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value      :: count 
     type (c_ptr),value          :: cstream ! cudaStream
   end subroutine kernelc_update_pc

   ! kernelc_update_x
   subroutine kernelc_update_xc( ph, sh, alpha, omega, x, count ) & 
    &              bind (C, name="kernelc_update_xc" )
     use iso_c_binding
     implicit none
     ! perform y = x + by, (complex)
     ! force convert a chunk of memory from *single* to *double*
     type (c_ptr),value          :: ph, sh, x ! pointers
     complex(c_double), value    :: alpha ! scaler alpha
     complex(c_double), value    :: omega ! scaler omega
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value      :: count 
   end subroutine kernelc_update_xc

   ! cf_free
   subroutine cf_free(ptr)  bind (C, name="free" )
     use iso_c_binding
     implicit none
     ! this calls free from the c side to deal with arrays
     ! associated by c_loc
     type (c_ptr),value          :: ptr 
   end subroutine cf_free

   ! cf_resetFlags
   integer(c_int) function cf_resetFlag(device_idx) &
                   bind (C, name="cf_resetFlag" )
     ! reset the computation flag of the GPU device
     ! see the c code in kernel_c for details
     ! can be "cudaDeviceScheduleSpin" or 
     ! "cudaDeviceScheduleYield"
     use iso_c_binding
     implicit none
     integer(c_int),value ::  device_idx
     ! this calls resetflags from the c side
   end function cf_resetFlag

   ! kernelc_hookDev
   integer(c_int) function cf_hookDev(device_idx) & 
    &              bind (C, name="cf_hookDev" )
     ! bind the current context to a given device
     ! need some extra tweaks to use the "internal" functions instead of 
     ! the NVML lib - this is actually fine for CUDA, but not possible 
     ! for other APIs (AMD ROCM or Intel oneAPI) as they don't have the 
     ! equivalents of NVML
     use iso_c_binding
     implicit none
     integer(c_int),value ::  device_idx
     ! get the current device idx
   end function cf_hookDev


   ! TODO
   ! currently do not need it here, since we already defined this at 
   ! Declarition_MPI.f90
   ! since currently the GPU version is bind to MPI - I don't see 
   ! why we need it here
   ! still - we may need this at some stage when I find the time to
   ! peel GPU from MPI...
   ! kernelc_getDevNum
   ! integer(c_int) function kernelc_getDevNum() & 
   !  &              bind (C, name="kernelc_getDevNum" )
   !   use iso_c_binding
   !   implicit none
   !   ! get the number of GPU devices  
   ! end function kernelc_getDevNum

   ! cf_releaseDev - release the GPU gate so other processes can hook on
   subroutine cf_releaseDev(device_idx) &
    &              bind(C, name="cf_releaseDev")
     use iso_c_binding
     implicit none
     integer(c_int), value :: device_idx
   end subroutine cf_releaseDev

   ! cf_cleanupLock - cleanup shared memory lock on program exit
   subroutine cf_cleanupLock() &
    &              bind(C, name="cf_cleanupLock")
     use iso_c_binding
     implicit none
   end subroutine cf_cleanupLock

   end interface  

end module
