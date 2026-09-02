module hipFortMap

   ! A BAND-AID MODULE TO MAP HIP C INTERFACES TO FORTRAN
   ! WITH ISO_C_BINDING
   ! ESSENTIALLY THE HIP IS A IMATATION OF THE CUDA INTERFACE BY THE AMD
   ! COMPANY. IT DOES HAVE SOME INCOMPATIBLE PLACES HERE AND THERE -  
   ! PLACES THAT NEED TO BE MANUALLY MODIFIED, PLEASE USE WITH CAUTION
   ! BUT OVERALL YOU ARE PROBABLY FINE JUST REPLACING SOME "cuda***" with 
   ! "hip***"
   ! 
   ! Hao 
   ! 2024.05
   use iso_c_binding
   use math_constants   ! math/physics constants
   implicit none
   save
   ! ======================== HIP enumerators ======================== !
   ! note here we explicitly setup the identical interface as CUDA
   ! (the real interface is through hip libs)
   ! for the compatibility in the solver module
   ! for C these are already setup in hipart/hipblas/hipsparse headers
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
   ! policy for hipsparse solve
   integer*4         :: CUSPARSE_SOLVE_POLICY_NO_LEVEL
   integer*4         :: CUSPARSE_SOLVE_POLICY_USE_LEVEL
   parameter (CUSPARSE_SOLVE_POLICY_NO_LEVEL=0)
   parameter (CUSPARSE_SOLVE_POLICY_USE_LEVEL=1)
   ! matrix attributes for SpMat
   integer*4         :: CUSPARSE_SPMAT_FILL_MODE
   integer*4         :: CUSPARSE_SPMAT_DIAG_TYPE
   parameter (CUSPARSE_SPMAT_FILL_MODE=0)
   parameter (CUSPARSE_SPMAT_DIAG_TYPE=1)
   ! hipStream flag    
   integer*4         :: cudaStreamNonBlocking
   integer*4         :: cudaStreamDefault
   parameter (cudaStreamNonBlocking=1)
   parameter (cudaStreamDefault=0)
   ! additional C_ONE in FP32 - useful in mixed precision calculations
   complex(kind=SP)  :: C_ONEs
   parameter (C_ONEs=(1.0, 0.0))
#ifdef FG
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
#endif
   ! ======================= RCCL C derived types ===================== !
#ifdef FG
   ! ncclUniqueId
   type, bind(c) :: ncclUniqueId
       character(c_char) :: internal(NCCL_UNIQUE_ID_BYTES)
   end type
   ! rcclComm
   type, bind(c) :: ncclComm
       type(c_ptr) :: member
   end type
   ! rcclResult
   type, bind(c) :: ncclResult
       integer(c_int) :: member
   end type
#endif
   ! ========================= HIP C interfaces ========================= !
   ! I only included some "might-be-useful" interfaces here...
   ! which ranges from basic hip memory manipulation to hipblas and hipsparse
   ! feel free to add more if you need other interfaces, but keep the 
   ! original naming convention.
   ! Hao
   ! 
   ! ==================================================================== !
   ! ======================= BASIC HIP INTERFACES ======================= !
   ! ==================================================================== !
   interface
   ! hipSetDevice
   integer (c_int) function cudaSetDevice( device_idx ) &
    &              bind (C, name="hipSetDevice" ) 
     ! set the current working device according to device_idx
     use iso_c_binding
     implicit none
     integer (c_int), value  :: device_idx
   end function cudaSetDevice

   ! hipGetErrorString(err));
   type (c_ptr) function cudaGetErrorString( ierr ) &
    &       bind (C, name="hipGetErrorString" )
     ! get the number of available devices
     use iso_c_binding
     implicit none
     integer (c_int), value           :: ierr
   end function cudaGetErrorString

   ! hipGetDeviceCount
   integer (c_int) function cudaGetDeviceCount( dev_count ) &
    &              bind (C, name="hipGetDeviceCount" )
     ! get the number of available devices
     use iso_c_binding
     implicit none
     type(c_ptr), value  :: dev_count
   end function cudaGetDeviceCount

   ! hipMemset
   integer (c_int) function cudaMemset( devPtr,value, count ) &
    &              bind (C, name="hipMemset" ) 
     ! fill the count bytes of device memory pointed to by devPr with value
     ! doesn't seems quite useful unless you are filling everything with 
     ! zeros
     use iso_c_binding
     implicit none
     ! devPtr is the *device* pointer to memory 
     type (c_ptr), value  :: devPtr
     integer(c_int), value :: value
     integer(c_size_t), value :: count
   end function cudaMemset

   ! hipMalloc
   integer (c_int) function cudaMalloc ( devPtr, count ) &
    &              bind (C, name="hipMalloc" ) 
     ! allocate count bytes of memory pointed by devPtr 
     use iso_c_binding
     implicit none
     ! devPtr is the *device* pointer to memory 
     type (c_ptr)  :: devPtr
     integer (c_size_t), value :: count
   end function cudaMalloc

   ! hipFree
   integer (c_int) function cudaFree(devPtr) & 
    &              bind(C, name="hipFree")
     ! free the chunk of memory used by buffer
     use iso_c_binding
     implicit none
     ! devPtr is the *device* pointer to memory that you want to free
     type (c_ptr),value :: devPtr
   end function cudaFree

   ! hipMallocHost
   integer (c_int) function cudaMallocHost(hostPtr, count ) &
    &              bind (C, name="hipMallocHost" ) 
     ! allocate count bytes of memory pointed by devPtr 
     use iso_c_binding
     implicit none
     ! devPtr is the *host* pointer to memory 
     type (c_ptr)  :: hostPtr
     integer (c_size_t), value :: count
   end function cudaMallocHost

   ! hipFreeHost
   integer (c_int) function cudaFreeHost(hostPtr) & 
    &              bind(C, name="hipFreeHost")
     ! free the chunk of memory used by buffer
     use iso_c_binding
     implicit none
     ! devPtr is the *host* pointer to memory that you want to free
     type (c_ptr),value :: hostPtr
   end function cudaFreeHost

   ! hipMemcpy
   integer (c_int) function cudaMemcpy ( dst, src, count, kind ) & 
    &              bind (C, name="hipMemcpy" )
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

   ! hipMemcpyAsync
   integer (c_int) function cudaMemcpyAsync( dst, src, count, kind, stream ) & 
    &              bind (C, name="hipMemcpyAsync" )
     use iso_c_binding
     implicit none
     ! copy a chunk of memory from *src* to *dst*
     type (c_ptr), value       :: dst, src
     ! count is the size of memory (with unit of size_t)
     integer (c_size_t), value :: count, kind 
     type(c_ptr),value         :: stream ! streamid
     ! kind specifies the direction 
     ! cudaMemcpyHostToDevice = 1
     ! cudaMemcpyDeviceToHost = 2
     ! cudaMemcpyDeviceToDevice = 3
   end function cudaMemcpyAsync

   ! hipMemGetInfo
   integer (c_int) function cudaMemGetInfo(free, total)  &
    &     bind(C, name="hipMemGetInfo")
     use iso_c_binding
     implicit none
     type(c_ptr),value :: free  ! free memory in bytes
     type(c_ptr),value :: total ! total memory in bytes
   end function cudaMemGetInfo

   ! cudaHostRegister
   integer (c_int) function cudaHostRegister(hostPtr, size, flag)  &
    &     bind(C, name="hipHostRegister")
     use iso_c_binding
     implicit none
     type(c_ptr),value       :: hostPtr  ! host pointer
     integer(c_size_t),value :: size  ! memory in bytes
     integer, value          :: flag ! mem type 
   end function cudaHostRegister

   ! cudaHostUnregister
   integer (c_int) function cudaHostUnregister(hostPtr)  &
    &     bind(C, name="hipHostUnregister")
     use iso_c_binding
     implicit none
     type(c_ptr),value       :: hostPtr  ! host pointer
   end function cudaHostUnregister

   ! hipDeviceSynchronize
   integer(c_int) function cudaDeviceSynchronize() &
    &            bind(C,name="hipDeviceSynchronize")
     use iso_c_binding
     implicit none
   end function cudaDeviceSynchronize

   ! hipStreamCreate
   integer(c_int) function cudaStreamCreate(stream) & 
    &            bind(C,name="hipStreamCreate")
     use iso_c_binding
     implicit none
     type(c_ptr)::stream ! streamid
   end function cudaStreamCreate

   ! hipStreamCreateWithFlags
   integer(c_int) function cudaStreamCreateWithFlags(stream,flag ) & 
    &            bind(C,name="hipStreamCreate")
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
    &            bind(C,name="hipStreamWaitEvent")
    ! this let the stream to wait for a event to finish  
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: cstream ! streamid, in
    type(c_ptr), value    :: event   ! event ID, in
    integer(c_int),value  :: flag    ! in
   end function cudaStreamWaitEvent

   ! hipStreamDestroy
   integer(c_int) function cudaStreamDestroy(stream) & 
    &            bind(C,name="hipStreamDestroy")
     ! this eliminates a cuda stream
     use iso_c_binding
     implicit none
     type(c_ptr), value ::stream ! streamid
   end function cudaStreamDestroy

   ! ==================================================================== !
   ! ======================= NCCL API INTERFACES ======================== !
   ! ==================================================================== !

#ifdef FG
   ! ncclGetErrorString
   character(c_char) function ncclGetErrorString(ierr) &
    &              bind (C, name="rcclGetErrorString" )
     ! getting error string
     use iso_c_binding
     implicit none
     integer(c_int), value :: ierr
   end function ncclGetErrorString

   ! ncclGetUniqueId
   integer (c_int) function ncclGetUniqueId(uid) &
    &              bind (C, name="rcclGetUniqueId" )
     ! get the unique id for a device
     use iso_c_binding
     import ncclUniqueId
     implicit none
     type(ncclUniqueID)          :: uid ! address of
   end function ncclGetUniqueId

   ! ncclCommInitRank
   integer (c_int) function ncclCommInitRank(comm, nDevice, uid, rank) &
    &              bind (C, name="rcclCommInitRank" )
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
    &              bind (C, name="rcclCommInitAll" )
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
    &              bind (C, name="rcclCommUserRank" )
     ! initialize the Communicator for a certain rank
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm), value  :: comm
     type(c_ptr)            :: rank
   end function ncclCommUserRank

   ! ncclCommCount
   integer (c_int) function ncclCommCount(comm, count) &
    &              bind (C, name="rcclCommCount" )
     ! initialize the Communicator for a certain rank
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm), value  :: comm
     type(c_ptr)            :: count
   end function ncclCommCount

   integer(c_int) function ncclSend(sendbuff, count, datatype, &
    &              peer, comm, stream)  bind (C, name="rcclSend" )
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
    &              peer, comm, stream)  bind (C, name="rcclRecv" )
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
    &              comm, stream)  bind (C, name="rcclAllReduce" )
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
    &              comm, stream)  bind (C, name="rcclAllGather" )
     use iso_c_binding
     import ncclComm
     implicit none
     ! return the device assoicated with the comm
     type(ncclComm),value      :: comm
     type (c_ptr), value       :: sendbuff, recvbuff
     type (c_ptr), value       :: stream! pointer
     integer (c_size_t), value :: count
     integer (c_int), value    :: datatype
   end function ncclAllGather

   ! ncclAllGatherV
   integer(c_int) function ncclAllGatherV(sendbuff, sendcount, datatype, &
    &              recvbuff, recvcounts, displs, root, n, comm, stream)  &
    &              bind (C, name="rcclAllGatherV" )
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
    &              bind (C, name="rcclCommFinalize" )
     ! destroy the Communicator
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm),value      :: comm
     ! type(c_ptr),value          :: comm
   end function ncclCommFinalize

   ! ncclCommDestroy
   integer (c_int) function ncclCommDestroy(comm) &
    &              bind (C, name="rcclCommDestroy" )
     ! destroy the Communicator
     use iso_c_binding
     import ncclComm
     implicit none
     type(ncclComm),value      :: comm
   end function ncclCommDestroy
#endif

   ! ==================================================================== !
   ! ===================== CUDA EVENT INTERFACES ======================== !
   ! ==================================================================== !
   ! not really useful for computation, needed for performance assessment

   ! hipEventCreate
   integer(c_int) function cudaEventCreate(event) &
    &            bind(C,name="hipEventCreate")
    ! this creates a cuda event 
    use iso_c_binding
    implicit none
    type(c_ptr)           :: event ! event ID
   end function cudaEventCreate

   ! hipEventCreateWithFlags
   integer(c_int) function cudaEventCreateWithFlags(event, flag) &
    &            bind(C,name="hipEventCreateWithFlags")
    ! this creates a cuda event with given flags
    use iso_c_binding
    implicit none
    type(c_ptr)           :: event ! event ID
    integer(c_int), value :: flag  ! flag 
    ! could be hipEventDefault, hipEventBlockingSync
    ! hipEventDisableTiming
   end function cudaEventCreateWithFlags

   ! hipEventRecord
   integer(c_int) function cudaEventRecord(event, stream) &
    &            bind(C,name="hipEventRecord")
    ! this records a cuda event (in a given stream)
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: event ! event ID
    type(c_ptr), value    :: stream  ! cuda stream ID
   end function cudaEventRecord

   ! hipEventElapsedTime
   integer(c_int) function cudaEventElapsedTime(time, starts, ends) &
    &            bind(C,name="hipEventElapsedTime")
    ! this calculates the time betweent two recoreded cuda events
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: time  ! (pointer) in milliseconds
    type(c_ptr), value    :: starts! starting event ID
    type(c_ptr), value    :: ends  ! ending event ID
   end function cudaEventElapsedTime

   ! hipEventQuery
   integer(c_int) function cudaEventQuery(event) &
    &            bind(C,name="hipEventQuery")
    ! this queries the status of the device preceding a cuda event 
    ! that is recorded by hipEventRecord
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: event ! event ID
   end function cudaEventQuery

   ! hipEventDestroy
   integer(c_int) function cudaEventDestroy(event) &
    &            bind(C,name="hipEventDestroy")
    ! this distroys the specified cuda event 
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: event ! event ID
   end function cudaEventDestroy

   ! ==================================================================== !
   ! ====================== HIPBLAS INTERFACES ========================== !
   ! ==================================================================== !
   ! hipblasCreate - need to be called to initialize hipblas
   integer(c_int) function cublasCreate(handle) &
    &            bind(C,name="hipblasCreate")
     use iso_c_binding
     implicit none
     type(c_ptr):: handle
   end function cublasCreate

   ! hipblasDestroy - need to be called after hipblas ends
   integer(c_int) function cublasDestroy(handle) &
    &            bind(C,name="hipblasDestroy")
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
   end function cublasDestroy

   ! hipblasSetStream
   integer(c_int) function cublasSetStream(handle,stream) &
    &            bind(C,name="hipblasSetStream")
     ! this sets the stream to be used by the hipblas lib
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     type(c_ptr),value::stream
   end function cublasSetStream

   ! hipblasGetStream
   integer(c_int) function cublasGetStream(handle) &
    &            bind(C,name="hipblasGetStream")
     ! this gets the stream to be used by the hipblas lib
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
   end function cublasGetStream

   ! hipblasGetVector
   integer(c_int) function cublasGetVector(n, elemSize, x, incx, y, incy) &
    &            bind(C,name="hipblasGetVector")
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

   ! hipblasSetVector
   integer(c_int) function cublasSetVector(n, elemSize, x, incx, y, incy) &
    &            bind(C,name="hipblasSetVector")
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

   ! hipblasSetMatrix
   integer(c_int) function cublasSetMatrix(nrow, ncol, elemSize, &
    &            A, lda, B, ldb) &
    &            bind(C,name="hipblasSetMatrix")
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

   ! hipblasGetMatrix
   integer(c_int) function cublasGetMatrix(nrow, ncol, elemSize, &
    &            A, lda, B, ldb) &
    &            bind(C,name="hipblasGetMatrix")
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

   ! hipblasDaxpy
   integer(c_int) function cublasDaxpy(handle,n,alpha,x,incx,y,incy) &
    &            bind(C,name="hipblasDaxpy_v2")
     ! compute y = y + a*x with double precision
     ! note that x and y should be located in GPU memory
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     real(c_double)        ::alpha        ! scaler a
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y in/out
     integer(c_int),value  :: incy ! strides between elements of x
   end function cublasDaxpy

   ! hipblasZaxpy
   integer(c_int) function cublasZaxpy(handle,n,alpha,x,incx,y,incy) &
    &             bind(C,name="hipblasZaxpy")
     ! compute y = y + a*x with complex double precision
     ! note that x and y should be located in GPU memory
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     complex(c_double)::alpha        ! complex scaler a
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y in/out
     integer(c_int),value  :: incy ! strides between elements of x
   end function cublasZaxpy

   ! hipblasDcopy
   integer(c_int) function cublasDcopy(handle,n,x,incx,y,incy) &
    &             bind(C,name="hipblasDcopy_v2")
     ! compute y = x with double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr), value    :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
   end function cublasDcopy

   ! hipblasZcopy
   integer(c_int) function cublasZcopy(handle,n,x,incx,y,incy) &
    &             bind(C,name="hipblasZcopy")
     ! compute y = x with complex double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
   end function cublasZcopy

   ! hipblasDdot
   integer(c_int) function cublasDdot(handle,n,x,incx,y,incy,result) &
    &            bind(C,name="hipblasDdot")
     ! compute result = x dot y with double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
     type(c_ptr),value     :: result  ! output result (can be in host mem)
   end function cublasDdot

   ! hipblasZdot
   integer(c_int) function cublasZdot(handle,n,x,incx,y,incy,result) &
    &            bind(C,name="hipblasZdotc")
     ! compute result = x dot y with complex double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
     type(c_ptr),value     :: result  ! output result (can be in host mem)
   end function cublasZdot

   ! hipblasDnrm2
   integer(c_int) function cublasDnrm2(handle,n,x,incx,norm) &
    &            bind(C,name="hipblasDnrm2_v2")
     ! compute result = norm(x) in double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: norm    ! output result (can be in host mem)
   end function cublasDnrm2

   ! hipblasZnrm2 there is no such thing like znrm2!
   integer(c_int) function cublasZnrm2(handle,n,x,incx,norm) &
    &            bind(C,name="hipblasDznrm2")
     ! compute result = norm(x) in complex double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: norm    ! output result (can be in host mem)
   end function cublasZnrm2


   ! hipblasDscal
   integer(c_int) function cublasDscal(handle,n,alpha,x,incx) &
    &            bind(C,name="hipblasDscal_v2")
     ! this scales the vector x by the scalar alpha 
     ! x = x/alpha in double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int), value :: n   ! n --> size of x
     real(c_double)        :: alpha   ! double scaler alpha
     type(c_ptr),value     :: x       ! vector x (in/out)
     integer(c_int),value  :: incx ! strides between elements of x
   end function cublasDscal

   ! hipblasZscal
   integer(c_int) function cublasZscal(handle,n,alpha,x,incx) &
    &            bind(C,name="hipblasZscal_v2")
     ! this scales the vector x by the scalar alpha 
     ! x = x/alpha in double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int), value :: n   ! n --> size of x
     complex(c_double):: alpha   ! complex double scaler alpha
     type(c_ptr),value     :: x       ! vector x (in/out)
     integer(c_int),value  :: incx ! strides between elements of x
   end function cublasZscal

   ! ==================================================================== !
   ! ====================== HIPSPARSE INTERFACES ======================== !
   ! ==================================================================== !

   ! hipsparseCreate - need to be called to initialize hipsparse
   integer(c_int) function cusparseCreate(cusparseHandle) & 
    &            bind(C,name="hipsparseCreate")
     use iso_c_binding
     implicit none
     type(c_ptr)::cusparseHandle  ! hipsparse context
   end function cusparseCreate

   ! hipsparseDestroy - need to be called after the hipsparse operation ends
   integer(c_int) function cusparseDestroy(cusparseHandle) &
    &            bind(C,name="hipsparseDestroy")
     use iso_c_binding
     implicit none
     type(c_ptr),value::cusparseHandle   ! hipsparse context
   end function cusparseDestroy

   ! hipsparseGetStream
   integer(c_int) function cusparseGetStream(cusparseHandle,stream) &
    &            bind(C,name="hipsparseGetStream")
     use iso_c_binding
     implicit none
     type(c_ptr),value :: cusparseHandle !handle to cuSPARSE context
     type(c_ptr),value :: stream !out, NULL if the stream is not set
   end function cusparseGetStream

   ! hipsparseSetStream
   integer(c_int) function cusparseSetStream(cusparseHandle,stream) &
    &            bind(C,name="hipsparseSetStream")
     ! this sets the stream to be used by the hipsparse lib
     use iso_c_binding
     implicit none
     type(c_ptr),value :: cusparseHandle !handle to cuSPARSE context
     type(c_ptr),value :: stream !in, stream id 
   end function cusparseSetStream

   ! hipsparseCreateMatDescr
   integer(c_int) function cusparseCreateMatDescr(descrA) &
    &            bind(C,name="hipsparseCreateMatDescr")
     ! this creates an empty sparse mat with descrA
     use iso_c_binding
     implicit none
     type(c_ptr):: descrA ! the matrix descr
   end function cusparseCreateMatDescr

   ! hipsparseDestroyMatDescr
   integer(c_int) function cusparseDestroyMatDescr(descrA) &
    &            bind(C,name="hipsparseDestroyMatDescr")
     ! this distroys the mat descrA 
     use iso_c_binding
     implicit none
     type(c_ptr), value :: descrA ! the matrix descr
   end function cusparseDestroyMatDescr

   ! hipsparseDestroySpMat
   integer(c_int) function cusparseDestroySpMat(SpMatA) &
    &            bind(C,name="hipsparseDestroySpMat")
     ! this distroys the mat descrA and releases the memory
     ! note this is different from the hipsparseDestroyMatDescr...
     use iso_c_binding
     implicit none
     type(c_ptr), value :: SpMatA ! the matrix descr
   end function cusparseDestroySpMat

   ! hipsparseGetMatType
   integer(c_int) function cusparseGetMatType(descrA) &
    &            bind(C,name="hipsparseGetMatType")
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

   ! hipsparseSetMatType
   integer(c_int) function cusparseSetMatType(descrA, type) &
    &            bind(C,name="hipsparseGetMatType")
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

   ! hipsparseGetMatIndexBase
   integer(c_int) function cusparseGetMatIndexBase(descrA)             &
    &            bind(C,name="hipsparseGetMatIndexBase")
     ! this gets the index base for the mat descrA
     ! could be zero (c convention) or one (fortran convention)
     ! CUSPARSE_INDEX_BASE_ZERO
     ! CUSPARSE_INDEX_BASE_ONE
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function cusparseGetMatIndexBase

   ! hipsparseSetMatIndexBase
   integer(c_int) function cusparseSetMatIndexBase(descrA,idxbase)      &
    &            bind(C,name="hipsparseSetMatIndexBase")
     ! this sets the index base for the mat descrA
     ! could be zero (c convention) or one (fortran convention)
     ! CUSPARSE_INDEX_BASE_ZERO
     ! CUSPARSE_INDEX_BASE_ONE
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int), value:: idxbase! the matrix index base
   end function cusparseSetMatIndexBase


   ! hipsparseGetMatFillMode
   integer(c_int) function cusparseGetMatFillMode(descrA) &
    &            bind(C,name="hipsparseGetMatFillMode")
     ! this gets the fill mode for the mat descrA
     ! could be 
     ! CUSPARSE_FILL_MODE_LOWER
     ! CUSPARSE_FILL_MODE_UPPER
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function cusparseGetMatFillMode

   ! hipsparseSetMatFillMode
   integer(c_int) function cusparseSetMatFillMode(descrA,fillmode) &
    &            bind(C,name="hipsparseSetMatFillMode")
     ! this sets the fill mode for the mat descrA
     ! could be 
     ! CUSPARSE_FILL_MODE_LOWER
     ! CUSPARSE_FILL_MODE_UPPER
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int),value :: fillmode !in
   end function cusparseSetMatFillMode

   ! hipsparseGetMatDiagType
   integer(c_int) function cusparseGetMatDiagType(descrA) &
    &            bind(C,name="hipsparseGetMatDiagType")
     ! this gets the Mat diag type for the mat descrA
     ! could be 
     ! CUSPARSE_DIAG_TYPE_NON_UNIT
     ! CUSPARSE_DIAG_TYPE_UNIT
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function cusparseGetMatDiagType

   ! hipsparseGetMatDiagType
   integer(c_int) function cusparseSetMatDiagType(descrA,diagtype) &
    &            bind(C,name="hipsparseSetMatDiagType")
     ! this sets the Mat diag type for the mat descrA
     ! could be 
     ! CUSPARSE_DIAG_TYPE_NON_UNIT
     ! CUSPARSE_DIAG_TYPE_UNIT
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int),value :: diagtype ! in
   end function cusparseSetMatDiagType

   ! hipsparseCreateCsr
   integer(c_int) function cusparseCreateCsr(SpMatA, &
    &           nRow,nCol,nnz,row,col,val,&
    &           rowIndType, ColIndType, idxbase, valueType) &
    &           bind(C,name="hipsparseCreateCsr")
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

   ! hipsparseCsrGet
   integer(c_int) function cusparseCsrGet(SpMatA, &
    &           nRow,nCol,nnz,row,col,val,&
    &           rowIndType, ColIndType, idxbase, valueType) &
    &           bind(C,name="hipsparseCsrGet")
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
   
   ! hipsparseCreateDnVec
   integer(c_int) function cusparseCreateDnVec(vec, n, &
    &           val, valueType) bind(C,name="hipsparseCreateDnVec")
     ! this creates a Vec datatype used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr)         ::vec    ! the vec descr (out)
     integer(c_int),value::n      ! size of the vec (in)
     type(c_ptr),value   ::val       ! values, size of n (in) on device
     integer(c_int),value::valueType   ! datatype of val (in)
   end function cusparseCreateDnVec

   ! hipsparseDestroyDnVec
   integer(c_int) function cusparseDestroyDnVec(vec) &
    &             bind(C,name="hipsparseDestroyDnVec")
     ! this destroys a Vec datatype used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr), value  ::  vec    ! the vec descr (in)
   end function cusparseDestroyDnVec

   ! hipsparseDnVecGet
   integer(c_int) function cusparseDnVecGet(vec, n, &
    &           val, valueType) bind(C,name="hipsparseDnVecGet")
     ! this gets all the field of the Vec datatype used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr)   ,value::vec    ! the vec descr (in)
     integer(c_int),value::n      ! size of the vec (out)
     type(c_ptr)         ::val       ! values, size of n (out)
     integer(c_int),value::valueType   ! datatype of val (out)
   end function cusparseDnVecGet

   ! hipsparseDnVecGetValues
   integer(c_int) function cusparseDnVecGetValues(vec,  &
    &           val) bind(C,name="hipsparseDnVecGetValues")
     ! this gets the Vec datatype values used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr), value::vec    ! the vec descr (in)
     type(c_ptr)       ::val    ! values, size of n (out)
   end function cusparseDnVecGetValues

   ! hipsparseDnVecSetValues
   integer(c_int) function cusparseDnVecSetValues(vec,  &
    &           val) bind(C,name="hipsparseDnVecSetValues")
     ! this sets the Vec datatype values used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr), value :: vec    ! the vec descr (in)
     type(c_ptr), value :: val    ! values, size of n (in)
   end function cusparseDnVecSetValues


   ! hipsparseSpMatGetSize
   integer(c_int) function cusparseSpMatGetSize(SpMatA, &
    &           nRow,nCol,nnz)  bind(C,name="hipsparseSpMatGetSize")
     ! this gets the size of a sparse matrix 
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     integer(c_int),value::nRow  ! number of rows 
     integer(c_int),value::nCol  ! number of cols
     integer(c_int),value::nnz   ! number of none zero elements
   end function cusparseSpMatGetSize

   ! hipsparseSpMatGetValues
   integer(c_int) function cusparseSpMatGetValues(SpMatA, &
    &           val) bind(C,name="hipsparseSpMatGetValues")
     ! this gets the val of a sparse matrix
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     type(c_ptr),value::val      ! values, size of nnz (output)
   end function cusparseSpMatGetValues

   ! hipsparseSpMatSetValues
   integer(c_int) function cusparseSpMatSetValues(SpMatA, &
    &           val) bind(C,name="hipsparseSpMatSetValues")
     ! this gets the val of a sparse matrix
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     type(c_ptr),value::val      ! values, size of nnz (input)
   end function cusparseSpMatSetValues

   ! hipsparseSpMatSetAttribute
   integer(c_int) function cusparseSpMatSetAttribute(SpMatA, attribute, &
    &           data,dataSize) bind(C,name="hipsparseSpMatSetAttribute")
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

   ! hipsparseSpMV_bufferSize
   integer(c_int) function cusparseSpMV_bufferSize(handle,opA,    &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, bufferSize) bind(C,name="hipsparseSpMV_bufferSize")
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

   ! hipsparseSpMV_bufferSize_c
   integer(c_int) function cusparseSpMV_bufferSize_cmplx(handle,opA, &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, bufferSize) bind(C,name="hipsparseSpMV_bufferSize")
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
   
   ! hipsparseSpMV
   integer(c_int) function cusparseSpMV(handle,opA, &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, pBuffer) bind(C,name="hipsparseSpMV")
     ! this calculates the matrix-vector multiplition
     ! y = alpha*A*x + beta*y (this is called Axpy in most libs)
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by hipsparseSpMV_bufferSize
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

   ! hipsparseSpMV_cmplx
   integer(c_int) function cusparseSpMV_cmplx(handle,opA, &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, pBuffer) bind(C,name="hipsparseSpMV")
     ! this calculates the matrix-vector multiplition
     ! y = alpha*A*x + beta*y (this is called Axpy in most libs)
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by hipsparseSpMV_bufferSize
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

   ! hipsparseSpSV_createDescr
   integer(c_int) function cusparseSpSV_createDescr(spsvDescr  &
    &            ) bind(C,name="hipsparseSpSV_createDescr")
     ! this creates a handle for the triangle solver hipsparseSpSV
     use iso_c_binding
     implicit none
     type(c_ptr):: spsvDescr ! spsv context handle
   end function cusparseSpSV_createDescr

   ! hipsparseSpSV_destroyDescr
   integer(c_int) function cusparseSpSV_destroyDescr(spsvDescr &
    &            ) bind(C,name="hipsparseSpSV_destroyDescr")
     ! this destroys a handle for the triangle solver hipsparseSpSV
     ! and releases the memory associated
     use iso_c_binding
     implicit none
     type(c_ptr), value :: spsvDescr ! spsv context handle
   end function cusparseSpSV_destroyDescr

   ! hipsparseSpSV_buffersize
   integer(c_int) function cusparseSpSV_bufferSize(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            bufferSize) bind(C,name="hipsparseSpSV_bufferSize")
     ! this computes the buffersize needed by the triangle solver 
     ! hipsparseSpSV
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
     type(c_ptr),value   :: spsvDescr ! spsv context handle
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function cusparseSpSV_bufferSize

   ! hipsparseSpSV_buffersize_cmplx
   integer(c_int) function cusparseSpSV_bufferSize_cmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            bufferSize) bind(C,name="hipsparseSpSV_bufferSize")
     ! this computes the buffersize needed by the triangle solver hipsparseSpSV
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
     type(c_ptr),value   :: spsvDescr ! spsv context handle
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function cusparseSpSV_bufferSize_cmplx

   ! hipsparseSpSV_buffersize_dcmplx
   integer(c_int) function cusparseSpSV_bufferSize_dcmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            bufferSize) bind(C,name="hipsparseSpSV_bufferSize")
     ! this computes the buffersize needed by the triangle solver hipsparseSpSV
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
     type(c_ptr),value   :: spsvDescr ! spsv context handle
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function cusparseSpSV_bufferSize_dcmplx

   ! hipsparseSpSV_analysis
   integer(c_int) function cusparseSpSV_analysis(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            pBuffer ) bind(C,name="hipsparseSpSV_analysis")
     ! this does the analysis phase needed by the triangle solver hipsparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by hipsparseSpMV_bufferSize
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

   ! hipsparseSpSV_analysis_cmplx
   integer(c_int) function cusparseSpSV_analysis_cmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            pBuffer ) bind(C,name="hipsparseSpSV_analysis")
     ! this does the analysis phase needed by the triangle solver hipsparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by hipsparseSpMV_bufferSize
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

   ! hipsparseSpSV_analysis_dcmplx
   integer(c_int) function cusparseSpSV_analysis_dcmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            pBuffer ) bind(C,name="hipsparseSpSV_analysis")
     ! this does the analysis phase needed by the triangle solver hipsparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by hipsparseSpSV_bufferSize
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

   ! hipsparseSpSV_solve
   integer(c_int) function cusparseSpSV_solve(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, &
    &            spsvDescr) bind(C,name="hipsparseSpSV_solve")
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

   ! hipsparseSpSV_solve_cmplx
   integer(c_int) function cusparseSpSV_solve_cmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, &
    &            spsvDescr) bind(C,name="hipsparseSpSV_solve")
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

   ! hipsparseSpSV_solve_dcmplx
   integer(c_int) function cusparseSpSV_solve_dcmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, &
    &            spsvDescr) bind(C,name="hipsparseSpSV_solve")
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
   
   
   ! hipsparseCreateCsrilu0Info
   integer(c_int) function cusparseCreateCsrilu0Info(info) &
    &           bind(C,name="hipsparseCreateCsrilu02Info")
     ! this create the info pointer for the ilu0 context
     use iso_c_binding
     implicit none
     type(c_ptr)::info
   end function cusparseCreateCsrilu0Info

   ! hipsparseDestroyCsrilu0Info
   integer(c_int) function cusparseDestroyCsrilu0Info(info) &
    &           bind(C,name="hipsparseDestroyCsrilu02Info")
     ! this destroys the info pointer for the ilu0 context
     ! and frees its memory
     use iso_c_binding
     implicit none
     type(c_ptr), value ::info
   end function cusparseDestroyCsrilu0Info

   ! hipsparseDcsrilu0_bufferSize
   integer(c_int) function cusparseDcsrilu0_bufferSize(handle,m,nnz,descrA,&
    &           valA,rowA,colA,info,bufferSize) &
    &           bind(C,name="hipsparseDcsrilu02_bufferSize")
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

   ! hipsparseZcsrilu0_bufferSize
   integer(c_int) function cusparseZcsrilu0_bufferSize(handle,m,nnz,descrA,&
    &           valA,rowA,colA,info,bufferSize) &
    &           bind(C,name="hipsparseZcsrilu02_bufferSize")
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

   ! hipsparseDcsrilu0_analysis
   integer(c_int) function cusparseDcsrilu0_analysis(handle,n,nnz,descrA,&
    &           valA,rowA,colA,info,policy,pBuffer) &
    &           bind(C,name="hipsparseDcsrilu02_analysis")
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

   ! hipsparseZcsrilu0_analysis
   integer(c_int) function cusparseZcsrilu0_analysis(handle,n,nnz,descrA,&
    &           valA,rowA,colA,info,policy,pBuffer) &
    &           bind(C,name="hipsparseZcsrilu02_analysis")
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

   ! hipsparseDcsrilu0
   integer(c_int) function cusparseDcsrilu0(handle,m,nnz, &
    &           descrA,valA,rowA,colA,info, policy, pBuffer) &
    &           bind(C,name="hipsparseDcsrilu02")
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

   ! hipsparseZcsrilu0
   integer(c_int) function cusparseZcsrilu0(handle,m,nnz, &
    &           descrA,valA,rowA,colA,info, policy, pBuffer) &
    &           bind(C,name="hipsparseZcsrilu02")
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

   ! hipsparseXcsrilu0_zeroPivot
   integer(c_int) function cusparseXcsrilu0_zeroPivot(handle,info,loc) &
    &           bind(C,name="hipsparseXcsrilu02_zeroPivot")
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
   subroutine kernelc_s2d( single, double, count, stream) & 
    &              bind (C, name="kernelc_s2d" )
     use iso_c_binding
     implicit none
     ! force convert a chunk of memory from *single* to *double*
     type (c_ptr), value :: single, double ! pointers 
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
     type (c_ptr), value :: stream
   end subroutine kernelc_s2d

   ! kernelc_d2s 
   subroutine kernelc_d2s( double, single, count, stream) & 
    &              bind (C, name="kernelc_d2s" )
     use iso_c_binding
     implicit none
     ! force convert a chunk of memory from *double* to *single*
     type (c_ptr),value   :: double, single ! pointers
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
     type (c_ptr), value :: stream
   end subroutine kernelc_d2s

   ! kernelc_hadar
   subroutine kernelc_hadar( a, b, c, count, cstream) & 
    &              bind (C, name="kernelc_hadar" )
     use iso_c_binding
     implicit none
     ! perform hadamard multiply c = a*b, (real)
     type (c_ptr),value   :: a, b, c ! pointers
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
     type (c_ptr),value     :: cstream ! hipStream
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
     type (c_ptr),value     :: cstream ! hipStream
   end subroutine kernelc_hadac

   ! kernelc_xpbyc
   subroutine kernelc_xpbyc( a, b, c, count, cstream) & 
    &              bind (C, name="kernelc_xpbyc" )
     use iso_c_binding
     implicit none
     ! perform y = x + by, (complex)
     type (c_ptr),value          :: a, c ! pointers
     complex(c_double), value    :: b ! scaler b
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value      :: count 
     type (c_ptr),value          :: cstream ! hipStream
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
     type (c_ptr),value          :: cstream ! hipStream
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

  ! cf_hookCtx
   integer(c_int) function cf_hookDev(device_idx) & 
    &              bind (C, name="cf_hookDev" )
     ! bind the current context to a given device
     ! need some extra tweaks to use the "internal" functions instead of 
     ! the NVML lib - this is actually fine for CUDA, but not possible 
     ! for other APIs (AMD ROCM or Intel oneAPI) as they don't have the 
     ! equivalents of NVML
     use iso_c_binding
     implicit none
     ! the current device idx
     integer(c_int),value ::  device_idx
   end function cf_hookDev

  ! cf_resetFlag
   integer(c_int) function cf_resetFlag(device_idx) & 
    &              bind (C, name="cf_resetFlag" )
     ! reset the computation flag of the GPU device
     ! can be "hipDeviceScheduleSpin" or 
     ! "hipDeviceScheduleYield"
     use iso_c_binding
     implicit none
     ! the current device idx
     integer(c_int),value ::  device_idx
   end function cf_resetFlag

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

end module hipFortMap
