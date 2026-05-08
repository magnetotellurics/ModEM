#include <new>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <atomic>
#if defined(FG) 
#include <nccl.h>
#endif

// simple kernel function that converts double vectors to single
__global__ void real64to32(const double *in, float *out, const int N)
{
    // a position every 64 bits
    // int pos = blockDim.x * blockIdx.x + threadIdx.x;
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	float sing = __double2float_rn(in[pos]);
	out[pos] = sing;
    }
}

// simple kernel function that converts single vectors to double
__global__ void real32to64(const uint32_t *in, uint64_t *out, const int N)
{
    // a position every 32 bits
    // int pos = blockDim.x * blockIdx.x + threadIdx.x;
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	// last sign bit (1)
	uint64_t s = in[pos] & 0x80000000;
	// exponent bits (8)
	uint64_t e = ((in[pos] & 0x7f800000) >> 23);
	e = e + 896;
	// mantissa bits (23)
	uint64_t m = in[pos] & 0x7fffff;
	// double through bitwise or
	uint64_t doub = (s<<32) | (e<<52) | (m<<29);
        // a new position every 64 bits
	out[pos] = doub;
    }
 }

__global__ void hada_real(const double *ina, const double *inb, double *out, const int N)
{
    //hadamard multiplication, as real 
    // c = a * b
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	out[pos] = ina[pos]*inb[pos];
    }
 }

__global__ void hada_cmplx(const double *ina, const double *inb, double *out, const int N)
{
    //hadamard multiplication, as complex
    //(a + bi) * (c + di) = (ac -bd) + (ad + bc)i
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
        if ( pos%2 == 0) // 0 2 4 
        {
	        out[pos] = ina[pos]*inb[pos]-ina[pos+1]*inb[pos+1];
        }
        else // 1 3 5 
        {
	        out[pos] = ina[pos-1]*inb[pos]+ina[pos]*inb[pos-1];
        }
    }
 }

__global__ void xpby_real(const double *x, const double b, double *y, const int N)
{
    // CUDA kernel implementing xpby:
    // y = x + b*y
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	        y[pos] = x[pos]+b*y[pos];
    }
}

__global__ void xpby_cmplx(const double *x, double br, double bi, double *y, const int N)
{
    // CUDA kernel implementing xpby:
    // y = x + b*y
    // complex makes this a little tricky - need a temp variable
    int pos = ((blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x)*2;
    if (pos<N) 
    {  
        // lower 64bit --> real
	    double lower = x[pos] + br*y[pos] - bi*y[pos+1];
        // higher 64bit --> imag
	    y[pos+1] = x[pos+1] + br*y[pos+1] + bi*y[pos];
        // and copy them back
        y[pos] = lower;
    }
}

__global__ void p_update_real(const double *v, const double *r, const double beta, const double omega, double *p, const int N)
{
    // CUDA kernel implementing p update:
    // p = r + beta * (p - omega*v)
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	    p[pos] = r[pos] + beta*(p[pos] - omega*v[pos]);
    }
}

__global__ void p_update_cmplx(const double *r, const double *v, const double br, const double bi, const double wr, const double wi, double *p, const int N)
{
    // CUDA kernel implementing p update:
    // p = r + beta * (p - omega*v)
    int pos = ((blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x)*2;
    if (pos<N) 
    {  
        // lower 64bit --> real p - omega*v
	    double lower = p[pos] - wr*v[pos] + wi*v[pos+1];
        // higher 64bit --> imag p - omega*v
	    double higher = p[pos+1] - wr*v[pos+1] - wi*v[pos];
        // lower 64bit --> real r + beta * (p - omega*v)
	    p[pos] = r[pos] + br*lower - bi*higher;
        // higher 64bit --> imag r + beta * (p - omega*v)
	    p[pos+1] = r[pos+1] + br*higher + bi*lower;
    }
}

__global__ void x_update_cmplx(const double *p, const double *s, const double ar, const double ai, const double wr, const double wi, double *x, const int N)
{
    // CUDA kernel implementing p update:
    int pos = ((blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x)*2;
    if (pos<N) 
    {  
        // x = x + alpha * ph + omega * sh
        // lower 64bit --> real x + alpha*ph
	    x[pos] = x[pos] + ar*p[pos] - ai*p[pos+1] + wr*s[pos] - wi*s[pos];
        // higher 64bit --> imag x + alpha*ph
	    x[pos+1]= x[pos+1] + ar*p[pos+1] + ai*p[pos] + wr*s[pos+1] + wi*s[pos];
    }
}

__global__ void reduce_real(double *a, double s, const int N)
{
    // CUDA kernel implementing reduce s = sum(a)
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    int number_of_threads = (blockDim.x*blockDim.y);
    int step = 1; //initial step
    while (number_of_threads > 0)
    {
        if (pos < number_of_threads) 
        {
            const int left = pos * step*2; 
            const int right = left + step;
            a[left] += a[right];
        }
        step <<= 1;
        number_of_threads >>= 1;
        __syncthreads();
    }
}

// function to wait for sometime (accepts fractional seconds)
void isleep(float seconds)
{
    // Converting time into milli_seconds
    int milli_seconds = int(1000.0 * seconds);
    // Storing start time
    clock_t start_time = clock();
    // looping till required time is not achieved
    while (clock() < start_time + milli_seconds)
        ;
}

// function called from main fortran program
// single to double conversion
extern "C" void kernelc_s2d(const uint32_t *a_d, uint64_t *b_d, int Np,
		cudaStream_t cstream)
{
    //uint32_t  *a_d;  // declare GPU vector double (but stored as uint)
    //uint64_t  *b_d;  // declare GPU vector float
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = Np*2;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // call function on GPU
    real32to64<<< grids, blocks, 0, cstream >>>( a_d, b_d, N);

    return;
}

// function called from main fortran program
// double to single conversion
extern "C" void kernelc_d2s(const double *a_d, float *b_d, int Np,
		cudaStream_t cstream)
{
    //double  *a_d;  // declare GPU vector double 
    //float  *b_d;  // declare GPU vector float
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = Np*2;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // call function on GPU
    real64to32<<< grids, blocks, 0, cstream >>>( a_d, b_d, N);

    return;
}

// function called from main fortran program
// hadamard multiply real version
extern "C" void kernelc_hadar(const double *a_d, const double *b_d, double *c_d,
		int Np, cudaStream_t cstream)
{
    //double  *a_d;  // declare GPU vector double 
    //double  *b_d;  // declare GPU vector double
    //double  *c_d;  // declare GPU vector double
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = Np;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // call function on GPU
    hada_real<<< grids, blocks, 0, cstream >>>( a_d, b_d, c_d, N);

    return;
}

// function called from main fortran program
// hadamard multiply complex version
extern "C" void kernelc_hadac(const double *a_d, const double *b_d, double *c_d,
		int Np, cudaStream_t cstream)
{
    //double  *a_d;  // declare GPU vector double 
    //double  *b_d;  // declare GPU vector double
    //double  *c_d;  // declare GPU vector double
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = Np*2;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // call function on GPU
    hada_cmplx<<< grids, blocks, 0, cstream >>>( a_d, b_d, c_d, N);

    return;
}

// function called from main fortran program
// xpby complex version
extern "C" void kernelc_xpbyc(const double *x_d, double _Complex b, double *y_d,
		int Np, cudaStream_t cstream)
{
    //double  *x_d;  // declare GPU vector double 
    //double  *y_d;  // declare GPU vector double
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = 2*Np;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // get upper 64 bit
    double br = creal(b);
    // get lower 64 bit
    double bi = cimag(b);
    // call function on GPU
    xpby_cmplx<<< grids, blocks, 0, cstream>>>( x_d, br, bi, y_d, N);

    return;
}

// function called from main fortran program
// update p complex version
extern "C" void kernelc_update_pc(const double *r_d, const double *v_d, double _Complex beta, double _Complex omega, double *p_d, int Np, cudaStream_t cstream)
{
    //double  *a_d;  // declare GPU vector double 
    //double  *b_d;  // declare GPU vector double
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = 2*Np;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // get upper 64 bit
    double br = creal(beta);
    // get lower 64 bit
    double bi = cimag(beta);
    // get upper 64 bit
    double wr = creal(omega);
    // get lower 64 bit
    double wi = cimag(omega);
    // call function on GPU
    p_update_cmplx<<< grids, blocks, 0, cstream >>>( r_d, v_d, br, bi, wr, wi, p_d, N);

    return;
}

// function called from main fortran program
// update x complex version
extern "C" void kernelc_update_xc(const double *ph_d, const double *sh_d, 
		double _Complex alpha, double _Complex omega, double *x_d, 
		int Np, cudaStream_t cstream)
{
    //double  *ph_d;  // declare GPU vector double 
    //double  *sh_d;  // declare GPU vector double
    //double  *x_d;  // declare GPU vector double
           
    int ngrid;              // number of grid
    int N;                  // length of N 

    N = 2*Np;
    ngrid = N/256 + 1;
    dim3 grids(ngrid,1,1);
    dim3 blocks(32,8,1);
    // get upper 64 bit
    double ar = creal(alpha);
    // get lower 64 bit
    double ai = cimag(alpha);
    // get upper 64 bit
    double wr = creal(omega);
    // get lower 64 bit
    double wi = cimag(omega);
    // call function on GPU
    x_update_cmplx<<< grids, blocks, 0, cstream >>>( ph_d, sh_d, ar, ai, wr, wi, x_d, N);

    return;
}

#if defined(FG) 
// This function is copied from nccl
static __inline__ int ncclTypeSize(ncclDataType_t type) {
  switch (type) {
    case ncclInt8:
    case ncclUint8:
      return 1;
    case ncclFloat16:
      return 2;
    case ncclInt32:
    case ncclUint32:
    case ncclFloat32:
      return 4;
    case ncclInt64:
    case ncclUint64:
    case ncclFloat64:
      return 8;
    default:
      return -1;
  }
}
// a homebrew allgatherv method by stitching the basic nccl apis togather  
extern "C" ncclResult_t ncclAllGatherV(void *sendbuff, size_t sendcount,
	ncclDataType_t senddatatype, void *recvbuff,
	size_t recvcounts[], size_t recvdispls[], 
	int root, int n, ncclComm_t comm, cudaStream_t stream)
{
    int rank_nccl;
    int size_nccl;
    cudaError_t cuErr;
    ncclResult_t Err;
    // need to check the nccl interface here
    ncclCommUserRank(comm, &rank_nccl);
    ncclCommCount(comm, &size_nccl);
    // for debug
    // printf("#C rank = %i, size = %i \n", rank_nccl, size_nccl);
    // for debug
    // printf("#C rank = %i, sizes = %ld \n", rank_nccl, recvcounts[rank_nccl]);
    // for debug
    // printf("#C rank = %i, displs = %ld \n",rank_nccl, recvdispls[rank_nccl]);
    // firstly gather to root
    // send to root (if you are not)
    if(rank_nccl!=root) {
	// non-root send the data
        // printf("#C sending data to %i @ %i\n", 0, rank_nccl);
        Err = ncclSend(sendbuff, sendcount, senddatatype,
        	root, comm, stream);
        if(Err){
            return Err;
	}
    }
    else {
	// rootrecieve the data from other processes
        for(int i=0;i<size_nccl;++i){
            if(i == root) continue; //skip the root
                // printf("#C waiting for data from %i @ %i\n", i, rank_nccl);
            Err=ncclRecv(static_cast<std::byte*>(recvbuff)+
	        ncclTypeSize(senddatatype)*recvdispls[i],
                recvcounts[i], senddatatype, i, comm, stream);
            if(Err){
        	 return Err;
	    }
	}
        // root copy sendbuff to recvbuff 
        size_t self_displ = recvdispls[rank_nccl];
        cuErr = cudaMemcpyAsync(static_cast<std::byte*>(recvbuff) + 
	    ncclTypeSize(senddatatype) * self_displ,
            sendbuff, sendcount * ncclTypeSize(senddatatype),
            cudaMemcpyDeviceToDevice, stream);
	if (cuErr != cudaSuccess) return ncclSystemError;
    }
    // now broadcast to all
    size_t total = 0;
    for (int i = 0; i < size_nccl; i++){
        total += recvcounts[i];
    }
    // for debug
    // now broadcast the full vector
    // printf("#C now broadcasting @ %i\n", rank_nccl);
    Err=ncclBcast(recvbuff, total, senddatatype, 0, comm, stream);
    if(Err){
        return Err;
    }
    return ncclSuccess;
}

// a homebrew allgatherv method by stitching the basic nccl apis togather  
extern "C" ncclResult_t ncclAllGatherV2(void *sendbuff, size_t sendcount,
	ncclDataType_t senddatatype, void *recvbuff,
        size_t recvcounts[], size_t recvdispls[],
        int root, int n, ncclComm_t comm, cudaStream_t stream)
{
    int rank_nccl;
    int size_nccl;
    cudaError_t cuErr;
    ncclResult_t Err;

    // Get the rank and size of the communicator
    Err = ncclCommUserRank(comm, &rank_nccl);
    if (Err != ncclSuccess) return Err;
    Err = ncclCommCount(comm, &size_nccl);
    if (Err != ncclSuccess) return Err;

    // for debug
    // printf("#C rank = %i, sizes = %ld \n", rank_nccl, recvcounts[rank_nccl]);
    // for debug
    // printf("#C rank = %i, displs = %ld \n",rank_nccl, recvdispls[rank_nccl]);

    // firstly everyone copies the local sendbuff to recvbuff
    size_t size_byte  = sendcount * ncclTypeSize(senddatatype);
    size_t displ_byte = recvdispls[rank_nccl] * ncclTypeSize(senddatatype);
    cuErr = cudaMemcpyAsync(static_cast<std::byte*>(recvbuff) + displ_byte,
            sendbuff, size_byte, cudaMemcpyDeviceToDevice,stream);
    if (cuErr != cudaSuccess) { 
        printf("Failed to copy the device Mem: %s\n",  
           cudaGetErrorString(cuErr));
        return ncclSystemError;
    }

    // Perform the allgather operation using a ring communication pattern
    // find the right and left neighbour
    int left = (rank_nccl - 1 + size_nccl) % size_nccl;
    int right = (rank_nccl + 1) % size_nccl;

    // we perform an in-place send and recv to avoid frequent alloc/dealloc
    // operations
    // For each rank, we'll send and receive data
    // except the local one (as we already have that) 
    for (int i = 1; i < size_nccl; ++i) {
	int crnt = (rank_nccl  - i + 1 + size_nccl) % size_nccl;
        // Synchronize the gather operation
        Err = ncclGroupStart();
        if (Err != ncclSuccess) return Err;
        // always Send data to the right neighbor, use the size 
	// and displs from the last iteration
        Err = ncclSend(static_cast<std::byte*>(recvbuff) + displ_byte,
                      recvcounts[crnt], senddatatype, right, comm, stream);
        if (Err) {
            return Err;
        }
	int next = (rank_nccl  - i + size_nccl) % size_nccl;
	displ_byte = recvdispls[next] * ncclTypeSize(senddatatype);
        // always Receive data from the left neighbor
        Err = ncclRecv(static_cast<std::byte*>(recvbuff) + displ_byte,
                      recvcounts[next], senddatatype, left, comm, stream);
        if (Err) {
            return Err;
        }
        Err = ncclGroupEnd();
        if (Err != ncclSuccess) return Err;
    }

    return ncclSuccess;
}

// function called from main fortran program 
#endif

extern "C" int cf_resetFlag(int dev_idx)
{
    cudaError_t err;
    unsigned int flag = 0;
    err = cudaGetDeviceFlags(&flag);
    if (err != cudaSuccess)
    {
        printf("Failed to get the device flags %i: %s\n", dev_idx, 
           cudaGetErrorString(err));
    }
    // printf("Dev = %i, flag = %i \n", dev_idx, flag);
    if (flag == cudaDeviceScheduleYield)
    {
        err = cudaSetDeviceFlags(cudaDeviceScheduleSpin);
        if (err != cudaSuccess)
        {
            printf("Failed to set the device flag %i: %s\n", dev_idx, 
            cudaGetErrorString(err));
            return 1;
        }
    }
    else
    {
        isleep(0.05);
        err = cudaSetDeviceFlags(cudaDeviceScheduleYield);
        if (err != cudaSuccess)
        {
            printf("Failed to set the device flag %i: %s\n", dev_idx, 
            cudaGetErrorString(err));
            return 1;
        }
    }
    return 0;
}

// ============================================================================
// Shared-memory GPU lock: cross-process atomic occupation flags.
// Multiple MPI ranks can target the same GPU. This uses a POSIX shared-memory
// region containing std::atomic<int> flags, one per device index.
// cf_hookDev() acquires (CAS 0->1), cf_releaseDev() releases (store 0).
// ============================================================================

static constexpr int LOCK_MAX_DEVICES = 64;

// Raw byte flag at offset 0 — no C++ object lifetime requirements.
// Used *only* to coordinate who runs placement-new on the atomics.
struct alignas(64) GpuLock {
    unsigned char constructed;            // 0 = raw storage, !=0 = atomics live
    std::atomic<int> initialized;         // cross-process init guard (0->1 CAS)
    std::atomic<int> occupied[LOCK_MAX_DEVICES];
};

static GpuLock* g_lock       = nullptr;
static bool     g_lock_inited = false;

static int init_gpu_lock()
{
    const char* name = "/ModEM_gpu_lock";
    // Try to open existing shared memory segment first
    int fd = shm_open(name, O_RDWR, 0666);
    if (fd < 0) {
        // Doesn't exist -- create it
        fd = shm_open(name, O_CREAT | O_RDWR, 0666);
        if (fd < 0) return 1;
        if (ftruncate(fd, sizeof(GpuLock)) < 0) { close(fd); return 1; }
    }

    g_lock = static_cast<GpuLock*>(mmap(nullptr, sizeof(GpuLock),
                                         PROT_READ | PROT_WRITE,
                                         MAP_SHARED, fd, 0));
    close(fd);
    if (g_lock == MAP_FAILED) { g_lock = nullptr; return 1; }

    // Begin object lifetimes once across all processes...
    if (__atomic_test_and_set(&g_lock->constructed, __ATOMIC_ACQ_REL)) {
        // we are not the first - just wait...
    } else {
        // we are the first - need to construct the atomics in shared memory
        // begin lifetimes of all std::atomic<int> via placement-new.
        new (&g_lock->initialized) std::atomic<int>(0);
        for (int i = 0; i < LOCK_MAX_DEVICES; i++)
            new (&g_lock->occupied[i]) std::atomic<int>(0);
    }

    g_lock_inited = true;
    return 0;
}

// function called from main fortran program 
extern "C" int cf_hookDev(int dev_idx)
{
    cudaDeviceProp dev_prop;
    cudaError_t err;
    size_t freeBytes, totalBytes, usedBytes;
    int nDevices;
    double memusage;
     
    // unsigned int flag = 0;
    // unsigned int i;
    err = cudaGetDeviceCount(&nDevices);
    if (err != cudaSuccess) 
    {
        printf("Error: Failed to inquire the number of devices %s\n", cudaGetErrorString(err));
        return 1;
    }
    // Validate dev_idx before any lock or CUDA access:
    if (dev_idx < 0 || dev_idx >= nDevices || dev_idx >= LOCK_MAX_DEVICES)
    {
        printf("Error: device idx out of range! \n");
        printf(" inquired device idx = %i, nDevices = %d, LOCK_MAX_DEVICES = %d \n",
               dev_idx, nDevices, LOCK_MAX_DEVICES);
        goto Error;
    }
    err = cudaGetDeviceProperties(&dev_prop, dev_idx);
    //printf("Device Idx: %d, name: %s \n", dev_idx, dev_prop.name);
    //printf("  Memory Clock Rate (KHz): %d\n",
    //         dev_prop.memoryClockRate);
    //printf("  Memory Bus Width (bits): %d\n",
    //         dev_prop.memoryBusWidth);
    //printf("  Peak Memory Bandwidth (GB/s): %f\n",
    //         2.0*dev_prop.memoryClockRate*(dev_prop.memoryBusWidth/8)/1.0e6);
    //printf("  PCI bus ID: %i\n",
    //         dev_prop.pciBusID);
    //printf("  is integrated GPU : %i\n",
    //         dev_prop.integrated);
    // now try to get 
    if (dev_prop.integrated == 1 ) // is integrated gpu
    {
	printf(" WARNING: integrated device found for device: %i ! \n", dev_idx);
	printf(" this may lead to unsupported computation error...\n");
    }
    // Init shared-memory lock (note that any rank can trigger)
    if (!g_lock_inited && init_gpu_lock() != 0)
    {
        printf("Error: failed to init GPU lock shared memory\n");
        goto Error;
    }
    //now try to hook on to that device
    err = cudaSetDevice(dev_idx);
	if( err != cudaSuccess )
    {
        printf("Failed to set device %i: %u \n", dev_idx, err);
        goto Error;
    }

    // Atomic lock loop: CAS 0->1 to claim the GPU
    while (true)
    {
        int expected = 0;
        if (g_lock->occupied[dev_idx].compare_exchange_strong(
                expected, 1, std::memory_order_acq_rel))
	{
            // CLAIMED -- we are the sole occupant
            // Safety: verify memory is actually available after acquiring lock
            err = cudaMemGetInfo(&freeBytes, &totalBytes);
            if (err != cudaSuccess)
            {
                g_lock->occupied[dev_idx].store(0, std::memory_order_release);
                printf("Failed to get usage for device %i: %s\n", dev_idx,
                       cudaGetErrorString(err));
                goto Error;
            }
            err = cudaSetDeviceFlags(cudaDeviceScheduleYield);
	    if (err != cudaSuccess)
	    {
		g_lock->occupied[dev_idx].store(0, std::memory_order_release);
		printf("Failed to set the device flag %i: %s\n", dev_idx, 
		       cudaGetErrorString(err));
		goto Error;
	    }
            printf(" # Dev Selected: %i. %s \n", dev_idx, dev_prop.name);
	    // just calculate the memory usage and report it 
            usedBytes = totalBytes - freeBytes;
            memusage = 100.0*usedBytes/totalBytes;
            printf(" # Dev Status  : GPU-mem = %f %%, PCI %i: %i\n", 
                   memusage, dev_prop.pciBusID, dev_prop.pciDeviceID );
	    return 0;
        }
        // else let it spin
        isleep(0.05);
    }
Error:
    printf("Error: Failed to attach to device: %i \n", dev_idx);
    return 1;
}

// Release the GPU lock so other processes can hook on
extern "C" void cf_releaseDev(int dev_idx)
{
    if (g_lock_inited && g_lock != nullptr &&
        dev_idx >= 0 && dev_idx < LOCK_MAX_DEVICES)
    {
        g_lock->occupied[dev_idx].store(0, std::memory_order_release);
    }
}

// Cleanup shared memory on program exit
extern "C" void cf_cleanupLock()
{
    if (g_lock != nullptr && g_lock != MAP_FAILED)
        munmap(g_lock, sizeof(GpuLock));
    shm_unlink("/ModEM_gpu_lock");
    g_lock = nullptr;
    g_lock_inited = false;
}
