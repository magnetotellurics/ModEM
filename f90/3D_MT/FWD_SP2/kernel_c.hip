#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <complex.h>
#include "hip/hip_runtime.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <atomic>
#if defined(FG)
#include "hip/rccl.h"
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
    // HIP kernel implementing xpby:
    // y = x + b*y
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	        y[pos] = x[pos]+b*y[pos];
    }
}

__global__ void xpby_cmplx(const double *x, double br, double bi, double *y, const int N)
{
    // HIP kernel implementing xpby:
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
    // HIP kernel implementing p update:
    // p = r + beta * (p - omega*v)
    int pos = (blockDim.x*blockDim.y)*blockIdx.x+blockDim.x*threadIdx.y+threadIdx.x;
    if (pos<N) 
    {  
	    p[pos] = r[pos] + beta*(p[pos] - omega*v[pos]);
    }
}

__global__ void p_update_cmplx(const double *r, const double *v, const double br, const double bi, const double wr, const double wi, double *p, const int N)
{
    // HIP kernel implementing p update:
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
    // HIP kernel implementing p update:
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
    // HIP kernel implementing reduce s = sum(a)
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
		hipStream_t cstream)
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
		hipStream_t cstream)
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
		int Np, hipStream_t cstream)
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
		int Np, hipStream_t cstream)
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
extern "C" void kernelc_xpbyc(const double *x_d, double _Complex b, double *y_d, int Np, hipStream_t cstream)
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
extern "C" void kernelc_update_pc(const double *r_d, const double *v_d, double _Complex beta, double _Complex omega, double *p_d, int Np, hipStream_t cstream)
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
		int Np, hipStream_t cstream)
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
// This function is copied from rccl
static __inline__ int rcclTypeSize(rcclDataType_t type) {
  switch (type) {
    case rcclInt8:
    case rcclUint8:
      return 1;
    case rcclFloat16:
      return 2;
    case rcclInt32:
    case rcclUint32:
    case rcclFloat32:
      return 4;
    case rcclInt64:
    case rcclUint64:
    case rcclFloat64:
      return 8;
    default:
      return -1;
  }
}
// a homebrew allgatherv method by stitching the basic rccl apis togather  
extern "C" rcclResult_t rcclAllGatherV0(void *sendbuff, size_t sendcount,
	rcclDataType_t senddatatype, void *recvbuff,
	size_t recvcounts[], size_t recvdispls[], 
	int root, int n, rcclComm_t comm, hipStream_t stream)
{
    int rank_rccl;
    int size_rccl;
    hipError_t hipErr;
    rcclResult_t Err;
    // need to check the rccl interface here
    rcclCommUserRank(comm, &rank_rccl);
    rcclCommCount(comm, &size_rccl);
    // for debug
    // printf("#C rank = %i, size = %i \n", rank_rccl, size_rccl);
    // for debug
    // printf("#C rank = %i, sizes = %ld \n", rank_rccl, recvcounts[rank_rccl]);
    // for debug
    // printf("#C rank = %i, displs = %ld \n",rank_rccl, recvdispls[rank_rccl]);
    // firstly gather to root
    // send to root (if you are not)
    if(rank_rccl!=root) {
	// non-root send the data
        // printf("#C sending data to %i @ %i\n", 0, rank_rccl);
        Err = rcclSend(sendbuff, sendcount, senddatatype,
        	root, comm, stream);
        if(Err){
            return Err;
	}
    }
    else {
        // root firstly copy sendbuff to recvbuff 
        size_t self_displ = recvdispls[rank_rccl];
        hipErr = hipMemcpyAsync(static_cast<std::byte*>(recvbuff) + 
	    rcclTypeSize(senddatatype) * self_displ,
            sendbuff, sendcount * rcclTypeSize(senddatatype),
            hipMemcpyDeviceToDevice, stream);
	if (hipErr != hipSuccess) return rcclSystemError;
	// then recieve the data from other processes
        for(int i=0;i<size_rccl;++i){
            if(i == root) continue; //skip the root
            // printf("#C waiting for data from %i @ %i\n", i, rank_rccl);
            Err=rcclRecv(static_cast<std::byte*>(recvbuff)+
	        rcclTypeSize(senddatatype)*recvdispls[i],
                recvcounts[i], senddatatype, i, comm, stream);
            if(Err){
        	 return Err;
	    }
	}
    }
    // now broadcast to all
    size_t total = 0;
    for (int i = 0; i < size_rccl; i++){
        total += recvcounts[i];
    }
    // for debug
    // now broadcast the full vector
    // printf("#C now broadcasting @ %i\n", rank_rccl);
    Err=rcclBcast(recvbuff, total, senddatatype, 0, comm, stream);
    if(Err){
        return Err;
    }
    return rcclSuccess;
}
#
// a homebrew allgatherv method by stitching the basic rccl apis togather  
extern "C" rcclResult_t rcclAllGatherV(void *sendbuff, size_t sendcount,
	rcclDataType_t senddatatype, void *recvbuff,
        size_t recvcounts[], size_t recvdispls[],
        int root, int n, rcclComm_t comm, hipStream_t stream)
{
    int rank_rccl;
    int size_rccl;
    hipError_t hipErr;
    rcclResult_t Err;

    // Get the rank and size of the communicator
    Err = rcclCommUserRank(comm, &rank_rccl);
    if (Err != rcclSuccess) return Err;
    Err = rcclCommCount(comm, &size_rccl);
    if (Err != rcclSuccess) return Err;

    // for debug
    // printf("#C rank = %i, size = %i \n", rank_rccl, size_rccl);
    // for debug
    // printf("#C rank = %i, sizes = %ld \n", rank_rccl, recvcounts[rank_rccl]);
    // for debug
    // printf("#C rank = %i, displs = %ld \n",rank_rccl, recvdispls[rank_rccl]);

    // firstly everyone copies the local sendbuff to recvbuff
    size_t size_byte  = sendcount * rcclTypeSize(senddatatype);
    size_t displ_byte = recvdispls[rank_rccl] * rcclTypeSize(senddatatype);
    // printf("Base: %p, Offset: %zu bytes, Final: %p\n", recvbuff, displ_byte, 
    // 	    static_cast<void*>(static_cast<std::byte*>(recvbuff) + displ_byte));
    hipErr = hipMemcpyAsync(static_cast<std::byte*>(recvbuff) + displ_byte,
            sendbuff, size_byte, hipMemcpyDeviceToDevice,stream);
    if (hipErr != hipSuccess) { 
        printf("Failed to copy the device Mem: %s\n",  
           hipGetErrorString(hipErr));
        return rcclSystemError;
    }

    // Perform the allgather operation using a ring communication pattern
    // find the right and left neighbour
    int left = (rank_rccl - 1 + size_rccl) % size_rccl;
    int right = (rank_rccl + 1) % size_rccl;

    // we perform an in-place send and recv to avoid frequent alloc/dealloc
    // operations
    // For each rank, we'll send and receive data
    // except the local one (as we already have that) 
    for (int i = 1; i < size_rccl; ++i) {
        // Synchronize the gather operation
        Err = rcclGroupStart();
        if (Err != rcclSuccess) return Err;
	int crnt = (rank_rccl  - i + 1 + size_rccl) % size_rccl;
        // always Send data to the right neighbor, use the size 
	// and displs from the last iteration
        Err = rcclSend(static_cast<std::byte*>(recvbuff) + displ_byte,
                      recvcounts[crnt], senddatatype, right, comm, stream);
        if (Err) {
            return Err;
        }
        // printf("Send to right:  %i, from Rank: %i\n", right, rank_rccl);
        // printf("Send: %p, Offset: %zu bytes, Final: %p\n",
        //	recvbuff, displ_byte,  
        //	static_cast<void*>(static_cast<std::byte*>(recvbuff) + displ_byte));
	// size and displs to receive from the rank to the left
	int next = (rank_rccl  - i + size_rccl) % size_rccl;
	displ_byte = recvdispls[next] * rcclTypeSize(senddatatype);
        // always Receive data from the left neighbor
        Err = rcclRecv(static_cast<std::byte*>(recvbuff) + displ_byte,
                      recvcounts[next], senddatatype, left, comm, stream);
        if (Err) {
            return Err;
        }
        // printf("Recv from left:  %i, to Rank: %i\n", left, rank_rccl);
        // printf("Recv: %p, Offset: %zu bytes, Final: %p\n",
	// 	recvbuff, displ_byte, 
	//	static_cast<void*>(static_cast<std::byte*>(recvbuff) + displ_byte));
        // Synchronize the gather operation
        Err = rcclGroupEnd();
        if (Err != rcclSuccess) return Err;
    }

    return rcclSuccess;
}
#endif

// function called from main fortran program 
extern "C" int cf_resetFlag(int dev_idx)
{
    hipError_t err;
    unsigned int flag = 0;
    err = hipGetDeviceFlags(&flag);
    if (err != hipSuccess)
    {
        printf("Failed to get the device flags %i: %s\n", dev_idx, 
           hipGetErrorString(err));
    }
    // printf("Dev = %i, flag = %i \n", dev_idx, flag);
    if (flag == hipDeviceScheduleYield)
    {
        err = hipSetDeviceFlags(hipDeviceScheduleSpin);
        if (err != hipSuccess)
        {
            printf("Failed to set the device flag %i: %s\n", dev_idx, 
            hipGetErrorString(err));
            return 1;
        }
    }
    else
    {
        isleep(0.05);
        err = hipSetDeviceFlags(hipDeviceScheduleYield);
        if (err != hipSuccess)
        {
            printf("Failed to set the device flag %i: %s\n", dev_idx, 
            hipGetErrorString(err));
            return 1;
        }
    }
    return 0;
}

// ============================================================================
// Shared-memory GPU lock: cross-process atomic occupation flags.
// Identical logic to CUDA version
// cf_hookDev() acquires (CAS 0->1), cf_releaseDev() releases (set to 0).
// ============================================================================

static constexpr int LOCK_MAX_DEVICES = 64;

struct alignas(64) GpuLock {
    std::atomic<int> occupied[LOCK_MAX_DEVICES];
};

static GpuLock* g_lock       = nullptr;
static bool     g_lock_inited = false;

static int init_gpu_lock()
{
    const char* name = "/ModEM_gpu_lock";

    int fd = shm_open(name, O_RDWR, 0666);
    if (fd < 0) {
        fd = shm_open(name, O_CREAT | O_RDWR, 0666);
        if (fd < 0) return 1;
        if (ftruncate(fd, sizeof(GpuLock)) < 0) { close(fd); return 1; }
    }

    g_lock = static_cast<GpuLock*>(mmap(nullptr, sizeof(GpuLock),
                                         PROT_READ | PROT_WRITE,
                                         MAP_SHARED, fd, 0));
    close(fd);
    if (g_lock == MAP_FAILED) { g_lock = nullptr; return 1; }

    static std::atomic<int> init_guard(0);
    int expected = 0;
    if (init_guard.compare_exchange_strong(expected, 1)) {
        for (int i = 0; i < LOCK_MAX_DEVICES; i++)
            g_lock->occupied[i].store(0, std::memory_order_relaxed);
    }

    g_lock_inited = true;
    return 0;
}

// function called from main fortran program, this is different from the original cuda version
// as we do not have nvml in hip, we will just use the hip api to check the free memory after claiming the lock
extern "C" int cf_hookDev(int dev_idx)
{
    hipDeviceProp_t dev_prop;
    hipError_t err;
    size_t freeBytes, totalBytes, usedBytes;
    int nDevices;
    double memusage;
     
    // unsigned int flag = 0;
    // unsigned int i;
    err = hipGetDeviceCount(&nDevices);
    if (err != hipSuccess) 
    {
        printf("Error: Failed to inquire the number of devices %s\n", hipGetErrorString(err));
        return 1;
    }
    // see if we have enough number of device
    if (dev_idx + 1 > nDevices) 
    { 
        printf(" inquired device idx = %i, nDevices = %i \n", dev_idx, nDevices);
        printf(" Error: device idx out of range! \n");
        goto Error;
    }
    // get the basic properties
    err = hipGetDeviceProperties(&dev_prop, dev_idx);
    if (err != hipSuccess )
    {
        printf("Failed to get device property %i: %i \n", dev_idx, err);
        goto Error;
    }
    //printf("Device Idx: %d, name: %s \n", dev_idx, dev_prop.name);
    //printf("  Memory Clock Rate (KHz): %d\n",
    //         dev_prop.memoryClockRate);
    //printf("  Memory Bus Width (bits): %d\n",
    //         dev_prop.memoryBusWidth);
    //printf("  Peak Memory Bandwidth (GB/s): %f\n",
    //         2.0*dev_prop.memoryClockRate*(dev_prop.memoryBusWidth/8)/1.0e6);
    //printf("  PCI bus ID: %d\n",
    //         dev_prop.pciBusID);
    //printf("  is integrated GPU : %i\n",
    //         dev_prop.integrated);
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
    // set the device context before claiming the lock.
    err = hipSetDevice(dev_idx);
    if (err != hipSuccess)
    {
        printf("Failed to set device %i: %u \n", dev_idx, err);
        return 1;
    }

    // now try to get affiliated device
    // Atomic lock loop: CAS 0->1 to claim the GPU
    while (true)
    {
        int expected = 0;
        if (g_lock->occupied[dev_idx].compare_exchange_strong(
                expected, 1, std::memory_order_acq_rel))
        {
            // CLAIMED -- we are the sole occupant
            // Safety: verify memory is actually available after acquiring lock
            err = hipMemGetInfo(&freeBytes, &totalBytes);
            if (err != hipSuccess)
            {
                g_lock->occupied[dev_idx].store(0, std::memory_order_release);
                printf("Failed to get usage for device %i: %s\n", dev_idx,
                       hipGetErrorString(err));
                goto Error;
            }
            err = hipSetDeviceFlags(hipDeviceScheduleYield);
	    if (err != hipSuccess)
	    {
		g_lock->occupied[dev_idx].store(0, std::memory_order_release);
		printf("Failed to set the device flag %i: %s\n", dev_idx, 
		       hipGetErrorString(err));
		goto Error;
	    }
            printf(" # Dev Selected: %i. %s \n", dev_idx, dev_prop.name);
	    // unlike the cuda version, we do not have the nvml to report the GPU memory usage
	    // just calculate the used memory percentage and report that
            usedBytes = totalBytes - freeBytes;
            memusage = 100.0*usedBytes/totalBytes;
            printf(" # Dev Status  : GPU-mem = %5.2f %%, PCI %i: %i\n", 
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

