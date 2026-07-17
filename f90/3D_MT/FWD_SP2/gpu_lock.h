// gpu_lock.h — Platform-independent, cross-process GPU locking via POSIX 
// shared memory.
// Multiple MPI ranks (should...) safely target the same GPU:

#ifndef MODEM_GPU_LOCK_H
#define MODEM_GPU_LOCK_H

#include <atomic>
#include <fcntl.h>
#include <new>
#include <sys/mman.h>
#include <unistd.h>

// Define device state flags                                                 

#define DEVICE_FREE  0
#define DEVICE_IN_USE 1

// Lock structure & global state 
// we probably don't have more than 64 devices on a single node
static constexpr int LOCK_MAX_DEVICES = 64;

// cnstr: compiler-built-in atomic on raw int — safe before object lifetime begins.
// All other members are std::atomic<int>, constructed by placement-new once,
// then used via normal C++ atomic operations.
struct alignas(64) GpuLock {
    int cnstr;      // 0 = not constructed, 1 = atomics live, 2 = init in progress
    std::atomic<int> occupied[LOCK_MAX_DEVICES];
};

static GpuLock* g_lock       = nullptr;
static bool     g_lock_inited = false;

// Internal: create / map shared-memory segment and construct atomics 
static inline int init_gpu_lock()
{
    const char* name = "/ModEM_gpu_lock";

    // Try to open existing shared memory segment first 
    int fd = shm_open(name, O_RDWR, 0600);
    if (fd < 0) {
        // Doesn't exist -- create it 
        fd = shm_open(name, O_CREAT | O_RDWR, 0600);
        if (fd < 0) return 1;
        if (ftruncate(fd, sizeof(GpuLock)) < 0) { close(fd); return 1; }
    }

    g_lock = static_cast<GpuLock*>(mmap(nullptr, sizeof(GpuLock),
                                         PROT_READ | PROT_WRITE,
                                         MAP_SHARED, fd, 0));
    close(fd);
    if (g_lock == MAP_FAILED) { g_lock = nullptr; return 1; }

    // Coordinate exactly-once construction of std::atomic members.   
    // cnstr uses __atomic_* built-ins: safe on raw storage       
    // cnstr: 0 = not constructed, 1 = atomics live, 2 = init in progress
    // the winner does placement-new on every std::atomic<int>     
    // the other (losers) spin-wait on 2, until it's 1, 
    // then proceed to use the atomics.
    // so CAS either succeeds, reveals cnstr=2 (spin), or reveals cnstr=1
    // (ready).
    // useful for stale lock state if a previous run crashed or didn't 
    // clean up properly. I only encountered this once in HIP, but not in 
    // CUDA, so it may be a HIP bug.

    int old = __atomic_load_n(&g_lock->cnstr, __ATOMIC_ACQUIRE);

    if (old == 2) {
        while (__atomic_load_n(&g_lock->cnstr, __ATOMIC_ACQUIRE) == 2)
            ;
    } else if (__atomic_compare_exchange_n(&g_lock->cnstr, &old, 2,
                          false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE))
    {
        for (int i = 0; i < LOCK_MAX_DEVICES; i++)
            new (&g_lock->occupied[i]) std::atomic<int>(DEVICE_FREE);
        for (int i = 0; i < LOCK_MAX_DEVICES; i++)
            g_lock->occupied[i].store(DEVICE_FREE, std::memory_order_release);
        __atomic_store_n(&g_lock->cnstr, 1, __ATOMIC_RELEASE);
    } else if (old == 2) {
        while (__atomic_load_n(&g_lock->cnstr, __ATOMIC_ACQUIRE) == 2)
            ;
    }
    // else old == 1: someone completed init before our CAS, no action needed

    g_lock_inited = true;
    return 0;
}

// Public C-bindings (called from Fortran)                            

extern "C" void cf_releaseDev(int dev_idx)
{
    if (g_lock_inited && g_lock != nullptr &&
        dev_idx >= 0 && dev_idx < LOCK_MAX_DEVICES)
    {
        g_lock->occupied[dev_idx].store(DEVICE_FREE, std::memory_order_release);
    }
}

extern "C" void cf_cleanupLock()
{
    if (g_lock != nullptr && g_lock != MAP_FAILED) {
        // also reset cnstr so the next init_gpu_lock starts cleanly
        __atomic_store_n(&g_lock->cnstr, 0, __ATOMIC_RELEASE);
        munmap(g_lock, sizeof(GpuLock));
    }
    shm_unlink("/ModEM_gpu_lock");
    g_lock = nullptr;
    g_lock_inited = false;
}

#endif // MODEM_GPU_LOCK_H
