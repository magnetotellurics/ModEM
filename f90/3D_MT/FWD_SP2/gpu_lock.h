// gpu_lock.h — Platform-independent, cross-process GPU locking via POSIX 
// shared memory.
// Multiple MPI ranks (should...) safely target the same GPU:

#ifndef MODEM_GPU_LOCK_H
#define MODEM_GPU_LOCK_H

#include <atomic>
#include <fcntl.h>
#include <new>
#include <sched.h>
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

    // Exactly-once (re-)initialization of std::atomic members.
    // cnstr: 0=uninit/stale, 2=init-in-progress, 1=ready.
    //
    // Both 0 (fresh shm) and 1 (stale from a previous crash — some
    // occupied[] may still be DEVICE_IN_USE) trigger CAS to 2. The
    // winner does placement-new + reset, clearing all stale lock bits.
    // Losers yield while cnstr == 2. If the initializer crashes, cnstr
    // may revert to 0; the next iteration retries.
    //
    // A cleanly-terminated previous job always calls shm_unlink, so
    // finding an existing segment implies the last job crashed — it is
    // safe to reconstruct atomics unconditionally.

    while (1) {
        int old = __atomic_load_n(&g_lock->cnstr, __ATOMIC_ACQUIRE);

        if (old == 2) {
	        // Another rank is initializing; yield until it finishes.
            do {
		        // yield, instead of busy-waiting...
                sched_yield();
                old = __atomic_load_n(&g_lock->cnstr, __ATOMIC_ACQUIRE);
            } while (old == 2);
            continue;
        }

        // old is 0 (fresh) or 1 (stale ready with possibly stuck locks).
        // Try CAS {0, 1} → 2 to become the (re-)initializer.
        if (__atomic_compare_exchange_n(&g_lock->cnstr, &old, 2,
                                         false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE))
        {
            for (int i = 0; i < LOCK_MAX_DEVICES; i++)
                new (&g_lock->occupied[i]) std::atomic<int>(DEVICE_FREE);
            for (int i = 0; i < LOCK_MAX_DEVICES; i++)
                g_lock->occupied[i].store(DEVICE_FREE, std::memory_order_release);
            __atomic_store_n(&g_lock->cnstr, 1, __ATOMIC_RELEASE);
            break;
        }
        // else
        // CAS failed — another rank changed cnstr; loop and re-evaluates.
    }

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
