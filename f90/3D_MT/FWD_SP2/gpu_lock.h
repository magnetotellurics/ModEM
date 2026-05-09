/*
 * gpu_lock.h — Platform-independent, cross-process GPU locking via POSIX shared memory.
 *
 * Multiple MPI ranks can safely target the same GPU:
 *   cf_hookDev()    — acquires GPU via CAS 0->1 on an atomic flag array (platform-specific)
 *   cf_releaseDev() — releases the lock so other processes can attach
 *   cf_cleanupLock()— unmaps and unlinks shared memory on exit
 */

#ifndef MODEM_GPU_LOCK_H
#define MODEM_GPU_LOCK_H

#include <atomic>
#include <fcntl.h>
#include <new>
#include <sys/mman.h>
#include <unistd.h>

/* ------------------------------------------------------------------ */
/* Device state flags                                                  */
/* ------------------------------------------------------------------ */

#define DEVICE_FREE  0
#define DEVICE_IN_USE 1

/* ------------------------------------------------------------------ */
/* Lock structure & global state                                       */
/* ------------------------------------------------------------------ */

static constexpr int LOCK_MAX_DEVICES = 64;

/*
 * Raw byte flag at offset 0 — no C++ object lifetime requirements.
 * Used *only* to coordinate who runs placement-new on the atomics.
 */
struct alignas(64) GpuLock {
    unsigned char constructed;             /* 0 = raw storage, !=0 = atomics live */
    std::atomic<int> initialized;          /* cross-process init guard (0->1 CAS) */
    std::atomic<int> occupied[LOCK_MAX_DEVICES];
};

static GpuLock* g_lock       = nullptr;
static bool     g_lock_inited = false;

/* ------------------------------------------------------------------ */
/* Internal: create / map shared-memory segment and construct atomics  */
/* ------------------------------------------------------------------ */

static inline int init_gpu_lock()
{
    const char* name = "/ModEM_gpu_lock";

    /* Try to open existing shared memory segment first */
    int fd = shm_open(name, O_RDWR, 0600);
    if (fd < 0) {
        /* Doesn't exist -- create it */
        fd = shm_open(name, O_CREAT | O_RDWR, 0600);
        if (fd < 0) return 1;
        if (ftruncate(fd, sizeof(GpuLock)) < 0) { close(fd); return 1; }
    }

    g_lock = static_cast<GpuLock*>(mmap(nullptr, sizeof(GpuLock),
                                         PROT_READ | PROT_WRITE,
                                         MAP_SHARED, fd, 0));
    close(fd);
    if (g_lock == MAP_FAILED) { g_lock = nullptr; return 1; }

    /* Begin object lifetimes once across all processes */
    if (__atomic_test_and_set(&g_lock->constructed, __ATOMIC_ACQ_REL)) {
        /* we are not the first -- just wait */
    } else {
        /* we are the first -- construct the atomics in shared memory */
        new (&g_lock->initialized) std::atomic<int>(0);
        for (int i = 0; i < LOCK_MAX_DEVICES; i++)
            new (&g_lock->occupied[i]) std::atomic<int>(DEVICE_FREE);
    }

    g_lock_inited = true;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Public C-bindings (called from Fortran)                              */
/* ------------------------------------------------------------------ */

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
    if (g_lock != nullptr && g_lock != MAP_FAILED)
        munmap(g_lock, sizeof(GpuLock));
    shm_unlink("/ModEM_gpu_lock");
    g_lock = nullptr;
    g_lock_inited = false;
}

#endif /* MODEM_GPU_LOCK_H */
