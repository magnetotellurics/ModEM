#include <sys/resource.h>
#include <stdio.h>

/* It should be noted that the below function might have different
 * results depending on the system you are on. POSIX doesn't guarantee
 * ru_maxrss even be implemented.
 *
 * It's on most systems, but your results may very.
 */

/*
* To call this C function from Fortran, use this interface:
*
* interface
*   subroutine get_maxrss(maxrss_bytes) bind(c)
*       use iso_c_binding, only : c_long
*       integer (c_long), intent(out) :: maxrss_bytes
*   end subroutine get_maxrss
* end interface
*/
void get_maxrss(long *maxrss_bytes) {

    long conversion;
    struct rusage usage;

#if defined(__APPLE__) && defined(__MACH__)
    // On macOS, ru_maxrss is already in bytes.
    conversion = 1L;
#elif defined(__linux__)
    // On Linux, ru_maxrss is reported in KiB.
    conversion = 1024L;
#else
    // Other BSDs (FreeBSD, OpenBSD, NetBSD) also report KiB.
    conversion = 1024L;
#endif

    getrusage(RUSAGE_SELF, &usage);
    *maxrss_bytes = usage.ru_maxrss * conversion;
}
