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
*   subroutine get_maxrss() result(maxrss) bind(c)
*       use iso_c_binding, only : c_long
*       integer (c_long), intent(out) :: maxrss
*   end subroutine get_maxrss
* end interface
*/
void get_maxrss(long *maxrss_bytes) {

    int conversion;
    struct rusage usage;

#if  defined(__APPLE_) || defined(__MACH__)
    // BSD's ru_maxrss (used by Apple) is in bytes
    conversion = 1.0;
#elif __linux__
    // Linux's ru_maxrss is in KiB
    conversion = 1000.0;
#else
    conversion = 1.0;
#endif

    getrusage(RUSAGE_SELF, &usage);
    *maxrss_bytes = usage.ru_maxrss / conversion;
}
