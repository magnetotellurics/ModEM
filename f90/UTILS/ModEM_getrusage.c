#include <sys/resource.h>
#include <stdio.h>

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

    printf("Inside getrusage!\n");

#if  defined(__APPLE_) || defined(__MACH__)
    conversion = 1.0;
#elif __linux__
    conversion = 1.0;
#else  
    conversion = 1.0;
#endif

    getrusage(RUSAGE_SELF, &usage);
    printf("getrusage before conversion: %ld\n", usage.ru_maxrss);
    *maxrss_bytes = usage.ru_maxrss / conversion;    

    printf("Inside getrusage - Maxrss_bytes: %ld - conversion: %d\n", *maxrss_bytes, conversion);
}
