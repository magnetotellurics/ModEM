! module ModEM_memory 
!
! Use getrusage to track memory high-water mark.
!
! This module contains helpful routines to call and
! log the results of ru_maxrss by calling the 
! getrusage system call.
!
! The public routines print ru_maxrss, which is
! a high-water mark of memory. More information on
! getrusage can be found by in `man 2 getrusage`.
!

module ModEM_memory

    use utilities
    use iso_c_binding, only : c_long
#ifdef MPI
    use mpi
#endif
    use ModEM_logger

    implicit none

    private

    integer :: task_id

    public ModEM_memory_print_report
    public ModEM_memory_log_report
#ifdef MPI
    public ModEM_memory_get_all
#endif

contains

    subroutine ModEM_memory_get_maxrss(maxrss_bytes)

        implicit none

        integer, intent(out) :: maxrss_bytes
        integer (c_long) :: maxrss_bytes_c

#ifdef MODEM_MEMORY
        interface
            subroutine get_maxrss(maxrss_bytes) bind(c)
                use iso_c_binding, only : c_long
                integer (c_long), intent(out) :: maxrss_bytes
            end subroutine get_maxrss
        end interface

        call get_maxrss(maxrss_bytes_c)
        maxrss_bytes = maxrss_bytes_c
#else
        write(0,*) "ERROR: ModEM was not built with the ModEM_getrusage.c file!"
        write(0,*) "ERROR: In order to use ModEM_memory routines you must"
        write(0,*) "ERROR: build ModEM with CMake and specify '-DMODEM_MEMORY=on'"
        call ModEM_abort()
#endif

    end subroutine ModEM_memory_get_maxrss

    subroutine ModEM_memory_print_report(message)

        implicit none

        integer :: maxrss
        character (len=*), intent(in) :: message
        real :: maxrss_bytes, maxrss_kb, maxrss_mb, maxrss_gb
        character (len=*), parameter :: LOG_MSG_FMT = "(I4.4, A, A, A, F18.1, A, F18.1, A, F18.1, A)"
        character (len=512) :: log_message

        call ModEM_memory_get_maxrss(maxrss)
        call ModEM_memory_convert_maxrss(maxrss, maxrss_kb, maxrss_mb, maxrss_gb)

        write(log_message, LOG_MSG_FMT) task_id, ' - ', trim(message), ", ", maxrss_kb, ' kb ', maxrss_mb, ' mb ', maxrss_gb, ' gb'
        write(6,*) trim(log_message)

    end subroutine ModEM_memory_print_report

    subroutine ModEM_memory_convert_maxrss(maxrss_bytes, maxrss_kb, maxrss_mb, maxrss_gb)

        implicit none

        integer :: maxrss_bytes
        real, intent(out) :: maxrss_kb
        real, intent(out) :: maxrss_mb
        real, intent(out) :: maxrss_gb

        maxrss_kb = maxrss_bytes / 1.0
        maxrss_mb = maxrss_kb / 1000.0
        maxrss_gb = maxrss_mb / 1000.1

    end subroutine ModEM_memory_convert_maxrss

    subroutine ModEM_memory_log_report(message)

        implicit none

        character (len=*), intent(in) :: message
        character (len=*), parameter :: LOG_MSG_FMT = "(A, A, F18.1, A, F18.1, A, F18.1, A)"
        integer :: maxrss
        real :: maxrss_kb, maxrss_mb, maxrss_gb
        character (len=512) :: log_message

        call ModEM_memory_get_maxrss(maxrss)

        call ModEM_memory_convert_maxrss(maxrss, maxrss_kb, maxrss_mb, maxrss_gb)
        write(log_message, LOG_MSG_FMT) trim(message), ", ", maxrss_kb, ' kb ', maxrss_mb, ' mb ', maxrss_gb, ' gb'
        call ModEM_log(log_message, mainOnly=.false., flush_log=.true.)

    end subroutine ModEM_memory_log_report

#ifdef MPI
    subroutine ModEM_memory_get_all(message)

        implicit none

        character (len=*), intent(in) :: message
        character (len=*), parameter :: LOG_MSG_FMT = "(A, A, F18.1, A, F18.1, A, F18.1, A)"
        character (len=512) :: log_message

        integer :: maxrss, global_maxrss
        real :: maxrss_kb, maxrss_mb, maxrss_gb

        call ModEM_memory_get_maxrss(maxrss)

        call MPI_reduce(maxrss, global_maxrss, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if (taskid == 0) then
            call ModEM_memory_convert_maxrss(global_maxrss, maxrss_kb, maxrss_mb, maxrss_gb)
            write(log_message, LOG_MSG_FMT) trim(message), ", ", maxrss_kb, ' kb ', maxrss_mb, ' mb ', maxrss_gb, ' gb'
            call ModEM_log(log_message, mainOnly=.false., flush_log=.true.)
        end if

    end subroutine ModEM_memory_get_all
#endif

end module ModEM_memory
