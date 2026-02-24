module ModEM_memory

! ModEM_memory - Call and log maxrss from getrusage (man getrusage)
!
! This module is intended as a diagnostic module. Because of the way ModEM is currently built,
! this module is not intended for release. However, if the build process is eventually upgraded
! so that the fmkmf.pl includes C files (or we do manual makefiles).
!
! To compile and be able to use this module with ModEM, perform the following steps:
!
! 1. Set the environment variable CC to your desired C compiler: `export CC=mpicc`
!
! 2. Add the following compilation steps to your makefile:
!
!   $(OBJDIR)/ModEM_getrusage.o:
!	    $(CC) -c UTILS/ModEM_getrusage.c -o $(OBJDIR)/ModEM_getrusage.o
!   $(OBJDIR)/ModEM_memory.o:UTILS/ModEM_memory.f90  $(OBJDIR)/ModEM_getrusage.o
!	    $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) UTILS/ModEM_memory.f90 -o $(OBJDIR)/ModEM_memory.o
!
! 3. Add both `$(OBJDIR)/ModEM_getrusage.o` and `$(OBJDIR)/ModEM_memory.o` to `OBJ` (the big linking line on line ~56).
!
! 4. You'll manually need to add `$(OBJDIR)/ModEM_getrusage.o` to the prerequisite of any modules/targets
! where you want to call ModEM_memory routines.
!
! 5. Ensure you call ModEM_memory_init and (less critically) ModEM_memory_finalize.
!
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
    public ModEM_memory_get_all

contains

    subroutine ModEM_memory_create_log_fname()

    end subroutine

    subroutine ModEM_memory_get_maxrss(maxrss_bytes)

        implicit none

        integer, intent(out) :: maxrss_bytes
        integer (c_long) :: maxrss_bytes_c

        interface
            subroutine get_maxrss(maxrss_bytes) bind(c)
                use iso_c_binding, only : c_long
                integer (c_long), intent(out) :: maxrss_bytes
            end subroutine get_maxrss
        end interface

        call get_maxrss(maxrss_bytes_c)
        maxrss_bytes = maxrss_bytes_c

    end subroutine ModEM_memory_get_maxrss

    subroutine ModEM_memory_print_report()

        implicit none

        integer :: maxrss
        character (len=512) :: message

        call ModEM_memory_get_maxrss(maxrss)

        write(message, "(A,i4.1,A,i16.1)") "Task: ", task_id, " Max RSS: ", maxrss

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

        write(0,*) "Maxrss_bytes: ", maxrss_bytes, maxrss_kb, maxrss_mb, maxrss_gb

    end subroutine ModEM_memory_convert_maxrss

    subroutine ModEM_memory_log_report(message)

        implicit none

        character (len=*), intent(in) :: message
        character (len=*), parameter :: LOG_MSG_FMT = "(A, A, F18.1, A, F18.1, A, F18.1, A)"
        integer :: maxrss
        real :: maxrss_bytes, maxrss_kb, maxrss_mb, maxrss_gb
        character (len=512) :: log_message

        call ModEM_memory_get_maxrss(maxrss)

        call ModEM_memory_convert_maxrss(maxrss, maxrss_kb, maxrss_mb, maxrss_gb)
        write(log_message, LOG_MSG_FMT) trim(message), ", ", maxrss_kb, ' kb ', maxrss_mb, ' mb ', maxrss_gb, ' gb'
        call ModEM_log(log_message, mainOnly=.false., flush_log=.true.)

    end subroutine ModEM_memory_log_report

    subroutine ModEM_memory_get_all(message)

        implicit none

        character (len=*), intent(in) :: message
        character (len=*), parameter :: LOG_MSG_FMT = "(A, A, F18.1, A, F18.1, A, F18.1, A)"
        character (len=512) :: log_message


        integer :: maxrss, global_maxrss
        real :: maxrss_bytes, maxrss_kb, maxrss_mb, maxrss_gb

        call ModEM_memory_get_maxrss(maxrss)

        call MPI_reduce(maxrss, global_maxrss, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if (taskid == 0) then
            call ModEM_memory_convert_maxrss(global_maxrss, maxrss_kb, maxrss_mb, maxrss_gb)
            write(log_message, LOG_MSG_FMT) trim(message), ", ", maxrss_kb, ' kb ', maxrss_mb, ' mb ', maxrss_gb, ' gb'
            call ModEM_log(log_message, mainOnly=.false., flush_log=.true.)
        end if

    end subroutine ModEM_memory_get_all

end module ModEM_memory
