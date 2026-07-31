module ModEM_timers

! By default, this module will use Fortran's SYSTEM_CLOCK as it's Monotonic; however, it
! has the capability to use C timers. Because ModEM's makefile generation tool doesn't
! automatically compile C files, in order to use the C timers you will need to do the 
! following:
! 
! 1. Compile ModEM_machine_timers.c - By manually editing your makefile, see instructions
!  in ModEM_machine_timers.c on how to do this and makefile rules to use.
! 2. Compile this module with USE_C_TIMERS defined (Add `-DUSE_C_TIMERS` to FFLAGS).
!
    use iso_fortran_env, only : int64
    use iso_c_binding, only : c_int

    implicit none

    private

    type ModEM_timer_t
        character(len=:), allocatable :: name
        integer (c_int) :: id

        logical :: is_running = .false.

        integer :: hours = 0
        integer :: mins = 0
        integer :: secs = 0
        integer(kind=int64) :: nsecs = 0

        integer (c_int) :: secs_c = 0
        integer (c_int) :: nsecs_c = 0
        integer (kind=int64) :: total_secs = 0
        integer (kind=int64) :: total_nsecs = 0

        type (ModEM_timer_t), pointer :: next => null()
    end type ModEM_timer_t

    type (ModEM_timer_t) :: timers

    integer, private :: id_counter = 0
    integer, parameter :: NSECS_TO_SECS = 1000000000

#ifndef USE_C_TIMERS
    integer, private, parameter :: MAX_FORTRAN_TIMERS = 100
    integer (int64), private :: fortran_clock_rate
    integer (int64), private, dimension(MAX_FORTRAN_TIMERS) :: fortran_timers_start
    integer (int64), private, dimension(MAX_FORTRAN_TIMERS) :: fortran_timers_end
#endif

    public ModEM_timer_t

    public ModEM_timers_create, ModEM_timers_destroy, ModEM_timers_destroy_all
    public ModEM_timers_start, ModEM_timers_stop, ModEM_timers_stop_all
    public ModEM_timers_exist, ModEM_timers_get
    public ModEM_timers_print, ModEM_timers_print_all, ModEM_timers_report

    ! The following will not be public in ModEM (only public for unit testing)
    public accumulate_time, convert_time

    contains

    function ModEM_timers_exist(timer_name) result(exists)

        implicit none

        character (len=*), intent(in) :: timer_name
        type (ModEM_timer_t), pointer :: timer => null()
        logical :: exists

        exists = .false.

        timer => ModEM_timers_get(timer_name)
        if (associated(timer)) then
            exists = .true.
        end if

    end function ModEM_timers_exist

    function get_new_timer_id() result(id)

        implicit none

        integer (c_int) :: id

        id_counter = id_counter + 1
        id = id_counter 

    end function get_new_timer_id

    subroutine ModEM_timers_create(timer_name, start_now)

        implicit none

        character (len=*), intent(in) :: timer_name
        logical, intent(in), optional :: start_now
        logical :: start_timer_now

        type (ModEM_timer_t), pointer :: new_timer

        if (present(start_now)) then
            start_timer_now = start_now
        else
            start_timer_now = .false.
        end if

        ! Don't allow timers of the same name
        if (ModEM_timers_exist(timer_name)) then
            return
        end if

        allocate(new_timer)

        new_timer % id = get_new_timer_id()
        new_timer % name = trim(timer_name)
        new_timer % next => timers % next
        timers % next => new_timer

        if (start_timer_now) then 
            call ModEM_timers_start(timer_name)
        end if

    end subroutine ModEM_timers_create

    function ModEM_timers_get(timer_name) result(timer)

        implicit none

        character (len=*), intent(in) :: timer_name
        type (ModEM_timer_t), pointer :: cur, timer

        timer => null()
        cur => timers % next 

        do while (associated(cur))
            if (cur % name == timer_name) then
                timer => cur
                exit
            end if 

            cur => cur % next
        end do

    end function ModEM_timers_get

    subroutine ModEM_timers_start_machine(timer)

        implicit none

#ifdef USE_C_TIMERS
        interface
            subroutine timer_start(timer_id) bind(C)
               use iso_c_binding, only : c_int
               integer (c_int), intent(in), value :: timer_id
            end subroutine timer_start
        end interface
#endif

        type (ModEM_timer_t) :: timer

#ifdef USE_C_TIMERS
        call timer_start(timer % id)
#else
        call system_clock(count=fortran_timers_start(timer % id))
#endif

    end subroutine ModEM_timers_start_machine

    subroutine ModEM_timers_start(timer_name)

        implicit none

        character (len=*), intent(in) :: timer_name

        type (ModEM_timer_t), pointer :: timer

        ! Create a timer if it doesn't exist
        if (ModEM_timers_exist(timer_name)) then
            timer => ModEM_timers_get(timer_name)
        else
            call ModEM_timers_create(timer_name)
            timer => ModEM_timers_get(timer_name)
        end if

        if (timer % is_running) then
            return
        end if

        timer % is_running = .true.
        call ModEM_timers_start_machine(timer)

    end subroutine ModEM_timers_start

    subroutine ModEM_timers_stop_machine(timer)

        implicit none

#ifdef USE_C_TIMERS
        interface
            subroutine timer_stop(timer_id, sec, nsec) bind(C)
               use iso_c_binding, only : c_int
               integer (c_int), intent(in), value :: timer_id
               integer (c_int), intent(out) :: sec, nsec
            end subroutine timer_stop
        end interface
#endif
 
        type (ModEM_timer_t), intent(inout) :: timer

#ifndef USE_C_TIMERS
        real :: seconds
#endif

#ifdef USE_C_TIMERS
        call timer_stop(timer % id, timer % secs_c, timer % nsecs_c)
#else
        call system_clock(count=fortran_timers_end(timer % id))
        call system_clock(count_rate=fortran_clock_rate)

        seconds = real(fortran_timers_end(timer % id) - fortran_timers_start(timer % id)) / real(fortran_clock_rate)
        timer % secs_c = int(seconds)
        timer % nsecs_c = (seconds - int(seconds)) * NSECS_TO_SECS
#endif

    end subroutine ModEM_timers_stop_machine

    subroutine ModEM_timers_stop(timer_name, accumulate)

        implicit none

        character (len=*), intent(in) :: timer_name
        logical, optional :: accumulate

        type (ModEM_timer_t), pointer :: timer
        logical :: accumulate_time_bool

        if (.not. ModEM_timers_exist(timer_name)) then
            return
        end if

        if (present(accumulate)) then
            accumulate_time_bool = accumulate
        else
            accumulate_time_bool = .true.
        end if

        timer => ModEM_timers_get(timer_name)

        if (.not. timer % is_running) then
            return
        end if

        timer % is_running = .false.

        call ModEM_timers_stop_machine(timer)

        if (accumulate_time_bool) then
            call accumulate_time(timer)
        else
            timer % total_secs = timer % secs_c
            timer % total_nsecs = timer % nsecs_c
        end if

        call convert_time(timer)

    end subroutine ModEM_timers_stop

    subroutine ModEM_timers_stop_all()

        implicit none

        type (ModEM_timer_t), pointer :: cur

        cur => timers % next

        do while (associated(cur))
            call ModEM_timers_stop(cur % name)
            cur => cur % next
        end do

    end subroutine ModEM_timers_stop_all

    subroutine ModEM_timers_destroy(timer_name)

        implicit none

        character (len=*), intent(in) :: timer_name
        type (ModEM_timer_t), pointer :: cur, prev

        if (.not. associated(timers % next)) then
            return
        end if

        cur => timers % next
        prev => timers % next

        do while (associated(cur))
            if (cur % name == timer_name) then

                if (associated(timers % next, cur)) then
                    ! Handle the first element
                    timers % next => cur % next
                else
                    prev % next => cur % next
                end if

                deallocate(cur % name)
                deallocate(cur)
                return
            end if
            
            prev => cur
            cur => cur % next
        end do

    end subroutine ModEM_timers_destroy

    subroutine ModEM_timers_destroy_all()

        implicit none

        type (ModEM_timer_t), pointer :: cur

        cur => timers % next

        do while (associated(cur))
            call ModEM_timers_destroy(cur % name)
            cur => cur % next
        end do

    end subroutine ModEM_timers_destroy_all

    subroutine ModEM_timers_print(timer_name, file)

        implicit none

        character (len=*), intent(in) :: timer_name
        integer, optional :: file
        integer :: file_descriptor

        type(ModEM_timer_t), pointer :: timer
        integer :: hours, mins, secs, nsecs
        character(len=80) :: ModEM_timer_str_format = "(A, A30, A, A, A, I2.2, A, I2.2, A, I2.2, A, I0.8)"

        if (present(file)) then
            file_descriptor = file
        else
            file_descriptor = 0
        end if

        ! "Timer 'Timer Name' - Elapsed Time: 00:00:00"
        timer => ModEM_timers_get(timer_name)

        hours = timer % hours
        mins = timer % mins
        secs = timer % secs
        nsecs = timer % nsecs

        write(file_descriptor, ModEM_timer_str_format) "Timer: '", trim(timer % name), "'", achar(9) // achar(9), " Elapsed Time: ", &
            hours, ":", mins, ":", secs, ":", nsecs

    end subroutine ModEM_timers_print

    subroutine ModEM_timers_print_all(file)

        implicit none

        integer, optional :: file
        integer :: file_descriptor
        type(ModEM_timer_t), pointer :: cur

        if (present(file)) then
            file_descriptor = file
        else
            file_descriptor = 0
        end if

        cur => timers % next

        do while (associated(cur))
            call ModEM_timers_print(cur % name, file_descriptor)
            cur => cur % next
        end do

    end subroutine ModEM_timers_print_all

    subroutine ModEM_timers_report(file)

        implicit none

        integer, intent(in), optional :: file
        integer :: file_descriptor

        if (present(file)) then
            file_descriptor = file
        else
            file_descriptor = 0
        end if

        write(file_descriptor,*) "============= ModEM Run Report ============"

        call ModEM_timers_print_all(file_descriptor)

        write(file_descriptor,*) "============================"

    end subroutine ModEM_timers_report

    subroutine accumulate_time(timer)

        implicit none

        type (ModEM_timer_t), pointer :: timer

        integer(c_int) :: remainder
        integer(kind=int64) :: converted_nsecs

        timer % total_nsecs = timer % total_nsecs + timer % nsecs_c

        converted_nsecs = timer % total_nsecs / NSECS_TO_SECS
        timer % total_secs = timer % total_secs + timer % secs_c + (converted_nsecs)
        timer % total_nsecs = timer % total_nsecs - (converted_nsecs * NSECS_TO_SECS)

    end subroutine accumulate_time

    subroutine convert_time(timer)

        implicit none

        type(ModEM_timer_t), intent(inout) :: timer

        timer % hours = timer % total_secs / 3600
        timer % mins = (timer % total_secs -(3600 * timer % hours)) / 60
        timer % secs = timer % total_secs -(3600 * timer % hours) - (timer % mins * 60) 
        timer % nsecs = timer % total_nsecs 

    end subroutine convert_time

end module ModEM_timers
