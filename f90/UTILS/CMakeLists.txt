# f90/UTILS

target_sources(${MODEM_EXE}
    PRIVATE elements.f90
    PRIVATE fields_orientation.f90
    PRIVATE file_units.f90
    PRIVATE math_constants.f90
    PRIVATE ModEM_timers.f90
    PRIVATE polpak.f90
    PRIVATE utilities.f90
    PRIVATE ModEM_logger.f90
)

if (USE_C_TIMERS)
    target_sources(${MODEM_EXE}
        PRIVATE ModEM_machine_timers.c
    )

    add_compile_definitions(USE_C_TIMERS)
endif()
