!
!> defines I/O file unit information, used throughout the code
!> 
!
module FileUnits
    !
    implicit none
    !
    integer, parameter :: ioStartup = 1111
    !
    integer, parameter :: ioPredData = 1112
    !
    integer, parameter :: ioESolution = 1113
    !
    integer, parameter :: ioModelParam = 1114
    !
    integer, parameter :: ioInvLog = 1115
    !
    integer, parameter :: ioFwdTmp = 1116
    !
    integer, parameter :: ioInvTmp = 1117
    !
end module FileUnits
