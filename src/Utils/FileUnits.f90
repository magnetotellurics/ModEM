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
    integer, parameter :: ioFwdTmp = 1115
    !
    integer, parameter :: ioCovariance = 1116
    !
    integer, parameter :: ioInvLog = 1117
    !
    integer, parameter :: ioInvTmp = 1118
    !
    integer, parameter :: ioGradLog = 1119
    !
    integer, parameter :: ioGradNorm = 1120
    !
    integer, parameter :: ioGradRMS = 1121
    !
    integer, parameter :: ioPlot = 1122
    !
end module FileUnits
