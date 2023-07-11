!
!> defines I/O file unit information, used throughout the code
!> 
!
module FileUnits
    !
    implicit none
    !
    integer, parameter :: ioStartup = 1
    !
    integer, parameter :: ioPredData = 2
    !
    integer, parameter :: ioESolution = 3
    !
    integer, parameter :: ioModelParam = 4
    !
    integer, parameter :: ioFwdTmp = 5
    !
    integer, parameter :: ioCovariance = 6
    !
    integer, parameter :: ioInvLog = 7
    !
    integer, parameter :: ioInvPlot = 8
    !
    integer, parameter :: ioFuncPlot = 9
    !
    integer, parameter :: ioInvTmp = 10
    !
    integer, parameter :: ioGradLog = 11
    !
    integer, parameter :: ioGradNorm = 12
    !
    integer, parameter :: ioGradRMS = 13
    !
    integer, parameter :: ioPlot = 14
    !
end module FileUnits
