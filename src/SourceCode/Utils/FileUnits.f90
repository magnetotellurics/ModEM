!
!> defines I/O file unit information, used throughout the code
!> 
!
module FileUnits
    !
    implicit none
    !
    integer, parameter :: ioStartup = 10018
    !
    integer, parameter :: ioPredData = 10019
    !
    integer, parameter :: ioESolution = 10020
    !
    integer, parameter :: ioModelParam = 10021
    !
    integer, parameter :: ioFwdTmp = 10022
    !
    integer, parameter :: ioCovariance = 10023
    !
    integer, parameter :: ioInvLog = 10024
    !
    integer, parameter :: ioInvTmp = 10025
    !
    integer, parameter :: ioGradLog = 10026
    !
    integer, parameter :: ioGradNorm = 10027
    !
    integer, parameter :: ioGradRMS = 10028
    !
    integer, parameter :: ioPlot = 10029
    !
end module FileUnits
