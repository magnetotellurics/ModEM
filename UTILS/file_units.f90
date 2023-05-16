! *****************************************************************************
module file_units
  ! This module defines I/O file unit information, used throughout the code

  implicit none

  ! User-defined parameter that sets the level of output
  integer, save									:: output_level=0
  character(12), save                           :: node_info=''

  ! Startup and control files
  integer, parameter							:: ioStartup=101
  integer, parameter							:: ioFwdCtrl=105
  integer, parameter							:: ioInvCtrl=106

  ! Log file
  integer, parameter							:: ioLog=111

  ! Weighting Log file
   integer, parameter							:: ioWeighting=222  
   
  ! Weighting Log file
   integer, parameter							:: ioConfig=999 

  ! MPI Log file
  integer, parameter                            :: ioMPI=2000

  ! Grid, model, data
  integer, parameter							:: ioPrm=1
  integer, parameter							:: ioGrd=2
  integer, parameter							:: ioDat=3
  integer, parameter							:: ioMdl=4

  ! Generic error, read and write units
  integer, parameter							:: ioERR=7
  integer, parameter							:: ioREAD=8
  integer, parameter							:: ioWRITE=9

  ! Dictionary files, if they exist
  integer, parameter							:: ioRX=17
  integer, parameter							:: ioTX=18
  integer, parameter							:: ioDT=19

  ! Fields, currents etc
  integer, parameter	  				        :: ioH=10
  integer, parameter	  				        :: ioE=11
  integer, parameter	  				        :: ioJ=12

  ! Fields in x,y,z directions
  integer, parameter	  				        :: ioHx=20
  integer, parameter	  				        :: ioHy=30
  integer, parameter	  				        :: ioHz=40
  integer, parameter	  				        :: ioEx=21
  integer, parameter	  				        :: ioEy=31
  integer, parameter	  				        :: ioEz=41

  ! Others, as necessary
  integer, parameter					        :: ioSens=27

  integer, parameter							:: ioShell=28
  integer, parameter							:: ioBound=29
  integer, parameter							:: ioRad=15
  integer, parameter							:: ioResp=25
  integer, parameter                            :: ioC=23, ioD=24
  integer, parameter                            :: ioAvgC=33, ioAvgD=34
  
  integer, parameter                            :: ioSolverStat=50

end module file_units
