! *****************************************************************************
module iospec
  ! for i/o, the file names have numbers; initialize them here before using

  use file_units
  implicit none

  !integer, parameter	  				:: ioHx=8, ioHy=9, ioHz=10
  !integer, parameter					:: ioJ=11, ioH=12, ioERR=7, ioC=13, ioD=14
  !integer, parameter					:: ioAvgC=33, ioAvgD=34, ioBV=2, ioJac=25
  !integer, parameter					:: ioMdl=20, ioResp=21, ioOut=22, ioEarth=23

  logical								:: opened ! for I/O inquiries

  ! indicates the required level of output
  ! integer, save							:: output_level

Contains

  ! ***************************************************************************
  ! * initFileWrite opens the chosen output file and writes the relevant header
  subroutine initFileWrite(fName, ioNum)

    implicit none
	character(80), intent(in)			:: fName
    integer, intent(in)					:: ioNum
    integer								:: ios

    open(unit=ioNum,file=fName,status='unknown',iostat=ios)

    if (ios/=0) then
       write(0,*) 'Error opening file in initFileWrite: ', fName
       stop
    endif

    ! write the header
	if (ioNum .eq. ioHx) then
	  write(ioNum,*)' phi theta r Re(Hx) Im(Hx)'
	else if (ioNum .eq. ioHy) then
	  write(ioNum,*)' phi theta r Re(Hy) Im(Hy)'
	else if (ioNum .eq. ioHz) then
	  write(ioNum,*)' phi theta r Re(Hz) Im(Hz)'
    else if (ioNum .eq. ioJ) then
	  write(ioNum,'(a)')' phi theta r'// &
		' Re(Hx) Im(Hx) Re(Hy) Im(Hy) Re(Hz) Im(Hz)'// &
		' Re(Ex) Im(Ex) Re(Ey) Im(Ey) Re(Ez) Im(Ez)'
	else if (ioNum .eq. ioERR) then
	  write(ioNum,*)' icount erre herr divH'
	end if

  end subroutine initFileWrite ! initFileWrite

  ! ***************************************************************************
  ! * closeOutFiles is called by the main program before exiting completely
  subroutine closeOutFiles()

	close(ioC); close(ioD); close(ioH); close(ioJ)
	close(ioHx); close(ioHy); close(ioHz)
	close(ioAvgC); close(ioAvgD)
	close(ioERR); close(ioWRITE)

  end subroutine closeOutFiles	! closeOutFiles


end module iospec
