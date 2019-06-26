
! **********************************************************************************

! Module for basic input and output of WS

module wsio
  use math_constants
  use modelparameter
  use soln2d
  use dataspace
  use emsolver
use utilities

  implicit none

  ! public routine

  public   ::  readModel


contains

! **********************************************************************************

  subroutine readModel(fid,cfile,grid,sigma)
 
  integer,           intent(in)    :: fid
  character*80,      intent(in)    :: cfile
  type(grid2d_t),   intent(inout) :: grid
  type(earthCond2D), intent(inout) :: sigma

  ! local variable 
  integer            ::  Ny, Nzair, Nzearth, Nz
  character*80       ::  paramType
  character*80       ::  ctmp
  character*256      ::  cindex
  logical            ::  endOFfile
  integer            ::  ifile, bwrd,ewrd ,iy,iz, len, iyt, nyt
  integer            ::  irmax,ir

  real(kind=8), allocatable, dimension(:)   :: Dy, Dz
  real(kind=8), allocatable, dimension(:,:) :: cond 
  integer, allocatable, dimension(:,:) :: irho
  real(kind=8), allocatable, dimension(:) :: rhoval


  endOFfile = .false.

  open(unit=fid,file=cfile,form='formatted',status='OLD')

  ifile = 0
  do
    read(fid,'(a80)') ctmp

!   write(6,*) ifile, ctmp

    if (FindStr(ctmp(1:10),'#') > 0) then
    endif
    if (FindStr(ctmp(1:10),'          ') > 0) then
    endif

    if (FindStr(ctmp,'title') > 0) then
       ifile = ifile + 1
    endif

    if (FindStr(ctmp,'ny') > 0) then
       ifile = ifile + 1
       bwrd = BegWrd(ctmp,2)
       ewrd = EndWrd(ctmp,2)
       read(ctmp(bwrd:ewrd),*) Ny

       allocate(Dy(Ny))
       read(fid,*) (Dy(iy), iy=1,Ny)
    endif

    if (FindStr(ctmp,'nzb') > 0) then
       ifile = ifile + 1
       bwrd = BegWrd(ctmp,2)
       ewrd = EndWrd(ctmp,2)
       read(ctmp(bwrd:ewrd),*) Nzearth

       Nzair = 10
       Nz = Nzair + Nzearth
       allocate(Dz(Nz))
 
       Dz(1) = 300000.
       Dz(2) = 100000.
       Dz(3) = 30000.
       Dz(4) = 10000.
       Dz(5) = 3000.
       Dz(6) = 1000.
       Dz(7) = 300.
       Dz(8) = 100.
       Dz(9) = 30.
       Dz(10)= 10.

       read(fid,*) (Dz(iz), iz=Nzair+1,Nz)
    endif

    if (FindStr(ctmp,'resistivity_model') > 0) then
       ifile = ifile + 1

       allocate(irho(Ny,Nzearth))

       do iz = 1,Nzearth
          nyt = 0
          read(fid,'(a256)') cindex
          bwrd = BegWrd(cindex,1)
          ewrd = EndWrd(cindex,1)
          len = ewrd - bwrd + 1
          iyt = 0
          do iy = bwrd,ewrd
	    iyt = iyt + 1
            read(cindex(iy:iy),'(i1)') irho(iyt+nyt,iz)
          enddo
          nyt = len+nyt
       enddo
   
       irmax = 0
       do iz = 1,Nzearth
          do iy = 1,Ny
            irmax = MAX(irho(iy,iz),irmax)
          enddo
       enddo

       allocate(rhoval(irmax))
       read(fid,*) (rhoval(ir),ir=1,irmax)

       allocate(cond(Ny,Nzearth))
       do iz = 1,Nzearth
          do iy = 1,Ny
            cond(iy,iz) = DLOG(1./rhoval(irho(iy,iz)))
!           write(14,*) iy,iz, cond(iy,iz), irho(iy,iz)
          enddo
       enddo
    endif


    if (ifile == 4) endOffile = .true.

    if (endOFfile) exit

  enddo 

  close(fid)


! for grid type 

  grid%Ny  = Ny
  grid%Nza = Nzair
  grid%Nz  = Nzearth + Nzair

  call allocateGrid2D(grid)
  
  grid%Dy = Dy
  grid%Dz = Dz


! for earthCond2D type
  call allocateEarthCond(grid,sigma)
  sigma%v = cond
  sigma%AirCond = SIGMA_AIR
  sigma%paramType = 'LOGE'


!  write(6,*) grid%Ny,grid%Nz
!  write(6,*) (grid%Dy(iy),iy=1,grid%Ny) 
!  write(6,*) (grid%Dz(iz),iz=1,grid%Nz)
!  do iz = 1,Nzearth
!     write(72,*) (irho(iy,iz),iy=1,Ny)
!  enddo
!  write(6,*) (rhoval(ir),ir=1,irmax)
   do iy = 1,Ny
     do iz = 1,Nzearth
       write(72,*) iy,iz,sigma%v(iy,iz)
     enddo
   enddo
   write(6,*) 'Air cond = ',sigma%AirCond, sigma%paramType

  end subroutine readModel


! **********************************************************************************


  subroutine readData(fid,cfile,im,nComp,nPer,nSites,periods, sites, allData)

  integer,           intent(in)    :: fid,im
  integer,           intent(in)    :: nComp
  character*80,      intent(in)    :: cfile

  integer, intent(inout)           :: nPer, nSites
  real (kind=8), intent(inout), allocatable, dimension(:)   :: periods
  real (kind=8), intent(inout), allocatable, dimension(:,:) :: sites
  type (dvecMultTX), intent(inout) :: allData


  ! local variables
  integer            :: iPer, ifile, nRes
  character*80       :: ctmp, cres
  integer            :: bwrd,ewrd ,iy,iz, len, iyt, nyt, is, ip, ic
  integer            :: irmax,ir
  logical            ::  endOFfile

  real (kind=8)      :: per, errtmp,zr,zi, absz
  real (kind=8), allocatable, dimension(:) :: dattemp

  integer Ndata, k , iDt

!!!!!!!!!!!!!!!!!!!!! read data file start here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  endOFfile = .false.

  ifile = 0
  open(unit=fid,file=cfile,form='formatted',status='OLD')
  
  do 
    read(fid,'(a80)') ctmp

  ! write(6,*) ifile, ctmp

    if (FindStr(ctmp(1:10),'#') > 0) then
    endif
    if (FindStr(ctmp(1:10),'          ') > 0) then
    endif

    if (FindStr(ctmp,'title') > 0) then
       ifile = ifile + 1
    endif

    if (FindStr(ctmp,'mode_type') > 0) then
       ifile = ifile + 1
    endif

    if (FindStr(ctmp,'number_of_response') > 0) then
       ifile = ifile + 1
       bwrd = BegWrd(ctmp,2)
       ewrd = EndWrd(ctmp,2)
       read(ctmp(bwrd:ewrd),*) nRes
       write(6,*) 'nRes = ',nRes,im

       if (nRes .eq. 0) then
          write(6,*) 'ERROR in nRes '
          stop
       endif
    endif

    if (FindStr(ctmp,'number_of_period') > 0) then
       ifile = ifile + 1
       bwrd = BegWrd(ctmp,2)
       ewrd = EndWrd(ctmp,2)
       read(ctmp(bwrd:ewrd),*) nPer

       if (im.eq.1) allocate(periods(nPer))
       read(fid,*) (periods(iPer), iPer=1,nPer)

       allocate(allData%d(nPer))
       allData%allocated = .true.
       allData%errorBar  = .true.
       allData%nTx = nPer
    endif

    if (FindStr(ctmp,'number_of_station') > 0) then
       ifile = ifile + 1
       bwrd = BegWrd(ctmp,2)
       ewrd = EndWrd(ctmp,2)
       read(ctmp(bwrd:ewrd),*) nSites

       if (im.eq.1) allocate(sites(2,nSites))
       read(fid,*) (sites(1,is), is=1,nSites)

       do is = 1,nSites
         sites(2,is) = 444440.
       enddo

       Ndata = 0
       do ip = 1,nPer
         allData%d(ip)%allocated = .true.
       ! allData%d(ip)%errorBar  = .true.
       ! allData%d(ip)%tx        = ip

         iDt = 1
         call createDvec(nComp,nSites,allData%d(ip),.true.,ip, iDt)
         Ndata = Ndata + nComp*nSites

       ! if (ip.eq.1) then
           do k = 1,nSites
             allData%d(ip)%rx(k) = k
           enddo
       ! endif

       enddo

       allData%Ndata = Ndata
       allocate(dattemp(nSites))

       write(6,*) 'Ndata = ',Ndata
    endif


    if (FindStr(ctmp,'data_response') > 0) then
       ifile = ifile + 1
       bwrd = BegWrd(ctmp,2)
       ewrd = EndWrd(ctmp,2)
       cres = ctmp(bwrd:ewrd)


       if (FindStr(cres,'realZ') > 0) then
         do ip = 1,nPer
	   read(fid,*) per, (dattemp(is),is=1,nSites)
           do is = 1,nSites
	     allData%d(ip)%data(1,is) = dattemp(is)
           enddo
         enddo
       endif 

       if (FindStr(cres,'imagZ') > 0) then
         do ip = 1,nPer
	   read(fid,*) per, (dattemp(is),is=1,nSites)
           do is = 1,nSites
	     allData%d(ip)%data(2,is) = dattemp(is)
           enddo
         enddo
       endif 
    endif

    if (FindStr(ctmp,'error_response') > 0) then
       ifile = ifile + 1
       bwrd = BegWrd(ctmp,2)
       ewrd = EndWrd(ctmp,2)
       cres = ctmp(bwrd:ewrd)

       if (FindStr(cres,'realZ') > 0) then
         do ip = 1,nPer
	   read(fid,*) per, (dattemp(is),is=1,nSites)
           do is = 1,nSites
	     allData%d(ip)%err(1,is) = dattemp(is)
           enddo
         enddo
       endif 

       if (FindStr(cres,'imagZ') > 0) then
         do ip = 1,nPer
	   read(fid,*) per, (dattemp(is),is=1,nSites)
           do is = 1,nSites
	     allData%d(ip)%err(2,is) = dattemp(is)
           enddo
         enddo
       endif 

    endif


    if (ifile == 8) endOffile = .true.

    if (endOFfile) exit

  enddo 
  close(fid)

  deallocate(dattemp)

!  extra  add 5% error bar to "synthetic data"

  do ip = 1,nPer
    do is = 1,nSites
      zr = allData%d(ip)%data(1,is)
      zi = allData%d(ip)%data(2,is)
      errtmp = dsqrt(abs(zr*zr + zi*zi))*0.05
      allData%d(ip)%err(1,is) = errtmp
      allData%d(ip)%err(2,is) = errtmp
    enddo
  enddo

! do ip = 1,nPer
!  do is = 1,nSites
!   do ic = 1,2
!    write(49,*) allData%d(ip)%err(ic,is)
!   enddo
!   enddo
! enddo


! write(6,*) 'nPer = ',nPer
! write(6,*) 'period ',periods
! write(6,*) 'nSite = ',nSites
! write(6,*) 'sites = ',(sites(1,is),is=1,nSites)

! write(13,*) 'im = ',im,' ---------------------------------'
! do ip = 1,nPer
!   write(13,*) ' data_response ' 
!   write(13,*)(allData%d(ip)%data(1,1,is),is=1,nSites+1)
!   write(13,*) ' data_response '
!   write(13,*)(allData%d(ip)%data(1,2,is),is=1,nSites+1)

!   write(13,*) ' err_response ' 
!   write(13,*)(allData%d(ip)%err(1,1,is),is=1,nSites+1)
!   write(13,*) ' err_response '
!   write(13,*)(allData%d(ip)%err(1,2,is),is=1,nSites+1)
! enddo  

  allData%nTx = nPer

  end subroutine readData


! **********************************************************************************


  subroutine writeData(fid,cfile,im,nPer,nSites,periods, sites, dat)

  integer,           intent(in)    :: fid, im
  character*80,      intent(in)    :: cfile

  integer, intent(inout)           :: nPer, nSites
  real (kind=8), intent(inout), allocatable, dimension(:) :: periods
  real (kind=8), intent(inout), allocatable, dimension(:,:) :: sites
  type (dvecMultTX), intent(inout) :: dat


  ! local variables
  integer            :: iPer, ifile, nRes
  character*80       :: ctmp, cres
  integer            :: bwrd,ewrd ,iy,iz, len, iyt, nyt, is, ip
  integer            :: irmax,ir
  logical            ::  endOFfile

  integer            :: ic, nrow, nr2, iss, isa, isb

  real (kind=8)      :: seven


  seven = 7.0

  ifile = 0
  open(unit=fid,file=cfile,form='formatted')
  
  nr2 = MOD(nSites, 7)
  if (nr2 .eq. 0) then
    nrow = idint(nSites/seven)
  else
    nrow = idint(nSites/seven) + 1
  endif

  write(fid,*) 'TITLE                TEST MODEL'
  if (im.eq.1) write(fid,*) 'MODE_TYPE            tm'
  if (im.eq.2) write(fid,*) 'MODE_TYPE            te'
  write(fid,*) 'NUMBER_OF_RESPONSE   2'
  write(fid,*) 'NUMBER_OF_PERIOD', nPer
  write(fid,100) (periods(iPer), iPer=1,nPer)
  write(fid,*) 'NUMBER_OF_STATION', nSites
  write(fid,100) (sites(1,is), is=1,nSites)

  write(fid,*) 'DATA_RESPONSE        realZ' 
  do ip = 1,nPer
    iss = 0
    do ic = 1, nrow
      isa = iss + 1
      isb = iss + 7

      if (isb >= nSites) isb = nSites
      if (ic .eq. 1) then
        write(fid,200) periods(ip),(dat%d(ip)%data(1,is),is=isa,isb)
      else
        write(fid,300) ' ',(dat%d(ip)%data(1,is),is=isa,isb)
      endif

      iss = iss + 7
    enddo
  enddo

  write(fid,*) 'ERROR_RESPONSE        realZ' 
  do ip = 1,nPer
    iss = 0
    do ic = 1, nrow
      isa = iss + 1
      isb = iss + 7
  
      if (isb >= nSites) isb = nSites
      if (ic .eq. 1) then
        write(fid,200) periods(ip),(dat%d(ip)%err(1,is),is=isa,isb)
      else
        write(fid,300) ' ',(dat%d(ip)%err(1,is),is=isa,isb)
      endif

      iss = iss + 7
    enddo
  enddo

  write(fid,*) 'DATA_RESPONSE        imagZ' 
  do ip = 1,nPer
    iss = 0
    do ic = 1, nrow
      isa = iss + 1
      isb = iss + 7
  
      if (isb >= nSites) isb = nSites
      if (ic .eq. 1) then
        write(fid,200) periods(ip),(dat%d(ip)%data(2,is),is=isa,isb)
      else
        write(fid,300) ' ',(dat%d(ip)%data(2,is),is=isa,isb)
      endif

      iss = iss + 7
    enddo
  enddo

  write(fid,*) 'ERROR_RESPONSE        imagZ' 
  do ip = 1,nPer
    iss = 0
    do ic = 1, nrow
      isa = iss + 1
      isb = iss + 7
  
      if (isb >= nSites) isb = nSites
      if (ic .eq. 1) then
        write(fid,200) periods(ip),(dat%d(ip)%err(2,is),is=isa,isb)
      else
        write(fid,300) ' ',(dat%d(ip)%err(2,is),is=isa,isb)
      endif

      iss = iss + 7
    enddo
  enddo

  close(fid)

100 format(7E12.4)
200 format(8E12.4)
300 format(a12,7E12.4)

  end subroutine writeData


! **********************************************************************************

end module wsio
