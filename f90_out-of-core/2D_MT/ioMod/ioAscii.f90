! *****************************************************************************
!  Module for basic input and output of standard data structures
!  for 2D MT modeling and inversion code
!
!  Here, we are using Randie Mackie's model format for input
!  and output of the model and the grid, and Naser Meqbel's
!  format for the data files.
!
!  The ascii input and output routines for cvector and solnVector
!  have not yet been created, so the binary versions are used
!  instead.

module ioAscii
   use math_constants
   use emfield
   use dataspace
   use ForwardSolver
   use utilities
   use transmitters
   use receivers
   use datatypes
   implicit none

   !  routines that are public
   private	::  read_cvector,write_cvector, &
		read_grid, write_grid

   public   :: write_Z, read_Z
   public   :: write_solnVectorMTX, read_solnVectorMTX


   Contains
     !******************************************************************
      subroutine read_cvector(fid,cfile,vec)

      !  open cfile on unit fid, read in object of
      !   type cvector in standard format
      !   vec must be allocated before calling


      integer, intent(in)		:: fid
      character*80, intent(in)	:: cfile
      type(cvector), intent(inout)		:: vec

      !  local variables
      integer 		:: N1, N2
      character*80	:: gridType


      if(vec%allocated) then
         open(unit=fid, file=cfile, form='unformatted')
         read(fid) gridType
         read(fid) N1,N2
         if((vec%N1 .NE. N1).OR.(vec%N2 .NE. N2)) then
            close(fid)
            call errStop('Size of vec does not agree with contents of file in read_cvector')
         else
            vec%gridType = gridType
            read(fid) vec%v
            close(fid)
         endif
      else
         call errStop('vec must be allocated before call to read_cvector')
      endif
      end subroutine read_cvector

     !******************************************************************
      subroutine openR_cvector(fid,cfile,header,nVec)

      !  open cfile on unit fid for reading multiple cvector objects
      !   reads header, number of vectors in file

      integer, intent(in)		:: fid
      character*80, intent(in)		:: cfile
      character*80, intent(out)		:: header
      integer, intent(out)		:: nVec

      open(unit=fid, file=cfile, form='unformatted',status = 'old')
      read(fid) header
      read(fid) nVec
      end subroutine openR_cvector

     !******************************************************************
      subroutine cvectorRead(fid,vec)

      !  reads one object of type cvector from unit fid
      !  by default, just reads next object ... might modify to
      !   allow skipping to read an arbitrary record number

      integer, intent(in)		:: fid
      type(cvector), intent(inout)	:: vec

      ! local variables:
      character*80	:: gridType, msg
      integer		:: N1,N2

      if(vec%allocated) then
         read(fid) gridType
         read(fid) N1,N2
         if((vec%N1 .NE. N1).OR.(vec%N2 .NE. N2)) then
            msg = 'Size of vec does not agree with contents of file'
            call errStop(msg)
         else
            vec%gridType = gridType
            read(fid) vec%v
         endif
      else
         msg = 'vec must be allocated before call to cvectorRead'
         call errStop(msg)
      endif

      end subroutine cvectorRead

     !******************************************************************
      subroutine openW_cvector(fid,cfile,header,nVec)

      !  open cfile on unit fid for writing multiple cvector objects
      !   writes header, number of vectors to be put in file

      integer, intent(in)		:: fid
      character*80, intent(in)		:: cfile
      character*80, intent(in)		:: header
      integer, intent(in)		:: nVec

      open(unit=fid, file=cfile, form='unformatted',status='unknown')
      write(fid) header
      write(fid) nVec
      end subroutine openW_cvector

     !******************************************************************
      subroutine cvectorWrite(fid,vec)

      !  writes one object of type cvector to unit fid, already opened

      integer, intent(in)		:: fid
      type(cvector), intent(in) 		:: vec

      write(fid) vec%gridType
      write(fid) vec%N1,vec%N2
      write(fid) vec%v

      end subroutine cvectorWrite

     !******************************************************************
      subroutine write_cvector(fid,cfile,vec)

      !  open cfile on unit fid, writes out object of
      !   type cvector in standard format (readable by matlab
      !   routine readcvector.m), closes file

      integer, intent(in)		:: fid
      character*80, intent(in)		:: cfile
      type(cvector), intent(in)		:: vec

      open(unit=fid, file=cfile, form='unformatted',status='unknown')
      write(fid) vec%gridType
      write(fid) vec%N1,vec%N2
      write(fid) vec%v
      close(fid)
      end subroutine write_cvector

     !******************************************************************
      subroutine write_solnVectorMTX(eAll,cfile)

      !  open cfile on unit ioE, writes out object of
      !   type cvector in standard format (readable by matlab
      !   routine readcvector.m), closes file
      !  NOT coded at present to specifically write out TE/TM
      !    solutions, periods, etc. (can get this infor from
      !    eAll%solns(j)%tx, but only with access to TXdict.

      character*80, intent(in)		:: cfile
      type(solnVectorMTX_t), intent(in)		:: eAll

      integer		:: j

      open(unit=ioE, file=cfile, form='unformatted',status='unknown')

      write(ioE) eAll%nTx
      do j = 1,eAll%nTx
         write(ioE) eAll%solns(j)%vec%gridType
         write(ioE) eAll%solns(j)%vec%N1,eAll%solns(j)%vec%N2
         write(ioE) eAll%solns(j)%vec%v
      enddo
      close(ioE)
      return
      end subroutine write_solnVectorMTX

     !******************************************************************
      subroutine read_solnVectorMTX(grid,eAll,cfile)

      !  open cfile on unit ioE, writes out object of
      !   type cvector in standard format (readable by matlab
      !   routine readcvector.m), closes file
      !  NOT coded at present to specifically write out TE/TM
      !    solutions, periods, etc. (can get this infor from
      !    eAll%solns(j)%tx, but only with access to TXdict.

      character*80, intent(in)      :: cfile
      type(grid_t), intent(in)      :: grid
      type(solnVectorMTX_t), intent(inout)     :: eAll

      integer       :: nTx,j

      open(unit=ioE, file=cfile, form='unformatted',status='unknown')
      if(eAll%allocated) then
        call deall_solnVectorMTX(eAll)
      endif
      read(ioE) nTx
      call create_solnVectorMTX(nTx,eAll)
      eAll%nTx = nTx
      do j = 1,eAll%nTx
         call create_solnVector(grid,j,eAll%solns(j))
         read(ioE) eAll%solns(j)%vec%gridType
         read(ioE) eAll%solns(j)%vec%N1,eAll%solns(j)%vec%N2
         read(ioE) eAll%solns(j)%vec%v
      enddo
      close(ioE)
      return
      end subroutine read_solnVectorMTX

      !******************************************************************
      subroutine read_rhsVectorMTX(grid,bAll,cfile,format)

          ! reads an array of solution vectors for all transmitters & subgrids
          ! currently uses the old binary format; will switch to NetCDF when
          ! time allows
          ! will allocate the grid from a full format file. for sparse files,
          ! this will only work if the grid is preallocated and compatible

          character(*), intent(in)                    :: cfile
          type(rhsVectorMTX_t), intent(inout)         :: bAll
          type(grid_t), intent(inout)                 :: grid
          character(*), optional, intent(in)          :: format

          !   local variables
          integer           :: j,k,jj,kk,nMode = 2, ios,istat,ig,cdot,filePer, nPol
          integer           :: iTx,nTx
          character (len=3)         :: igchar
          character (len=20)        :: version = '',source_type,tx_type
          character (len=2)         :: tmp,mode
          character (len=30)        :: str
          character (len=200)       :: fn_input
          logical                   :: sparse
          real (kind=prec)   :: omega, period

          sparse = .true.

          if (present(format)) then
            if (trim(format) .eq. 'full') then
                sparse = .false.
            end if
          end if

          if (.not. sparse) then

              fn_input = trim(cfile)

              write(*,*) 'Reading RHS E-fields from file: ',trim(fn_input)
              call setup_grid(grid)

              read(ioE) nTx
              call create_rhsVectorMTX(nTx,bAll)
              do iTx=1,nTx
                  bAll%combs(iTx)%nonzero_BC = .false.
                  bAll%combs(iTx)%nonzero_Source = .true.
                  call create_rhsVector(grid,iTx,bAll%combs(iTx))
              end do

              do j = 1,bAll%nTx
                  read(ioE) bAll%combs(j)%source%gridType
                  read(ioE) bAll%combs(j)%source%N1,bAll%combs(j)%source%N2
                  read(ioE) bAll%combs(j)%source%v
              enddo
              close(ioE)


          else

              open (unit=ioE,file=cfile,STATUS='unknown', form ='formatted', iostat=ios)

              if( ios/=0) then
                  write(0,*) 'Error opening sparse vector output file in read_rhsVectorMTX: ', cfile
              else
                  read(ioE,'(a30,i5)',iostat=istat) str,nTx
                  call create_rhsVectorMTX(nTx,bAll)
                  do iTx=1,nTx
                      bAll%combs(iTx)%nonzero_BC = .true.
                      bAll%combs(iTx)%nonzero_Source = .false.
                      call create_rhsVector(grid,iTx,bAll%combs(iTx))
                  end do


                  do j = 1,bAll%nTx

                          ! now that everything is allocated, set them both to false
                          ! because one of them could be missing from file...
                          bAll%combs(j)%nonzero_BC = .false.
                          bAll%combs(j)%nonzero_Source = .false.

                          read(ioE,'(a30,i5,es14.6,a2,a2)',iostat=istat) str, jj, period, tmp, mode
                          write(*,*) str,jj,period,mode
                          if (abs((period - txDict(bAll%combs(j)%tx)%period)/period) > TOL6) then
                              write(0,*) 'Warning: periods don''t match on RHS E-field input ',j
                              write(0,*) period, txDict(bAll%combs(j)%tx)%period
                          end if

                          read(ioE,*) bAll%combs(j)%bc
                          bAll%combs(j)%nonzero_BC = .true.
                          write(*,'(a34,i5,a6,a20)') 'Completed reading BC for period #',jj,' mode ',mode

                  enddo
                  close(ioE)
              end if

          end if

      end subroutine read_rhsVectorMTX

      !******************************************************************
      subroutine write_rhsVectorMTX(bAll,cfile,format)

          !  open cfile on unit fid, writes out object of
          !   type cvector in standard format (readable by matlab
          !   routine readcvector.m), closes file
          !  NOT coded at present to specifically write out TE/TM
          !    solutions, periods, etc. (can get this infor from
          !    eAll%solns(j)%tx, but only with access to TXdict.

          character(*), intent(in)                    :: cfile
          type(rhsVectorMTX_t), intent(in)            :: bAll
          character(*), intent(in), optional          :: format

          !   local variables
          type (cvector)    :: rhsFull
          type (sparsevecc) :: rhsSparse
          integer           :: j,k,nMode = 2, ios,istat,ig,cdot
          character (len=3)         :: igchar
          character (len=2)         :: mode
          character (len=20)        :: version = ''
          character (len=200)       :: fn_output
          logical            :: sparse
          real (kind=prec)   :: omega,period

          sparse = .true.

          if (present(format)) then
            if (trim(format) .eq. 'full') then
                sparse = .false.
            end if
          end if

          fn_output = trim(cfile)

          if (.not. sparse) then

              open(unit=ioE, file=cfile, form='unformatted',status='unknown')

              write(ioE) bAll%nTx
              do j = 1,bAll%nTx
                  if (bAll%combs(j)%nonzero_Source) then
                    write(ioE) bAll%combs(j)%source%gridType
                    write(ioE) bAll%combs(j)%source%N1,bAll%combs(j)%source%N2
                    write(ioE) bAll%combs(j)%source%v
                  end if
              enddo
              close(ioE)
              return

          else

              open (unit=ioE,file=cfile,STATUS='unknown', form ='formatted', iostat=ios)

              if( ios/=0) then
                  write(0,*) 'Error opening sparse vector output file in write_rhsVectorMTX: ', cfile
              else
                  write(ioE,'(a30,i5)',iostat=istat) 'Number of transmitters      : ',bAll%nTx
                  do j = 1,bAll%nTx
                          omega = txDict(bAll%combs(j)%tx)%omega
                          period = txDict(bAll%combs(j)%tx)%period
                          mode = txDict(bAll%combs(j)%tx)%mode
                          if (bAll%combs(j)%nonzero_BC) then
                            write(ioE,'(a30,i5,es14.6,a2,a2)',iostat=istat) 'iTx, period (secs) and mode : ',j,period,'  ',mode
                            write(ioE,*) bAll%combs(j)%bc
                          end if
                  enddo
                  close(ioE)
              end if

          end if

          write(*,*) 'RHS E-fields written to ',trim(fn_output)

      end subroutine write_rhsVectorMTX

     !**********************************************************************
      subroutine write_Z(fid,cfile,nTx,periods,modes,nSites,sites,units,allData)
     ! writes impedance file, including list of periods, siteLocations
     !   NOTE: this assumes that the arrays "sites" and "periods" are
     !    essentially identical to the receiver and transmitter dictionaries
     !    (in which case, why have both these arrays and the dicts?)
		 !   NOTE ALSO that for the 2D case the number of components = 2

      integer, intent(in)       :: fid
      character(*), intent(in)  :: cfile
      integer, intent(in)       :: nTx,nSites
      real(kind=8),intent(in)   :: periods(nTx)
      character*2, intent(in)	:: modes(nTx)
      real(kind=8),intent(in)   :: sites(2,nSites)
      type(dataVectorMTX_t), intent(in)      :: allData
      real(kind=8), dimension(:,:), allocatable :: siteTemp
      character(80)   :: info,units
      integer   :: nComp = 2
      complex*8  temp

      ! local variables
      type(dataBlock_t) :: data
      integer   :: nSite,iTx,iDt,k,j,isite,icomp
      character(10) :: siteid, tab='           '
      character(80) :: header
      real(kind=8) :: SI_factor

      info = 'Impedance responses from ModEM'
      if (index(units,'[V/m]/[T]')>0) then
         ! SI units for E/B
         SI_factor = 1.0
      else if (index(units,'[mV/km]/[nT]')>0) then
         ! practical units for E/B
         SI_factor = 1000.0
      else if ((index(units,'[V/m]/[A/m]')>0) .or. (index(units,'Ohm')>0)) then
         ! SI units for E/H
         SI_factor = 1000.0 * 10000.0/(4*PI) ! approx. 796000.0
      else
         call errStop('Unknown units in output data file '//cfile)
      end if

      open(unit=fid,file=cfile,form='formatted',status='unknown')
      write(fid,'(a12)',advance='no') 'Description:'
      write(fid,*) trim(info)
      write(fid,'(a6)',advance='no') 'Units:'
      write(fid,*) trim(units)
      write(fid,'(a17,i3)') 'Sign convention: ',ISIGN
      write(fid,*)

      write(fid,'(i5)') allData%nTx
      ! loop over periods
      do iTx = 1,allData%nTx
         do iDt = 1,allData%d(iTx)%nDt
             data = allData%d(iTx)%data(iDt)
	         header = tab//'  Re   '//tab//'   Im    '
	         nSite = data%nSite
	         nComp = data%nComp
	         ! write period, number of sites for this period
	         write(fid,'(es13.6,a5,i5)') periods(iTx),modes(iTx),nSite
	         ! write site locations for this period
	         allocate(siteTemp(2,nSite))
	         do k = 1,nSite
	            siteTemp(:,k) = sites(:,data%rx(k))
	         enddo
			! write latitude, longitude and elevation
			 do j = 1,2
			   do k = 1,nSite
	        	  write(fid,'(f14.3)',advance='no') siteTemp(j,k)
			   enddo
	           write(fid,*)
			 enddo
		    ! write the comment line for the data block
		     write(fid,'(a80)') header
	         do isite = 1,nSite
	         	! Note: temporarily, we write site id's according to their number;
	         	! in the future, they will be stored in the receiver dictionary
	         	write(siteid,'(i5)') isite
	         	write(fid,'(a10)',advance='no') trim(siteid)
	         	!  note that each data field in a dataVec (i.e., allData%d(iTx)%data(iDt))
	         	!   is a real array of size (nComp,nSite)
	         	do icomp = 1,nComp
	         		write(fid,'(es15.6)',advance='no') data%value(icomp,isite)/SI_factor
	         	enddo
	         	write(fid,*)
	         	write(fid,'(a10)',advance='no') tab
	         	do icomp = 1,nComp
	         	    if (data%errorBar) then
	         		   write(fid,'(es15.6)',advance='no') data%error(icomp,isite)/SI_factor
	         		else
	         		   write(fid,'(es15.6)',advance='no') R_ZERO
	         		endif
	         	enddo
	         	write(fid,*)
	         enddo
	         deallocate(siteTemp)
         enddo
      enddo
      close(fid)

      call deall_dataBlock(data)

      ! Calculating apparent resistivities and phases
      !temp=dcmplx(data%value(1,k),data%value(2,k))
      !write(fid,*)tab, ((periods(iTx)*MU_0)/(2.0*PI))*abs(temp)**2, abs(atan(data%value(2,k)/data%value(1,k)))*(180/PI)

      end subroutine write_Z


     !**********************************************************************
      subroutine read_Z(fid,cfile,nTx,periods,modes,nSites,sites,units,allData)
     ! reads in data file, returns list of periods, modes, siteLocations, and
     !   sets up data vector structure, including data and error bars
     ! First four lines are assumed to be:
     ! Description: (up to 80 char)
     ! Units: [V/m]/[A/m] OR [mV/km]/[nT]
     ! Sign Convention: -1 OR 1
     ! and an empty line
      integer, intent(in)       :: fid
      character(*), intent(in)  :: cfile
      integer, intent(out)      :: nTx,nSites
      real(kind=8),dimension(:), pointer     :: periods
      real(kind=8),dimension(:,:), pointer   :: sites
      character*2, dimension(:), pointer	:: modes
      type(dataVectorMTX_t), intent(inout)   :: allData

     ! local variables
      integer   :: nComp = 2
      integer   :: nSite,nDt,iTx,iDt,k,l,j,Ndata
      character(10)siteid
      real(kind=8), allocatable, dimension(:,:) :: siteTemp,siteTempAll
      logical   :: newSite, conjugate
      logical   :: isComplex, errorBar
      character(80) temp, description, units
      integer   :: sign_in_file
      real(kind=8) :: SI_factor

      open(unit=fid,file=cfile,status='old')
      read(fid,'(a13,a80)') temp,description
      read(fid,'(a7,a80)') temp,units
      read(fid,'(a17,i3)') temp,sign_in_file
      read(fid,*)

      if (index(units,'[V/m]/[T]')>0) then
         ! SI units for E/B
         SI_factor = 1.0
      else if (index(units,'[mV/km]/[nT]')>0) then
         ! practical units for E/B
         SI_factor = 1000.0
      else if ((index(units,'[V/m]/[A/m]')>0) .or. (index(units,'Ohm')>0)) then
         ! SI units for E/H
         SI_factor = 1000.0 * 10000.0/(4*PI) ! approx. 796000.0
      else
         call errStop('Unknown units in input data file '//cfile)
      end if

      if (sign_in_file == ISIGN) then
        conjugate = .false.
      else if (abs(sign_in_file) == 1) then
        conjugate = .true.
      else
        call errStop('Unknown sign convention in the data file '//cfile)
      end if

      read(fid,*) nTx
      call create(nTx, allData)
      allocate(periods(nTx))
      allocate(modes(nTx))


     ! loop over transmitters (periods/modes)
      Ndata = 0
      do iTx = 1,nTx

         ! hey, for now we have to assume that there's only one data type per period
         ! in the data file - otherwise the reading would be too tedious... (no need
         ! for this assumption anywhere else in the code, including write_Z routine)
         ! the structure of the data file should really be different.
         nDt = 1
         call create(nDt, allData%d(iTx))
         allData%d(iTx)%tx = iTx
         allData%d(iTx)%allocated = .true.

         do iDt = 1,nDt
	         ! read in number of sites for this period
	         read(fid,*) periods(iTx),modes(iTx),nSite
	         ! read in site locations
	         allocate(siteTemp(2,nSite))

	         read(fid,*) (siteTemp(1,k),k=1,nSite)
	         read(fid,*) (siteTemp(2,k),k=1,nSite)

	         ! read comment line just before the data block
	         read(fid,*)

	         ! create dataVec object, read in data
	         isComplex = .true.
	         errorBar = .true.
	         call create_dataBlock(nComp,nSite,allData%d(iTx)%data(iDt),isComplex,errorBar)
	         Ndata  = Ndata + nComp*nSite

	         allData%d(iTx)%data(iDt)%tx = iTx

		     if(modes(iTx) .eq. 'TM') then
		         allData%d(iTx)%data(iDt)%datatype = 2
		     else
		         allData%d(iTx)%data(iDt)%datatype = 1
		     endif

	         do k=1,nSite
	             read(fid,*)siteid, (allData%d(iTx)%data(iDt)%value(j,k),j=1,nComp)
	             read(fid,*)        (allData%d(iTx)%data(iDt)%error(j,k),j=1,nComp)
	         end do

	         ! convert data to SI units
	         allData%d(iTx)%data(iDt)%value = SI_factor * allData%d(iTx)%data(iDt)%value
	         allData%d(iTx)%data(iDt)%error = SI_factor * allData%d(iTx)%data(iDt)%error

	         ! conjugate data as necessary
	         if (conjugate) then
	           do j=2,nComp,2
	              allData%d(iTx)%data(iDt)%value(j,:) = - allData%d(iTx)%data(iDt)%value(j,:)
	           end do
	         end if

	         if(iTx .eq. 1) then
	           ! allocate temporary storage for full sites list
	           ! (this might not always work ... I am assuming that the
	           !        max number of sites is ns from period one * nTx)
	            allocate(siteTempAll(2,nSite*nTx))
	            nSites = nSite
	            do k = 1,nSite
	               siteTempAll(1,k) = siteTemp(1,k)
	               siteTempAll(2,k) = siteTemp(2,k)
	               allData%d(iTx)%data(iDt)%rx(k) = k
	            enddo
	         else
	            ! check to see if site locations are already in list
	            !  if not add to list; in any event set reciever "pointer" rx
	            do k = 1,nSite
	               newSite = .true.
	               do l = 1,nSites
	                  if((siteTemp(1,k).eq.siteTempAll(1,l)).and.  &
			  	(siteTemp(2,k).eq.siteTempAll(2,l))) then
	                     newSite = .false.
	                     allData%d(iTx)%data(iDt)%rx(k) = l
	                     exit
	                  endif
	               enddo
	               if(newSite) then
	                  nSites = nSites+1
	                  siteTempAll(1,nSites) = siteTemp(1,k)
	                  siteTempAll(2,nSites) = siteTemp(2,k)
	                  allData%d(iTx)%data(iDt)%rx(k) = nSites
	               endif
	            enddo
	         endif
	         deallocate(siteTemp)
         enddo
      enddo

      allData%allocated = .true.

      ! copy list of unique sites into "sites" array
      allocate(sites(2,nSites))
      do k = 1,nSites
         do j = 1,2
            sites(j,k) = siteTempAll(j,k)
         enddo
      enddo
      close(fid)
      end subroutine read_Z

     !**********************************************************************
      subroutine read_grid(fid,cfile,grid)
     !  reads in basic grid, allocating for Dy, Dz
      integer, intent(in)		:: fid
      character(*), intent(in)		:: cfile
      type(grid_t), intent(inout)	:: grid

      ! local variables:
      integer				:: Ny, Nz, Nza, NzEarth, j
      logical               :: newFile

      ! We are using Randie Mackie's format, which does not have information
      ! about the air layers. So we make it equal 10 in this routine.
      Nza = 10

      if (len_trim(cfile)>0) then
         newFile = .true.
         open(unit=fid,file=cfile,status='OLD')
      end if

      !  Read in grid geometry definitions, store in structure TEgrid
      !    first grid dimensions ...
      read(fid,*) Ny,NzEarth
      Nz = NzEarth + Nza
      ! then allocate for grid
      call create_grid(Ny,Nz,Nza,grid)

      !    read in grid spacings

        read(fid,*) (grid%Dy(j),j=1,Ny)

        read(fid,*) (grid%Dz(j),j=Nza+1,Nz)

        ! set the air layers spacing to that of the top 10 Earth layers
        if (NzEarth <= Nza) then
        	close(fid)
            call errStop('Too few Earth layers in the Mackie input model file in read_grid')
        else
        	do j = 1,Nza
        		grid%Dz(Nza-j+1) = grid%Dz(Nza+j)
        	end do
        end if

      if (newFile) then
         close(fid)
      end if

      ! complete grid definition
      call setup_grid(grid)

      end subroutine read_grid

     !**********************************************************************
      subroutine write_grid(fid,cfile,grid)
     !  writes basic grid, if cfile is empty, assumes file already open
      integer, intent(in)		    :: fid
      character(*), intent(in)		:: cfile
      type(grid_t), intent(in)	:: grid

      ! local variables:
      integer				:: Ny, Nz, Nza, NzEarth, j
      logical               :: newFile

      if (len_trim(cfile)>0) then
         newFile = .true.
         open(unit=fid,file=cfile,status='unknown')
      end if
      !  Read in grid geometry definitions, store in structure TEgrid
      !    first grid dimensions ...
      Ny=grid%ny
      Nz=grid%nz
      NzEarth=grid%nz - grid%nza
      Nza = grid%nza
      write(fid,'(2i5)') Ny,NzEarth

      !    write grid spacings: NOTE that Randie Mackie inserts an empty
      !    line after every 100 values in the grid (i.e. 10 lines);
      !    we do not do that here to save extra coding

        write(fid,'(10g11.4)') (grid%Dy(j),j=1,Ny)

        write(fid,'(10g11.4)') (grid%Dz(j),j=Nza+1,Nz)

      if (newFile) then
         close(fid)
      end if

      end subroutine write_grid

end module ioAscii
