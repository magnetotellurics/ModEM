! *****************************************************************************
module DataIO
  ! This module contains io routines for reading and writing the data vectors
  ! Generic interface required to call this from Level I inversion routines...
  ! Version: Global 3D

  use UserData
  use transmitters
  use receivers
  use dataTypes

  implicit none

  private

  complex(8), dimension(:), allocatable, save   :: avgHth

  public     :: write_dataVectorMTX, initData, outputResponses


Contains

  ! ***************************************************************************
  ! * initData initializes the data according to the information stored in
  ! * the list of transfer functions TFList; then merges all data types into
  ! * a single data vector allData
  subroutine initData(cUserDef,allData,obsList,freqList,TFList)

    implicit none
    type (userdef_control), intent(in)                          :: cUserDef
    type (dataVectorMTX_t), intent(inout)                   :: allData
    type (Obs_List), target, intent(inout)                  :: obsList
    type (Freq_List), target, intent(inout)                 :: freqList
    type (TF_List), target, intent(in)                      :: TFList
    ! local
    type (dataVectorMTX_t)     :: cdata,ddata,hdata
    real(8)                    :: theta
    complex(8)                 :: H
    integer                    :: i,j,k,iTx,iDt,iRx,numHth,istat

    ! Read the data from files into data blocks, and save in the data vector
    do iDt=1,TFList%n

      select case ( TFList%info(iDt)%name )

      case ('C')
        ! Read C responses into a data vector with a single data type (cdata)
        call initDataList(cdata,iDt,obsList,freqList,cUserDef%fn_cdata)
        if (.not. cdata%allocated) then
          write(0,*) node_info,'Warning: Data misfit will not be calculated for C responses: data file not available'
          call initDataList(cdata,iDt,obsList,freqList)
        end if
        do i=1,cdata%ntx
          do k=1,cdata%d(i)%data(1)%nSite
                iRx = cdata%d(i)%data(1)%rx(k)
                theta = obsList%info(iRx)%colat*d2r;
                if ((theta >= pi/2-EPS_GRID).and.(theta <= pi/2+EPS_GRID)) then
                  write(0,*) node_info,'Error: (initData) receiver ',trim(obsList%info(iRx)%code),' located too close to equator;'
                  write(0,*) node_info,'Error: (initData) please delete this receiver from the list and try again.'
                    exit
                end if
                ! convert from C responses to C ratios
                cdata%d(i)%data(1)%value(:,k) = cdata%d(i)%data(1)%value(:,k)/dtan(theta)
                cdata%d(i)%data(1)%error(:,k) = cdata%d(i)%data(1)%error(:,k)/abs(dtan(theta))
          end do
        end do

      case ('D')
        ! Read D responses into a data vector with a single data type (ddata)
        call initDataList(ddata,iDt,obsList,freqList,cUserDef%fn_ddata)
        if (.not. ddata%allocated) then
          write(0,*) node_info,'Warning: Data misfit will not be calculated for D responses: data file not available'
          call initDataList(ddata,iDt,obsList,freqList)
        end if
        do i=1,ddata%ntx
          do k=1,ddata%d(i)%data(1)%nSite
                iRx = ddata%d(i)%data(1)%rx(k)
                theta = obsList%info(iRx)%colat*d2r;
                if ((theta >= pi/2-EPS_GRID).and.(theta <= pi/2+EPS_GRID)) then
                  write(0,*) node_info,'Error: (initData) receiver ',trim(obsList%info(iRx)%code),' located too close to equator;'
                  write(0,*) node_info,'Error: (initData) please delete this receiver from the list and try again.'
                    exit
                end if
                ! convert from D responses to D ratios
                ddata%d(i)%data(1)%value(:,k) = ddata%d(i)%data(1)%value(:,k)/dsin(theta)
                ddata%d(i)%data(1)%error(:,k) = ddata%d(i)%data(1)%error(:,k)/abs(dsin(theta))
          end do
        end do

      case ('T') ! generalised magnetic field transfer functions (i.e., scaled magnetic fields)
        ! Read raw H fields into a data vector with a single data type (hdata)
        call initDataList(hdata,iDt,obsList,freqList,cUserDef%fn_hdata)
        if (.not. hdata%allocated) then
          write(0,*) node_info,'Warning: Data misfit will not be calculated for H responses: data file not available'
          call initDataList(hdata,iDt,obsList,freqList)
        end if
        ! once frequency list is defined, allocate space for avgHth and compute them
        allocate(avgHth(freqList%n),STAT=istat)
        do i=1,hdata%ntx
          iTx = freqList%info(i)%i
          avgHth(iTx) = Dst_index(hdata%d(i)%data(1))
          ! convert from raw H fields to H responses
          do k=1,hdata%d(i)%data(1)%nSite
            do j=1,3
                H    = cmplx(hdata%d(i)%data(1)%value(2*j-1,k),hdata%d(i)%data(1)%value(2*j,k)) / avgHth(iTx)
                hdata%d(i)%data(1)%value(2*j-1,k) = dreal(H)
                hdata%d(i)%data(1)%value(2*j,k)   = dimag(H)
                hdata%d(i)%data(1)%error(2*j-1,k) = hdata%d(i)%data(1)%error(2*j-1,k) / abs(avgHth(iTx))
                hdata%d(i)%data(1)%error(2*j,k)   = hdata%d(i)%data(1)%error(2*j,k) / abs(avgHth(iTx))
             end do
           end do
        end do

      case ('H')
        ! Read raw H fields into a data vector with a single data type (hdata)
        call initDataList(hdata,iDt,obsList,freqList,cUserDef%fn_hdata)
        if (.not. hdata%allocated) then
          write(0,*) node_info,'Warning: Data misfit will not be calculated for H responses: data file not available'
          call initDataList(hdata,iDt,obsList,freqList)
        end if

      case default

        write(0,*) node_info,'Please specify correct transfer functions in the file ',trim(cUserDef%fn_func)
        stop

      end select

    end do

    ! Now, merge the data types into a single data vector
    if (cdata%allocated .or. ddata%allocated) then
        call merge_dataVectorMTX(cdata,ddata,allData)
    end if

    ! Add H fields, if present
    if (hdata%allocated) then
        call merge_dataVectorMTX(hdata,allData,allData)
    end if

    ! Clean up
    call deall_dataVectorMTX(cdata)
    call deall_dataVectorMTX(ddata)
    call deall_dataVectorMTX(hdata)

  end subroutine initData ! initData

  ! ***************************************************************************
  ! * initData reads the file fn_data (in future could be several files) that
  ! * contain all the information about the available data: the values of the
  ! * responses, data errors corresponding to an observatory and a frequency.

  subroutine initDataList(mydat,itype,obsList,freqList,fname)

    implicit none
    type (dataVectorMTX_t), intent(inout)      :: mydat
    integer, intent(in)                        :: itype
    type (Obs_List), target, intent(in)        :: obsList
    type (Freq_List), target, intent(in)       :: freqList
    character(*), intent(in), optional         :: fname
    ! local
    integer                 :: ifreq,iobs
    character(100)          :: label
    integer                 :: i,j,k,ios=0,istat=0
    real(8)                 :: const,large
    real(8)                 :: freq,days,prev
    character(80)           :: code
    real(8)                 :: lon,lat
    real(8)                 :: creal(3),cimag(3),cerr(3)
    integer                 :: nfreq,nobs,ncomp,nsite
    integer                 :: countData,countFreq
    complex(8), allocatable :: value(:,:,:) ! (nfreq,nobs,ncomp)
    real(8), allocatable    :: error(:,:,:) ! (nfreq,nobs,ncomp)
    logical, allocatable    :: exist(:,:,:) ! (nfreq,nobs,ncomp)
    logical                 :: new,isComplex,errorBar

    large = 2.0e15
    const = 1.0d-2
    prev = 0.0d0

    ! Initialize the data functionals
    nfreq = freqList%n
    nobs  =  obsList%n
    ncomp =   TFList%info(itype)%nComp
    allocate(value(nfreq,nobs,ncomp),STAT=istat)
    allocate(error(nfreq,nobs,ncomp),STAT=istat)
    allocate(exist(nfreq,nobs,ncomp),STAT=istat)
    value(:,:,:) = dcmplx(0.0d0,0.0d0)
    error(:,:,:) = large
    exist(:,:,:) = .FALSE.

    ! If called without a data file, use dictionaries to create data vector
    if (.not. present(fname)) then
        call create_dataVectorMTX(nfreq,mydat)
        mydat%allocated = .TRUE.
        do i = 1,nfreq
           isComplex = .TRUE.
           errorBar = .TRUE.
           call create_dataVector(1,mydat%d(i))
           mydat%d(i)%tx = i
           mydat%d(i)%allocated = .TRUE.
           call create_dataBlock(2*ncomp,nobs,mydat%d(i)%data(1),isComplex,errorBar)
           do j = 1,nobs
               do k = 1,ncomp
                   exist(i,j,k) = .TRUE.
                   if(exist(i,j,k)) then
                       mydat%d(i)%data(1)%value(2*k-1,j) = real(value(i,j,k))
                       mydat%d(i)%data(1)%value(2*k,j) = imag(value(i,j,k))
                       mydat%d(i)%data(1)%error(2*k-1,j) = error(i,j,k)
                       mydat%d(i)%data(1)%error(2*k,j) = error(i,j,k)
                       mydat%d(i)%data(1)%exist(2*k-1,j) = exist(i,j,k)
                       mydat%d(i)%data(1)%exist(2*k,j) = exist(i,j,k)
                       mydat%d(i)%data(1)%rx(j) = j
                   end if
               end do
           end do
           mydat%d(i)%data(1)%dataType = itype
           mydat%d(i)%data(1)%tx = i
           mydat%d(i)%data(1)%allocated = .TRUE.
        end do
        deallocate(value,STAT=istat)
        deallocate(error,STAT=istat)
        deallocate(exist,STAT=istat)
        return
    end if

    ! Otherwise, if data file does not exist, exit
    inquire(FILE=trim(fname),EXIST=exists)
    if (.not.exists) then
      deallocate(value,STAT=istat)
      deallocate(error,STAT=istat)
      deallocate(exist,STAT=istat)
      return
    end if

    ! If data file is present, initialize data, functionals and residuals
    open(ioDat,file=trim(fname),status='old',form='formatted',iostat=ios)

    write(6,*) node_info,'Reading from the data file ',trim(fname)
    read(ioDat,'(a)') label
    write(6,*) label
    read(ioDat,'(a)') label

    countData = 0
    countFreq = 0


    READ_DATA: do

      !read(ioDat,'(a12,1x)',advance='no',iostat=ios)  label !'# freq, obs: '
      !read(ioDat,*,iostat=ios) freq, num

      if (ios /= 0) exit

      ! Assume either 1 component (C/D responses) or 3 components (H-fields)
      if (ncomp == 1) then
        read(ioDat,*,iostat=ios) days,code,lon,lat,creal(1),cimag(1),cerr(1)
        creal(1) = km2m * creal(1)
        cimag(1) = km2m * cimag(1)
        cerr(1)  = km2m * cerr(1)
      else
        read(ioDat,*,iostat=ios) days,code,lon,lat,creal(1),cimag(1),cerr(1), &
                      creal(2),cimag(2),cerr(2),creal(3),cimag(3),cerr(3)
      end if

      ! Find the relevant frequency in the list
      freq  = 1/(days*24*60*60)
      new = .FALSE.
      if (abs((prev-freq)/freq) > const) then
        new = .TRUE.
        prev = freq
      end if
      if (new) then
        ifreq = getFreq(freqList,freq)
        if (ifreq > 0) then
           write(6,'(a12,a32,i6,a2,es12.6,a5)') node_info,'Reading data for the period ',ifreq,': ',days,' days'!freqList%info(i)%value
           countFreq = countFreq + 1
        end if
      end if
      if (ifreq == 0) then
        print *, 'Frequency ',freq,' is not in the frequency list; ignore'
        !do j=1,num
          !read(ioDat,*,iostat=ios)
        !end do
        cycle READ_DATA
      end if
      ! Find the relevant observatory in the list
      iobs = getObs(obsList,code)
      if (iobs > 0) then
        countData = countData + ncomp
      else
        print *, 'Observatory ',trim(code),' in file ',trim(fname),' is not in our list'
        !print *, 'This is an error; exiting...'
        !stop
        cycle READ_DATA
      end if

      ! If everything is correct, save the data into the 2D arrays
      do k = 1,ncomp
          value(ifreq,iobs,k) = dcmplx(creal(k),cimag(k))
          error(ifreq,iobs,k) = cerr(k)
          exist(ifreq,iobs,k) = .TRUE.
      end do

      ! Now account for "non-existent" observatories, e.g., too close to poles
      if (.not. obsList%info(iobs)%defined) then
          exist(ifreq,iobs,:) = .FALSE.
      end if

    end do READ_DATA

    close(ioDat)

    if (countData==0) then
      write(0,*) node_info,'Warning: No data matches given frequencies and observatories'
    else
      write(0,*) node_info,'Number of complex data values in total: ',countData
    end if

    ! Finally, store the data in the data vector
    call create_dataVectorMTX(countFreq,mydat)
    mydat%allocated = .TRUE.
    i = 0
    ifreq = 0
    SAVE_DATA: do ifreq = 1,nfreq

       nSite = count(exist(ifreq,:,1)) ! may need to be generalized
       isComplex = .TRUE.
       errorBar = .TRUE.
       if(nSite == 0) then
           ! no data for this frequency ... try next one
           cycle SAVE_DATA
       end if
       i = i + 1
       call create_dataVector(1,mydat%d(i))
       mydat%d(i)%tx = ifreq
       mydat%d(i)%allocated = .TRUE.
       call create_dataBlock(2*ncomp,nSite,mydat%d(i)%data(1),isComplex,errorBar)
       j = 1
       do iobs = 1,nobs
           if(count(exist(ifreq,iobs,:))>0) then
               do k = 1,ncomp
                   mydat%d(i)%data(1)%value(2*k-1,j) = real(value(ifreq,iobs,k))
                   mydat%d(i)%data(1)%value(2*k,j) = imag(value(ifreq,iobs,k))
                   mydat%d(i)%data(1)%error(2*k-1,j) = error(ifreq,iobs,k)
                   mydat%d(i)%data(1)%error(2*k,j) = error(ifreq,iobs,k)
                   mydat%d(i)%data(1)%exist(2*k-1,j) = exist(ifreq,iobs,k)
                   mydat%d(i)%data(1)%exist(2*k,j) = exist(ifreq,iobs,k)
               end do
               mydat%d(i)%data(1)%rx(j) = iobs
               j = j + 1
          end if
       end do
       mydat%d(i)%data(1)%dataType = itype
       mydat%d(i)%data(1)%tx = ifreq
       mydat%d(i)%data(1)%allocated = .TRUE.

    end do SAVE_DATA

    deallocate(value,STAT=istat)
    deallocate(error,STAT=istat)
    deallocate(exist,STAT=istat)

    return

  end subroutine initDataList ! initData

!**********************************************************************
! writes global responses file in ASCII format

   subroutine write_dataVectorMTX(allResp,cfile)

      character(*), intent(in)                      :: cfile
      type(dataVectorMTX_t), intent(in)             :: allResp
      ! local variables
      integer									:: i
      character(200)								:: strtmp

      do i=1,allData%nTx

		write(strtmp,*) trim(cfile),'.cout'
		outFiles%fn_cdat = strtmp

		write(strtmp,*) trim(cfile),'.dout'
		outFiles%fn_ddat = strtmp

        write(strtmp,*) trim(cfile),'.hout'
        outFiles%fn_hdat = strtmp

		call outputResponses(allResp%d(i),outFiles,allData%d(i))

      end do

   end subroutine write_dataVectorMTX

  ! ***************************************************************************
  ! * OutputResponses writes the chosen kind of responses calculated at every
  ! * observatory location to an output file

  subroutine outputResponses(psi,outFiles,dat)

    type (dataVector_t), intent(in)                 :: psi
    type (dataVector_t), intent(in)                 :: dat
    type (output_info), intent(in)                  :: outFiles
    character(200)                                  :: fn_response
    character(200)                                  :: comment
    complex(8), dimension(:,:), allocatable         :: Resp,RespRatio
    complex(8), dimension(:,:), allocatable         :: FieldData
    real(8), dimension(:,:), allocatable            :: FieldError
    real(8)                                         :: rval,ival,err,rms
    integer                                         :: i,j,k,icomp,itype,iobs,nSite,nComp
    integer                                         :: ios,istat

    i = psi%tx

    do j=1,psi%ndt

      nSite = psi%data(j)%nSite
      nComp = psi%data(j)%nComp/2
      allocate(Resp(nSite,nComp),RespRatio(nSite,nComp),FieldData(nSite,nComp),FieldError(nSite,nComp), STAT=istat)
      do icomp=1,nComp
        RespRatio(:,icomp) = cmplx(psi%data(j)%value(2*icomp-1,:),psi%data(j)%value(2*icomp,:))
        FieldData(:,icomp) = cmplx(dat%data(j)%value(2*icomp-1,:),dat%data(j)%value(2*icomp,:))
        FieldError(:,icomp) = dat%data(j)%error(2*icomp-1,:)
      end do
      itype = dat%data(j)%dataType

      select case ( trim(TFList%info(itype)%name) )

      case ('C')
        fn_response = outFiles%fn_cdat
        do k=1,nSite
            iobs = psi%data(j)%rx(k)
            Resp(k,:) = RespRatio(k,:) * dtan(obsList%info(iobs)%colat*d2r) * m2km
        end do
        comment = "#Period Code GM_Lon GM_Lat Real(km) Imag(km) Error(km)"

      case ('D')
        fn_response = outFiles%fn_ddat
        do k=1,nSite
            iobs = psi%data(j)%rx(k)
            Resp(k,:) = RespRatio(k,:) * dsin(obsList%info(iobs)%colat*d2r) * m2km
        end do
        comment = "#Period Code GM_Lon GM_Lat Real(km) Imag(km) Error(km)"

      case ('T') ! generalised magnetic field transfer functions computed; scale them to obtain magnetic fields
        fn_response = outFiles%fn_hdat
        Resp(:,:) = RespRatio(:,:) * avgHth(i)
        comment = "#Period Code GM_Lon GM_Lat Real(Hp) Imag(Hp) Error(Hp) Real(Ht) Imag(Ht) Error(Ht) Real(Hr) Imag(Hr) Error(Hr)"

      case ('H') ! raw magnetic fields computed; no scaling required
        fn_response = outFiles%fn_hdat
        Resp(:,:) = RespRatio(:,:)
        comment = "#Period Code GM_Lon GM_Lat Real(Hp) Imag(Hp) Error(Hp) Real(Ht) Imag(Ht) Error(Ht) Real(Hr) Imag(Hr) Error(Hr)"

      case default
        write(0,*) 'Warning: unknown transfer function: ',&
          trim(TFList%info(itype)%name)
        cycle
      end select

!          inquire(FILE=fn_response,EXIST=exists)
!     if ((.not.exists).and.(i==1)) then
!       open(ioResp,file=fn_response,status='unknown',form='formatted',iostat=ios)
!       write(ioResp,*) "#Output of earth3d (for ",freqList%n," frequency values)";
!       write(ioResp,*) "#Period Code GM_Lon GM_Lat Real(km) Imag(km) Error(km)";
!     else if(exists .and.(i==1)) then
!                write(0,*) 'Warning: Response file already exists. Appending...'
!       open(ioResp,file=fn_response,position='append', form='formatted',iostat=ios)
!          else
!       open(ioResp,file=fn_response,position='append', form='formatted',iostat=ios)
!     end if


        ! NB: No appending for now to make it easier to debug and cleaner. Jan 19, 2007
      if (i==1) then
          open(ioResp,file=fn_response,status='unknown',form='formatted',iostat=ios)
          write(ioResp,*) "#Output of earth3d (for ",freqList%n," frequency values)"
          write(ioResp,*) trim(comment)
      else
          open(ioResp,file=fn_response,position='append', form='formatted',iostat=ios)
      end if
        !write(ioResp,*) "freq = ",freq%value
        do k=1,nSite
          iobs = psi%data(j)%rx(k)
          if (.not.obsList%info(iobs)%defined) then
            cycle
          end if
          write(ioResp,'(f8.3,a12,2g15.7)',advance='no') &
              freqList%info(i)%period,&
              trim(obsList%info(iobs)%code),&
              obsList%info(iobs)%lon,obsList%info(iobs)%lat

          do icomp = 1,nComp
            write(ioResp,'(2g15.7)',advance='no') Resp(k,icomp)
            if (count(dat%data(j)%exist(:,k))==2*nComp) then
              rval = dreal(FieldData(k,icomp)-RespRatio(k,icomp))
              ival = dimag(FieldData(k,icomp)-RespRatio(k,icomp))
              err = FieldError(k,icomp)
              rms = ((rval/err)**2 + (ival/err)**2)/2
            else
              rms = 999999.9
            end if
            write(ioResp,'(g15.7)',advance='no') rms
          end do
          write(ioResp,*)
        end do
      close(ioResp)

      deallocate(Resp,RespRatio,FieldData,FieldError, STAT=istat)

!     if (fn_response == '') then
!       do k=1,size(Resp)
!         write(*,*) trim(obsList%info(k)%code),Resp(k) * m2km
!       end do
!     end if

    end do

    return

  end subroutine outputResponses  ! outputResponses

  ! ***************************************************************************
  ! * Computes a quantity similar to the Dst index: the average value of Hth
  ! * for all reference observatories in the refObsList; assume that the input
  ! * is a data vector for raw magnetic fields
  function Dst_index(data) result (avgHth)

    type (dataBlock_t), intent(in)      :: data
    complex(8)                          :: avgHth
    integer                             :: numHth,iRx,k

    if (data%allocated) then
      avgHth = C_ZERO
      numHth = 0
      do k=1,data%nSite
            iRx = data%rx(k)
            ! if one of the reference observatories, add to the "Dst index"
            if (getObs(refObsList,obsList%info(iRx)%code)>0) then
                avgHth = avgHth + cmplx(data%value(3,k),data%value(4,k))
                numHth = numHth + 1
            end if
      end do
      if (numHth > 0) then
        avgHth = avgHth / numHth
      else
        avgHth = C_ONE
      end if
    else
      avgHth = C_ONE
    end if

  end function Dst_index


end module DataIO
