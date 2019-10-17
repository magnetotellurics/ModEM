! *****************************************************************************
module DataIO
  ! This module contains io routines for reading and writing the data vectors
  ! Version: 3D MT

  use math_constants
  use file_units
  use utilities
  use dataspace
  use transmitters
  use receivers
  use datatypes

  implicit none

  private

  ! switch between data formats by leaving uncommented one of the options below
  interface read_dataVectorMTX
	MODULE PROCEDURE read_Z_list
  end interface

  interface write_dataVectorMTX
	MODULE PROCEDURE write_Z_list
  end interface


  public     :: read_dataVectorMTX, write_dataVectorMTX

  type :: data_file_block

      ! this block of information constitutes user preferences about the data format;
      ! there is one entry per each transmitter type and data type... (iTxt,iDt)
      ! if there are multiple data blocks of the same transmitter & data types,
      ! the last value is used.
      character(200) :: info_in_file
      character(20)  :: sign_info_in_file
      integer        :: sign_in_file
      character(20)  :: units_in_file
      real           :: origin_in_file(2)
      real           :: geographic_orientation

     ! these lists contain the indices into the data vector for each data type;
     ! they make it possible to sort the data by receiver for output.
     ! no data denoted by zero index; dimensions (nTx) and (nTx,nRx).
     ! these indices are typically allocated as we read the data file
     integer, pointer, dimension(:)   :: tx_index
     integer, pointer, dimension(:)   :: dt_index
     integer, pointer, dimension(:,:) :: rx_index

     ! some transmitter types and data types don't go together
     logical         :: defined

  end type data_file_block

  ! private dictionary of data block info dimension (nTxt,nDt)
  ! where nTxt = number of all possible transmitter types
  !       nDt  = number of all possible data types
  ! number of transmitter types comes from the DICT/txTypes module
  ! and defines the number of conceptually different types of sources
  type (data_file_block), pointer, save, private, dimension(:,:) :: fileInfo


Contains

!**********************************************************************
! Sorts out the data block header

  function DataBlockHeader(txType,dataType) result (header)

    integer, intent(in)         :: txType
    integer, intent(in)         :: dataType
    character(200)              :: header

    select case (dataType)

       case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)
          header = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error'

       case(Full_Interstation_TF)
          header = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Ref_Code Ref_Lat '// &
                   'Ref_Lon Ref_X(m) Ref_Y(m) Ref_Z(m) Component Real Imag Error'

       case(Off_Diagonal_Rho_Phase,Phase_Tensor)
          header = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Value Error'

    end select

  end function DataBlockHeader

  ! **************************************************************************
  ! Cleans up and deletes type dictionary at end of program execution
  subroutine init_fileInfo(nTxt,nDt,nTx,nRx)

    integer, intent(in) :: nTxt,nDt
    integer, intent(in), optional :: nTx,nRx
    integer     :: istat,iTxt,iDt

    allocate(fileInfo(nTxt,nDt),STAT=istat)


     do iTxt = 1,nTxt
       do iDt = 1,nDt
         fileInfo(iTxt,iDt)%defined = .false.
         if (present(nTx) .and. present(nRx)) then
           allocate(fileInfo(iTxt,iDt)%tx_index(nTx),STAT=istat)
           allocate(fileInfo(iTxt,iDt)%dt_index(nTx),STAT=istat)
           allocate(fileInfo(iTxt,iDt)%rx_index(nTx,nRx),STAT=istat)
         end if
       end do
     end do

  end subroutine init_fileInfo

  ! **************************************************************************
  ! Cleans up and deletes type dictionary at end of program execution
  subroutine deall_fileInfo()

    integer     :: i,j, istat

    if (associated(fileInfo)) then

     do i = 1,size(fileInfo,1)
       do j = 1,size(fileInfo,2)
          if (associated(fileInfo(i,j)%tx_index)) then
             deallocate(fileInfo(i,j)%tx_index,STAT=istat)
          end if
          if (associated(fileInfo(i,j)%dt_index)) then
             deallocate(fileInfo(i,j)%dt_index,STAT=istat)
          end if
          if (associated(fileInfo(i,j)%rx_index)) then
             deallocate(fileInfo(i,j)%rx_index,STAT=istat)
          end if
       end do
     end do

     deallocate(fileInfo,STAT=istat)

    end if

  end subroutine deall_fileInfo

!**********************************************************************
! writes data in the ASCII list data file; it is convenient to work
! with data ordered by site (NOT by frequency), so we are going into
! some pains here to write them out in this order ...

   subroutine write_Z_list(allData,cfile)

    character(*), intent(in)                  :: cfile
    type(dataVectorMTX_t), intent(in)         :: allData
    ! local variables
    integer                         :: nTx,nRx,nDt,ncomp
    integer                         :: countData
    real(8), allocatable            :: value(:) ! (ncomp)
    real(8), allocatable            :: error(:) ! (ncomp)
    logical, allocatable            :: exist(:) ! (ncomp)
    character(2)                    :: temp = '> '
    character(40)                   :: siteid,ref_siteid,compid
    character(1000)                 :: strtemp
    integer                         :: iTxt,iTx,iRx,iDt,icomp,i,j,k,istat,ios,nBlocks
    real(8)                         :: x(3),ref_x(3), Period,SI_factor,large
    real(8)                         :: lat,lon,ref_lat,ref_lon
    logical                         :: conjugate, isComplex

    ! 2019.04.23, Liu Zhongyin, add local varible azimu
    real(8), allocatable            :: azimu(:) ! (ncomp)

    iTxt = 1

    open(unit=ioDat,file=cfile,form='formatted',status='unknown')

    ! For each data type in dictionary, if data of this type exists, write it out.
    WRITE_DATA_TYPE: do iDt = 1,size(typeDict)

      nBlocks = countDataBlock(allData,iDt,iTxt)
      if (nBlocks == 0) then
	! no data for this data type; skip it - this shouldn't happen anymore
	! since the "defined" logical deals with this on input
        cycle WRITE_DATA_TYPE
      else
        ! count the number of transmitters and receivers
        nTx = 0
        nRx = 0
        do iTx = 1,size(txDict)
            if (fileInfo(iTxt,iDt)%tx_index(iTx) > 0) then
                nTx = nTx + 1
            end if
        end do
        do iRx = 1,size(rxDict)
            do iTx = 1,size(txDict)
                if (fileInfo(iTxt,iDt)%rx_index(iTx,iRx) > 0) then
                    nRx = nRx + 1
                    exit
                end if
            end do
        end do
      end if

      ! write the data type header
      call compact(fileInfo(iTxt,iDt)%info_in_file)
      write(strtemp,*) adjustl(trim(fileInfo(iTxt,iDt)%info_in_file))
      write(ioDat,'(a32)',advance='no') '# ModEM impedance responses for '
      write(ioDat,*,iostat=ios) strtemp(1:100)
      write(strtemp,*) adjustl(trim(DataBlockHeader(iTxt,iDt)))
      write(ioDat,'(a2)',advance='no') '# '
      write(ioDat,*,iostat=ios) strtemp(1:100)
      !if (.not. (tx_type_name(iTxt) .eq. 'MT')) then
      !    write(ioDat,'(a2)',advance='no') '+ '
      !    write(ioDat,*,iostat=ios) trim(tx_type_name(iTxt))
      !end if
      call compact(typeDict(iDt)%name)
      write(ioDat,'(a2)',advance='no') temp
      write(ioDat,*,iostat=ios) trim(typeDict(iDt)%name)
      call compact(fileInfo(iTxt,iDt)%sign_info_in_file)
      write(ioDat,'(a2)',advance='no') temp
      write(ioDat,*,iostat=ios) trim(fileInfo(iTxt,iDt)%sign_info_in_file)
      call compact(fileInfo(iTxt,iDt)%units_in_file)
      write(ioDat,'(a2)',advance='no') temp
      write(ioDat,*,iostat=ios) trim(fileInfo(iTxt,iDt)%units_in_file)
      write(ioDat,'(a2,f8.2)',iostat=ios) temp,fileInfo(iTxt,iDt)%geographic_orientation
      write(ioDat,'(a2,2f9.3)',iostat=ios) temp,fileInfo(iTxt,iDt)%origin_in_file(1),fileInfo(iTxt,iDt)%origin_in_file(2)
      write(ioDat,'(a2,2i6)',iostat=ios) temp,nTx,nRx

      if (fileInfo(iTxt,iDt)%sign_in_file == ISIGN) then
          conjugate = .false.
      else if (abs(fileInfo(iTxt,iDt)%sign_in_file) == 1) then
          conjugate = .true.
      end if
      SI_factor = ImpUnits(typeDict(iDt)%units,fileInfo(iTxt,iDt)%units_in_file)

      ncomp = typeDict(iDt)%nComp
      allocate(value(ncomp),error(ncomp),exist(ncomp),STAT=istat)

      ! 2019.04.23, Liu Zhongyin, add azimu allocation
      allocate(azimu(ncomp),stat=istat)

      isComplex = typeDict(iDt)%isComplex
      countData = 0

      ! write data
      do iRx = 1,size(rxDict)
        do iTx = 1,size(txDict)

            k = fileInfo(iTxt,iDt)%rx_index(iTx,iRx)
            i = fileInfo(iTxt,iDt)%dt_index(iTx)
            j = fileInfo(iTxt,iDt)%tx_index(iTx)
            if (k == 0) then
                cycle
            end if
            value = SI_factor * allData%d(j)%data(i)%value(:,k)
            if (allData%d(j)%data(i)%errorBar) then
                error = SI_factor * allData%d(j)%data(i)%error(:,k)
            else
                error = LARGE_REAL
            end if
            exist = allData%d(j)%data(i)%exist(:,k)
            Period = txDict(iTx)%period
            siteid = rxDict(iRx)%id
            x = rxDict(iRx)%x

            ! 2019.04.23, Liu Zhongyin, assign azimuth to azimu
            azimu = allData%d(j)%data(i)%Azimuth(:,k)

            select case (iDt)

                case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)

                    do icomp = 1,ncomp/2
                        if (.not. exist(2*icomp-1)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        write(ioDat,'(es12.6)',    iostat=ios,advance='no') Period
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(a40,3f15.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        !if (conjugate) then
                        !    write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),-value(2*icomp),error(2*icomp)
                        !else
                        !    write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp),error(2*icomp)
                        !end if
                        ! 2019.04.23, Liu Zhongyin, add azimu while writing
                        if (conjugate) then
                            write(ioDat,'(a8,4es15.6)',iostat=ios) trim(compid),value(2*icomp-1),-value(2*icomp),error(2*icomp),azimu(2*icomp)
                        else
                            write(ioDat,'(a8,4es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp),error(2*icomp),azimu(2*icomp)
                        end if
                        countData = countData + 1
                    end do

                case(Full_Interstation_TF)

                    do icomp = 1,ncomp/2
                        if (.not. exist(2*icomp-1)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        ref_siteid = rxDict(iRx)%id_ref
                        ref_x = rxDict(iRx)%r
                        write(ioDat,'(es12.6)',    iostat=ios,advance='no') Period
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(a40,3f15.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        write(ioDat,'(a40,3f15.3)',iostat=ios,advance='no') trim(ref_siteid),ref_x(:)
                        !if (conjugate) then
                        !    write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),-value(2*icomp),error(2*icomp)
                        !else
                        !    write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp),error(2*icomp)
                        !end if
                        ! 2019.04.23, Liu Zhongyin, add azimu while writing
                        if (conjugate) then
                            write(ioDat,'(a8,4es15.6)',iostat=ios) trim(compid),value(2*icomp-1),-value(2*icomp),error(2*icomp),azimu(2*icomp)
                        else
                            write(ioDat,'(a8,4es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp),error(2*icomp),azimu(2*icomp)
                        end if
                        countData = countData + 1
                    end do

                case(Off_Diagonal_Rho_Phase,Phase_Tensor)

                    do icomp = 1,ncomp
                        if (.not. exist(icomp)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        ! For apparent resistivities only, log10 of the values was used internally in the program;
                        ! writing out the linear apparent resistivity
                        if (index(compid,'RHO')>0) then
                            value(icomp) = 10**value(icomp)
                            ! Avoid Inf for FWD calculation
                            if (error(icomp) .ge. LARGE_REAL) then
                                error(icomp) = LARGE_REAL
                            else
                                ! error(icomp) = 10**error(icomp)
                                ! Liu Zhongyin, 2019.10.08, change the error
                                error(icomp) = value(icomp)*error(icomp)*Ln10
                            endif
                        endif

                        ! For Phase only, Rad Used, Writen by LiuZhongyin
                        ! 2017.05.27
                        if (index(compid,'PHS')>0) then
                            if (conjugate) then
                             value(icomp) = value(icomp)*R2D
                            else
                                value(icomp) = -value(icomp)*R2D
                            endif
                        end if
                        write(ioDat,'(es12.6)',    iostat=ios,advance='no') Period
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(a40,3f15.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        !write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(icomp),error(icomp)
                        ! 2019.04.23, Liu Zhongyin, add azimu while writing
                        write(ioDat,'(a8,4es15.6)',iostat=ios) trim(compid),value(icomp),error(icomp),azimu(icomp)
                        countData = countData + 1
                    end do
            end select

        end do  ! transmitters
      end do  ! receivers

      if (output_level > 4) then
        write(0,*) 'Written ',countData,' data values of type MT: ',trim(typeDict(iDt)%name),' to file'
      end if
      deallocate(value, error, exist, STAT=istat)
      ! 2019.04.23, Liu Zhongyin, deallocate azimu
      deallocate(azimu,stat=istat)

    end do WRITE_DATA_TYPE ! data types

    close(ioDat)

   end subroutine write_Z_list


!**********************************************************************
! reads in the ASCII list data file, sets up all dictionaries
! and the allData structure, including data and error bars.
! logic here is quite complicated, but once written can be used
! to read any kind of data, by adding a new case statement.

   subroutine read_Z_list(allData,cfile)

    character(*), intent(in)               :: cfile
    type(dataVectorMTX_t), intent(inout)   :: allData
    ! local variables
    type(dataVectorMTX_t)           :: newData
    integer                         :: nTx,nRx,nDt,ncomp,iRx,iTx,icomp
    integer                         :: countData,countRx
    complex(8), allocatable         :: value(:,:,:) ! (nTx,nRx,ncomp)
    real(8), allocatable            :: error(:,:,:) ! (nTx,nRx,ncomp)
    logical, allocatable            :: exist(:,:,:) ! (nTx,nRx,ncomp)
    integer, allocatable            :: new_TxType(:) ! contains txType indices (nTx)
    integer, allocatable            :: new_Tx(:) ! contains txDict indices (nTx)
    integer, allocatable            :: new_Rx(:) ! contains rxDict indices (nRx)
    character(2)                    :: temp
    character(200)                  :: txTypeName,typeName,typeInfo,typeHeader
    character(40)                   :: siteid,ref_siteid,compid
    integer                         :: nTxt,iTxt,iDt,i,j,k,istat,ios
    character(40)                   :: code,ref_code
    real(8)                         :: x(3),ref_x(3), Period,SI_factor,large
    real(8)                         :: lat,lon,ref_lat,ref_lon
    real(8)                         :: Zreal, Zimag, Zerr
    logical                         :: conjugate, errorBar, isComplex

    ! 2019.03.18, Liu Zhongyin, add azimu variable
    real(8), allocatable            :: azimu(:,:,:) ! (nTx,nRx,ncomp)
    real(8)                         :: angle
    integer                         :: ncount
    character(1000)                 :: tmpline    

    ! First, set up the data type dictionary, if it's not in existence yet
    call setup_typeDict()

    ! Save the user preferences
    nDt = size(typeDict)
    call init_fileInfo(1,nDt)

    ! Now, read the data file
    open(unit=ioDat,file=cfile,form='formatted',status='old')
      
    ! Read the data blocks for each data type
    READ_DATA_TYPE: do
      
    	read(ioDat,'(a2,a200)',iostat=ios) temp,typeInfo
    	read(ioDat,'(a2,a200)',iostat=ios) temp,typeHeader
    	read(ioDat,'(a2,a100)',iostat=ios) temp,typeName

        ! If transmitter name exists, it precedes the typeName
        if (temp(1:1) == '+') then
            txTypeName = typeName
            read(ioDat,'(a2,a100)',iostat=ios) temp,typeName
        else
            txTypeName = 'MT'
        end if
    	iTxt = 1
    	if (ios /= 0) exit
    
    	! Read new data type
    	call compact(typeName)
    	iDt = ImpType(typeName)
    	ncomp = typeDict(iDt)%nComp
    	if (typeDict(iDt)%isComplex) then
        	ncomp = ncomp/2
    	end if

    	fileInfo(iTxt,iDt)%defined = .true.
    	fileInfo(iTxt,iDt)%info_in_file = typeInfo
    	
    	! Sort out the sign convention
		read(ioDat,'(a2,a20)',iostat=ios) temp,fileInfo(iTxt,iDt)%sign_info_in_file
    	if(index(fileInfo(iTxt,iDt)%sign_info_in_file,'-')>0) then
      		fileInfo(iTxt,iDt)%sign_in_file = - 1
    	else
      		fileInfo(iTxt,iDt)%sign_in_file = 1
    	end if
    	if (fileInfo(iTxt,iDt)%sign_in_file == ISIGN) then
      		conjugate = .false.
    	else
      		conjugate = .true.
    	end if

        read(ioDat,'(a2,a20)',iostat=ios) temp,fileInfo(iTxt,iDt)%units_in_file
        SI_factor = ImpUnits(fileInfo(iTxt,iDt)%units_in_file,typeDict(iDt)%units)

        read(ioDat,*,iostat=ios) temp,fileInfo(iTxt,iDt)%geographic_orientation
        read(ioDat,*,iostat=ios) temp,fileInfo(iTxt,iDt)%origin_in_file(1),fileInfo(iTxt,iDt)%origin_in_file(2)
        read(ioDat,*,iostat=ios) temp,nTx,nRx


        if (output_level > 3) then
            write(0,*) node_info,'Reading data type: ',trim(typeName)
            write(0,*) node_info,'Sign convention in file: ',trim(fileInfo(iTxt,iDt)%sign_info_in_file)
            write(0,*) node_info,'Units in file: ',trim(fileInfo(iTxt,iDt)%units_in_file)
            write(0,*) node_info,'Number of transmitters: ',nTx
            write(0,*) node_info,'Number of receivers: ',nRx
        end if


        ! Allocate temporary data arrays
        allocate(new_TxType(nTx),new_Tx(nTx),new_Rx(nRx),STAT=istat)
        allocate(value(nTx,nRx,ncomp),error(nTx,nRx,ncomp),exist(nTx,nRx,ncomp),STAT=istat)

        ! 2019.04.23, Liu Zhongyin, add azimu allocation
        allocate(azimu(nTx,nRx,ncomp),stat=istat)

        new_TxType(:) = 0
        new_Tx(:) = 0
        new_Rx(:) = 0
        value(:,:,:) = dcmplx(0.0d0,0.0d0)
        error(:,:,:) = large
        exist(:,:,:) = .FALSE.
        countData = 0

        ! 2019.04.23, Liu Zhongyin, add azimu initial
        azimu(:,:,:) = R_ZERO


        READ_DATA_LINE: Do

            select case (iDt)

            case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)
                !read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr
                ! 2019.04.23, Liu Zhongyin, add angle while reading line
                !read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr,angle
                read(ioDat,'(a)',iostat=ios) tmpline

                if (ios /= 0) then
                    backspace(ioDat)
                    exit
                end if

                ! Liu Zhongyin, 2019.08.27, add new codes for reading data
                backspace(ioDat)
                call strcount(tmpline, ' ', ncount)
                select case (ncount)
                case(11)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr
                    angle = 0.0
                case(12)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr,angle
                case default
                    call errStop('read datafile error')
                end select

                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! Update the transmitter dictionary and the index (sets up if necessary)
                iTx = update_txDict(Period,2)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a20,2f9.3)') code,lat,lon
                iRx = update_rxDict(x,siteid)

            case(Full_Interstation_TF)
                !read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                !    ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr
                ! 2019.04.23, Liu Zhongyin, add angle while reading line
                ! read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                !     ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr,angle
                read(ioDat,'(a)',iostat=ios) tmpline

                if (ios /= 0) then
                    backspace(ioDat)
                    exit
                end if

                ! Liu Zhongyin, 2019.08.27, add new codes for reading data
                backspace(ioDat)
                call strcount(tmpline, ' ', ncount)
                select case (ncount)
                case(17)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr
                    angle = 0.0
                case(18)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr,angle
                case default
                    call errStop('read datafile error')
                end select

                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! Update the transmitter dictionary and the index (sets up if necessary)
                iTx = update_txDict(Period,2)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a22,2f9.3)') code,lat,lon
                write(ref_siteid,'(a22,2f9.3)') ref_code,ref_lat,ref_lon
                iRx = update_rxDict(x,siteid,ref_x,ref_siteid)


            case(Off_Diagonal_Rho_Phase,Phase_Tensor)
                !read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr
                ! 2019.04.23, Liu Zhongyin, add angle while reading line
                ! read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr,angle
                read(ioDat,'(a)',iostat=ios) tmpline

                if (ios /= 0) then
                    backspace(ioDat)
                    exit
                end if

                ! Liu Zhongyin, 2019.08.27, add new codes for reading data
                backspace(ioDat)
                call strcount(tmpline, ' ', ncount)
                select case (ncount)
                case(10)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr
                    angle = 0.0
                case(11)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr,angle
                case default
                    call errStop('read datafile error')
                end select

                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! For apparent resistivities only, use log10 of the values
                if (index(compid,'RHO')>0) then
                    Zerr  = Zerr/Zreal/Ln10
                    Zreal = log10(Zreal)
            !   Zerr  = log10(Zerr)
            end if

            ! For Phase only, use rad, Writen by LiuZhongyin 2017.05.27
            if (index(compid,'PHS')>0) then
                    if (conjugate) then
                Zreal = Zreal*D2R
                    else
                        Zreal = -Zreal*D2R
                    endif
                Zerr  = Zerr*D2R
                end if

                ! Update the transmitter dictionary and the index (sets up if necessary)
                iTx = update_txDict(Period,2)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a22,2f9.3)') code,lat,lon
                iRx = update_rxDict(x,siteid)

            end select

            ! complete transmitter dictionary update
            do i = 1,nTx
                if ((new_Tx(i) == iTx) .or. (new_Tx(i) == 0)) then
                    exit
                end if
            end do
            new_Tx(i) = iTx
            new_TxType(i) = iTxt

            ! complete receiver dictionary update
            do j = 1,nRx
                if ((new_Rx(j) == iRx) .or. (new_Rx(j) == 0)) then
                    exit
                end if
            end do
            new_Rx(j) = iRx

            ! record the value for storage in the data vector
            if (typeDict(iDt)%isComplex) then
                if (conjugate) then
                    value(i,j,icomp) = SI_factor * dcmplx(Zreal,-Zimag)
                else
                    value(i,j,icomp) = SI_factor * dcmplx(Zreal,Zimag)
                end if
            else
                value(i,j,icomp) = SI_factor * Zreal
            end if
            error(i,j,icomp) = SI_factor * Zerr
            exist(i,j,icomp) = .TRUE.

            ! 2019.04.23, Liu Zhongyin, assign angle to azimu
            ! azimu(i,j,icomp) = angle
            azimu(i,j,:) = angle

            countData = countData + 1

        end do READ_DATA_LINE

        write(0,*) 'Read ',countData,' data values of MT type ',trim(typeDict(iDt)%name),' from file'

        ! Create a single-type data vector from the new values
	call create_dataVectorMTX(nTx,newData)
	newData%allocated = .TRUE.
	errorBar = .TRUE.
        SAVE_DATA: do i = 1,nTx

	       ! Count how many receivers we really have for this transmitter
	       countRx = 0
	       do j = 1,nRx
	        if(count(exist(i,j,:))>0) then
	            countRx = countRx + 1
	        end if
	       end do

	       ! Create a data vector for this transmitter and data type
	       call create_dataVector(1,newData%d(i))
	       newData%d(i)%tx = new_Tx(i)
	       newData%d(i)%txType = new_TxType(i)
	       newData%d(i)%allocated = .TRUE.
	       call create_dataBlock(typeDict(iDt)%nComp,countRx,newData%d(i)%data(1),typeDict(iDt)%isComplex,errorBar)
	       k = 1
	       do j = 1,nRx
	           ! If no data for this receiver, skip it
	           if(count(exist(i,j,:))==0) then
	            cycle
	           end if
	           ! Otherwise, write all components to data vector
	           do icomp = 1,ncomp
	            if(typeDict(iDt)%isComplex) then
	               newData%d(i)%data(1)%value(2*icomp-1,k) = real(value(i,j,icomp))
	               newData%d(i)%data(1)%value(2*icomp  ,k) = imag(value(i,j,icomp))
	               newData%d(i)%data(1)%error(2*icomp-1,k) = error(i,j,icomp)
	               newData%d(i)%data(1)%error(2*icomp  ,k) = error(i,j,icomp)
	               newData%d(i)%data(1)%exist(2*icomp-1,k) = exist(i,j,icomp)
                   newData%d(i)%data(1)%exist(2*icomp  ,k) = exist(i,j,icomp)
                   
                   ! 2019.04.23, Liu Zhongyin, add azimuth
                   newData%d(i)%data(1)%Azimuth(2*icomp-1,k) = azimu(i,j,icomp)
                   newData%d(i)%data(1)%Azimuth(2*icomp  ,k) = azimu(i,j,icomp)
	            else
	               newData%d(i)%data(1)%value(icomp,k) = real(value(i,j,icomp))
	               newData%d(i)%data(1)%error(icomp,k) = error(i,j,icomp)
                   newData%d(i)%data(1)%exist(icomp,k) = exist(i,j,icomp)
                   
                   ! 2019.04.23, Liu Zhongyin, add azimuth
                   newData%d(i)%data(1)%Azimuth(icomp,k) = azimu(i,j,icomp)
	            end if
	           end do
	           newData%d(i)%data(1)%rx(k) = new_Rx(j)
	           k = k+1
	       end do
	       newData%d(i)%data(1)%dataType = iDt
	       newData%d(i)%data(1)%tx = new_Tx(i)
	       newData%d(i)%data(1)%txType = new_TxType(i)
	       newData%d(i)%data(1)%allocated = .TRUE.

        end do SAVE_DATA

	! Merge the new data into the main data vector
    call merge_dataVectorMTX(allData,newData,allData)
    
    ! 2019.04.23, Liu Zhongyin, deallocate azimu
    deallocate(azimu,stat=istat)

	deallocate(value,error,exist,STAT=istat)
	deallocate(new_TxType,new_Tx,new_Rx,STAT=istat)
	call deall_dataVectorMTX(newData)

    end do READ_DATA_TYPE

    close(ioDat)

    ! Finished reading the data: write an empty line to screen
    write(0,*)

    ! Finally, set up the index vectors in the data type dictionary - used for output
    nTxt = 1
    nTx = size(txDict)
    nRx = size(rxDict)
    do iTxt = 1,nTxt
    	do iDt = 1,nDt
		allocate(fileInfo(iTxt,iDt)%tx_index(nTx),STAT=istat)
	        allocate(fileInfo(iTxt,iDt)%dt_index(nTx),STAT=istat)
	        allocate(fileInfo(iTxt,iDt)%rx_index(nTx,nRx),STAT=istat)
	        call index_dataVectorMTX(allData,iTxt,iDt,fileInfo(iTxt,iDt)%tx_index,fileInfo(iTxt,iDt)%dt_index,fileInfo(iTxt,iDt)%rx_index)
	end do
    end do

   end subroutine read_Z_list

end module DataIO
