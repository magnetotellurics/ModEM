! *****************************************************************************
module DataIO
  ! This module contains io routines for reading and writing the data vectors
  ! Version: 3D MT-CSEM-DC-Geoelectric

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

  ! this block of information constitutes user preferences about the data format
  character(200),allocatable, private, save :: info_in_file(:)
  character(20), allocatable, private, save :: sign_info_in_file(:)
  integer,       allocatable, private, save :: sign_in_file(:)
  character(20), allocatable, private, save :: units_in_file(:)
  real,          allocatable, private, save :: origin_in_file(:,:) ! (nDt,2)
  real,          allocatable, private, save :: geographic_orientation(:)

Contains

!**********************************************************************
! writes data in the ASCII list data file; it is convenient to work
! with data ordered by site (NOT by frequency), so we are going into
! some pains here to write them out in this order ...

subroutine write_Z_list(allData,cfile)

!================================================================================
!================================================================ Aihua's version
!================================================================================
!      character(*), intent(in)					:: cfile
!      type(dataVectorMTX_t), intent(in)			:: allData
!      ! local variables
!      type(dataBlock_t)                           :: data
!      character(400)                            :: info
!      integer									:: i,j,k,istat
!      integer                                   :: nBlocks,iTx,iDt
!      real(kind=prec)    						:: SI_factor
!      logical                                   :: conjugate
!      character(10)                             :: tab='           '
!================================================================================
!================================================================ Modular version
!================================================================================
	character(*), intent(in)                  :: cfile
	type(dataVectorMTX_t), intent(in)         :: allData
	! local variables
	integer                         :: nTx,nRx,nDt,ncomp
	integer                         :: countData
	real(8), allocatable            :: value(:) ! (ncomp)
	real(8), allocatable            :: error(:) ! (ncomp)
	logical, allocatable            :: exist(:) ! (ncomp)
	character(2)                    :: temp = '> '
	character(100)                  :: siteid,ref_siteid,compid
	integer                         :: iTx,iRx,iDt,icomp,i,j,k,istat,ios,nBlocks
	real(8)                         :: x(3),ref_x(3), Period,SI_factor,large
	real(8)                         :: lat,lon,ref_lat,ref_lon
	logical                         :: conjugate, isComplex

	!===========================================================================
    !======================================================== New Local Variable
	!===========================================================================
	character(8)                    :: Dipole
	character(40)              		:: Txid=''
    real(8) 						:: Moment, Azi, Dip, LatTx, LongTx, Tx(3)
    type(transmitter_t)					:: aTx
    
    open(unit=ioDat,file=cfile,form='formatted',status='unknown')

    ! For each data type in dictionary, if data of this type exists, write it out.
    WRITE_DATA_TYPE: do iDt = 1,size(typeDict)
    
		nBlocks = countDataBlock(allData,iDt)
		if (nBlocks == 0) then
			! no data for this data type; skip it
			cycle WRITE_DATA_TYPE
		else
			! count the number of transmitters and receivers
			nTx = 0
			nRx = 0
			do iTx = 1,size(txDict)
				if (typeDict(iDt)%tx_index(iTx) > 0) then
					nTx = nTx + 1
				end if
			end do
			do iRx = 1,size(rxDict)
				do iTx = 1,size(txDict)
					if (typeDict(iDt)%rx_index(iTx,iRx) > 0) then
						nRx = nRx + 1
						exit
					end if
				end do
			end do
		end if
		
		! write the data type header
		call compact(info_in_file(iDt))
		write(ioDat,'(a32)',advance='no') '# ModEM impedance responses for '
		write(ioDat,'(a100)',iostat=ios) info_in_file(iDt)
		write(ioDat,'(a200)',iostat=ios) ImpHeader(iDt)
		call compact(typeDict(iDt)%name)
		write(ioDat,'(a2)',advance='no') temp
		write(ioDat,*,iostat=ios) trim(typeDict(iDt)%name)
		call compact(sign_info_in_file(iDt))
		write(ioDat,'(a2)',advance='no') temp
		write(ioDat,*,iostat=ios) trim(sign_info_in_file(iDt))
		call compact(units_in_file(iDt))
		write(ioDat,'(a2)',advance='no') temp
		write(ioDat,*,iostat=ios) trim(units_in_file(iDt))
		write(ioDat,'(a2,f8.2)',iostat=ios) temp,geographic_orientation(iDt)
		write(ioDat,'(a2,2f8.3)',iostat=ios) temp,origin_in_file(iDt,1),origin_in_file(iDt,2)
		write(ioDat,'(a2,2i6)',iostat=ios) temp,nTx,nRx
		
		if (sign_in_file(iDt) == ISIGN) then
			conjugate = .false.
		else if (abs(sign_in_file(iDt)) == 1) then
			conjugate = .true.
		end if
 		SI_factor = ImpUnits(typeDict(iDt)%units,units_in_file(iDt))
		!SI_factor = 1.0d0
		
		ncomp = typeDict(iDt)%nComp
		allocate(value(ncomp),error(ncomp),exist(ncomp),STAT=istat)
		isComplex = typeDict(iDt)%isComplex
		countData = 0
		
		! write data
			do iTx = 1,size(txDict)
		        do iRx = 1,size(rxDict)

		
				k = typeDict(iDt)%rx_index(iTx,iRx)
				i = typeDict(iDt)%dt_index(iTx)
				j = typeDict(iDt)%tx_index(iTx)
				if (k == 0) then
					cycle
				end if
				value = SI_factor * allData%d(j)%data(i)%value(:,k)
				if (allData%d(j)%data(i)%errorBar) then
					error = SI_factor * allData%d(j)%data(i)%error(:,k)
				else
					error = large
				end if
				exist = allData%d(j)%data(i)%exist(:,k)
				Dipole = txDict(iTx)%Dipole
				Period = txDict(iTx)%period
				Moment = txDict(iTx)%moment
				Azi = txDict(iTx)%AzimuthTx
				Dip = txDict(iTx)%dipTx
				Tx = txDict(iTx)%xyzTx
				!Txid = txDict(iTx)%id
				
				siteid = rxDict(iRx)%id
				x = rxDict(iRx)%x
				
				select case (iDt)
case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)

                    do icomp = 1,ncomp/2
                        if (.not. exist(2*icomp-1)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        write(ioDat,'(es12.6)',    iostat=ios,advance='no') Period
                        write(ioDat,'(a40,3f12.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        if (conjugate) then
                            write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),-value(2*icomp),error(2*icomp)
                        else
                            write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp),error(2*icomp)
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
                        write(ioDat,'(a40,3f12.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        write(ioDat,'(a40,3f12.3)',iostat=ios,advance='no') trim(ref_siteid),ref_x(:)
                        if (conjugate) then
                            write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),-value(2*icomp),error(2*icomp)
                        else
                            write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp),error(2*icomp)
                        end if
                        countData = countData + 1
                    end do

case(Off_Diagonal_Rho_Phase,Phase_Tensor)

                    do icomp = 1,ncomp
                        if (.not. exist(icomp)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        write(ioDat,'(es12.6)',    iostat=ios,advance='no') Period
                        write(ioDat,'(a40,3f12.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(icomp),error(icomp)
                        countData = countData + 1
                    end do			
case(Ex_Field,Ey_Field,Bx_Field,By_Field,Bz_Field)
						do icomp = 1,ncomp/2
							if (.not. exist(2*icomp-1)) then
								cycle
							end if
							compid = typeDict(iDt)%id(icomp)
							write(ioDat,'(a8)',iostat=ios,advance='no') trim(Dipole)
							write(ioDat, '(a1)', iostat=ios,advance='no') ' '
							write(ioDat,'(2es12.6)', iostat=ios, advance='no') Period
							write(ioDat, '(a1)', iostat=ios,advance='no') ' '
							write(ioDat,'(2es12.6)', iostat=ios, advance='no') Moment
							write(ioDat,'(2f9.3, 3f12.3)',iostat=ios,advance='no') Azi, Dip, Tx
							write(ioDat,'(a15,3f12.3)',iostat=ios,advance='no') trim(siteid),x(:)
							write(ioDat, '(a1)', iostat=ios,advance='no') ' '
!							write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp), error(2*icomp)
!							error(2*icomp) = sqrt(value(2*icomp-1)*value(2*icomp-1) + value(2*icomp)*value(2*icomp))/100.0;
							write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp), error(2*icomp)
							countData = countData + 1
                        end do
case(Pole_Pole_DC_Rho)
						do icomp = 1,1
							if (.not. exist(2*icomp-1)) then
								cycle
							end if
							compid = typeDict(iDt)%id(icomp)
							write(ioDat,'(3f12.3)',iostat=ios,advance='no') Tx
							write(ioDat,'(a15,3f12.3)',iostat=ios,advance='no') trim(siteid),x(:)
							write(ioDat, '(a1)', iostat=ios,advance='no') ' '
!							write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp), error(2*icomp)
!							error(2*icomp) = sqrt(value(2*icomp-1)*value(2*icomp-1) + value(2*icomp)*value(2*icomp))/100.0;
							write(ioDat,'(a16,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1), value(2*icomp-1)*0.05 !error(2*icomp)
							countData = countData + 1
						end do
end select
		
			end do  ! transmitters
		end do  ! receivers
		
		if (output_level > 4) then
		write(0,*) 'Written ',countData,' data values of type ',trim(typeDict(iDt)%name),' to file'
		end if
		deallocate(value, error, exist, STAT=istat)
		
    end do WRITE_DATA_TYPE

	close(ioDat)

end subroutine write_Z_list

!******************************************************************
! reads in the ASCII data file, sets up all dictionaries
! and the allData structure, including data and error bars

subroutine read_Z_list(allData,cfile)

!================================================================================
!================================================================ Aihua's version
!================================================================================
!     character(*), intent(in)  				:: cfile
!     type(dataVectorMTX_t), intent(inout)   	:: allData
!     ! local variables
!     type(dataBlock_t), dimension(:), pointer :: Data
!     integer      							:: nBlocks,nTx,nDt,nSite,nComp,nData
!     integer                                :: iTx,iRx,iDt,i,j,k,istat
!     character(10), pointer,dimension(:)    :: siteIDs
!     real(kind=prec), pointer, dimension(:,:) :: siteLoc
!     real(kind=prec)                        :: SI_factor, Period, tx, ty, tz
!     logical      							:: conjugate, errorBar, isComplex
!     character(400) 						:: temp,header
!     type(MTtx)                            :: aTx

!================================================================================
!================================================================ Modular version
!================================================================================
    character(*), intent(in)               :: cfile
    type(dataVectorMTX_t), intent(inout)   :: allData
    ! local variables
    type(dataVectorMTX_t)           :: newData
    integer                         :: nTx,nRx,nDt,ncomp,iRx,iTx,icomp
    integer                         :: countData,countRx
    complex(8), allocatable         :: value(:,:,:) ! (nTx,nRx,ncomp)
    real(8), allocatable            :: error(:,:,:) ! (nTx,nRx,ncomp)
    logical, allocatable            :: exist(:,:,:) ! (nTx,nRx,ncomp)
    integer, allocatable            :: new_Tx(:) ! contains txDict indices (nTx)
    integer, allocatable            :: new_Rx(:) ! contains rxDict indices (nRx)
    character(2)                    :: temp
    character(200)                  :: typeName,typeInfo,typeHeader
    character(40)                   :: siteid,ref_siteid,compid,Txid
    integer                         :: iDt,i,j,k,istat,ios
    character(12)                   :: code,ref_code
    real(8)                         :: x(3),ref_x(3), Period,SI_factor,large
    real(8)                         :: lat,lon,ref_lat,ref_lon
    real(8)                         :: Zreal, Zimag, Zerr
    logical                         :: conjugate, errorBar, isComplex
    
    !===========================================================================
    !======================================================== New Local Variable
    !===========================================================================
    character(8)                    :: Dipole
    real(8) 						:: Moment, Azi, Dip, LatTx, LongTx, Tx(3)
    type(transmitter_t)					:: aTx
    
	! First, set up the data type dictionary, if it's not in existence yet
	call setup_typeDict()
	
	! Save the user preferences
	nDt = size(typeDict)
	allocate(info_in_file(nDt),STAT=istat)
	allocate(sign_info_in_file(nDt),STAT=istat)
	allocate(sign_in_file(nDt),STAT=istat)
	allocate(units_in_file(nDt),STAT=istat)
	allocate(origin_in_file(nDt,2),STAT=istat)
	allocate(geographic_orientation(nDt),STAT=istat)
		
    ! Now, read the data file
    open(unit=ioDat,file=cfile,form='formatted',status='old')
      
    ! Read the data blocks for each data type
    READ_DATA_TYPE: do
      
    	read(ioDat,'(a2,a200)',iostat=ios) temp,typeInfo
    	read(ioDat,'(a2,a200)',iostat=ios) temp,typeHeader
    	read(ioDat,'(a2,a100)',iostat=ios) temp,typeName
    	if (ios /= 0) exit
    
    	! Read new data type
    	call compact(typeName)
    	iDt = ImpType(typeName)
    	ncomp = typeDict(iDt)%nComp
    	if (typeDict(iDt)%isComplex) then
        	ncomp = ncomp/2
    	end if
    	info_in_file(iDt) = typeInfo
    	
    	! Sort out the sign convention
		read(ioDat,'(a2,a20)',iostat=ios) temp,sign_info_in_file(iDt)
    	if(index(sign_info_in_file(iDt),'-')>0) then
      		sign_in_file(iDt) = - 1
    	else
      		sign_in_file(iDt) = 1
    	end if
    	if (sign_in_file(iDt) == ISIGN) then
      		conjugate = .false.
    	else
      		conjugate = .true.
    	end if

   	read(ioDat,'(a2,a20)',iostat=ios) temp,units_in_file(iDt)
   	SI_factor = ImpUnits(units_in_file(iDt),typeDict(iDt)%units)

	read(ioDat,*,iostat=ios) temp,geographic_orientation(iDt)
	read(ioDat,*,iostat=ios) temp,origin_in_file(iDt,1),origin_in_file(iDt,2)
    read(ioDat,*,iostat=ios) temp,nTx,nRx


    	
    	if (output_level > 3) then
        	write(0,*) node_info,'Reading data type: ',trim(typeName)
        	write(0,*) node_info,'Sign convention in file: ',trim(sign_info_in_file(iDt))
        	write(0,*) node_info,'Units in file: ',trim(units_in_file(iDt))
        	write(0,*) node_info,'Number of transmitters: ',nTx
        	write(0,*) node_info,'Number of receivers: ',nRx
    	end if
    	

    	! Allocate temporary data arrays
    	allocate(new_Tx(nTx),new_Rx(nRx),STAT=istat)
    	allocate(value(nTx,nRx,ncomp),error(nTx,nRx,ncomp),exist(nTx,nRx,ncomp),STAT=istat)
    	
    	new_Tx(:) = 0
    	new_Rx(:) = 0
    	value(:,:,:) = dcmplx(0.0d0,0.0d0)
    	error(:,:,:) = large
	    exist(:,:,:) = .FALSE.
	countData = 0

	
      	select case (iDt)
		    case (Pole_Pole_DC_Rho)
               Do
				read(ioDat,*,iostat=ios) Tx(1), Tx(2), Tx(3), code, x(1), x(2), x(3), compid, Zreal,Zerr
				if (ios /= 0 ) then
					backspace(ioDat)
					exit
				end if
				! Find component id for this value
            			icomp = ImpComp(compid,iDt)
            		    aTx%Tx_type='DC'
						aTx%nPol=1
            			aTx%Dipole = ''
				        aTx%period = 1
            			aTx%omega = 1		
            			aTx%xyzTx = Tx
            			aTx%azimuthTx = 0
            			aTx%dipTx = 0
            			atx%moment = 0
            		
!            			write(Txid,'(2f9.3)') LatTx,LongTx
!            			aTx%id = Txid
            		
            			iTx = update_txDict(aTx)
            	do i = 1,nTx
					if ((new_Tx(i) == iTx) .or. (new_Tx(i) == 0)) then
						exit
					end if
				end do
				new_Tx(i) = iTx

				! Update the receiver dictionary and index (sets up if necessary)
				! For now, make lat & lon part of site ID; could use directly in the future
				write(siteid,'(a12,2f9.3)') code
				iRx = update_rxDict(x,siteid)
				do j = 1,nRx
					if ((new_Rx(j) == iRx) .or. (new_Rx(j) == 0)) then
						exit
					end if
				end do
				new_Rx(j) = iRx					
				value(i,j,icomp) = SI_factor * Zreal
				error(i,j,icomp) = SI_factor * Zerr
				exist(i,j,icomp) = .TRUE.
				countData = countData + 1
			end Do
				write(0,*) nTx, nRx
				
      		case(Ex_Field,Ey_Field,Bx_Field,By_Field,Bz_Field)
      			Do
				read(ioDat,*,iostat=ios) Dipole, Period, Moment, Azi, Dip, Tx(1), Tx(2), Tx(3), code, x(1), x(2), x(3), compid, Zreal, Zimag, Zerr
				if (ios /= 0 ) then
					backspace(ioDat)
					exit
				end if
				! Find component id for this value
            			icomp = ImpComp(compid,iDt)
            		    aTx%Tx_type='CSEM'
						aTx%nPol=1
            			aTx%Dipole = Dipole
				        aTx%period = Period
            			aTx%omega = 2.0d0*PI/Period            		
            			aTx%xyzTx = Tx
            			aTx%azimuthTx = Azi
            			aTx%dipTx = Dip
            			atx%moment = Moment
            		
!            			write(Txid,'(2f9.3)') LatTx,LongTx
!            			aTx%id = Txid
            		
            			iTx = update_txDict(aTx)
            			do i = 1,nTx
					if ((new_Tx(i) == iTx) .or. (new_Tx(i) == 0)) then
						exit
					end if
				end do
				new_Tx(i) = iTx

				! Update the receiver dictionary and index (sets up if necessary)
				! For now, make lat & lon part of site ID; could use directly in the future
				write(siteid,'(a12,2f9.3)') code
				iRx = update_rxDict(x,siteid)
				do j = 1,nRx
					if ((new_Rx(j) == iRx) .or. (new_Rx(j) == 0)) then
						exit
					end if
				end do
				new_Rx(j) = iRx					
				value(i,j,icomp) = SI_factor * dcmplx(Zreal,Zimag)
				error(i,j,icomp) = SI_factor * Zerr
				exist(i,j,icomp) = .TRUE.
				countData = countData + 1
			end Do
				write(*,*) nTx, nRx
       case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)

        do

            read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr

            if (ios /= 0) then
                backspace(ioDat)
                exit
            end if

            ! Find component id for this value
            icomp = ImpComp(compid,iDt)
            			icomp = ImpComp(compid,iDt)
            		    aTx%Tx_type='MT'
						aTx%nPol=2
            			aTx%Dipole =''
				        aTx%period = Period
            			aTx%omega = 2.0d0*PI/Period            		
            			aTx%xyzTx(1:3) = 0
            			aTx%azimuthTx = 0
            			aTx%dipTx = 0
            			atx%moment = 0
            		
!            			write(Txid,'(2f9.3)') LatTx,LongTx
!            			aTx%id = Txid
            		
            			iTx = update_txDict(aTx)
            do i = 1,nTx
                if ((new_Tx(i) == iTx) .or. (new_Tx(i) == 0)) then
                    exit
                end if
            end do
            new_Tx(i) = iTx

            ! Update the receiver dictionary and index (sets up if necessary)
            ! For now, make lat & lon part of site ID; could use directly in the future
            write(siteid,'(a12,2f9.3)') code,lat,lon
            iRx = update_rxDict(x,siteid)
            do j = 1,nRx
                if ((new_Rx(j) == iRx) .or. (new_Rx(j) == 0)) then
                    exit
                end if
            end do
            new_Rx(j) = iRx

            if (conjugate) then
                value(i,j,icomp) = SI_factor * dcmplx(Zreal,-Zimag)
            else
                value(i,j,icomp) = SI_factor * dcmplx(Zreal,Zimag)
            end if
            error(i,j,icomp) = SI_factor * Zerr
            exist(i,j,icomp) = .TRUE.

            countData = countData + 1

        end do


       case(Full_Interstation_TF)

        do

            read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr

            if (ios /= 0) then
                backspace(ioDat)
                exit
            end if

            ! Find component id for this value
            icomp = ImpComp(compid,iDt)
            		    aTx%Tx_type='MT'
						aTx%nPol=2
            			aTx%Dipole =''
				        aTx%period = Period
            			aTx%omega = 2.0d0*PI/Period            		
            			aTx%xyzTx(1:3) = 0
            			aTx%azimuthTx = 0
            			aTx%dipTx = 0
            			atx%moment = 0
            		
!            			write(Txid,'(2f9.3)') LatTx,LongTx
!            			aTx%id = Txid
            		
            			iTx = update_txDict(aTx)
            do i = 1,nTx
                if ((new_Tx(i) == iTx) .or. (new_Tx(i) == 0)) then
                    exit
                end if
            end do
            new_Tx(i) = iTx

            ! Update the receiver dictionary and index (sets up if necessary)
            ! For now, make lat & lon part of site ID; could use directly in the future
            write(siteid,'(a12,2f9.3)') code,lat,lon
            write(ref_siteid,'(a12,2f9.3)') ref_code,ref_lat,ref_lon
            iRx = update_rxDict(x,siteid,ref_x,ref_siteid)
            do j = 1,nRx
                if ((new_Rx(j) == iRx) .or. (new_Rx(j) == 0)) then
                    exit
                end if
            end do
            new_Rx(j) = iRx

            if (conjugate) then
                value(i,j,icomp) = SI_factor * dcmplx(Zreal,-Zimag)
            else
                value(i,j,icomp) = SI_factor * dcmplx(Zreal,Zimag)
            end if
            error(i,j,icomp) = SI_factor * Zerr
            exist(i,j,icomp) = .TRUE.

            countData = countData + 1

        end do


       case(Off_Diagonal_Rho_Phase,Phase_Tensor)

        do

            read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr

            if (ios /= 0) then
                backspace(ioDat)
                exit
            end if

            ! Find component id for this value
            icomp = ImpComp(compid,iDt)
            		    aTx%Tx_type='MT'
						aTx%nPol=2
            			aTx%Dipole =''
				        aTx%period = Period
            			aTx%omega = 2.0d0*PI/Period            		
            			aTx%xyzTx(1:3) = 0
            			aTx%azimuthTx = 0
            			aTx%dipTx = 0
            			atx%moment = 0
            		
!            			write(Txid,'(2f9.3)') LatTx,LongTx
!            			aTx%id = Txid
            		
            			iTx = update_txDict(aTx)
            do i = 1,nTx
                if ((new_Tx(i) == iTx) .or. (new_Tx(i) == 0)) then
                    exit
                end if
            end do
            new_Tx(i) = iTx

            ! Update the receiver dictionary and index (sets up if necessary)
            ! For now, make lat & lon part of site ID; could use directly in the future
            write(siteid,'(a12,2f9.3)') code,lat,lon
            iRx = update_rxDict(x,siteid)
            do j = 1,nRx
                if ((new_Rx(j) == iRx) .or. (new_Rx(j) == 0)) then
                    exit
                end if
            end do
            new_Rx(j) = iRx

            value(i,j,icomp) = SI_factor * Zreal
            error(i,j,icomp) = SI_factor * Zerr
            exist(i,j,icomp) = .TRUE.
       
            countData = countData + 1

        end do			
		end select
      
		write(0,*) 'Read ',countData,' data values of type ',trim(typeDict(iDt)%name),' from file'
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
				else
				   newData%d(i)%data(1)%value(icomp,k) = real(value(i,j,icomp))
				   newData%d(i)%data(1)%error(icomp,k) = error(i,j,icomp)
				   newData%d(i)%data(1)%exist(icomp,k) = exist(i,j,icomp)
				end if
			   end do
			   newData%d(i)%data(1)%rx(k) = new_Rx(j)
			   k = k+1
		   end do
		   newData%d(i)%data(1)%dataType = iDt
		   newData%d(i)%data(1)%tx = new_Tx(i)
		   newData%d(i)%data(1)%allocated = .TRUE.
	
		end do SAVE_DATA
		
		! Merge the new data into the main data vector
		call merge_dataVectorMTX(allData,newData,allData)

		deallocate(value,error,exist,STAT=istat)
		deallocate(new_Tx,new_Rx,STAT=istat)
		call deall_dataVectorMTX(newData)
	
	end do READ_DATA_TYPE
	
	close(ioDat)

	! Finally, set up the index vectors in the data type dictionary - used for output
	nTx = size(txDict)
	nRx = size(rxDict)
	do iDt = 1,nDt
		allocate(typeDict(iDt)%tx_index(nTx),STAT=istat)
		allocate(typeDict(iDt)%dt_index(nTx),STAT=istat)
		allocate(typeDict(iDt)%rx_index(nTx,nRx),STAT=istat)
		call index_dataVectorMTX(allData,iDt,typeDict(iDt)%tx_index,typeDict(iDt)%dt_index,typeDict(iDt)%rx_index)
	end do
	
end subroutine read_Z_list	  ! Finally, set up allData



end module DataIO
