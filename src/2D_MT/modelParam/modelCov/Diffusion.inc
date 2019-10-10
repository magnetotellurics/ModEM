!-------------------------------------------------------------------!
! 2D_MT model covariance by Weerachai Siripunvaraporn: procedures.  !
!-------------------------------------------------------------------!
  function multBy_Cm(m1) result (m2)

   ! Multiplies by the square root of the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).

    type (modelParam_t), intent(in)     :: m1
    type (modelParam_t)				    :: m2

    if (.not. CmSqrt%allocated) then
    	call create_CmSqrt(m1)
    end if

    m2 = m1
    call solveDiff(CmSqrt,m2)
    m2%temporary = .true.

  end function multBy_Cm
! *******************************************************************
  function multBy_CmSqrt(m1) result (m2)

   ! Multiplies by the square root of the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).

    type (modelParam_t), intent(in)     :: m1
    type (modelParam_t)				    :: m2

    if (.not. CmSqrt%allocated) then
    	call create_CmSqrt(m1)
    end if

    m2 = m1
    call solveDiff(CmSqrt,m2)
    m2%temporary = .true.

  end function multBy_CmSqrt

! *******************************************************************
  subroutine create_CmSqrt(m,cfile)

    type (modelParam_t), intent(in)     :: m
    character(*), intent(in), optional  :: cfile

    if (present(cfile)) then
      ! do nothing: a placeholder for other covariance operators
    end if

    call setup1DCM(m%grid,CmSqrt,cfile)

    CmSqrt%allocated = .true.

  end subroutine create_CmSqrt

! *******************************************************************
  subroutine deall_CmSqrt()

    integer                             :: istat

    deallocate(CmSqrt%YDIF, STAT = istat)
    deallocate(CmSqrt%ZDIF, STAT = istat)

  end subroutine deall_CmSqrt

! *******************************************************************
!
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.
!
!----------------------------------------------------------------------
!
!     SOLVING DIFFUSION EQU. IN EXPLICIT WAY.
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.


  subroutine setup1DCM(gd,Cm,cfile)

  type (grid_t),   intent(in)    :: gd
  type (modelCov_t), intent(out)   :: Cm
  character(*), intent(in), optional  :: cfile

  real (kind=prec), pointer, dimension(:,:)   ::  YDF
  real (kind=prec), pointer, dimension(:,:)   ::  ZDF

  integer  :: iy, iz, NzEarth ,Ny, istat
  real (kind=prec) :: ylen, zlen

  Ny      = gd%Ny
  NzEarth = gd%Nz - gd%Nza
  allocate(Cm%YDIF(2,gd%Ny,NzEarth),STAT = istat)
  allocate(Cm%ZDIF(2,gd%Ny,NzEarth),STAT = istat)

  	allocate(Cm%Sx(NzEarth), STAT=istat)
	allocate(Cm%Sy(NzEarth), STAT=istat)
	allocate(Cm%mask(Ny,NzEarth), STAT=istat)

  Cm%ylen = 0.05
  Cm%zlen = 0.05
  Cm%N=10
  !By default all model parameters as free for changes
  Cm%mask=ONE
  
  Cm%allocated=.true.
  
  ! If model Cov file is present ---> read and overwrite Cm%ylen and Cm%zlen from the file and also the masked cells  
  if (present(cfile)) then
    call read_Cm(cfile,CM)
  end if

  ylen = Cm%ylen
  zlen = Cm%zlen

  ! horizontal diffusion
  do iz = 1,NzEarth
    call setupYCM(Ny,ylen,YDF,cm,iz,gd)
    do iy = 1,Ny
      Cm%YDIF(1,iy,iz) = YDF(1,iy)
      Cm%YDIF(2,iy,iz) = YDF(2,iy)
    enddo
  	deallocate(YDF,STAT=istat)
  enddo

  ! vertical diffusion
  do iy = 1,Ny
    call setupZCM(NzEarth,zlen,ZDF,cm,iy,gd)
    do iz = 1,NzEarth
      Cm%ZDIF(1,iy,iz) = ZDF(1,iz)
      Cm%ZDIF(2,iy,iz) = ZDF(2,iz)
    enddo
    deallocate(ZDF,STAT=istat)
  enddo

  end subroutine setup1DCM

!------------------------------------------------------------------------------
 subroutine read_Cm(cfile,Cm)

	! The minimal covariance information includes the AR parameters
	! alpha(k), beta(k) for smoothing in x, y directions and gamma for
	! the vertical smoothing. Both alpha and beta could depend on the
	! vertical layer. The scaling is the identity when not specified.
	! This information is read from a file. Also, we read an integer mask
	! array that subdivides the model grid into different regions
	! (AIR, OCEAN, EARTH) and a set of rules that overrides the default
	! smoothing parameters across a particular surface between two
	! distinct regions. We use this to set up the covariance CmSqrt.
	!
	! Strictly speaking, to define the smoothing across surfaces in
	! full generality while maintaining efficiency, it has to be a sparse
	! real vector defined on FACES (sparsevecr). We only have a complex
	! sparse vector implemented (sparsevecc). We could either use that,
	! or imitate the structure.

    character(*), intent(in)                         :: cfile
	type (modelCov_t), intent(inout)   :: Cm
    ! Exception rules
    integer, pointer, dimension(:)                   :: mask1, mask2, ii, jj, kk, xyz
    real (kind=prec), pointer, dimension(:)  :: smoothing, S

    integer       			 :: Nx, Ny, NzEarth, nrules, nS, i, j, k, n, istat
    integer                  :: fid = 30
	integer                  :: ios

	if (.not. Cm%allocated) then
		call errStop('Model covariance must be allocated before reading from file in read_CmSqrt')
	end if

    open(unit=fid,file=cfile,form='formatted',status='old')

    ! skip the 16 lines header
    do j = 1,16
       read(fid,*)
    end do

	! read grid dimensions
	read(fid,*) Ny,NzEarth
	Cm%Ny = Ny
	Cm%NzEarth = NzEarth

	! read smoothing parameters
    read(fid,*) Cm%Sy
    read(fid,*) Cm%Sz
    
 Cm%ylen=Cm%Sy(1)
 Cm%zlen=Cm%Sz
 
	! read number of times to apply the smoothing
	read(fid,*) Cm%N

	! read exception rules for smoothing across surfaces
    read(fid,*) nrules
    allocate(mask1(nrules),mask2(nrules),smoothing(nrules),STAT=istat)
    do n = 1,nrules
       read(fid,*) mask1(n), mask2(n), smoothing(n)
    end do

	! read the mask array
      do j=1,NzEarth
      	read(fid,*,iostat=ios) (Cm%mask(k,j), k=1,Ny)
      end do

	close(fid)



       
  end subroutine read_Cm
! *******************************************************************
  subroutine setupYCM(Ny,ylen,YDF,Cm,iz,gd)

  !use wsLAPACK
  real (kind=prec), pointer, dimension(:,:)   ::  YDF
  integer  :: Ny,iz
  real (kind=prec) :: ylen
  type (modelCov_t), intent(in)   :: Cm
    type (grid_t),   intent(in)    :: gd

  integer  :: iy,info,istat
  real (kind=prec) :: cr, cl,min_scale
  real (kind=prec), allocatable, dimension(:)   ::  dfyc, dfyr


  allocate(YDF(2,Ny),STAT=istat)
  allocate(dfyc(Ny))
  allocate(dfyr(Ny))
  do iy = 1,Ny
    dfyc(iy) = 0.
    dfyr(iy) = 0.
  enddo

  do iy = 1,Ny
    cr = ylen
    cl = ylen
	    
		IF (iy.LT.Ny) THEN
          IF ((cm%mask(iy,iz) .eq. OCEAN).AND.(cm%mask(iy+1,iz).eq. OCEAN)) THEN
            cr = R_ZERO
          ENDIF
        ENDIF

		IF (iy.LT.Ny) THEN
          IF ((cm%mask(iy,iz) .eq. AIR).AND.(cm%mask(iy+1,iz).eq. AIR)) THEN
            cr = R_ZERO
          ENDIF
        ENDIF




        IF (iy.GT.1) THEN
          IF ((cm%mask(iy,iz) .eq. OCEAN).AND.(cm%mask(iy-1,iz).eq. OCEAN)) THEN
            cl = R_ZERO
          ENDIF
        ENDIF


        IF (iy.GT.1) THEN
          IF ((cm%mask(iy,iz) .eq. AIR).AND.(cm%mask(iy-1,iz).eq. AIR)) THEN
            cl = R_ZERO
          ENDIF
        ENDIF





    dfyc(iy) = 1. + cr + cl
    dfyr(iy) = -cr
  enddo

  do iy =1,Ny
    YDF(1,iy) = 0.
    YDF(2,iy) = 0.
  enddo

  do iy =2,Ny
    YDF(1,iy) = dfyr(iy-1)
  enddo
  do iy =1,Ny
    YDF(2,iy) = dfyc(iy)
  enddo

  call dpbtrf('U',Ny,1,YDF,2,info)

  deallocate(dfyc)
  deallocate(dfyr)

  end subroutine setupYCM

!------------------------------------------------------------------------------

  subroutine setupZCM(Nzb,zlen,ZDF,cm,iy,gd)

  !use wsLAPACK
  real (kind=prec), pointer,  dimension(:,:)   ::  ZDF
  integer  :: Nzb,iy
  real (kind=prec) :: zlen
    type (modelCov_t), intent(in)   :: Cm
	type (grid_t),   intent(in)    :: gd

  integer  :: iz,info,istat
  real (kind=prec) :: cr, cl,min_scale
  real (kind=prec), allocatable, dimension(:)   ::  dfzc, dfzr


  allocate(ZDF(2,Nzb),STAT=istat)
  allocate(dfzc(Nzb))
  allocate(dfzr(Nzb))
  do iz = 1,Nzb
    dfzc(iz) = 0.
    dfzr(iz) = 0.
  enddo

  do iz = 1,Nzb
    cr = zlen
    cl = zlen

        IF (iz.LT.Nzb) THEN
          IF ((cm%mask(iy,iz) .eq. OCEAN).AND.(cm%mask(iy,iz+1).eq. OCEAN)) THEN
            cr = R_ZERO
          ENDIF
        ENDIF

        IF (iz.LT.Nzb) THEN
          IF ((cm%mask(iy,iz) .eq. AIR).AND.(cm%mask(iy,iz+1).eq. AIR)) THEN
            cr = R_ZERO
          ENDIF
        ENDIF

        IF (iz.GT.1) THEN
          IF ((cm%mask(iy,iz) .eq. OCEAN).AND.(cm%mask(iy,iz-1).eq. OCEAN)) THEN
            cl = R_ZERO
          ENDIF
        ENDIF

        IF (iz.GT.1) THEN
          IF ((cm%mask(iy,iz) .eq. AIR).AND.(cm%mask(iy,iz-1).eq. AIR)) THEN
            cl = R_ZERO
          ENDIF
        ENDIF


    dfzc(iz) = 1. + cr + cl
    dfzr(iz) = -cr

  enddo

  do iz =1,Nzb
    ZDF(1,iz) = 0.
    ZDF(2,iz) = 0.
  enddo

  do iz =2,Nzb
    ZDF(1,iz) = dfzr(iz-1)
  enddo
  do iz =1,Nzb
    ZDF(2,iz) = dfzc(iz)
  enddo

  call dpbtrf('U',Nzb,1,ZDF,2,info)

  deallocate(dfzc)
  deallocate(dfzr)

  end subroutine setupZCM

!------------------------------------------------------------------------------

  subroutine solveDiff(Cm,CmJ)

  type (modelCov_t),   intent(in)    :: Cm
  type (modelParam_t), intent(inout) :: CmJ

  integer            :: Ny,Nz,it, dTime,iy,iz

  Ny = CmJ%Ny
  Nz = CmJ%NzEarth

  dTime = Cm%N

 !do iz = 1,Nz
 ! do iy = 1,Ny
 !  write(45,*) CmJ%v(iy,iz)
 ! enddo
 !enddo

  do it = 1,dTime
    call yDiff(Cm,CmJ)
    call zDiff(Cm,CmJ)
  enddo
  do it = 1,dTime
    call zDiff(Cm,CmJ)
    call yDiff(Cm,CmJ)
  enddo

 do iz = 1,Nz
   do iy = 1,Ny
     if (cm%mask(iy,iz) == OCEAN .or. cm%mask(iy,iz) == AIR) then
	   CmJ%v(iy,iz)=R_ZERO
     end if
   enddo
 enddo

  end subroutine solveDiff

!------------------------------------------------------------------------------

  subroutine yDiff(Cm,CmJ)

  !use wsLAPACK
  type (modelCov_t),   intent(in)    :: Cm
  type (modelParam_t), intent(inout) :: CmJ

  real (kind=prec), allocatable, dimension(:,:)   ::  YDF
  real (kind=prec), allocatable, dimension(:)     ::  uy

  integer            :: Ny,Nz,iy,iz,info,istat

  Ny = CmJ%Ny
  Nz = CmJ%NzEarth

  allocate(YDF(2,Ny))
  allocate(uy(Ny))

  do iz = 1,Nz
    do iy = 1,Ny
      uy(iy)    = CmJ%v(iy,iz)
      YDF(1,iy) = Cm%YDIF(1,iy,iz)
      YDF(2,iy) = Cm%YDIF(2,iy,iz)
    enddo

    call dpbtrs('U',Ny,1,1,YDF,2,uy,Ny,info)

    do iy = 1,Ny
      CmJ%v(iy,iz) = uy(iy)
    enddo

  enddo

  deallocate(YDF,STAT=istat)
  deallocate(uy,STAT=istat)

  end subroutine yDiff

!------------------------------------------------------------------------------

  subroutine zDiff(Cm,CmJ)

  !use wsLAPACK
  type (modelCov_t),   intent(in)    :: Cm
  type (modelParam_t), intent(inout) :: CmJ

  real (kind=prec), allocatable, dimension(:,:)   ::  ZDF
  real (kind=prec), allocatable, dimension(:)     ::  uz


  integer            :: Ny,Nz,iy,iz,info,istat

  Ny = CmJ%Ny
  Nz = CmJ%NzEarth

  allocate(ZDF(2,Nz))
  allocate(uz(Nz))

  do iy = 1,Ny
    do iz = 1,Nz
      uz(iz)    = CmJ%v(iy,iz)
      ZDF(1,iz) = Cm%ZDIF(1,iy,iz)
      ZDF(2,iz) = Cm%ZDIF(2,iy,iz)
    enddo

    call dpbtrs('U',Nz,1,1,ZDF,2,uz,Nz,info)

    do iz = 1,Nz
      CmJ%v(iy,iz) = uz(iz)
    enddo

  enddo

  deallocate(ZDF,STAT=istat)
  deallocate(uz,STAT=istat)

  end subroutine zDiff
