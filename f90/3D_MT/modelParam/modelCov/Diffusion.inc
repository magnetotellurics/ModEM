!-------------------------------------------------------------------!
! 3D_MT model covariance by Weerachai Siripunvaraporn: procedures.  !
!                                                                   !
! Converted to Fortran 90 and adapted to the modular structure by   !
! Anna Kelbert  (March 24, 2008).                                   !
!                                                                   !
! Weerachai uses the prior model (integer format) to switch off     !
! the model covariance for frozen regions in Setup1DCM.             !
! Judging from the parameters in the Setup subroutines, he has      !
! replaced the explicit time forwards scheme with the implicit time !
! backwards scheme (see solution of the diffusion equation and      !
! Crank-Nicholson scheme in the Numerical Recipes for a discussion).!
!-------------------------------------------------------------------!
! *******************************************************************
  function multBy_Cm(m1) result (m2)
   ! Multiplies by the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m \tilde{m} + m_0.
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).

    type (modelParam_t), intent(in)     :: m1
    type (modelParam_t)                 :: m2
    type (rscalar)                      :: logcond
    character(80)                       :: paramType
    Integer dtime



    if (.not. CmSqrt%allocated) then
    	call create_CmSqrt(m1)
    end if
    m2=m1
    paramType = LOGE

    call getValue_modelParam(m1,paramType,logcond)

    dtime=5

    call solveDiff(dtime,CmSqrt,logcond)

    call setValue_modelParam(m2,paramType,logcond)

 end function multBy_Cm

! *******************************************************************

! *******************************************************************
  function multBy_CmSqrt(m1) result (m2)
   ! Multiplies by the square root of the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).

    type (modelParam_t), intent(in)     :: m1
    type (modelParam_t)                 :: m2
    type (rscalar)                      :: logcond
    character(80)                       :: paramType
    Integer dtime

    if (.not. CmSqrt%allocated) then
    	call create_CmSqrt(m1)
    end if
    m2=m1

    paramType = LOGE
    call getValue_modelParam(m1,paramType,logcond)

    dtime=1
    call solveDiff(dtime,CmSqrt,logcond)
    call setValue_modelParam(m2,paramType,logcond)

 end function multBy_CmSqrt

! *******************************************************************
  function multBy_CmSqrtInv(dm) result (mhat)

   ! Multiplies by the inverse square root of the model covariance,
   ! which is viewed as a roughening operator. Intended
   ! to be used to compute \tilde{m} = C_m^{-1/2} ( m - m_0 ).
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).

    type (modelParam_t), intent(in)     :: dm
    type (modelParam_t)                 :: mhat

    if (.not. CmSqrt%allocated) then
        call create_CmSqrt(dm)
    end if

    call errStop('The inverse has not been implemented for the diffusive covariance')

  end function multBy_CmSqrtInv

! *******************************************************************
  subroutine create_CmSqrt(m,cfile)
	! Initializes CmSqrt variable stored in RecursiveAR.hd. If cfile
  	! is specified, gets this information from file.

    type (modelParam_t), intent(in)     :: m
    character(*), intent(in), optional  :: cfile
    integer                             :: istat
    logical                             :: exists

	! initializing CmSqrt
	CmSqrt%Nx = m%Nx
	CmSqrt%Ny = m%Ny
	CmSqrt%NzEarth = m%NzEarth
	allocate(CmSqrt%Sx(CmSqrt%NzEarth), STAT=istat)
	allocate(CmSqrt%Sy(CmSqrt%NzEarth), STAT=istat)
	CmSqrt%Sx = 0.5
	CmSqrt%Sy = 0.5
	CmSqrt%Sz = 0.5
	 CmSqrt%xlen=0.5
     CmSqrt%ylen=0.5
     CmSqrt%zlen=0.5

	call create_iscalar(m%grid,CmSqrt%mask,CELL_EARTH)
	CmSqrt%mask%v = FREE
	CmSqrt%N = 5
    CmSqrt%allocated = .true.

    if (.not. present(cfile)) then
    	! CmSqrt already initialized
    else
		! attempt to read CmSqrt from cfile
	    inquire(FILE=cfile,EXIST=exists)
	    if (exists) then
	      call read_Cm(cfile)
	    else
	      call errStop('Unable to find the input covariance file '//trim(cfile)//' in create_CmSqrt')
	    end if

		if ((CmSqrt%Nx /= m%Nx) .or. (CmSqrt%Ny /= m%Ny) .or. (CmSqrt%NzEarth /= m%NzEarth)) then
	       call errStop('Grid dimensions do not match in input covariance file '//cfile)
	    end if
    end if


    call setup1DCM(m%grid,CmSqrt)



  end subroutine create_CmSqrt
  ! *******************************************************************
 subroutine read_Cm(cfile)

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

    ! Exception rules
    integer, pointer, dimension(:)                   :: mask1, mask2, ii, jj, kk, xyz
    real (kind=prec), pointer, dimension(:)  :: smoothing, S

    integer       			 :: Nx, Ny, NzEarth, nrules, nS, i, j, k, n, istat
    integer                  :: fid = 30

	if (.not. CmSqrt%allocated) then
		call errStop('Model covariance must be allocated before reading from file in read_CmSqrt')
	end if

    open(unit=fid,file=cfile,form='formatted',status='old')

    ! skip the 16 lines header
    do j = 1,16
       read(fid,*)
    end do

	! read grid dimensions
	read(fid,*) Nx,Ny,NzEarth
	CmSqrt%Nx = Nx
	CmSqrt%Ny = Ny
	CmSqrt%NzEarth = NzEarth

	! read smoothing parameters
    read(fid,*) CmSqrt%Sx
    read(fid,*) CmSqrt%Sy
    read(fid,*) CmSqrt%Sz

 CmSqrt%xlen=CmSqrt%Sx(1)
 CmSqrt%ylen=CmSqrt%Sy(1)
 CmSqrt%zlen=CmSqrt%Sz

	! read number of times to apply the smoothing
	read(fid,*) CmSqrt%N

	! read exception rules for smoothing across surfaces
    read(fid,*) nrules
    allocate(mask1(nrules),mask2(nrules),smoothing(nrules),STAT=istat)
    do n = 1,nrules
       read(fid,*) mask1(n), mask2(n), smoothing(n)
    end do

	! create and read the mask array
	call read_iscalar(fid,CmSqrt%mask)

	close(fid)




  end subroutine read_Cm
! *******************************************************************
  subroutine deall_CmSqrt()

    integer                             :: istat

    deallocate(CmSqrt%XDIF, STAT = istat)
    deallocate(CmSqrt%YDIF, STAT = istat)
    deallocate(CmSqrt%ZDIF, STAT = istat)

  end subroutine deall_CmSqrt



!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     SOLVING DIFFUSION EQU. IN EXPLICIT WAY.
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.
!...  Adapted from Setup1DCM_NEW by Anna Kelbert

  SUBROUTINE Setup1DCM(gd,Cm)

     !integer, dimension(:,:,:), intent(in)   :: PMod
     type (grid_t), intent(in)             :: gd
     type (modelCov_t), intent(inout)              :: Cm


     real (kind=prec), pointer, dimension(:,:)   ::  XDF
     real (kind=prec), pointer, dimension(:,:)   ::  YDF
     real (kind=prec), pointer, dimension(:,:)   ::  ZDF
     Integer         , pointer, dimension(:)     ::  xgd,ygd,zgd

     integer                  :: jj, ix, iy, iz, Nx, Ny, NzEarth
     real (kind=prec) :: xlen, ylen, zlen

	  Nx      = gd%Nx
	  Ny      = gd%Ny
	  NzEarth = gd%Nz - gd%nzAir

	  allocate(CmSqrt%XDIF(2,Nx,Ny,NzEarth))
	  allocate(CmSqrt%YDIF(2,Nx,Ny,NzEarth))
	  allocate(CmSqrt%ZDIF(2,Nx,Ny,NzEarth))

	  allocate(XDF(2,Nx))
	  allocate(YDF(2,Ny))
	  allocate(ZDF(2,NzEarth))

	  allocate (xgd(Nx),ygd(Ny),zgd(NzEarth))


	  XLen = CmSqrt%xlen
	  YLen = CmSqrt%ylen
	  ZLen = CmSqrt%zlen

!>>>  X : Horizontal Diffusion
      DO iz = 1,NzEarth
        DO iy = 1,Ny

          DO ix = 1,Nx
            xgd(ix) = CmSqrt%mask%v(ix,iy,iz)
          ENDDO ! ix

          CALL SetupXCM(Nx,xgd,XLen,XDF)
          DO ix = 1,Nx
            CmSqrt%XDIF(1,ix,iy,iz) = XDF(1,ix)
            CmSqrt%XDIF(2,ix,iy,iz) = XDF(2,ix)
          ENDDO ! ix
        ENDDO ! iy
      ENDDO ! iz


!>>>  Y : Horizontal Diffusion
      DO iz = 1,NzEarth
        DO ix = 1,Nx

          DO iy = 1,Ny
            ygd(iy) = CmSqrt%mask%v(ix,iy,iz)
          ENDDO ! iy

          CALL SetupYCM(Ny,ygd,YLen,YDF)
          DO iy = 1,Ny
            CmSqrt%YDIF(1,ix,iy,iz) = YDF(1,iy)
            CmSqrt%YDIF(2,ix,iy,iz) = YDF(2,iy)
          ENDDO ! iy
        ENDDO ! ix
      ENDDO ! iz


!>>>  Vertical Diffusion
      DO ix = 1,Nx
        DO iy = 1,Ny

          DO iz = 1,NzEarth
            zgd(iz) = CmSqrt%mask%v(ix,iy,iz)
          ENDDO ! iz


          CALL SetupZCM(NzEarth,zgd,ZLen,ZDF)
          DO iz = 1,NzEarth
            CmSqrt%ZDIF(1,ix,iy,iz) = ZDF(1,iz)
            CmSqrt%ZDIF(2,ix,iy,iz) = ZDF(2,iz)
          ENDDO ! iy
        ENDDO ! iy
      ENDDO ! ix

     deallocate(XDF,YDF,ZDF)

  END SUBROUTINE Setup1DCM ! Setup1DCM

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!>>>  Horizontal Diffusion : X
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE SetupXCM(Nx,xgd,xlen,XDF)


      Integer,intent(in), dimension(:)     ::  xgd
      integer Nx,info
      real (kind=prec)  XDF(2,Nx),xlen
      real (kind=prec)  cr,cl,difxc(Nx),difxr(Nx)






      integer ix
      real (kind=prec)  D0,D1

      D0 = 0.
      D1 = 1.

      DO ix = 1,Nx
        cr = xlen
        cl = xlen


        IF (ix.LT.Nx) THEN
          IF ((xgd(ix) .eq. OCEAN).AND.(xgd(ix+1).eq. OCEAN)) THEN
            cr = D0
          ENDIF
        ENDIF

        IF (ix.GT.1) THEN
          IF ((xgd(ix) .eq. OCEAN).AND.(xgd(ix-1).eq. OCEAN)) THEN
            cl = D0
          ENDIF
        ENDIF

        difxc(ix) = D1 + cr + cl
        difxr(ix) = -cr
      ENDDO ! ix

      !CALL ConstantMatrixR8(XDF,2,Nx,2,Nx,D0)
      XDF=D0
      DO ix = 2,Nx
        XDF(1,ix) = difxr(ix-1)
      ENDDO
      DO ix = 1,Nx
        XDF(2,ix) = difxc(ix)
      ENDDO

      CALL DPBTRF('U',Nx,1,XDF,2,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : XDiff',info
        WRITE(12,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(12,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : XDiff',info
        STOP
      ENDIF

      RETURN

  END SUBROUTINE SetupXCM ! SetupXCM


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!>>>  Horizontal Diffusion : Y
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE SetupYCM(Ny,ygd,ylen,YDF)

      Integer,intent(in), dimension(:)     ::  ygd
      integer Ny,info
      real (kind=prec)  YDF(2,Ny),ylen
      real (kind=prec)  cr,cl,difyc(Ny),difyr(Ny)
      integer iy
      real (kind=prec)  D0,D1

      D0 = 0.
      D1 = 1.


      DO iy = 1,Ny
        cr = ylen
        cl = ylen

        IF (iy.LT.Ny) THEN
          IF ((ygd(iy) .eq. OCEAN).AND.(ygd(iy+1).eq. OCEAN)) THEN
            cr = D0
          ENDIF
        ENDIF

        IF (iy.GT.1) THEN
          IF ((ygd(iy) .eq. OCEAN).AND.(ygd(iy-1).eq. OCEAN)) THEN
            cl = D0
          ENDIF
        ENDIF

        difyc(iy) = D1 + cr + cl
        difyr(iy) = -cr
      ENDDO ! iy

      !CALL ConstantMatrixR8(YDF,2,Ny,2,Ny,D0)
   YDF=D0
      DO iy = 2,Ny
        YDF(1,iy) = difyr(iy-1)
      ENDDO
      DO iy = 1,Ny
        YDF(2,iy) = difyc(iy)
      ENDDO

      CALL DPBTRF('U',Ny,1,YDF,2,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : YDiff',info
        WRITE(12,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(12,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : YDiff',info
        STOP
      ENDIF

      RETURN

  END SUBROUTINE SetupYCM ! SetupYCM


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!>>>  Vertical Diffusion : Z
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE SetupZCM(Nz,zgd,zlen,ZDF)

      Integer,intent(in), dimension(:)     ::  zgd
      integer Nz,info
      real (kind=prec)  ZDF(2,Nz),zlen
      real (kind=prec)  cr,cl,difzc(Nz),difzr(Nz)
      integer iz
      real (kind=prec)  D0,D1

      D0 = 0.
      D1 = 1.

      DO iz = 1,Nz
        cr = zlen
        cl = zlen


        IF (iz.LT.Nz) THEN
          IF ((zgd(iz) .eq. OCEAN).AND.(zgd(iz+1).eq. OCEAN)) THEN
            cr = D0
          ENDIF
        ENDIF

        IF (iz.GT.1) THEN
          IF ((zgd(iz) .eq. OCEAN).AND.(zgd(iz-1).eq. OCEAN)) THEN
            cl = D0
          ENDIF
        ENDIF


        difzc(iz) = D1 + cr + cl
        difzr(iz) = -cr
      ENDDO ! iz

      !CALL ConstantMatrixR8(ZDF,2,Nz,2,Nz,D0)
ZDF=D0
      DO iz = 2,Nz
        ZDF(1,iz) = difzr(iz-1)
      ENDDO
      DO iz = 1,Nz
        ZDF(2,iz) = difzc(iz)
      ENDDO
!     DO iz = 1,Nz
!       write(6,*) 'ZDF BEFORe ',iz,ZDF(1,iz),ZDF(2,iz)
!     ENDDO

      CALL DPBTRF('U',Nz,1,ZDF,2,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : ZDiff',info
        STOP
      ENDIF

!     DO iz = 1,Nz
!       write(6,*) 'ZDF AFTER ',iz,ZDF(1,iz),ZDF(2,iz)
!     ENDDO
!     stop

      RETURN

  END SUBROUTINE SetupZCM ! SetupZCM


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE SolveDiff(DTime,Cm,u2)

      integer, intent(in)                     :: DTime
      type (modelCov_t), intent(in)               :: Cm
      type (rscalar), intent(inout)           :: u2


      integer  idt,ix,iy,iz
	  integer            :: Nx,Ny,NzEarth

      Nx = u2%nx
      Ny = u2%ny
      NzEarth = u2%nz


      DO idt = 1,DTime
        CALL X_Diff(CmSqrt,u2)
        CALL Y_Diff(CmSqrt,u2)
        CALL Z_Diff(CmSqrt,u2)
      ENDDO ! idt

      DO idt = 1,DTime
        CALL Z_Diff(CmSqrt,u2)
        CALL Y_Diff(CmSqrt,u2)
        CALL X_Diff(CmSqrt,u2)
      ENDDO ! idt

  	! apply the scaling operator C
    do iz = 1,NzEarth
    	do iy = 1,Ny
    		do ix = 1,Nx
				u2%v(ix,iy,iz) = Scaling(ix,iy,iz) * u2%v(ix,iy,iz)
    		end do
    	end do
    end do


  END SUBROUTINE SolveDiff ! SolveDiff


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE X_Diff(Cm,u2)

	  type (modelCov_t), intent(in)      :: Cm
	  type (rscalar), intent(inout)  :: u2


	  real (kind=prec), allocatable, dimension(:,:)   ::  XDF
	  real (kind=prec), allocatable, dimension(:)     ::  ux

	  integer            :: Nx,Ny,NzEarth,ix,iy,iz,info

      Nx = u2%nx
      Ny = u2%ny
      NzEarth = u2%nz
      allocate(XDF(2,Nx))
      allocate(ux(Nx))

      DO iz = 1,NzEarth
        DO iy = 1,Ny
          DO ix = 1,Nx
            ux(ix) = u2%v(ix,iy,iz)
            XDF(1,ix) = CmSqrt%XDIF(1,ix,iy,iz)
            XDF(2,ix) = CmSqrt%XDIF(2,ix,iy,iz)
          ENDDO ! ix

          CALL DPBTRS('U',Nx,1,1,XDF,2,ux,Nx,info)
          IF (info.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
            WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : XDiff',info
            STOP
          ENDIF

          DO ix = 1,Nx
            u2%v(ix,iy,iz) = ux(ix)
          ENDDO
        ENDDO ! iy
      ENDDO ! iz

      deallocate(XDF,ux)

  END SUBROUTINE X_DIFF ! X_DIFF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE Y_Diff(Cm,u2)

	  type (modelCov_t), intent(in)      :: Cm
	  type (rscalar), intent(inout)  :: u2

	  real (kind=prec), allocatable, dimension(:,:)   ::  YDF
	  real (kind=prec), allocatable, dimension(:)     ::  uy

	  integer            :: Nx,Ny,NzEarth,ix,iy,iz,info

      Nx = u2%nx
      Ny = u2%ny
      NzEarth = u2%nz
      allocate(YDF(2,Ny))
      allocate(uy(Ny))

      DO iz = 1,NzEarth
        DO ix = 1,Nx
          DO iy = 1,Ny
            uy(iy) = u2%v(ix,iy,iz)
            YDF(1,iy) = CmSqrt%YDIF(1,ix,iy,iz)
            YDF(2,iy) = CmSqrt%YDIF(2,ix,iy,iz)
          ENDDO ! iy

          CALL DPBTRS('U',Ny,1,1,YDF,2,uy,Ny,info)
          IF (info.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
            WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : YDiff',info
            STOP
          ENDIF

          DO iy = 1,Ny
            u2%v(ix,iy,iz) = uy(iy)
          ENDDO ! iy
        ENDDO ! ix
      ENDDO ! iz

      deallocate(YDF,uy)

  END SUBROUTINE Y_Diff ! Y_DIFF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE Z_Diff(Cm,u2)

	  type (modelCov_t), intent(in)      :: Cm
	  type (rscalar), intent(inout)  :: u2

	  real (kind=prec), allocatable, dimension(:,:)   ::  ZDF
	  real (kind=prec), allocatable, dimension(:)     ::  uz

	  integer            :: Nx,Ny,NzEarth,ix,iy,iz,info


      Nx = u2%nx
      Ny = u2%ny
      NzEarth = u2%nz
      allocate(ZDF(2,NzEarth))
      allocate(uz(NzEarth))

      DO ix = 1,Nx
        DO iy = 1,Ny
          DO iz = 1,NzEarth
            uz(iz) = u2%v(ix,iy,iz)
            ZDF(1,iz) = CmSqrt%ZDIF(1,ix,iy,iz)
            ZDF(2,iz) = CmSqrt%ZDIF(2,ix,iy,iz)
          ENDDO ! iy

          CALL DPBTRS('U',NzEarth,1,1,ZDF,2,uz,NzEarth,info)
          IF (info.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
            WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : ZDiff',info
            STOP
          ENDIF

          DO iz = 1,NzEarth
            u2%v(ix,iy,iz) = uz(iz)
          ENDDO ! iz
        ENDDO ! iy
      ENDDO ! ix







      deallocate(ZDF,uz)

  END SUBROUTINE Z_Diff ! Z_DIFF
! ****************************************************************************
  function Scaling(i,j,k) result(c)

  	! computes the scaling coefficient based on CmSqrt

    implicit none
    integer, intent(in)                      :: i,j,k
    real (kind=prec)                 :: c

    if ((i < 1) .or. (i > CmSqrt%Nx)) then
    	call errStop('index i out of bounds in Scaling(i,j,k)')
    else if ((j < 1) .or. (j > CmSqrt%Ny)) then
    	call errStop('index j out of bounds in Scaling(i,j,k)')
    else if ((k < 1) .or. (k > CmSqrt%NzEarth)) then
    	call errStop('index k out of bounds in Scaling(i,j,k)')
    end if

	if (CmSqrt%mask%v(i,j,k) == AIR) then
		c = 0.0
	else if (CmSqrt%mask%v(i,j,k) == OCEAN) then
		c = 0.0
	else
		c = 1.0
	end if

  end function Scaling
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

!end module modelCov
