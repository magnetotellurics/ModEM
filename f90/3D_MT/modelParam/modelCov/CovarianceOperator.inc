!----------------------------------------------------------------------!
!  Procedures of model covariance operators.                           !
!  Unified interfaces to all operators.                                !
!----------------------------------------------------------------------!


!**********************************************************************!
!   Interface functions.                                               !
!**********************************************************************!
    function multBy_Cm(mhat) result (dm)
!
! Multiplies by the full model covariance,
! which is viewed as a smoothing operator. Intended
! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
! For efficiency, CmSqrt is a saved, private variable inside
! the modelParam module. Before this routine can be called,
! it has to be initialized by calling create_CmSqrt(m).
!
    type (modelParam_t), intent(in)     :: mhat
    type (modelParam_t)                 :: dm
    character(80)                       :: paramType
    Integer                             :: dtime

    if (.not. CovAllocated) then
        call create_CmSqrt(mhat)
    end if

    dm = mhat
    
    select case (flagCov)
    case (1) ! 1D AR
        call RecursiveAR(mhat%cellCond%v, dm%cellCond%v, 2)
    case (2) ! L2(L^TL)
        logcond%v = mhat%cellCond%v
        logcondSmth%v = R_ZERO
        call PCG_Poisson(logcond,logCondSmth)
        logcond%v = R_ZERO
        call PCG_Poisson(logcondSmth,logCond)
        !  can probably just pass dm%cellCond to PCG_Poisson ...
        dm%cellCond%v = logCond%v
    case (3) ! L1
        logcondSmth%v = mhat%cellCond%v
        logcond%v = R_ZERO
        call PCG_Poisson(logcondSmth,logCond)
        dm%cellCond%v = logCond%v
    case default
        write(0,*) 'Unsupported model covariance operator!'
        stop
    end select

    dm%temporary = .true.

    end function multBy_Cm

!----------------------------------------------------------------------!
    function multBy_CmSqrt(mhat) result (dm)
!
! Multiplies by the square root of the model covariance,
! which is viewed as a smoothing operator. Intended
! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
! For efficiency, CmSqrt is a saved, private variable inside
! the modelParam module. Before this routine can be called,
! it has to be initialized by calling create_CmSqrt(m).
!
    type (modelParam_t), intent(in)     :: mhat
    type (modelParam_t)                 :: dm
    character(80)                       :: paramType
    Integer                             :: dtime

    if (.not. CovAllocated) then
        call create_CmSqrt(mhat)
    end if
    dm = mhat
    paramType = LOGE
        call getValue_modelParam(mhat,paramType,logcond)
    select case (flagCov)
    case (1)
        call RecursiveAR(mhat%cellCond%v, dm%cellCond%v, N)
    case (2)
        logcond%v = mhat%cellCond%v
        !call writeRscalarBin(logcond,'x.bin')
        !call writeRscalarBin(logcondSmth,'y.bin')
        logcondSmth%v = R_ZERO
        LaplaceNiter = 0
        LaplaceRerr = R_Zero
        call PCG_Poisson(logcond,logcondSmth)
        ! write(6,*) 'Number of iterations = ', LaplaceNiter
        ! write(6,*) 'relative error = ', LaplaceRerr(LaplaceNiter)
        !  can probably just pass dm%cellCond to PCG_Poisson ...
        !    ... maybe after zeroing
        dm%cellCond%v = logCondSmth%v
    case (3)
        call errStop('cmsqrt is not yet implemented for L1 covariance')
    case default
        write(0,*) 'Unsupported model covariance operator!'
        stop
    end select

    dm%temporary = .true.

    end function multBy_CmSqrt

!----------------------------------------------------------------------!
    function multBy_CmSqrtInv(dm) result (mhat)
!
! Multiplies by the inverse square root of the model covariance,
! which is viewed as a roughening operator. Intended
! to be used to compute \tilde{m} = C_m^{-1/2} ( m - m_0 ).
! For efficiency, CmSqrt is a saved, private variable inside
! the modelParam module. Before this routine can be called,
! it has to be initialized by calling create_CmSqrt(m).
!
    type (modelParam_t), intent(in)     :: dm
    type (modelParam_t)                 :: mhat
    character(80)                       :: paramType
    Integer                             :: dtime

    if (.not. CovAllocated) then
        call create_CmSqrt(dm)
    end if

    mhat = dm
    
    select case (flagCov)
    case (1)
        call RecursiveARInv(dm%cellCond%v, mhat%cellCond%v, N)
    case (2)
        !logcond%v = mhat%cellCond%v
        !call writeRscalarBin(logcond,'x.bin')
        !call writeRscalarBin(logcondSmth,'y.bin')
        !logcondSmth%v = R_ZERO
        mHat%cellcond%v = R_ZERO
        call A(dm%cellCond,mhat%cellCond)
    case (3)
         call errStop('cmsqrtinv is not yet implemented for L1 covariance')
    case default
        write(0,*) 'Unsupported model covariance operator!'
        stop
    end select

    mhat%temporary = .true.

    end function multBy_CmSqrtInv

!----------------------------------------------------------------------!
    function multBy_CmInv(dm) result (mhat)
!
! Multiplies by the inverse square root of the model covariance,
! which is viewed as a roughening operator. Intended
! to be used to compute \tilde{m} = C_m^{-1/2} ( m - m_0 ).
! For efficiency, CmSqrt is a saved, private variable inside
! the modelParam module. Before this routine can be called,
! it has to be initialized by calling create_CmSqrt(m).
!
    type (modelParam_t), intent(in)     :: dm
    type (modelParam_t)                 :: mhat
    character(80)                       :: paramType
    Integer                             :: dtime

    if (.not. CovAllocated) then
        call create_CmSqrt(dm)
    end if

    mhat = dm
    
    select case (flagCov)
    case (1)
        call RecursiveARInv(dm%cellCond%v, mhat%cellCond%v, 2)
    case (2)
        logcond%v = R_ZERO
        call A(dm%cellCond,logcond)
        mHat%cellcond%v = R_ZERO
        call A(logcond,mHat%cellCond)
    case (3)
        mHat%cellcond%v = R_ZERO
        call A(dm%cellCond,mHat%cellCond)
    case default
        write(0,*) 'Unsupported model covariance operator!'
        stop
    end select

    mhat%temporary = .true.

    end function multBy_CmInv
!----------------------------------------------------------------------!
    subroutine create_CmSqrt(m,covType,cfile)
!
! Initializes CmSqrt variable stored in CovarianceOperator.hd. If cfile
! is specified, gets this information from file.
!
    type (modelParam_t), intent(in)     :: m
    integer, intent(in), optional       :: covType
    character(*), intent(in), optional  :: cfile
    integer                             :: istat
    logical                             :: exists
! initializing CmSqrt. ... 
! I am changing all variables to be module veriables, eliminating the
! covariance type ... 
! this is just used to store all these variables used in the
! model_space module for covariance operations
!
    Nx = m%Nx
    Ny = m%Ny
    NzEarth = m%NzEarth
    flagCov = covType

    allocate(Sx(NzEarth), STAT=istat)
    allocate(Sy(NzEarth), STAT=istat)
    Sx(:) = 0.3
    Sy(:) = 0.3
    select case (flagCov)
    case (1) ! 1D AR
        allocate(Sz(1), STAT=istat)
        Sz(:) = 0.3
        N  = 1
    case (2) ! L2
        allocate(Sz(NzEarth), STAT=istat)
        Sz(:) = 0.3
        !   this sets module variable modelGrid
        call setLaplaceGrid(m%grid)
        !    this sest Laplace operator
        call initializeLaplace()
        !call writeRvectorBin(Lm,'Lm.bin')
        !call writeRvectorBin(Lp,'Lp.bin')
        !call writeRscalarBin(Lc,'Lc.bin')
        !     initialize rscalar arrays for smoothing wih Laplacian
        call create_rscalar(modelGrid,logcond,CELL_EARTH)
        call create_rscalar(modelGrid,logcondSmth,CELL_EARTH)
    case (3) ! L1, do the same as in L2
        allocate(Sz(NzEarth), STAT=istat)
        Sz(:) = 0.3
        !   this sets module variable modelGrid
        call setLaplaceGrid(m%grid)
        !   this sest Laplace operator
        call initializeLaplace()
        !   initialize rscalar arrays for smoothing wih Laplacian
        call create_rscalar(modelGrid,logcond,CELL_EARTH)
        call create_rscalar(modelGrid,logcondSmth,CELL_EARTH)
    case default
        call errStop('Unsupported covariance operator!')
    end select

    call create_iscalar(m%grid,mask,CELL_EARTH)
    mask%v = FREE
    CovAllocated = .true.
!
! read file if needed.
!
    if (.not. present(cfile)) then
        ! CmSqrt already initialized
    else
        ! attempt to read CmSqrt from cfile
        inquire(FILE=cfile,EXIST=exists)
        if (exists) then
            call read_Cm(cfile,flagCov)
        else
            call errStop('Unable to find the input covariance file '//trim(cfile)//' in create_CmSqrt')
        end if

        if ((Nx /= m%Nx) .or. (Ny /= m%Ny) .or. (NzEarth /= m%NzEarth)) then
            call errStop('Grid dimensions do not match in input covariance file '//cfile)
        end if
    end if

    end subroutine create_CmSqrt
!----------------------------------------------------------------------!
    subroutine deall_CmSqrt()

    integer :: istat

    deallocate(Sx,Sy,Sz, STAT=istat)
    call deall_iscalar(mask)
    select case (flagCov)
    case (1)
        call deall_sparsevecc(S)
    case(2)
        call deall_Laplace()
        call deall_rscalar(logcond)
        call deall_rscalar(logcondSmth)
    case (3)
        call deall_Laplace()
        call deall_rscalar(logcond)
        call deall_rscalar(logcondSmth)
    case default
    end select

    CovAllocated = .false.

    end subroutine deall_CmSqrt
!----------------------------------------------------------------------!
    subroutine read_Cm(cfile, CovType)
    implicit none
    character(*), intent(in)             :: cfile
    integer, intent(in),optional         :: CovType

    ! local variables 
    integer                              :: flagCov
    logical                              :: exists

    if (.not.present(CovType)) then
        flagCov = 1
    else
        flagCov = CovType
    endif

    inquire(FILE=cfile,EXIST=exists)
    if (exists) then
        select case (flagCov)
        case (1)
            call read_Cm_RecursiveAR(cfile)
        case (2)
            call read_Cm_Laplacian(cfile)
        case (3)
            call read_Cm_Laplacian(cfile)
        case default
            call errStop('Unsupported covariance operator!')
        end select
    else
        call errStop('Unable to find the input covariance file '//trim(cfile)//' in create_CmSqrt')
    end if

    end subroutine read_Cm
!----------------------------------------------------------------------!

!**********************************************************************!
!   Implementations for RecursiveAR operator.                          !
!**********************************************************************!
!----------------------------------------------------------------------!
!...  Copyright (C) 2008 Anna Kelbert. All rights reserved.
!----------------------------------------------------------------------!
!+-----------------------------------------------------------------------------+
!| This file defines model covariance for a recursive autoregression scheme.   |
!| The model space may be divided into distinct areas using integer masks.     |
!| Mask 0 is reserved for air; mask 1 is reserved for ocean. Smoothing between |
!| air, ocean and the rest of the model is turned off automatically. You can   |
!| also define exceptions to override smoothing between any two model areas.   |
!| To turn off smoothing set it to zero. This header is 16 lines long.         |
!| 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |
!| 2. Smoothing in the X direction (NzEarth real values)                       |
!| 3. Smoothing in the Y direction (NzEarth real values)                       |
!| 4. Vertical smoothing (1 real value)                                        |
!| 5. Number of times the smoothing should be applied (1 integer >= 0)         |
!| 6. Number of exceptions (1 integer >= 0)                                    |
!| 7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 3 & 4) |
!| 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|
!+-----------------------------------------------------------------------------+

    subroutine read_Cm_RecursiveAR(cfile)
!
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
!
    character(*), intent(in)             :: cfile
!
! Exception rules
!
    integer,pointer,dimension(:)         :: mask1, mask2, ii, jj, kk, xyz
    real(kind=prec),pointer,dimension(:) :: smoothing

    integer :: nrules, nS, i, j, k, n, istat, NxFile, NyFile, NzFile
    integer :: fid = 30

    if (.not. CovAllocated) then
        call errStop('Model covariance must be allocated before reading from file in read_CmSqrt')
    end if

    open(unit=fid,file=cfile,form='formatted',status='old')

    ! skip the 16 lines header
    do j = 1,16
        read(fid,*)
    end do

    ! read grid dimensions: check consistency with grid
    read(fid,*) NxFile,NyFile,NzFile    
    if ((Nx.ne.NxFile) .or. (Ny.ne.NyFile).or.(NzEarth.ne.NzFile)) then
            call errStop('Grid dimensions do not match in input covariance file '//cfile)
        end if
!!   this makes no sense  ... we allocated from model parameter, and we check consistency with model
!            parameter AFTER reading Sx etc.  But if things are not consistent, we will have an error before
!!       we check!!

    ! read smoothing parameters
    read(fid,*) Sx
    read(fid,*) Sy
    read(fid,*) Sz

    ! read number of times to apply the smoothing
    read(fid,*) N

    ! read exception rules for smoothing across surfaces
    read(fid,*) nrules
    allocate(mask1(nrules),mask2(nrules),smoothing(nrules),STAT=istat)
    do n = 1,nrules
        read(fid,*) mask1(n), mask2(n), smoothing(n)
    end do

    ! create and read the mask array
    call read_iscalar(fid,mask)

    close(fid)

    ! create a huge sparse vector to make sure we accommodate all smoothing exceptions
    call create_sparsevecc(Nx*Ny*NzEarth, S, FACE)

    ! now, parse the exceptions
    nS = 0
    do k = 2,NzEarth
        do j = 2,Ny
            do i = 2,Nx
                do n = 1,nrules
                    ! look back in the X-direction
                    if (((mask%v(i-1,j,k) == mask1(n)) .and. (mask%v(i,j,k) == mask2(n))) &
                        .or. ((mask%v(i-1,j,k) == mask2(n)) .and. (mask%v(i,j,k) == mask1(n)))) &
                        then
                        nS = nS+1
                        S%i(nS) = i-1
                        S%j(nS) = j
                        S%k(nS) = k
                        S%xyz(nS) = 1
                        S%c(nS) = smoothing(n)
                    end if
                    ! look back in the Y-direction
                    if (((mask%v(i,j-1,k) == mask1(n)) .and. (mask%v(i,j,k) == mask2(n))) &
                        .or. ((mask%v(i,j-1,k) == mask2(n)) .and. (mask%v(i,j,k) == mask1(n)))) &
                        then
                        nS = nS+1
                        S%i(nS) = i
                        S%j(nS) = j-1
                        S%k(nS) = k
                        S%xyz(nS) = 2
                        S%c(nS) = smoothing(n)
                    end if
                    ! look back in the Z-direction
                    if (((mask%v(i,j,k-1) == mask1(n)) .and. (mask%v(i,j,k) == mask2(n))) &
                        .or. ((mask%v(i,j,k-1) == mask2(n)) .and. (mask%v(i,j,k) == mask1(n)))) &
                        then
                        nS = nS+1
                        S%i(nS) = i
                        S%j(nS) = j
                        S%k(nS) = k-1
                        S%xyz(nS) = 3
                        S%c(nS) = smoothing(n)
                    end if
                end do
            end do
        end do
    end do
    deallocate(mask1,mask2,smoothing)

    ! now, truncate the smoothing vector to the correct number of components
    call reall_sparsevecc(nS, S)

    end subroutine read_Cm_RecursiveAR

!----------------------------------------------------------------------!
  subroutine RecursiveAR(w,v,n)
!
! Implements the recursive autoregression algorithm for a 3D real array.
! In our case, the assumed-shape array would be e.g. conductivity
! in each cell of the Nx x Ny x NzEarth grid.
!
    real (kind=prec), intent(in)     :: w(:,:,:)
    real (kind=prec), intent(out)    :: v(:,:,:)
    integer, intent(in)              :: n
    integer                          :: Nx, Ny, NzEarth, i, j, k, iSmooth

    Nx      = size(w,1)
    Ny      = size(w,2)
    NzEarth = size(w,3)

    if (maxval(abs(shape(w) - shape(v)))>0) then
        call errStop('The input arrays should be of the same shapes in RecursiveAR')
    end if

    v = w

    do iSmooth = 1,n

        ! smooth in the X-direction (Sx)
        do k = 1,NzEarth
            do j = 1,Ny
                !v(1,j,k) = v(1,j,k)
                do i = 2,Nx
                    v(i,j,k) = SmoothX(i-1,j,k) * v(i-1,j,k) + v(i,j,k)
                end do
            end do
        end do

        ! smooth in the Y-direction (Sy)
        do k = 1,NzEarth
            do i = 1,Nx
                ! v(i,1,k) = v(i,1,k)
                do j = 2,Ny
                    v(i,j,k) = SmoothY(i,j-1,k) * v(i,j-1,k) + v(i,j,k)
                end do
            end do
        end do

        ! smooth in the Z-direction (Sz)
        do j = 1,Ny
            do i = 1,Nx
                ! v(i,j,1) = v(i,j,1)
                do k = 2,NzEarth
                    v(i,j,k) = SmoothZ(i,j,k-1) * v(i,j,k-1) + v(i,j,k)
                end do
            end do
        end do

        ! smooth in the Z-direction (Sz^T)
        do j = Ny,1,-1
            do i = Nx,1,-1
                ! v(i,j,NzEarth) = v(i,j,NzEarth)
                do k = NzEarth,2,-1
                    v(i,j,k-1) = v(i,j,k-1) + SmoothZ(i,j,k-1) * v(i,j,k)
                end do
            end do
        end do

        ! smooth in the Y-direction (Sy^T)
        do k = NzEarth,1,-1
            do i = Nx,1,-1
                ! v(i,Ny,k) = v(i,Ny,k)
                do j = Ny,2,-1
                    v(i,j-1,k) = v(i,j-1,k) + SmoothY(i,j-1,k) * v(i,j,k)
                end do
            end do
        end do

        ! smooth in the X-direction (Sx^T)
        do k = NzEarth,1,-1
            do j = Ny,1,-1
                ! v(Nx,j,k) = v(Nx,j,k)
                do i = Nx,2,-1
                    v(i-1,j,k) = v(i-1,j,k) + SmoothX(i-1,j,k) * v(i,j,k)
                end do
            end do
        end do

    end do

    ! apply the scaling operator C
    do k = 1,NzEarth
        do j = 1,Ny
            do i = 1,Nx
                v(i,j,k) = (Scaling(i,j,k)**n) * v(i,j,k)
            end do
        end do
    end do

    end subroutine RecursiveAR

!----------------------------------------------------------------------!
    subroutine RecursiveARInv(w,v,n)
!
! ... and the inverse "roughening" operator useful for starting
! the inversion with an arbitrary model: \tilde{m} = C_m^{-1/2} (m - m_0).
! In our case, the assumed-shape array would be e.g. conductivity
! in each cell of the Nx x Ny x NzEarth grid.
! NOTE: the inverse covariance operator is poorly conditioned!!!
! e.g., at alpha=0.3 n=4 white noise completely overwhelmes the
! inverse model. Be extra careful when you use this function and
! always look at the result before using it to start the inversion.
! In the future, may want to stabilize this.
!
    real (kind=prec), intent(in)     :: w(:,:,:)
    real (kind=prec), intent(out)    :: v(:,:,:)
    real (kind=prec), allocatable    :: u(:,:,:)
    integer, intent(in)              :: n
    integer                          :: Nx, Ny, NzEarth, i, j, k, iSmooth, istat

    Nx      = size(w,1)
    Ny      = size(w,2)
    NzEarth = size(w,3)

    if (maxval(abs(shape(w) - shape(v)))>0) then
        call errStop('The input arrays should be of the same shapes in RecursiveARInv')
    end if

    allocate(u(Nx,Ny,NzEarth),stat=istat)
    v = w

    do iSmooth = 1,n

        u = v

        ! invert smoothing in the X-direction (Sx^T)
        do k = NzEarth,1,-1
            do j = Ny,1,-1
                v(Nx,j,k) = u(Nx,j,k)
                do i = Nx,2,-1
                    v(i-1,j,k) = u(i-1,j,k) - SmoothX(i-1,j,k) * u(i,j,k)
                end do
            end do
        end do
        
        u = v
        
        ! invert smoothing in the Y-direction (Sy^T)
        do k = NzEarth,1,-1
            do i = Nx,1,-1
                v(i,Ny,k) = u(i,Ny,k)
                do j = Ny,2,-1
                    v(i,j-1,k) = u(i,j-1,k) - SmoothY(i,j-1,k) * u(i,j,k)
                end do
            end do
        end do
        
        u = v
        
        ! invert smoothing in the Z-direction (Sz^T)
        do j = Ny,1,-1
            do i = Nx,1,-1
                v(i,j,NzEarth) = u(i,j,NzEarth)
                do k = NzEarth,2,-1
                    v(i,j,k-1) = u(i,j,k-1) - SmoothZ(i,j,k-1) * u(i,j,k)
                end do
            end do
        end do
        
        u = v
        
        ! invert smoothing in the Z-direction (Sz)
        do j = 1,Ny
            do i = 1,Nx
                v(i,j,1) = u(i,j,1)
                do k = 2,NzEarth
                    v(i,j,k) = - SmoothZ(i,j,k-1) * u(i,j,k-1) +u(i,j,k)
                end do
            end do
        end do
        
        u = v
        
        ! invert smoothing in the Y-direction (Sy)
        do k = 1,NzEarth
            do i = 1,Nx
                v(i,1,k) = u(i,1,k)
                do j = 2,Ny
                    v(i,j,k) = - SmoothY(i,j-1,k) * u(i,j-1,k) + u(i,j,k)
                end do
            end do
        end do
        
        u = v
        
        ! invert smoothing in the X-direction (Sx)
        do k = 1,NzEarth
            do j = 1,Ny
                v(1,j,k) = u(1,j,k)
                do i = 2,Nx
                    v(i,j,k) = - SmoothX(i-1,j,k) * u(i-1,j,k) + u(i,j,k)
                end do
            end do
        end do
    
    end do
    
    ! apply the inverse of the scaling operator C
    do k = 1,NzEarth
        do j = 1,Ny
            do i = 1,Nx
                if (abs(Scaling(i,j,k)) < R_TINY) then
                    v(i,j,k) = 0
                else
                    v(i,j,k) =  v(i,j,k) / (Scaling(i,j,k)**n)
                end if
            end do
        end do
    end do
    
    deallocate(u,stat=istat)

    end subroutine RecursiveARInv

!----------------------------------------------------------------------!
    function SmoothX(i,j,k) result(alpha)
!
! computes the smoothing coefficient in the x-direction based on CmSqrt
!
    implicit none
    integer, intent(in)              :: i,j,k
    real (kind=prec)                 :: alpha
    
    ! local variables
    integer                          :: n
    
    if (.not. associated(Sx)) then
        call errStop('Sx has to be allocated before calling SmoothX')
    end if
    
    alpha = Sx(k)
    
    if ((i < 1) .or. (i > Nx)) then
        call errStop('index i out of bounds in SmoothX(i,j,k)')
    else if ((j < 1) .or. (j > Ny)) then
        call errStop('index j out of bounds in SmoothX(i,j,k)')
    else if ((k < 1) .or. (k > NzEarth)) then
        call errStop('index k out of bounds in SmoothX(i,j,k)')
    end if
    
    if (S%allocated .and. (S%nCoeff > 0)) then
        ! scan through the special rules and possibly update the result
        do n = 1,S%nCoeff
            if (S%xyz(n) == 1) then
                if ((S%i(n) == i) .and. (S%j(n) == j) .and. (S%k(n) == k)) then
                    alpha = S%c(n)
                    exit
                end if
            end if
        end do
    end if

    end function SmoothX

!----------------------------------------------------------------------!
    function SmoothY(i,j,k) result(beta)
!
! computes the smoothing coefficient in the y-direction based on CmSqrt
!
    implicit none
    integer, intent(in)              :: i,j,k
    real (kind=prec)                 :: beta
    
    ! local variables
    integer                          :: n
    
    if (.not. associated(Sy)) then
        call errStop('Sy has to be allocated before calling SmoothY')
    end if
    
    beta = Sy(k)
    
    if ((i < 1) .or. (i > Nx)) then
        call errStop('index i out of bounds in SmoothY(i,j,k)')
    else if ((j < 1) .or. (j > Ny)) then
        call errStop('index j out of bounds in SmoothY(i,j,k)')
    else if ((k < 1) .or. (k > NzEarth)) then
        call errStop('index k out of bounds in SmoothY(i,j,k)')
    end if
    
    if (S%allocated .and. (S%nCoeff > 0)) then
        ! scan through the special rules and possibly update the result
        do n = 1,S%nCoeff
            if (S%xyz(n) == 2) then
                if ((S%i(n) == i) .and. (S%j(n) == j) .and. (S%k(n) == k)) then
                    beta = S%c(n)
                    exit
                end if
            end if
        end do
    end if

    end function SmoothY

!----------------------------------------------------------------------!
    function SmoothZ(i,j,k) result(gamma)
!
! computes the smoothing coefficient in the z-direction based on CmSqrt
!
    implicit none
    integer, intent(in)              :: i,j,k
    real (kind=prec)                 :: gamma
    
    ! local variables
    integer                          :: n
    
    gamma = Sz(1)
    
    if ((i < 1) .or. (i > Nx)) then
        call errStop('index i out of bounds in SmoothZ(i,j,k)')
    else if ((j < 1) .or. (j >Ny)) then
        call errStop('index j out of bounds in SmoothZ(i,j,k)')
    else if ((k < 1) .or. (k > NzEarth)) then
        call errStop('index k out of bounds in SmoothZ(i,j,k)')
    end if
    
    if (S%allocated .and. (S%nCoeff > 0)) then
        ! scan through the special rules and possibly update the result
        do n = 1,S%nCoeff
            if (S%xyz(n) == 3) then
                if ((S%i(n) == i) .and. (S%j(n) == j) .and. (S%k(n) == k)) then
                    gamma = S%c(n)
                    exit
                end if
            end if
        end do
    end if

    end function SmoothZ

!----------------------------------------------------------------------!
    function Scaling(i,j,k) result(c)
!
! computes the scaling coefficient based on CmSqrt
!
    implicit none
    integer, intent(in)              :: i,j,k
    real (kind=prec)                 :: c, alpha, beta, gamma
    
    if ((i < 1) .or. (i > Nx)) then
        call errStop('index i out of bounds in Scaling(i,j,k)')
    else if ((j < 1) .or. (j > Ny)) then
        call errStop('index j out of bounds in Scaling(i,j,k)')
    else if ((k < 1) .or. (k > NzEarth)) then
        call errStop('index k out of bounds in Scaling(i,j,k)')
    end if
    
    alpha = SmoothX(i,j,k)
    beta  = SmoothY(i,j,k)
    gamma = SmoothZ(i,j,k)

    if (mask%v(i,j,k) == AIR) then
        c = 0.0
    else if (mask%v(i,j,k) == OCEAN) then
        c = 0.0
    else
        c = (1 - alpha)**2 * (1 - beta)**2 * (1 - gamma)**2
    end if

    end function Scaling

!----------------------------------------------------------------------!

!**********************************************************************!
!   Implementations for Diffusion operator.                            !
!**********************************************************************!
!----------------------------------------------------------------------!

!    this is eliminated for now -- not really used for a long time, not
!             clear it is worth supporing at all

!----------------------------------------------------------------------!

!**********************************************************************!
!   Implementations for Laplacian operator.                            !
!**********************************************************************!

 subroutine setLaplaceCov(gridType,constantScale, verticalScale, horizontalScale)
   integer, intent(in)   :: gridType
   real (kind=prec)      :: verticalScale, constantScale, horizontalScale

   !   just sets some module variables for Laplacian smoothing
   laplaceGridType = gridType
   laplaceConstantScale = constantScale
   laplaceVerticalScale = verticalScale
   laplaceHorizontalScale = horizontalScale
   ! write(6,*) 'vertical scale =', laplaceVerticalScale
   ! write(6,*) 'horizontal scale =', laplaceHorizontalScale
   ! write(6,*) 'constant scale =', laplaceConstantScale
   end subroutine
!----------------------------------------------------------------------!

  !**********************************************************************
  ! subroutine to setup the value of the laplacian operator
  subroutine setLaplaceGrid(fullGrid)
  ! This routine makes a grid with no air layers; a uniform grid (with
  ! all cell dimensions identical in all directions) is optional,
  ! controlled by "LaplaceGridType"
  ! I am now making this part of the Laplace covariance; perhaps should
  ! be a utility in one of the grid modules -- would need to change some
  ! things then (now using module variables)
  ! create grid (no air layers) for model parameter
    implicit none
    type(grid_t), intent(in)    :: fullGrid
    real (kind=prec)            :: origin(3)
    integer                     :: Nz, NzSurface, xCum, yCum, zCum
    integer                     :: ix, iy, iz
    real (kind=prec)            :: ax, ay, az

  !
    call create_grid(Nx,Ny,0,NzEarth,modelGrid)
    Nz = NzEarth
    select case (LaplaceGridType)
        case (UNIFORM_GRID)
            modelGrid%dx(1:Nx) = 1.0
            modelGrid%dy(1:Ny) = 1.0
            modelGrid%dz(1:Nz) = 1.0
        case (MODEL_GRID)
            nzSurface =  fullGrid%nzAir+1
            ax = fullGrid%dx(1)
            ay = fullGrid%dy(1)
            az = fullGrid%dz(1)
            ! write(6,*) 'ax, ay, az =', ax, ay, az
            ! write(6,*) 'Nx, Ny, Nz =', Nx, Ny, NzEarth
            modelGrid%dx(1:Nx) = fullGrid%dx(1:Nx)/ax
            modelGrid%dy(1:Ny) = fullGrid%dy(1:Ny)/ay
            modelGrid%dz(1:Nz) = fullGrid%dz(1:NzSurface)/az
    end select
    origin(1) = fullGrid%ox
    origin(2) = fullGrid%oy
    origin(3) = fullGrid%oz
    ! call setup_grid(modelGrid,origin)
    ! NOTE: here we do the grid setup manually, as the setup_grid 
    ! subroutine will only work with an airlayer
    ! 
    modelGrid%dxinv = 1/ modelGrid%dx
    modelGrid%dyinv = 1/ modelGrid%dy
    modelGrid%dzinv = 1/ modelGrid%dz

    modelGrid%rotdeg = modelGrid%rotdeg
    modelGrid%ox = origin(1)
    modelGrid%oy = origin(2)
    modelGrid%oz = origin(3)
    modelGrid%xEdge(1) = modelGrid%ox
    modelGrid%yEdge(1) = modelGrid%oy
    modelGrid%zEdge(1) = modelGrid%oz
    xCum = 0.0
    yCum = 0.0
    zCum = 0.0
    do ix = 1, modelGrid%nx
       xCum = xCum + modelGrid%dx(ix)
       modelGrid%xEdge(ix+1) = xCum + modelGrid%ox
    enddo
    do iy = 1, modelGrid%ny
       yCum = yCum + modelGrid%dy(iy)
       modelGrid%yEdge(iy+1) = yCum + modelGrid%oy
    enddo
    !  NOTE: adjust for origin later to get airthickness, reference to origin
    !    at Earth's surface correct!
    do iz = 1, modelGrid%nz
       zCum = zCum + modelGrid%dz(iz)
       modelGrid%zEdge(iz+1) = zCum
    enddo
    modelGrid%zAirThick = 0.0

    ! distance between center of the grids
    modelGrid%delX(1) = modelGrid%dx(1)
    DO ix = 2,modelGrid%nx
       modelGrid%delX(ix) = modelGrid%dx(ix-1) + modelGrid%dx(ix)
    ENDDO
    modelGrid%delX(modelGrid%nx+1) = modelGrid%dx(modelGrid%nx)
    modelGrid%delX = modelGrid%delX/2.0

    modelGrid%delY(1)    = modelGrid%dy(1)
    DO iy = 2,modelGrid%ny
       modelGrid%delY(iy) = modelGrid%dy(iy-1) + modelGrid%dy(iy)
    ENDDO
    modelGrid%delY(modelGrid%ny+1) = modelGrid%dy(modelGrid%ny)
    modelGrid%delY = modelGrid%delY/2.0

    modelGrid%delZ(1)    = modelGrid%dz(1)
    DO iz = 2,modelGrid%nz
       modelGrid%delZ(iz) = modelGrid%dz(iz-1) + modelGrid%dz(iz)
    ENDDO
    modelGrid%delZ(modelGrid%nz+1) = modelGrid%dz(modelGrid%nz)
    modelGrid%delZ = modelGrid%delZ/ 2.0

    modelGrid%delXinv = 1/ modelGrid%delX
    modelGrid%delYinv = 1/ modelGrid%delY
    modelGrid%delZinv = 1/ modelGrid%delZ

    xCum = 0.0
    yCum = 0.0
    zCum = 0.0
    ! cumulative distance between the centers, adjusted to model origin
    do ix = 1, modelGrid%nx
       xCum = xCum + modelGrid%delX(ix)
       modelGrid%xCenter(ix) = xCum + modelGrid%ox
    enddo
    do iy = 1, modelGrid%ny
       yCum = yCum + modelGrid%delY(iy)
       modelGrid%yCenter(iy) = yCum + modelGrid%oy
    enddo
    do iz = 1, modelGrid%nz
       zCum = zCum + modelGrid%delZ(iz)
       modelGrid%zCenter(iz) = zCum
    enddo

    !  need to be careful here ... grid origin is given at Earth's surface,
    !   not top of model domain!
    do iz = 1, modelGrid%nz
       modelGrid%zCenter(iz) = modelGrid%zCenter(iz)-modelGrid%zAirThick+modelGrid%oz
       modelGrid%zEdge(iz) = modelGrid%zEdge(iz)-modelGrid%zAirThick+modelGrid%oz
    enddo
    modelGrid%zEdge(modelGrid%nz+1) = modelGrid%zEdge(modelGrid%nz+1)-modelGrid%zAirThick+modelGrid%oz
    ! finished setting up the grid here...

  end subroutine !   setLaplaceGrid

  subroutine initializeLaplace()
!   again, this is written to be part of the module where the covariance
!   implementation sits  ... modelGrid is thus a saved module variable
    
       implicit none
       integer                           :: ix,iy,iz, status
       type (rscalar)                    :: Lt

       real (kind=prec), pointer, dimension(:)        :: dx, dy, dz
       real (kind=prec)                               :: s, s1

       Call create_rvector(modelGrid, Lm, FACE)
       Call create_rvector(modelGrid, Lp, FACE)
       Call create_rscalar(modelGrid, Lc, CENTER)
       ! might need these if a preconditoner is created
       ! Call create_rscalar(modelGrid, Ld, CENTER)
       ! Call create_rscalar(modelGrid, Lt, CORNER)
       
       Call zero(Lm)
       Call zero(Lp)
       do iz = 1,NzEarth
           do iy = 1,Ny 
              do ix = 1,Nx-1 
                  Lm%x(ix+1,iy,iz) = MinusOne*laplaceHorizontalScale/   &
                      (modelGrid%dx(ix+1)*modelGrid%delX(ix+1))
                  Lp%x(ix+1,iy,iz) = MinusOne*laplaceHorizontalScale/   &
                      (modelGrid%dx(ix)*modelGrid%delX(ix+1))
             enddo 
          enddo 
       enddo 
       do iz = 1,NzEarth
           do iy = 1,Ny-1
              do ix = 1,Nx 
                  Lm%y(ix,iy+1,iz) = MinusOne*laplaceHorizontalScale/   &
                      (modelGrid%dy(iy+1)*modelGrid%delY(iy+1))
                  Lp%y(ix,iy+1,iz) = MinusOne*laplaceHorizontalScale/   &
                      (modelGrid%dy(iy)*modelGrid%delY(iy+1))
             enddo 
          enddo 
       enddo 
       do iz = 1,NzEarth-1
           do iy = 1,Ny
              do ix = 1,Nx 
                  Lm%z(ix,iy,iz+1) = MinusOne*laplaceVerticalScale/   &
                      (modelGrid%dz(iz+1)*modelGrid%delZ(iz+1))
                  Lp%z(ix,iy,iz+1) = MinusOne*laplaceVerticalScale/   &
                      (modelGrid%dz(iz)*modelGrid%delZ(iz+1))
             enddo 
         enddo 
      enddo 
      
      Lc%v(1:Nx,1:Ny,1:NzEarth) = -(Lm%x(1:Nx,1:Ny,1:NzEarth) &
                           +Lp%x(2:Nx+1,1:Ny,1:NzEarth) &
                           +Lm%y(1:Nx,1:Ny,1:NzEarth)   &
                           +Lp%y(1:Nx,2:Ny+1,1:NzEarth) &
                           +Lm%z(1:Nx,1:Ny,1:NzEarth)   &
                           +Lp%z(1:Nx,1:Ny,2:NzEarth+1))
      ! try to normalize the laplacian by the mean of Lc
      s = R_ZERO
      s1 = R_ZERO
      do iz = 1,NzEarth
          do iy = 1,Ny
             do ix = 1,Nx 
                 s = s + Lc%v(ix,iy,iz) 
                 s1 = s1 + ONE
             enddo
          enddo
      enddo
      s = s/s1/6.0
      ! write(6,*) 'scaling for laplacian is', s 

      ! now scale the Laplacian
      Lc%v(1:Nx,1:Ny,1:NzEarth) = Lc%v(1:Nx,1:Ny,1:NzEarth)/s
      Lm%x(1:Nx+1,1:Ny,1:NzEarth) = Lm%x(1:Nx+1,1:Ny,1:NzEarth)/s
      Lm%y(1:Nx,1:Ny+1,1:NzEarth) = Lm%y(1:Nx,1:Ny+1,1:NzEarth)/s
      Lm%z(1:Nx,1:Ny,1:NzEarth+1) = Lm%z(1:Nx,1:Ny,1:NzEarth+1)/s
      Lp%x(1:Nx+1,1:Ny,1:NzEarth) = Lp%x(1:Nx+1,1:Ny,1:NzEarth)/s
      Lp%y(1:Nx,1:Ny+1,1:NzEarth) = Lp%y(1:Nx,1:Ny+1,1:NzEarth)/s
      Lp%z(1:Nx,1:Ny,1:NzEarth+1) = Lp%z(1:Nx,1:Ny,1:NzEarth+1)/s

    ! TODO: setup the preconditioner - if I can understand the matrix free 
    ! convention some day!
    ! first set the top nodes to 1.0
    ! Ld%v(1,:,:) = 1.0
    ! Ld%v(:,1,:) = 1.0
    ! Ld%v(:,:,1) = 1.0
    ! 
    ! 
    ! compute the inverse diagonal elements for D-ILU (for all nodes)
    ! i suppose this could be simplified - 
    ! Call ZERO(Lt)
    ! do iz = 2, NzEarth+1
    !    do iy = 2, Ny+1
    !       do ix = 2, Nx+1

           !  Ld(i,j,k) should only be modified by Ld from the immediate
           !  previous adjacent cells, i.e. Ld(i-1,j,k), Ld(i,j-1,k) and 
           !  Ld(i,j,k-1), these should all be zero for the first cell (1,1,1)
    !          Lt%v(ix, iy, iz) = Lc%v(ix-1, iy-1, iz-1) - &
    !               Lm%x(ix-1,iy-1,iz-1)*Lp%x(ix,iy-1,iz-1)*Lt%v(ix-1,iy,iz)-&
    !               Lm%y(ix-1,iy-1,iz-1)*Lp%y(ix-1,iy,iz-1)*Lt%v(ix,iy-1,iz)-&
    !               Lm%z(ix-1,iy-1,iz-1)*Lp%z(ix-1,iy-1,iz)*Lt%v(ix,iy,iz-1)
    !          Lt%v(ix, iy, iz) = 1.0/ Lt%v(ix, iy, iz)

    !       enddo
    !    enddo
    ! enddo

    ! Ld%v(1:Nx,1:Ny,1:NzEarth) = Lt%v(2:Nx+1,2:Ny+1,2:NzEarth+1)
    !   allocate for PCG rerr
    allocate(LaplaceRerr(LaplaceMaxIt), STAT=status)
    LaplaceRerr = 0.0
    Call deall_rscalar(Lt)
   end subroutine  ! intializeLaplace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_Cm_Laplacian(cfile)
!
! this mimics the read_Cm_AR subroutine that reads the cov parameters
! the actual information is read from a .cov file. 
! Also, we read an integer mask
! array that subdivides the model grid into different regions
! (AIR, OCEAN, EARTH) and a set of rules that overrides the default
! smoothing parameters across a particular surface between two
! distinct regions. We use this to set up the covariance CmSqrt.
!
! 
! I don't see the use of setting different number of smoothings with L1/L2,
! so - 
! essentially the mask idea is like turning on/off the smoothing at boundaries 
! between different numbered regions -
! 
!      *-------*-------*-------*
!      |       |       |       |
!      |   2   X   1   X   2   |
!      |       |       |       |
!      *---O---*---O---*---X---*
!      |       |       |       |
!      |   2   X   1   O   1   |
!      |       |       |       |
!      *---O---*---X---*---O---*
!      |       |       |       |
!      |   2   O   2   X   1   |
!      |       |       |       |
!      *-------*-------*-------*
!
! e.g. the boundary smoothing between the above two regions 1 and 2 are 
! turned off (X), while the smoothing for inner boundaries between cells are
! left on (O)
! 
! in practice, this is done by multiplying the corresponding stencil by zero 
! or one, which becomes awkward as the Lc should be the sum of Lps and Lms(?)
! in 2D we have something like:
!
!              *-------*                         *-------*
!              |       |                         |       |
!              |   1   |                         |Lm%y=1 |
!              |       |                         |       |
!      *-------*---O---*-------*         *-------*-------*-------*
!      |       |       |       |         |       |       |       |
!      |   2   X   1   O   1   |  ==>    |Lm%x=0 !Lc = 2 |Lp%x=1 |     
!      |       |       |       |         |       |       |       |
!      *-------*---X---*-------*         *-------*-------*-------*
!              |       |                         |       |
!              |   2   |                         |Lp%y=0 |
!              |       |                         |       |
!              *-------*                         *-------*
! 
! consider the centre cell - Lm%y and Lp%x should be one, while the other 
! Lm/Lps are set to zero. but in that case the Lc should become 2 instead 
! of 4 (?). 
! 
! The only good part (?) is we just need to do it once for all. so, we don't 
! using the complex sparse vector (sparcevecc) to store the exeption - 
! but just modify the laplace operator on-the-fly. 
! 
! I will leave the discussion of whether we will need a sparcevecr class to 
! Gary and Anna. 
!
! NOTE: we keep all the parameters, to keep the file compatible with the AR
! one - some of the following parameters are not used
!
    implicit none
    character(*), intent(in)             :: cfile
!   Exception rules
    integer,pointer,dimension(:)         :: mask1, mask2, ii, jj, kk, xyz
    real(kind=prec),pointer,dimension(:) :: smoothing

    integer :: nrules, nS, i, j, k, n, istat, NxFile, NyFile, NzFile
    integer :: fid = 30
    real(kind=prec) :: tmp

    if (.not. CovAllocated) then
        call errStop('Model covariance must be allocated before reading from file in read_CmSqrt')
    end if

    open(unit=fid,file=cfile,form='formatted',status='old')

    ! skip the 16 lines header
    do j = 1,16
        read(fid,*)
    end do

    ! read grid dimensions: check consistency with grid
    read(fid,*) NxFile,NyFile,NzFile    
    if ((Nx /= NxFile) .or. (Ny /= NyFile) .or. (NzEarth /= NzFile)) then
        call errStop('Grid dimensions do not match in input covariance file '&
            //cfile)
    end if
        
!! this makes no sense  ... we allocated from model parameter, 
!! and we check consistency with model parameter AFTER reading Sx etc.  
!! But if things are not consistent, we will have an error before we check!!
!! NOTE: the above is corrected already 

    ! read smoothing parameters
    read(fid,*) Sx
    read(fid,*) Sy
    read(fid,*) tmp

    Sx = Sx/laplaceHorizontalScale
    Sy = Sy/laplaceHorizontalScale
    Sz(:) = tmp/laplaceVerticalScale
    ! read number of times to apply the smoothing
    ! NOTE: this is not used for laplacian 
    ! or maybe we still want it to switch between L1 and L2? 
    read(fid,*) N

    ! read exception rules for smoothing across surfaces
    read(fid,*) nrules
    allocate(mask1(nrules),mask2(nrules),smoothing(nrules),STAT=istat)
    do n = 1,nrules
        read(fid,*) mask1(n), mask2(n), smoothing(n)
    end do

    ! create and read the mask array
    ! NOTE: this is not used for laplacian 
    call read_iscalar(fid,mask)
    close(fid)

    ! create a *huge* sparse vector to make sure we accommodate all 
    ! smoothing exceptions
    ! NOTE: not used anymore for laplacian
    ! call create_sparsevecc(Nx*Ny*NzEarth, S, FACE)

    ! now, parse the exceptions

    ! the boundary is irretating, so I just stop thinking and make the 
    ! subroutine to do it for multiple times (!)
    ! this might still be simplified...
    do k = 1,NzEarth
        do j = 1,Ny
            do i = 1,Nx
                ! modify the smoothing before applying mask
                Lp%x(i+1,j,k) = Lp%x(i+1,j,k)*Sx(k)
                Lp%y(i,j+1,k) = Lp%y(i,j+1,k)*Sy(k)
                Lp%z(i,j,k+1) = Lp%z(i,j,k+1)*Sz(k)
                Lm%x(i,j,k) = Lm%x(i,j,k)*Sx(k)
                Lm%y(i,j,k) = Lm%y(i,j,k)*Sy(k)
                Lm%z(i,j,k) = Lm%z(i,j,k)*Sz(k)
            end do
        end do
    end do
    ! look ahead
    do k = 1,NzEarth
        do j = 1,Ny
            do i = 1,Nx-1
                do n = 1,nrules
                    ! look ahead in the X-direction
                    if (((mask%v(i+1,j,k) .eq. mask1(n)) .and. &
                        (mask%v(i,j,k) .eq. mask2(n))) .or. &
                        ((mask%v(i+1,j,k) .eq. mask2(n)) .and. &
                        (mask%v(i,j,k) .eq. mask1(n)))) then
                        if (smoothing(n) .eq. 0) then
                            Lp%x(i+1,j,k) = R_ZERO
                        else 
                            ! do nothing currently 
                        endif
                    endif
                end do
            end do
        end do
    end do
    do k = 1,NzEarth
        do j = 1,Ny
            do i = 2,Nx
                do n = 1,nrules
                    ! look back in the X-direction
                    if (((mask%v(i-1,j,k) .eq. mask1(n)) .and. &
                        (mask%v(i,j,k) .eq. mask2(n))) .or. &
                        ((mask%v(i-1,j,k) .eq. mask2(n)) .and. &
                        (mask%v(i,j,k) .eq. mask1(n)))) then
                        if (smoothing(n) .eq. 0) then
                            Lm%x(i,j,k) = R_ZERO
                        else 
                            ! do nothing currently 
                        endif
                    endif
                end do
            end do
        end do
    end do
    do k = 1,NzEarth
        do j = 1,Ny-1
            do i = 1,Nx
                do n = 1,nrules
                    ! look ahead in the Y-direction
                    if (((mask%v(i,j+1,k) .eq. mask1(n)) .and. &
                        (mask%v(i,j,k) .eq. mask2(n))) .or. &
                        ((mask%v(i,j+1,k) .eq. mask2(n)) .and. &
                        (mask%v(i,j,k) .eq. mask1(n)))) then
                        if (smoothing(n) .eq. 0) then
                            Lp%y(i,j+1,k) = R_ZERO
                        else 
                            ! do nothing currently 
                        endif
                    end if
                end do
            end do
        end do
    end do
    do k = 1,NzEarth
        do j = 2,Ny
            do i = 1,Nx
                do n = 1,nrules
                    ! look back in the Y-direction
                    if (((mask%v(i,j-1,k) == mask1(n)) .and. &
                        (mask%v(i,j,k) == mask2(n))) .or. &
                        ((mask%v(i,j-1,k) == mask2(n)) .and. &
                        (mask%v(i,j,k) == mask1(n))))  then
                        if (smoothing(n) .eq. 0) then
                            Lm%y(i,j,k) = R_ZERO
                        else 
                            ! do nothing currently 
                        endif
                    end if
                end do
            end do
        end do
    end do
    ! look back
    do k = 1,NzEarth-1
        do j = 1,Ny
            do i = 1,Nx
                do n = 1,nrules
                    ! look ahead in the Z-direction
                    if (((mask%v(i,j,k+1) .eq. mask1(n)) .and. &
                        (mask%v(i,j,k) .eq. mask2(n))) .or. &
                        ((mask%v(i,j,k+1) .eq. mask2(n)) .and. &
                        (mask%v(i,j,k) .eq. mask1(n)))) then
                        if (smoothing(n) .eq. 0) then
                            Lp%z(i,j,k+1) = R_ZERO
                        else 
                            ! do nothing currently 
                        endif
                    end if
                end do
            end do
        end do
    end do
    do k = 2,NzEarth
        do j = 1,Ny
            do i = 1,Nx
                do n = 1,nrules
                    ! look back in the Z-direction
                    if (((mask%v(i,j,k-1) == mask1(n)) .and. &
                        (mask%v(i,j,k) == mask2(n))) .or. &
                        ((mask%v(i,j,k-1) == mask2(n)) .and. &
                        (mask%v(i,j,k) == mask1(n)))) then
                        if (smoothing(n) .eq. 0) then
                            Lm%z(i,j,k) = R_ZERO
                        else 
                            ! do nothing currently 
                        endif
                    end if
                end do
            end do
        end do
    end do
    ! re-calculate the Lc to compensate 
    ! now behold the might of fortran (?)
    Lc%v(1:Nx,1:Ny,1:NzEarth) = -(Lm%x(1:Nx,1:Ny,1:NzEarth) &
                           +Lp%x(2:Nx+1,1:Ny,1:NzEarth) &
                           +Lm%y(1:Nx,1:Ny,1:NzEarth)   &
                           +Lp%y(1:Nx,2:Ny+1,1:NzEarth) &
                           +Lm%z(1:Nx,1:Ny,1:NzEarth)   &
                           +Lp%z(1:Nx,1:Ny,2:NzEarth+1))

    deallocate(mask1,mask2,smoothing)

    end subroutine read_Cm_Laplacian

   subroutine A(x,y)
   !     computes y = -DelSq*x   ... y needs to be created/allocated
   !            before calling ; let's just assume compatible
   !    again assuming this is all part of the module containing model
   !          covariance
       type(rscalar), intent(in)            :: x
       type(rscalar), intent(inout)         :: y
       real (kind = prec)                   :: s, s2
       integer ix,iy,iz,Nz
       Nz = NzEarth
       y%v(:,:,:) = Lc%v(:,:,:)*x%v(:,:,:)

       !  the boundaries complicate things in Fortran!
       !    maybe this could be made more efficient ...
       do ix = 1,Nx-1 
          y%v(ix,:,:) = y%v(ix,:,:) + Lp%x(ix+1,:,:)*x%v(ix+1,:,:) 
          y%v(ix+1,:,:) = y%v(ix+1,:,:) + Lm%x(ix+1,:,:)*x%v(ix,:,:)
       enddo
       do iy = 1,Ny-1 
          y%v(:,iy,:) = y%v(:,iy,:) + Lp%y(:,iy+1,:)*x%v(:,iy+1,:) 
          y%v(:,iy+1,:) = y%v(:,iy+1,:) + Lm%y(:,iy+1,:)*x%v(:,iy,:)
       enddo
       do iz = 1,Nz-1 
          y%v(:,:,iz) = y%v(:,:,iz) + Lp%z(:,:,iz+1)*x%v(:,:,iz+1) 
          y%v(:,:,iz+1) = y%v(:,:,iz+1) + Lm%z(:,:,iz+1)*x%v(:,:,iz)
       enddo
       !  add constant
       s = 0.0
       s2 = 0.0
       do ix = 1,Nx
          do iy = 1,Ny
             do iz = 1,Nz
                 if (mask%v(ix,iy,iz) == AIR) then
                     s = s + 0.0
                 else if (mask%v(ix,iy,iz) == OCEAN) then
                     s = s + 0.0
                 else
                     s = s + x%v(ix,iy,iz)
                     s2 = s2 + 1.0
                 end if
             enddo
          enddo
      enddo
      ! s = s*LaplaceConstantScale/s2
      ! y%v(:,:,:) = y%v(:,:,:)+s

   end subroutine  ! Laplace
  !**********************************************************************
  ! to deallocate the coefficients used for divergence correction
  subroutine deall_Laplace()

    implicit none

    Call deall_rvector(Lm)
    Call deall_rvector(Lp)
    Call deall_rscalar(Lc)

  end subroutine ! deall_Laplace 
!   for now I am just including the CG solver within the Laplace routines.
!         this is almost identical to PCG (from module solver) but is for
!           real scalar fields, and the interface is setup for the Poisson
!           equation; also, so far there is no preconditioner; this could be
!             simplified, or a preconditioner developed

  subroutine PCG_Poisson(b,x)!  ,PCG_iter)   note that PCG_iter is now a module
  !                                           variable and does not need to be
  !                                           passed
  ! Purpose: a quasi-generic pre-conditioned conjugate gradient
  ! routine, set up for solving DivCgrad phiSOL = phiEMrhs using
  ! routines in divCor.  Actual code is generic, but the interface
  ! is not

   ! inteface..............
   ! redefining some of the interfaces for our convenience (locally)
   ! generic routines for scalar operations on corner/ center nodes
   ! on a staggered grid
    use sg_scalar

  ! routines for divergence correction
  ! use modeloperator3D, only:  A => Laplacian, Minv => Ident

  implicit none
  type (rscalar), intent(in)            :: b
  type (rscalar), intent(inout)         :: x
  ! type (solverControl_t), intent(inout)         :: PCGiter

  ! local variables
  type (rscalar)        :: r,s,p,q
  real(kind=prec)    :: beta,alpha,delta,deltaOld
  real(kind=prec)    :: bnorm, rnorm
  integer               :: i

  if (.not.b%allocated) then
      write(0,*) 'b in PCG not allocated yet'
      stop
  end if

  if (.not.x%allocated) then
      write(0,*) 'x in PCG not allocated yet'
      stop
end if

  ! Allocation of r, z, p, q
  Call create(x%grid, r, x%gridType)
  Call create(x%grid, s, x%gridType)
  Call create(x%grid, p, x%gridType)
  Call create(x%grid, q, x%gridType)

  bnorm = dotProd(b,b)
  if (isnan(bnorm)) then
      ! this usually means an bad RHS, in which case the pcg fails
      write(6,*) 'ERROR b in PCG contains NaNs; exiting...'
      stop
  else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
      write(6,*) 'WARNING: b in PCG has all zeros, returning zero solution'
      x = b
      LaplaceNiter = 1
      Laplacefailed = .false.
      LaplaceRerr(1) = 0.0
      return
  endif
  Call A(x,r)
  Call linComb(one,b,MinusOne,r,r)
  rnorm = dotProd(r,r)
  i = 1
  LaplaceRerr(i) = real(rnorm/bnorm)

  loop: do while ((LaplaceRerr(i).gt.LaplaceTol).and.(i.lt.LaplaceMaxIt))
     !Call Minv(r,s)
     !delta = dotProd(r,s)
     !    at present no preconditioner -- should be able to simplify ... or
     !    develop   a preconditioner
     s  = r 
     delta = rnorm
     if(i.eq.1) then
        p = s
     else
        beta = delta/deltaOld
        Call linComb(one,s,beta,p,p)
     end if

     Call A(p,q)
     alpha = delta/dotProd(p,q)
     Call scMultAdd(alpha,p,x)
     Call scMultAdd(-alpha,q,r)
     deltaOld = delta
     i = i + 1
     rnorm = dotProd(r,r)
     LaplaceRerr(i) = rnorm/bnorm

  end do loop
  LaplaceNiter = i

  ! deallocate all the work arrays
  Call deall(r)
  Call deall(s)
  Call deall(p)
  Call deall(q)

end subroutine PCG_Poisson

  !**********************************************************************
  ! apply pre-conditioner, solving lower and upper triangular systems using
  ! coefficients in Lm, Lp, and Ld.
  subroutine Minv(r, rPrime)
  ! essentially this calculates rPrime = M^-1 * r
  ! by solving M*rPrime=r, where M=LU, 
  ! 1.  solves L * y = r
  ! 2.  solves U * rPrime = y

    implicit none
    type (rscalar), intent(in)                :: r
    type (rscalar), intent(inout)             :: rPrime
    integer                                   :: ix, iy, iz

    IF(.not.r%allocated) THEN
         WRITE(0,*) 'r not allocated in Minv'
         STOP
    ENDIF

    IF(.not.rPrime%allocated) THEN
         WRITE(0,*) 'rPrime not allocated in Minv'
         STOP
    ENDIF

    if ((r%nx == rPrime%nx).and.(r%ny == rPrime%ny).and.&
         (r%nz == rPrime%nz)) then
        ! check if the r and rPrime are of the same size
        if (r%gridType == rPrime%gridType) then
           rPrime%v = 0.0
          ! forward substitution (Solve lower triangular system)
          ! the coefficients are only for the interior nodes
          do iz = 2, r%nz
             do iy = 2, r%ny
                do ix = 2, r%nx
                    rPrime%v(ix, iy, iz) = r%v(ix, iy, iz) &
                        - rPrime%v(ix-1,iy,iz)*Lm%x(ix,iy,iz)&
                        *Ld%v(ix-1,iy,iz) &
                        - rPrime%v(ix,iy-1,iz)*Lm%y(ix,iy,iz)&
                        *Ld%v(ix,iy-1,iz) &
                        - rPrime%v(ix,iy,iz-1)*Lm%z(ix,iy,iz)&
                        *Ld%v(ix,iy,iz-1)
                enddo
             enddo
          enddo
          ! backward substitution (Solve upper triangular system)
          ! the coefficients are only for the interior nodes
          do iz = r%nz,2,-1
             do iy = r%ny,2,-1
                do ix = r%nx,2,-1
                    rPrime%v(ix, iy, iz) = (rPrime%v(ix, iy, iz)  &
                        - rPrime%v(ix+1, iy, iz)*Lp%x(ix, iy, iz)  &
                        - rPrime%v(ix, iy+1, iz)*Lp%y(ix, iy, iz)  &
                        - rPrime%v(ix, iy, iz+1)*Lp%z(ix, iy, iz)) &
                        *Ld%v(ix, iy, iz)
                enddo
             enddo
          enddo
        else
          write (0, *) 'Minv: not compatible existing data types'
        end if
    else
        write(0, *) 'Error: Minv: scalars not even of the same size'
    end if
    end subroutine Minv          ! Minv

!**********************************************************************!

    subroutine writeRvectorBin(r,cfile)
    character(*),intent(in)        :: cfile
    type(rvector), intent(in)        :: r
    integer        :: fid

    open(unit=fid,file=cfile,form='unformatted',status='unknown')
    write(fid) r%nx,r%ny,r%nz
    write(fid) r%x
    write(fid) r%y
    write(fid) r%z
    close(fid)
    end subroutine writeRvectorBin

    subroutine writeRscalarBin(r,cfile)

    character(*),intent(in)        :: cfile
    type(rscalar), intent(in)      :: r
    integer                        :: fid

    open(unit=fid,file=cfile,form='unformatted',status='unknown')
    write(fid) r%nx,r%ny,r%nz
    write(fid) r%v
    close(fid)
    end subroutine writeRscalarBin

     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine H(x,y,lambda1,gamma1)
   !   computes y = H*x   ... y needs to be created/allocated
   !   before calling ; let's just assume compatible
   !   again assuming this is all part of the module containing model
   !   covariance
   !   
   !   in a smoothed (untransfored) space, Hessian should be something like: 
   !           2*(J^T*Cd^-1*J + \lambda*Cm^-1)
   !   to avoid computations with J, we want to approximate Hessian with:
   !           H = (\gamma*I + \lambda*Cm^-1)
   ! 
   !   in Gary's unsmoothed model space, we want to similarly approximate:
   !           2*(Cm^0.5*J^T*Cd^-1*J + \lambda*Cm^-0.5)
   !   with: 
   !           H = (\gamma*I + \lambda*Cm^-0.5)
   ! 
   !   so firstly we calculate yPrime = Cm^-0.5*x
   !   then we compute y = \gamma*x + \lambda*yPrime
      type(rscalar), intent(in)            :: x
      type(rscalar), intent(inout)         :: y
      real(kind=prec), intent(in)          :: lambda1, gamma1

      ! local variables
      type(rscalar)                        :: yPrime
      integer                              :: Nmodel

      ! initialize yPrime
      yPrime = x
      y = x
      select case (flagCov)
          case (1)
              ! firstly we calculate y = L*x
              call RecursiveARInv(x%v, yPrime%v, N)
              ! then we compute yPrime = L^T*y
              call RecursiveARInv(yPrime%v, y%v, N)
          case (2)
              ! firstly we calculate y = L*x
              call A(x,yPrime)
              ! then we compute yPrime = L^T*y
              call A(yPrime,y)
          case (3)
              call A(x,y)
          case default
              write(0,*) 'Unsupported model covariance operator!'
              stop
      end select

      Nmodel = x%Nx * x%Ny * x%Nz
      call linComb(gamma1,x,lambda1/Nmodel,y,y)
      ! deallocate all the work arrays
      call deall(yPrime)
      return
  end subroutine H ! Hessian

!**********************************************************************!

   subroutine Hinv(bbox, xbox, gamma1, lambda1)
    ! this subroutine is used to apply preconditioner to iteration-based 
    ! optimization methods like CG or Quasi-Newton
    ! essentially this constructs a preconditioner P that acts like the 
    ! inverse of Hessian (H^-1) and apply it to vector r:
    !     x = H^-1*b 
    ! here we use Randy's rational, where
    !     H^-1 = (\gamma*I + \lambda*C_m^-1)^-1 = C^-1
    ! we don't have to construct H^-1 explicitly, but just apply it 
    ! by solving the system  
    !     H*x = b
    ! apparently, the operator: 
    !     H = (\gamma*I + \lambda*L^TL)
    ! will be needed, but also not explicitly, we just need y = H(x)
    ! by a (preconditioned) conjugate gradient solver, here gamma is selected 
    ! to be a real value that approximate the matrix of JTCd^-1J
    ! we might just want to do the same as Randy, to use an value that is
    ! based on the Jacobian of a homogeneous half space (or layered model)
    ! 
    ! the solver itself is provided in the Model Covariance module for 
    ! laplacian smoother/rougher

     type (modelParam_t), intent(in)           :: bbox
     type (modelParam_t), intent(inout)        :: xbox
     real (kind=prec),intent(in)               :: gamma1
     real (kind=prec),intent(in)               :: lambda1

     ! local variables 
     type (rscalar)                            :: b,x,r,s,p,q
     real(kind=prec)                           :: beta,alpha,delta,deltaOld
     real(kind=prec)                           :: bnorm, rnorm
     integer                                   :: i

     if (.not.bbox%allocated) then
         write(0,*) 'b in PCG (H^-1) not allocated yet'
         stop
     end if

     if (.not.xbox%allocated) then
         write(0,*) 'x in PCG (H^-1) not allocated yet'
         stop
     end if
     ! Allocation of x and b
     Call create(bbox%cellCond%grid, b, bbox%cellCond%gridType)
     b%v = bbox%cellCond%v
     Call create(b%grid, x, b%gridType)
     x%v = R_ZERO
     ! Allocation of r, z, p, q
     Call create(x%grid, r, x%gridType)
     Call create(x%grid, s, x%gridType)
     Call create(x%grid, p, x%gridType)
     Call create(x%grid, q, x%gridType)
     bnorm = dotProd(b,b)
     if (isnan(bnorm)) then
     ! this usually means an bad RHS, in which case the pcg fails
         write(6,*) 'ERROR b in PCG(H^-1) contains NaNs; exiting...'
         stop  
     else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
         write(6,*) 'WARNING: b in PCG has all zeros, returning zero solution'
         x = b
         LaplaceNiter = 1
         Laplacefailed = .false.
         LaplaceRerr(1) = 0.0
         return
     endif
     Call H(x,r,gamma1,lambda1)
     Call linComb(one,b,MinusOne,r,r)
     rnorm = dotProd(r,r)
     i = 1
     LaplaceRerr(i) = real(rnorm/bnorm)
     ! now start looping 
     loop: do while ((LaplaceRerr(i).gt.LaplaceTol).and.(i.lt.LaplaceMaxIt))
         ! Call Minv(r,s)
         ! delta = dotProd(r,s)
         ! at present no preconditioner -- should be able to simplify ... or
         ! develop a preconditioner
         s = r 
         delta = rnorm
         if(i.eq.1) then
             p = s
         else
             beta = delta/deltaOld
             Call linComb(one,s,beta,p,p)
         end if
         Call H(p,q,gamma1,lambda1)
         alpha = delta/dotProd(p,q)
         Call scMultAdd(alpha,p,x)
         Call scMultAdd(-alpha,q,r)
         deltaOld = delta
         i = i + 1
         rnorm = dotProd(r,r)
         LaplaceRerr(i) = rnorm/bnorm
     end do loop
     LaplaceNiter = i
     
     ! write(6,*) 'number of Hinv iterations: ', LaplaceNiter
     xbox%cellCond%v = x%v
     ! deallocate all the work arrays
     Call deall(x)
     Call deall(b)
     Call deall(r)
     Call deall(s)
     Call deall(p)
     Call deall(q)

end subroutine Hinv !H^-1

