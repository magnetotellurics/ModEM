! *****************************************************************************
! This module implements boundary conditions and first guess as in WS3D
! subroutine bound3d (in boundary.f)
module boundary_ws

  use sg_vector
  use sg_scalar
  use sg_boundary
  use fwdtemod

  implicit none

  ! Workhorse
  public		:: BC_x0_WS

Contains


  !****************************************************************************
  ! Generates boundary conditions and initial guess for the iterative solution.
  subroutine BC_x0_WS(imode,period,grid3D,Cond3D,E0,BC)
    !  sets up boundary condtions in BC and initial electric
    !  field solution vector in E0
    !  E0 should be complex vector, allocated, gridType = EDGE
    !  BC should be complex BC, allocated

    !  Input mode, period
    integer, intent(in)		:: imode
    real(kind=prec)	:: period
    !  Input 3D grid
    type(grid_t), intent(in)	:: grid3D
    !  Input 3D conductivity in cells
    type(rscalar), intent(in)	:: Cond3D

    ! Output electric field first guess (for iterative solver)
    type(cvector), intent(inout)	:: E0
    ! Output boundary conditions
    type(cboundary), optional, intent(inout)	:: BC


    ! local variables
    ! 2D grid definitions
    type(grid2d_t)		        :: grid2D
    integer				:: iSlice,ih,iz,nEXB,nv,nh,nSlice
    integer				:: IER,istat

    ! These are temporary work arrays
    ! Array for 2D solutions
    complex (kind=prec),allocatable, dimension(:,:)	:: Esol
    ! Array for 2D BC
    complex (kind=prec), allocatable, dimension(:)      :: EXB
    ! Array for 2D Conductivity
    real (kind=prec), allocatable, dimension(:,:)	:: Cond2D

    !  nv is number of vertical levels in 2D
    nv = grid3D%nz
    grid2D%nz = nv
    grid2D%nza = grid3D%nzAir
    allocate(grid2D%Dz(nv))
    do iz = 1,nv
       grid2D%Dz(iz) = grid3D%Dz(iz)
    enddo

    ! allocate local 2D arrays ... size depends on mode
    ! nh is number of horizontal grid cells
    ! nSlice is number of slices to do 2D computations for
    ! Definitions: imode = 1 Ey polarization and imode = 2 Ex polarization
    if(imode .eq. 2) then
       ! Ex-polarization
       nh = grid3D%ny
       nSlice = grid3D%nx
    else
       !  just assume imode.eq. 1 (Ey-polarization)
       nh = grid3D%nx
       nSlice = grid3D%ny
    endif
    grid2D%Ny = nh
    nEXB = 2*(nh+1)+2*(nv+1)
    allocate(grid2D%Dy(nh))
    allocate(Cond2D(nh,nv))
    allocate(Esol(nh+1,nv+1))
    allocate(EXB(nEXB))

    !  copy Dy, Cond3D for first slice into 2D arrays
    ! Definitions: imode = 1 Ey polarization and imode = 2 Ex polarization
    if(imode .eq. 2) then
       do ih = 1,nh
          grid2D%Dy(ih) = grid3D%dy(ih)
       enddo
    else
       !  just assume imode.eq.1
       do ih = 1,nh
          grid2D%Dy(ih) = grid3D%dx(ih)
       enddo
    endif

    !  just start with dummy initialization of conductivity
    Cond2D(:,:) = 1.d+0

    Call setWSparams(nh,nv)
    Call FWD2DSetupTE(grid2D,IER)
    Call zero(E0)

    do iSlice = 1,nSlice
       ! now set conductivity for eacy slice ... differently for
       ! x/y modes
       ! Definitions: imode = 1 Ey polarization and imode = 2 Ex polarization
       if(imode .eq. 2) then
          do ih = 1,nh
             do iz = 1,nv
                Cond2D(ih,iz) = Cond3D%v(iSlice,ih,iz)
             enddo
          enddo
       else
          !  just assume imode.eq.1
          do ih = 1,nh
             do iz = 1,nv
                Cond2D(ih,iz) = Cond3D%v(ih,iSlice,iz)
             enddo
          enddo
       endif
       ! using conductivity for this slice set coefficients
       Call UpdateCondTE(Cond2D)
       ! finish setting up operator and factor matrix
       Call UpdateFreqTE(period)
       ! set BC for 2D problem
       Call SetBoundTE(period,EXB)
       ! solve 2D TE equations
       Call Fwd2DsolveTE(EXB,Esol,IER)

       ! copy Esol into a slice of E0
       ! Definitions: imode = 1 Ey polarization and imode = 2 Ex polarization
       if(imode .eq. 2) then
          !  only x components are non-zero ...
          do ih = 1,nh+1
             do iz = 1,nv+1
                E0%x(iSlice,ih,iz) = eSol(ih,iz)
             enddo
          enddo
       else
          !  just assume imode.eq.1 ...
          !    now only y components are non-zero ...
          do ih = 1,nh+1
             do iz = 1,nv+1
                E0%y(ih,iSlice,iz) = eSol(ih,iz)
             enddo
          enddo
       endif
    enddo

    !  Extract BC data from E0 ...
    if (present(BC)) then
        Call getBC(E0,BC)
    end if
    Call Fwd2DdeallTE()
    BC%xZMax = 0.0
    BC%yZMax = 0.0

    ! clean up of the temporary work arrays
    deallocate(grid2D%Dy, STAT=istat)
    deallocate(grid2D%Dz, STAT=istat)
    deallocate(Cond2D, STAT=istat)
    deallocate(Esol, STAT=istat)
    deallocate(EXB, STAT=istat)

    return

  end subroutine BC_x0_WS


end module boundary_ws
