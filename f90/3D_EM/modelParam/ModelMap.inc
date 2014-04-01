!    routines which define mappings between the "natural" representation
!        of conductivity/resistivity on the model grid and the formal
!        model parameter structure: in general (allowing for a non-linear
!        mapping, as will occur when log condutivity is used for inversion;
!        in the implementation here,we combine the linear and non-linear
!        mappings in a single routine)
!        one needs the non-linear mapping, a linearized mapping, and the
!        adjoint of the linearized mapping.  If data functionals depend
!        on conductivity or resistivity additional routines are needed--
!        i.e., a real function (SigC here), with input arguments a conductivity
!        parameter and cell indicies, which returns the conductivity of
!        this cell (or resistivity), as well as what is essentially the
!        adjoint of the linearization of this function (sort of ...
!         see below; this is called QtoModelParam here.)
!       Note that the functionality is always needed, but names can be
!        flexible ... these routines are called by the interpolation
!        routines in EMfieldInterp* (for the latter two only), by the forward
!        modeling codes (the non-linear mapping only), and by SolverSens/DataSens
!        (the linearized mappings and adjoint) and thus names/interfaces
!        need to be kept consistent with what is used in these routines.
!        Note that for 2D MT different mapping routines are used for
!        TE and TM modes, as the natural representations for the two cases
!        are different (conductivity vs. resistivity)

   !**********************************************************************
   subroutine ModelParamToCell(m,cCond,paramType,grid,AirCond)

     type(modelParam_t), intent(in)        :: m
     type(rscalar), intent(inout)	       :: cCond
     character(80), intent(out), optional  :: paramType
     type(grid_t), intent(out), optional :: grid
     real(kind=prec), intent(out), optional :: AirCond

     if(cCond%allocated) then
        if((cCond%Ny .ne. m%Ny).or. (cCond%Nx .ne. m%Nx) .or. &
		(cCond%Nz .ne. m%grid%Nz)) then
           call deall_rscalar(cCond)
           call create_rscalar(m%grid,cCond,CENTER)
        endif
     else
        call create_rscalar(m%grid,cCond,CENTER)
     endif
     if(m%paramType .eq. LOGE) then
        cCond%v(:,:,1:m%grid%NzAir) = exp(m%AirCond)
        cCond%v(:,:,m%grid%NzAir+1:m%grid%Nz) = exp(m%cellCond%v)
     else
        cCond%v(:,:,1:m%grid%NzAir) = m%AirCond
        cCond%v(:,:,m%grid%NzAir+1:m%grid%Nz) = m%cellCond%v
     endif

     if (present(paramType)) then
        paramType = m%paramType
     end if

     if (present(grid)) then
        if (grid%allocated) then
           call deall_grid(grid)
        end if
        grid = m%grid
     end if

     if (present(AirCond)) then
        AirCond = m%AirCond
     end if

   end subroutine ModelParamToCell
  !**********************************************************************

  subroutine dModelParamToEdge(m, Econd,m0)
  !  Maps conductivity defined on grid cells (prisms) to edge-nodes ...
  !  Linear interface to ModelParamToEdge which makes the code more readable.

    implicit none
    !  INPUTS:  model parameter
    type(modelParam_t), intent(in)            :: m
    !  OUTPUTS: conductivities on edge nodes, as needed for FD calc.
    !    stored as rvector; allocate Econd before calling
    type(rvector), intent(inout)         :: Econd
    !  INPUT: background model parameter for linearization
    type(modelParam_t),  intent(in)	     :: m0

    call ModelParamToEdge(m, Econd,m0)

  end subroutine dModelParamToEdge

  !**********************************************************************

  subroutine ModelParamToEdge(m, Econd,m0)
  !  Maps conductivity defined on grid cells (prisms) to edge-nodes ... this
  !   can be used for both non-linear and linearized mappings when
  !   paramtype = LOGE.  In this case, if optional argument m0 is present
  !   the linearized mapping is used; otherwise the non-linear mappring is used.
  !   If paramtype = LINEAR (linear conductivity) then the same linear mapping
  !   is done in whether m0 is present or not (and in fact, m0 is not referenced
  !   in this case)
  !  NOTE: there is a subtlty (and possible source of error) associated with
  !   treatement of air conductivity.  This is (in the present implementation)
  !   viewed as a non-adjustable parameter, but when this routine is called for
  !   the full (i.e., in general non-linear) mapping it must set the edges
  !   in the air to airCond.  On the other hand, when called to implement the
  !   linearization (i.e., to map perturbations in the model parameter) it
  !   should leave air edges set to zero.  WHen the conductivity mapping is
  !   linear, there is no way to distinguish between whether this is a linear
  !   mapping of the perturbation, or the full linear mapping.
  !   One should thus set airCond in perturbations to zero (or they will
  !   be propagated to perturbation in the air edges, only for the linear
  !   model parameter case).

    implicit none
    !  INPUTS:  model parameter
    type(modelParam_t), intent(in)            :: m
    !  OUTPUTS: conductivities on edge nodes, as needed for FD calc.
    !    stored as rvector; allocate Econd before calling
    type(rvector), intent(inout)         :: Econd
    !   OPTIONAL INPUT: background model parameter for linearization
    type(modelParam_t), optional, intent(in)	:: m0

    !   LOCAL Variables
    type(rscalar)  :: temp,Ctemp,V_C
    type(rvector)  :: Etemp,V_E
    integer		:: ix,iy,iz,ize,nx,ny,nz,nzEarth,nzAir
    logical		:: linearizedMapping
    real(kind=prec)   :: airCond
 
    linearizedMapping = present(m0)

    if (.not. m%allocated) then
      write(0,*) 'input model parameter in ModelParamToEdge not allocated yet'
      stop
    end if

    if (.not. Econd%allocated) then
      call create_rvector(m%grid,Econd,EDGE)
    end if

    ! Could add code to check whether the bounds are the same ...

    nx = m%Nx
    ny = m%Ny
    nz = m%grid%nz
    nzEarth = m%nzEarth
    nzAir = nz-nzEarth

    airCond = m%airCond
    call create_rscalar(m%grid,temp,CENTER)
    call setup_grid(temp%grid)

    if(m%paramType .EQ. LOGE) then
       if(linearizedMapping) then
          temp%v(:,:,nzAir+1:nz) = m%cellCond%v*exp(m0%cellCond%v)
          airCond = R_ZERO
          temp%v(:,:,1:nzAir) = airCond
       else
          temp%v(:,:,nzAir+1:nz) = exp(m%cellCond%v)
          airCond = exp(m%airCond)
          temp%v(:,:,1:nzAir) = airCond
       endif
    endif

	call CellVolume(m%grid,V_C)
	call EdgeVolume(m%grid,V_E)	
	Ctemp = temp
	call diagMult(V_C,temp,Ctemp)
	call Cell2Edge(m%grid,Ctemp,Etemp)
	Econd = Etemp
	call diagDiv(Etemp,V_E,Econd)

	call deall_rscalar(V_C)
	call deall_rvector(V_E)
	call deall_rscalar(temp)
	call deall_rscalar(Ctemp)
	call deall_rvector(Etemp)

  end subroutine ModelParamToEdge

  !**********************************************************************

  subroutine dEdgeToModelParam(Econd,m,m0)
  !  Maps from a real vector (defined on edges) to modelParam;
  !    the adjoint of linear mapping implemented by ModelParamToEdge
  !  Interface to EdgeToModelParam that makes the code more readable.

    implicit none
    !  INPUTS:  real vector defined on edges
    type(rvector), intent(in)            :: Econd
    !  OUTPUTS: model parameter
    type(modelParam_t), intent(inout)         :: m
    !  INPUT background model parameter
    type(modelParam_t), intent(in)	    :: m0

    call EdgeToModelParam(Econd,m,m0)

    m%zeroValued = .false.

  end subroutine dEdgeToModelParam

  !**********************************************************************

  subroutine EdgeToModelParam(Econd,m,m0)
  !  Maps from a real vector (defined on edges) to modelParam;
  !    the adjoint of linear mapping implemented by ModelParamToEdge
  !
    implicit none
    !  INPUTS:  real vector defined on edges
    type(rvector), intent(in)            :: Econd
    !  OUTPUTS: model parameter
    type(modelParam_t), intent(inout)         :: m
    !  INPUT (OPTIONAL) background model parameter,
    !         required if m%paramtype=LOGE
    type(modelParam_t), optional, intent(in)	:: m0

    !   LOCAL Variables
    type(rscalar)  :: temp,Ctemp,V_C
    type(rvector)  :: Etemp,V_E
    integer		:: ix,iy,iz,ize,nx,ny,nz,nzEarth,nzAir

    if ((.not.m%allocated).or.(.not.Econd%allocated)) then
      call errStop('m or Econd not allocated yet in EdgeToModelParam')
      stop
    end if

    if((m%paramtype.eq.LOGE).and. (.not.present(m0))) then
       call errStop('Background conductivity required for paramType LOGE in EdgeToModelParam')
    endif

    ! Could add code to check whether the bounds are the same ...

    nx = m%Nx
    ny = m%Ny
    nz = m%grid%nz
    nzEarth = m%nzEarth
    nzAir = nz-nzEarth

    call create_rscalar(Econd%grid,temp,CENTER)
    
	call CellVolume(m%grid,V_C)
	call EdgeVolume(m%grid,V_E)	
	Etemp = Econd
	call diagDiv(Econd,V_E,Etemp)
	call Edge2Cell(m%grid,Etemp,Ctemp)
	temp = Ctemp
	call diagMult(V_C,Ctemp,temp)


    if(m%paramType .EQ. LOGE) then
       m%cellCond%v = temp%v(:,:,nzAir+1:nz)*exp(m0%cellCond%v)
    else
       m%cellCond%v = temp%v(:,:,nzAir+1:nz)
    endif

    m%zeroValued = .false.

	call deall_rscalar(V_C)
	call deall_rvector(V_E)
	call deall_rscalar(temp)
	call deall_rscalar(Ctemp)
	call deall_rvector(Etemp)

  end subroutine EdgeToModelParam

  !**********************************************************************
  function ModelParamToOneEdge(m,xyz,ix,iy,iz) result(r)
  !   computes conductivity for edge xyz/i,j,k using input modelParam
  !    structure m.  This function defines how the abstract
  !     conductivity parameter is mapped to edge conductivities needed
  !     for more accurate electric field interpolation.  The derivative
  !     of this function is required for evaluation of linearized
  !     data functionals, and for construction of the direct
  !     parameter space part of the comb

    type (modelParam_t), intent(in)     :: m
    integer, intent(in)               :: xyz,ix,iy,iz
    real(kind=prec)           :: r

    ! local variables
    integer		:: izE,nx,ny,nz,nzEarth,nzAir
    real(kind=prec)           :: w11,w21,w12,w22,S,temp(2,2)

    nx = m%Nx
    ny = m%Ny
    nz = m%grid%nz
    nzEarth = m%nzEarth
    nzAir = nz-nzEarth

    izE = iz-nzAir
    selectcase(xyz)
    case(1)		 ! Ex edge
       w11 = m%grid%dy(iy)*m%grid%dz(iz)
       temp(1,1) = m%cellCond%v(ix,iy,izE)
       w21 = m%grid%dy(iy-1)*m%grid%dz(iz)
       temp(2,1) = m%cellCond%v(ix,iy-1,izE)
       w12 = m%grid%dy(iy)*m%grid%dz(iz-1)
       temp(1,2) = m%cellCond%v(ix,iy,izE-1)
       w22 = m%grid%dy(iy-1)*m%grid%dz(iz-1)
       temp(2,2) = m%cellCond%v(ix,iy-1,izE-1)
       S = w11+w21+w12+w22
       if(m%paramtype .eq. LOGE) then
          temp = exp(temp)
       endif
       if(izE.eq.1) then
          r = (temp(1,1)*w11+temp(2,1)*w21)/S
       else
          r = (temp(1,1)*w11 + temp(2,1)*w21+    &
               temp(1,2)*w12 + temp(2,2)*w22)/S
       endif

    case(2)		 ! Ey edge
       w11 = m%grid%dx(ix)*m%grid%dz(iz)
       temp(1,1) = m%cellCond%v(ix,iy,izE)
       w21 = m%grid%dx(ix-1)*m%grid%dz(iz)
       temp(2,1) = m%cellCond%v(ix-1,iy,izE)
       w12 = m%grid%dx(ix)*m%grid%dz(iz-1)
       temp(1,2) = m%cellCond%v(ix,iy,izE-1)
       w22 = m%grid%dx(ix-1)*m%grid%dz(iz-1)
       temp(2,2) = m%cellCond%v(ix-1,iy,izE-1)
       S = w11+w21+w12+w22
       if(m%paramtype .eq. LOGE) then
          temp = exp(temp)
       endif
       if(izE.eq.1) then
          r = (temp(1,1)*w11+temp(2,1)*w21)/S
       else
          r = (temp(1,1)*w11 + temp(2,1)*w21+ &
               temp(1,2)*w12 + temp(2,2)*w22)/S
       endif

    case(3)		 ! Ez edge
       w11 = m%grid%dx(ix)*m%grid%dy(iy)
       temp(1,1) = m%cellCond%v(ix,iy,izE)
       w21 = m%grid%dx(ix)*m%grid%dy(iy-1)
       temp(2,1) = m%cellCond%v(ix,iy-1,izE)
       w12 = m%grid%dx(ix-1)*m%grid%dy(iy)
       temp(1,2) = m%cellCond%v(ix-1,iy,izE)
       w22 = m%grid%dx(ix-1)*m%grid%dy(iy-1)
       temp(2,2) = m%cellCond%v(ix-1,iy-1,izE)
       S = w11+w21+w12+w22
       if(m%paramtype .eq. LOGE) then
          temp = exp(temp)
       endif
       r = (temp(1,1)*w11 + temp(2,1)*w21+  &
            temp(1,2)*w12 + temp(2,2)*w22)/S

    endselect
  end function ModelParamToOneEdge
