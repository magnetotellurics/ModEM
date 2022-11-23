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
     type(rscalar), intent(inout)	       :: cCond(3)
     character(80), intent(out), optional  :: paramType
     type(grid_t), intent(out), optional :: grid
     real(kind=prec), intent(out), optional :: AirCond
     integer ia

     do ia = 1,3
       if(cCond(ia)%allocated) then
          if((cCond(ia)%Ny .ne. m%Ny).or. (cCond(ia)%Nx .ne. m%Nx) .or. &
		      (cCond(ia)%Nz .ne. m%grid%Nz)) then
             call deall_rscalar(cCond(ia))
             call create_rscalar(m%grid,cCond(ia),CENTER)
          endif
       else
          call create_rscalar(m%grid,cCond(ia),CENTER)
       endif
       if(m%paramType .eq. LOGE) then
          cCond(ia)%v(:,:,1:m%grid%NzAir) = exp(m%AirCond)
          cCond(ia)%v(:,:,m%grid%NzAir+1:m%grid%Nz) = exp(m%cellCond(ia)%v)
       else
          cCond(ia)%v(:,:,1:m%grid%NzAir) = m%AirCond
          cCond(ia)%v(:,:,m%grid%NzAir+1:m%grid%Nz) = m%cellCond(ia)%v
       endif
     enddo

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
    type(rscalar)  :: temp(3),Ctemp(3)
    type(rvector)  :: Etemp,V_E_ghost
    integer		:: ix,iy,iz,ize,nx,ny,nz,nzEarth,nzAir,ioNum,ios
    logical		:: linearizedMapping
    real(kind=prec)   :: airCond
    integer ia
 
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
    do ia = 1,3
      call create_rscalar(m%grid,temp(ia),CENTER)
      call setup_grid(temp(ia)%grid)

      if(m%paramType .EQ. LOGE) then
         if(linearizedMapping) then
            temp(ia)%v(:,:,nzAir+1:nz) = m%cellCond(ia)%v*exp(m0%cellCond(ia)%v)
            airCond = R_ZERO
            temp(ia)%v(:,:,1:nzAir) = airCond
         else
            temp(ia)%v(:,:,nzAir+1:nz) = exp(m%cellCond(ia)%v)
            airCond = exp(m%airCond)
            temp(ia)%v(:,:,1:nzAir) = airCond
         endif
      endif
    enddo

	! V_C and V_E are global variables stored in GridCalc, initialized in initSolver;
	! for efficiency, do not recompute unless called outside of solver for testing.
	if (.not. V_C%allocated) then
		call CellVolume(m%grid,V_C) ! global variable
	end if
	if (.not. V_E%allocated) then
		call EdgeVolume(m%grid,V_E) ! global variable
	end if
	
	! for spherical coords to be strictly correct, need to also average the volumes
	! instead of using the edge volumes V_E, use the "ghost" averaged version...
	! which also has full volumes on boundary edges as if the grid extended further
	call Cell2Edge(m%grid,V_C,V_E_ghost) 	
	
	do ia = 1,3
	  Ctemp(ia) = temp(ia)
	  call diagMult(V_C,temp(ia),Ctemp(ia))
	enddo
	call Cell2EdgeArray(m%grid,Ctemp,Etemp)
	
	Econd = Etemp
	call diagDiv(Etemp,V_E_ghost,Econd)
	
	call deall_rvector(V_E_ghost)
	do ia = 1,3
	  call deall_rscalar(temp(ia))
	  call deall_rscalar(Ctemp(ia))
	enddo
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
    type(rscalar)  :: temp(3),Ctemp(3)
    type(rvector)  :: Etemp,V_E_ghost
    integer		:: ix,iy,iz,ize,nx,ny,nz,nzEarth,nzAir,ioNum,ios
    integer ia

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
    
    do ia = 1,3
      call create_rscalar(Econd%grid,temp(ia),CENTER)
    enddo
    
	! V_C and V_E are global variables stored in GridCalc, initialized in initSolver;
	! for efficiency, do not recompute unless called outside of solver for testing.
	if (.not. V_C%allocated) then
		call CellVolume(m%grid,V_C) ! global variable
	end if
	if (.not. V_E%allocated) then
		call EdgeVolume(m%grid,V_E) ! global variable
	end if
	
	! for spherical coords to be strictly correct, need to also average the volumes
	! instead of using the edge volumes V_E, use the "ghost" averaged version...
	! which also has full volumes on boundary edges as if the grid extended further
	call Cell2Edge(m%grid,V_C,V_E_ghost) 
	
	Etemp = Econd
	call diagDiv(Econd,V_E_ghost,Etemp)
	call Edge2CellArray(m%grid,Etemp,Ctemp)
	do ia = 1,3
	  temp(ia) = Ctemp(ia)
	  call diagMult(V_C,Ctemp(ia),temp(ia))
	  
    if(m%paramType .EQ. LOGE) then
       m%cellCond(ia)%v = temp(ia)%v(:,:,nzAir+1:nz)*exp(m0%cellCond(ia)%v)
    else
       m%cellCond(ia)%v = temp(ia)%v(:,:,nzAir+1:nz)
    endif
    	  
	enddo


    m%zeroValued = .false.

	call deall_rvector(V_E_ghost)
	do ia = 1,3
	  call deall_rscalar(temp(ia))
	  call deall_rscalar(Ctemp(ia))
	enddo
	call deall_rvector(Etemp)

  end subroutine EdgeToModelParam
