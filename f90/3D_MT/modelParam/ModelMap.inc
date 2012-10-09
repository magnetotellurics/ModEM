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
   subroutine ModelParamToCell(m ,cCond,paramType,grid,AirCond)

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

  ! **********************************************************************
  subroutine dModelParamToEdge(m, Econd,m0)
  !  Maps conductivity defined on grid cells (prisms) to edge-nodes ...
  !  Linear interface to ModelParamToEdge which makes the code more readable.

    implicit none
    !  INPUTS:  model parameter
    type(modelParam_t), intent(in)            :: m
    !  OUTPUTS: conductivities on edge nodes, as needed for FD calc.
    !    stored as rvector; allocate Econd before calling
    type(rvector_mg), intent(inout)         :: Econd
    !  INPUT: background model parameter for linearization
    type(modelParam_t),  intent(in)	     :: m0

    call ModelParamToEdgeMG(m, Econd,m0)

  end subroutine dModelParamToEdge

  !**********************************************************************
  ! subroutine for destribution the model parameters to edges of multigrid
  ! arithmetic averaging of the m.p. for coarser grid
  ! LINEARIZED MAPPING case need to be checked

    subroutine ModelParamtoEdgeMG(m, Econd,m0)

    implicit none
    !  INPUTS:  model parameter
    type(modelParam_t), intent(in)  :: m
    !  OUTPUTS: conductivities on edge nodes of multigrids, as needed for FD calc.
    ! stored as rvector_mg; allocate Econd before calling; Econd on multigrid
    type(rvector_mg), intent(inout)  :: Econd
    !   OPTIONAL INPUT: background model parameter for linearization
    type(modelParam_t), optional, intent(in)  :: m0


    !   LOCAL Variables
    integer  :: imgrid, ix, iy, ic, icx, icy, izE, iz, izv
    integer  :: nx,ny,nz,nzAir, nzCum
    integer  :: ccoeff_current,mgridSize
    integer  :: status
    logical		:: linearizedMapping
    ! the prefix V are the volumes S is total over surrounding cells
    !   common sides are ommited from volume calculation
    real(kind=prec)  :: Vdr,Vdl,Vur,Vul,S
    real(kind=prec)  :: Vdn,Vds,Vun,Vus
    real(kind=prec)  :: Vnr,Vnl,Vsr,Vsl
    real(kind=prec)  ::airCond
    real(kind=prec),pointer,dimension(:,:,:)	:: temp ! x,y,z
    real(kind=prec),pointer,dimension(:,:,:)    :: v ! x,y,z

    linearizedMapping = present(m0)

    ! model parameters must be allocated
      if (.not. m%allocated) then
        write(0,*) 'input model parameter in ModelParamToEdgeMG not allocated yet'
        stop
      end if

    ! allocate Econd for multigrid
      if (.not. Econd%allocated) then
        call create(m%grid,Econd,EDGE) !create_rvector_mg
      end if

      allocate(v(m%grid%nx, m%grid%ny, 0:m%grid%nzEarth+1))
      allocate(temp(m%grid%nx, m%grid%ny, 0:m%grid%nzEarth+1))
      airCond = m%airCond
      nzCum = 0
      v(:,:,1:m%grid%nzEarth) = m%cellCond%v(:,:,:)
      v(:,:,0) = R_ZERO
      v(:,:,m%grid%nzEarth+1) = R_ZERO

      do imgrid = 1, m%grid%mgridSize ! Global loop over sub-grids
        temp = R_ZERO
        nx = m%grid%gridArray(imgrid)%nx ! nx of the current subgrid
        ny = m%grid%gridArray(imgrid)%ny ! ny
        nz = m%grid%gridArray(imgrid)%nz ! nz = nzEarth+nzAir
        nzAir = m%grid%gridArray(imgrid)%nzAir
        ccoeff_current = 2**m%grid%coarseness(imgrid)
        do iz = nzAir,nz+1
           if(nzAir == nz) cycle
           izE = iz - nzAir
           izv = izE + nzCum
          do iy = 1, ny
            do ix = 1, nx
              do icy = 1, ccoeff_current
                do icx = 1, ccoeff_current
                  temp(ix,iy,izE) = temp(ix,iy,izE) + v(ccoeff_current*(ix-1)+icx, ccoeff_current*(iy-1)+icy,izv)
                enddo ! icx
              enddo ! icy
              temp(ix,iy,izE) = temp(ix,iy,izE)/ccoeff_current**2
            enddo ! nx
          enddo ! ny
        enddo ! nz

        if(m%paramType .EQ. LOGE) then
        !  if(linearizedMapping) then   ! NEED TO BE CHECKED
        !    do iz = nzAir,nz+1
        !       izE = iz - nzAir
        !       izv = izE + nzCum
        !      do iy = 1, ny
        !        do ix = 1, nx
        !          do ic = 1, ccoeff_current
        !            temp(ix,iy,izE)= temp(ix,iy,izE)*exp(m0%cellCond%v(ccoeff_current*(ix-1)+1, ccoeff_current*(iy-1)+1,izv))
        !          enddo
        !        enddo
        !      enddo
        !    enddo
        !    airCond = R_ZERO
        !  else
            temp = exp(temp)
            airCond = exp(m%airCond)
        !  endif
        endif

        ! Ex edges
        Econd%rvArray(imgrid)%x = 0.0
        do ix =1 ,nx
          do iy = 2,ny
            do iz = nzAir+1, nz+1
               izE = iz-nzAir
               izv = izE + nzCum+m%grid%nzAir
              Vdr = m%grid%gridArray(imgrid)%dy(iy)*m%grid%dz(izv)
              Vdl = m%grid%gridArray(imgrid)%dy(iy-1)*m%grid%dz(izv)
              Vur = m%grid%gridArray(imgrid)%dy(iy)*m%grid%dz(izv-1)
              Vul = m%grid%gridArray(imgrid)%dy(iy-1)*m%grid%dz(izv-1)
              S = Vdr+Vdl+Vur+Vul
              if(izE.eq.1.and.nzAir.ne.0) then
                Econd%rvArray(imgrid)%x(ix,iy,iz) = (temp(ix,iy,izE)*Vdr + &
                                    temp(ix,iy-1,izE)*Vdl)/S
              else
                Econd%rvArray(imgrid)%x(ix,iy,iz) = (temp(ix,iy,izE)*Vdr + &
                                   temp(ix,iy-1,izE)*Vdl + &
                                   temp(ix,iy  ,izE-1)*Vur + &
                                   temp(ix,iy-1,izE-1)*Vul)/S
                endif
            enddo
            do iz = 1,nzAir
              Econd%rvArray(imgrid)%x(ix,iy,iz) = airCond
            enddo
          enddo
        enddo

        ! Ey edges
        Econd%rvArray(imgrid)%y = 0.0
        do ix = 2, nx
          do iy = 1, ny
             do iz = nzAir+1,nz+1
                izE = iz-nzAir
                izv = izE + nzCum + m%grid%nzAir
               Vdn = m%grid%gridArray(imgrid)%dx(ix)*m%grid%dz(izv)
               Vds = m%grid%gridArray(imgrid)%dx(ix-1)*m%grid%dz(izv)
               Vun = m%grid%gridArray(imgrid)%dx(ix)*m%grid%dz(izv-1)
               Vus = m%grid%gridArray(imgrid)%dx(ix-1)*m%grid%dz(izv-1)
               S = Vdn+Vds+Vun+Vus
               if(izE.eq.1.and. nzAir.ne.0) then
                  Econd%rvArray(imgrid)%y(ix,iy,iz) = (temp(ix,iy,izE)*Vdn + &
                                     temp(ix-1,iy,izE)*Vds)/S
               else
                  Econd%rvArray(imgrid)%y(ix,iy,iz) = (temp(ix,iy,izE)*Vdn + &
                                      temp(ix-1,iy,izE)*Vds + &
                                      temp(ix,iy  ,izE-1)*Vun + &
                                      temp(ix-1,iy,izE-1)*Vus)/S
               endif
             enddo
             do iz = 1,nzAir
                  Econd%rvArray(imgrid)%y(ix,iy,iz) = airCond
             enddo
           enddo
        enddo

         ! Ez edges
         Econd%rvArray(imgrid)%z = 0.0
         do ix = 2, nx
            do iy = 2, ny
               do iz = nzAir+1,nz
                  izE = iz-nzAir
                  Vnr = m%grid%gridArray(imgrid)%dx(ix)*m%grid%gridArray(imgrid)%dy(iy)
                  Vnl = m%grid%gridArray(imgrid)%dx(ix)*m%grid%gridArray(imgrid)%dy(iy-1)
                  Vsr = m%grid%gridArray(imgrid)%dx(ix-1)*m%grid%gridArray(imgrid)%dy(iy)
                  Vsl = m%grid%gridArray(imgrid)%dx(ix-1)*m%grid%gridArray(imgrid)%dy(iy-1)
                  S =  Vnr+Vnl+Vsr+Vsl
                  Econd%rvArray(imgrid)%z(ix,iy,iz) = (temp(ix,iy,izE)*Vnr + &
                                      temp(ix,iy-1,izE)*Vnl + &
                                      temp(ix-1,iy  ,izE)*Vsr + &
                                      temp(ix-1,iy-1,izE)*Vsl)/S
               enddo
              do iz = 1,nzAir
                  Econd%rvArray(imgrid)%z(ix,iy,iz) = airCond
               enddo
            enddo
          enddo
       nzCum = nzCum + m%grid%gridArray(imgrid)%nzEarth
      enddo ! Global loop over sub-grids

        Econd%rvArray(m%grid%mgridSize)%x(:,:,m%grid%gridArray(m%grid%mgridSize)%nz+1) = &
        Econd%rvArray(m%grid%mgridSize)%x(:,:,m%grid%gridArray(m%grid%mgridSize)%nz)
        Econd%rvArray(m%grid%mgridSize)%y(:,:,m%grid%gridArray(m%grid%mgridSize)%nz+1) = &
        Econd%rvArray(m%grid%mgridSize)%y(:,:,m%grid%gridArray(m%grid%mgridSize)%nz)


    deallocate(temp,stat=status)
    deallocate(v,stat=status)
    end subroutine ModelParamtoEdgeMG

! ***********************************************************************
  subroutine dEdgeToModelParam(Econd,m,m0)
  !  Maps from a real vector (defined on edges) to modelParam;
  !    the adjoint of linear mapping implemented by ModelParamToEdge
  !  Interface to EdgeToModelParam that makes the code more readable.

    implicit none
    !  INPUTS:  real vector defined on edges
    type(rvector_mg), intent(in)            :: Econd
    !  OUTPUTS: model parameter
    type(modelParam_t), intent(inout)         :: m
    !  INPUT background model parameter
    type(modelParam_t), intent(in)	    :: m0

    call EdgeToModelParamMG(Econd,m,m0)

    m%zeroValued = .false.

  end subroutine dEdgeToModelParam

    !**********************************************************************
    subroutine EdgeToModelParamMG(Econd,m,m0)
    !  Maps from a real vector (defined on edges of the MULTIGRID) to modelParam;
    !    the adjoint of linear mapping implemented by ModelParamToEdgeMG
    !
    implicit none
    !  INPUTS:  real vector defined on edges
    type(rvector_mg), intent(in)  :: Econd
    !  OUTPUTS: model parameter
    type(modelParam_t), intent(inout)  :: m
    !  INPUT (OPTIONAL) background model parameter,
    !         required if m%paramtype=LOGE
    type(modelParam_t), optional, intent(in)  :: m0

    !   LOCAL Variables
    integer  :: ix,iy,iz,ize,izv,imgrid,ic
    integer  :: nx,ny,nz,nzE,nzAir,nzCum, ccoeff_current
    integer  :: status
    ! the prefix V are the volumes S is total over surrounding cells
    !   common sides are ommited from volume calculation
    real(kind=prec)  :: Vdr,Vdl,Vur,Vul,S
    real(kind=prec)  :: Vdn,Vds,Vun,Vus
    real(kind=prec)  :: Vnr,Vnl,Vsr,Vsl
    real(kind=prec),pointer,dimension(:,:,:)  :: temp


      if ((.not.m%allocated).or.(.not.Econd%allocated)) then
        call errStop('m or Econd not allocated yet in EdgeToModelParamMG')
        stop
      end if

      if((m%paramtype.eq.LOGE).and. (.not.present(m0))) then
        call errStop('Background conductivity required for paramType LOGE in EdgeToModelParamMG')
      endif

      ! Could add code to check whether the bounds are the same ...

      allocate(temp(m%grid%nx,m%grid%ny,0:m%grid%nzEarth+1))
      ! Now map onto edges, x, y, then z ... first  zero output
      nzCum = 0
      do imgrid = 1,m%grid%mgridSize  ! Global loop over subgrids
        temp = R_ZERO
        nx = m%grid%gridArray(imgrid)%nx ! nx of the current subgrid
        ny = m%grid%gridArray(imgrid)%ny ! ny
        nz = m%grid%gridArray(imgrid)%nz ! nz = nzEarth+nzAir
        nzAir = m%grid%gridArray(imgrid)%nzAir
        ccoeff_current = 2**m%grid%coarseness(imgrid)
        if(nzAir /= nz) then ! Do nothing for Air
          ! Ex edges: (leave boundary edges set to 0)
          do ix = 1, nx
            do iy = 2, ny
              do iz = nzAir+1,nz+1
                 izE = iz-nzAir
                 izv = izE + nzCum+ m%grid%nzAir
                 Vdr = m%grid%gridArray(imgrid)%dy(iy)*m%grid%dz(izv)
                 Vdl = m%grid%gridArray(imgrid)%dy(iy-1)*m%grid%dz(izv)
                 Vur = m%grid%gridArray(imgrid)%dy(iy)*m%grid%dz(izv-1)
                 Vul = m%grid%gridArray(imgrid)%dy(iy-1)*m%grid%dz(izv-1)
                 S = Vdr+Vdl+Vur+Vul
                 !if(izE.gt.1) then
                    temp(ix,iy,izE-1) = temp(ix,iy,izE-1)+ Econd%rvArray(imgrid)%x(ix,iy,iz)*Vur/S
                    temp(ix,iy-1,izE-1) = temp(ix,iy-1,izE-1)+Econd%rvArray(imgrid)%x(ix,iy,iz)*Vul/S
                 !endif
                 temp(ix,iy,izE) = temp(ix,iy,izE)+ Econd%rvArray(imgrid)%x(ix,iy,iz)*Vdr/S
                 temp(ix,iy-1,izE) = temp(ix,iy-1,izE)+Econd%rvArray(imgrid)%x(ix,iy,iz)*Vdl/S
              enddo ! iz
            enddo ! iy
          enddo ! ix

          ! Ey edges
          do ix = 2, nx
            do iy = 1, ny
               do iz = nzAir+1,nz+1
                  izE = iz-nzAir
                  izv = izE + nzCum+ m%grid%nzAir
                  Vdn = m%grid%gridArray(imgrid)%dx(ix)*m%grid%dz(izv)
                  Vds = m%grid%gridArray(imgrid)%dx(ix-1)*m%grid%dz(izv)
                  Vun = m%grid%gridArray(imgrid)%dx(ix)*m%grid%dz(izv-1)
                  Vus = m%grid%gridArray(imgrid)%dx(ix-1)*m%grid%dz(izv-1)
                  S = Vdn+Vds+Vun+Vus
                  !if(izE.gt.1) then
                    temp(ix,iy,izE-1) = temp(ix,iy,izE-1)+ Econd%rvArray(imgrid)%y(ix,iy,iz)*Vun/S
                    temp(ix-1,iy,izE-1) = temp(ix-1,iy,izE-1)+Econd%rvArray(imgrid)%y(ix,iy,iz)*Vus/S
                  !endif
                  temp(ix,iy,izE) = temp(ix,iy,izE)+ Econd%rvArray(imgrid)%y(ix,iy,iz)*Vdn/S
                  temp(ix-1,iy,izE) = temp(ix-1,iy,izE)+Econd%rvArray(imgrid)%y(ix,iy,iz)*Vds/S
               enddo
             enddo
           enddo

           ! Ez edges
           do ix = 2, nx
            do iy = 2, ny
               do iz = nzAir+1,nz
                  izE = iz-nzAir
                  Vnr = m%grid%gridArray(imgrid)%dx(ix)*m%grid%gridArray(imgrid)%dy(iy)
                  Vnl = m%grid%gridArray(imgrid)%dx(ix)*m%grid%gridArray(imgrid)%dy(iy-1)
                  Vsr = m%grid%gridArray(imgrid)%dx(ix-1)*m%grid%gridArray(imgrid)%dy(iy)
                  Vsl = m%grid%gridArray(imgrid)%dx(ix-1)*m%grid%gridArray(imgrid)%dy(iy-1)
                  S =  Vnr+Vnl+Vsr+Vsl
                  temp(ix,iy,izE) = temp(ix,iy,izE)+ Econd%rvArray(imgrid)%z(ix,iy,iz)*Vnr/S
                  temp(ix,iy-1,izE) = temp(ix,iy-1,izE)+Econd%rvArray(imgrid)%z(ix,iy,iz)*Vnl/S
                  temp(ix-1,iy,izE) = temp(ix-1,iy,izE)+ Econd%rvArray(imgrid)%z(ix,iy,iz)*Vsr/S
                  temp(ix-1,iy-1,izE) = temp(ix-1,iy-1,izE)+Econd%rvArray(imgrid)%z(ix,iy,iz)*Vsl/S
               enddo
            enddo
           enddo

          !Store temp on subgrids to model parameters on the finnest grid, where m is defined
           do iz = 1,nz
              izv = iz+nzCum
             do iy =1 ,ny
               do ix = 1, nx
                do ic = 1, ccoeff_current
                  if(m%paramType .eq. LOGE) then
                    m%cellCond%v(ccoeff_current*(ix-1)+ic,ccoeff_current*(iy-1)+ic,izv) = temp(ix,iy,iz) &
                    *exp(m0%cellCond%v(ccoeff_current*(ix-1)+ic,ccoeff_current*(iy-1)+ic,izv))
                  else
                    m%cellCond%v(ccoeff_current*(ix-1)+ic,ccoeff_current*(iy-1)+ic,izv) = temp(ix,iy,iz)
                  endif
                enddo ! ic
               enddo !ix
              enddo ! iy
            enddo !iz
        endif ! Do nothing for Air
        nzCum = nzCum + m%grid%gridArray(imgrid)%nzEarth

      enddo ! Global loop over subgrids

        m%zeroValued = .false.

        deallocate(temp, stat=status)
      end subroutine EdgeToModelParamMG

  ! *********************************************************************
  function ModelParamToOneEdgeMG(m,xyz,is,ix,iy,iz) result(r)

  ! This function is needed for the interpolation in Multigrid case
  !______________________________________________________________________
  !   computes conductivity for edge xyz/i,j,k using input modelParam
  !    structure m.  This function defines how the abstract
  !     conductivity parameter is mapped to edge conductivities needed
  !     for more accurate electric field interpolation.  The derivative
  !     of this function is required for evaluation of linearized
  !     data functionals, and for construction of the direct
  !     parameter space part of the comb
  !______________________________________________________________________

    type (modelParam_t), intent(in)     :: m
    integer, intent(in)                 :: xyz,ix,iy,iz,is
    real(kind=prec)                     :: r

    ! local variables
    integer                             :: izE,nx,ny,nz,nZE,nZa
    real(kind=prec)                     :: w11,w21,w12,w22,S,temp(2,2)

    !local
    integer :: ccoeff_current

    nx = m%grid%gridArray(is)%Nx
    ny = m%grid%gridArray(is)%Ny
    nz = m%grid%gridArray(is)%nZ
    nZE = m%grid%gridArray(is)%nZEarth
    nZa = Nz-nZE

    izE = iz-Nza
    ccoeff_current = 2**m%grid%coarseness(is)
    selectcase(xyz)
    case(1)      ! Ex edge
       w11 = m%grid%dy(iy)*m%grid%dz(iz)
       temp(1,1) = m%cellCond%v(ccoeff_current*(ix-1)+1, ccoeff_current*(iy-1)+1,izE)
       w21 = m%grid%dy(iy-1)*m%grid%dz(iz)
       temp(2,1) = m%cellCond%v(ccoeff_current*(ix-1)+1, ccoeff_current*(iy-2)+1,izE)
       w12 = m%grid%dy(iy)*m%grid%dz(iz-1)
       temp(1,2) = m%cellCond%v(ccoeff_current*(ix-1)+1, ccoeff_current*(iy-1)+1,izE-1)
       w22 = m%grid%dy(iy-1)*m%grid%dz(iz-1)
       temp(2,2) = m%cellCond%v(ccoeff_current*(ix-1)+1, ccoeff_current*(iy-2)+1,izE-1)
       S = w11+w21+w12+w22
       if(m%paramtype .eq. LOGE) then
          temp = exp(temp)
       endif
       if(izE.eq.1) then
          r = (temp(1,1)*w11+temp(1,2)*w21)/S
       else
          r = (temp(1,1)*w11 + temp(2,1)*w21+    &
               temp(1,2)*w12 + temp(2,2)*w22)/S
       endif

    case(2)      ! Ey edge
       w11 = m%grid%dx(ix)*m%grid%dz(iz)
       temp(1,1) = m%cellCond%v(ccoeff_current*(ix-1)+1, ccoeff_current*(iy-1)+1,izE)
       w21 = m%grid%dx(ix-1)*m%grid%dz(iz)
       temp(2,1) = m%cellCond%v(ccoeff_current*(ix-2)+1, ccoeff_current*(iy-1)+1,izE)
       w12 = m%grid%dx(ix)*m%grid%dz(iz-1)
       temp(1,2) = m%cellCond%v(ccoeff_current*(ix-1)+1, ccoeff_current*(iy-1)+1,izE-1)
       w22 = m%grid%dx(ix-1)*m%grid%dz(iz-1)
       temp(2,2) = m%cellCond%v(ccoeff_current*(ix-2)+1, ccoeff_current*(iy-1)+1,izE-1)
       S = w11+w21+w12+w22
       if(m%paramtype .eq. LOGE) then
          temp = exp(temp)
       endif
       if(izE.eq.1) then
          r = (temp(1,1)*w11+temp(1,2)*w21)/S
       else
          r = (temp(1,1)*w11 + temp(2,1)*w21+ &
               temp(1,2)*w12 + temp(2,2)*w22)/S
       endif

    case(3)      ! Ez edge
       w11 = m%grid%dx(ix)*m%grid%dy(iy)
       temp(1,1) = m%cellCond%v(ccoeff_current*(ix-1)+1, ccoeff_current*(iy-1)+1,izE)
       w21 = m%grid%dx(ix)*m%grid%dy(iy-1)
       temp(2,1) = m%cellCond%v(ccoeff_current*(ix-1)+1, ccoeff_current*(iy-2)+1,izE)
       w12 = m%grid%dx(ix-1)*m%grid%dy(iy)
       temp(1,2) = m%cellCond%v(ccoeff_current*(ix-2)+1, ccoeff_current*(iy-1)+1,izE)
       w22 = m%grid%dx(ix-1)*m%grid%dy(iy-1)
       temp(2,2) = m%cellCond%v(ccoeff_current*(ix-2)+1, ccoeff_current*(iy-2)+1,izE)
       S = w11+w21+w12+w22
       if(m%paramtype .eq. LOGE) then
          temp = exp(temp)
       endif
       r = (temp(1,1)*w11 + temp(2,1)*w21+  &
            temp(1,2)*w12 + temp(2,2)*w22)/S

    endselect
  end function ModelParamToOneEdgeMG

