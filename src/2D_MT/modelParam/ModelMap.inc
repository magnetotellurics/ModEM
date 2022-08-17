!!!!!!  >>>>>>>>>>>>>>>>>>>>   TE mode routines <<<<<<<<<<<<<<<<<<
! *****************************************************************************
      subroutine dModelParamToNode(CellCond,NodeCond,Sigma0)
      !  LINEARIZED mapping from conductivity or log conductivity,
      !      defined on Earth cells to interior Earth node conductivities,
      !      averaging over adjacent cells.  Overwrites inputs on
      !      these cells, leaving  other nodes (in air,
      !      on grid boundaries) unaltered
      !   Sigma0 is background conductivity for linearization
      !        -- this is an optional argument,
      !           required when CellCond%paramType == 'LOGE'
      !
      !   NOTE: CellCond array is real, while NodeCond array is
      !       complex  !!!
      !
      type (modelParam_t), intent(in)   :: CellCond
      type (modelParam_t), intent(in), optional   :: Sigma0

      !  OUTPUT is of type vec2D, of gridType NODE
      !   i.e., output should be defined on all grid nodes, though
      !   this routine only modifies interior nodes within the
      !   Earth.
      type (cvector), intent(inout)   :: NodeCond

      ! local variables
      integer iy,iz,Ny,Nz,Nza,izc, NzE
      real(kind=prec) w00,w01,w10,w11,wsum
      real(kind=prec), allocatable, dimension(:,:)	:: temp

      Ny = cellCond%grid%Ny
      Nz = cellCond%grid%Nz
      Nza = cellCond%grid%Nza
      NzE = Nz-Nza
      allocate(temp(Ny,NzE))
      if(NodeCond%gridType .NE. NODE) then
         call errStop('NodeCond must be of gridType NODE in CellToNode')
      endif
      if(CellCond%paramType .EQ. LOGE) then
         if(present(Sigma0)) then
            temp = CellCond%v*exp(Sigma0%v)
         else
            call errStop('Sigma0 input to CellToNode is required for case LOGE')
         endif
      else
         temp = CellCond%v
      endif

      !  loop over interior Earth nodes
      !        (including air/earth interface)
      do iz = Nza+1,Nz
         izc = iz-Nza
         do iy = 2,Ny
            w00 = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz-1)
            w10 = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz-1)
            w01 = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz)
            w11 = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz)
            wsum = w00+w10+w01+w11
            w00 = w00/wsum
            w10 = w10/wsum
            w01 = w01/wsum
            w11 = w11/wsum
            if(izc .EQ. 1) then
               NodeCond%v(iy,iz) = w01*temp(iy-1,izc) + &
			      w11*temp(iy,izc)
            else
               NodeCond%v(iy,iz) = w00*temp(iy-1,izc-1)+ &
			      w10*temp(iy,izc-1) + &
			      w01*temp(iy-1,izc) + &
			      w11*temp(iy,izc)
            endif
         enddo
      enddo
      deallocate(temp)
      end subroutine dModelParamToNode

   ! *****************************************************************************
      subroutine dNodeToModelParam(NodeCond,CellCond,Sigma0)
      !  Maps from interior Earth nodes to Earth cells;
      !   transpose (adjoint) of EarthCellToNode
      !  INPUT is of type vec2D, of gridType NODE
      !   NodeCond is assumed defined on all nodes, but only
      !   interior Earth nodes are used
      type (cvector), intent(in)   :: NodeCond


      !  OUTPUT is linear cell conductivity parameter structure
      !    (paramType should be LINEAR)
      !  CellCond is overwritten by call to this routine
      type (modelParam_t), intent(inout)   :: CellCond
      !   Optional input Sigma0 is background conductivity, required
      !        for log conductivity parameterization
      type (modelParam_t), intent(in), optional   :: Sigma0

      !   NOTE: CellCond array is real, while NodeCond array is
      !       complex  !!!

      ! local variables
      integer iy,iz,izc,Ny,Nz,Nza
      real(kind=prec) w00,w01,w10,w11,wsum

      if(NodeCond%gridType .NE. NODE) then
         call errStop('NodeCond must be of gridType NODE in NodeToCell')
      endif
      if((CellCond%paramType .eq. LOGE) .and. (.not.present(Sigma0))) then
         call errStop('Background conductivity required for paramType LOGE in NodeToCell')
      endif

      Ny = cellCond%grid%Ny
      Nz = cellCond%grid%Nz
      Nza = cellCond%grid%Nza

      !   zero cell conductivity values
      CellCond%v = R_ZERO

      !  loop over interior Earth nodes
      !        (including air/earth interface)
      do iz = Nza+1,Nz
         izc = iz-Nza
         do iy = 2,Ny
            w00 = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz-1)
            w10 = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz-1)
            w01 = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz)
            w11 = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz)
            wsum = w00+w10+w01+w11
            w00 = w00/wsum
            w10 = w10/wsum
            w01 = w01/wsum
            w11 = w11/wsum
            if(izc .GT. 1) then
               CellCond%v(iy-1,izc-1) = &
		 CellCond%v(iy-1,izc-1)+w00*real(NodeCond%v(iy,iz))
               CellCond%v(iy,izc-1) = &
		 CellCond%v(iy,izc-1)+w10*real(NodeCond%v(iy,iz))
            endif
            CellCond%v(iy-1,izc) = &
		 CellCond%v(iy-1,izc)+w01*real(NodeCond%v(iy,iz))
            CellCond%v(iy,izc) = &
		 CellCond%v(iy,izc)+w11*real(NodeCond%v(iy,iz))
         enddo
      enddo
      if(CellCond%paramType .eq. LOGE) then
         CellCond%v = CellCond%v*exp(Sigma0%v)
      endif

      CellCond%zeroValued = .false.


      end subroutine dNodeToModelParam

      !**********************************************************************
      subroutine ModelParamToCell(CellCond,Ny,Nz,Cond)
      !  Copies from formal Earth cell conductivity structure
      !    to a standard real array defined for the whole grid.
      !     Sets air cells to AirCond
      !
      !  INPUT is cell conductivity parameter structure
      !    (i.e., paramType should be LINEAR)
      type (modelParam_t), intent(in)   :: CellCond
      integer, intent(in)	:: Ny,Nz

      real(kind=prec), intent(inout)   :: Cond(Ny,Nz)

      ! local variables
      integer iy,iz,izc,NyC,NzC,Nza

      if (.not. cellCond%allocated) then
        call errStop('CellCond not allocated on input to ModelParamToCell')
      endif

      NyC = cellCond%grid%Ny
      NzC = cellCond%grid%Nz
      Nza = cellCond%grid%Nza

      if((NyC .NE. Ny).or.(NzC .NE. Nz)) then
        call errStop('dimensions of CellCond, Cond disagree in ModelParamToCell')
      endif

      Cond = CellCond%AirCond
      do iz = Nza+1,Nz
         izc = iz-Nza
         do iy = 1,Ny
            Cond(iy,iz) = CellCond%v(iy,izc)
         enddo
      enddo
      if(CellCond%paramType .eq. LOGE) then
         Cond = exp(Cond)
      endif

      end subroutine ModelParamToCell

!!!!!!  >>>>>>>>>>>>>>>>>>>>   TM mode routines <<<<<<<<<<<<<<<<<<
      !**********************************************************************
      !   computes resistivity for cell j,k using input parameter
      !     structure sigma.  This function defines how the abstract
      !     conductivity parameter is mapped to cell resistivites needed
      !     for TM electric field interpolation functions.  The derivative
      !     of this function is required for evaluation of linearized
      !     data functionals, and for construction of the direct
      !     parameter space part of the comb
      function ModelParamToOneCell(sigma,j,k) result(r)
      type (modelParam_t), intent(in)	:: sigma
      integer, intent(in)		:: j,k
      real(kind=prec)	                :: r

      ! local variables
      if(sigma%paramType .eq. LOGE) then
         r = exp(-sigma%v(j,k))
      else
         r = 1./sigma%v(j,k)
      endif

      end function ModelParamToOneCell

   !**********************************************************************
   subroutine dModelParamToEdge(CellCond, EdgeCond_y,EdgeCond_z, sigma0)
   ! map from Cell to Edge on Y and Z axis

   type (modelParam_t), intent(in)    :: CellCond
   type (modelParam_t), intent(in)    :: sigma0
   type (cvector), intent(inout)   :: EdgeCond_y, EdgeCond_z

   real (kind=prec), allocatable, dimension(:,:) :: temp


   integer :: iz, iy, Nz, Nza, Ny, Nzb
   real (kind=prec)    ::  wll, wrr, wdd, wuu, wsum

   Ny  = CellCond%grid%Ny
   Nz  = CellCond%grid%Nz
   Nza = CellCond%grid%Nza
   Nzb = Nz - Nza

   allocate(temp(Ny,Nzb))
   if(CellCond%paramType .EQ. LOGE) then
      temp = - CellCond%v/exp(sigma0%v)
   else
      temp = - CellCond%v/(sigma0%v*sigma0%v)
   endif
   do iy = 2, Ny
     do iz = 1, Nzb
       wll = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz+Nza)
       wrr = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz+Nza)
       wsum = wll + wrr
       wll = wll/wsum
       wrr = wrr/wsum
       EdgeCond_y%v(iy,iz) = wll*(temp(iy-1,iz)) + &
                             wrr*(temp(iy,iz))
     enddo
   enddo

   do iz = 2, Nzb
     do iy = 1, Ny
       wuu = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz-1+Nza)
       wdd = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz+Nza)
       wsum = wuu + wdd
       wuu = wuu/wsum
       wdd = wdd/wsum
       EdgeCond_z%v(iy,iz) = wuu*(temp(iy,iz-1)) + &
                             wdd*(temp(iy,iz))
     enddo
   enddo

   deallocate(temp)

   end subroutine dModelParamToEdge

   !**********************************************************************

   subroutine dEdgeToModelParam(EdgeCond_y,EdgeCond_z, CellCond, sigma0)
   ! map from Edge to Cell on Y and Z axis

   type (cvector), intent(in)   :: EdgeCond_y, EdgeCond_z
   type (modelParam_t), intent(in) :: sigma0
   type (modelParam_t), intent(inout)    :: CellCond


   integer :: iz, iy, Nz, Nza, Ny, Nzb
   real (kind=prec)    ::  wll, wrr, wdd, wuu, wsum

   Ny  = CellCond%grid%Ny
   Nz  = CellCond%grid%Nz
   Nza = CellCond%grid%Nza
   Nzb = Nz - Nza

   CellCond%v = 0.
   do iy = 2, Ny
     do iz = 1, Nzb
       wll = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz+Nza)
       wrr = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz+Nza)
       wsum = wll + wrr
       wll = wll/wsum
       wrr = wrr/wsum
       CellCond%v(iy-1,iz) = CellCond%v(iy-1,iz) + wll*real(EdgeCond_y%v(iy,iz))
       CellCond%v(iy,iz)   = CellCond%v(iy,iz)   + wrr*real(EdgeCond_y%v(iy,iz))
     enddo
   enddo

   do iz = 2, Nzb
     do iy = 1, Ny
       wuu = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz-1+Nza)
       wdd = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz+Nza)
       wsum = wuu + wdd
       wuu = wuu/wsum
       wdd = wdd/wsum
       CellCond%v(iy,iz-1) = CellCond%v(iy,iz-1) + wuu*real(EdgeCond_z%v(iy,iz))
       CellCond%v(iy,iz)   = CellCond%v(iy,iz)   + wdd*real(EdgeCond_z%v(iy,iz))
     enddo
   enddo

   if (CellCond%paramType .eq. LOGE) then
      CellCond%v = - CellCond%v/exp(sigma0%v)
   else
      CellCond%v = - CellCond%v/(sigma0%v*sigma0%v)
   endif

   CellCond%zeroValued = .false.

   end subroutine dEdgeToModelParam

   !**********************************************************************
   subroutine SparseCellToModelParam(Qj,sigma0,dsigmaReal,dSigmaImag)

   !  given input sparse vector defined on cells,
   !   compute vector defined on parameter space
   !   output result as dsigmaReal, dsigmaImag
   !   i.e. multiplies a sparse vector by (d\pi/dm)^T
   !   used to compute and multiply by Q^T

   type (sparsevecc), intent(in)   :: Qj
   type (modelParam_t), intent(in) :: sigma0
   type (modelParam_t), intent(out)    :: dsigmaReal
   type (modelParam_t), intent(out),optional    :: dsigmaImag

   ! local variables
   integer		:: jk,j,k

   dsigmaReal = sigma0
   call zero(dsigmaReal)

   if(present(dSigmaImag)) then
      dsigmaImag = sigma0
      call zero(dsigmaImag)
   endif

   do jk = 1,Qj%nCoeff
      j = Qj%J(jk)
      k = Qj%K(jk)
      if(sigma0%paramType .eq. LOGE) then
         dsigmaReal%v(j,k) = dSigmaReal%v(j,k) &
			-real(Qj%C(jk))/exp(sigma0%v(j,k))
         if(present(dSigmaImag)) then
            dsigmaImag%v(j,k) = dSigmaImag%v(j,k) &
			-imag(Qj%C(jk))/exp(sigma0%v(j,k))
         endif
      else
         dsigmaReal%v(j,k) = dSigmaReal%v(j,k) &
			-real(Qj%C(jk))/(sigma0%v(j,k)*sigma0%v(j,k))
         if(present(dSigmaImag)) then
            dsigmaImag%v(j,k) = dSigmaImag%v(j,k) &
		-imag(Qj%C(jk))/(sigma0%v(j,k)*sigma0%v(j,k))
         endif
      endif
   enddo

   dsigmaReal%zeroValued = .false.
   if(present(dsigmaImag)) then
      dsigmaImag%zeroValued = .false.
   endif

   end subroutine SparseCellToModelParam
