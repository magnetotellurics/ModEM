!**
! Basic model parameter class --
! cell conductivities defined on numerical grid
! in Earth cells only -- as coded this is specific
! to usual CSG.
!*
module ModelParameterCell_SG
   !
   use Constants
   use Grid
   use Grid3D_SG
   use rScalar
   use rVector
   use rScalar3D_SG
   use rVector3D_SG 
   use ModelParameter
   use Grid1D
   use Grid2D
   use ModelParameter1D
   use ModelParameter2D
   !
   type, extends( ModelParameter_t ) :: ModelParameterCell_SG_t
       
       ! This will be slighty modified from numerical grid --
       ! NzAir = 0 -- more generally this might be a
       ! completely different grid.
       type( Grid3D_SG_t ) :: paramGrid
       !
       ! Pointer to the original grid
       class( Grid3D_SG_t ), pointer :: grid
       
       type( rScalar3D_SG_t ) :: cellCond
       type( rVector3D_SG_t ) :: V_Ex4    ! edge volume times 4 -- also provides an
                              ! instance of an edge vector.
       type( rScalar3D_SG_t ) :: V_C       ! cell volume -- also an

    contains
	   !
       final :: ModelParameterCell_SG_dtor
       !
       procedure, public :: Length
       
       procedure, public :: Zeros
       procedure, public :: CopyFrom
    
       procedure, public :: EdgeLength
       procedure, public :: EdgeArea
       
       ! Model mapping methods
       procedure, public :: PDEmapping
       procedure, public :: dPDEmapping
       procedure, public :: dPDEmappingT
       !
       procedure, public :: Slice1D => Slice1DModelParameterCell
       procedure, public :: Slice2D => Slice2DModelParameterCell
	   !
	   procedure, public :: AvgModel1D => AvgModel1DModelParameterCell
       !
       procedure, public :: SetType
       
   end type ModelParameterCell_SG_t

   interface ModelParameterCell_SG_t
       module procedure ModelParameterCell_SG_ctor
   end interface ModelParameterCell_SG_t
   
contains
   
   !**
   !
   !*
   function ModelParameterCell_SG_ctor( grid, ccond, p_paramType ) result( self )
      class( Grid3D_SG_t ), target, intent( in ) :: grid
      class( rScalar3D_SG_t ), intent( in )      :: ccond
      character(*), optional, intent(in)         :: p_paramType
	  !
	  type( ModelParameterCell_SG_t ) :: self
      !
      integer :: nx, ny, nz, nzAir
      character(30) :: paramType
      !
      write(*,*) "Constructor ModelParameterCell_SG_t"
      !
      if (.not.present(p_paramType)) then
         paramType = LOGE;
      else
         paramType = trim(p_paramType)
      end if
      !
      ! Point to the original grid
      self%grid => grid
      !
      nx = grid%nx
      ny = grid%ny
      nz = grid%nz - grid%nzAir

      nzAir = 0
      self%ParamGrid = Grid3D_SG_t(nx, ny, nzAir, nz, &
        grid%dx, grid%dy, &
        grid%dz( grid%nzAir+1:grid%nz ) )
      
      self%cellCond = rScalar3D_SG_t( self%paramGrid, CELL_EARTH )
      self%cellCond = ccond
      
      if ( present( p_paramType ) ) then
     call self%SetSigMap(p_paramType)
     ! We initially specify airCond as linear conductivity!
     self%AirCond = self%SigMap(self%airCond, 'inverse')
      end if
      
      self%isAllocated = .true.
      
   end function ModelParameterCell_SG_ctor
   !
   ! ModelParameterCell_SG_ destructor
   subroutine ModelParameterCell_SG_dtor( self )
      !
      implicit none
      !
      type( ModelParameterCell_SG_t ), intent( inout ) :: self
      !
      write(*,*) "Destructor ModelParameterCell_SG_t"
      !
	  !call self%deallocate()
	  !
   end subroutine ModelParameterCell_SG_dtor
   !
  !    **********
  !
  function Slice1DModelParameterCell( self, ix, iy ) result( m1D )
      !   extracts slice corresponding to column j of model parameter
      implicit none
      ! Arguments
      class(ModelParameterCell_SG_t), intent(in) :: self
	  !
      integer, intent(in)      :: ix, iy
      type(ModelParameter1D_t) ::  m1D 
      !   local variables
      type(Grid1D_t) :: grid1
      real(kind=prec), allocatable, dimension(:) :: CondSlice

      !   extract 1D grid
      grid1 = self%grid%Slice1D()
      !   create 1D model parameter
      m1D = ModelParameter1D_t(grid1)
      !   comnductivity slice
      allocate(CondSlice(grid1%nz))
      !
      CondSlice = self%CellCond%v( ix, iy, : )
      !
      call m1d%SetConductivity(CondSlice, self%AirCond, &
      self%paramType, self%mKey)
      !
   end function Slice1DModelParameterCell
   !
   function AvgModel1DModelParameterCell( self ) result( m1D )
      !   extracts slice corresponding to column j of model parameter
      implicit none
      ! Arguments
      class(ModelParameterCell_SG_t), intent(in) :: self
	  !
      type(ModelParameter1D_t) ::  m1D 
      !   local variables
      type(Grid1D_t) :: grid1
      real(kind=prec), allocatable, dimension(:) :: CondSlice
	  real(kind=prec) :: wt, temp_sigma_value
	  integer :: i, j, k
      !   extract 1D grid
      grid1 = self%grid%Slice1D()
      !   create 1D model parameter
      m1D = ModelParameter1D_t( grid1 )
      !   comnductivity slice
      allocate( CondSlice( grid1%nzEarth ) )
	  !
	  do k = 1, self%grid%nzEarth
		wt = R_ZERO
		temp_sigma_value = R_ZERO
		do i = 1, self%grid%Nx
		   do j = 1, self%grid%Ny
		      wt = wt + self%grid%dx(i) * self%grid%dy(j)
			  temp_sigma_value = temp_sigma_value + self%CellCond%v( i, j, k ) * &
			  self%grid%dx(i) * self%grid%dy(j)
		   end do
		end do
		!
		CondSlice( k ) = exp( temp_sigma_value / wt )
		!
		!write(*,*) k, self%CellCond%v( 1, 1, k ), CondSlice( k ), 1.0 / CondSlice( k )
	  end do
      !
      call m1d%SetConductivity( CondSlice, self%AirCond, &
      self%paramType, self%mKey )
      !
   end function AvgModel1DModelParameterCell
   !
   function Slice2DModelParameterCell( self, axis, j ) result( m2D )
     !   extracts slice corresponding to column j of model parameter
     implicit none
     ! Arguments
     class( ModelParameterCell_SG_t ), intent( in ) :: self
     integer, intent( in )                          :: axis, j
     !
     type( ModelParameter2D_t ) ::  m2D 
     !   local variables
     character(:), allocatable :: paramType
     type(Grid2D_t)   :: grid2
     real(kind=prec), allocatable, dimension(:,:) :: CondSlice
     !
     paramType = LINEAR
     !   extract 2D grid
     grid2 = self%grid%Slice2D()
     !   create 2D model parameter
     m2D = ModelParameter2D_t( grid2 )
     !   comnductivity slice
     allocate( CondSlice( grid2%ny, grid2%nzEarth ) )
     !
     if( axis == 1 ) then
       CondSlice = self%CellCond%v(j,:,:)
     else if( axis == 2 ) then
        CondSlice = self%CellCond%v(:,j,:)
     else if( axis == 3 ) then
        CondSlice = self%CellCond%v(:,:,j)
     else
        stop "ModelParameter:Slice2D: wrong axis"
     endif
     !
     call m2D%SetConductivity(CondSlice, self%AirCond, &
       paramType, self%mKey)
   !
   end function Slice2DModelParameterCell
   !
   function Length(self) result(nParam)
      implicit none
      class(ModelParameterCell_SG_t), intent(in) :: self
      integer :: nParam
      
      nParam = self%ParamGrid%Length()
   end function Length
   
   !**
   ! Zeros
   ! Zero model parameter
   !*
   subroutine Zeros(self)
      implicit none
      ! Arguments
      class(ModelParameterCell_SG_t), intent(inout) :: self
      
      call self%cellCond%Zeros()
   end subroutine Zeros

   !**
   ! Copy rhs to self.
   !*
   subroutine CopyFrom(self, rhs)
      implicit none
      ! Arguments
      class(ModelParameterCell_SG_t), intent(inout) :: self
      class(ModelParameter_t)     , intent(in)      :: rhs
      
      select type(rhs)
      class is(ModelParameterCell_SG_t)
     self%ParamGrid = rhs%ParamGrid
     self%cellCond   = rhs%cellCond
     self%airCond    = rhs%airCond
     self%paramType = rhs%paramType
     self%mKey         = rhs%mKey
      class default
     write(*, *) 'ERROR:ModelParameterCell:CopyFrom'
     write(*, *) '         Incompatible input. Exiting.'
     STOP
      end select
   end subroutine CopyFrom
   !
   function EdgeLength( self ) result ( length )

         class(ModelParameterCell_SG_t), intent( in ) :: self
         type( rVector3D_SG_t )                   :: length
      !
         ! local variables
         integer :: ix, iy, iz

         length = rVector3D_SG_t( self%grid, EDGE )

         ! x-component edge length elements
         do ix = 1, self%grid%nx
        length%x(ix, :, :) = self%grid%dx(ix)
         enddo

         ! y-component edge length elements
         do iy = 1, self%grid%ny
        length%y(:, iy, :) = self%grid%dy(iy)
         enddo

         ! z-component edge length elements
         do iz = 1, self%grid%nz
        length%z(:, :, iz) = self%grid%dz(iz)
         enddo

   end function EdgeLength
   !
   function EdgeArea( self ) result ( area )

         class(ModelParameterCell_SG_t), intent( in ) :: self
         type( rVector3D_SG_t )                   :: area
      ! local variables
      integer                  :: ix,iy,iz

      area = rVector3D_SG_t( self%grid, EDGE )

      ! edge areas are made for all the edges
      ! for x-components
      do ix = 1, self%grid%nx
       do iy = 1, self%grid%ny+1
         do iz = 1, self%grid%nz+1

        ! area%x values are centered within dx.
        area%x(ix, iy, iz) = self%grid%delY(iy)*self%grid%delZ(iz)

         enddo
       enddo
      enddo

      ! edge areas are made for all the edges
      ! for y-components
      do ix = 1, self%grid%nx+1
       do iy = 1, self%grid%ny
         do iz = 1, self%grid%nz+1

        ! area%y values are centered within dy.
        area%y(ix, iy, iz) = self%grid%delX(ix) * self%grid%delZ(iz)

         enddo
       enddo
      enddo

      ! edge areas are made for all the edges
      ! for z-components
      do ix = 1, self%grid%nx+1
       do iy = 1, self%grid%ny+1
         do iz = 1, self%grid%nz

        ! area%z values are centered within dz.
        area%z(ix, iy, iz) = self%grid%delX(ix) * self%grid%delY(iy)

         enddo
       enddo
      enddo

   end function EdgeArea

   function PDEmapping( self ) result( eVec )
      implicit none
      ! Arguments
      class(ModelParameterCell_SG_t), intent(inout)   :: self
      ! Local variables
      type(rScalar3D_SG_t) :: SigmaCell
      integer :: i, j, k, k0, k1, k2
      class(rVector_t), allocatable :: eVec
      !
      type( rVector3D_SG_t ) :: length, area
      
      allocate( eVec, source = rVector3D_SG_t( self%grid, EDGE ) )
      
      SigmaCell = rScalar3D_SG_t( self%grid, CELL )
      
      k0 =   self%ParamGrid%nzAir
      k1 = k0 + 1
      k2 = self%ParamGrid%Nz
      SigmaCell%v(:, :, 1:k0) = self%airCond
      
      ! Note: AirCond should always be in linear domain, but conductivity
      ! in cells is generally transformed -- SigMap converts to linear
      SigmaCell%v(:, :, k1:k2) = self%SigMap(self%cellCond%v)
      !
      self%V_C = rScalar3D_SG_t( self%grid, CENTER )
      !
      do i = 1, self%grid%nx
       do j = 1, self%grid%ny
       do k = 1, self%grid%nz
           self%V_C%v(i, j, k) = self%grid%dx(i)*self%grid%dy(j)*self%grid%dz(k)
       enddo
       enddo
      enddo
      !
      ! Form Conductivity--cell volume product
      sigmaCell = sigmaCell * self%V_C
      ! Sum onto edges
      call eVec%SumCells( SigmaCell )
      !
      length = self%EdgeLength()
      !
      area = self%EdgeArea()
      !
      self%V_Ex4 = length%diagMult( area )
      !
      ! Divide by total volume -- sum of 4 cells
      ! surrounding edge -- just 4*V_E      
      eVec = eVec / self%V_Ex4
      
   end function PDEmapping
   
   !**
   ! PDE mapping linearized at background model
   ! parameter m0, applied to dm result is an edge-vector eVec.
   !*
   function dPDEmapping(self, dm) result(eVec)
      implicit none
      class(ModelParameterCell_SG_t), intent(in)   :: self
      class(ModelParameter_t)     , intent(in)   :: dm
      ! Local variables
      type(rScalar3D_SG_t) :: SigmaCell
      character(len = 5), parameter :: JOB = 'DERIV'
      integer :: k0, k1, k2
      class(rVector_t), allocatable :: eVec

      select type(dm)
      class is(modelParameterCell_SG_t)
     allocate(eVec, source = rVector3D_SG_t(self%paramGrid, EDGE))
     SigmaCell = rScalar3D_SG_t(self%ParamGrid, CELL)
     
     ! Set Earth cells using m0, SigMap and dm
     ! I am doing this explicitly -- could make SigmaCell on ParamGrid
     ! then move Earth part to a Vector on ModelGrid (this is how we
     ! would do this more generally, when model space was really different
     ! from modeling grid.
     k0 = self%ParamGrid%NzAir
     k1 = k0 + 1
     k2 = self%ParamGrid%Nz
     
     SigmaCell%v(:,:,k1:k2) = self%SigMap(self%cellCond%v, JOB)
     
     SigmaCell%v(:,:,k1:k2) = SigmaCell%v(:,:,k1:k2)*dm%cellCond%v
     
     ! Average onto edges, as in PDEmapping ...
     sigmaCell = sigmaCell * self%V_C
     call eVec%SumCells(SigmaCell)
     eVec = eVec / self%V_Ex4
     
      class default
     write(*, *) 'ERROR:ModelParameterCell_SG:dPDEmapping:'
     write(*, *) '         Incompatible input [dm]. Exiting.'
     STOP
      end select
   end function dPDEmapping
   
   !**
   ! Transpose (adjoint) of dPDEmapping, applied to an edge-vector eVec
   ! result is a model parameter dm.
   !*
   function dPDEmappingT(self, eVec) result(dm)
      implicit none
      class(ModelParameterCell_SG_t), intent(in)      :: self
      class(rVector_t)           , intent(in)      :: eVec
      ! Local variables
      class(rScalar_t), allocatable :: sigmaCell
      type(rVector3D_SG_t) :: vTemp
      character(len = 5), parameter :: JOB = 'DERIV'
      integer :: k0, k1, k2
      class(Grid3D_SG_t), pointer :: grid
      class(ModelParameter_t), allocatable :: dm

      select type(eVec)
      class is(rVector3D_SG_t)
     allocate(dm, source = ModelParameterCell_SG_t(self%paramGrid, self%cellCond, self%paramType))

     select type(dm)
     class is(ModelParameterCell_SG_t)
          ! Create local temporary scalar and vector
          vTemp = rVector3D_SG_t(grid, EDGE)
          
          vTemp = eVec / self%V_Ex4
          call vTemp%SetAllBoundary(R_ZERO) 

          sigmaCell = vTemp%SumEdges()
          sigmaCell = sigmaCell * self%V_C
          dm%cellCond%v = self%SigMap(self%cellCond%v) 

          k0 = self%ParamGrid%NzAir
          k1 = k0 + 1
          k2 = self%ParamGrid%Nz
          
          select type(sigmaCell)
          class is(rScalar3D_SG_t)
              dm%cellCond%v = dm%cellCond%v * SigmaCell%v(:,:,k1:k2)
          end select
     end select
      class default
     write(*, *) 'ERROR:ModelParameterCell_SG:dPDEmappingT:'
     write(*, *) '         Incompatible input [eVec]. Exiting.'
     
     STOP
      end select
      
   end function dPDEmappingT
   
   !**
   ! SetType
   !*
   subroutine SetType(self, paramType)
      ! Arguments
      class(ModelParameterCell_SG_t), intent(inout) :: self
      character(*)                 , intent(in) :: paramType
      
      if (.not.(self%isAllocated)) then
     write(*, *) 'ERROR:ModelPArameterCell_t:SetType:'
     write(*, *) '         Not allocated.'
     STOP
      end if
      
      if (trim(paramType) .eq. trim(self%paramType)) then
     ! We are done
      else if (self%paramType == "") then
           self%paramType = trim(paramType)
      else if (self%paramType == LINEAR) then
     ! Convert to log
     if (paramType == LOGE) then
          self%cellCond%v = log(self%cellCond%v)
          self%airCond = log(self%AirCond)
          
     else if (paramType == LOG_10) then
          self%cellCond%v = log10(self%cellCond%v)
          self%airCond = log10(self%airCond)
     end if
      else if (paramType == LINEAR) then
     ! Convert from log to linear
     if (self%paramType == LOGE) then
          self%cellCond%v = exp(self%cellCond%v)
          self%airCond = exp(self%airCond)
     else if(self%paramType == LOG_10) then
          self%cellCond%v = exp(self%cellCond%v * log(10.))
          self%airCond = exp(self%AirCond * log(10.))
     end if
      else if ((self%paramType == LOGE) .and. (paramType == LOG_10)) then
     ! Convert from natural log to log10
     self%cellCond%v = self%cellCond%v / log(10.)
     self%airCond = self%airCond / log(10.)
      else if ((self%paramType == LOG_10) .and. (paramType == LOGE)) then
     ! Convert from log10 to natural log
     self%cellCond%v = self%cellCond%v * log(10.)
     self%airCond = self%airCond * log(10.)
      else
     write(*, *) 'ERROR:ModelParameterCell_t:SetType:'
     write(*, *) '         Unknown paramType.'
     STOP
      end if
      
      self%paramType = paramType
   end subroutine SetType
   
end Module ModelParameterCell_SG
