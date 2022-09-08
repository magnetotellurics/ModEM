module GridDef
   !  Defines grid data structure, basic grid methods

   use math_constants

   ! Don't forget to overload the '=' sign: depending on the compiler, might
   ! run into trouble with the default assignment, since that sometimes doesn't
   ! copy allocatable or pointer arrays
   INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE copy_grid
   END INTERFACE

   ! Possible grid types, on which cvector is defined. Viewing the grid
   ! as 2D, NODES and EDGES correspond to the corners and the edges of
   ! all square grid elements, and CELLS corresponds to the centers of
   ! these squares. EARTH means the air is excluded from the grid.
   character(len=80), parameter		:: CELL = 'CELL'
   character(len=80), parameter		:: NODE = 'NODE'
   character(len=80), parameter		:: EDGE = 'EDGE'
   character(len=80), parameter		:: NODE_EARTH = 'NODE EARTH'
   character(len=80), parameter		:: CELL_EARTH = 'CELL EARTH'
   character(len=80), parameter		:: EDGE_EARTH = 'EDGE EARTH'

   type ::  grid_t
      !  grid_t is derived data type used to store grid geometry
      integer   :: Nz = 0
      integer   :: Ny = 0
      integer   :: Nza = 0
      real (kind=prec), pointer, dimension(:) :: Dy,Dz
      real (kind=prec), pointer, dimension(:) :: Dely,Delz
      real (kind=prec), pointer, dimension(:) :: yNode,zNode
      real (kind=prec), pointer, dimension(:) :: yCenter,zCenter
      real (kind=prec)  :: zAir = R_ZERO
      logical	:: allocated = .false.
   end type grid_t

   public         :: create_grid, deall_grid, setup_grid, copy_grid

   Contains

     !************************************************************************
     subroutine create_grid(Ny,Nz,Nza,grid)
       !  creates finite differences grid_t structure of
       !  size Nz x Ny, allocates arrays
       !
       implicit none
       integer, intent(in)		        :: Nz,Ny,Nza
       type (grid_t), intent(inout)	:: grid
       integer                          :: istat

       grid%Nza = Nza
       grid%Nz = Nz
       grid%Ny = Ny
       grid%zAir = R_ZERO
       allocate(grid%Dz(Nz), STAT=istat)
       allocate(grid%Dy(Ny), STAT=istat)
       allocate(grid%Delz(Nz+1), STAT=istat)
       allocate(grid%Dely(Ny+1), STAT=istat)
       allocate(grid%zNode(Nz+1), STAT=istat)
       allocate(grid%yNode(Ny+1), STAT=istat)
       allocate(grid%zCenter(Nz), STAT=istat)
       allocate(grid%yCenter(Ny), STAT=istat)
       grid%allocated = .true.

     end subroutine create_grid

  ! **************************************************************************
  logical function valid_grid(grid)
    ! basic sanity check for grid manipulations - can be expanded in the future [AK]
  
    type (grid_t), intent(inout) 	:: grid
    
    valid_grid = .true.
    
    if ((grid%Ny <= 0) .or. (grid%Nz <= 0)) then
    	write(0,*) 'Grid information: Ny=',grid%Ny,' Nz=',grid%Nz
    	valid_grid = .false.
    end if
    
  end function valid_grid
  
     !************************************************************************
     subroutine deall_grid(grid)
       !  deallocates finite differences grid_t structure
       !
       implicit none
       type (grid_t), intent(inout)	:: grid
       integer                          :: istat

       if (grid%allocated) then
       if (associated(grid%Dz)) deallocate(grid%Dz, STAT=istat)
       if (associated(grid%Dy)) deallocate(grid%Dy, STAT=istat)
       if (associated(grid%Delz)) deallocate(grid%Delz, STAT=istat)
       if (associated(grid%Dely)) deallocate(grid%Dely, STAT=istat)
       if (associated(grid%zNode)) deallocate(grid%zNode, STAT=istat)
       if (associated(grid%yNode)) deallocate(grid%yNode, STAT=istat)
       if (associated(grid%zCenter)) deallocate(grid%zCenter, STAT=istat)
       if (associated(grid%yCenter)) deallocate(grid%yCenter, STAT=istat)
       end if

     end subroutine deall_grid

     !***********************************************************************
     ! setup_grid: after allocation, read in Dy, Dz, then call:
     subroutine setup_grid(grid)
       implicit none
       type (grid_t) , intent(inout)	:: grid
       !  local variables
       integer ::	Nz,Ny,iy,iz,Nza

       Nz = grid%Nz
       Nza = grid%Nza
       Ny = grid%Ny
       do iy = 2,Ny
          grid%Dely(iy) = (grid%Dy(iy-1)+grid%Dy(iy))/2.
       enddo
       grid%Dely(1) = grid%Dely(1)/2.
       grid%Dely(Ny+1) = grid%Dely(Ny)/2.
       do iz = 2,Nz
            grid%Delz(iz) = (grid%Dz(iz-1)+grid%Dz(iz))/2.
        enddo
        grid%Delz(1) = grid%Delz(1)/2.
        grid%Delz(Nz+1) = grid%Delz(Nz)/2.

        !  for now assume y0 = 0, z0 = 0 (left edge, Earth surface)
        grid%yNode(1) = 0
        grid%zNode(1) = 0
        do iy = 1,Ny
           grid%yNode(iy+1) = grid%yNode(iy) + grid%Dy(iy)
        enddo
        ! start summing from top of domain (includes air)
        grid%zAir = 0
        do iz = 1,Nz
           grid%zNode(iz+1) = grid%zNode(iz) + grid%Dz(iz)
           ! also save the total thickness of the air layers
           if(iz .le. Nza) then
             grid%zAir = grid%zAir + grid%Dz(iz)
           endif
        enddo
        !  next construct positions of cell centers
        do iy = 1,Ny
           grid%yCenter(iy) = (grid%yNode(iy)+grid%yNode(iy+1))/2.;
        enddo
        do iz = 1,grid%Nz
           grid%zCenter(iz) = (grid%zNode(iz)+grid%zNode(iz+1))/2.;
        enddo
     end subroutine setup_grid

     !************************************************************************
     subroutine copy_grid(gridOut,gridIn)
       !  overloads the '=' sign
       !
       implicit none
       type (grid_t), intent(in)		:: gridIn
       type (grid_t), intent(inout)		:: gridOut
       integer		        			:: Nz,Ny,Nza

       Ny = gridIn%Ny
       Nz = gridIn%Nz
       Nza = gridIn%Nza

       call deall_grid(gridOut)
       call create_grid(Ny,Nz,Nza,gridOut)

       gridOut%Dy = gridIn%Dy
       gridOut%Dz = gridIn%Dz
       gridOut%Dely = gridIn%Dely
       gridOut%Delz = gridIn%Delz
       gridOut%yNode = gridIn%yNode
       gridOut%zNode = gridIn%zNode
       gridOut%yCenter = gridIn%yCenter
       gridOut%zCenter = gridIn%zCenter
       gridOut%zAir = gridIn%zAir
       gridOut%allocated = .true.

     end subroutine copy_grid

end module GridDef
