Module nestedEM
  use EMfieldInterp
  use sg_boundary			! work between different data types
  					! (between boundary conditions and
					! complex vectors)
  use sg_sparse_vector !, only: add_scvector

  implicit none

!################################ saved and used for the model operators ##################################

  !Vector to hold the BC interpolated from a file E solution file for all transmitters and polarizations
  type(cboundary), pointer, dimension(:), save    :: BC_from_file
  type(cvector), pointer, dimension(:), save      :: E0_from_file


!###################################### used for nested modelling ############################################  
    type :: location
     ! p is the three dimensional coordinate for a location
     real(kind=prec), dimension(3)                ::  p

     ! in some instances, it is prudent to store for what component of
     ! electrical or magnetic field is the location being used for
     ! xyz = component where xyz = 1 for x-component, xyz = 2 for y-
     ! component, and xyz = 3 for z-component. 0 means that that location
     ! is either for mutliple components or neutral
     integer                            :: xyz = 0

     ! Remember, ix, iy, and iz (array locations) are not needed in
     ! every instance of location only 
     ! So that we always know where they come from in the array of interest 
     ! during the nested like intepolation. Why are they needed?
     ! Refer to NestedE: that is how we make sure that the intepolated 
     ! values are  put in the right position for eOut cvector by carrying the
     ! ix, iy, and iz values of the slices all along
     integer                            :: ix, iy, iz

  end type location
  
  ! local variables set up by CreateSlices; saved and used for nesting calculations
  type(location), pointer, dimension(:), save     :: xTslice,xBslice
  type(location), pointer, dimension(:), save     :: yTslice,yBslice
  type(location), pointer, dimension(:), save     :: zTslice,zBslice
  integer, save                                   :: xslSites, yslSites, zslSites
  logical, save                                   :: slicesSetUp = .false.

  !  the following are private data types that are just used internally
  ! defining a electrical field interpolation site and then using them
  ! for interpolation
  type :: einterp_site

     ! x gives location of electrical field calculation
     type (location)                    ::  l

     ! these sparse vectors contain coefficients of measurement
     ! functionals for electric field components at this location
     !type(sparsevecc)                   ::  e


  end type einterp_site

type(sparsevecc), pointer, dimension(:)                         ::  e_sparse_vector

  ! a dictionary that stores electrical fields that need to be intepolated.
  ! one possible use is for interpolating electrical fields for nested modeling
  type (einterp_site), pointer, save, dimension (:)		::enestDict
  
  ! total number of sites
  integer, save							:: totalSites  

!##########################################################################################################


Contains

  !**********************************************************************
  subroutine setup_BC_from_file(grid,nTx,nPol)

  type(grid_t),intent(in)                   ::  grid
  integer, intent(in)                       ::  nTx, nPol
  ! local
  type(cboundary)  :: BC
  integer     :: j,status

    allocate (BC_from_file(nTx*nPol),E0_from_file(nTx*nPol), STAT=status)

    call create_cboundary(grid,BC)
    do j=1,nTx*nPol
      BC_from_file(j)=BC
    end do

    call CreateSlices(grid)

  end subroutine setup_BC_from_file


  !***********************************************************************
  ! computes a single value in the nTx * nPol array of boundary conditions
  ! stored in the BC_from_file variable and saved in this module

  subroutine compute_BC_from_file(Large_Grid,e_large,Grid,counter)

    type(grid_t)  ,intent(in)                 ::  Large_Grid
    type(cvector),intent(in)                  ::  e_large
    type(grid_t),intent(in)                   ::  Grid
    integer, intent(in)                       ::  counter

    ! local variables needed for nesting calculations
    character (len=80)           :: whichFace = ''
    type(cvector)                :: elecSoln


    if (.not. slicesSetUp) then
        call CreateSlices(Grid)
    end if

    call create_cvector(Grid, elecSoln, EDGE)

    whichFace = 'zFace'
   ! top face
       Call NestingSetUp(Large_Grid, zslSites, zTslice, whichFace)
       Call NestedE(e_large, elecSoln)

   ! bottom face
       Call NestingSetUp(Large_Grid, zslSites, zBslice, whichFace)
       Call NestedE(e_large, elecSoln)

       ! the two x-faces
       whichFace = 'xFace'
   ! top face
       Call NestingSetUp(Large_Grid, xslSites, xTslice, whichFace)
       Call NestedE(e_large, elecSoln)
   ! bottom face
       Call NestingSetUp(Large_Grid, xslSites, xBslice, whichFace)
       Call NestedE(e_large, elecSoln)

       ! the two y-faces
       whichFace = 'yFace'
   ! top face
       Call NestingSetUp(Large_Grid, yslSites, yTslice, whichFace)
       Call NestedE(e_large, elecSoln)
   ! bottom face
       Call NestingSetUp(Large_Grid, yslSites, yBslice, whichFace)
       Call NestedE(e_large, elecSoln)


    Call getBC(elecSoln,  BC_from_file(counter))

    Call deall_cvector(elecSoln)


  end subroutine compute_BC_from_file


  !**********************************************************************
  subroutine deall_BC_from_file()

  ! local
  integer     :: j, status

    deallocate(BC_from_file,E0_from_file, STAT=status)

    deallocate(xTslice,xBslice, STAT=status)
    deallocate(yTslice,yBslice, STAT=status)
    deallocate(zTslice,zBslice, STAT=status)
    slicesSetUp = .false.

  end subroutine deall_BC_from_file


  !**********************************************************************
  subroutine CreateSlices(Grid)

  type(grid_t),intent(in)                   ::  Grid

  ! local variables needed for nesting calculations
  character (len=80)           :: whichFace = '', whatType = ''
  integer                      :: whereAt

     ! we are creating the slices (the outer layers of the mGrid cube)
     ! Values are extracted at the center of the edge
      whatType = 'Edge'

     ! the two z-faces
     whichFace = 'zFace'
     whereAt = 1                ! top
     Call OneSliceAtCenter(Grid, zTslice, zslSites, whichFace, whatType, whereAt)
     whereAt = Grid%nz+1      ! bottom
     Call OneSliceAtCenter(Grid, zBslice, zslSites,whichFace, whatType, whereAt)

     ! the two x-faces
     whichFace = 'xFace'
     whereAt = 1               ! nearest
     Call OneSliceAtCenter(Grid, xTslice, xslSites, whichFace, whatType, whereAt)
     whereAt = Grid%nx+1      ! farthest
     Call OneSliceAtCenter(Grid, xBslice, xslSites, whichFace, whatType, whereAt)

     ! the two y-faces
     whichFace = 'yFace'
     whereAt = 1               ! nearest
     Call OneSliceAtCenter(Grid, yTslice, yslSites,whichFace, whatType, whereAt)
     whereAt = Grid%ny+1      ! farthest
     Call OneSliceAtCenter(Grid, yBslice, yslSites,whichFace, whatType, whereAt)

     slicesSetUp = .true.

    end subroutine CreateSlices


  ! ***************************************************************************
  ! * BOP
  ! creates an array that has the coordinates for a given slice at the center.
  ! the slice is described by "whichFace" i.e. whichFace = zFace 
  ! (horizontal), xFace, and yFace. whatType = 'Edge' or 'Face'. 
  ! whereAt = position of the slice. OneSliceAtCenter could be on the edge or
  ! face. Remember, the absolute locations are global but ix, iy, iz are wrt
  ! the grid inputted 
  subroutine OneSliceAtCenter(inGrid, loc, num, whichFace, whatType, whereAt)

    implicit none
    type (grid_t), intent(in) 	        			:: inGrid
    type (location), pointer, dimension(:)	:: loc
    integer, intent(out)					:: num
    ! whichFace = 'xFace', 'yFace', 'zFace' 
    character (len=80), intent(in)     				:: whichFace
    ! whatType = 'Edge', 'Face'
    character (len=80), intent(in)				:: whatType
    ! whereAt = 1, nx+1; 1, ny+1; 1, nz+1 
    integer, intent(in)						:: whereAt
    ! * EOP

    integer 							:: i, ix, iy
    integer							:: iz
    integer							:: status    
    
   ! no baggage
    deallocate(loc, STAT = status)

    if (whichFace == 'zFace') then
       if (whereAt.gt.(inGrid%nz+1)) then
          write(0, *) 'whereAt: Out of bounds in OneSlice' 
          stop
       end if
       if (whatType == 'Edge') then
          ! total number for a given slice	
          num = inGrid%nx*(inGrid%ny+1)+(inGrid%nx+1)*inGrid%ny
          allocate(loc(num), STAT=status)
          i = 0
          do ix = 1, inGrid%nx
             do iy = 1, inGrid%ny+1
                ! x-nodes
                i = i + 1
                loc(i)%p(1) = inGrid%xCenter(ix)
                loc(i)%p(2) = inGrid%yEdge(iy)
                loc(i)%p(3) = inGrid%zEdge(whereAt)
		loc(i)%xyz = 1
		loc(i)%ix = ix
		loc(i)%iy = iy
		loc(i)%iz = whereAt
             end do
          end do
          do ix = 1, inGrid%nx+1
             do iy = 1, inGrid%ny	
                ! y-nodes 
                i = i + 1
                loc(i)%p(1) = inGrid%xEdge(ix)
                loc(i)%p(2) = inGrid%yCenter(iy)
                loc(i)%p(3) = inGrid%zEdge(whereAt)
		loc(i)%xyz = 2
		loc(i)%ix = ix
		loc(i)%iy = iy
		loc(i)%iz = whereAt
             end do
          end do

       else if (whatType == 'Face') then
          ! total number for a given slice	
          num = (inGrid%nx)*(inGrid%ny)
          allocate(loc(num), STAT=status)
          i = 0
          do ix = 1, inGrid%nx
             do iy = 1, inGrid%ny
                ! z-nodes
                i = i + 1
                loc(i)%p(1) = inGrid%xCenter(ix)
                loc(i)%p(2) = inGrid%yCenter(iy)
                loc(i)%p(3) = inGrid%zEdge(whereAt)
		loc(i)%ix = ix
		loc(i)%iy = iy
		loc(i)%iz = whereAt 
             end do
          end do

       else 
          write (0, *) 'Redefine the type of slice needed: OneSlice'
       end if

    else if (whichFace == 'xFace') then
       if (whereAt.gt.(inGrid%nx+1)) then
          write(0, *) 'whereAt: Out of bounds in OneSlice' 
          stop
       end if
       if (whatType == 'Edge') then	
	  ! total number for a given slice	
          num = inGrid%ny*(inGrid%nz+1)+(inGrid%ny+1)*inGrid%nz
          allocate(loc(num), STAT=status)
	  i = 0
          do iy = 1, inGrid%ny
             do iz = 1, inGrid%nz+1
                ! y-nodes
                i = i + 1
                loc(i)%p(1) = inGrid%xEdge(whereAt)
                loc(i)%p(2) = inGrid%yCenter(iy)
                loc(i)%p(3) = inGrid%zEdge(iz)
		loc(i)%xyz = 2
		loc(i)%ix = whereAt
		loc(i)%iy = iy
		loc(i)%iz = iz
             end do
          end do
          do iy = 1, inGrid%ny+1
             do iz = 1, inGrid%nz	
                ! z-nodes 
                i = i + 1
                loc(i)%p(1) = inGrid%xEdge(whereAt)
                loc(i)%p(2) = inGrid%yEdge(iy)
                loc(i)%p(3) = inGrid%zCenter(iz)
		loc(i)%xyz = 3
		loc(i)%ix = whereAt
		loc(i)%iy = iy
		loc(i)%iz = iz
             end do
          end do

       else if (whatType == 'Face') then 
          ! total number for a given slice	
          num = (inGrid%ny)*(inGrid%nz)
          allocate(loc(num), STAT=status)
	  i = 0
          do iy = 1, inGrid%ny
             do iz = 1, inGrid%nz
                i = i + 1
                ! x-nodes
                loc(i)%p(1) = inGrid%xEdge(whereAt)
                loc(i)%p(2) = inGrid%yCenter(iy)
                loc(i)%p(3) = inGrid%zCenter(iz)
		loc(i)%ix = whereAt
		loc(i)%iy = iy
		loc(i)%iz = iz 
             end do
          end do

       else 
          write (0, *) 'Redefine the type of slice needed: OneSlice'
       end if

    else if (whichFace == 'yFace') then
       if (whereAt.gt.(inGrid%ny+1)) then
          write(0, *) 'whereAt: Out of bounds in OneSlice' 
          stop
       end if
       if (whatType == 'Edge') then
	  ! total number for a given slice	
          num = inGrid%nx*(inGrid%nz+1)+(inGrid%nx+1)*inGrid%nz
          allocate(loc(num), STAT=status)
	  i = 0
	  do ix = 1, inGrid%nx+1
             do iz = 1, inGrid%nz	
                ! z-nodes 
                i = i + 1
                loc(i)%p(1) = inGrid%xEdge(ix)
                loc(i)%p(2) = inGrid%yEdge(whereAt)
                loc(i)%p(3) = inGrid%zCenter(iz)
		loc(i)%xyz = 3 
		loc(i)%ix = ix
		loc(i)%iy = whereAt
		loc(i)%iz = iz
             end do
          end do
          do ix = 1, inGrid%nx
             do iz = 1, inGrid%nz+1
                ! x-nodes
                i = i + 1
                loc(i)%p(1) = inGrid%xCenter(ix)
                loc(i)%p(2) = inGrid%yEdge(whereAt)
                loc(i)%p(3) = inGrid%zEdge(iz)
		loc(i)%xyz = 1
		loc(i)%ix = ix
		loc(i)%iy = whereAt
		loc(i)%iz = iz
             end do
          end do

       else if (whatType == 'Face') then
          ! total number for a given slice	
          num = (inGrid%nx)*(inGrid%nz)
          allocate(loc(num), STAT=status)
	  i = 0
          do ix = 1, inGrid%ny
             do iz = 1, inGrid%nz
                i = i + 1
                ! y-nodes
                loc(i)%p(1) = inGrid%xCenter(ix)
                loc(i)%p(2) = inGrid%yEdge(whereAt)
                loc(i)%p(3) = inGrid%zCenter(iz)
		loc(i)%ix = ix
		loc(i)%iy = whereAt
		loc(i)%iz = iz 
             end do
          end do

       else 
          write (0, *) 'Redefine the type of slice needed: OneSlice'
       end if

    else 

       write(0, *) 'Redefine the face of slice needed: OneSlice'

    end if

  end subroutine OneSliceAtCenter
  
!##############################################################  
  ! **************************************************************************
  ! * BOP
  ! Sets up the electrical fields interpolation for nested models. 
  ! both the init and set up are to done together. inGrid = larger grid.
  ! the inGrid could be an actual grid or pointer to the grid with all the
  ! information. But remember the inGrid has to be the larger grid
  ! siteLocations are from a smaller grid but corrected for the origin or an
  ! arbitrary set of locations within the same large location. 
  ! Remember, the absolute location are global but ix, iy, iz are wrt 
  ! the grid inputted.
  ! Call NestingSetUp(inGrid,nSites,siteLocations, whichFace) first
  ! and then call NestedE(ef, eOut) for complete intepolation together.
  ! However, the siteLocations(nSites) needs to be calculated before hand.
  subroutine NestingSetUp(inGrid,nSites,siteLocations, whichFace)

    implicit none
    type (grid_t), target, intent(in) 	:: inGrid
    integer, intent(in)	 		:: nSites
    type (location), intent(in)		:: siteLocations(nSites)
    ! whichFace = 'xFace', 'yFace', 'zFace' 
    character (len=80), intent(in)     	:: whichFace
    ! * EOP

    type (sparsevecc)			:: LC
    integer				:: i, xyz
    integer				:: status

    ! in most cases, the allocation is not done in set up. this is an
    ! exception to the rule
    allocate(e_sparse_vector(nSites), STAT=status)
    allocate(enestDict(nSites), STAT=status)
    
    do i = 1,nSites
       enestDict(i)%l%p = siteLocations(i)%p
       enestDict(i)%l%xyz = siteLocations(i)%xyz
       enestDict(i)%l%ix = siteLocations(i)%ix
       enestDict(i)%l%iy = siteLocations(i)%iy
       enestDict(i)%l%iz = siteLocations(i)%iz

       totalSites = nSites


       if (whichFace == 'zFace') then       
          if (siteLocations(i)%xyz == 1) then
	  ! sparse vector for ex measurement
             xyz = 1 
             Call EinterpSetUp(inGrid,siteLocations(i)%p,xyz,LC)      
             e_sparse_vector(i) = LC
          else if (siteLocations(i)%xyz == 2) then
          ! sparse vector for ey measurement
             xyz = 2 
             Call EinterpSetUp(inGrid,siteLocations(i)%p,xyz,LC)      
             e_sparse_vector(i)= LC
	  else 
	      write (0, *) 'wrong component in NestingSetUp for ', whichFace  
	      stop
          end if
       end if

       if (whichFace == 'xFace') then       
          if (siteLocations(i)%xyz == 2) then
             xyz = 2 
             Call EinterpSetUp(inGrid,siteLocations(i)%p,xyz,LC)      
             e_sparse_vector(i) = LC
          else if (siteLocations(i)%xyz == 3) then
             ! sparse vector for ey measurement
             xyz = 3
             Call EinterpSetUp(inGrid,siteLocations(i)%p,xyz,LC)      
             e_sparse_vector(i) = LC
	  else 
	      write (0, *) 'wrong component in NestingSetUp for ', whichFace 
	      stop
          end if
       end if

       if (whichFace == 'yFace') then
          ! remember, same circular positioning of nodes, 
          ! as in OneSliceAtCenter
	  ! sparse vector for ez measurement
          if (siteLocations(i)%xyz == 3) then
             xyz = 3
             Call EinterpSetUp(inGrid,siteLocations(i)%p,xyz,LC)      
             e_sparse_vector(i) = LC
          else if (siteLocations(i)%xyz == 1) then
             ! sparse vector for ex measurement
             xyz = 1
             Call EinterpSetUp(inGrid,siteLocations(i)%p,xyz,LC)      
             e_sparse_vector(i) = LC
	  else 
	      write (0, *) 'wrong component in NestingSetUp for ', whichFace  
	      stop
          end if
       end if

    enddo

  end subroutine NestingSetUp
    ! **************************************************************************
  ! * BOP
  ! Intepolates the electrical fields for a slice used for nested modeling. 
  ! The input is in cvector data type, eIn and the output eOut is also a 
  ! cvector. This is a very generic routine that has cvector as an output.
  ! Call NestingSetUp(inGrid,nSites,siteLocations, whichFace) first
  ! and then call NestedE(ef, eOut) for complete intepolation together.
  ! However, the siteLocations(nSites) needs to be calculated before hand.
  ! eOut has to be contained within ef
  subroutine NestedE(ef, eOut)

    implicit none
    type (cvector), intent(in) 				:: ef
    type (cvector), intent(inout)			:: eOut 
    ! * EOP 

    complex (kind=prec)					:: c
    integer 						:: i, ix, iy
    integer						:: iz
    integer						:: xyz
    integer						:: status
  
    if(.not.ef%allocated) then
       write(0,*) 'input not allocated yet for NestedE_E'
       stop
    end if

    if(.not.eOut%allocated) then
       write(0,*) 'input not allocated yet for NestedE_E'
       stop
    end if

    do i = 1, totalSites
       c =dotProd(e_sparse_vector(i),ef)
       ! that is how we make sure that the intepolated values are 
       ! put in the right position eOut cvector by carrying the
       ! ix, iy, and iz values of the slices all along
       ix = enestDict(i)%l%ix
       iy = enestDict(i)%l%iy
       iz = enestDict(i)%l%iz
       xyz = enestDict(i)%l%xyz
       if (xyz == 1) then 
       	   eOut%x(ix, iy, iz) = c
       else if (xyz == 2) then
           eOut%y(ix, iy, iz) = c
       else if (xyz == 3) then
           eOut%z(ix, iy, iz) = c
	else
	   write (0, *) 'xyz not provided: NestedE'
	   stop
       end if
    end do 

    ! no baggage left behind
    deallocate(enestDict, STAT=status)

    do i=1,totalSites
       call deall_sparsevecc (e_sparse_vector(i))
    end do
    deallocate(e_sparse_vector, STAT=status)
    
     totalSites = 0
     
end subroutine NestedE






end module nestedEM
