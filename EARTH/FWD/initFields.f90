! *****************************************************************************
module initFields
  ! Module containing the subroutines to initialize the field vectors.
  ! No memory allocation is done here. Vectors are assumed to already exist.

  use field_vectors
  use dimensions
  use ringcurrent
  use coreFwd
  use wrapper
  use sg_sparse_vector
  use griddef
  implicit none

  private
      ! store the grid dimensions in the original format
  integer                                       :: nx, ny, nz, nzEarth, nzAir
      ! storing the (spherical) grid in Randie Mackie's format
      ! when allocated, dimensions will be x(nx), y(ny+1), z(nz+1)
  real(8), allocatable, dimension(:)            :: x,y,z


  public        :: initialize_fields
  public        :: initialize_fields_vec
  public        :: initialize_div_corr

Contains


  ! ***************************************************************************
  ! * initialize_fields is the routine to compute the initial values of the
  ! * magnetic fields based on the boundary conditions. Currently there is an
  ! * ambiguity that has to be dealt with, in that the _vertical_ outer boundary
  ! * information is not saved in the transformation Hx/Hy/Hz -> vectorh.
  ! * Hence the converse transformation vectorh -> Hx/Hy/Hz does not reconstruct
  ! * the original boundary values, but rather sets them to zero.
  ! * If Hx/Hy/Hz have not had the 'save' attribute, the boundary information in
  ! * them is soon lost. To avoid this problem we apply a patch by saving the
  ! * upper and lower boundary fields in bvH.
  ! * Then Hx/Hy/Hz vectors should be identical to their initialized values after
  ! * the following two-line transformation:
  !	* call copyd1_d3_b(nx,ny,nz,Hx,Hy,Hz,hvec,x,y,z)
  ! * call insertBoundaryValues(bvH,Hx,Hy,Hz)

  subroutine initialize_fields_vec(hvec,bvec,grid,rho,bvH)


	complex(8), dimension(:), intent(inout)	:: hvec,bvec
    real(8), dimension(:,:,:), intent(in)     :: rho
	type (sparsevecc), intent(inout), optional :: bvH
	type (grid_t), intent(in)               :: grid
	integer									:: istat

    nx = grid%nx; ny = grid%ny; nz = grid%nz
    nzEarth = grid%nzEarth; nzAir = grid%nzAir
    allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
    x = grid%x; y = grid%y; z = grid%z

	!allocate(Hx(np1),Hy(np1),Hz(np1),STAT=istat)

	! Initialise the H field at the outer boundary, then make it
	! depend linearly on radius throughout the computational domain
	! such that all field components are zero at the CMB. The P10
	! source at the boundary has been modified to the more accurate
	! initialization by Ikuko Fujii allowing for the ringcurrent.
	call set_initial_p10(nx,ny,nz,x,y,z,Hx,Hy,Hz,np1)
	! call set_initial_ringcurrent(nx,ny,nz,x,y,z,Hx,Hy,Hz,np1)
	if (present(bvH)) then
	   call initBoundaryValues(bvH,Hx,Hy,Hz,grid)
	end if

	! save the boundary values of Hz
	! call insertBoundaryValues(bvH,Hx,Hy,Hz)

	! initialize vector <h> of A <h> = <b>
	call copyd3_d1_b(nx,ny,nz,Hx,Hy,Hz,hvec,x,y,z)

	! compute <b> of A <h> = <b>
	call calcb(nx,ny,nz,Hx,Hy,Hz,bvec,rho,x,y,z)

	!deallocate(Hx,Hy,Hz)

	!do i=1,bvH%ncoeff
	  !print *,'bvH: ',bvH%i(i),bvH%j(i),bvH%k(i),bvH%xyz(i),bvH%c(i)
	!end do
	!call create_cvector(grid, Bzero, EDGE)
	!call compare_fields(B, Bzero)

    deallocate(x,y,z,STAT=istat)

  end subroutine initialize_fields_vec



  ! ***************************************************************************
  ! * initialize_fields is the routine to compute the initial values of the
  ! * magnetic fields based on the boundary conditions
  subroutine initialize_fields(grid,H,Hb)

	!use field_vectors

	type (grid_t), intent(in)               :: grid
	type (cvector), intent(inout)			:: H
	type (sparsevecc), intent(inout)	    :: Hb
    type (cvector)                          :: F
	complex(8), dimension(:), allocatable	:: hvec
	!complex(8), dimension(:), allocatable   :: Hx,Hy,Hz
	integer									:: istat,i

    nx = grid%nx; ny = grid%ny; nz = grid%nz
    nzEarth = grid%nzEarth; nzAir = grid%nzAir
    allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
    x = grid%x; y = grid%y; z = grid%z

	! Local allocations; perhaps not required if the variables become global
	allocate(hvec(np2),STAT=istat)
	!allocate(hx(np1),hy(np1),hz(np1),STAT=istat)

	!---- Initialise the H field at the outer boundary, then make it
	!     depend linearly with radius throughout the computational domain
	!     such that all field components are zero at the CMB. The P10
	!		source at the boundary has been modified to the more accurate
	!		initialization by Ikuko Fujii allowing for the ringcurrent.
	call set_initial_p10(nx,ny,nz,x,y,z,Hx,Hy,Hz,np1)
	! NB --- This is now done in InitGlobalArrays
	! call set_initial_ringcurrent(nx,ny,nz,x,y,z,hx,hy,hz,np1)
	call copyd3_d1_b(nx,ny,nz,hx,hy,hz,hvec,x,y,z)	  ! extract interior components
	call initBoundaryValues(Hb,hx,hy,hz,grid)		  ! extract boundary components

	! Initialize the magnetic field
	call copyd1_d3_d(hvec,H,grid,bc=Hb)

	! Initialize the forcing with the boundary magnetic fields
	!call create_cvector(grid,F,EDGE)
	!call add_scvector(C_ONE,Hb,F)

	!deallocate(hx,hy,hz)
	deallocate(hvec)
	!call deall_sparsevecc(Hb)
    deallocate(x,y,z,STAT=istat)

  end subroutine initialize_fields

  ! ***************************************************************************
  ! * initialize_div_corr is the routine that initializes vectors required
  ! * for divergence correction. Input has to be the true H fields and true
  ! * interior forcing for the system s (no boundary conditions included).
  ! * Optionally, insert the boundary conditions in hx,hy,hz - this should
  ! * only be done when initialized for FWD (NOT for secondary field formulation
  ! * and NOT for sensitivity - use zero BC in these cases).

  subroutine initialize_div_corr(hvec,svec,grid,bvH)

    complex(8), dimension(:), intent(inout) :: hvec,svec
    type (sparsevecc), intent(in), optional  :: bvH
    type (grid_t), intent(in)               :: grid
    integer                                 :: istat

    nx = grid%nx; ny = grid%ny; nz = grid%nz
    nzEarth = grid%nzEarth; nzAir = grid%nzAir
    allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
    x = grid%x; y = grid%y; z = grid%z

    ! Set initial conditions for the divergence correction - interior source
    sx = C_ZERO
    sy = C_ZERO
    sz = C_ZERO
    call copyd1_d3_b(nx,ny,nz,sx,sy,sz,svec,x,y,z)

    ! Set initial conditions for the divergence correction - magnetic fields
    hx = C_ZERO
    hy = C_ZERO
    hz = C_ZERO
    call copyd1_d3_b(nx,ny,nz,hx,hy,hz,hvec,x,y,z)

    ! This step is only needed while the divergence correction uses b.c.
    if (present(bvH)) then
        call insertBoundaryValues(bvH,hx,hy,hz,grid)
    end if

    deallocate(x,y,z,STAT=istat)

  end subroutine initialize_div_corr


end module initFields
