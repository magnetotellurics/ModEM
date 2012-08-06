! *****************************************************************************
module jacobian
  ! Module containing the operators required for the Jacobian and penalty
  ! functional derivative calculations

  use iotypes
  use modeldef
  use model_operators
  use modelmap
  use elements
  use dimensions
  use griddef
  use sg_vector
  use sg_scalar
  use sg_spherical
  use wrapper
  implicit none


Contains

  ! ***************************************************************************
  ! * This is the operator Mii (for Maxwell's equations): E_i -> E_i.

  subroutine operatorMii(h,s,omega,rho,grid,fwdCtrls,errflag,adjoint,BC)

	use maxwells
	use initFields
	implicit none

    type (cvector), intent(inout)             :: h
    type (cvector), intent(in)                :: s
    !complex(8),dimension(np2),intent(inout)   :: vectorh
    !complex(8),dimension(np2),intent(in)      :: vectors
	real(8), intent(in)						  :: omega
	type (rscalar), intent(in)	              :: rho
	type (grid_t), intent(in)			      :: grid
	type (fwdCtrl_t), intent(in)			  :: fwdCtrls
	integer, intent(out)					  :: errflag
	logical, intent(in)		                  :: adjoint
	type (sparsevecc), intent(in), optional   :: BC
	! local
	real(8)									  :: om
    complex(8),dimension(:),allocatable       :: vectorx,vectory,vectorb,vectors,vectorh
    logical                                   :: sens
	integer									  :: istat

	! Indicator of whether the forward solver or the adjoint has been called
	sens = .true.
	if(adjoint .and. present(BC)) then
	  write(0,*) 'Warning: (operatorMii) adjoint solver called; BC will not be used'
	elseif(present(BC)) then
	  sens = .false.
	else
	  sens = .true.
	end if

    ! Initialize
    allocate(vectorx(np2),vectory(np2),STAT=istat)
    allocate(vectorb(np2),vectors(np2),vectorh(np2),STAT=istat)

    ! Check that input and output are initialized
    if(.not. s%allocated) then
       write(0,*) 'Error: (operatorMii) input vector not allocated, exiting...'
       stop
    endif
    call copyd3_d1_d(vectors,s,grid)  ! extract interior source components

    if(.not. h%allocated) then
       ! Output vector not allocated, initializing now...
       call create_cvector(grid,h,EDGE)
    endif
    call copyd3_d1_d(vectorh,h,grid)  ! extract interior field components

	! If we are computing the data sensitivities start from zero
	if (sens) then
	  vectorh = C_ZERO
	end if

    ! Initialize all temporary vectors
	vectorx = C_ZERO
	vectory = C_ZERO
	vectorb = C_ZERO

	! Initialize the internal grid for SolveMaxwells
    !write(0,*) 'Preparing SolveMaxwells - initializing the grid; nzEarth =',grid%nzEarth
    nx = grid%nx; ny = grid%ny; nz = grid%nz
    nzEarth = grid%nzEarth; nzAir = grid%nzAir
    allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
    x = grid%x; y = grid%y; z = grid%z

    !-----------------------------------
    ! Set up the RHS <y> for A <x> = <y>
    !-----------------------------------
    vectory = vectors
    if (adjoint) then
      call divide_vec_by_l(nx,ny,nz,vectory,grid%x,grid%y,grid%z)
    else
      call mult_vec_by_S(nx,ny,nz,vectory,grid%x,grid%y,grid%z)
    end if
    ! For the forward solver, compute forcing from BCs and add to the RHS
    ! (DO NOT ADD for data sensitivities or secondary field formulation)
    if (.not. sens) then
        call calcb_from_bc(nx,ny,nz,BC,vectorb,rho%v,grid%x,grid%y,grid%z,grid)
        vectory = vectory + vectorb
    end if

    !------------------------------------------------
    ! Set up the initial value of <x> for A <x> = <y>
    !------------------------------------------------
    vectorx = vectorh
    if (adjoint) then
      call divide_vec_by_S(nx,ny,nz,vectorx,grid%x,grid%y,grid%z)
    else
      call mult_vec_by_l(nx,ny,nz,vectorx,grid%x,grid%y,grid%z)
    end if

    !-------------------------------------------------------------------
    ! Perform all preliminary computations for the divergence correction
    !-------------------------------------------------------------------
    ! Regardless of whether adjoint or not, recompute vectors for div(s)
    ! Need to divide by unit areas in both cases. Failure to do so results
    ! in broken code ... I forget why.
    vectors = vectory - vectorb
    call divide_vec_by_S(nx,ny,nz,vectors,grid%x,grid%y,grid%z)
    call copyd1_d3_b(nx,ny,nz,sx,sy,sz,vectors,grid%x,grid%y,grid%z)
    ! Set initial conditions for the divergence correction (stored in h)
    call copyd1_d3_b(nx,ny,nz,hx,hy,hz,vectorh,grid%x,grid%y,grid%z)
    ! This step is only needed while the divergence correction uses b.c.
    if (.not. sens) then
      ! Then use original boundary conditions for hx,hy,hz
      call insertBoundaryValues(BC,hx,hy,hz,grid)
    end if

    !----------------------
    ! Set the sign of omega
    !----------------------
    ! FOR TRN ALWAYS +omega
    !----------------------
    !if (adjoint) then
    !  om = - omega
    !else
      om = omega
    !end if

    !------------------------------------------------
    ! Call the forward solver $A_{\rho,\omega} x = y$
    !------------------------------------------------
    write(0,*) 'Starting SolveMaxwells'

    call SolveMaxwells(vectorx,vectory,om,rho%v,fwdCtrls,errflag)

    write(0,*) 'Completed SolveMaxwells'

    !-----------------------------------------------
    ! From <x>, compute the interior components of h
    !-----------------------------------------------
    vectorh = vectorx
    if (adjoint) then
      call mult_vec_by_S(nx,ny,nz,vectorh,grid%x,grid%y,grid%z)
    else
      call divide_vec_by_l(nx,ny,nz,vectorh,grid%x,grid%y,grid%z)
    end if

    !-------------------------------------------------
    ! Just making sure the output BCs are also correct
    !-------------------------------------------------
    call zero_cvector(h)
    if (.not. sens) then
        call copyd1_d3_d(vectorh,h,grid,bc=BC)
    else
        call copyd1_d3_d(vectorh,h,grid)
    end if

	!--------------
	! Done; exiting
	!--------------
    deallocate(vectorx,vectory,vectorh,vectorb,vectors)
	deallocate(x,y,z,STAT=istat)
	return

  end subroutine operatorMii	! operatorMii


  ! ***************************************************************************
  ! * This is the operator M (for Maxwell's equations). It is more general than
  ! * operator $A_{\rho,\omega}$ : E_i -> E_i (implemented as the subroutine
  ! * SolveMaxwells()). Operator $A_{\rho,\omega}$ acts on a vector defined on
  ! * internal edges, to produce an output also defined on internal edges.
  ! * It would not need to know anything about the input and output vectors,
  ! * if it didn't do the divergence correction, for which the pure values of
  ! * the internal & external magnetic fields and the forcing are required.
  ! * If the adjoint of this operator is required, it is obtained by setting
  ! * the omega which is passed on to the operator to -omega. Therefore, we pass
  ! * the information to operator $A_{\rho,\omega}$ whether it is an adjoint or
  ! * the forward solution that we would like to obtain. If it is in fact an
  ! * adjoint, two things happen inside SolveMaxwells(). First, omega is set to
  ! * minus omega; second, divergence correction takes into account that to obtain
  ! * the pure fields different transformations are required.
  ! *
  ! * Operator M here is the operator the solves the system of equations.
  ! * The input to this operator is assumed to be the boundary conditions (pure
  ! * fields on boundaries), and pure forcing. The output is the internal and
  ! * boundary magnetic fields (both are saved in the same vector H). For the
  ! * forward operator the boundary magnetic fields will always be the same as
  ! * the input boundary conditions. However, for the adjoint this will not always
  ! * be the case. In fact, the boundary magnetic fields will be computed using
  ! * both the input boundary conditions and the computed internal magnetic field.
  ! * This has not yet been implemented. Currently, for the adjoint, the boundary
  ! * magnetic fields will always be assumed zero; since for the adjoint the
  ! * solution on the internal nodes does not depend on them, except through the
  ! * divergence correction.

  subroutine operatorM(Xfull,Yfull,omega,rho,grid,fwdCtrls,errflag,adjflag,sens)

    use maxwells
    use field_vectors
    use boundaries
    implicit none

    type (cvector), intent(inout)             :: Xfull
    type (cvector), intent(in)                :: Yfull
    type (sparsevecc)                         :: Xb,Yb
    type (grid_t), intent(in)             :: grid
    real(8), intent(in)                       :: omega
    real(8)                                   :: om
    real(8), dimension(:,:,:), intent(in)     :: rho
    type (fwdCtrl_t), intent(in)              :: fwdCtrls
    integer, intent(inout)                    :: errflag
    logical, intent(inout), optional          :: adjflag
    logical, intent(inout), optional          :: sens
    logical                                   :: adjoint,delta
    complex(8),dimension(:),allocatable       :: vectorx,vectory,vectorb,vectors,vectorh
    integer                                   :: istat

    ! Indicator of whether the forward solver or the adjoint has been called
    if(.not.present(adjflag)) then
      adjoint = .FALSE.
    else
      adjoint = adjflag
    end if
    ! Indicator of whether we are computing the data sensitivities;
    ! if so, starting solution and boundary conditions are zero.
    if(.not.present(sens)) then
      delta = .FALSE.
    else
      delta = sens
    end if

    ! Check that input and output are initialized
    if(.not.Yfull%allocated) then
       write(0,*) 'Error: (operatorM) Input vector not allocated, exiting...'
       stop
    endif
    if(.not.Xfull%allocated) then
       ! Output vector not allocated, initializing now...
       call create_cvector(grid,Xfull,EDGE)
    endif

    allocate(vectorx(np2),vectory(np2),STAT=istat)
    allocate(vectorb(np2),vectors(np2),vectorh(np2),STAT=istat)
    vectorx = C_ZERO
    vectory = C_ZERO
    vectorb = C_ZERO
    vectors = C_ZERO
    vectorh = C_ZERO

    ! Initialize the internal grid for SolveMaxwells
    nx = grid%nx; ny = grid%ny; nz = grid%nz
    nzEarth = grid%nzEarth; nzAir = grid%nzAir
    allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
    x = grid%x; y = grid%y; z = grid%z

    !-----------------------------------
    ! Set up the RHS <y> for A <x> = <y>
    !-----------------------------------
    call copyd3_d1_d(vectors,Yfull,grid)  ! extract interior components
    call extractBC(Yb,Yfull)              ! extract boundary components
    call calcb_from_bc(nx,ny,nz,Yb,vectorb,rho,grid%x,grid%y,grid%z,grid)
    vectory = vectors
    if (adjoint) then
      call divide_vec_by_l(nx,ny,nz,vectory,grid%x,grid%y,grid%z)
    else
      call mult_vec_by_S(nx,ny,nz,vectory,grid%x,grid%y,grid%z)
      vectory = vectory + vectorb
    end if

    !------------------------------------------------
    ! Set up the initial value of <x> for A <x> = <y>
    !------------------------------------------------
    call copyd3_d1_d(vectorh,Xfull,grid)  ! extract interior components
    call extractBC(Xb,Xfull)              ! extract boundary components
    vectorx = vectorh
    if (adjoint) then
      call divide_vec_by_S(nx,ny,nz,vectorx,grid%x,grid%y,grid%z)
    else
      call mult_vec_by_l(nx,ny,nz,vectorx,grid%x,grid%y,grid%z)
    end if

    !-------------------------------------------------------------------
    ! Perform all preliminary computations for the divergence correction
    !-------------------------------------------------------------------
    ! Regardless of whether adjoint or not, recompute vectors for div(s)
    vectors = vectory - vectorb
    call divide_vec_by_S(nx,ny,nz,vectors,grid%x,grid%y,grid%z)
    call copyd1_d3_b(nx,ny,nz,sx,sy,sz,vectors,grid%x,grid%y,grid%z)
    ! Set initial conditions for the divergence correction (stored in Xfull)
    call copyd1_d3_b(nx,ny,nz,hx,hy,hz,vectorh,grid%x,grid%y,grid%z)
    ! This step is only needed while the divergence correction uses b.c.
    if (.not.delta) then
      Xb = Yb ! then use original boundary conditions for H
    end if
    call insertBoundaryValues(Xb,hx,hy,hz,grid)


    !----------------------
    ! Set the sign of omega
    !----------------------
    ! FOR TRN ALWAYS +omega
    !----------------------
    !if (adjoint) then
    !  om = - omega
    !else
      om = omega
    !end if

    !------------------------------------------------
    ! Call the forward solver $A_{\rho,\omega} x = y$
    !------------------------------------------------
    call SolveMaxwells(vectorx,vectory,om,rho,fwdCtrls,errflag)


    !-----------------------------------------------
    ! From <x>, compute the interior components of X
    !-----------------------------------------------
    vectorh = vectorx
    if (adjoint) then
      call mult_vec_by_S(nx,ny,nz,vectorh,grid%x,grid%y,grid%z)
    else
      call divide_vec_by_l(nx,ny,nz,vectorh,grid%x,grid%y,grid%z)
    end if

    !--------------------------------------------------------------------
    ! Boundary conditions: in the future, a calculation will be performed
    !--------------------------------------------------------------------
    if (adjoint) then
      ! Xb will be computed from X_i and Y_b using the expression
      ! $X_b = - \length{b} R_\rho^T \area{i}^{-1} X_i + Y_b$
      ! Since $R_\rho^T$ has not yet been implemented, and since we
      ! do not use X_b for computing data sensitivities, assume zero.
      Xb%c = C_ZERO
      if (.not.delta) then
        write (0, *) 'Warning: (operatorM) Xb has been set to zero; it is not zero'
      end if
    else
      Xb = Yb
    end if

    !---------------------
    ! Copy everything to X
    !---------------------
    call copyd1_d3_d(vectorh,Xfull,grid,bc=Xb)

    !--------------
    ! Done; exiting
    !--------------
    call deall_sparsevecc(Xb)
    call deall_sparsevecc(Yb)
    deallocate(vectorx,vectory,vectorh,vectorb,vectors)
    deallocate(x,y,z,STAT=istat)
    return

  end subroutine operatorM  ! operatorM

  ! ***************************************************************************
  ! * operator D_{S_i}^{-1}: E -> E_i represents the diagonal operator that
  ! * pre-divides the interior components of the input vector X by unit area,
  ! * and nullifies the boundary components to obtain a vector of interior
  ! * components of S^{-1} X. Here the output will include the boundary values,
  ! * which will be zero.

  subroutine operatorD_Si_divide(Xfull,grid)

	type (cvector), intent(inout)			  :: Xfull
	type (grid_t), intent(in)			  :: grid
	complex(8),dimension(:),allocatable		  :: vectorx
	integer									  :: istat,nx,ny,nz

	nx = grid%nx
	ny = grid%ny
	nz = grid%nz

	allocate(vectorx(np2),STAT=istat)
	call copyd3_d1_d(vectorx,Xfull,grid)  ! extract interior components
	call divide_vec_by_S(nx,ny,nz,vectorx,grid%x,grid%y,grid%z)
	call copyd1_d3_d(vectorx,Xfull,grid)  ! map back with zero b.c.
	deallocate(vectorx)

  end subroutine operatorD_Si_divide


  ! ***************************************************************************
  ! * operator D_{l}^{-1}: E -> E represents the diagonal operator that
  ! * pre-divides all components of the input vector X by unit lengths,
  ! * including the boundary conditions

  subroutine operatorD_l_divide(Xfull,grid)

	use boundaries

	type (cvector), intent(inout)			  :: Xfull
	type (sparsevecc)						  :: Xb,bc
	type (grid_t), intent(in)			  :: grid
	complex(8),dimension(:),allocatable		  :: vectorx
	integer									  :: istat,nx,ny,nz

	nx = grid%nx
	ny = grid%ny
	nz = grid%nz

	call extractBC(Xb,Xfull)
	call divide_bc_by_l(Xb,bc,grid)
	allocate(vectorx(np2),STAT=istat)
	call copyd3_d1_d(vectorx,Xfull,grid)  ! extract interior components
	call divide_vec_by_l(nx,ny,nz,vectorx,grid%x,grid%y,grid%z)
	call copyd1_d3_d(vectorx,Xfull,grid,bc=bc)  ! map back with b.c.
	deallocate(vectorx)
	call deall_sparsevecc(Xb)
	call deall_sparsevecc(bc)

  end subroutine operatorD_l_divide

! *****************************************************************************
! * multiplies E on FACES by corresponding length elements

  subroutine lFmult_cvector(E)

      integer                   :: nx,ny,nz
      real(8),dimension(:),allocatable      :: x,y,z
      complex(8),dimension(np2) :: rvec
      real(8)                   :: xlen,ylen,zlen
      integer                   :: ic,i,j,k,ii,istat
      type(cvector), intent(inout) :: E
      type(cvector)             :: E1
      type(rvector)             :: lF

      ! Check whether input vector is defined on edges
      if (E%gridType /= FACE) then
        write(0, *) 'Error: (lFmult_cvector) input vector not defined on faces'
        stop
      end if

      E1 = E

      call create_rvector(E%grid,lF,FACE)

      nx = E%grid%nx
      ny = E%grid%ny
      nz = E%grid%nz

      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x = E%grid%x
      y = E%grid%y
      z = E%grid%z

      do i=1,nx
        do k=1,nz
          do j=2,ny
            call leng_xijk2_shifted(i,j,k,x,y,z,xlen)
            lF%x(i,j,k)=xlen
          end do
        end do
      end do

      do j=2,ny
        do k=1,nz
          call leng_yijk2_shifted(j,k,y,z,ylen)
          do i=1,nx
            lF%y(i,j,k)=ylen
          end do
        end do
      end do

      do k=1,nz
        call leng_zijk2_shifted(k,z,zlen)
        lF%y(:,1,k)=zlen
        do j=2,ny
          do i=1,nx
            lF%y(i,j,k)=zlen
          end do
        end do
        lF%y(:,ny+1,k)=zlen
      end do

      call validate_rvector(lF,.true.)
      call diagMult_rcvector(lF,E1,E)

      call deall_rvector(lF)
      call deall_cvector(E1)
      deallocate(x,y,z)

      return
      end subroutine lFmult_cvector

! *****************************************************************************
! * divides E on FACES by corresponding area elements

  subroutine SFdiv_cvector(E)

      integer                   :: nx,ny,nz
      real(8),dimension(:),allocatable      :: x,y,z
      complex(8),dimension(np2) :: rvec
      real(8)                   :: sijk2,sjki2,skij2
      real(8)                   :: ym,yp,zp
      integer                   :: ic,i,j,k,ii,istat
      type(cvector), intent(inout) :: E
      type(cvector)             :: E1
      type(rvector)             :: SFinv

      ! Check whether input vector is defined on edges
      if (E%gridType /= FACE) then
        write(0, *) 'Error: (SFdiv_cvector) input vector not defined on faces'
        stop
      end if

      E1 = E

      call create_rvector(E%grid,SFinv,FACE)

      nx = E%grid%nx
      ny = E%grid%ny
      nz = E%grid%nz

      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x = E%grid%x
      y = E%grid%y
      z = E%grid%z

      ic=0
      do i=1,nx
        do k=2,nz
          do j=2,ny
            ic=ic+1
            call area_sijk2(j-1,k-1,y,z,sijk2)
            SFinv%x(i,j,k) = ONE/sijk2
          end do
        end do
      end do

      do j=1,ny
        do k=2,nz
          ! Zero longitude
          ic=ic+1
          call area_sjki2(nx,1,j,k-1,x,y,z,sjki2)
          SFinv%y(i,j,k) = ONE/sjki2
          do i=2,nx
            ic=ic+1
            call area_sjki2(i-1,i,j,k-1,x,y,z,sjki2)
            SFinv%y(i,j,k) = ONE/sjki2
          end do
        end do
      end do

      do k=2,nz
        ic=ic+1
        ! North pole cap
        j=1
        zp=(z(k)+z(k+1))/2.0d0
        yp=(y(j)+y(j+1))/2.0d0
        skij2=2.0d0*pi*(zp**2)*(1.0d0-dcos(yp))
        SFinv%z(i,j,k) = ONE/skij2
        do j=2,ny
          ! Zero longitude
          ic=ic+1
          call area_skij2(nx,1,j-1,k,x,y,z,skij2)
          SFinv%z(i,j,k) = ONE/skij2
          do i=2,nx
            ic=ic+1
            call area_skij2(i-1,i,j-1,k,x,y,z,skij2)
            SFinv%z(i,j,k) = ONE/skij2
          end do
        end do
        ic=ic+1
        ! South pole cap
        j=ny+1
        zp=(z(k)+z(k+1))/2.0d0
        ym=(y(j)+y(j-1))/2.0d0
        skij2=2.0d0*pi*(zp**2)*(dcos(ym)+1.0d0)
        SFinv%z(i,j,k) = ONE/skij2
      end do

      call validate_rvector(SFinv,.true.)
      call diagMult_rcvector(SFinv,E1,E)

      call deall_rvector(SFinv)
      call deall_cvector(E1)
      deallocate(x,y,z)

      return
      end subroutine SFdiv_cvector

  ! ***************************************************************************
  ! * operator D_{l}: E -> E represents the diagonal operator that
  ! * pre-multiplies all components of the input vector X by unit lengths,
  ! * including the boundary conditions

  subroutine operatorD_l_mult(Xfull,grid)

	use boundaries

	type (cvector), intent(inout)			  :: Xfull
	type (sparsevecc)						  :: Xb,bc
	type (grid_t), intent(in)			  :: grid
	complex(8),dimension(:),allocatable		  :: vectorx
	integer									  :: istat,nx,ny,nz

	nx = grid%nx
	ny = grid%ny
	nz = grid%nz

	call extractBC(Xb,Xfull)
	call mult_bc_by_l(Xb,bc,grid)
	allocate(vectorx(np2),STAT=istat)
	call copyd3_d1_d(vectorx,Xfull,grid)  ! extract interior components
	call mult_vec_by_l(nx,ny,nz,vectorx,grid%x,grid%y,grid%z)
	call copyd1_d3_d(vectorx,Xfull,grid,bc=bc)  ! map back with b.c.
	deallocate(vectorx)
	call deall_sparsevecc(Xb)
	call deall_sparsevecc(bc)

  end subroutine operatorD_l_mult


  ! ***************************************************************************
  ! * operator lC is the full curl operator, lC: E -> F
  ! * It acts on a complex vector defined on edges,
  ! * producing the respective vector defined on faces.
  ! * It is an implementation of the following computations:
  ! $e_\phi(i,j,k)=h_r(i,j,k) + h_\th(i,j,k+1) - h_r(i,j+1,k) - h_\th(i,j,k)$
  ! $e_\th(i,j,k)= - h_r(i,j,k) + h_\phi(i,j,k) + h_r(i+1,j,k) - h_\phi(i,j,k+1)$
  ! $e_r(i,j,k)=h_\th(i,j,k) + h_\phi(i,j+1,k) - h_\th(i+1,j,k) - h_\phi(i,j,k)$
  ! * where h is the vector field defined on edges, pre-multiplied by edge lengths.

  subroutine operatorlC(vecE,vecF,grid)

    implicit none
	type (cvector), intent(in)					 :: vecE
	type (cvector), intent(inout)					 :: vecF
	type (grid_t), intent(in)				 :: grid
    real(8)										 :: xlen,xlenjp,xlenkp
    real(8)										 :: ylen,ylenip,ylenkp
	real(8)										 :: zlen,zlenip,zlenjp
    integer										 :: i,j,k,istat
	integer										 :: nx,ny,nz
	real(8), dimension(:), allocatable			 :: x,y,z

    if (.not.vecE%allocated) then

 	  write(0,*) 'Error: (operatorC) vecE not allocated'
 	  stop

	else if (vecE%gridType /= EDGE) then

 	  write(0,*) 'Error: (operatorC) input vector is not defined on edges'
 	  stop

    endif

    if (.not.vecF%allocated) then

 	  call create_cvector(grid,vecF,FACE)

	else if (vecF%gridType /= FACE) then

 	  write(0,*) 'Error: (operatorC) output vector is not defined on faces'
 	  stop

    endif

	nx = grid%nx
	ny = grid%ny
	nz = grid%nz

	allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
	x = grid%x
	y = grid%y
	z = grid%z

	! Boundary values of h at k=1 and k=nz+1 are fixed during initialization;
	! loop over radii; do not loop over k=nz+1 (fields undefined at nz+2).
    do k=1,nz
	  call leng_zijk(k,z,zlen)
	  zlenip=zlen
	  zlenjp=zlen
	  ! loop over co-latitudes from the North pole to the South
      do j=1,ny
		call leng_yijk(j,k,y,z,ylen)
		call leng_yijk(j,k+1,y,z,ylenkp)
		ylenip=ylen
		! loop over longitudes including i=nx (nx+1'th fields same as i=1)
        do i=1,nx
		  call leng_xijk(i,j,k,x,y,z,xlen)
		  call leng_xijk(i,j+1,k,x,y,z,xlenjp)
		  call leng_xijk(i,j,k+1,x,y,z,xlenkp)

		  vecF%x(i,j,k) = zlen*vecE%z(i,j,k) + ylenkp*vecE%y(i,j,k+1) &
						- zlenjp*vecE%z(i,j+1,k) - ylen*vecE%y(i,j,k)

		  vecF%y(i,j,k) = - zlen*vecE%z(i,j,k) + xlen*vecE%x(i,j,k) &
						+ zlenip*vecE%z(i+1,j,k) - xlenkp*vecE%x(i,j,k+1)

		  vecF%z(i,j,k) = ylen*vecE%y(i,j,k) + xlenjp*vecE%x(i,j+1,k) &
						- ylenip*vecE%y(i+1,j,k) - xlen*vecE%x(i,j,k)

         end do
      end do
    end do

	! x-component at nx+1 (zero longitude)
	vecF%x(nx+1,:,:) = vecF%x(1,:,:)

	! y-component at ny+1 (South pole): undefined
	vecF%y(:,ny+1,:) = C_ZERO

	! z-component at nz+1 (lower domain boundary)
	k=nz+1
      do j=1,ny
		call leng_yijk(j,k,y,z,ylen)
		ylenip=ylen
        do i=1,nx
		  call leng_xijk(i,j,k,x,y,z,xlen)
		  call leng_xijk(i,j+1,k,x,y,z,xlenjp)
		  vecF%z(i,j,k) = ylen*vecE%y(i,j,k) + xlenjp*vecE%x(i,j+1,k) &
						- ylenip*vecE%y(i+1,j,k) - xlen*vecE%x(i,j,k)
         end do
      end do


    deallocate(x,y,z)

  end subroutine operatorlC	! lC


  ! ***************************************************************************
  ! * operator C, C: E -> F, represents the multiplication by a sparse matrix
  ! * with non-zero values equal to +/- 1.
  ! * It acts on a complex vector defined on edges,
  ! * producing the respective vector defined on faces.
  ! * It is an implementation of the following computations:
  ! $e_\phi(i,j,k)=h_\th(i,j,k) + h_r(i,j+1,k) - h_\th(i,j,k+1) - h_r(i,j,k)$
  ! $e_\th(i,j,k)= h_r(i,j,k) + h_\phi(i,j,k+1) - h_r(i+1,j,k) - h_\phi(i,j,k)$
  ! $e_r(i,j,k)=h_\phi(i,j,k) + h_\th(i+1,j,k) - h_\phi(i,j+1,k) - h_\th(i,j,k)$
  ! * where h is the input vector (normally, edge lengths time the magnetic field
  ! * defined on edges).

  subroutine operatorC(vecE,vecF,grid)

    implicit none
	type (cvector), intent(inout)				 :: vecE
	type (cvector), intent(inout)					 :: vecF
	type (grid_t), intent(in)				 :: grid
	logical										 :: verbose
    integer										 :: i,j,k,istat
	integer										 :: nx,ny,nz
	real(8), dimension(:), allocatable			 :: x,y,z

    if (.not.vecE%allocated) then

 	  write(0,*) 'Error: (operatorC) vecE not allocated'
 	  stop

	else if (vecE%gridType /= EDGE) then

 	  write(0,*) 'Error: (operatorC) input vector is not defined on edges'
 	  stop

    endif

    if (.not.vecF%allocated) then

 	  call create_cvector(grid,vecF,FACE)

	else if (vecF%gridType /= FACE) then

 	  write(0,*) 'Error: (operatorC) output vector is not defined on faces'
 	  stop

    endif

	nx = grid%nx
	ny = grid%ny
	nz = grid%nz

	allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
	x = grid%x
	y = grid%y
	z = grid%z

	! Set the undefined values to zero just in case they are non-zero
	! (which should never happen! - unless a test vector is used)
	!vecE%x(:,1,:) = C_ZERO
	!vecE%x(:,ny+1,:) = C_ZERO

	! Set the repetitious values in case they are not already correct
	!vecE%y(nx+1,:,:) = vecE%y(1,:,:)
	!vecE%z(nx+1,:,:) = vecE%z(1,:,:)
	call validate_cvector(vecE,verbose)


  ! $e_\phi(i,j,k)=h_\th(i,j,k) + h_r(i,j+1,k) - h_\th(i,j,k+1) - h_r(i,j,k)$
  ! $e_\th(i,j,k)= h_r(i,j,k) + h_\phi(i,j,k+1) - h_r(i+1,j,k) - h_\phi(i,j,k)$
  ! $e_r(i,j,k)=h_\phi(i,j,k) + h_\th(i+1,j,k) - h_\phi(i,j+1,k) - h_\th(i,j,k)$

	! Boundary values of h at k=1 and k=nz+1 are fixed during initialization;
	! loop over radii; do not loop over k=nz+1 (fields undefined at nz+2).
    do k=1,nz
	  ! loop over co-latitudes from the North pole to the South
      do j=1,ny
		! loop over longitudes including i=nx (nx+1'th fields same as i=1)
        do i=1,nx

		  vecF%x(i,j,k) = vecE%y(i,j,k) + vecE%z(i,j+1,k) &
						  - vecE%y(i,j,k+1) - vecE%z(i,j,k)

		  vecF%y(i,j,k) = vecE%z(i,j,k) + vecE%x(i,j,k+1) &
						  - vecE%z(i+1,j,k) - vecE%x(i,j,k)

		  vecF%z(i,j,k) = vecE%x(i,j,k) + vecE%y(i+1,j,k) &
						  - vecE%x(i,j+1,k) - vecE%y(i,j,k)

         end do
      end do
    end do

	! x-component at nx+1 (zero longitude)
	vecF%x(nx+1,:,:) = vecF%x(1,:,:)

	! y-component at ny+1 (South pole): undefined
	vecF%y(:,ny+1,:) = C_ZERO

	! z-component at nz+1 (lower domain boundary)
	k=nz+1
      do j=1,ny
        do i=1,nx
		  vecF%z(i,j,k) = vecE%x(i,j,k) + vecE%y(i+1,j,k) &
						  - vecE%x(i,j+1,k) - vecE%y(i,j,k)
         end do
      end do

	! z-component at North pole: undefined
	!vecF%z(:,1,:) = C_ZERO

	call validate_cvector(vecF)

    deallocate(x,y,z)

  end subroutine operatorC	! C


  ! ***************************************************************************
  ! * operator Ct, Ct: F -> E, maps onto center edge of the "paddle wheel".
  ! * It acts on a complex vector defined on faces,
  ! * producing the respective vector defined on edges.
  ! * It is an implementation of the following computations:
  ! $h_\phi(i,j,k)=e_r(i,j,k) - e_\th(i,j,k) - e_r(i,j-1,k) + e_\th(i,j,k-1)$
  ! $h_\th(i,j,k)= e_\phi(i,j,k) - e_r(i,j,k) - e_\phi(i,j,k-1) + e_r(i-1,j,k)$
  ! $h_r(i,j,k)=e_\th(i,j,k) - e_\phi(i,j,k) - e_\th(i-1,j,k) + e_\phi(i,j-1,k)$
  ! * where e is the input vector defined on faces.
  ! * This operator is not required in the Jacobian computations. Currently
  ! * used for testing purposes only.

  subroutine operatorCt(vecF,vecE,grid)

    implicit none
	type (cvector), intent(inout)				 :: vecF
	type (cvector), intent(inout)					 :: vecE
	type (grid_t), intent(in)				 :: grid
	logical										 :: verbose
    integer										 :: i,j,k,istat
	integer										 :: nx,ny,nz
	complex(8)									 :: total
	real(8), dimension(:), allocatable			 :: x,y,z

    if (.not.vecF%allocated) then

 	  write(0,*) 'Error: (operatorCt) vecF not allocated'
 	  stop

	else if (vecF%gridType /= FACE) then

 	  write(0,*) 'Error: (operatorCt) input vector is not defined on faces'
 	  stop

    endif

    if (.not.vecE%allocated) then

 	  call create_cvector(grid,vecE,EDGE)

	else if (vecE%gridType /= EDGE) then

 	  write(0,*) 'Error: (operatorCt) output vector is not defined on edges'
 	  stop

    endif

	nx = grid%nx
	ny = grid%ny
	nz = grid%nz

	allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
	x = grid%x
	y = grid%y
	z = grid%z

	! Set the undefined values to zero just in case they are non-zero
	! (which should never happen! - unless a test vector is used)
	!vecF%y(:,1,:) = C_ZERO
	!vecF%y(:,ny+1,:) = C_ZERO

	! Set the repetitious values in case they are not already correct
	!vecF%x(nx+1,:,:) = vecF%x(1,:,:)
	call validate_cvector(vecF)

	! x-component

  ! $h_\phi(i,j,k)=e_r(i,j,k) - e_\th(i,j,k) - e_r(i,j-1,k) + e_\th(i,j,k-1)$
  ! $h_\th(i,j,k)= e_\phi(i,j,k) - e_r(i,j,k) - e_\phi(i,j,k-1) + e_r(i-1,j,k)$
  ! $h_r(i,j,k)=e_\th(i,j,k) - e_\phi(i,j,k) - e_\th(i-1,j,k) + e_\phi(i,j-1,k)$

	! loop over radii; do not loop over k=1 (fields undefined at k-1).
    do k=2,nz
	  ! loop over co-latitudes from the North pole to the South
      do j=2,ny
		! zero longitude (nx+1'th fields same as i=1)
        do i=1,nx
		  vecE%x(i,j,k) = vecF%z(i,j,k) + vecF%y(i,j,k-1) &
						- vecF%z(i,j-1,k) - vecF%y(i,j,k)
		end do
	  end do
	end do
	k = 1
	  do j=2,ny
		do i=1,nx
		  vecE%x(i,j,k) = vecF%z(i,j,k) &
						- vecF%z(i,j-1,k) - vecF%y(i,j,k)
		end do
	  end do
	k = nz+1
	  do j=2,ny
		do i=1,nx
		  vecE%x(i,j,k) = vecF%z(i,j,k) + vecF%y(i,j,k-1) &
						- vecF%z(i,j-1,k)
		end do
	  end do
	vecE%x(:,1,:) = C_ZERO	!undefined
	vecE%x(:,ny+1,:) = C_ZERO	!undefined

	! y-component

    do k=2,nz
      do j=1,ny
		vecE%y(1,j,k) = - vecF%z(1,j,k) + vecF%x(1,j,k) &
					  - vecF%x(1,j,k-1)
        do i=2,nx
		  vecE%y(i,j,k) = - vecF%z(i,j,k) + vecF%x(i,j,k) &
						+ vecF%z(i-1,j,k) - vecF%x(i,j,k-1)
		end do
		vecE%y(nx+1,j,k) = vecF%z(nx,j,k)
	  end do
	end do
	k = 1
	  do j=1,ny
		vecE%y(1,j,k) = - vecF%z(1,j,k) + vecF%x(1,j,k)
		do i=2,nx
		  vecE%y(i,j,k) = - vecF%z(i,j,k) + vecF%x(i,j,k) &
						+ vecF%z(i-1,j,k)
		end do
		vecE%y(nx+1,j,k) = vecF%z(nx,j,k)
	  end do
	k = nz+1
	  do j=1,ny
		vecE%y(1,j,k) = - vecF%z(1,j,k) - vecF%x(1,j,k-1)
		do i=2,nx
		  vecE%y(i,j,k) = - vecF%z(i,j,k) &
						+ vecF%z(i-1,j,k) - vecF%x(i,j,k-1)
		end do
		vecE%y(nx+1,j,k) = vecF%z(nx,j,k)
	  end do

	! z-component

    do j=2,ny
      do k=1,nz
		vecE%z(1,j,k) = vecF%y(1,j,k) + vecF%x(1,j-1,k) &
					   - vecF%x(1,j,k)
        do i=2,nx
		  vecE%z(i,j,k) = vecF%y(i,j,k) + vecF%x(i,j-1,k) &
						- vecF%y(i-1,j,k) - vecF%x(i,j,k)
		end do
		vecE%z(nx+1,j,k) = - vecF%y(nx,j,k)
	  end do
	end do
	j = 1
	  do k=1,nz
		do i=1,nx
		  vecE%z(i,j,k) = - vecF%x(i,j,k)
		end do
		vecE%z(nx+1,j,k) = C_ZERO
	  end do
	j = ny+1
	  do k=1,nz
		do i=1,nx
		  vecE%z(i,j,k) = vecF%x(i,j-1,k)
		end do
		vecE%z(nx+1,j,k) = C_ZERO
	  end do

	! Map back from full set of edges to unique edges

	do k=1,nz
	  vecE%z(1,1,k) = sum(vecE%z(:,1,k))
	  vecE%z(1,2:ny,k) = vecE%z(1,2:ny,k) + vecE%z(nx+1,2:ny,k)
	  vecE%z(1,ny+1,k) = sum(vecE%z(:,ny+1,k))
	end do

	do k=1,nz+1
	  vecE%y(1,:,k) = vecE%y(1,:,k) + vecE%y(nx+1,:,k)
	end do

	call validate_cvector(vecE)


    deallocate(x,y,z)

  end subroutine operatorCt	! This operator is the transpose of C


  ! ***************************************************************************
  ! * operator L is the mapping from cells onto faces, L: G -> F
  ! * This operator is the implementation of the following mapping:
  ! * $\rho -> (1/2*S)({l^+ \rho^+} + {l^- \rho^-})$

  subroutine operatorL(resist,vecF,grid)

    implicit none
	type (rscalar), intent(in)			         :: resist	 !(nx,ny,nz)
	type (rvector), intent(inout)				 :: vecF
	type (grid_t), intent(in)				     :: grid
    real(8)										 :: lm,lp,S
	integer										 :: i,j,k,istat
	integer										 :: nx,ny,nz
	real(8), dimension(:), allocatable			 :: x,y,z

    if (.not.vecF%allocated) then

 	  call create_rvector(grid,vecF,FACE)

	else if (vecF%gridType /= FACE) then

 	  write(0,*) 'Error: (operatorL) output vector is not defined on faces'
 	  stop

	else

	  vecF%x = R_ZERO
	  vecF%y = R_ZERO
	  vecF%z = R_ZERO

    endif

	nx = grid%nx
	ny = grid%ny
	nz = grid%nz

	allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
	x = grid%x
	y = grid%y
	z = grid%z


    do k=1,nz
	  ! loop over co-latitudes from the North pole to the South
      do j=1,ny
		! loop over longitudes including i=nx (nx+1'th fields same as i=1)
        do i=1,nx
		  call leng_xijk2_shifted(i,j,k,x,y,z,lp)
		  if (i==1) then
			call leng_xijk2_shifted(nx,j,k,x,y,z,lm)
		  else
			call leng_xijk2_shifted(i-1,j,k,x,y,z,lm)
		  end if
          call area_sijk(j,k,y,z,S)

		  if (i==1) then
			vecF%x(i,j,k)=(lp*resist%v(i,j,k)+lm*resist%v(nx,j,k))/(2*S)
		  else
			vecF%x(i,j,k)=(lp*resist%v(i,j,k)+lm*resist%v(i-1,j,k))/(2*S)
		  end if

		  call leng_yijk2_shifted(j,k,y,z,lp)
		  if (j==1) then
		  else
			call leng_yijk2_shifted(j-1,k,y,z,lm)
		  end if
          call area_sjki(i,j,k,x,y,z,S)

		  if (j==1) then
			vecF%y(i,j,k)=R_ZERO  !undefined
		  else
			vecF%y(i,j,k)=(lp*resist%v(i,j,k)+lm*resist%v(i,j-1,k))/(2*S)
		  end if

		  call leng_zijk2_shifted(k,z,lp)
		  if (k==1) then
			lm=lp
		  else
			call leng_zijk2_shifted(k-1,z,lm)
		  end if
          call area_skij(i,j,k,x,y,z,S)

		  if (k==1) then
			vecF%z(i,j,k)=(lp*resist%v(i,j,k))/(2*S)
		  else
			vecF%z(i,j,k)=(lp*resist%v(i,j,k)+lm*resist%v(i,j,k-1))/(2*S)
		  end if

		end do
	  end do
	end do

	vecF%x(nx+1,:,:) = vecF%x(1,:,:)
	vecF%y(:,ny+1,:) = R_ZERO !undefined

	call leng_zijk2_shifted(nz,z,lm)
    do j=1,ny
      do i=1,nx
        call area_skij(i,j,nz+1,x,y,z,S)
		vecF%z(i,j,nz+1) = lm*resist%v(i,j,nz)/(2*S)
	  end do
	end do

	deallocate(x,y,z)

  end subroutine operatorL	! L


  ! ***************************************************************************
  ! * operator Lt is the mapping from faces onto cells, L^t: F -> G

  subroutine operatorLt(resist,vecF,grid)

    implicit none
	type (rscalar), intent(inout)		         :: resist	 !(nx,ny,nz)
	type (rvector), intent(inout)				 :: vecF
	type (grid_t), intent(in)				 :: grid
    real(8)										 :: l,S,Sp
	integer										 :: i,j,k,istat
	integer										 :: nx,ny,nz
	real(8), dimension(:), allocatable			 :: x,y,z

    if (.not.vecF%allocated) then

 	  write(0,*) 'Error: (operatorLt) input vector is not allocated'
 	  stop

	else if (vecF%gridType /= FACE) then

 	  write(0,*) 'Error: (operatorLt) input vector is not defined on faces'
 	  stop

    endif


	nx = grid%nx
	ny = grid%ny
	nz = grid%nz

	allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
	x = grid%x
	y = grid%y
	z = grid%z

    ! Allocate the resistivity vector
	! resist = R_ZERO
    call create_rscalar(grid,resist,CENTER)

	! making sure the input vector makes scientific sense...
	vecF%x(nx+1,:,:) = vecF%x(1,:,:)
	vecF%y(:,ny+1,:) = R_ZERO !undefined

    do k=1,nz
	  ! loop over co-latitudes from the North pole to the South
      do j=1,ny
		! loop over longitudes including i=nx (nx+1'th fields same as i=1)
        do i=1,nx
		  call leng_xijk2_shifted(i,j,k,x,y,z,l)
          call area_sijk(j,k,y,z,S)
		  Sp=S
		  resist%v(i,j,k)=vecF%x(i,j,k)*l/(2*S)+vecF%x(i+1,j,k)*l/(2*Sp)

		  call leng_yijk2_shifted(j,k,y,z,l)
          call area_sjki(i,j,k,x,y,z,S)
          call area_sjki(i,j+1,k,x,y,z,Sp)
		  if (j==1) then
			resist%v(i,j,k)=resist%v(i,j,k)+vecF%y(i,j+1,k)*l/(2*Sp)
		  else if (j==ny) then
			resist%v(i,j,k)=resist%v(i,j,k)+vecF%y(i,j,k)*l/(2*S)
		  else
  			resist%v(i,j,k)=resist%v(i,j,k)+vecF%y(i,j,k)*l/(2*S)+vecF%y(i,j+1,k)*l/(2*Sp)
		  end if

		  call leng_zijk2_shifted(k,z,l)
          call area_skij(i,j,k,x,y,z,S)
          call area_skij(i,j,k+1,x,y,z,Sp)
		  resist%v(i,j,k)=resist%v(i,j,k)+vecF%z(i,j,k)*l/(2*S)+vecF%z(i,j,k+1)*l/(2*Sp)

		end do
	  end do
	end do

	! North pole
	j = 1
	do k=2,nz
	  do i=1,nx
	  end do
	end do

	! South pole
	j = ny
	do k=2,nz
	  do i=1,nx
	  end do
	end do

	! Uppermost grid layer
	k = 1
	do j=1,ny
	  do i=1,nx
		if (j==1) then
		else if (j==ny) then
		else
		end if
	  end do
	end do

	! Lowermost grid layer
	k = nz
	do j=1,ny
	  do i=1,nx
		if (j==1) then
		else if (j==ny) then
		else
		end if
	  end do
	end do


	deallocate(x,y,z)

  end subroutine operatorLt	! This operator is the transpose of L


  ! ***************************************************************************
  ! * operator P is the mapping from the model parameters to cell centres,
  ! * P: R^n -> G; by definition, \drho = P \da.
  ! * vector p_j is the j'th column of P, that corresponds to \pd{\rho}{a_j}.
  ! * p_j \In G.

  function vectorPj(m0,da,n,grid) result (dm)

    type (modelParam_t), intent(in)                 :: m0
    type (modelParam_t), intent(inout)				:: da
	integer, intent(in)								:: n
    type (modelParam_t)                             :: m  !current resistivity on the grid
    type (rscalar)                                  :: resist, resist0 !resist0 is fixed background
    type (rscalar)                                  :: dresist
    type (grid_t), intent(in)                       :: grid
	type (rscalar)									:: dm !(nx,ny,nz)
	integer											:: i,j,k,l,istat
	real(8)											:: value
	type (modelLayer_t), pointer						:: this_layer
	type (modelPoint_t)								:: point
	type (modelFunc_t)								:: func
	type (modelCoeff_t)								:: coeff

    ! Preliminary computations: only needed if arctan parametrization is involved
    resist0 = m0%rho0
    call mapToGrid(m0,resist)
    dresist = resist
    dresist%v = dtan((PI/2)*(resist%v/resist0%v - 1.0d0))

	! Create a zero-valued output in G
	call create_rscalar(grid,dm,CENTER)

	coeff = getCoeff_modelParam(da,n)

	if (coeff%frozen) then
	  return
	end if

	this_layer => coeff%L

	! Resistivity is constant at air layers, hence derivative is zero
	forall (i=1:grid%nx, j=1:grid%ny, k=1:grid%nzAir)
	  dm%v(i,j,k) = 0.0d0
	end forall

	do k=grid%nzAir+1,grid%nz

	  if (.not.in_layer(grid%r(k),coeff%L)) then
		cycle
	  end if

	  do i=1,grid%nx
		do j=1,grid%ny

		  point%phi   = (grid%ph(i) + grid%ph(i+1))/2
		  point%theta = (grid%th(j) + grid%th(j+1))/2
		  point%r     = (grid%r(k) + grid%r(k+1))/2

		  func = coeff%F
		  value = F_at_point(func,point)

		  if (this_layer%if_log) then
			dm%v(i,j,k) = value*resist%v(i,j,k)*log(10.)
          elseif (this_layer%if_tan) then
            dm%v(i,j,k) = (2.0d0/PI)*value*(resist0%v(i,j,k)/(1.0d0+dresist%v(i,j,k)**2))
          elseif (this_layer%if_exp) then
            dm%v(i,j,k) = value*(-resist0%v(i,j,k)*resist%v(i,j,k))
		  else
			dm%v(i,j,k) = value
		  end if

		end do
	  end do
	end do

	call deall_modelParam(m)
	call deall_rscalar(resist)
    call deall_rscalar(resist0)
    call deall_rscalar(dresist)
    call deall_rscalar(dm)

  end function vectorPj	! j'th column of P


  ! ***************************************************************************
  ! * operator P is the mapping from the model parameters to cell centres,
  ! * P: R^n -> G; by definition, \drho = P \da.
  ! * We define $\tau^l_{ijk}=\drho_{ijk}/\da_l$. Hence $\tau_l$ is the
  ! * l'th column of operator P, and therefore the l'th row of Pt.
  ! * Therefore, we obtain $\da_l=\sum_{i,j,k} \tau^l_{ijk}\drho_{ijk}$.
  ! * We compute $\tau_l$ for the following two cases:
  ! * 1)		 $\rho_{ijk}  = \sum_l a_l F_l(\phi,\th,r)$,
  ! * 2) $\log_10(\rho_{ijk}) = \sum_l a_l F_l(\phi,\th,r)$,
  ! * where (phi,th,r) is the (i,j,k)'th cell centre.
  ! * In the first case, $\tau^l_{ijk}=F_l(\phi,\th,r)$;
  ! * in the second case, $\tau^l_{ijk}=\rho_{ijk} F_l(\phi,\th,r) \log(10)$.
  ! * The rho, drho, dprm are global variables defined in module modeldef
  ! * Functions to substitute for F are declared in module paramfunc
  ! *
  ! * mapToGrid is the routine that generates the 3-D resistivity map on the
  ! * grid using the information stored in the model parametrization,
  ! * including the grid and background resistivity pointers.

  subroutine operatorP(da,dm,grid,m0)

    type (modelParam_t), intent(in)                 :: m0
    type (modelParam_t), intent(in)					:: da
	type (rscalar), intent(inout)						:: dm !(nx,ny,nz)
    type (modelParam_t)                              :: m  !current resistivity on the grid
    type (rscalar)                                   :: resist, resist0 !resist0 is fixed background
    type (grid_t), intent(in)                       :: grid
	integer											:: i,j,k,l,istat
	integer											:: iL,ip
	real(8)											:: value,coeff
	type (modelLayer_t), pointer						:: this_layer
	type (modelPoint_t)								:: point
	type (modelFunc_t)								:: func
    type (rscalar)                                  :: dresist !(nx,ny,nz)

    ! Preliminary computations: only needed if arctan parametrization is involved
    resist0 = m0%rho0
    call mapToGrid(m0,resist)
    dresist = resist
    dresist%v = dtan((PI/2)*(resist%v/resist0%v - 1.0d0))

	! Create the output resistivity model on the grid
	call create_rscalar(grid,dm,CENTER)

	! Resistivity is constant at air layers, hence derivative is zero
	forall (i=1:grid%nx, j=1:grid%ny, k=1:grid%nzAir)
	  dm%v(i,j,k) = 0.0d0
	end forall

	do k=grid%nzAir+1,grid%nz

	  ! Find current layer by locating the upper boundary of a cell
	  do l=1,m0%nL
		if (in_layer(grid%r(k),m0%L(l))) then
		  this_layer => m0%L(l)
		  exit
		end if
	  end do

	  iL = this_layer%num

	  do i=1,grid%nx
		do j=1,grid%ny

		  point%phi   = (grid%ph(i) + grid%ph(i+1))/2
		  point%theta = (grid%th(j) + grid%th(j+1))/2
		  point%r     = (grid%r(k) + grid%r(k+1))/2

		  dm%v(i,j,k) = 0.0d0

		  ! Sum up the coeffs * F_at_point in the given layer
		  value = 0.0d0
		  do ip=1,m0%nF

			if(da%c(iL,ip)%frozen) then
			  cycle
			end if
			coeff = da%c(iL,ip)%value
			func = da%F(ip)
			value = F_at_point(func,point)

			if (this_layer%if_log) then
			  dm%v(i,j,k) = dm%v(i,j,k) + value*coeff*resist%v(i,j,k)*log(10.)
            elseif (this_layer%if_tan) then
              dm%v(i,j,k) = dm%v(i,j,k) + (2.0d0/PI)*value*coeff*(resist0%v(i,j,k)/(1.0d0+dresist%v(i,j,k)**2))
            elseif (this_layer%if_exp) then
              dm%v(i,j,k) = dm%v(i,j,k) + value*(-resist0%v(i,j,k)*resist%v(i,j,k))
			else
			  dm%v(i,j,k) = dm%v(i,j,k) + value*coeff
			end if

		  end do

		end do
	  end do
	end do

    call deall_modelParam(m)
    call deall_rscalar(resist)
    call deall_rscalar(resist0)
    call deall_rscalar(dresist)

  end subroutine operatorP	! P


  ! ***************************************************************************
  ! * operator Pt is the mapping from cells to the model parameters,
  ! * P^t: G -> R^n
  ! * This is an implementation of \da^t * P^t = \drho^t
  ! * In the general case, P^t is a real matrix $n\times \mod{\mathbb{G}}}$
  ! * defined as the transpose of $P = \drho/\da$.
  ! * Its' element $\tau_l(ijk) = \drho(ijk)/\da_l$.
  ! *
  ! * Output vector size equals to the number of parametrization coefficients.
  ! * The values are only computed for variable parameters, otherwize they are
  ! * set to zero. This is done so that it would be easy to map back to the
  ! * original parametrization structure when we output and otherwise use these
  ! * values. It is easy to extract the components corresponding to the variable
  ! * parameters only:
  ! * count = 1
  ! * do index = 1,size(vecN)
  ! *   if (.not.param%p(index)%frozen) then
  ! *	  value(count) = vecN(index)
  ! *	  count = count + 1
  ! *	end if
  ! * end do

  subroutine operatorPt(dm,da,grid,m0)

    type (modelParam_t), intent(in)                 :: m0 !(nx,ny,nz)
	type (rscalar), intent(in)						:: dm !(nx,ny,nz)
    type (modelParam_t), intent(inout)				:: da
    type (grid_t), intent(in)                        :: grid
	type (modelParam_t)				                 :: m  !current resistivity on the grid
    type (rscalar)                                   :: resist, resist0 !resist0 is fixed background
	!real(8),dimension(:),intent(out)				 :: da  !ncoeff
	integer											 :: i,j,k,l
	integer											 :: iL,ip
	real(8)											 :: tau ! single entry of matrix P^t
    type (rscalar)                                   :: dresist !(nx,ny,nz)
	type (modelPoint_t)								 :: point
	type (modelFunc_t)								 :: func
	type (modelLayer_t)								 :: this_layer

!	if(size(dm) /= grid%nx * grid%ny * grid%nz) then
!		write(0, *) 'Error: (operatorPt) input vector should be defined at every grid cell'
!		stop
!	end if

!	if(size(da) /= param%np) then
!		write(0, *) 'Error: (operatorPt) output vector size should = # of parametrization coefficients'
!		stop
!	end if

!	if(size(da) /= count(.not.param%p(1:param%np)%frozen)) then
!		write(0, *) 'Error: (operatorPt) output vector size should = # of variable parameters'
!		stop
!	end if

    ! Preliminary computations: only needed if arctan parametrization is involved
    resist0 = m0%rho0
    call mapToGrid(m0,resist)
    dresist = resist
    dresist%v = dtan((PI/2)*(resist%v/resist0%v - 1.0d0))

	da = m0
	call zero(da)

	do iL = 1,m0%nL

	  this_layer = m0%L(iL)

	  do ip = 1,m0%nF

		! If this parameter is frozen, the derivative is zero, cycle
		if (m0%c(iL,ip)%frozen) then
		  cycle
		end if

		! Otherwise find the functional
		func = m0%c(iL,ip)%F

		! Going vertically down through the chosen layer
		do k=grid%nzAir+1,grid%nz

		  ! If grid radius is not in this layer, ignore
		  if(.not.in_layer(grid%r(k),this_layer)) then
			cycle
		  end if
		  ! Otherwise, compute an expression for each cell
		  do i=1,grid%nx
			do j=1,grid%ny
			  point%phi   = (grid%ph(i) + grid%ph(i+1))/2
			  point%theta = (grid%th(j) + grid%th(j+1))/2
			  point%r     = (grid%r(k) + grid%r(k+1))/2

			  tau = F_at_point(func,point)
			  if (this_layer%if_log) then
				tau = tau * resist%v(i,j,k) * log(10.0d0)
              elseif (this_layer%if_tan) then
                tau = tau * (2.0d0/PI) * (resist0%v(i,j,k)/(1.0d0+dresist%v(i,j,k)**2))
              elseif (this_layer%if_exp) then
                tau = tau * (-resist0%v(i,j,k)*resist%v(i,j,k))
			  end if
			  ! Add this expression to the output vector component
			  da%c(iL,ip)%value = da%c(iL,ip)%value + tau * dm%v(i,j,k)
			end do
		  end do

		end do
	  end do
	end do

    call deall_modelParam(m)
    call deall_rscalar(resist)
    call deall_rscalar(resist0)
    call deall_rscalar(dresist)
	da%zeroValued = .false.

  end subroutine operatorPt	! This operator is the transpose of P



end module jacobian
