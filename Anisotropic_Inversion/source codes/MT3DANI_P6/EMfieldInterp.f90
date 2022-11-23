! *****************************************************************************
module EMfieldInterp
  ! Generic data functionals for interpolation of electric and magnetic
  ! fields to an arbitrary point within the model domain for 3D staggered
  ! grid finite difference solutions

  use utilities
  use sg_sparse_vector
  use sg_scalar
  use ModelSpace

  implicit none

  public				:: HBinterpSetUp
  public				:: HEinterpSetUp, EfromHESetUp
  type(rscalar), pointer, dimension(:), private :: aRes,dxx,dxy,dxz,dyy,dyz,dzz
  
Contains


  ! **************************************************************************
  ! HEdgeCenter,magnetic field coeffcients in the sparse vector
  subroutine HBinterpSetUp(inGrid,x,xyz,LC)

    implicit none
    type(grid_t), target, intent(in)              :: inGrid
    real(kind=prec), dimension(3), intent(in)     :: x
    integer, intent(in)                         :: xyz
    type(sparsevecc), intent(inout)             :: LC

    ! Local Variables
    integer                                     :: i0,j0,k0,ii,n,m,p,ic,jc,kc
    integer, dimension(8)                       :: I,J,K
    integer                                     :: nxMax, nyMax, nzMax
    integer                                     :: status
    complex(kind=prec), dimension(8)    :: C
    real(kind=prec), dimension(3,2)     :: w
    real(kind=prec)                     :: wadd
    character (len = 80)                        :: gridType = ''

    integer, parameter		:: IX = 1, IY = 2, IZ = 3

    if(LC%allocated) then
       deallocate(LC%i,STAT=status)
       deallocate(LC%j, STAT=status)
       deallocate(LC%k, STAT=status)
       deallocate(LC%xyz, STAT=status)
       deallocate(LC%c, STAT=status)
       LC%gridType = ''
       LC%allocated = .false.
    endif

    ! maximum number of edge nodes
    nxMax = inGrid%nx+1
    nyMax = inGrid%ny+1
    nzMax = inGrid%nz+1
    i0 = minNode(x(1),inGrid%xEdge)
    j0 = minNode(x(2),inGrid%yEdge)
    k0 = minNode(x(3),inGrid%zEdge)
    if (xyz .eq. 1) then
       ! Evaluation of x component wrt magnetic vectors on cubic edges
       ic = i0
       i0 = minNode(x(1),inGrid%xCenter)
       ! maximum number of center nodes
       nxMax = nxMax -1
    elseif (xyz .eq. 2) then
       ! Evaluation of y component wrt magnetic vectors on cubic edges
       jc = j0
       j0 = minNode(x(2),inGrid%yCenter)
       ! maximum number of center nodes
       nyMax = nyMax -1
    elseif (xyz .eq. 3) then
       ! Evaluation of z component wrt magnetic vectors on cubic edges
       kc = k0
       k0 = minNode(x(3),inGrid%zCenter)
       ! maximum number of center nodes
       nzMax = nzMax -1
    else
       write(0,*) 'Error: component # out of range in EinterpSetUp'
    endif

    if((i0.gt.0).and.(i0.lt.nxMax)) then
       if(xyz.eq.1) then
          w(1,2) = (x(1) - inGrid%xCenter(i0))/(inGrid%delX(i0+1))
       else
          w(1,2) = (x(1) - inGrid%xEdge(i0))/(inGrid%dx(i0))
       endif
    elseif(i0.le.0) then
       w(1,2) = 1  
    else
       w(1,2) = 0
    endif

    if((j0.gt.0).and.(j0.lt.nyMax)) then
       if(xyz.eq.2) then
          w(2,2) = (x(2) - inGrid%yCenter(j0))/(inGrid%delY(j0+1))
       else
          w(2,2) = (x(2) - inGrid%yEdge(j0))/(inGrid%dy(j0))
       endif
    elseif(j0.le.0) then
       w(2,2) = 1
    else
       w(2,2) = 0
    endif

    if((k0.gt.0).and.(k0.lt.nzMax)) then
       if(xyz.eq.3) then
          w(3,2) = (x(3) - inGrid%zCenter(k0))/(inGrid%delZ(k0+1))
       else
          w(3,2) = (x(3) - inGrid%zEdge(k0))/(inGrid%dz(k0))
       endif
    elseif(k0.le.0) then
       w(3,2) = 1
    else
       w(3,2) = 0
    endif

    w(1,1) = 1-w(1,2)
    w(2,1) = 1-w(2,2)
    w(3,1) = 1-w(3,2)

    ii = 0
    do n = 1,2
       do m = 1,2
          do p = 1,2
             wadd = w(1,n)*w(2,m)*w(3,p)
             if(wadd .gt. 0) then
                ii = ii + 1
                I(ii) = i0+n-1
                J(ii) = j0+m-1
                K(ii) = k0+p-1
                C(ii) = wadd * MU_0 !!!!!!!!!!!!!!!!!!
             endif
          enddo
       enddo
    enddo

    gridType = EDGE
    Call create_sparsevecc(ii,LC,gridType)
    ! I, J, K and C inside LC have a pointer.
    ! In some compiler (e.g. gfortran), it is not allowed
    ! to copy implicitly a dimensioned array (I) into a pointer.  
    ! Thus, copy explicitly
     do n=1,ii
      LC%i(n)=I(n)
      LC%j(n)=J(n)
      LC%k(n)=K(n)
      LC%c(n)=C(n)
     end do
    !   assuming xyz will be assigned to all elements of LC%xyz
    LC%xyz = xyz
  end subroutine HBinterpSetUp


  ! **************************************************************************
  ! HEdgeCenter, electric field coefficients in the sparse vector
  ! electric field is assumed to be on the earth surface, and no z component
  subroutine HEinterpSetUp(inGrid,x,xyz,LC)

    implicit none
    type (grid_t), intent(in)    		:: inGrid
    real (kind=prec), dimension(3), intent(in) 	:: x
    integer, intent(in)        			:: xyz
    type (sparsevecc), intent(inout) 		:: LC   

    integer 					:: i0,j0,k0,ii,n,m,p
    integer                                     :: nxMax, nyMax, nzMax
    integer, dimension(8)			:: I,J,K
    integer					:: status
    real (kind=prec), dimension(2,2)	:: w
    real (kind=prec)			:: wadd
    complex (kind=prec), dimension(8)  	:: C
    character (len=80)				:: gridType = ''

    if(LC%allocated) then
       deallocate(LC%i, STAT=status)
       deallocate(LC%j, STAT=status)
       deallocate(LC%k, STAT=status)
       deallocate(LC%xyz, STAT=status)
       deallocate(LC%c, STAT=status)
       LC%gridType = ''
       LC%allocated = .false.
    endif

    ! maximum number of center nodes
    nxMax = inGrid%nx
    nyMax = inGrid%ny
    nzMax = inGrid%nz

    if (xyz .eq. 1 .or. xyz .eq. 2) then
       ! Evaluation of x component wrt magnetic vectors on cubic faces
       i0 = minNode(x(1),inGrid%xCenter)
       j0 = minNode(x(2),inGrid%yCenter)
       k0 = minNode(x(3),inGrid%zEdge)
    else
       write(0,*) 'Error: component # out of range in HEinterpSetUp'
    endif

    if((i0.gt.0).and.(i0.lt.nxMax)) then
       w(1,2) = (x(1) - inGrid%xCenter(i0))/(inGrid%delX(i0+1))
    elseif(i0.le.0) then
       w(1,2) = 1
    else
       w(1,2) = 0
    endif

    if((j0.gt.0).and.(j0.lt.nyMax)) then
       w(2,2) = (x(2) - inGrid%yCenter(j0))/(inGrid%delY(j0+1))
    elseif(j0.le.0) then
       w(2,2) = 1
    else
       w(2,2) = 0
    endif


    w(1,1) = 1-w(1,2)
    w(2,1) = 1-w(2,2)

    ii = 0
    do n = 1,2
       do m = 1,2
          wadd = w(1,n)*w(2,m)
          if(wadd .gt. 0) then
             ii = ii + 1
             I(ii) = i0+n-1
             J(ii) = j0+m-1
             K(ii) = k0
             C(ii) = wadd
          endif
       enddo
    enddo

    ! we are dealing with electric fields (therefore, gridtype = FACE)
    gridType = EDGE
    Call create_sparsevecc(ii,LC,gridType)
    ! The same as for H
     do n=1,ii
      LC%i(n)=I(n)
      LC%j(n)=J(n)
      LC%k(n)=K(n)
      LC%c(n)=C(n)
     end do    
    !   assuming xyz will be assigned to all elements of LC%xyz
    LC%xyz = xyz

  end subroutine HEinterpSetUp


  ! determine the limitation of model parameters' cycle
  subroutine M_Index(inGrid,x,ilm,jlm)
    implicit none
    type (grid_t), target, intent(in) 		:: inGrid
    real(kind=prec), dimension(3), intent(in)	:: x
    integer :: ilm(2),jlm(2)
    ! local variables
    integer :: imin,imax,jmin,jmax
    integer :: nxMax,nyMax
    integer :: ii,xyz,status
    type (sparsevecc) :: LCE
    
    nxMax = inGrid%nx
    nyMax = inGrid%ny
    
    
    xyz = 1 ! Esx
    Call HEinterpSetUp(inGrid,x,xyz,LCE)
    
    imin = LCE%i(1)
    imax = LCE%i(1)
    jmin = LCE%j(1)
    jmax = LCE%j(1)
    do ii = 1,LCE%nCoeff
       if(imin.gt.LCE%i(ii)) imin = LCE%i(ii)
       if(imax.lt.LCE%i(ii)) imax = LCE%i(ii)
       if(jmin.gt.LCE%j(ii)) jmin = LCE%j(ii)
       if(jmax.lt.LCE%j(ii)) jmax = LCE%j(ii)
    enddo    
    
    xyz = 2 ! Esy
    Call HEinterpSetUp(inGrid,x,xyz,LCE)
    
    do ii = 1,LCE%nCoeff
       if(imin.gt.LCE%i(ii)) imin = LCE%i(ii)
       if(imax.lt.LCE%i(ii)) imax = LCE%i(ii)
       if(jmin.gt.LCE%j(ii)) jmin = LCE%j(ii)
       if(jmax.lt.LCE%j(ii)) jmax = LCE%j(ii)
    enddo    
       
    
    if(imin.eq.1) then
      ilm(1) = imin
    else
      ilm(1) = imin - 1
    endif
    if(imax.eq.nxMax) then
      ilm(2) = imax
    else
      ilm(2) = imax + 1
    endif
 
    if(jmin.eq.1) then
      jlm(1) = jmin
    else
      jlm(1) = jmin - 1
    endif
    if(jmax.eq.nyMax) then
      jlm(2) = jmax
    else
      jlm(2) = jmax + 1
    endif  
    
    ! release memory
    deallocate(LCE%i, STAT=status)
    deallocate(LCE%j, STAT=status)
    deallocate(LCE%k, STAT=status)
    deallocate(LCE%xyz, STAT=status)
    deallocate(LCE%c, STAT=status)
    LCE%gridType = ''
    LCE%allocated = .false.
  
  end subroutine M_Index

  ! **************************************************************************
  !  electric field from magnetic field in a sparse vector data structures
  subroutine EfromHESetUp(inGrid,omega,x,xyz,LC,rho0,iPar,mx)

    implicit none
    type (grid_t), target, intent(in) 		:: inGrid
    real(kind=prec), 	intent(in)	:: omega
    real(kind=prec), dimension(3), intent(in)	:: x
    integer, intent(in)              		:: xyz
    type (sparsevecc), intent(inout) 		:: LC
    type(modelParam_t), intent(in) :: rho0
    integer, intent(in),optional  :: iPar,mx(3)

    ! local variables
    integer 					    :: ii
    complex (kind=prec)   	       	:: c1,c2
    type (sparsevecc)				:: LCE, LCEH, LCtemp
    character (len=80)				:: gridType = ''
    integer					:: status
    

    if(LC%allocated) then 
       deallocate(LC%i, STAT=status)
       deallocate(LC%j, STAT=status)
       deallocate(LC%k, STAT=status)
       deallocate(LC%xyz, STAT=status)
       deallocate(LC%c, STAT=status)
       LC%gridType = ''
       LC%allocated = .false.
    endif

    gridType = EDGE
    Call HEinterpSetUp(inGrid,x,xyz,LCE) ! weight coefficients

    !print*,'LCE%nCoeff',LCE%nCoeff !debug

    do ii = 1,LCE%nCoeff

       if(ii.eq.1) then
         
         if(present(mx)) then
           Call E_HinterpSetUp(ii,omega,inGrid,LCE,LCEH,rho0,iPar,mx)
         else
           Call E_HinterpSetUp(ii,omega,inGrid,LCE,LCEH,rho0)
         endif
		 LC = LCEH ! i,j,k,xyz
		 LC%c = LCE%c(ii) * LCEH%c

       else

         if(present(mx)) then
           Call E_HinterpSetUp(ii,omega,inGrid,LCE,LCEH,rho0,iPar,mx)
         else
           Call E_HinterpSetUp(ii,omega,inGrid,LCE,LCEH,rho0)
         endif
		 LCtemp = LC
         c1 = C_ONE
         c2 = LCE%c(ii)
         Call linComb_sparsevecc(LCtemp,c1,LCEH,c2,LC)

       endif

    enddo
    
    !pause !debug

    deallocate(LCE%i, STAT=status)
    deallocate(LCE%j, STAT=status)
    deallocate(LCE%k, STAT=status)
    deallocate(LCE%xyz, STAT=status)
    deallocate(LCE%c, STAT=status)
    LCE%gridType = ''
    LCE%allocated = .false.

    deallocate(LCEH%i, STAT=status)
    deallocate(LCEH%j, STAT=status)
    deallocate(LCEH%k, STAT=status)
    deallocate(LCEH%xyz, STAT=status)
    deallocate(LCEH%c, STAT=status)
    LCEH%gridType = ''
    LCEH%allocated = .false.
    
    deallocate(LCtemp%i, STAT=status)
    deallocate(LCtemp%j, STAT=status)
    deallocate(LCtemp%k, STAT=status)
    deallocate(LCtemp%xyz, STAT=status)
    deallocate(LCtemp%c, STAT=status)
    LCtemp%gridType = ''
    LCtemp%allocated = .false.

  end subroutine EfromHESetUp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine E_HinterpSetUp(ii,omega,inGrid,LCE,LCEH,rho0,iPar,mx)
    
    implicit none
	  integer ii
    real(kind=prec), 	intent(in)	        :: omega
    type (grid_t), target, intent(in) 		:: inGrid
	  type (sparsevecc), intent(in) 		    :: LCE
    type (sparsevecc), intent(inout) 		:: LCEH
    type(modelParam_t), intent(in) :: rho0
    integer, intent(in), optional  :: iPar,mx(3)
	  ! local variables
	  integer jj
    integer					:: num
    integer, dimension(4)			:: I,J,K,AXES
    complex (kind=prec), dimension(4)  	:: C
	  complex (kind=prec)   	       	:: c1,c2,i_omega_mu0
	  type (sparsevecc)  :: LCEJ,LCtemp,LCJH,LCH
	  integer					        :: status
    character (len=80)				:: gridType = ''

    num = 4
    gridType = EDGE
    if(present(mx)) then
      Call SensE_JinterpSetUp(ii,inGrid,LCE,LCEJ,iPar,mx) 
    else
      Call E_JinterpSetUp(ii,inGrid,LCE,LCEJ,rho0) 
    endif
	Call create_sparsevecc(num,LCEH,gridType)
	
	!print*,'LCEJ%nCoeff',LCEJ%nCoeff ! debug

    do jj = 1,LCEJ%nCoeff

       if(LCEJ%xyz(jj).eq.1) then  ! Jx

          AXES(1) = 2  ! Hy(i,j,k)
          AXES(2) = 2  ! Hy(i,j,k+1)
          AXES(3) = 3  ! Hz(i,j,k)
          AXES(4) = 3  ! Hz(i,j+1,k)
          I(1) = LCEJ%i(jj)
          I(2) = LCEJ%i(jj)
          I(3) = LCEJ%i(jj)
          I(4) = LCEJ%i(jj)
          J(1) = LCEJ%j(jj)
          J(2) = LCEJ%j(jj)
          J(3) = LCEJ%j(jj)
          J(4) = LCEJ%j(jj)+1
          K(1) = LCEJ%k(jj)
          K(2) = LCEJ%k(jj)+1
          K(3) = LCEJ%k(jj)
          K(4) = LCEJ%k(jj)
          C(1) = 1./inGrid%dz(K(1))
	      C(2) = -1./inGrid%dz(K(1))
          C(3) = -1./inGrid%dy(J(1))
          C(4) = 1./inGrid%dy(J(1))

       elseif(LCEJ%xyz(jj).eq.2) then ! Jy

          AXES(1) = 3 ! Hz(i,j,k)
          AXES(2) = 3 ! Hz(i+1,j,k)
          AXES(3) = 1 ! Hx(i,j,k)
          AXES(4) = 1 ! Hx(i,j,k+1)
          J(1) = LCEJ%j(jj)
          J(2) = LCEJ%j(jj)
          J(3) = LCEJ%j(jj)
          J(4) = LCEJ%j(jj)
          K(1) = LCEJ%k(jj)
          K(2) = LCEJ%k(jj)
          K(3) = LCEJ%k(jj)
          K(4) = LCEJ%k(jj)+1
          I(1) = LCEJ%i(jj)
          I(2) = LCEJ%i(jj)+1
          I(3) = LCEJ%i(jj)
          I(4) = LCEJ%i(jj)
          C(1) = 1./inGrid%dx(I(1))
          C(2) = -1./inGrid%dx(I(1))
          C(3) = -1./inGrid%dz(K(1))
          C(4) = 1./inGrid%dz(K(1))

       else ! Jz

          AXES(1) = 1 ! Hx(i,j,k)
          AXES(2) = 1 ! Hx(i,j+1,k)
          AXES(3) = 2 ! Hy(i,j,k)
          AXES(4) = 2 ! Hy(i+1,j,k)
          K(1) = LCEJ%k(jj)
          K(2) = LCEJ%k(jj)
          K(3) = LCEJ%k(jj)
          K(4) = LCEJ%k(jj)
          I(1) = LCEJ%i(jj)
          I(2) = LCEJ%i(jj)
          I(3) = LCEJ%i(jj)
          I(4) = LCEJ%i(jj)+1
          J(1) = LCEJ%j(jj)
          J(2) = LCEJ%j(jj)+1
          J(3) = LCEJ%j(jj)
          J(4) = LCEJ%j(jj)
          C(1) = 1./inGrid%dy(J(1))
          C(2) = -1./inGrid%dy(J(1))
          C(3) = -1./inGrid%dx(I(1))
          C(4) = 1./inGrid%dx(I(1))

       endif


       if(jj.eq.1) then

          Call create_sparsevecc(num,LCJH,gridType)
          LCEH%i = I
          LCEH%j = J
          LCEH%k = K
          LCEH%xyz = AXES
          LCEH%c= C*LCEJ%c(jj)
       else
          ! add coefficients for next face
          ! first store cumulative sum so far in LCtemp
          ! Call create_sparsevecc(LCEH%nCoeff,LCtemp,gridType)
	      LCtemp = LCEH
          LCJH%i = I
          LCJH%j = J
          LCJH%k = K
          LCJH%xyz = AXES
          LCJH%c = C

          c1 = C_ONE
          c2 = LCEJ%c(jj)
          Call linComb_sparsevecc(LCtemp,c1,LCJH,c2,LCEH)

       endif

    end do

    deallocate(LCEJ%i, STAT=status)
    deallocate(LCEJ%j, STAT=status)
    deallocate(LCEJ%k, STAT=status)
    deallocate(LCEJ%xyz, STAT=status)
    deallocate(LCEJ%c, STAT=status)
    LCEJ%gridType = ''
    LCEJ%allocated = .false.

    if(.not.present(mx)) then
        i_omega_mu0 = ISIGN * cmplx(0.0 ,omega*MU_0,kind=prec)

	    num = 2
        Call create_sparsevecc(num,LCH,gridType)
	    ! the remaining magnetic fields mutiplied by i*omega*mu_0
	    if(LCE%xyz(ii).eq.1) then ! Esx
	      ! Hy(i,j,k)
	      LCH%i(1) = LCE%I(ii)
	      LCH%j(1) = LCE%J(ii)
	      LCH%k(1) = LCE%K(ii)
	      LCH%xyz(1) = 2
	      LCH%c(1) = 0.25d0*i_omega_mu0*inGrid%dz(LCH%k(1))
	      ! Hy(i+1,j,k)
	      LCH%i(2) = LCE%I(ii) + 1
	      LCH%j(2) = LCE%J(ii)
	      LCH%k(2) = LCE%K(ii)
	      LCH%xyz(2) = 2
	      LCH%c(2) = 0.25d0*i_omega_mu0*inGrid%dz(LCH%k(2))
      elseif(LCE%xyz(ii).eq.2) then ! Esy
	      ! Hx(i,j,k)
	      LCH%i(1) = LCE%I(ii)
	      LCH%j(1) = LCE%J(ii)
	      LCH%k(1) = LCE%K(ii)
	      LCH%xyz(1) = 1
	      LCH%c(1) = -0.25d0*i_omega_mu0*inGrid%dz(LCH%k(1))
	      ! Hx(i,j+1,k)
	      LCH%i(2) = LCE%I(ii)
	      LCH%j(2) = LCE%J(ii) + 1
	      LCH%k(2) = LCE%K(ii)
	      LCH%xyz(2) = 1
	      LCH%c(2) = -0.25d0*i_omega_mu0*inGrid%dz(LCH%k(2))
      endif
        
        ! final construction
        LCtemp = LCEH
        c1 = C_ONE
        c2 = C_ONE
        Call linComb_sparsevecc(LCtemp,c1,LCH,c2,LCEH)
        
        deallocate(LCH%i, STAT=status)
        deallocate(LCH%j, STAT=status)
        deallocate(LCH%k, STAT=status)
        deallocate(LCH%xyz, STAT=status)
        deallocate(LCH%c, STAT=status)
        LCH%gridType = ''
        LCH%allocated = .false.
            
    endif
  
    deallocate(LCJH%i, STAT=status)
    deallocate(LCJH%j, STAT=status)
    deallocate(LCJH%k, STAT=status)
    deallocate(LCJH%xyz, STAT=status)
    deallocate(LCJH%c, STAT=status)
    LCJH%gridType = ''
    LCJH%allocated = .false.

    deallocate(LCtemp%i, STAT=status)
    deallocate(LCtemp%j, STAT=status)
    deallocate(LCtemp%k, STAT=status)
    deallocate(LCtemp%xyz, STAT=status)
    deallocate(LCtemp%c, STAT=status)
    LCtemp%gridType = ''
    LCtemp%allocated = .false.

  end subroutine E_HinterpSetUp


  !!!!!!!!!!!!!!!E function with variable J!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine E_JinterpSetUp(ii,inGrid,LCE,LCEJ,rho0)

    use modelOperator3D,only: getResTensor

    implicit none
		integer ii
    type (grid_t), target, intent(in) 		:: inGrid
		type (sparsevecc), intent(in) 		    :: LCE
    type (sparsevecc), intent(inout) 		:: LCEJ
    type(modelParam_t), intent(in) :: rho0
    ! local variables
		integer :: nxMax, nyMax, nzMax
		integer :: num,l,m,n,ia
		integer					:: status
		type(rscalar),allocatable  :: RhoE(:)
    character (len=80)	    :: gridType = ''
    
    if(LCEJ%allocated) then
       deallocate(LCEJ%i, STAT=status)
       deallocate(LCEJ%j, STAT=status)
       deallocate(LCEJ%k, STAT=status)
       deallocate(LCEJ%xyz, STAT=status)
       deallocate(LCEJ%c, STAT=status)
       LCEJ%gridType = ''
       LCEJ%allocated = .false.
    endif

    ! maximum number of center nodes
    nxMax = inGrid%nx
    nyMax = inGrid%ny
    nzMax = inGrid%nz

		l = LCE%I(ii)
		m = LCE%J(ii)
		n = LCE%K(ii)

    ! read model parameters form module modelOperator3D_ANI
    allocate(RhoE(6))
		call getResTensor(RhoE) 

    if(LCE%xyz(ii).eq.1) then ! Esx = Fx(Jx,Jy,Jz)

      if(l.eq.1) then ! left

				num = 10
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
        ! 1 - Jx(i,j,k)
				LCEJ%i(1) = l
				LCEJ%j(1) = m
				LCEJ%k(1) = n
				LCEJ%xyz(1) = 1  ! 1 represents x component
				LCEJ%c(1) = RhoE(1)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(5)%v(l,m,n)
        ! 2 - Jx(i+1,j,k)
				LCEJ%i(2) = l + 1
				LCEJ%j(2) = m
				LCEJ%k(2) = n
				LCEJ%xyz(2) = 1
				LCEJ%c(2) = RhoE(1)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(5)%v(l,m,n)  &
				   - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(5)%v(l+1,m,n)
        ! 3 - Jx(i+2,j,k)
 				LCEJ%i(3) = l + 2
				LCEJ%j(3) = m
				LCEJ%k(3) = n
				LCEJ%xyz(3) = 1
				LCEJ%c(3) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(5)%v(l+1,m,n)				   
! ----------------------------------Jy-----------------------------------
        ! 4 - Jy(i,j,k)
 				LCEJ%i(4) = l
				LCEJ%j(4) = m
				LCEJ%k(4) = n
				LCEJ%xyz(4) = 2 ! 2 represents y component
				LCEJ%c(4) = RhoE(4)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(6)%v(l,m,n)		
        ! 5 - Jy(i,j+1,k)
 				LCEJ%i(5) = l
				LCEJ%j(5) = m + 1
				LCEJ%k(5) = n
				LCEJ%xyz(5) = 2
				LCEJ%c(5) = RhoE(4)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(6)%v(l,m,n)		
        ! 6 - Jy(i+1,j,k)
 				LCEJ%i(6) = l + 1
				LCEJ%j(6) = m
				LCEJ%k(6) = n
				LCEJ%xyz(6) = 2
				LCEJ%c(6) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(6)%v(l+1,m,n)			
        ! 7 - Jy(i+1,j+1,k)
 				LCEJ%i(7) = l + 1
				LCEJ%j(7) = m + 1
				LCEJ%k(7) = n
				LCEJ%xyz(7) = 2
				LCEJ%c(7) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(6)%v(l+1,m,n)		
! ----------------------------------Jz-----------------------------------
        ! 8 - Jz(i,j,k)
 				LCEJ%i(8) = l
				LCEJ%j(8) = m
				LCEJ%k(8) = n
				LCEJ%xyz(8) = 3 ! 3 represents z component
				LCEJ%c(8) = RhoE(5)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(3)%v(l,m,n)	
        ! 9 - Jz(i+1,j,k)
 				LCEJ%i(9) = l + 1
				LCEJ%j(9) = m
				LCEJ%k(9) = n
				LCEJ%xyz(9) = 3
				LCEJ%c(9) =  - (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(3)%v(l+1,m,n)
        ! 10 - Jz(i,j,k+1)
 				LCEJ%i(10) = l
				LCEJ%j(10) = m
				LCEJ%k(10) = n + 1
				LCEJ%xyz(10) = 3
				LCEJ%c(10) = RhoE(5)%v(l,m,n)/2.d0
		  								          
	  elseif(l.eq.nxMax) then ! right

				num = 10
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
				! 1 - Jx(i-1,j,k)
				LCEJ%i(1) = l - 1
				LCEJ%j(1) = m
				LCEJ%k(1) = n
				LCEJ%xyz(1) = 1
				LCEJ%c(1) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(5)%v(l-1,m,n)
				! 2 - Jx(i,j,k)
				LCEJ%i(2) = l
				LCEJ%j(2) = m
				LCEJ%k(2) = n
				LCEJ%xyz(2) = 1
				LCEJ%c(2) = RhoE(1)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(5)%v(l-1,m,n)  &
				   + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
				   (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(5)%v(l,m,n)
				! 3 - Jx(i+1,j,k)
				LCEJ%i(3) = l + 1
				LCEJ%j(3) = m
				LCEJ%k(3) = n
				LCEJ%xyz(3) = 1
				LCEJ%c(3) = RhoE(1)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(5)%v(l,m,n)
! ----------------------------------Jy-----------------------------------
				! 4 - Jy(i-1,j,k)
				LCEJ%i(4) = l - 1
				LCEJ%j(4) = m
				LCEJ%k(4) = n
				LCEJ%xyz(4) = 2
				LCEJ%c(4) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(6)%v(l-1,m,n)
				! 5 - Jy(i-1,j+1,k)
				LCEJ%i(5) = l - 1
				LCEJ%j(5) = m + 1
				LCEJ%k(5) = n
				LCEJ%xyz(5) = 2
				LCEJ%c(5) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(6)%v(l-1,m,n)
				! 6 - Jy(i,j,k)
				LCEJ%i(6) = l
				LCEJ%j(6) = m
				LCEJ%k(6) = n
				LCEJ%xyz(6) = 2
				LCEJ%c(6) = RhoE(4)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(6)%v(l,m,n)
				! 7 - Jy(i,j+1,k)
				LCEJ%i(7) = l
				LCEJ%j(7) = m + 1
				LCEJ%k(7) = n
				LCEJ%xyz(7) = 2
				LCEJ%c(7) = RhoE(4)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(6)%v(l,m,n)
! ----------------------------------Jz-----------------------------------
				! 8 - Jz(i-1,j,k)
				LCEJ%i(8) = l - 1
				LCEJ%j(8) = m
				LCEJ%k(8) = n
				LCEJ%xyz(8) = 3
				LCEJ%c(8) = (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(3)%v(l-1,m,n)
				! 9 - Jz(i,j,k)
				LCEJ%i(9) = l
				LCEJ%j(9) = m
				LCEJ%k(9) = n
				LCEJ%xyz(9) = 3
				LCEJ%c(9) = RhoE(5)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(3)%v(l,m,n)
				! 10 - Jz(i,j,k+1)
				LCEJ%i(10) = l
				LCEJ%j(10) = m
				LCEJ%k(10) = n + 1
				LCEJ%xyz(10) = 3
				LCEJ%c(10) = RhoE(5)%v(l,m,n)/2.d0

	  else

				num = 14
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
				! 1 - Jx(i-1,j,k)
				LCEJ%i(1) = l - 1
				LCEJ%j(1) = m
				LCEJ%k(1) = n
				LCEJ%xyz(1) = 1
				LCEJ%c(1) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(5)%v(l-1,m,n)
				! 2 - Jx(i,j,k)
				LCEJ%i(2) = l
				LCEJ%j(2) = m
				LCEJ%k(2) = n
				LCEJ%xyz(2) = 1
				LCEJ%c(2) = RhoE(1)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(5)%v(l-1,m,n)  &
				   + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*(inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1))  &
				   - inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(5)%v(l,m,n)
        ! 3 - Jx(i+1,j,k)
				LCEJ%i(3) = l + 1
				LCEJ%j(3) = m
				LCEJ%k(3) = n
				LCEJ%xyz(3) = 1
				LCEJ%c(3) = RhoE(1)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(5)%v(l+1,m,n)  &
				   + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*(inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1))  &
				   - inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(5)%v(l,m,n)
        ! 4 - Jx(i+2,j,k)
 				LCEJ%i(4) = l + 2
				LCEJ%j(4) = m
				LCEJ%k(4) = n
				LCEJ%xyz(4) = 1
				LCEJ%c(4) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(5)%v(l+1,m,n)
! ----------------------------------Jy-----------------------------------
        ! 5 - Jy(i-1,j,k)
				LCEJ%i(5) = l - 1
				LCEJ%j(5) = m
				LCEJ%k(5) = n
				LCEJ%xyz(5) = 2
				LCEJ%c(5) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(6)%v(l-1,m,n)
				! 6 - Jy(i-1,j+1,k)
				LCEJ%i(6) = l - 1
				LCEJ%j(6) = m + 1
				LCEJ%k(6) = n
				LCEJ%xyz(6) = 2
				LCEJ%c(6) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(6)%v(l-1,m,n)
				! 7 - Jy(i,j,k)
				LCEJ%i(7) = l
				LCEJ%j(7) = m
				LCEJ%k(7) = n
				LCEJ%xyz(7) = 2
				LCEJ%c(7) = RhoE(4)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1))  &
				   - inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(6)%v(l,m,n)
				! 8 - Jy(i,j+1,k)
				LCEJ%i(8) = l
				LCEJ%j(8) = m + 1
				LCEJ%k(8) = n
				LCEJ%xyz(8) = 2
				LCEJ%c(8) = RhoE(4)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1))  &
				   - inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(6)%v(l,m,n)
        ! 9 - Jy(i+1,j,k)
 				LCEJ%i(9) = l + 1
				LCEJ%j(9) = m
				LCEJ%k(9) = n
				LCEJ%xyz(9) = 2
				LCEJ%c(9) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(6)%v(l+1,m,n)			
        ! 10 - Jy(i+1,j+1,k)
 				LCEJ%i(10) = l + 1
				LCEJ%j(10) = m + 1
				LCEJ%k(10) = n
				LCEJ%xyz(10) = 2
				LCEJ%c(10) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(6)%v(l+1,m,n)	
! ----------------------------------Jz-----------------------------------
				! 11 - Jz(i-1,j,k)
				LCEJ%i(11) = l - 1
				LCEJ%j(11) = m
				LCEJ%k(11) = n
				LCEJ%xyz(11) = 3
				LCEJ%c(11) = (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*RhoE(3)%v(l-1,m,n)
				! 12 - Jz(i,j,k)
				LCEJ%i(12) = l
				LCEJ%j(12) = m
				LCEJ%k(12) = n
				LCEJ%xyz(12) = 3
				LCEJ%c(12) = RhoE(5)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1))  &
				   - inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(3)%v(l,m,n)
        ! 13 - Jz(i+1,j,k)
 				LCEJ%i(13) = l + 1
				LCEJ%j(13) = m
				LCEJ%k(13) = n
				LCEJ%xyz(13) = 3
				LCEJ%c(13) =  - (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*RhoE(3)%v(l+1,m,n)
				! 14 - Jz(i,j,k+1)
				LCEJ%i(14) = l
				LCEJ%j(14) = m
				LCEJ%k(14) = n + 1
				LCEJ%xyz(14) = 3
				LCEJ%c(14) = RhoE(5)%v(l,m,n)/2.d0

	  endif


	elseif(LCE%xyz(ii).eq.2) then ! Esy = Fy(Jx,Jy,Jz)


      if(m.eq.1) then ! behind

				num = 10
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
        ! 1 - Jx(i,j,k)
				LCEJ%i(1) = l
				LCEJ%j(1) = m
				LCEJ%k(1) = n
				LCEJ%xyz(1) = 1  ! 1 represents x component
				LCEJ%c(1) = RhoE(4)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(5)%v(l,m,n)
        ! 2 - Jx(i,j+1,k)
				LCEJ%i(2) = l
				LCEJ%j(2) = m + 1
				LCEJ%k(2) = n
				LCEJ%xyz(2) = 1
				LCEJ%c(2) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(5)%v(l,m+1,n)
        ! 3 - Jx(i+1,j,k)
				LCEJ%i(3) = l + 1
				LCEJ%j(3) = m
				LCEJ%k(3) = n
				LCEJ%xyz(3) = 1
				LCEJ%c(3) = RhoE(4)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(5)%v(l,m,n) 
        ! 4 - Jx(i+1,j+1,k)
				LCEJ%i(4) = l + 1
				LCEJ%j(4) = m + 1
				LCEJ%k(4) = n
				LCEJ%xyz(4) = 1
				LCEJ%c(4) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(5)%v(l,m+1,n)
! ----------------------------------Jy-----------------------------------
        ! 5 - Jy(i,j,k)				   
 				LCEJ%i(5) = l
				LCEJ%j(5) = m
				LCEJ%k(5) = n
				LCEJ%xyz(5) = 2 ! 2 represents y component
				LCEJ%c(5) = RhoE(2)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(6)%v(l,m,n)  				   
        ! 6 - Jy(i,j+1,k)				   
 				LCEJ%i(6) = l
				LCEJ%j(6) = m + 1
				LCEJ%k(6) = n
				LCEJ%xyz(6) = 2
				LCEJ%c(6) = RhoE(2)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(6)%v(l,m,n)  &
				   - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(6)%v(l,m+1,n)
        ! 7 - Jy(i,j+2,k)				   
 				LCEJ%i(7) = l
				LCEJ%j(7) = m + 2
				LCEJ%k(7) = n
				LCEJ%xyz(7) = 2
				LCEJ%c(7) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(6)%v(l,m+1,n)						   
! ----------------------------------Jz-----------------------------------
				! 8 - Jz(i,j,k)
				LCEJ%i(8) = l
				LCEJ%j(8) = m
				LCEJ%k(8) = n
				LCEJ%xyz(8) = 3 ! 3 represents z component
				LCEJ%c(8) = RhoE(6)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(3)%v(l,m,n)	
				! 9 - Jz(i,j+1,k)
				LCEJ%i(9) = l
				LCEJ%j(9) = m + 1
				LCEJ%k(9) = n
				LCEJ%xyz(9) = 3
				LCEJ%c(9) = - (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(3)%v(l,m+1,n)
				! 10 - Jz(i,j,k+1)
				LCEJ%i(10) = l
				LCEJ%j(10) = m
				LCEJ%k(10) = n + 1
				LCEJ%xyz(10) = 3
				LCEJ%c(10) = RhoE(6)%v(l,m,n)/2.d0
									           
	  elseif(m.eq.nyMax) then ! front

				num = 10
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
        ! 1 - Jx(i,j-1,k)
				LCEJ%i(1) = l
				LCEJ%j(1) = m - 1
				LCEJ%k(1) = n
				LCEJ%xyz(1) = 1
				LCEJ%c(1) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(5)%v(l,m-1,n)
        ! 2 - Jx(i+1,j-1,k)
				LCEJ%i(2) = l + 1
				LCEJ%j(2) = m - 1
				LCEJ%k(2) = n
				LCEJ%xyz(2) = 1
				LCEJ%c(2) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(5)%v(l,m-1,n)
        ! 3 - Jx(i,j,k)
				LCEJ%i(3) = l
				LCEJ%j(3) = m
				LCEJ%k(3) = n
				LCEJ%xyz(3) = 1
				LCEJ%c(3) = RhoE(4)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(5)%v(l,m,n)
        ! 4 - Jx(i+1,j,k)
				LCEJ%i(4) = l + 1
				LCEJ%j(4) = m
				LCEJ%k(4) = n
				LCEJ%xyz(4) = 1
				LCEJ%c(4) = RhoE(4)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(5)%v(l,m,n)
! ----------------------------------Jy-----------------------------------
        ! 5 - Jy(i,j-1,k)				   
 				LCEJ%i(5) = l
				LCEJ%j(5) = m - 1
				LCEJ%k(5) = n
				LCEJ%xyz(5) = 2
				LCEJ%c(5) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(6)%v(l,m-1,n)
        ! 6 - Jy(i,j,k)				   
 				LCEJ%i(6) = l
				LCEJ%j(6) = m
				LCEJ%k(6) = n
				LCEJ%xyz(6) = 2
				LCEJ%c(6) = RhoE(2)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(6)%v(l,m-1,n)  &
				   + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(6)%v(l,m,n)
        ! 7 - Jy(i,j+1,k)				   
 				LCEJ%i(7) = l
				LCEJ%j(7) = m + 1
				LCEJ%k(7) = n
				LCEJ%xyz(7) = 2
				LCEJ%c(7) = RhoE(2)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(6)%v(l,m,n)
! ----------------------------------Jz-----------------------------------
				! 8 - Jz(i,j-1,k)
				LCEJ%i(8) = l
				LCEJ%j(8) = m - 1
				LCEJ%k(8) = n
				LCEJ%xyz(8) = 3
				LCEJ%c(8) = (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(3)%v(l,m-1,n)
				! 9 - Jz(i,j,k)
				LCEJ%i(9) = l
				LCEJ%j(9) = m
				LCEJ%k(9) = n
				LCEJ%xyz(9) = 3
				LCEJ%c(9) = RhoE(6)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(3)%v(l,m,n)
				! 10 - Jz(i,j,k+1)
				LCEJ%i(10) = l
				LCEJ%j(10) = m
				LCEJ%k(10) = n + 1
				LCEJ%xyz(10) = 3
				LCEJ%c(10) = RhoE(6)%v(l,m,n)/2.d0

	  else

				num = 14
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
        ! 1 - Jx(i,j-1,k)
				LCEJ%i(1) = l
				LCEJ%j(1) = m - 1
				LCEJ%k(1) = n
				LCEJ%xyz(1) = 1
				LCEJ%c(1) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(5)%v(l,m-1,n)
        ! 2 - Jx(i+1,j-1,k)
				LCEJ%i(2) = l + 1
				LCEJ%j(2) = m - 1
				LCEJ%k(2) = n
				LCEJ%xyz(2) = 1
				LCEJ%c(2) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(5)%v(l,m-1,n)
        ! 3 - Jx(i,j,k)
				LCEJ%i(3) = l
				LCEJ%j(3) = m
				LCEJ%k(3) = n
				LCEJ%xyz(3) = 1
				LCEJ%c(3) = RhoE(4)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1))  &
				   - inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(5)%v(l,m,n)
        ! 4 - Jx(i,j+1,k)
				LCEJ%i(4) = l
				LCEJ%j(4) = m + 1
				LCEJ%k(4) = n
				LCEJ%xyz(4) = 1
				LCEJ%c(4) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(5)%v(l,m+1,n)
        ! 5 - Jx(i+1,j,k)
				LCEJ%i(5) = l + 1
				LCEJ%j(5) = m
				LCEJ%k(5) = n
				LCEJ%xyz(5) = 1
				LCEJ%c(5) = RhoE(4)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1))  &
				   - inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(5)%v(l,m,n)
        ! 6 - Jx(i+1,j+1,k)
				LCEJ%i(6) = l + 1
				LCEJ%j(6) = m + 1
				LCEJ%k(6) = n
				LCEJ%xyz(6) = 1
				LCEJ%c(6) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(5)%v(l,m+1,n)
! ----------------------------------Jy-----------------------------------
        ! 7 - Jy(i,j-1,k)				   
 				LCEJ%i(7) = l
				LCEJ%j(7) = m - 1
				LCEJ%k(7) = n
				LCEJ%xyz(7) = 2
				LCEJ%c(7) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(6)%v(l,m-1,n)
        ! 8 - Jy(i,j,k)				   
 				LCEJ%i(8) = l
				LCEJ%j(8) = m
				LCEJ%k(8) = n
				LCEJ%xyz(8) = 2
				LCEJ%c(8) = RhoE(2)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(6)%v(l,m-1,n)  &
				   + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*(inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1))  &
		           - inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(6)%v(l,m,n)
        ! 9 - Jy(i,j+1,k)				   
 				LCEJ%i(9) = l
				LCEJ%j(9) = m + 1
				LCEJ%k(9) = n
				LCEJ%xyz(9) = 2
				LCEJ%c(9) = RhoE(2)%v(l,m,n)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(6)%v(l,m+1,n)  &
				   + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*(inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1))  &
		           - inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(6)%v(l,m,n)
        ! 10 - Jy(i,j+2,k)				   
 				LCEJ%i(10) = l
				LCEJ%j(10) = m + 2
				LCEJ%k(10) = n
				LCEJ%xyz(10) = 2
				LCEJ%c(10) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(6)%v(l,m+1,n)
! ----------------------------------Jz-----------------------------------
				! 11 - Jz(i,j-1,k)
				LCEJ%i(11) = l
				LCEJ%j(11) = m - 1
				LCEJ%k(11) = n
				LCEJ%xyz(11) = 3
				LCEJ%c(11) = (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*RhoE(3)%v(l,m-1,n)
				! 12 - Jz(i,j,k)
				LCEJ%i(12) = l
				LCEJ%j(12) = m
				LCEJ%k(12) = n
				LCEJ%xyz(12) = 3
				LCEJ%c(12) = RhoE(6)%v(l,m,n)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1))  &
				   - inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(3)%v(l,m,n)
				! 13 - Jz(i,j+1,k)
				LCEJ%i(13) = l
				LCEJ%j(13) = m + 1
				LCEJ%k(13) = n
				LCEJ%xyz(13) = 3
				LCEJ%c(13) = - (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*RhoE(3)%v(l,m+1,n)
				! 14 - Jz(i,j,k+1)
				LCEJ%i(14) = l
				LCEJ%j(14) = m
				LCEJ%k(14) = n + 1
				LCEJ%xyz(14) = 3
				LCEJ%c(14) = RhoE(6)%v(l,m,n)/2.d0
						
	  endif

	else

      write(0,*) 'Error: component # out of range in E_JinterpSetUp'

	endif
    
    do ia = 1, 6
      call deall_rscalar(RhoE(ia))
    enddo
    deallocate(RhoE)
    
  end subroutine E_JinterpSetUp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine SensE_JinterpSetUp(ii,inGrid,LCE,LCEJ,iPar,mx)
    implicit none
	  integer ii
    type (grid_t), target, intent(in) 		:: inGrid
	  type (sparsevecc), intent(in) 		    :: LCE
    type (sparsevecc), intent(inout) 		:: LCEJ
    integer, intent(in) :: iPar,mx(3)
    ! local variables
	  integer :: nxMax, nyMax, nzMax
	  integer :: num,l,m,n,km,ia
	  integer	 :: status
    character (len=80)	:: gridType = ''
    real(kind=prec) :: mpar
    
    if(LCEJ%allocated) then
       deallocate(LCEJ%i, STAT=status)
       deallocate(LCEJ%j, STAT=status)
       deallocate(LCEJ%k, STAT=status)
       deallocate(LCEJ%xyz, STAT=status)
       deallocate(LCEJ%c, STAT=status)
       LCEJ%gridType = ''
       LCEJ%allocated = .false.
    endif

    ! maximum number of center nodes
    nxMax = inGrid%nx
    nyMax = inGrid%ny
    nzMax = inGrid%nz

	  l = LCE%I(ii)
	  m = LCE%J(ii)
	  n = LCE%K(ii)
	  km = n - inGrid%nzAir

	  if(iPar.le.3) then
	    mpar = aRes(iPar)%v(mx(1),mx(2),mx(3))
	  else
	    mpar = D2R
	  endif
	   
	!print*,'ix-iy-iz',l,m,n  !debug

    if(LCE%xyz(ii).eq.1) then ! Esx = Fx(Jx,Jy,Jz)

      if(l.eq.1) then ! left

	    num = 10
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
        ! 1 - Jx(i,j,k)
		LCEJ%i(1) = l
		LCEJ%j(1) = m
		LCEJ%k(1) = n
		LCEJ%xyz(1) = 1  ! 1 represents x component
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then 
		  LCEJ%c(1) = (dxx(ipar)%v(l,m,km)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(1) = C_ZERO
		endif
        ! 2 - Jx(i+1,j,k)
		LCEJ%i(2) = l + 1
		LCEJ%j(2) = m
		LCEJ%k(2) = n
		LCEJ%xyz(2) = 1
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(2) = (dxx(ipar)%v(l,m,km)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*dxz(ipar)%v(l,m,km))*mpar
		elseif((l+1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(2) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*dxz(ipar)%v(l+1,m,km)*mpar
		else
		  LCEJ%c(2) = C_ZERO
		endif
        ! 3 - Jx(i+2,j,k)
 		LCEJ%i(3) = l + 2
		LCEJ%j(3) = m
		LCEJ%k(3) = n
		LCEJ%xyz(3) = 1
		if((l+1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(3) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*dxz(ipar)%v(l+1,m,km)*mpar
		else
		  LCEJ%c(3) = C_ZERO 
		endif				   
! ----------------------------------Jy-----------------------------------
        ! 4 - Jy(i,j,k)
 		LCEJ%i(4) = l
		LCEJ%j(4) = m
		LCEJ%k(4) = n
		LCEJ%xyz(4) = 2 ! 2 represents y component
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(4) = (dxy(ipar)%v(l,m,km)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(4) = C_ZERO
		endif
        ! 5 - Jy(i,j+1,k)
 		LCEJ%i(5) = l
		LCEJ%j(5) = m + 1
		LCEJ%k(5) = n
		LCEJ%xyz(5) = 2
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(5) = (dxy(ipar)%v(l,m,km)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(5) = C_ZERO
		endif		
        ! 6 - Jy(i+1,j,k)
 		LCEJ%i(6) = l + 1
		LCEJ%j(6) = m
		LCEJ%k(6) = n
		LCEJ%xyz(6) = 2
		if((l+1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(6) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*dyz(ipar)%v(l+1,m,km)*mpar	
		else
		  LCEJ%c(6) = C_ZERO
		endif			
        ! 7 - Jy(i+1,j+1,k)
 		LCEJ%i(7) = l + 1
		LCEJ%j(7) = m + 1
		LCEJ%k(7) = n
		LCEJ%xyz(7) = 2
		if((l+1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(7) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*dyz(ipar)%v(l+1,m,km)*mpar
		else
		  LCEJ%c(7) = C_ZERO
		endif				
! ----------------------------------Jz-----------------------------------
        ! 8 - Jz(i,j,k)
 		LCEJ%i(8) = l
		LCEJ%j(8) = m
		LCEJ%k(8) = n
		LCEJ%xyz(8) = 3 ! 3 represents z component
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(8) = (dxz(ipar)%v(l,m,km)/2.d0 - (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*dzz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(8) = C_ZERO
		endif		
        ! 9 - Jz(i+1,j,k)
 		LCEJ%i(9) = l + 1
		LCEJ%j(9) = m
		LCEJ%k(9) = n
		LCEJ%xyz(9) = 3
		if((l+1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(9) =  - (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*dzz(ipar)%v(l+1,m,km)*mpar
		else
		  LCEJ%c(9) = C_ZERO
		endif		
        ! 10 - Jz(i,j,k+1)
 		LCEJ%i(10) = l
		LCEJ%j(10) = m
		LCEJ%k(10) = n + 1
		LCEJ%xyz(10) = 3
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(10) = dxz(ipar)%v(l,m,km)/2.d0*mpar
		else
		  LCEJ%c(10) = C_ZERO
		endif
		  								          
	  elseif(l.eq.nxMax) then ! right

	    num = 10
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
		! 1 - Jx(i-1,j,k)
		LCEJ%i(1) = l - 1
		LCEJ%j(1) = m
		LCEJ%k(1) = n
		LCEJ%xyz(1) = 1
		if((l-1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(1) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*dxz(ipar)%v(l-1,m,km)*mpar
		else
		  LCEJ%c(1) = C_ZERO
		endif		
		! 2 - Jx(i,j,k)
		LCEJ%i(2) = l
		LCEJ%j(2) = m
		LCEJ%k(2) = n
		LCEJ%xyz(2) = 1
		if((l-1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(2) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*dxz(ipar)%v(l-1,m,km)*mpar
		elseif(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(2) = (dxx(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
				   (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(2) = C_ZERO
		endif		
		! 3 - Jx(i+1,j,k)
		LCEJ%i(3) = l + 1
		LCEJ%j(3) = m
		LCEJ%k(3) = n
		LCEJ%xyz(3) = 1
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(3) = (dxx(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(3) = C_ZERO
		endif				
! ----------------------------------Jy-----------------------------------
		! 4 - Jy(i-1,j,k)
		LCEJ%i(4) = l - 1
		LCEJ%j(4) = m
		LCEJ%k(4) = n
		LCEJ%xyz(4) = 2
		if((l-1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(4) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*dyz(ipar)%v(l-1,m,km)*mpar
		else
		  LCEJ%c(4) = C_ZERO
		endif		
		! 5 - Jy(i-1,j+1,k)
		LCEJ%i(5) = l - 1
		LCEJ%j(5) = m + 1
		LCEJ%k(5) = n
		LCEJ%xyz(5) = 2
		if((l-1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(5) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*dyz(ipar)%v(l-1,m,km)*mpar
		else
		  LCEJ%c(5) = C_ZERO
		endif		
		! 6 - Jy(i,j,k)
		LCEJ%i(6) = l
		LCEJ%j(6) = m
		LCEJ%k(6) = n
		LCEJ%xyz(6) = 2
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(6) = (dxy(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(6) = C_ZERO
		endif		
		! 7 - Jy(i,j+1,k)
		LCEJ%i(7) = l
		LCEJ%j(7) = m + 1
		LCEJ%k(7) = n
		LCEJ%xyz(7) = 2
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(7) = (dxy(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(7) = C_ZERO
		endif		        
! ----------------------------------Jz-----------------------------------
		! 8 - Jz(i-1,j,k)
		LCEJ%i(8) = l - 1
		LCEJ%j(8) = m
		LCEJ%k(8) = n
		LCEJ%xyz(8) = 3
		if((l-1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(8) = (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*dzz(ipar)%v(l-1,m,km)*mpar
		else
		  LCEJ%c(8) = C_ZERO
		endif		    
		! 9 - Jz(i,j,k)
		LCEJ%i(9) = l
		LCEJ%j(9) = m
		LCEJ%k(9) = n
		LCEJ%xyz(9) = 3
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(9) = (dxz(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1)))*dzz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(9) = C_ZERO
		endif		
		! 10 - Jz(i,j,k+1)
		LCEJ%i(10) = l
		LCEJ%j(10) = m
		LCEJ%k(10) = n + 1
		LCEJ%xyz(10) = 3
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(10) = dxz(ipar)%v(l,m,km)/2.d0*mpar
		else
		  LCEJ%c(10) = C_ZERO
		endif		

	  else

	    num = 14
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
		! 1 - Jx(i-1,j,k)
		LCEJ%i(1) = l - 1
		LCEJ%j(1) = m
		LCEJ%k(1) = n
		LCEJ%xyz(1) = 1
		if((l-1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(1) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*dxz(ipar)%v(l-1,m,km)*mpar
		else
		  LCEJ%c(1) = C_ZERO
		endif		
		! 2 - Jx(i,j,k)
		LCEJ%i(2) = l
		LCEJ%j(2) = m
		LCEJ%k(2) = n
		LCEJ%xyz(2) = 1
		if((l-1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(2) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*dxz(ipar)%v(l-1,m,km)*mpar
		elseif(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(2) = (dxx(ipar)%v(l,m,km)/2.d0  &
				   + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*(inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1))  &
				   - inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(2) = C_ZERO
		endif		
        ! 3 - Jx(i+1,j,k)
		LCEJ%i(3) = l + 1
		LCEJ%j(3) = m
		LCEJ%k(3) = n
		LCEJ%xyz(3) = 1
		if((l+1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(3) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*dxz(ipar)%v(l+1,m,km)*mpar
		elseif(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(3) = (dxx(ipar)%v(l,m,km)/2.d0  &
				   + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*(inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1))  &
				   - inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(3) = C_ZERO
		endif		
        ! 4 - Jx(i+2,j,k)
 		LCEJ%i(4) = l + 2
		LCEJ%j(4) = m
		LCEJ%k(4) = n
		LCEJ%xyz(4) = 1
		if((l+1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(4) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*dxz(ipar)%v(l+1,m,km)*mpar
		else
		  LCEJ%c(4) = C_ZERO
		endif		
! ----------------------------------Jy-----------------------------------
        ! 5 - Jy(i-1,j,k)
		LCEJ%i(5) = l - 1
		LCEJ%j(5) = m
		LCEJ%k(5) = n
		LCEJ%xyz(5) = 2
		if((l-1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(5) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*dyz(ipar)%v(l-1,m,km)*mpar
		else
		  LCEJ%c(5) = C_ZERO
		endif		
		! 6 - Jy(i-1,j+1,k)
		LCEJ%i(6) = l - 1
		LCEJ%j(6) = m + 1
		LCEJ%k(6) = n
		LCEJ%xyz(6) = 2
		if((l-1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(6) = (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*dyz(ipar)%v(l-1,m,km)*mpar
		else
		  LCEJ%c(6) = C_ZERO
		endif		
		! 7 - Jy(i,j,k)
		LCEJ%i(7) = l
		LCEJ%j(7) = m
		LCEJ%k(7) = n
		LCEJ%xyz(7) = 2
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(7) = (dxy(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1))  &
				   - inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(7) = C_ZERO
		endif		
		! 8 - Jy(i,j+1,k)
		LCEJ%i(8) = l
		LCEJ%j(8) = m + 1
		LCEJ%k(8) = n
		LCEJ%xyz(8) = 2
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(8) = (dxy(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1))  &
				   - inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(8) = C_ZERO
		endif		
        ! 9 - Jy(i+1,j,k)
 		LCEJ%i(9) = l + 1
		LCEJ%j(9) = m
		LCEJ%k(9) = n
		LCEJ%xyz(9) = 2
		if((l+1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(9) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*dyz(ipar)%v(l+1,m,km)*mpar
		else
		  LCEJ%c(9) = C_ZERO
		endif					
        ! 10 - Jy(i+1,j+1,k)
 		LCEJ%i(10) = l + 1
		LCEJ%j(10) = m + 1
		LCEJ%k(10) = n
		LCEJ%xyz(10) = 2
		if((l+1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(10) = - (inGrid%dz(n)/inGrid%dx(l)/4.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*dyz(ipar)%v(l+1,m,km)*mpar
		else
		  LCEJ%c(10) = C_ZERO
		endif			
! ----------------------------------Jz-----------------------------------
		! 11 - Jz(i-1,j,k)
		LCEJ%i(11) = l - 1
		LCEJ%j(11) = m
		LCEJ%k(11) = n
		LCEJ%xyz(11) = 3
		if((l-1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(11) = (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l-1)))*dzz(ipar)%v(l-1,m,km)*mpar
		else
		  LCEJ%c(11) = C_ZERO
		endif		
		! 12 - Jz(i,j,k)
		LCEJ%i(12) = l
		LCEJ%j(12) = m
		LCEJ%k(12) = n
		LCEJ%xyz(12) = 3
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(12) = (dxz(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l-1)/(inGrid%dx(l)+inGrid%dx(l-1))  &
				   - inGrid%dx(l+1)/(inGrid%dx(l)+inGrid%dx(l+1)))*dzz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(12) = C_ZERO
		endif		
        ! 13 - Jz(i+1,j,k)
 		LCEJ%i(13) = l + 1
		LCEJ%j(13) = m
		LCEJ%k(13) = n
		LCEJ%xyz(13) = 3
		if((l+1).eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(13) =  - (inGrid%dz(n)/inGrid%dx(l)/2.d0)*  &
		           (inGrid%dx(l)/(inGrid%dx(l)+inGrid%dx(l+1)))*dzz(ipar)%v(l+1,m,km)*mpar
		else
		  LCEJ%c(13) = C_ZERO
		endif		
		! 14 - Jz(i,j,k+1)
		LCEJ%i(14) = l
		LCEJ%j(14) = m
		LCEJ%k(14) = n + 1
		LCEJ%xyz(14) = 3
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(14) = dxz(ipar)%v(l,m,km)/2.d0*mpar
		else
		  LCEJ%c(14) = C_ZERO
		endif		

	  endif


	elseif(LCE%xyz(ii).eq.2) then ! Esy = Fy(Jx,Jy,Jz)


      if(m.eq.1) then ! behind

	    num = 10
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
        ! 1 - Jx(i,j,k)
		LCEJ%i(1) = l
		LCEJ%j(1) = m
		LCEJ%k(1) = n
		LCEJ%xyz(1) = 1  ! 1 represents x component
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(1) = (dxy(ipar)%v(l,m,km)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(1) = C_ZERO
		endif		
        ! 2 - Jx(i,j+1,k)
		LCEJ%i(2) = l
		LCEJ%j(2) = m + 1
		LCEJ%k(2) = n
		LCEJ%xyz(2) = 1
		if(l.eq.mx(1).and.(m+1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(2) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*dxz(ipar)%v(l,m+1,km)*mpar
		else
		  LCEJ%c(2) = C_ZERO
		endif		
        ! 3 - Jx(i+1,j,k)
		LCEJ%i(3) = l + 1
		LCEJ%j(3) = m
		LCEJ%k(3) = n
		LCEJ%xyz(3) = 1
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(3) = (dxy(ipar)%v(l,m,km)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(3) = C_ZERO
		endif		
        ! 4 - Jx(i+1,j+1,k)
		LCEJ%i(4) = l + 1
		LCEJ%j(4) = m + 1
		LCEJ%k(4) = n
		LCEJ%xyz(4) = 1
		if(l.eq.mx(1).and.(m+1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(4) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*dxz(ipar)%v(l,m+1,km)*mpar
		else
		  LCEJ%c(4) = C_ZERO
		endif		
! ----------------------------------Jy-----------------------------------
        ! 5 - Jy(i,j,k)				   
 		LCEJ%i(5) = l
		LCEJ%j(5) = m
		LCEJ%k(5) = n
		LCEJ%xyz(5) = 2 ! 2 represents y component
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(5) = (dyy(ipar)%v(l,m,km)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(5) = C_ZERO
		endif		
        ! 6 - Jy(i,j+1,k)				   
 		LCEJ%i(6) = l
		LCEJ%j(6) = m + 1
		LCEJ%k(6) = n
		LCEJ%xyz(6) = 2
		if(l.eq.mx(1).and.(m+1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(6) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*dyz(ipar)%v(l,m+1,km)*mpar
		elseif(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(6) = (dyy(ipar)%v(l,m,km)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(6) = C_ZERO
		endif		
        ! 7 - Jy(i,j+2,k)				   
 		LCEJ%i(7) = l
		LCEJ%j(7) = m + 2
		LCEJ%k(7) = n
		LCEJ%xyz(7) = 2
		if(l.eq.mx(1).and.(m+1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(7) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*dyz(ipar)%v(l,m+1,km)*mpar
		else
		  LCEJ%c(7) = C_ZERO
		endif								   
! ----------------------------------Jz-----------------------------------
		! 8 - Jz(i,j,k)
		LCEJ%i(8) = l
		LCEJ%j(8) = m
		LCEJ%k(8) = n
		LCEJ%xyz(8) = 3 ! 3 represents z component
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(8) = (dyz(ipar)%v(l,m,km)/2.d0 - (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*dzz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(8) = C_ZERO
		endif		
		! 9 - Jz(i,j+1,k)
		LCEJ%i(9) = l
		LCEJ%j(9) = m + 1
		LCEJ%k(9) = n
		LCEJ%xyz(9) = 3
		if(l.eq.mx(1).and.(m+1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(9) = - (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*dzz(ipar)%v(l,m+1,km)*mpar
		else
		  LCEJ%c(9) = C_ZERO
		endif		
		! 10 - Jz(i,j,k+1)
		LCEJ%i(10) = l
		LCEJ%j(10) = m
		LCEJ%k(10) = n + 1
		LCEJ%xyz(10) = 3
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(10) = dyz(ipar)%v(l,m,km)/2.d0*mpar
		else
		  LCEJ%c(10) = C_ZERO
		endif		
									           
	  elseif(m.eq.nyMax) then ! front

	    num = 10
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
        ! 1 - Jx(i,j-1,k)
		LCEJ%i(1) = l
		LCEJ%j(1) = m - 1
		LCEJ%k(1) = n
		LCEJ%xyz(1) = 1
		if(l.eq.mx(1).and.(m-1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(1) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*dxz(ipar)%v(l,m-1,km)*mpar
		else
		  LCEJ%c(1) = C_ZERO
		endif		
        ! 2 - Jx(i+1,j-1,k)
		LCEJ%i(2) = l + 1
		LCEJ%j(2) = m - 1
		LCEJ%k(2) = n
		LCEJ%xyz(2) = 1
		if(l.eq.mx(1).and.(m-1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(2) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*dxz(ipar)%v(l,m-1,km)*mpar
		else
		  LCEJ%c(2) = C_ZERO
		endif		
        ! 3 - Jx(i,j,k)
		LCEJ%i(3) = l
		LCEJ%j(3) = m
		LCEJ%k(3) = n
		LCEJ%xyz(3) = 1
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(3) = (dxy(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(3) = C_ZERO
		endif		
        ! 4 - Jx(i+1,j,k)
		LCEJ%i(4) = l + 1
		LCEJ%j(4) = m
		LCEJ%k(4) = n
		LCEJ%xyz(4) = 1
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(4) = (dxy(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(4) = C_ZERO
		endif		
! ----------------------------------Jy-----------------------------------
        ! 5 - Jy(i,j-1,k)				   
 		LCEJ%i(5) = l
		LCEJ%j(5) = m - 1
		LCEJ%k(5) = n
		LCEJ%xyz(5) = 2
		if(l.eq.mx(1).and.(m-1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(5) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*dyz(ipar)%v(l,m-1,km)*mpar
		else
		  LCEJ%c(5) = C_ZERO
		endif		
        ! 6 - Jy(i,j,k)				   
 		LCEJ%i(6) = l
		LCEJ%j(6) = m
		LCEJ%k(6) = n
		LCEJ%xyz(6) = 2
		if(l.eq.mx(1).and.(m-1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(6) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*dyz(ipar)%v(l,m-1,km)*mpar
		elseif(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(6) = (dyy(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(6) = C_ZERO
		endif		
        ! 7 - Jy(i,j+1,k)				   
 		LCEJ%i(7) = l
		LCEJ%j(7) = m + 1
		LCEJ%k(7) = n
		LCEJ%xyz(7) = 2
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(7) = (dyy(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(7) = C_ZERO
		endif		
! ----------------------------------Jz-----------------------------------
		! 8 - Jz(i,j-1,k)
		LCEJ%i(8) = l
		LCEJ%j(8) = m - 1
		LCEJ%k(8) = n
		LCEJ%xyz(8) = 3
		if(l.eq.mx(1).and.(m-1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(8) = (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*dzz(ipar)%v(l,m-1,km)*mpar
		else
		  LCEJ%c(8) = C_ZERO
		endif		
		! 9 - Jz(i,j,k)
		LCEJ%i(9) = l
		LCEJ%j(9) = m
		LCEJ%k(9) = n
		LCEJ%xyz(9) = 3
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(9) = (dyz(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1)))*dzz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(9) = C_ZERO
		endif		
		! 10 - Jz(i,j,k+1)
		LCEJ%i(10) = l
		LCEJ%j(10) = m
		LCEJ%k(10) = n + 1
		LCEJ%xyz(10) = 3
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(10) = dyz(ipar)%v(l,m,km)/2.d0*mpar
		else
		  LCEJ%c(10) = C_ZERO
		endif		

	  else

	    num = 14
        gridType = EDGE
        Call create_sparsevecc(num,LCEJ,gridType)
! ----------------------------------Jx-----------------------------------
        ! 1 - Jx(i,j-1,k)
		LCEJ%i(1) = l
		LCEJ%j(1) = m - 1
		LCEJ%k(1) = n
		LCEJ%xyz(1) = 1
		if(l.eq.mx(1).and.(m-1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(1) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*dxz(ipar)%v(l,m-1,km)*mpar
		else
		  LCEJ%c(1) = C_ZERO
		endif		
        ! 2 - Jx(i+1,j-1,k)
		LCEJ%i(2) = l + 1
		LCEJ%j(2) = m - 1
		LCEJ%k(2) = n
		LCEJ%xyz(2) = 1
		if(l.eq.mx(1).and.(m-1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(2) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*dxz(ipar)%v(l,m-1,km)*mpar
		else
		  LCEJ%c(2) = C_ZERO
		endif		
        ! 3 - Jx(i,j,k)
		LCEJ%i(3) = l
		LCEJ%j(3) = m
		LCEJ%k(3) = n
		LCEJ%xyz(3) = 1
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(3) = (dxy(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1))  &
				   - inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(3) = C_ZERO
		endif		
        ! 4 - Jx(i,j+1,k)
		LCEJ%i(4) = l
		LCEJ%j(4) = m + 1
		LCEJ%k(4) = n
		LCEJ%xyz(4) = 1
		if(l.eq.mx(1).and.(m+1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(4) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*dxz(ipar)%v(l,m+1,km)*mpar
		else
		  LCEJ%c(4) = C_ZERO
		endif		
        ! 5 - Jx(i+1,j,k)
		LCEJ%i(5) = l + 1
		LCEJ%j(5) = m
		LCEJ%k(5) = n
		LCEJ%xyz(5) = 1
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(5) = (dxy(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1))  &
				   - inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*dxz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(5) = C_ZERO
		endif		
        ! 6 - Jx(i+1,j+1,k)
		LCEJ%i(6) = l + 1
		LCEJ%j(6) = m + 1
		LCEJ%k(6) = n
		LCEJ%xyz(6) = 1
		if(l.eq.mx(1).and.(m+1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(6) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*dxz(ipar)%v(l,m+1,km)*mpar
		else
		  LCEJ%c(6) = C_ZERO
		endif		
! ----------------------------------Jy-----------------------------------
        ! 7 - Jy(i,j-1,k)				   
 		LCEJ%i(7) = l
		LCEJ%j(7) = m - 1
		LCEJ%k(7) = n
		LCEJ%xyz(7) = 2
		if(l.eq.mx(1).and.(m-1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(7) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*dyz(ipar)%v(l,m-1,km)*mpar
		else
		  LCEJ%c(7) = C_ZERO
		endif		
        ! 8 - Jy(i,j,k)				   
 		LCEJ%i(8) = l
		LCEJ%j(8) = m
		LCEJ%k(8) = n
		LCEJ%xyz(8) = 2
		if(l.eq.mx(1).and.(m-1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(8) = (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*dyz(ipar)%v(l,m-1,km)*mpar
		elseif(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(8) = (dyy(ipar)%v(l,m,km)/2.d0  &
				   + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*(inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1))  &
		           - inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(8) = C_ZERO
		endif		
        ! 9 - Jy(i,j+1,k)				   
 		LCEJ%i(9) = l
		LCEJ%j(9) = m + 1
		LCEJ%k(9) = n
		LCEJ%xyz(9) = 2
		if(l.eq.mx(1).and.(m+1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(9) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*dyz(ipar)%v(l,m+1,km)*mpar
		elseif(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(9) = (dyy(ipar)%v(l,m,km)/2.d0  &
				   + (inGrid%dz(n)/inGrid%dy(m)/4.d0)*(inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1))  &
		           - inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*dyz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(9) = C_ZERO
		endif		
        ! 10 - Jy(i,j+2,k)				   
 		LCEJ%i(10) = l
		LCEJ%j(10) = m + 2
		LCEJ%k(10) = n
		LCEJ%xyz(10) = 2
		if(l.eq.mx(1).and.(m+1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(10) = - (inGrid%dz(n)/inGrid%dy(m)/4.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*dyz(ipar)%v(l,m+1,km)*mpar
		else
		  LCEJ%c(10) = C_ZERO
		endif		
! ----------------------------------Jz-----------------------------------
		! 11 - Jz(i,j-1,k)
		LCEJ%i(11) = l
		LCEJ%j(11) = m - 1
		LCEJ%k(11) = n
		LCEJ%xyz(11) = 3
		if(l.eq.mx(1).and.(m-1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(11) = (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m-1)))*dzz(ipar)%v(l,m-1,km)*mpar
		else
		  LCEJ%c(11) = C_ZERO
		endif		
		! 12 - Jz(i,j,k)
		LCEJ%i(12) = l
		LCEJ%j(12) = m
		LCEJ%k(12) = n
		LCEJ%xyz(12) = 3
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(12) = (dyz(ipar)%v(l,m,km)/2.d0 + (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m-1)/(inGrid%dy(m)+inGrid%dy(m-1))  &
				   - inGrid%dy(m+1)/(inGrid%dy(m)+inGrid%dy(m+1)))*dzz(ipar)%v(l,m,km))*mpar
		else
		  LCEJ%c(12) = C_ZERO
		endif		
		! 13 - Jz(i,j+1,k)
		LCEJ%i(13) = l
		LCEJ%j(13) = m + 1
		LCEJ%k(13) = n
		LCEJ%xyz(13) = 3
		if(l.eq.mx(1).and.(m+1).eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(13) = - (inGrid%dz(n)/inGrid%dy(m)/2.d0)*  &
		           (inGrid%dy(m)/(inGrid%dy(m)+inGrid%dy(m+1)))*dzz(ipar)%v(l,m+1,km)*mpar
		else
		  LCEJ%c(13) = C_ZERO
		endif		
		! 14 - Jz(i,j,k+1)
		LCEJ%i(14) = l
		LCEJ%j(14) = m
		LCEJ%k(14) = n + 1
		LCEJ%xyz(14) = 3
		if(l.eq.mx(1).and.m.eq.mx(2).and.km.eq.mx(3)) then
		  LCEJ%c(14) = dyz(ipar)%v(l,m,km)/2.d0*mpar
		else
		  LCEJ%c(14) = C_ZERO
		endif		
						
	  endif

	else

      write(0,*) 'Error: component # out of range in E_JinterpSetUp'

	endif
          
  end subroutine SensE_JinterpSetUp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine create_SensE(sigma0)
      implicit none
      type(modelParam_t), intent(in) :: sigma0
      character(80)  :: paramType = 'LINEAR'  !!!!!
      
      allocate(aRes(3),dxx(6),dxy(6),dxz(6),dyy(6),dyz(6),dzz(6))
      call getResValue_modelParam(sigma0,paramType,aRes)
      call SensModelParamToTensor(sigma0,dxx,dxy,dxz,dyy,dyz,dzz)
      
  end subroutine create_SensE
   
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine deall_SensE
    implicit none
    integer ia
     
    ! release memory
    do ia = 1,3
       call deall_rscalar(aRes(ia))
    enddo        
    do ia=1,6
       call deall_rscalar(dxx(ia))
       call deall_rscalar(dxy(ia))
       call deall_rscalar(dxz(ia))
       call deall_rscalar(dyy(ia))
       call deall_rscalar(dyz(ia))
       call deall_rscalar(dzz(ia))
    end do
    
    deallocate(aRes,dxx,dxy,dxz,dyy,dyz,dzz)
    
  end subroutine deall_SensE

end module EMfieldInterp
