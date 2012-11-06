! *****************************************************************************
module EMfieldInterp
  ! Generic data functionals for interpolation of electric and magnetic
  ! fields to an arbitrary point within the model domain for 3D staggered
  ! grid finite difference solutions

  use utilities
  use sg_sparse_vector
  use ModelSpace, sigC => ModelParamToOneEdge

  implicit none

  public				:: EinterpSetUp
  public				:: BinterpSetUp, BfromESetUp

Contains

  ! **************************************************************************
  ! electrical field coeffcients in the sparse vector
  subroutine EinterpSetUp(inGrid,x,xyz,LC,CondE)

  ! _________________________________________________________________________________
    ! sets up coefficients in sparse vector LC for evaluation/interpolation
    ! at location x of electric field vector defined on edges of staggered
    ! grid. xyz (1 = x; 2 = y; 3 = z) gives the component
    ! INTERPOLATION METHOD: tri-linear spline
    ! Allocates arrays in sparse vector LC, deallocating first if already
    ! allocated   ...  note that in general the number of non-zero
    ! coefficients in LC is in general computed within this routine
    ! at execution time
    !   Oct 11 2004:  GE : optional edge conductivity added to allow interpolation
    !    of currents instead of electric field in direction of component
  !___________________________________________________________________________________

    implicit none
    type(grid_t), target, intent(in)              :: inGrid
    real(kind=prec), dimension(3), intent(in)     :: x
    integer, intent(in)                           :: xyz
    type(sparsevecc), intent(inout)               :: LC ! this is mg sparsevector
    type(modelParam_t), intent(in), optional      :: CondE

    ! Local Variables
    integer                                     :: i0,j0,k0,ii,n,m,p,ic,jc,kc
    integer, dimension(8)                       :: I,J,K
    integer                                     :: nxMax, nyMax, nzMax
    integer                                     :: izz,izc, imgrid,nzCum, k0SG,currSG
    integer                                     :: status
    complex(kind=prec), dimension(8)    :: C
    real(kind=prec), dimension(3,2)     :: w
    real(kind=prec)                     :: wadd
    real(kind=prec)                     :: x3SG
    character (len = 80)                        :: gridType = ''
    logical                     		:: UseCond, found

    integer, parameter		:: IX = 1, IY = 2, IZ = 3

    UseCond = present(CondE)
    if(LC%allocated) then
      call  deall_sparsevecc(LC)
    endif

    ! which sub-grid x(3) belong to?  ! x3SD
    ! derive the number of this sub-grid  ! currSG
    nzCum = 0
    Subgrids: do imgrid = 1, inGrid%mgridSize
      Zlayers: do izz = 1, inGrid%gridArray(imgrid)%nz
        izc = izz+nzCum
        if(clean(inGrid%zEdge(izc)) .gt. clean(x(3)) ) then
          currSG = imgrid
          if(izz == 1) then
            x3SG = inGrid%gridArray(imgrid-1)%zEdge(inGrid%gridArray(imgrid-1)%nz+1)
          else
            x3SG = inGrid%gridArray(imgrid)%zEdge(izz-1)
          endif
          found = .true.
        exit
        endif
      enddo Zlayers
      if(found == .true.) exit
      nzCum = nzCum + inGrid%gridArray(imgrid)%nz
    enddo Subgrids

    ! Do interpolation within current sub-grid, first
    ! if UseCond = .false. then do not need to worry about the interface layer
    ! maximum number of edge nodes
    nxMax = inGrid%gridArray(currSG)%nx+1
    nyMax = inGrid%gridArray(currSG)%ny+1
    nzMax = inGrid%nz+1
    i0 = minNode(x(1),inGrid%gridArray(currSG)%xEdge)
    j0 = minNode(x(2),inGrid%gridArray(currSG)%yEdge)
    k0 = minNode(x(3),inGrid%zEdge)
    ! we need to store k0 in sparse vector with respect to current sub-grid
    k0SG = minNode(x3SG,inGrid%gridArray(currSG)%zEdge)
    if (xyz .eq. 1) then
       ! Evaluation of x component wrt electrical vectors on cubic edges
       ic = i0
       i0 = minNode(x(1),inGrid%gridArray(currSG)%xCenter)
       ! maximum number of center nodes
       nxMax = nxMax -1
    elseif (xyz .eq. 2) then
       ! Evaluation of y component wrt electrical vectors on cubic edges
       jc = j0
       j0 = minNode(x(2),inGrid%gridArray(currSG)%yCenter)
       ! maximum number of center nodes
       nyMax = nyMax -1
    elseif (xyz .eq. 3) then
       ! Evaluation of z component wrt electrical vectors on cubic edges
       kc = k0
       k0 = minNode(x(3),inGrid%zCenter)
       ! we need to store k0 in sparse vector with respect to current sub-grid
       k0SG = minNode(x3SG,inGrid%gridArray(currSG)%zCenter)
       ! maximum number of center nodes
       nzMax = nzMax -1
    else
       write(0,*) 'Error: component # out of range in EinterpSetUp'
    endif

    if((i0.gt.0).and.(i0.lt.nxMax)) then
       if(xyz.eq.1) then
          w(1,2) = (x(1) - inGrid%gridArray(currSG)%xCenter(i0))/(inGrid%gridArray(currSG)%delX(i0+1))
       else
          w(1,2) = (x(1) - inGrid%gridArray(currSG)%xEdge(i0))/(inGrid%gridArray(currSG)%dx(i0))
       endif
    elseif(i0.le.0) then
       w(1,2) = 1
       UseCOnd = .false.
    else
       w(1,2) = 0
       UseCond = .false.
    endif

    if((j0.gt.0).and.(j0.lt.nyMax)) then
       if(xyz.eq.2) then
          w(2,2) = (x(2) - inGrid%gridArray(currSG)%yCenter(j0))/(inGrid%gridArray(currSG)%delY(j0+1))
       else
          w(2,2) = (x(2) - inGrid%gridArray(currSG)%yEdge(j0))/(inGrid%gridArray(currSG)%dy(j0))
       endif
    elseif(j0.le.0) then
       w(2,2) = 1
       UseCond = .false.
    else
       w(2,2) = 0
       UseCond = .false.
    endif

    if((k0.gt.0).and.(k0.lt.nzMax)) then
       if(xyz.eq.3) then
          w(3,2) = (x(3) - inGrid%zCenter(k0))/(inGrid%delZ(k0+1))
       else
          w(3,2) = (x(3) - inGrid%zEdge(k0))/(inGrid%dz(k0))
       endif
    elseif(k0.le.0) then
       w(3,2) = 1
       UseCond = .false.
    else
       w(3,2) = 0
       UseCond = .false.
    endif

    w(1,1) = 1-w(1,2)
    w(2,1) = 1-w(2,2)
    w(3,1) = 1-w(3,2)

    ! UseCond then have to think about interface layers (m to EDGE)
    if(UseCond) then
        !  modify weights to allow for discontinuities in conductivity
        !  in the direction of the interpolated E component
        if(xyz.eq.1) then
           if(ic.eq.i0) then
           !    interpolation location is in cell i0,j0,k0
              w(1,2) = w(1,2)*sigC(CondE,IX,currSG,i0+1,j0,k0)/    &
			sigC(CondE,IX,currSG,i0,j0,k0)
           else
           ! interpolation location is in cell i0+1,j0,k0
              w(1,1) = w(1,1)*sigC(CondE,IX,currSG,i0,j0,k0)/    &
			sigC(CondE,IX,currSG,i0+1,j0,k0)
           endif
        elseif(xyz.eq.2) then
           if(jc.eq.j0) then
           !    interpolation location is in cell i0,j0,k0
              w(2,2) = w(2,2)*sigC(CondE,IY,currSG,i0,j0+1,k0)/    &
			sigC(CondE,IY,currSG,i0,j0,k0)
           else
           ! interpolation location is in cell i0,j0+1,k0
              w(2,1) = w(2,1)*sigC(CondE,IY,currSG,i0,j0,k0)/   &
			sigC(CondE,IY,currSG,i0,j0+1,k0)
           endif
        elseif(xyz.eq.3) then
           if(kc.eq.k0) then
           !    interpolation location is in cell i0,j0,k0
              w(3,2) = w(3,2)*sigC(CondE,IZ,currSG,i0,j0,k0+1)/  &
			sigC(CondE,IZ,currSG,i0,j0,k0)
           else
           ! interpolation location is in cell i0,j0,k0+1
              w(3,1) = w(3,1)*sigC(CondE,IZ,currSG,i0,j0,k0)/   &
			sigC(CondE,IZ,currSG,i0,j0,k0+1)
           endif
        endif
    endif
    ii = 0
    do n = 1,2
       do m = 1,2
          do p = 1,2
             wadd = w(1,n)*w(2,m)*w(3,p)
             if(wadd .gt. 0) then
                ii = ii + 1
                I(ii) = i0+n-1
                J(ii) = j0+m-1
                K(ii) = k0SG+p-1
                C(ii) = wadd
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
    LC%currSG = currSG

  end subroutine EinterpSetUp


  ! **************************************************************************
  ! magnetic field coefficients in the sparse vector
  subroutine BinterpSetUp(inGrid,x,xyz,LC)

    ! 24.09.2012   Maria
    ! This routine is re-written for Multigrid
! ___________________________________________________________________________________
    ! sets up coefficients in multi-grid sparse vector LC for evaluation/interpolation
    ! of magnetic field vector component xyz (1 = x; 2 = y; 3 = z) at
    ! location given by x, using MAGNETIC field vector defined on staggered
    ! grid faces.  (For direct application, magnetic fields defined on faces
    ! would be required; can be used to construct an interpolator to compute
    ! H at arbitrary points directly from electric fields defined on
    ! staggered grid edges; see BfromESetup)
    ! INTERPOLATION METHOD:  bilinear spline
    ! NOTE: this is essentially identical to EinterpSetUp,
    ! except usage of xEdge and xCener are reversed
!_____________________________________________________________________________________

    implicit none
    type (grid_t), intent(in)    		:: inGrid
    real (kind=prec), dimension(3), intent(in) 	:: x
    integer, intent(in)        			:: xyz
    type (sparsevecc), intent(inout) 		:: LC

    integer 					:: i0,j0,k0,ii,n,m,p
    integer                                     :: nxMax, nyMax, nzMax
    integer, dimension(8)			:: I,J,K
    integer					:: status
    real (kind=prec), dimension(3,2)	:: w
    real (kind=prec)			:: wadd
    real (kind=prec)            :: x3SG
    complex (kind=prec), dimension(8)  	:: C
    character (len=80)				:: gridType = ''
    integer                         :: izz,izc,imgrid,nzCum
    integer                         ::currSG
    logical                         ::found

    if(LC%allocated) then
      call deall_sparsevecc(LC)
    endif

  ! which sub-grid x(3) belong to?  ! x3SD
  ! derive the number of this sub-grid  ! currSG
    nzCum = 0
    Subgrids: do imgrid = 1, inGrid%mgridSize
      Zlayers: do izz = 1, inGrid%gridArray(imgrid)%nz
        izc = izz+nzCum
        if(clean(inGrid%zEdge(izc)) .gt. clean(x(3)) ) then
          currSG = imgrid
          if(izz == 1) then
            x3SG = inGrid%gridArray(imgrid-1)%zEdge(inGrid%gridArray(imgrid-1)%nz+1)
          else
            x3SG = inGrid%gridArray(imgrid)%zEdge(izz-1)
          endif
          found = .true.
        exit
        endif
      enddo Zlayers
      if(found == .true.) exit
      nzCum = nzCum + inGrid%gridArray(imgrid)%nz
    enddo Subgrids

    ! maximum number of center nodes
    nxMax = inGrid%gridArray(currSG)%nx
    nyMax = inGrid%gridArray(currSG)%ny
    nzMax = inGrid%gridArray(currSG)%nz

    if (xyz .eq. 1) then
       ! Evaluation of x component wrt magnetic vectors on cubic faces
       i0 = minNode(x(1),inGrid%gridArray(currSG)%xEdge)
       j0 = minNode(x(2),inGrid%gridArray(currSG)%yCenter)
       k0 = minNode(x3SG,inGrid%gridArray(currSG)%zCenter)
       ! maximum number of edge nodes
       nxMax = nxMax + 1
    elseif (xyz .eq. 2) then
       ! Evaluation of y component wrt magnetic vectors on cubic faces
       i0 = minNode(x(1),inGrid%gridArray(currSG)%xCenter)
       j0 = minNode(x(2),inGrid%gridArray(currSG)%yEdge)
       k0 = minNode(x3SG,inGrid%gridArray(currSG)%zCenter)
       ! maximum number of edge nodes
       nyMax = nyMax + 1
    elseif (xyz .eq. 3) then
       ! Evaluation of z component wrt magnetic vectors on cubic faces
       i0 = minNode(x(1),inGrid%gridArray(currSG)%xCenter)
       j0 = minNode(x(2),inGrid%gridArray(currSG)%yCenter)
       k0 = minNode(x3SG,inGrid%gridArray(currSG)%zEdge)
       ! maximum number of edge nodes
       nzMax = nzMax + 1
    else
       write(0,*) 'Error: component # out of range in BinterpSetUp'
    endif

    if((i0.gt.0).and.(i0.lt.nxMax)) then
       if(xyz.eq.1) then
          w(1,2) = (x(1) - inGrid%gridArray(currSG)%xEdge(i0))/(inGrid%gridArray(currSG)%dx(i0))
       else
          w(1,2) = (x(1) - inGrid%gridArray(currSG)%xCenter(i0))/(inGrid%gridArray(currSG)%delX(i0+1))
       endif
    elseif(i0.le.0) then
       w(1,2) = 1
    else
       w(1,2) = 0
    endif

    if((j0.gt.0).and.(j0.lt.nyMax)) then
       if(xyz.eq.2) then
          w(2,2) = (x(2) - inGrid%gridArray(currSG)%yEdge(j0))/(inGrid%gridArray(currSG)%dy(j0))
       else
          w(2,2) = (x(2) - inGrid%gridArray(currSG)%yCenter(j0))/(inGrid%gridArray(currSG)%delY(j0+1))
       endif
    elseif(j0.le.0) then
       w(2,2) = 1
    else
       w(2,2) = 0
    endif

    if((k0.gt.0).and.(k0.lt.nzMax)) then
       if(xyz.eq.3) then
          w(3,2) = (x3SG - inGrid%gridArray(currSG)%zEdge(k0))/(inGrid%gridArray(currSG)%dz(k0))
       else
          w(3,2) = (x3SG - inGrid%gridArray(currSG)%zCenter(k0))/(inGrid%gridArray(currSG)%delZ(k0+1))
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
                C(ii) = wadd
             endif
          enddo
       enddo
    enddo

    ! we are dealing with magnetic fields (therefore, gridtype = FACE)
    gridType = FACE
    Call create_sparsevecc(ii,LC,gridType)

    ! The same as for E
     do n=1,ii
      LC%i(n)=I(n)
      LC%j(n)=J(n)
      LC%k(n)=K(n)
      LC%c(n)=C(n)
     end do    
    !   assuming xyz will be assigned to all elements of LC%xyz
    LC%xyz = xyz
    LC%currSG = currSG

  end subroutine BinterpSetUp


  ! **************************************************************************
  ! magnetic field from electrical field in a sparse vector data structures
  subroutine BfromESetUp(inGrid,omega,x,xyz,LC)
  ! Multi-grid version of the routine
  ! 24.09.2012 Maria
! ________________________________________________________________________________
    !  sets up coefficients in sparse vector LC for evaluation/interpolation
    !  of magnetic field vector component xyz (1 = x; 2 = y; 3 = z) at
    !  location given by x, using ELECTRIC field vector defined on staggered
    !  grid edges.  Calls BinterpSetUp, and various sparse_vector routines
!_________________________________________________________________________________

    implicit none
    type (grid_t), target, intent(in) 		:: inGrid
    real(kind=prec), 	intent(in)	:: omega
    real(kind=prec), dimension(3), intent(in)	:: x
    integer, intent(in)              		:: xyz
    type (sparsevecc), intent(inout) 		:: LC

    ! local variables
    integer 					:: ii,n,m,p
    integer					:: num = 4
    integer, dimension(4)			:: I,J,K,AXES
    real (kind=prec), dimension(3,2)   	:: w
    complex (kind=prec), dimension(4)  	:: C
    complex (kind=prec)   	       	:: c1,c2,i_omega_inv
    type (sparsevecc)				:: LConeH, LCtemp
    type (sparsevecc)            :: LCH
    character (len=80)				:: gridType = ''
    integer					:: status

    !   NOTE: Careful about minus sign here ...
    i_omega_inv = ISIGN*cmplx(0.0 ,1./omega,kind=prec)

    if(LC%allocated) then
      call deall(LC) ! deall_sparsevecc
    endif

    ! we work with electrical fields here, therefore gridType = EDGE
    gridType = EDGE
    ! compute coefficients for bi-linear spline interpolation of
    ! magnetic fields from grid cell faces to location x
    Call BinterpSetUp(inGrid,x,xyz,LCH)
    Call create(num,LC,gridType) !create_sparsevecc

    !   loop over coefficients for mag field interpolation
    do ii = 1,LCH%nCoeff

       if(LCH%xyz(ii).eq.1) then
          AXES(1) = 2
          AXES(2) = 2
          AXES(3) = 3
          AXES(4) = 3
          I(1) = LCH%i(ii)
          I(2) = LCH%i(ii)
          I(3) = LCH%i(ii)
          I(4) = LCH%i(ii)
          J(1) = LCH%j(ii)
          J(2) = LCH%j(ii)
          J(3) = LCH%j(ii)
          J(4) = LCH%j(ii)+1
          K(1) = LCH%k(ii)
          K(2) = LCH%k(ii)+1
          K(3) = LCH%k(ii)
          K(4) = LCH%k(ii)
          C(1) = 1./inGrid%gridArray(LCH%currSG)%dz(K(1))
	  C(2) = -1./inGrid%gridArray(LCH%currSG)%dz(K(1))
          C(3) = -1./inGrid%gridArray(LCH%currSG)%dy(J(1))
          C(4) = 1./inGrid%gridArray(LCH%currSG)%dy(J(1))

       elseif(LCH%xyz(ii).eq.2) then
          AXES(1) = 3
          AXES(2) = 3
          AXES(3) = 1
          AXES(4) = 1
          J(1) = LCH%j(ii)
          J(2) = LCH%j(ii)
          J(3) = LCH%j(ii)
          J(4) = LCH%j(ii)
          K(1) = LCH%k(ii)
          K(2) = LCH%k(ii)
          K(3) = LCH%k(ii)
          K(4) = LCH%k(ii)+1
          I(1) = LCH%i(ii)
          I(2) = LCH%i(ii)+1
          I(3) = LCH%i(ii)
          I(4) = LCH%i(ii)
          C(1) = 1./inGrid%gridArray(LCH%currSG)%dx(I(1))
          C(2) = -1./inGrid%gridArray(LCH%currSG)%dx(I(1))
          C(3) = -1./inGrid%gridArray(LCH%currSG)%dz(K(1))
          C(4) = 1./inGrid%gridArray(LCH%currSG)%dz(K(1))

       else
          AXES(1) = 1
          AXES(2) = 1
          AXES(3) = 2
          AXES(4) = 2
          K(1) = LCH%k(ii)
          K(2) = LCH%k(ii)
          K(3) = LCH%k(ii)
          K(4) = LCH%k(ii)
          I(1) = LCH%i(ii)
          I(2) = LCH%i(ii)
          I(3) = LCH%i(ii)
          I(4) = LCH%i(ii)+1
          J(1) = LCH%j(ii)
          J(2) = LCH%j(ii)+1
          J(3) = LCH%j(ii)
          J(4) = LCH%j(ii)
          C(1) = 1./inGrid%gridArray(LCH%currSG)%dy(J(1))
          C(2) = -1./inGrid%gridArray(LCH%currSG)%dy(J(1))
          C(3) = -1./inGrid%gridArray(LCH%currSG)%dx(I(1))
          C(4) = 1./inGrid%gridArray(LCH%currSG)%dx(I(1))
       endif

       if(ii.eq.1) then
          ! initialize LC (coefficients for H measurement functional:
          ! coefficients for first face
          Call create(num,LConeH,gridType) !create_sparsevecc
          Call create(num,LCtemp,gridType) !create_sparsevecc
          LC%i = I
          LC%j = J
          LC%k = K
          LC%xyz = AXES
          LC%c= C*LCH%c(1)*i_omega_inv
       else
          ! add coefficients for next face
          ! first store cumulative sum so far in LCtemp
	      LCtemp = LC
          LConeH%i = I
          LConeH%j = J
          LConeH%k = K
          LConeH%xyz = AXES
          LConeH%c = C

          c1 = C_ONE
          c2 = LCH%c(ii)*i_omega_inv
          Call linComb(LCtemp,c1,LConeH,c2,LC) !linComb_sparsevecc
       endif

    enddo
    LC%currSG = LCH%currSG
    ! do all the clean-up while exiting
    call deall(LConeH) !deall_sparsevecc
    call deall(LCH)
    call deall(LCtemp)

  end subroutine BfromESetUp

end module EMfieldInterp
