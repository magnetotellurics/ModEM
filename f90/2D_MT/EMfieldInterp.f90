
! *****************************************************************************
module EMfieldInterp
  ! Generic interpolation functionals for both electric and magnetic
  ! fields to an arbitrary point within the model domain,
  ! including routines for TE electric field, and TM magnetic field
  ! solutions.  These routines just create sparse vector representations
  ! of these basic interpolation functionals.  Routines in this module
  ! are called by higer level routines to actually apply data functionals,
  ! construct "combs" for inversion, etc.

  use math_constants
  use utilities
  use EMfield
  use ModelSpace, rhoC => ModelParamToOneCell

  implicit none
  public				:: NodeInterpSetup2D
  public				:: BinterpSetUp_TE, BfromESetUp_TE
  public				:: EinterpSetUp_TM, EfromBSetUp_TM

Contains

  ! **************************************************************************
  subroutine NodeInterpSetup2D(inGrid,x,mode,LC)
    ! sets up coefficients in sparse vector LC for evaluation/interpolation
    ! at location x of a field defined on TE or TM mode grid
    ! zero vertical location corresponds to the Earth's surface; positive down
    ! INPUT PARAMETER mode (character*2) is used to define which case
    !   to do the interpolation for (Only difference is that there is
    !   no air layer in the TM mode solution, so vertical node numbering
    !   is slightly different)
    ! INTERPOLATION METHOD: bi-linear spline
    ! Allocates arrays in sparse vector LC, deallocating first if necessary
    !   This can be used to setup for interpolation of E for TE solution
    !     or for B for TM solution

    !  on at least some compilers need to use intent(inout) for pointer
    !    targets????
    type (grid_t), target, intent(in)	:: inGrid
    real (kind=prec), dimension(2), intent(in)     :: x
    character (len=2), intent(in)     		:: mode  ! = 'TE' or 'TM'
    type (sparsevecc), intent(inout)          :: LC

    ! Local Variables
    integer                                     :: j0,k0,ii,n,m
    integer, dimension(4)                       :: J,K
    integer                                     :: nyMax, nzMax
    complex(kind=prec), dimension(4)	:: C
    real(kind=prec), dimension(2,2)	:: w
    real(kind=prec)			:: wadd

     if(LC%allocated) then
        call deall_sparsevecc(LC)
     endif

    ! maximum number of edge nodes
    nyMax = inGrid%ny+1
    nzMax = inGrid%nz+1

    j0 = minNode(x(1),inGrid%yNode)
    k0 = minNode(x(2)+inGrid%zAir,inGrid%zNode)
    if((j0.gt.0).and.(j0.lt.nyMax)) then
       w(1,2) = (x(1) - inGrid%yNode(j0))/(inGrid%Dy(j0))
    elseif(j0.le.0) then
       w(1,2) = ONE
    else
       w(1,2) = R_ZERO
    endif

    if((k0.gt.0).and.(k0.lt.nzMax)) then
       w(2,2) = (x(2)+inGrid%zAir - inGrid%zNode(k0))/(inGrid%dz(k0))
    elseif(k0.le.0) then
       w(2,2) = ONE
    else
       w(2,2) = R_ZERO
    endif

    w(1,1) = ONE-w(1,2)
    w(2,1) = ONE-w(2,2)

    ii = 0
    do n = 1,2
       do m = 1,2
          wadd = w(1,n)*w(2,m)
          if(wadd .gt. 0) then
             ii = ii + 1
             J(ii) = j0+n-1
             K(ii) = k0+m-1
             C(ii) = wadd
          endif
       enddo
    enddo

    if(mode == 'TE') then
       call create_sparsevecc(inGrid,NODE,ii,LC)
    else
       call create_sparsevecc(inGrid,NODE_EARTH,ii,LC)
    endif

    LC%j = J(1:ii)
    LC%k = K(1:ii)
    LC%c = C(1:ii)
    if(mode == 'TM') then
       LC%k = K(1:ii) - inGrid%Nza
    endif

  end subroutine NodeInterpSetup2D
  ! **************************************************************************
  subroutine BinterpSetUp_TE(inGrid,x,LC)
    ! sets up coefficients in sparse vector LC for evaluation/interpolation
    ! of HORIZONTAL magnetic field component at location given by x,
    ! using MAGNETIC field vector defined on face
    ! zero vertical location corresponds to the Earth's surface; positive down
    ! (For direct application, magnetic fields defined on faces
    ! would be required; can be used to construct an interpolator to compute
    ! H at arbitrary points directly from electric fields defined on
    ! staggered grid edges; see BfromESetup_TE)
    ! INTERPOLATION METHOD:  bilinear spline
    !   NEED TO MODIFY TO ALLOW INTERPOLATION OF VERTICAL

    type (grid_t), target, intent(in)	:: inGrid
    real (kind=prec), dimension(2), intent(in) 	:: x
    type (sparsevecc), intent(inout) 		:: LC

    !   local variables
    integer 					:: j0,k0,ii,n,m
    integer                                     :: nyMax, nzMax
    integer, dimension(4)			:: J,K
    real (kind=prec), dimension(2,2)	:: w
    real (kind=prec)			:: wadd
    complex (kind=prec), dimension(4)	:: C

    if(LC%allocated) then
       call deall_sparsevecc(LC)
    endif

    ! this part needs to change for vertical mag field
    nyMax = inGrid%ny+1
    nzMax = inGrid%nz

    j0 = minNode(x(1),inGrid%yNode)
    if((j0.gt.0).and.(j0.lt.nyMax)) then
       w(1,2) = (x(1) - inGrid%yNode(j0))/(inGrid%dy(j0))
    elseif(j0.le.0) then
       w(1,2) = ONE
    else
       w(1,2) = R_ZERO
    endif

    k0 = minNode(x(2)+inGrid%zAir,inGrid%zCenter)
    if((k0.gt.0).and.(k0.lt.nzMax)) then
       w(2,2) = (x(2)+inGrid%zAir - inGrid%zCenter(k0))/(inGrid%Delz(k0+1))
    elseif(k0.le.0) then
       w(2,2) = ONE
    else
       w(2,2) = R_ZERO
    endif

    w(1,1) = ONE-w(1,2)
    w(2,1) = ONE-w(2,2)

    ii = 0
    do n = 1,2
       do m = 1,2
          wadd = w(1,n)*w(2,m)
          if(wadd .gt. 0) then
             ii = ii + 1
             J(ii) = j0+n-1
             K(ii) = k0+m-1
             C(ii) = wadd
          endif
       enddo
    enddo

    call create_sparsevecc(inGrid,EDGE,ii,LC)
    LC%j = J(1:ii)
    LC%k = K(1:ii)
    LC%c = C(1:ii)

  end subroutine BinterpSetUp_TE
  ! **************************************************************************
  ! magnetic field from electrical field in a sparse vector data structures
  subroutine BfromESetUp_TE(inGrid,x,omega,LC)
    !  sets up coefficients in sparse vector LC for evaluation/interpolation
    !  of HORIZONTAL TE magnetic field xyz at
    !  location given by x, using TE ELECTRIC field defined on grid nodes
    !  zero vertical location corresponds to the Earth's surface; positive down
    !  Calls BinterpSetUp, and various sparse_vector routines
    !   NEED TO MODIFY TO ALLOW INTERPOLATION OF VERTICAL

    type (grid_t), target, intent(in)	:: inGrid
    real (kind=prec), dimension(2), intent(in)	:: x
    real (kind=prec), intent(in)	:: omega
    type (sparsevecc), intent(inout) 		:: LC

    integer 					:: ii
    integer					:: num = 2
    integer, dimension(2)			:: J,K
    complex (kind=prec), dimension(2)	:: C
    complex (kind=prec)			:: i_omega
    type (sparsevecc)				:: LConeH, LCH, LCtemp

    i_omega = MinusONE*ISIGN*cmplx(0.0 ,omega,kind=prec)

    if(LC%allocated) then
       call deall_sparsevecc(LC)
    endif
    call create_sparsevecc(inGrid,NODE,num,LC)

    !  first setup sparsevector for interpolation of mag fields from
    !  mag fields at center of edge
    call BinterpSetUp_TE(inGrid,x,LCH)

    !   loop over coefficients for mag field interpolation
    do ii = 1,LCH%nCoeff
       J(1) = LCH%j(ii)
       J(2) = LCH%j(ii)
       K(1) = LCH%k(ii)
       K(2) = LCH%k(ii)+1
       C(1) = MinusONE/(inGrid%dz(K(1))*i_omega)
       C(2) = ONE/(inGrid%dz(K(1))*i_omega)

       if(ii.eq.1) then
          ! initialize LC (coefficients for H measurement functional:
          ! coefficients for first face
          call create_sparsevecc(inGrid,NODE,num,LConeH)
          LC%j = J
          LC%k = K
          LC%c= C*LCH%c(1)
       else
          ! add coefficients for next face
          ! first store cumulative sum so far in LCtemp
	  call copy_sparsevecc(LCtemp,LC)
          LConeH%j = J
          LConeH%k = K
          LConeH%c = C
          call linComb_sparsevecc(LCtemp,C_ONE,LConeH,LCH%c(ii),LC)
       endif

    enddo

    ! clean-up
    call deall_sparsevecc(LConeH)
    call deall_sparsevecc(LCH)
    call deall_sparsevecc(LCtemp)

  end subroutine BfromESetUp_TE
  ! **************************************************************************
  subroutine EinterpSetUp_TM(inGrid,x,sigma,LC,Q)
    ! sets up coefficients in sparse vector LC for evaluation/interpolation
    ! of HORIZONTAL electric field component at location given by x,
    ! using ELECTRIC field vector defined on faces,  or
    !    on surface magnetic node
    ! zero vertical location corresponds to the Earth's surface; positive down
    ! (For direct application, electric fields defined on faces/surface nodes
    ! would be required; can be used to construct an interpolator to compute
    ! E at arbitrary points directly from magnetic fields defined on
    ! staggered grid edges; see EfromBsetup_TM)
    ! ALSO, optionally returns array of sparse vectors Q defined on cells for
    !   construction of derivatives of data functional derivaties with
    !   respect to resistivity parameters; one sparse vector is returned
    !    for each interpolation coefficient in LC
    !
    ! INTERPOLATION METHOD:  bilinear spline

    type (grid_t), target, intent(in)        :: inGrid
    real (kind=prec), dimension(2), intent(in)    :: x
    type (modelParam_t), intent(in)		:: sigma
    type (sparsevecc), intent(inout)		:: LC
    type (sparsevecc), intent(inout),optional	:: Q(4)

    !   local variables
    integer                                     :: j0,k0,ii,n,m,j0m1,j0p1
    integer                                     :: nyMax, nzMax
    integer, dimension(4)                       :: J,K
    complex (kind=prec), dimension(4)	:: C
    integer, dimension(2)                       :: JQ,KQ
    complex (kind=prec), dimension(2)	:: CQ
    real (kind=prec), dimension(2,2)	:: w
    real (kind=prec), dimension(2)	:: r
    real (kind=prec)			:: wadd
    real (kind=prec), dimension(:), allocatable	:: eNode,DelZ
    logical 					:: returnQ


    returnQ = present(Q)
    if(LC%allocated) then
       call deall_sparsevecc(LC)
    endif
    if(returnQ) then
       do ii = 1,4
          if(Q(ii)%allocated) then
             call deall_sparsevecc(Q(ii))
          endif
       enddo
    endif

    nyMax = inGrid%ny+1
    nzMax = inGrid%nz - inGrid%Nza+1

    !   construct a new "eNode" list:  first vertical node is at surface,
    !    subsequent vertical nodes are at face centers
    allocate(eNode(nzMax))
    allocate(DelZ(nzMax-1))
    eNode(1) = inGrid%zNode(inGrid%Nza+1)
    eNode(2:nzMax) = inGrid%zCenter(inGrid%Nza+1:inGrid%nz)
    DelZ = eNode(2:nzMax)-eNode(1:nzMax-1)

    !  vertical interpolation, using eNode and DelZ
    k0 = minNode(x(2)+inGrid%zAir,eNode)
    if((k0.gt.0).and.(k0.lt.nzMax)) then
       w(2,2) = (x(2)+inGrid%zAir - eNode(k0))/(DelZ(k0))
    elseif(k0.le.0) then
       w(2,2) = ONE
    else
       w(2,2) = R_ZERO
    endif
    w(2,1) = ONE-w(2,2)

    ! lateral interpolation, coefficients depend on resistivity
    j0 = minNode(x(1),inGrid%yNode)
    if((j0.gt.0).and.(j0.lt.nyMax)) then
       w(1,2) = (x(1) - inGrid%yNode(j0))/(inGrid%dy(j0))
    elseif(j0.le.0) then
       w(1,2) = ONE
    else
       w(1,2) = R_ZERO
    endif
    w(1,1) = ONE-w(1,2)

    !  vertical cell index for interpolated point
    k0 = minNode(x(2)+inGrid%zAir,inGrid%zNode) - inGrid%Nza
    j0m1 = max(1,j0-1)
    j0p1 = min(j0+1,ingrid%ny)
    r(1) = TWO*rhoC(sigma,j0,k0)/(rhoC(sigma,j0,k0)+rhoC(sigma,j0m1,k0))
    r(2) = TWO*rhoC(sigma,j0,k0)/(rhoC(sigma,j0,k0)+rhoC(sigma,j0p1,k0))

    ii = 0
    do n = 1,2
       do m = 1,2
          wadd = w(1,n)*w(2,m)*r(n)
          if(wadd .gt. 0) then
             ii = ii + 1
             J(ii) = j0+n-1
             !  NOTE: k0=0 corresponds to surface node
             K(ii) = k0+m-2
             C(ii) = wadd

             if(returnQ) then
                if((j0m1 .lt. j0).and.(j0p1.gt.j0)) then
                   call create_sparsevecc(inGrid,CELL,2,Q(ii))
                   if(n.eq.1) then
                      JQ(1) = j0m1
                      JQ(2) = j0
                   else
                      JQ(1) = j0p1
                      JQ(2) = j0
                   endif
                   KQ(1) = max(K(ii),1)
                   KQ(2) = KQ(1)
                   CQ(1) = -r(n)*w(1,n)*w(2,m)/&
			(rhoC(sigma,JQ(1),KQ(1))+rhoC(sigma,JQ(2),KQ(1)))
                   CQ(2) = (TWO-r(n))*w(1,n)*w(2,m)/&
			(rhoC(sigma,JQ(1),KQ(1))+rhoC(sigma,JQ(2),KQ(1)))
                   Q(ii)%j = JQ
                   Q(ii)%k = KQ
                   Q(ii)%c = CQ
                else
                   call errStop('obs location out of bounds in EinterpSetUp_TM')
                   !  need to come up with something simple for this case
                   !   (trying to evaluate E in the outer edge cells)
                endif
             endif
          endif
       enddo
    enddo

    call create_sparsevecc(inGrid,EDGE_EARTH,ii,LC)
    LC%j = J(1:ii)
    LC%k = K(1:ii)
    LC%c = C(1:ii)
    LC%grid => inGrid

    deallocate(DelZ)
    deallocate(eNode)

  end subroutine EinterpSetUp_TM

! **************************************************************************
  ! magnetic field from electrical field in a sparse vector data structures
  subroutine EfromBSetUp_TM(inGrid,x,omega,sigma,LC,b0,Q)
    !  sets up coefficients in sparse vector LC for evaluation/interpolation
    !  of HORIZONTAL TM electric field at
    !  location given by x, using TM Magetic field defined on grid nodes
    !  zero vertical location corresponds to the Earth's surface; positive down
    !  Calls EinterpSetUp_TM, and various sparse_vector routines

    !  Input: model grid data structure
    type (grid_t), target, intent(in)	:: inGrid
    !  Input: location of point observation of electric field
    real (kind=prec), dimension(2), intent(in)	:: x
    !  Input: frequency
    real (kind=prec), intent(in)	:: omega
    type (modelParam_t), intent(in)		:: sigma
    ! Output: sparse vector defined on magnetic field solution space
    type (sparsevecc), intent(inout) 		:: LC
    ! Opitional arguments, used for derivative of data due
    !     to dependence of data functional coefficients on model
    !     parameters.  BOTH ARE REQUIRED IF ONE IS PROVIDED
    !  Optional input: background magnetic field
    type (cvector), intent(in),optional 		:: b0
    !  Optional output: sparse vector
    !      defined on conductivity parameter space
    type (sparsevecc), intent(inout),optional :: Q

    ! local variables
    integer 					:: ii,Nza,nCoeff
    integer					:: num = 6
    integer, dimension(4)			:: J,K
    complex (kind=prec), dimension(4)	:: C
    complex (kind=prec)			:: e_ii,i_omega
    type (sparsevecc)				:: LCE, LC_ii, LCtemp, Qj
    type (sparsevecc)				:: QE(4),Qtemp
    logical 					:: returnQ = .false.
    logical,parameter          			:: Conj_Case = .false.
    real (kind=prec)			:: dz1,dy1,dy2

    i_omega = MinusONE*ISIGN*cmplx(0.,omega,kind=prec)

    if(LC%allocated) then
       call deall_sparsevecc(LC)
    endif

    !  first setup sparse vectors for interpolation of electric fields from
    !  electric fields defined on faces + surface B node
    returnQ = present(Q)
    if(returnQ) then
       call EinterpSetUp_TM(inGrid,x,sigma,LCE,QE)
    else
       call EinterpSetUp_TM(inGrid,x,sigma,LCE)
    endif

    ! then loop over coefficients in electric field interpolation
    Nza = inGrid%Nza
    do ii = 1,LCE%nCoeff
       if(LCE%K(ii).eq.0) then
          !  surface E node
          !   (only this case occurs for standard TM mode)
          dz1 = inGrid%dz(Nza+1)
	  nCoeff = 4
          J(1:2) = LCE%J(ii)
          J(3) = LCE%J(ii)-1
          J(4) = LCE%J(ii)+1
          K(1) = 1
          K(2:4) = 2
          C(1) = -(rhoC(sigma,J(1)-1,1)+rhoC(sigma,J(1),1))/(TWO*dz1*MU_0) + &
			i_omega*THREE*dz1/EIGHT
          dy1 = inGrid%dy(J(1)-1)*(inGrid%dy(J(1)-1)+inGrid%dy(J(1)))*MU_0
          dy2 = inGrid%dy(J(1))*(inGrid%dy(J(1)-1)+inGrid%dy(J(1)))*MU_0
          C(2) = (rhoC(sigma,J(1)-1,1)+rhoC(sigma,J(1),1))/(TWO*dz1*MU_0) &
			+ i_omega*dz1/EIGHT  &
          		- dz1*(rhoC(sigma,J(1)-1,1)&
			+rhoC(sigma,J(1)-1,2))/(EIGHT*dy1) &
          		- dz1*(rhoC(sigma,J(1),1)+ &
			rhoC(sigma,J(1),2))/(EIGHT*dy2)
          C(3) = (rhoC(sigma,J(1)-1,1)+ &
		rhoC(sigma,J(1)-1,2))*dz1/(EIGHT*dy1)
          C(4) = (rhoC(sigma,J(1),1)+ &
		rhoC(sigma,J(1),2))*dz1/(EIGHT*dy2)

          if(returnQ) then
             ! make sparse vector Qj defined on cells
             call create_sparsevecc(inGrid,CELL,4,Qj)

             ! First cell
             Qj%J(1) = J(1) - 1
             Qj%K(1) =  1
             call create_sparsevecc(inGrid,NODE_EARTH,3,LCtemp)
             LCtemp%J = J(1:3)
             LCtemp%K = K(1:3)
             LCtemp%C(1) = MinusONE/(TWO*dz1*MU_0)
             LCtemp%C(2) = ONE/(TWO*dz1*MU_0) - dz1/(EIGHT*dy1)
             LCtemp%C(3) = dz1/(EIGHT*dy1)
             Qj%C(1) = dotProd_scvector(LCtemp,b0,Conj_Case)

             ! Second cell
             Qj%J(2) = J(1)
             Qj%K(2) =  1
             LCtemp%J(1:2) = J(1:2)
             LCtemp%J(3) = J(4)
             LCtemp%K(1:2) = K(1:2)
             LCtemp%K(3) = K(4)
             LCtemp%C(1) = MinusONE/(TWO*dz1*MU_0)
             LCtemp%C(2) = ONE/(TWO*dz1*MU_0) - dz1/(EIGHT*dy2)
             LCtemp%C(3) = dz1/(EIGHT*dy2)
             Qj%C(2) = dotProd_scvector(LCtemp,b0,Conj_Case)

             ! Third cell
             Qj%J(3) = J(1)-1
             Qj%K(3) =  2
             call create_sparsevecc(inGrid,NODE_EARTH,2,LCtemp)
             LCtemp%J(1:2) = J(2:3)
             LCtemp%K(1:2) = K(2:3)
             LCtemp%C(1) = - dz1/(EIGHT*dy1)
             LCtemp%C(2) =  dz1/(EIGHT*dy1)
             Qj%C(3) = dotProd_scvector(LCtemp,b0,Conj_Case)

             ! Fourth cell
             Qj%J(4) = J(1)
             Qj%K(4) =  2
             LCtemp%J(1) = J(2)
             LCtemp%J(2) = J(4)
             LCtemp%K(1) = K(2)
             LCtemp%K(2) = K(4)
             LCtemp%C(1) = - dz1/(EIGHT*dy2)
             LCtemp%C(2) =  dz1/(EIGHT*dy2)
             Qj%C(4) = dotProd_scvector(LCtemp,b0,Conj_Case)
          endif

       else
          !  not a surface E node
          dz1 = inGrid%dz(K(1)+ingrid%Nza)
          nCoeff = 2
          J(1:2) = LCE%J(ii)
          K(1) = LCE%K(ii)
          K(2) = LCE%K(ii) + 1
          C(2) = (rhoC(sigma,J(1)-1,K(1))+rhoC(sigma,J(1),K(1)))&
				/(TWO*dz1*MU_0)
          C(1) = -C(2)

          if(returnQ) then
             ! make sparse vector Qj defined on cells
             call create_sparsevecc(inGrid,CELL,2,Qj)

             ! First cell
             Qj%J(1) = J(1) - 1
             Qj%K(1) =  K(1)
             call create_sparsevecc(inGrid,NODE_EARTH,2,LCtemp)
             LCtemp%J = J
             LCtemp%K = K
             LCtemp%C(1) = MinusONE/(TWO*dz1*MU_0)
             LCtemp%C(2) = ONE/(TWO*dz1*MU_0)
             Qj%C(1) = dotProd_scvector(LCtemp,b0,Conj_Case)

             ! Second cell
             Qj%J(2) = J(1)
             Qj%K(2) =  K(1)
             Qj%C(2) = Qj%C(1)
          endif

       endif
       !  Now merge into one sparse vector L for the solution space
       !    (and optionally a sparse vector Q for the parameter space)
       if(ii.eq.1) then
          call create_sparsevecc(inGrid,NODE_EARTH,nCoeff,LC)
          LC%J = J(1:nCoeff)
          LC%K = K(1:nCoeff)
          LC%C = C(1:nCoeff)
          LC%grid => inGrid
          if(returnQ) then
             call copy_sparsevecc(Qtemp,Qj)
             !   evaluate electric field at node ii
             e_ii = dotProd_scvector(LC,b0,Conj_Case)
             call linComb_sparsevecc(Qtemp,LCE%C(1),QE(1),e_ii,Q)
             call deall_sparsevecc(QE(1))
          endif
	  LC%C = LC%C*LCE%C(1)
       else
         !   WHAT SHOULD gridType be HERE??????  Guessing NODE_EARTH
         !   NEED TO CHECK !!!!!!!
          call create_sparsevecc(inGrid,NODE_EARTH,nCoeff,LC_ii)
          ! add coefficients from next face to LC
          call copy_sparsevecc(LCtemp,LC)
          LC_ii%J  = J(1:nCoeff)
          LC_ii%K  = K(1:nCoeff)
          LC_ii%C = C(1:nCoeff)
          call linComb_sparsevecc(LCtemp,C_ONE,LC_ii,LCE%C(ii),LC)
          if(returnQ) then
	     call copy_sparsevecc(Qtemp,Q)
             call linComb_sparsevecc(Qtemp,C_ONE,Qj,LCE%C(ii),Q)
	     call copy_sparsevecc(Qtemp,Q)
             e_ii = dotProd_scvector(LC_ii,b0,Conj_Case)
             call linComb_sparsevecc(Qtemp,C_ONE,QE(ii),e_ii,Q)
             call deall_sparsevecc(QE(ii))
          endif
       endif
    enddo
    call deall_sparsevecc(LC_ii)
    call deall_sparsevecc(LCE)
    call deall_sparsevecc(LCtemp)
    call deall_sparsevecc(Qj)
    call deall_sparsevecc(Qtemp)

  end subroutine EfromBSetUp_TM
end module EMfieldInterp
