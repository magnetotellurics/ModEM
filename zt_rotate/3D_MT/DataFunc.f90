! *****************************************************************************
module dataFunc
  ! 3D MT data functionals

  ! This module contains
  !  (1) routines for evaluation of impedances, and ultimately other
  !       interpretation parameters
  !  (2) routines to compute data functionals for linearized
  !       impedances,  and ultimately other interpretation paramters
  !   The idea:
  !     -> first the dictionaries txDict, typeDict, rxDict are initialized
  !         by calling appropriate initialization/setup routines
  !     -> data are stored in structures (defined in module DataSpace)
  !        which contain indices into transmitter and receiver dictionaries
  !        in addition to actual data values.  These indices are used by
  !        the data functional computation routines to compute predicted data.
  !
  !
  !  This module is specific to 3D MT; similar modules will need to be written
  !     to implement data functionals for other problems

  use EMfieldInterp
  use SolnSpace
  use receivers
  use transmitters
  use dataTypes

  implicit none

  !   Names of these routines must be as here, as these are called by
  !    top-level inversion routines
  public                        :: dataResp, Lrows, Qrows


  !Keep the model responses as complex numbers (Z) which are required in Lrows subroutine.
  complex(kind=prec),save, private	:: Z(6)


Contains

!******************************************************************************
  subroutine dataResp(ef,Sigma,iDT,iRX,Resp,Angle,Binv,Eout,Bout)
  ! given electric field solutions (both modes--and note
  !    that the solution knows about the transmitter used),
  ! and indices into data types and receiver dictionaries for one
  ! data vector compute the complex impedance tensor.
  ! Binv is optional output argument, used needed for linearized
  ! impedance calculation in this module (not used by higher levels)

  implicit none
  type (solnVector_t), intent(in)		:: ef
  type (modelParam_t), intent(in) :: Sigma ! used to compute ef
  integer, intent(in)			:: iDT
  integer, intent(in) 			:: iRX
  real(kind=prec), intent(inout)	:: Resp(:)

  ! 2019.03.20, Liu Zhongyin, add Angle
  real(kind=prec), intent(in)    :: Angle(:)




  ! Definition of the impedance elements:
  !   iDT=Full_Impedance
  ! 			Z(1) = Zxx; Z(2) = Zxy; Z(3) = Zyx; Z(4) = Zyy
  !   iDT=Off_Diagonal_Impedance
  !  			Z(1) = Zxy, Z(2) = Zyx
  !   iDT=Full_Vertical_Components
  !  			Z(1) = Tx, Z(2) = Ty
  !   iDT=Full_Interstation_TF
  ! 		 	Z(1) = Mxx; Z(2) = Mxy; Z(3) = Myx; Z(4) = Myy
  !   iDT=Off_Diagonal_Rho_Phase
  !  			Z(1) = log(Rhoxy) , Z(2) = Phixy, Z(3) = log(Rhoyx), Z(4) = Phiyx
  !   iDT=Phase_Tensor
  !  			Z(1) = PhiXX , Z(2) = PhiXY, Z(3) = PhiYX, Z(4) = PhiYY

  !  optional argument, useful for linearized impedance
  complex(kind=prec), optional, intent(out)	:: Binv(2,2),Eout(2,2),Bout(3,2)

  !  local variables
  integer			:: iMode, i,j,xyz,ij, iComp,ncomp,iFunc,nFunc
  real(kind=prec)	:: omega,x(3),x_ref(3),detX
  complex(kind=prec)    :: tempZ(4)
  complex(kind=prec)	:: BB(3,2),EE(2,2),RR(2,2)
  complex(kind=prec)	:: det,i_omega,ctemp
  type(sparsevecC)		:: Lex,Ley,Lbx,Lby,Lbz,Lrx,Lry
  logical			:: ComputeHz,ComputeE

  !  probably should dependence on omega into BinterpSetup, as in 2D!
  omega = txDict(ef%tx)%omega

  ncomp = typeDict(iDT)%ncomp
  if(typeDict(iDT)%isComplex) then
     !  data are complex; one sensitivity calculation can be
     !   used for both real and imaginary parts
     if(mod(ncomp,2).ne.0) then
        call errStop('for complex data # of components must be even in dataResp')
     endif
     nFunc = ncomp/2
  else
     !  data are treated as real
     nFunc = ncomp
  endif
  !allocate(Z(nFunc))

 selectcase (iDT)
     case(Full_Impedance)
               x     = rxDict(iRX)%x         !Local site position (x,y,z)
		     ! First set up interpolation functionals for Ex, Ey
			  xyz = 1
			  call EinterpSetUp(ef%grid,x,xyz,Lex)
			  xyz = 2
			  call EinterpSetUp(ef%grid,x,xyz,Ley)
			 ! Then set up interpolation functionals for Bx, By
			  xyz = 1
			  call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
			  xyz = 2
			  call BfromESetUp(ef%grid,omega,x,xyz,Lby)
			  ! loop over modes
			  do iMode = 1,2
			      ! electric fields
			      EE(1,iMode) =  dotProd_noConj_scvector_f(Lex,ef%pol(iMode))
			      EE(2,iMode) =  dotProd_noConj_scvector_f(Ley,ef%pol(iMode))
			      ! magnetic fields
			      BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
			      BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
			 end do

          ! 2019.05.16, Liu Zhongyin, Add Bout
          if (present(Bout)) then
               Bout = BB
          endif

			 !invert horizontal B matrix using Kramer's rule.
			  det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
			  ctemp = BB(1,1)
			  BB(1,1) =  BB(2,2)/det
			  BB(2,2) =  ctemp/det
			  BB(1,2) = -BB(1,2)/det
			  BB(2,1) = -BB(2,1)/det

			         do j = 1,2
			           do i = 1,2
			              ij = 2*(i-1)+j
		              Z(ij) = EE(i,1)*BB(1,j)+EE(i,2)*BB(2,j)
			           enddo
			        enddo

         ! 2019.03.20 Liu Zhongyin, Add rotateZ
          call rotateZ(Z,iDT,Angle(1))  

     case(Off_Diagonal_Impedance)
               x     = rxDict(iRX)%x         !Local site position (x,y,z)
		     ! First set up interpolation functionals for Ex, Ey
			  xyz = 1
			  call EinterpSetUp(ef%grid,x,xyz,Lex)
			  xyz = 2
			  call EinterpSetUp(ef%grid,x,xyz,Ley)
			 ! Then set up interpolation functionals for Bx, By
			  xyz = 1
			  call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
			  xyz = 2
			  call BfromESetUp(ef%grid,omega,x,xyz,Lby)
			  ! loop over modes
			  do iMode = 1,2
			      ! electric fields
			      EE(1,iMode) =  dotProd_noConj_scvector_f(Lex,ef%pol(iMode))
			      EE(2,iMode) =  dotProd_noConj_scvector_f(Ley,ef%pol(iMode))
			      ! magnetic fields
			      BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
			      BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
			 end do
 
          ! 2019.05.16, Liu Zhongyin, Add Bout
          if (present(Bout)) then
            Bout = BB
   endif
 
           !invert horizontal B matrix using Kramer's rule.
            det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
            ctemp = BB(1,1)
            BB(1,1) =  BB(2,2)/det
            BB(2,2) =  ctemp/det
            BB(1,2) = -BB(1,2)/det
            BB(2,1) = -BB(2,1)/det
 
    			!Z(1) = EE(1,1)*BB(1,2)+EE(1,2)*BB(2,2)
            !Z(2) = EE(2,1)*BB(1,1)+EE(2,2)*BB(2,1)

               do j = 1,2
                  do i = 1,2
                     ij = 2*(i-1)+j
                     Z(ij) = EE(i,1)*BB(1,j)+EE(i,2)*BB(2,j)
                  enddo
               enddo
 
             ! 2019.03.20 Liu Zhongyin, Add rotateZ
             call rotateZ(Z,iDT,Angle(1))
 
      case(Full_Vertical_Components)
                x     = rxDict(iRX)%x          !Local site position (x,y,z)
               !  Vertical field TF
           ! First set up interpolation functionals for Bx, By, Bz
            xyz = 1
            call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
            xyz = 2
            call BfromESetUp(ef%grid,omega,x,xyz,Lby)
            xyz = 3
              call BfromESetUp(ef%grid,omega,x,xyz,Lbz)
            ! loop over modes
            do iMode = 1,2
                ! magnetic fields
                BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
                BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
                BB(3,iMode) = dotProd_noConj_scvector_f(Lbz,ef%pol(iMode))
           end do
          
          ! 2019.05.16, Liu Zhongyin, Add Bout
          if (present(Bout)) then
            Bout = BB
          endif

           !invert horizontal B matrix using Kramer's rule.
            det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
            ctemp = BB(1,1)
            BB(1,1) =  BB(2,2)/det
            BB(2,2) =  ctemp/det
            BB(1,2) = -BB(1,2)/det
            BB(2,1) = -BB(2,1)/det
 
 
               Z(1) = BB(3,1)*BB(1,1)+BB(3,2)*BB(2,1)
               Z(2) = BB(3,1)*BB(1,2)+BB(3,2)*BB(2,2)
 
             ! 2019.03.20 Liu Zhongyin, Add rotateZ
               call rotateZ(Z,iDT,Angle(1))
 
      case(Full_Interstation_TF)
               x     = rxDict(iRX)%x          !Local site position (x,y,z)
               x_ref = rxDict(iRX)%r          !Reference site position (x,y,z)
             ! First set up interpolation functionals for Bx, By at local site
            xyz = 1
            call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
            xyz = 2
            call BfromESetUp(ef%grid,omega,x,xyz,Lby)
            !Then set up interpolation functionals for Bx, By at the referance site
            xyz = 1
            call BfromESetUp(ef%grid,omega,x_ref,xyz,Lrx)
            xyz = 2
            call BfromESetUp(ef%grid,omega,x_ref,xyz,Lry)
              do iMode = 1,2
                ! magnetic fields at local station
                BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
                BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
                ! magnetic fields, at the REFERANCE station
                RR(1,iMode) = dotProd_noConj_scvector_f(Lrx,ef%pol(iMode))
                RR(2,iMode) = dotProd_noConj_scvector_f(Lry,ef%pol(iMode))
              end do
             
          ! 2019.05.16, Liu Zhongyin, Add Bout
             if (present(Bout)) then
               Bout(1:2,:) = RR
            endif

            ! Compute the inverse of RR using Kramer's rule
            det = RR(1,1)*RR(2,2)-RR(1,2)*RR(2,1)
            ctemp = RR(1,1)
            RR(1,1) =  RR(2,2)/det
            RR(2,2) =  ctemp/det
            RR(1,2) = -RR(1,2)/det
            RR(2,1) = -RR(2,1)/det
            ! Z = BB * RR^-1
                   do j = 1,2
                     do i = 1,2
                        ij = 2*(i-1)+j
                        Z(ij) = BB(i,1)*RR(1,j)+BB(i,2)*RR(2,j)
                     enddo
                  enddo
                  Z(1)= Z(1)-ONE
                       Z(4)= Z(4)-ONE
 
           case(Off_Diagonal_Rho_Phase)
                 x     = rxDict(iRX)%x          !Local site position (x,y,z)
            ! First set up interpolation functionals for Ex, Ey
            xyz = 1
            call EinterpSetUp(ef%grid,x,xyz,Lex)
            xyz = 2
            call EinterpSetUp(ef%grid,x,xyz,Ley)
           ! Then set up interpolation functionals for Bx, By
            xyz = 1
            call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
            xyz = 2
            call BfromESetUp(ef%grid,omega,x,xyz,Lby)
            ! loop over modes
            do iMode = 1,2
                ! electric fields
                EE(1,iMode) =  dotProd_noConj_scvector_f(Lex,ef%pol(iMode))
                EE(2,iMode) =  dotProd_noConj_scvector_f(Ley,ef%pol(iMode))
                ! magnetic fields
                BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
                BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
           end do
          
          ! 2019.05.16, Liu Zhongyin, Add Bout
          if (present(Bout)) then
            Bout = BB
         endif

           !invert horizontal B matrix using Kramer's rule.
            det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
            ctemp = BB(1,1)
            BB(1,1) =  BB(2,2)/det
            BB(2,2) =  ctemp/det
            BB(1,2) = -BB(1,2)/det
            BB(2,1) = -BB(2,1)/det
 
		       !tempZ(1) = EE(1,1)*BB(1,2)+EE(1,2)*BB(2,2)
             !tempZ(2) = EE(2,1)*BB(1,1)+EE(2,2)*BB(2,1)

             do j = 1,2
               do i = 1,2
                  ij = 2*(i-1)+j
                  tempZ(ij) = EE(i,1)*BB(1,j)+EE(i,2)*BB(2,j)
               enddo
            enddo
 
             ! 2019.03.20 Liu Zhongyin, Add rotateZ
              call rotateZ(tempZ,iDT,Angle(1))
 
              ! For Phase only, use rad, changed By LiuZhongyin 2017.05.27
             Z(1)  = log10(abs(tempZ(1))**2*MU_0/omega)
              Z(2)  = atan2(ISIGN*dimag(tempZ(1)),real(tempZ(1)))
              Z(3)  = log10(abs(tempZ(2))**2*MU_0/omega)
              Z(4)  = atan2(ISIGN*dimag(tempZ(2)),real(tempZ(2)))
 
            case(Phase_Tensor)
             ! First calculate full impedance tensor
                x     = rxDict(iRX)%x         !Local site position (x,y,z)
            ! First set up interpolation functionals for Ex, Ey
            xyz = 1
            call EinterpSetUp(ef%grid,x,xyz,Lex)
            xyz = 2
            call EinterpSetUp(ef%grid,x,xyz,Ley)
           ! Then set up interpolation functionals for Bx, By
            xyz = 1
            call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
            xyz = 2
            call BfromESetUp(ef%grid,omega,x,xyz,Lby)
            ! loop over modes
            do iMode = 1,2
                ! electric fields
                EE(1,iMode) =  dotProd_noConj_scvector_f(Lex,ef%pol(iMode))
                EE(2,iMode) =  dotProd_noConj_scvector_f(Ley,ef%pol(iMode))
                ! magnetic fields
                BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
                BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
           end do
          
          ! 2019.05.16, Liu Zhongyin, Add Bout
          if (present(Bout)) then
            Bout = BB
         endif

           !invert horizontal B matrix using Kramer's rule.
            det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
            ctemp = BB(1,1)
            BB(1,1) =  BB(2,2)/det
            BB(2,2) =  ctemp/det
            BB(1,2) = -BB(1,2)/det
            BB(2,1) = -BB(2,1)/det
 
               do j = 1,2
                  do i = 1,2
                     ij = 2*(i-1)+j
                     tempZ(ij) = EE(i,1)*BB(1,j)+EE(i,2)*BB(2,j)
                  enddo
               enddo
 
             ! 2019.03.20 Liu Zhongyin, Add rotateZ
             call rotateZ(tempZ,iDT,Angle(1))
 
             detX = dreal(tempZ(1))*dreal(tempZ(4))-dreal(tempZ(2))*dreal(tempZ(3))
 
          Z(1) = ISIGN*(dreal(tempZ(4))*dimag(tempZ(1))-dreal(tempZ(2))*dimag(tempZ(3)))/detX
          Z(2) = ISIGN*(dreal(tempZ(4))*dimag(tempZ(2))-dreal(tempZ(2))*dimag(tempZ(4)))/detX
          Z(3) = ISIGN*(dreal(tempZ(1))*dimag(tempZ(3))-dreal(tempZ(3))*dimag(tempZ(1)))/detX
          Z(4) = ISIGN*(dreal(tempZ(1))*dimag(tempZ(4))-dreal(tempZ(3))*dimag(tempZ(2)))/detX
  end select
 
   !  copy responses in Z (possibly complex) into real output vector Resp
   !  Loop over components
   iComp = 0
   do iFunc  = 1, nFunc
        if(typeDict(iDT)%isComplex) then
           iComp = iComp + 1
           Resp(iComp) = real(Z(iFunc))
           iComp = iComp + 1
           Resp(iComp) = dimag(Z(iFunc))
        else
           iComp = iComp + 1
           Resp(iComp) = real(Z(iFunc))
        endif
    enddo
 
  if(present(Binv)) then
      if(typeDict(iDT)%tfType .eq. Full_Interstation_TF) then
         Binv = RR(1:2,:)
      else
         Binv = BB(1:2,:)
      end if
  endif

  ! 2019.05.15, Liu Zhongyin, Add Einv
  if(present(Eout)) then
      Eout = EE
  endif

   ! clean up
   call deall_sparsevecc(Lex)
   call deall_sparsevecc(Ley)
   call deall_sparsevecc(Lbx)
   call deall_sparsevecc(Lby)
   call deall_sparsevecc(Lbz)
   call deall_sparsevecc(Lrx)
   call deall_sparsevecc(Lry)
   !deallocate(Z)
 
  end subroutine dataResp

!****************************************************************************
  subroutine Lrows(e0,Sigma0,iDT,iRX,L,Angle)
  !subroutine Lrows(e0,Sigma0,iDT,iRX,L)
  !  given input background electric field solution (both modes; e0),
  !  indices into data type/receiver dictionaries
  !  compute array of sparse complex vectors giving coefficients
  !  of linearized impedance functional (complex representation)

  type (solnVector_t), intent(in)		   :: e0
  type (modelParam_t), intent(in)  :: Sigma0
  integer, intent(in)			   :: iDT, iRX
  !   NOTE: Lz and Qz have to be declared as arrays for
  !     consistency with calling program (in general the
  !     number nFunc of complex data functionals that will
  !     be returned by data functional routines could be > 1)
  !   NOTE:  in principal the comparable input arguments in
  !        the 2D program should also be of type sparseVector!
  type(sparseVector_t), intent(inout)		:: L(:)

  ! 2019.05.08, Liu Zhongyin, add para Angle
  real(kind=prec), intent(in) :: Angle(:)

  !  local variables
  complex(kind=prec)	:: Binv(2,2)
  complex (kind=prec)	:: i_omega,c1,c2
  real(kind=prec)	:: Resp(12),x(3),x_ref(3),omega,detX,PT(2,2)
  type(sparsevecc)		:: L1,L2,L3,Lp11,Lp12,Lp21,Lp22
  integer			:: i,j,k,nComp,IJ(3,6),xyz,n, iComp,predictedComp
  type(sparsevecC)		:: Lex,Ley,Lbx,Lby,Lbz,Lrx,Lry
  logical			:: ComputeHz

  ! 2019.05.15, Liu Zhongyin, Add Eout,Bout
  complex(kind=prec) :: Eout(2,2),Bout(3,2),c3,c4
  real(kind=prec) :: cos2a,sin2a,cosa,sina,det,det2


  omega = txDict(e0%tx)%omega
  	 x     = rxDict(iRX)%x
     x_ref = rxDict(iRX)%r          !Reference site position (x,y,z)

  !  set up which components are needed,  ... and ! evaluate
  !   impedance, Binv for background solution
  !          ... appear to need full impedance for offdiagonal

  !   Some modifications to allow for other TFs ... e.g., the case
  !     of Hz TFs only (changes also slightly simplify generalization
  !      to interstation TFs) :  increase first dimension of array IJ
  !      to 3: (IJ(1,:) = row index in TF matrix Z;
  !             IJ(2,:) = column index in TF martrix X;
  !             IJ(3,:) = predicted field component ...
  !                     Ex = 1; Ey =2; Bz = 3; (Bx = 4; By = 5,  at referance site)
  !						(can add more cases)
  !


  select case(iDT)
     case(Full_Impedance)
        nComp = 4
        ComputeHz = .false.
        do j = 1,2
           do i = 1,2
              IJ(1,2*(i-1)+j) = i
              IJ(2,2*(i-1)+j) = j
              IJ(3,2*(i-1)+j) = i
           enddo
        enddo
        ! Call dataResp(e0,Sigma0,Full_Impedance,iRX,Resp,Binv)
        ! 2019.05.08, Liu Zhongyin, add Angle
        Call dataResp(e0,Sigma0,Full_Impedance,iRX,Resp,Angle,Binv,Eout,Bout)
     case(Off_Diagonal_Impedance)
        nComp = 2
        ComputeHz = .false.
        IJ(1,1) = 1
        IJ(2,1) = 2
        IJ(1,2) = 2
        IJ(2,2) = 1
        IJ(3,1) = 1
        IJ(3,2) = 2
        ! Call dataResp(e0,Sigma0,Full_Impedance,iRX,Resp,Binv)
        ! 2019.05.08, Liu Zhongyin, add Angle
        Call dataResp(e0,Sigma0,Full_Impedance,iRX,Resp,Angle,Binv,Eout,Bout)
      case(Full_Vertical_Components)
        nComp = 2
        ComputeHz = .true.
        IJ(1,1) = 1
        IJ(1,2) = 1
        IJ(2,1) = 1
        IJ(2,2) = 2
        IJ(3,1) = 3
        IJ(3,2) = 3
        ! Call dataResp(e0,Sigma0,Full_Vertical_Components,iRX,Resp,Binv)
        ! 2019.05.08, Liu Zhongyin, add Angle
        Call dataResp(e0,Sigma0,Full_Vertical_Components,iRX,Resp,Angle,Binv,Eout,Bout)
     case(Full_Interstation_TF)
        nComp = 4
        ComputeHz = .false.
        do j = 1,2
           do i = 1,2
              IJ(1,2*(i-1)+j) = i
              IJ(2,2*(i-1)+j) = j
              IJ(3,2*(i-1)+j) = i+3
           enddo
        enddo
        ! Call dataResp(e0,Sigma0,Full_Interstation_TF,iRX,Resp,Binv)
        ! 2019.05.08, Liu Zhongyin, add Angle
        Call dataResp(e0,Sigma0,Full_Interstation_TF,iRX,Resp,Angle,Binv,Eout,Bout)
     case(Off_Diagonal_Rho_Phase)
        ! First calculate Off_Diagonal_Impedance Ls
        ! Rho_Phase actually has 4 (real) components, but nComp here refers to the
        ! two complex off-diagonal impedance tensor components
        nComp = 2
        ComputeHz = .false.
        IJ(1,1) = 1
        IJ(2,1) = 2
        IJ(1,2) = 2
        IJ(2,2) = 1
        IJ(3,1) = 1
        IJ(3,2) = 2
        ! Call dataResp(e0,Sigma0,Full_Impedance,iRX,Resp,Binv)
        ! 2019.05.08, Liu Zhongyin, add Angle
        Call dataResp(e0,Sigma0,Full_Impedance,iRX,Resp,Angle,Binv,Eout,Bout)
     case(Phase_Tensor)
	  ! First calculate Full_Impedances Ls
        nComp = 4
        ComputeHz = .false.
        do j = 1,2
           do i = 1,2
              IJ(1,2*(i-1)+j) = i
              IJ(2,2*(i-1)+j) = j
              IJ(3,2*(i-1)+j) = i
           enddo
        enddo
        ! Call dataResp(e0,Sigma0,Full_Impedance,iRX,Resp,Binv)
        ! 2019.05.08, Liu Zhongyin, add Angle
        Call dataResp(e0,Sigma0,Full_Impedance,iRX,Resp,Angle,Binv,Eout,Bout)

     endselect

  ! Then set up interpolation functionals for Ex, Ey, Bx, By, Bz, (Bx,By at referance station)
  xyz = 1
  call EinterpSetUp(e0%grid,x,xyz,Lex)
  xyz = 2
  call EinterpSetUp(e0%grid,x,xyz,Ley)
  xyz = 1
  call BfromESetUp(e0%grid,omega,x,xyz,Lbx)
  xyz = 2
  call BfromESetUp(e0%grid,omega,x,xyz,Lby)
  xyz = 1
  call BfromESetUp(e0%grid,omega,x_ref,xyz,Lrx)
  xyz = 2
  call BfromESetUp(e0%grid,omega,x_ref,xyz,Lry)
  if(ComputeHz) then
     xyz = 3
     call BfromESetUp(e0%grid,omega,x,xyz,Lbz)
  endif

  cos2a = cos(2*Angle(1)*D2R)
  sin2a = sin(2*Angle(1)*D2R)
  cosa  = cos(Angle(1)*D2R)
  sina  = sin(Angle(1)*D2R)
  det   = Bout(1,1)*Bout(2,2) - Bout(1,2)*Bout(2,1)
  det2  = det*det
  ! 2019.05.16, Liu Zhongyin, Add new methods for rotated impedance and tipper, replace the old one
select case(iDT)
case(Full_Impedance,Phase_Tensor)
   ! zxx
   ! mode1
   ! Zxx/Ex1
   c1 =  Bout(2,2) + cos2a*Bout(2,2) - Bout(1,2)*sin2a
   c1 =  c1/2./det
   ! Zxx/Ey1
   c2 =  (cos2a-1)*Bout(1,2) + Bout(2,2)*sin2a
   c2 =  c2/2./det
   ! Zxx/Hx1
   c3 =  (Eout(1,2)*Bout(2,1) - Eout(1,1)*Bout(2,2))*((1+cos2a)*Bout(2,2) - Bout(1,2)*sin2a)
   c3 =  c3 + Eout(2,1)*Bout(2,2)*((1-cos2a)*Bout(1,2) - Bout(2,2)*sin2a)
   c3 =  c3 + Eout(2,2)*Bout(2,1)*((cos2a-1)*Bout(1,2) + Bout(2,2)*sin2a)
   c3 =  c3/2./det2
   ! Zxx/Hy1
   c4 =  (Eout(1,2)*Bout(1,1) - Eout(1,1)*Bout(1,2))*((1+cos2a)*Bout(2,2) - Bout(1,2)*sin2a)
   c4 =  c4 + Eout(2,1)*Bout(1,2)*((1-cos2a)*Bout(1,2) - Bout(2,2)*sin2a)
   c4 =  c4 + Eout(2,2)*Bout(1,1)*((cos2a-1)*Bout(1,2) + Bout(2,2)*sin2a)
   c4 = -c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(1)%L(1))
   ! mode2
   ! Zxx/Ex2
   c1 =  Bout(2,1) + cos2a*Bout(2,1) - Bout(1,1)*sin2a
   c1 = -c1/2./det
   ! Zxx/Ey2
   c2 =  Bout(1,1) - cos2a*Bout(1,1) - Bout(2,1)*sin2a
   c2 =  c2/2./det
   ! Zxx/Hx2
   c3 =  Eout(1,2)*Bout(2,1)*((1+cos2a)*Bout(2,1) - Bout(1,1)*sin2a)
   c3 =  c3 + Eout(2,2)*Bout(2,1)*((cos2a-1)*Bout(1,1) + Bout(2,1)*sin2a)
   c3 =  c3 - Bout(2,2)*((cos2a-1)*Eout(2,1)*Bout(1,1) + Eout(2,1)*Bout(2,1)*sin2a + Eout(1,1)*((1+cos2a)*Bout(2,1) - Bout(1,1)*sin2a))
   c3 = -c3/2./det2
   ! Zxx/Hy2
   c4 =  (Eout(1,2)*Bout(1,1) - Eout(1,1)*Bout(1,2))*(-(cos2a+1)*Bout(2,1) + Bout(1,1)*sin2a)
   c4 =  c4 + Eout(2,2)*Bout(1,1)*((1-cos2a)*Bout(1,1) - Bout(2,1)*sin2a)
   c4 =  c4 + Eout(2,1)*Bout(1,2)*((cos2a-1)*Bout(1,1) + Bout(2,1)*sin2a)
   c4 = -c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(1)%L(2))

   ! zxy
   ! mode1
   ! Zxy/Ex1
   c1 =  Bout(1,2) + cos2a*Bout(1,2) + Bout(2,2)*sin2a
   c1 = -c1/2./det
   ! Zxy/Ey1
   c2 =  Bout(2,2) - cos2a*Bout(2,2) + Bout(1,2)*sin2a
   c2 = -c2/2./det
   ! Zxy/Hx1
   c3 =  (Eout(2,2)*Bout(2,1) - Eout(2,1)*Bout(2,2))*((cos2a-1)*Bout(2,2) - Bout(1,2)*sin2a)
   c3 =  c3 - Eout(1,2)*Bout(2,1)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c3 =  c3 + Eout(1,1)*Bout(2,2)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c3 =  c3/2./det2
   ! Zxy/Hy1
   c4 =  (Eout(2,2)*Bout(1,1) - Eout(2,1)*Bout(1,2))*((1-cos2a)*Bout(2,2) + Bout(1,2)*sin2a)
   c4 =  c4 + Eout(1,2)*Bout(1,1)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c4 =  c4 - Eout(1,1)*Bout(1,2)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c4 =  c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(2)%L(1))
   ! mode2
   ! Zxy/Ex2
   c1 =  Bout(1,1) + cos2a*Bout(1,1) + Bout(2,1)*sin2a
   c1 =  c1/2./det
   ! Zxy/Ey2
   c2 =  Bout(2,1) - cos2a*Bout(2,1) + Bout(1,1)*sin2a
   c2 =  c2/2./det
   ! Zxy/Hx2
   c3 =  Eout(2,2)*Bout(2,1)*((1-cos2a)*Bout(2,1) + Bout(1,1)*sin2a)
   c3 =  c3 + Eout(1,2)*Bout(2,1)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a)
   c3 =  c3 - Bout(2,2)*(Eout(2,1)*((1-cos2a)*Bout(2,1) + Bout(1,1)*sin2a) + Eout(1,1)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a))
   c3 =  c3/2./det2
   ! Zxy/Hy2
   c4 =  (Eout(2,2)*Bout(1,1) - Eout(2,1)*Bout(1,2))*((1-cos2a)*Bout(2,1) + Bout(1,1)*sin2a)
   c4 =  c4 + Eout(1,2)*Bout(1,1)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a)
   c4 =  c4 - Eout(1,1)*Bout(1,2)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a)
   c4 = -c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(2)%L(2))

   ! zyx
   ! mode1
   ! Zyx/Ex1
   c1 =  Bout(1,2) - cos2a*Bout(1,2) - Bout(2,2)*sin2a
   c1 =  c1/2./det
   ! Zyx/Ey1
   c2 =  Bout(2,2) + cos2a*Bout(2,2) - Bout(1,2)*sin2a
   c2 =  c2/2./det
   ! Zyx/Hx1
   c3 =  (Eout(2,2)*Bout(2,1) - Eout(2,1)*Bout(2,2))*((1+cos2a)*Bout(2,2) - Bout(1,2)*sin2a)
   c3 =  c3 + Eout(1,2)*Bout(2,1)*((1-cos2a)*Bout(1,2) - Bout(2,2)*sin2a)
   c3 =  c3 + Eout(1,1)*Bout(2,2)*((cos2a-1)*Bout(1,2) + Bout(2,2)*sin2a)
   c3 =  c3/2./det2
   ! Zyx/Hy1
   c4 = -(Eout(2,2)*Bout(1,1) - Eout(2,1)*Bout(1,2))*((1+cos2a)*Bout(2,2) - Bout(1,2)*sin2a)
   c4 =  c4 + Eout(1,1)*Bout(1,2)*((1-cos2a)*Bout(1,2) - Bout(2,2)*sin2a)
   c4 =  c4 + Eout(1,2)*Bout(1,1)*((cos2a-1)*Bout(1,2) + Bout(2,2)*sin2a) 
   c4 =  c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(3)%L(1))
   ! mode2
   ! Zyx/Ex2
   c1 =  (cos2a-1)*Bout(1,1) + Bout(2,1)*sin2a
   c1 =  c1/2./det
   ! Zyx/Ey2
   c2 =  Bout(2,1) + cos2a*Bout(2,1) - Bout(1,1)*sin2a
   c2 = -c2/2./det
   ! Zyx/Hx2
   c3 = -Eout(2,2)*Bout(2,1)*((1+cos2a)*Bout(2,1) - Bout(1,1)*sin2a) 
   c3 =  c3 + Eout(2,1)*Bout(2,2)*((1+cos2a)*Bout(2,1) - Bout(1,1)*sin2a)
   c3 =  c3 + Eout(1,1)*Bout(2,2)*((1-cos2a)*Bout(1,1) - Bout(2,1)*sin2a)
   c3 =  c3 + Eout(1,2)*Bout(2,1)*((cos2a-1)*Bout(1,1) + Bout(2,1)*sin2a)
   c3 =  c3/2./det2
   ! Zyx/Hy2
   c4 = -(Eout(2,2)*Bout(1,1) - Eout(2,1)*Bout(1,2))*(-(1+cos2a)*Bout(2,1) + Bout(1,1)*sin2a) 
   c4 =  c4 + Eout(1,2)*Bout(1,1)*((1-cos2a)*Bout(1,1) - Bout(2,1)*sin2a)
   c4 =  c4 + Eout(1,1)*Bout(1,2)*((cos2a-1)*Bout(1,1) + Bout(2,1)*sin2a)
   c4 =  c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(3)%L(2))

   ! zyy
   ! mode1
   ! Zyy/Ex1
   c1 =  Bout(2,2) - cos2a*Bout(2,2) + Bout(1,2)*sin2a
   c1 =  c1/2./det
   ! Zyy/Ey1
   c2 =  Bout(1,2) +cos2a*Bout(1,2) + Bout(2,2)*sin2a
   c2 = -c2/2./det
   ! Zyy/Hx1
   c3 =  (Eout(1,2)*Bout(2,1) - Eout(1,1)*Bout(2,2))*((cos2a-1)*Bout(2,2) - Bout(1,2)*sin2a)
   c3 =  c3 + Eout(2,2)*Bout(2,1)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c3 =  c3 - Eout(2,1)*Bout(2,2)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c3 = -c3/2./det2
   ! Zyy/Hy1
   c4 =  (Eout(1,2)*Bout(1,1) - Eout(1,1)*Bout(1,2))*((1-cos2a)*Bout(2,2) + Bout(1,2)*sin2a)
   c4 =  c4 - Eout(2,2)*Bout(1,1)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c4 =  c4 + Eout(2,1)*Bout(1,2)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c4 = -c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(4)%L(1))
   ! mode2
   ! Zyy/Ex2
   c1 =  Bout(2,1) - cos2a*Bout(2,1) + Bout(1,1)*sin2a
   c1 = -c1/2./det
   ! Zyy/Ey2
   c2 =  Bout(1,1) + cos2a*Bout(1,1) + Bout(2,1)*sin2a
   c2 =  c2/2./det
   ! Zyy/Hx2
   c3 =  Eout(1,2)*Bout(2,1)*((cos2a-1)*Bout(2,1) - Bout(1,1)*sin2a)
   c3 =  c3 + Eout(2,2)*Bout(2,1)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a)
   c3 =  c3 - Bout(2,2)*((cos2a-1)*Eout(1,1)*Bout(2,1) - Eout(1,1)*Bout(1,1)*sin2a + Eout(2,1)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a))
   c3 =  c3/2./det2
   ! Zyy/Hy2
   c4 =  (Eout(1,2)*Bout(1,1) - Eout(1,1)*Bout(1,2))*((1-cos2a)*Bout(2,1) + Bout(1,1)*sin2a)
   c4 =  c4 - Eout(2,2)*Bout(1,1)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a)
   c4 =  c4 + Eout(2,1)*Bout(1,2)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a) 
   c4 =  c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(4)%L(2))

case(Off_Diagonal_Impedance,Off_Diagonal_Rho_Phase)
   ! zxy
   ! mode1
   ! Zxy/Ex1
   c1 =  Bout(1,2) + cos2a*Bout(1,2) + Bout(2,2)*sin2a
   c1 = -c1/2./det
   ! Zxy/Ey1
   c2 =  Bout(2,2) - cos2a*Bout(2,2) + Bout(1,2)*sin2a
   c2 = -c2/2./det
   ! Zxy/Hx1
   c3 =  (Eout(2,2)*Bout(2,1) - Eout(2,1)*Bout(2,2))*((cos2a-1)*Bout(2,2) - Bout(1,2)*sin2a)
   c3 =  c3 - Eout(1,2)*Bout(2,1)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c3 =  c3 + Eout(1,1)*Bout(2,2)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c3 =  c3/2./det2
   ! Zxy/Hy1
   c4 =  (Eout(2,2)*Bout(1,1) - Eout(2,1)*Bout(1,2))*((1-cos2a)*Bout(2,2) + Bout(1,2)*sin2a)
   c4 =  c4 + Eout(1,2)*Bout(1,1)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c4 =  c4 - Eout(1,1)*Bout(1,2)*((1+cos2a)*Bout(1,2) + Bout(2,2)*sin2a)
   c4 =  c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(1)%L(1))
   ! mode2
   ! Zxy/Ex2
   c1 =  Bout(1,1) + cos2a*Bout(1,1) + Bout(2,1)*sin2a
   c1 =  c1/2./det
   ! Zxy/Ey2
   c2 =  Bout(2,1) - cos2a*Bout(2,1) + Bout(1,1)*sin2a
   c2 =  c2/2./det
   ! Zxy/Hx2
   c3 =  Eout(2,2)*Bout(2,1)*((1-cos2a)*Bout(2,1) + Bout(1,1)*sin2a)
   c3 =  c3 + Eout(1,2)*Bout(2,1)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a)
   c3 =  c3 - Bout(2,2)*(Eout(2,1)*((1-cos2a)*Bout(2,1) + Bout(1,1)*sin2a) + Eout(1,1)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a))
   c3 =  c3/2./det2
   ! Zxy/Hy2
   c4 =  (Eout(2,2)*Bout(1,1) - Eout(2,1)*Bout(1,2))*((1-cos2a)*Bout(2,1) + Bout(1,1)*sin2a)
   c4 =  c4 + Eout(1,2)*Bout(1,1)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a)
   c4 =  c4 - Eout(1,1)*Bout(1,2)*((1+cos2a)*Bout(1,1) + Bout(2,1)*sin2a)
   c4 = -c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(1)%L(2))

   ! zyx
   ! mode1
   ! Zyx/Ex1
   c1 =  Bout(1,2) - cos2a*Bout(1,2) - Bout(2,2)*sin2a
   c1 =  c1/2./det
   ! Zyx/Ey1
   c2 =  Bout(2,2) + cos2a*Bout(2,2) - Bout(1,2)*sin2a
   c2 =  c2/2./det
   ! Zyx/Hx1
   c3 =  (Eout(2,2)*Bout(2,1) - Eout(2,1)*Bout(2,2))*((1+cos2a)*Bout(2,2) - Bout(1,2)*sin2a)
   c3 =  c3 + Eout(1,2)*Bout(2,1)*((1-cos2a)*Bout(1,2) - Bout(2,2)*sin2a)
   c3 =  c3 + Eout(1,1)*Bout(2,2)*((cos2a-1)*Bout(1,2) + Bout(2,2)*sin2a)
   c3 =  c3/2./det2
   ! Zyx/Hy1
   c4 = -(Eout(2,2)*Bout(1,1) - Eout(2,1)*Bout(1,2))*((1+cos2a)*Bout(2,2) - Bout(1,2)*sin2a)
   c4 =  c4 + Eout(1,1)*Bout(1,2)*((1-cos2a)*Bout(1,2) - Bout(2,2)*sin2a)
   c4 =  c4 + Eout(1,2)*Bout(1,1)*((cos2a-1)*Bout(1,2) + Bout(2,2)*sin2a) 
   c4 =  c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(2)%L(1))
   ! mode2
   ! Zyx/Ex2
   c1 =  (cos2a-1)*Bout(1,1) + Bout(2,1)*sin2a
   c1 =  c1/2./det
   ! Zyx/Ey2
   c2 =  Bout(2,1) + cos2a*Bout(2,1) - Bout(1,1)*sin2a
   c2 = -c2/2./det
   ! Zyx/Hx2
   c3 = -Eout(2,2)*Bout(2,1)*((1+cos2a)*Bout(2,1) - Bout(1,1)*sin2a) 
   c3 =  c3 + Eout(2,1)*Bout(2,2)*((1+cos2a)*Bout(2,1) - Bout(1,1)*sin2a)
   c3 =  c3 + Eout(1,1)*Bout(2,2)*((1-cos2a)*Bout(1,1) - Bout(2,1)*sin2a)
   c3 =  c3 + Eout(1,2)*Bout(2,1)*((cos2a-1)*Bout(1,1) + Bout(2,1)*sin2a)
   c3 =  c3/2./det2
   ! Zyx/Hy2
   c4 = -(Eout(2,2)*Bout(1,1) - Eout(2,1)*Bout(1,2))*(-(1+cos2a)*Bout(2,1) + Bout(1,1)*sin2a) 
   c4 =  c4 + Eout(1,2)*Bout(1,1)*((1-cos2a)*Bout(1,1) - Bout(2,1)*sin2a)
   c4 =  c4 + Eout(1,1)*Bout(1,2)*((cos2a-1)*Bout(1,1) + Bout(2,1)*sin2a)
   c4 =  c4/2./det2
   Call linComb_sparsevecc(Lex,c1,Ley,c2,L1)
   Call linComb_sparsevecc(Lbx,c3,Lby,c4,L2)
   Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L(2)%L(2))

case(Full_Vertical_Components)
   ! tzx
   ! mode1
   ! Tzx/Hz1
   c1 =  cosa*Bout(2,2) - Bout(1,2)*sina
   c1 =  c1/det
   ! Tzx/Hx1
   c2 =  (-Bout(2,2)*Bout(3,1) + Bout(2,1)*Bout(3,2))*(-cosa*Bout(2,2) + Bout(1,2)*sina)
   c2 =  c2/det2
   ! Tzx/Hy1
   c3 =  (Bout(1,2)*Bout(3,1) - Bout(1,1)*Bout(3,2))*(cosa*Bout(2,2) - Bout(1,2)*sina)
   c3 =  c3/det2
   Call linComb_sparsevecc(Lbz,c1,Lbx,c2,L1)
   Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(1)%L(1))
   ! mode2
   ! Tzx/Hz2
   c1 =  cosa*Bout(2,1) - Bout(1,1)*sina
   c1 = -c1/det
   ! Tzx/Hx2
   c2 =  (Bout(2,2)*Bout(3,1) - Bout(2,1)*Bout(3,2))*(cosa*Bout(2,1) - Bout(1,1)*sina)
   c2 =  c2/det2
   ! Tzx/Hy2
   c3 =  (Bout(1,2)*Bout(3,1) - Bout(1,1)*Bout(3,2))*(-cosa*Bout(2,1) + Bout(1,1)*sina)
   c3 =  c3/det2
   Call linComb_sparsevecc(Lbz,c1,Lbx,c2,L1)
   Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(1)%L(2))

   ! tzy
   ! mode1
   ! Tzy/Hz1
   c1 =  cosa*Bout(1,2) + Bout(2,2)*sina
   c1 = -c1/det
   ! Tzy/Hx1
   c2 =  (Bout(2,2)*Bout(3,1) - Bout(2,1)*Bout(3,2))*(cosa*Bout(1,2) + Bout(2,2)*sina)
   c2 =  c2/det2
   ! Tzy/Hy1
   c3 =  (Bout(1,2)*Bout(3,1) - Bout(1,1)*Bout(3,2))*(cosa*Bout(1,2) + Bout(2,2)*sina)
   c3 = -c3/det2
   Call linComb_sparsevecc(Lbz,c1,Lbx,c2,L1)
   Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(2)%L(1))
   ! mode2
   ! Tzy/Hz2
   c1 =  cosa*Bout(1,1) + Bout(2,1)*sina
   c1 =  c1/det
   ! Tzy/Hx2
   c2 =  (-Bout(2,2)*Bout(3,1) + Bout(2,1)*Bout(3,2))*(cosa*Bout(1,1) + Bout(2,1)*sina)
   c2 =  c2/det2
   ! Tzy/Hy2
   c3 =  (Bout(1,2)*Bout(3,1) - Bout(1,1)*Bout(3,2))*(cosa*Bout(1,1) + Bout(2,1)*sina)
   c3 =  c3/det2
   Call linComb_sparsevecc(Lbz,c1,Lbx,c2,L1)
   Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(2)%L(2))
  

case(Full_Interstation_TF)

endselect
   ! select case(iDT)
   ! case(Full_Impedance,Phase_Tensor)
   !    ! zxx
   !    ! mode1
   !    ! Zxx/Ex1
   !    c1 = Binv(1,1)*(1+cos2a)/2. + Binv(1,2)*sin2a/2.
   !    ! Zxx/Hx1
   !    c2 = -z(1) + (1-cos2a)/2.*Eout(2,2)/Bout(2,2) + sin2a/2.*Eout(1,2)/Bout(2,2)
   !    c2 = Binv(1,1)*c2
   !    ! Zxx/Hy1
   !    c3 = -z(1) + (1+cos2a)/2.*Eout(1,2)/Bout(1,2) + sin2a/2.*Eout(2,2)/Bout(1,2)
   !    c3 = Binv(1,2)*c3
   !    Call linComb_sparsevecc(Lex,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(1)%L(1))
   !    ! mode2
   !    ! Zxx/Ex2
   !    c1 = Binv(2,1)*(1+cos2a)/2. + Binv(2,2)*sin2a/2.
   !    ! Zxx/Hx2
   !    c2 = -z(1) + (1-cos2a)/2.*Eout(2,1)/Bout(2,1) + sin2a/2.*Eout(1,1)/Bout(2,1)
   !    c2 = Binv(2,1)*c2
   !    ! Zxx/Hy2
   !    c3 = -z(1) + (1+cos2a)/2.*Eout(1,1)/Bout(1,1) + sin2a/2.*Eout(2,1)/Bout(1,1)
   !    c3 = Binv(2,2)*c3
   !    Call linComb_sparsevecc(Lex,c1,Lbx,c2,L1)
   !    call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(1)%L(2))

   !    ! zxy
   !    ! mode1
   !    ! Zxy/Ex1
   !    c1 = Binv(1,2)*(1+cos2a)/2. - Binv(1,1)*sin2a/2.
   !    ! Zxy/Hx1
   !    c2 = -z(2) + (1+cos2a)/2.*Eout(1,2)/Bout(2,2) + sin2a/2.*Eout(2,2)/Bout(2,2)
   !    c2 = Binv(1,1)*c2
   !    ! Zxy/Hy1
   !    c3 = -z(2) + (cos2a-1)/2.*Eout(2,2)/Bout(1,2) - sin2a/2.*Eout(1,2)/Bout(1,2)
   !    c3 = Binv(1,2)*c3
   !    Call linComb_sparsevecc(Lex,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(2)%L(1))
   !    ! mode2
   !    ! Zxy/Ex2
   !    c1 = Binv(2,2)*(1+cos2a)/2. - Binv(2,1)*sin2a/2.
   !    ! Zxy/Hx2
   !    c2 = -z(2) +(1+cos2a)/2.*Eout(1,1)/Bout(2,1) + sin2a/2.*Eout(2,1)/Bout(2,1)
   !    c2 = Binv(2,1)*c2
   !    ! Zxy/Hy2
   !    c3 = -z(2) +(cos2a-1)/2.*Eout(2,1)/Bout(1,1) - sin2a/2.*Eout(1,1)/Bout(1,1)
   !    c3 = Binv(2,2)*c3
   !    Call linComb_sparsevecc(Lex,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(2)%L(2))

   !    ! zyx
   !    ! mode1
   !    ! Zyx/Ey1
   !    c1 = Binv(1,1)*(1+cos2a)/2. + Binv(1,2)*sin2a/2.
   !    ! Zyx/Hx1
   !    c2 = -z(3) + (cos2a-1)/2.*Eout(1,2)/Bout(2,2) + sin2a/2.*Eout(2,2)/Bout(2,2)
   !    c2 = Binv(1,1)*c2
   !    ! Zyx/Hy1
   !    c3 = -z(3) + (1+cos2a)/2.*Eout(2,2)/Bout(1,2) - sin2a/2.*Eout(1,2)/Bout(1,2)
   !    c3 = Binv(1,2)*c3
   !    Call linComb_sparsevecc(Ley,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(3)%L(1))
   !    ! mode2
   !    ! Zyx/Ey2
   !    c1 = Binv(2,1)*(1+cos2a)/2. + Binv(2,2)*sin2a/2.
   !    ! Zyx/Hx2
   !    c2 = -z(3) + (cos2a-1)/2.*Eout(1,1)/Bout(2,1) + sin2a/2.*Eout(2,1)/Bout(2,1)
   !    c2 = Binv(2,1)*c2
   !    ! Zyx/Hy2
   !    c3 = -z(3) + (1+cos2a)/2.*Eout(2,1)/Bout(1,1) - sin2a/2.*Eout(1,1)/Bout(1,1)
   !    c3 = Binv(2,2)*c3
   !    Call linComb_sparsevecc(Ley,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(3)%L(2))

   !    ! zyy
   !    ! mode1
   !    ! Zyy/Ey1
   !    c1 = Binv(1,2)*(1+cos2a)/2. - Binv(1,1)*sin2a/2.
   !    ! Zyy/Hx1
   !    c2 = -z(4) + (1+cos2a)/2.*Eout(2,2)/Bout(2,2) - sin2a/2.*Eout(1,2)/Bout(2,2)
   !    c2 = Binv(1,1)*c2
   !    ! Zyy/Hy1
   !    c3 = -z(4) + (1-cos2a)/2.*Eout(1,2)/Bout(1,2) - sin2a/2.*Eout(2,2)/Bout(1,2)
   !    c3 = Binv(1,2)*c3
   !    Call linComb_sparsevecc(Ley,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(4)%L(1))
   !    ! mode2
   !    ! Zyy/Ey2
   !    c1 = Binv(2,2)*(1+cos2a)/2. - Binv(2,1)*sin2a/2.
   !    ! Zyy/Hx2
   !    c2 = -z(4) + (1+cos2a)/2.*Eout(2,1)/Bout(2,1) - sin2a/2.*Eout(1,1)/Bout(2,1)
   !    c2 = Binv(2,1)*c2
   !    ! Zyy/Hyx
   !    c3 = -z(4) + (1-cos2a)/2.*Eout(1,1)/Bout(1,1) -sin2a/2.*Eout(2,1)/Bout(1,1)
   !    c3 = Binv(2,2)*c3
   !    Call linComb_sparsevecc(Ley,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(4)%L(2))

   ! case(Off_Diagonal_Impedance,Off_Diagonal_Rho_Phase)
   !    ! zxy
   !    ! mode1
   !    ! Zxy/Ex1
   !    c1 = Binv(1,2)*(1+cos2a)/2. - Binv(1,1)*sin2a/2.
   !    ! Zxy/Hx1
   !    c2 = -z(2) + (1+cos2a)/2.*Eout(1,2)/Bout(2,2) + sin2a/2.*Eout(2,2)/Bout(2,2)
   !    c2 = Binv(1,1)*c2
   !    ! Zxy/Hy1
   !    c3 = -z(2) + (cos2a-1)/2.*Eout(2,2)/Bout(1,2) - sin2a/2.*Eout(1,2)/Bout(1,2)
   !    c3 = Binv(1,2)*c3
   !    Call linComb_sparsevecc(Lex,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(1)%L(1))
   !    ! mode2
   !    ! Zxy/Ex2
   !    c1 = Binv(2,2)*(1+cos2a)/2. - Binv(2,1)*sin2a/2.
   !    ! Zxy/Hx2
   !    c2 = -z(2) +(1+cos2a)/2.*Eout(1,1)/Bout(2,1) + sin2a/2.*Eout(2,1)/Bout(2,1)
   !    c2 = Binv(2,1)*c2
   !    ! Zxy/Hy2
   !    c3 = -z(2) +(cos2a-1)/2.*Eout(2,1)/Bout(1,1) - sin2a/2.*Eout(1,1)/Bout(1,1)
   !    c3 = Binv(2,2)*c3
   !    Call linComb_sparsevecc(Lex,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(1)%L(2))

   !    ! zyx
   !    ! mode1
   !    ! Zyx/Ey1
   !    c1 = Binv(1,1)*(1+cos2a)/2. + Binv(1,2)*sin2a/2.
   !    ! Zyx/Hx1
   !    c2 = -z(3) + (cos2a-1)/2.*Eout(1,2)/Bout(2,2) + sin2a/2.*Eout(2,2)/Bout(2,2)
   !    c2 = Binv(1,1)*c2
   !    ! Zyx/Hy1
   !    c3 = -z(3) + (1+cos2a)/2.*Eout(2,2)/Bout(1,2) - sin2a/2.*Eout(1,2)/Bout(1,2)
   !    c3 = Binv(1,2)*c3
   !    Call linComb_sparsevecc(Ley,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(2)%L(1))
   !    ! mode2
   !    ! Zyx/Ey2
   !    c1 = Binv(2,1)*(1+cos2a)/2. + Binv(2,2)*sin2a/2.
   !    ! Zyx/Hx2
   !    c2 = -z(3) + (cos2a-1)/2.*Eout(1,1)/Bout(2,1) + sin2a/2.*Eout(2,1)/Bout(2,1)
   !    c2 = Binv(2,1)*c2
   !    ! Zyx/Hy2
   !    c3 = -z(3) + (1+cos2a)/2.*Eout(2,1)/Bout(1,1) - sin2a/2.*Eout(1,1)/Bout(1,1)
   !    c3 = Binv(2,2)*c3
   !    Call linComb_sparsevecc(Ley,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(2)%L(2))

   ! case(Full_Vertical_Components)
   !    ! tzx
   !    ! mode1
   !    ! Tzx/Hz1
   !    c1 = Binv(1,1)*cosa + Binv(1,2)*sina
   !    ! Tzx/Hx1
   !    c2 = -z(1) + sina*Bout(3,2)/Bout(2,2)
   !    c2 = Binv(1,1)*c2
   !    ! Tzx/Hy1
   !    c3 = -z(1) + cosa*Bout(3,2)/Bout(1,2)
   !    c3 = Binv(1,2)*c3
   !    Call linComb_sparsevecc(Lbz,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(1)%L(1))
   !    ! mode2
   !    ! Tzx/Hz2
   !    c1 = Binv(2,1)*cosa + Binv(2,2)*sina
   !    ! Tzx/Hx2
   !    c2 = -z(1) + sina*Bout(3,1)/Bout(2,1)
   !    c2 = Binv(2,1)*c2
   !    ! Tzx/Hy2
   !    c3 = -z(1) + cosa*Bout(3,1)/Bout(1,1)
   !    c3 = Binv(2,2)*c3
   !    Call linComb_sparsevecc(Lbz,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(1)%L(2))

   !    ! tzy
   !    ! mode1
   !    ! Tzy/Hz1
   !    c1 = Binv(1,2)*cosa - Binv(1,1)*sina
   !    ! Tzy/Hx1
   !    c2 = -z(2) + cosa*Bout(3,2)/Bout(2,2)
   !    c2 = Binv(1,1)*c2
   !    ! Tzy/Hy1
   !    c3 = -z(2) -sina*Bout(3,2)/Bout(1,2)
   !    c3 = Binv(1,2)*c3
   !    Call linComb_sparsevecc(Lbz,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(2)%L(1))
   !    ! mode2
   !    ! Tzy/Hz2
   !    c1 = Binv(2,2)*cosa - Binv(2,1)*sina
   !    ! Tzy/Hx2
   !    c2 = -z(2) + cosa*Bout(3,1)/Bout(2,1)
   !    c2 = Binv(2,1)*c2
   !    ! Tzy/Hy2
   !    c3 = -z(2) - sina*Bout(3,1)/Bout(1,1)
   !    c3 = Binv(2,2)*c3
   !    Call linComb_sparsevecc(Lbz,c1,Lbx,c2,L1)
   !    Call linComb_sparsevecc(L1,C_ONE,Lby,c3,L(2)%L(2))
     

   ! case(Full_Interstation_TF)

   ! endselect





!   !  compute sparse vector representations of linearized functionals
!   do n = 1,nComp
!      !  i runs over rows of TF matrix, j runs over columns of TF
!      i = IJ(1,n)
!      j = IJ(2,n)
!      predictedComp = IJ(3,n)
!      c1 = Z(2*(i-1)+1)
!      c2 = Z(2*(i-1)+2)
!     if(typeDict(iDT)%tfType .eq. Full_Interstation_TF) then
!       Call linComb_sparsevecc(Lrx,c1,Lry,c2,L1)
!     else
!       Call linComb_sparsevecc(Lbx,c1,Lby,c2,L1)
!     end if
!      do k = 1,2
!         !  k defines which mode the linearized functional is
!         !   to be applied to
!         c1 = Binv(k,j)  !In case of interstaion TF, Binv = RRinv.
!         c2 = -c1
!         if(predictedComp.eq.1) then
!            !  component in x row of impedance tensor
!            Call linComb_sparsevecc(Lex,c1,L1,c2,L(n)%L(k))
!         elseif(predictedComp.eq.2) then
!            !  component in y row of impedance tensor
!            Call linComb_sparsevecc(Ley,c1,L1,c2,L(n)%L(k))
!         elseif(predictedComp.eq.3) then
!            !  component in Bz row (vertical field TF)
!            Call linComb_sparsevecc(Lbz,c1,L1,c2,L(n)%L(k))
!         elseif(predictedComp.eq.4) then
!            !  component in x row (interstation TF)
!            Call linComb_sparsevecc(Lbx,c1,L1,c2,L(n)%L(k))
!         elseif(predictedComp.eq.5) then
!            !  component in y row (interstation TF)
!            Call linComb_sparsevecc(Lby,c1,L1,c2,L(n)%L(k))
!         endif
!      enddo
!   enddo
if (typeDict(iDT)%tfType .eq. Off_Diagonal_Rho_Phase) then
       do k=1,2 ! 2 modes
        ! PHSYX
        c1 =dcmplx(0.0d0,1.0d0)*conjg(Z(3)) / (abs(Z(3))**TWO)
	     Call linComb_sparsevecc(L(2)%L(k),c1,L(2)%L(k),C_ZERO,L(4)%L(k))

		 !log RHOYX
	    ! c1 =  TWO*conjg(Z(3))/(abs(Z(3))**TWO)*dlog(10.0d0)
		 ! devided by Ln10, Liuzhongyin 2017.06.04
	     c1 =  TWO*conjg(Z(3))/(abs(Z(3))**TWO)/dlog(10.0d0)
        Call linComb_sparsevecc(L(2)%L(k),c1,L(2)%L(k),C_ZERO,L(3)%L(k))

        ! PHSXY
        c1 =dcmplx(0.0d0,1.0d0)*conjg(Z(2))/(abs(Z(2))**TWO)
		Call linComb_sparsevecc(L(1)%L(k),c1,L(1)%L(k),C_ZERO,L(2)%L(k))

          !log(RHOXY)
          ! c1 =  TWO*conjg(Z(2))  /(abs(Z(2))**TWO)*dlog(10.0d0)
          ! devided by Ln10, Liuzhongyin 2017.06.04
          c1 =  TWO*conjg(Z(2))  /(abs(Z(2))**TWO)/dlog(10.0d0)
	     Call linComb_sparsevecc(L(1)%L(k),c1,L(1)%L(k),C_ZERO,L1)
        L(1)%L(k) = L1

     enddo
  end if
   if (typeDict(iDT)%tfType .eq. Phase_Tensor) then
       do k=1,2 ! 2 modes
	    !calculate Phase Tensor Elements
			detX = dreal(Z(1))*dreal(Z(4))-dreal(Z(2))*dreal(Z(3))

		PT(1,1) = ISIGN*(dreal(Z(4))*dimag(Z(1))-dreal(Z(2))*dimag(Z(3)))/detX
		PT(1,2) = ISIGN*(dreal(Z(4))*dimag(Z(2))-dreal(Z(2))*dimag(Z(4)))/detX
		PT(2,1) = ISIGN*(dreal(Z(1))*dimag(Z(3))-dreal(Z(3))*dimag(Z(1)))/detX
		PT(2,2) = ISIGN*(dreal(Z(1))*dimag(Z(4))-dreal(Z(3))*dimag(Z(2)))/detX

		 !PTXX
		 !dx11
	     c1 =  dcmplx(MinusONE*PT(1,1) * dreal(Z(4)) / detX, R_ZERO)
		 !dx12
		  c2 =  dcmplx((PT(1,1) * dreal(Z(3)) - ISIGN*dimag(Z(3))) / detX ,R_ZERO)
         Call linComb_sparsevecc(L(1)%L(k),c1,L(2)%L(k),c2,L1)
		 !dx21
	     c1 =  dcmplx(PT(1,1) * dreal(Z(2)) / detX , R_ZERO)
		 !dx22
		  c2 =  dcmplx((MinusONE * PT(1,1) * dreal(Z(1)) + ISIGN*dimag(Z(1)))/ detX,R_ZERO)
         Call linComb_sparsevecc(L(3)%L(k),c1,L(4)%L(k),c2,L2)
         Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L3)
		 !dy11
	     c1 =  dcmplx(R_ZERO,dreal(Z(4)) / detX)
		 !dy21
		 c2 = dcmplx(R_ZERO,MinusONE* dreal(Z(2)) / detX)
		 Call linComb_sparsevecc(L(1)%L(k),c1,L(3)%L(k),c2,L1)
		 Call linComb_sparsevecc(L3,C_ONE,L1,C_ONE,Lp11)

		 !PTXY
		 !dx11
	     c1 =  dcmplx(MinusONE*PT(1,2) * dreal(Z(4)) / detX, R_ZERO)
		 !dx12
		  c2 =  dcmplx((PT(1,2) * dreal(Z(3)) - ISIGN*dimag(Z(4))) / detX, R_ZERO)
         Call linComb_sparsevecc(L(1)%L(k),c1,L(2)%L(k),c2,L1)
		 !dx21
	     c1 =  dcmplx(PT(1,2) * dreal(Z(2)) / detX, R_ZERO)
		 !dx22
		  c2 =  dcmplx((MinusONE*PT(1,2) * dreal(Z(1)) + ISIGN*dimag(Z(2)))/ detX, R_ZERO)
         Call linComb_sparsevecc(L(3)%L(k),c1,L(4)%L(k),c2,L2)
         Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L3)
		 
		 !dy12
	     c1 =  dcmplx(R_ZERO, dreal(Z(4)) / detX)
		 !dy22
		 c2 = dcmplx(R_ZERO, MinusONE* dreal(Z(2))/ detX)
         Call linComb_sparsevecc(L(2)%L(k),c1,L(4)%L(k),c2,L1)
		 Call linComb_sparsevecc(L3,C_ONE,L1,C_ONE,Lp12)

		 !PTYX
		 !dx11
	     c1 =  dcmplx((MinusONE*PT(2,1) * dreal(Z(4)) + ISIGN*dimag(Z(3)))/ detX, R_ZERO)
		 !dx12
		  c2 =  dcmplx(PT(2,1) * dreal(Z(3)) / detX, R_ZERO)
         Call linComb_sparsevecc(L(1)%L(k),c1,L(2)%L(k),c2,L1)
		 !dx21
	     c1 = dcmplx(( PT(2,1) * dreal(Z(2)) - ISIGN*dimag(Z(1)))/ detX, R_ZERO)
		 !dx22
		  c2 =  dcmplx(MinusONE*PT(2,1) * dreal(Z(1)) / detX, R_ZERO)
         Call linComb_sparsevecc(L(3)%L(k),c1,L(4)%L(k),c2,L2)
         Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L3)
		 
		 !dy11
	     c1 =  dcmplx(R_ZERO, MinusONE*dreal(Z(3)) / detX)
		 !dy21
		 c2 = dcmplx(R_ZERO,  dreal(Z(1)) / detX)
         Call linComb_sparsevecc(L(1)%L(k),c1,L(3)%L(k),c2,L1)
		 Call linComb_sparsevecc(L3,C_ONE,L1,C_ONE,Lp21)


		 !PTYY
		 !dx11
	     c1 =  dcmplx((MinusONE*PT(2,2) * dreal(Z(4)) + ISIGN*dimag(Z(4)))/ detX, R_ZERO)
		 !dx12
		  c2 =  dcmplx(PT(2,2) * dreal(Z(3)) / detX, R_ZERO)
         Call linComb_sparsevecc(L(1)%L(k),c1,L(2)%L(k),c2,L1)
		 !dx21
	     c1 = dcmplx(( PT(2,2) * dreal(Z(2)) - ISIGN*dimag(Z(2)))/ detX, R_ZERO)
		 !dx22
		  c2 =  dcmplx(MinusONE*PT(2,2) * dreal(Z(1)) / detX, R_ZERO)
         Call linComb_sparsevecc(L(3)%L(k),c1,L(4)%L(k),c2,L2)
         Call linComb_sparsevecc(L1,C_ONE,L2,C_ONE,L3)
		 
		 !dy12
	     c1 =  dcmplx(R_ZERO, MinusONE*dreal(Z(3)) / detX)
		 !dy22
		 c2 = dcmplx(R_ZERO,dreal(Z(1)) / detX)
         Call linComb_sparsevecc(L(2)%L(k),c1,L(4)%L(k),c2,L1)
		 Call linComb_sparsevecc(L3,C_ONE,L1,C_ONE,Lp22)


    	 !Finally overwrite Impedance Ls of this mode with Phase Tensor Ls
		 Call linComb_sparsevecc(Lp11,C_ONE,Lp11,C_ZERO,L(1)%L(k))
 		 Call linComb_sparsevecc(Lp12,C_ONE,Lp12,C_ZERO,L(2)%L(k))
 		 Call linComb_sparsevecc(Lp21,C_ONE,Lp21,C_ZERO,L(3)%L(k))
 		 Call linComb_sparsevecc(Lp22,C_ONE,Lp22,C_ZERO,L(4)%L(k))
    enddo
  end if


 ! clean up
 if (typeDict(iDT)%tfType .eq. Phase_Tensor) then
   call deall_sparsevecc(L2)
  call deall_sparsevecc(L3)
  call deall_sparsevecc(Lp11)
  call deall_sparsevecc(Lp12)
  call deall_sparsevecc(Lp21)
  call deall_sparsevecc(Lp22)
  end if

   call deall_sparsevecc(L1)
  call deall_sparsevecc(Lex)
  call deall_sparsevecc(Ley)
  call deall_sparsevecc(Lbx)
  call deall_sparsevecc(Lby)
  call deall_sparsevecc(Lbz)
  call deall_sparsevecc(Lrx)
  call deall_sparsevecc(Lry)

  end subroutine Lrows
!
!****************************************************************************
  subroutine Qrows(e0,Sigma0,iDT,iRX,zeroValued,Qreal,Qimag)
  !  given input background solution vector (e0) and model parameter (Sigma0)
  !  and indices into data type and receiver dictionaries
  !  compute derivative of data functional with respect to model parameters
  !  for all components of the data type ...
  !             (ZERO VECTORS FOR 3D MT!!!!)

  type (solnVector_t), intent(in)		    :: e0
  type (modelParam_t), intent(in)       :: Sigma0
  integer, intent(in)			              :: iDT, iRX
  logical, intent(out)                      :: zeroValued
  !   NOTE: Qreal and Qimag have to be declared as arrays for
  !     consistency with calling program (in general the
  !     number nFunc of complex data functionals that will
  !     be returned by data functional routines could be > 1)
  !   NOTE: Qreal and Qimag both exist regardless of whether the data
  !     are real or complex, since Q itself is complex
  type(modelParam_t), intent(inout)     :: Qreal(:), Qimag(:)

  ! local variables
  integer       :: istat,ncomp,nFunc,iFunc
  logical       :: isComplex

  ncomp = typeDict(iDT)%nComp
  isComplex = typeDict(iDT)%isComplex
  if(isComplex) then
     !  data are complex; one sensitivity calculation can be
     !   used for both real and imaginary parts
     if(mod(ncomp,2).ne.0) then
        call errStop('for complex data # of components must be even in Qrows')
     endif
     nFunc = ncomp/2
  else
     !  data are treated as real: full sensitivity computation is required
     !   for each component
     nFunc = ncomp
  endif

  ! set the rows of Q to zero
  if((size(Qreal) .ne. nFunc) .or. (size(Qimag) .ne. nFunc)) then
    call errStop('incorrect output size in Qrows')
  endif

  ! for efficiency, for zero vectors, just set the logical to true and exit
  zeroValued = .true.

  !do iFunc = 1, nFunc
  !  Qreal(iFunc) = Sigma0
  !  call zero(Qreal(iFunc))
  !  Qimag(iFunc) = Sigma0
  !  call zero(Qimag(iFunc))
  !enddo

  end subroutine Qrows

  ! **********************************************************************
  ! 2019.03.20 Liu Zhongyin, Add rotate for Z
  ! Note: For the Off_Diagonal_Rho_Phase, this subroutine should be 
  ! called before the calculation of Rho and Phase
  subroutine rotateZ(d,iDT,angle)
   implicit none
   complex(kind=prec), intent(inout)   :: d(4)
   integer, intent(in)                 :: iDT
   real(kind=prec), intent(in)         :: angle

   ! local varialbles
   complex(kind=prec)        :: d1, d2, d3, d4
   real(kind=prec)           :: cos2a, sin2a, sina, cosa

  selectcase (iDT)
      case(Full_Impedance, Phase_Tensor)
         cos2a = cos(2*angle*D2R)
         sin2a = sin(2*angle*D2R)

         d1 = 1/2.*((1+cos2a)*d(1)+sin2a*d(2)+sin2a*d(3)+(1-cos2a)*d(4))
         d2 = 1/2.*(-sin2a*d(1)+(1+cos2a)*d(2)-(1-cos2a)*d(3)+sin2a*d(4))
         d3 = 1/2.*(-sin2a*d(1)-(1-cos2a)*d(2)+(1+cos2a)*d(3)+sin2a*d(4))
         d4 = 1/2.*((1-cos2a)*d(1)-sin2a*d(2)-sin2a*d(3)+(1+cos2a)*d(4))

         d(1) = d1
         d(2) = d2
         d(3) = d3
         d(4) = d4

      case(Off_Diagonal_Impedance, Off_Diagonal_Rho_Phase)
         cos2a = cos(2*angle*D2R)
         sin2a = sin(2*angle*D2R)

         d1 = 1/2.*((1+cos2a)*d(1)+sin2a*d(2)+sin2a*d(3)+(1-cos2a)*d(4))
         d2 = 1/2.*(-sin2a*d(1)+(1+cos2a)*d(2)-(1-cos2a)*d(3)+sin2a*d(4))
         d3 = 1/2.*(-sin2a*d(1)-(1-cos2a)*d(2)+(1+cos2a)*d(3)+sin2a*d(4))
         d4 = 1/2.*((1-cos2a)*d(1)-sin2a*d(2)-sin2a*d(3)+(1+cos2a)*d(4))

         d(1) = d2
         d(2) = d3
         d(3) = cmplx(0,0)
         d(4) = cmplx(0,0)

      case(Full_Vertical_Components)
         sina = sin(angle*D2R)
         cosa = cos(angle*D2R)

         d1 =  d(1)*cosa+d(2)*sina
         d2 = -d(1)*sina+d(2)*cosa

         d(1) = d1
         d(2) = d2

      case(Full_Interstation_TF)

  end select

  end subroutine rotateZ

end module dataFunc
