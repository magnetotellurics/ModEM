
! *****************************************************************************
! model_data is used for sharing the data for the joint forward modeling-
! inversion scheme. Data module where the current model definition (grid,
! conductivity, frequency) is stored. This module is sued by SetUp routines
! (for initialization, modification of any parameter values), and PDE coefficient
! initialization routines.
module modelOperator3D
  !  This module merges old modules model_data, model_data_update, multA,
  !      preconditioner, divcorr
  !   ... into a single module that does all operations involving model
  !     operators including the full EM operator and operators needed for
  !     divergence correction and preconditioning.  Everything that
  !     initializes or uses the equation coefficients are now in this module
  !     allowing arrays of equation coefficients, etc. to be private
  !       and essentially global to this module

  use math_constants
  use utilities
  use gridcalc             ! staggered grid definitions
  use sg_vector            ! generic routines for vector operations on the
  use sg_boundary
  use ModelSpace
  use boundary_Pek          ! sets the boundary conditions
  implicit none

  ! * These variables are used by model equation
  ! * and preconditioner modules;
  ! * All variables are saved until changed/deallocated

  save

  !!!!!!!>>>>>>>>> Tesnor resistivity and conductivity in cells(private)
  type(rscalar), private    :: rho_C(6),sigma_C(6) ! cells - needed for boundary conditions

  !!!!!!!>>>>>>>>> FROM model_data
  type(grid_t), private, target 	::	mGrid   ! THE model grid
  !type(rvector), public			::	volE    ! THE volume elements
  !type(rvector), private		::	condE   ! THE edge conductivities
  real(kind=prec),private	::      omega   ! THE (active) frequency

  ! Mackie's A Matrix
  complex (kind=prec), pointer, dimension(:), private    :: spa
  integer, pointer, dimension(:), private :: ija
  ! Mackie's ILU Operator
  integer, pointer, dimension(:), private :: jlu,ju
  complex (kind=prec), pointer, dimension(:), private :: alu
  integer, pointer, dimension(:), private :: ijapot
  real (kind=prec), pointer, dimension(:), private :: sppot,sppot1


Contains

  subroutine ModelDataInit(inGrid)
  !**********************************************************************
  ! *   Copies grid to mGrid
  !      and/or compute variables stored in model_data module:
  !         create: allocate for edge conductivity
  !              volume weights;
  !         EdgeVolume:  compute volume weights for edge-centered prisms
  !
  !**********************************************************************

    implicit none
    !  INPUTS:
    type (grid_t), intent(in)		  :: inGrid

    !   copy inGrid to mGrid
    call copy_grid(mGrid,inGrid) !感觉就需要这个复制grid的，后面创建的体积好像用不到

    ! Allocate data structure for volume elements, and compute these
    Call create(mGrid, V_E, EDGE)

    ! Use the grid (which, potentially, maybe have been updated!) to set up
    !   all the grid length, surface and volume elements stored in GridCalc.
    ! Want to initialize them here in case the grid gets updated along the way.
    ! The reason for storing them in GridCalc is that they are also used
    !   by ModelMap, EMfieldInterp, nestedEM
    Call EdgeVolume(mGrid, V_E)
    Call NodeVolume(mGrid, V_N) ! used for divergence correction

    ! set a default omega
    omega = 0.0

  end subroutine ModelDataInit


  subroutine ModelDataCleanUp
  
    integer :: ia
    ! Deallocated the grid
    call deall_grid(mGrid)

    ! and the grid elements stored in GridCalc
    call deall_rvector(V_E)
    call deall_rscalar(V_N)

    do ia = 1,6
      call deall_rscalar(rho_C(ia))
      call deall_rscalar(sigma_C(ia))
    enddo
    
  end subroutine ModelDataCleanUp


  ! ***************************************************************************
  ! * UpdateFreqCond updates the frequency that is currently being use and
  !*  conductivity values on the edges
  subroutine UpdateFreqRho(inOmega, RhoParam)

    implicit none
    real(kind=prec)                 :: inOmega
    type(modelParam_t), intent(in)            :: RhoParam      ! input conductivity
    ! structure on the center of the grid

    omega = inOmega

    ! TEMPORARY; REQUIRED FOR BOUNDARY CONDITIONS
    !  set static array for cell conductivities
    call ModelParamToTensor(RhoParam,sigma_C,rho_C)

  end subroutine UpdateFreqRho  ! UpdateFreqCond


  ! ***************************************************************************
  subroutine UpdateRho(RhoParam)

    implicit none
    type(modelParam_t), intent(in)            :: RhoParam      ! input conductivity

    ! TEMPORARY; REQUIRED FOR BOUNDARY CONDITIONS
    !  set static array for cell conductivities
    call ModelParamToTensor(RhoParam,sigma_C,rho_C)

  end subroutine UpdateRho 
 
  
!**********************************************************************
! Sets boundary conditions for HH scheme.
  Subroutine SetBound(imode,period,E0,BC,iTx)
  !call SetBound(e0%Pol_index(iMode),period,e0%pol(imode),b0%bc,iTx)
  !type(RHS_t), save, private		:: b0
    implicit none
    !  Input mode, period
    integer, intent(in)		:: imode
    integer, intent(in)		:: iTx
    real(kind=prec)	:: period
    ! Output electric field first guess (for iterative solver)
    type(cvector), intent(inout)	:: E0
    ! Output boundary conditions
    type(cboundary), intent(inout)	:: BC

    ! top boundary
    if(imode.eq.1) then
      BC%xZMin=C_ZERO
      BC%yZMin=C_ONE	
		elseif(imode.eq.2) then
      BC%xZMin=C_ONE
      BC%yZMin=C_ZERO
		else
			call errStop('Polarization mode does not exist!')	
		end if
	 
    ! side boundaries
    call BC_x0_Pek(imode,period,mGrid,sigma_C,E0,BC)
    call BC_y0_Pek(imode,period,mGrid,sigma_C,E0,BC)

		call setBC(BC, E0) ! E0 (=) BC 

  end subroutine SetBound
  

  ! **************************************
  subroutine calcb(inE, adjt, outE)
    implicit none
    type (cvector), intent (in)              :: inE
    logical, intent (in)                     :: adjt
    type (cvector), intent (inout)           :: outE
	 ! local variables
    integer                     :: status     ! for dynamic memory allocation
    integer                     :: ia,i,j,k,l,m,n,nzn ! dummy variables
    integer                     :: nhx,nhy,nhz,np,np5
	  integer                     :: ii,jj,kk
!   for Hx>>
    real(kind=prec) :: ribcyx,ribcyz,ribkyx,ribkyz,rijcyx,rijcyz,rijkyx,rijkyz !total 8
	  real(kind=prec) :: ribczx,ribczy,ribkzx,ribkzy,rijczx,rijczy,rijkzx,rijkzy 
	  real(kind=prec) :: rijcy,rijky
	  real(kind=prec) :: ribkz,rijkz
!   for Hy>>
	  real(kind=prec) :: rajcxy,rajcxz,rajkxy,rajkxz,rijcxy,rijcxz,rijkxy,rijkxz
	  real(kind=prec) :: rajczx,rajczy,rajkzx,rajkzy !,rijczx,rijczy,rijkzx,rijkzy
	  real(kind=prec) :: rijcx,rijkx
	  real(kind=prec) :: rajkz !,rijkz
!   for Hz>>
	  real(kind=prec) :: rabkxy,rabkxz,ribkxy,ribkxz !,rajkxy,rajkxz,rijkxy,rijkxz
    real(kind=prec) :: rabkyx,rabkyz,rajkyx,rajkyz !,ribkyx,ribkyz,rijkyx,rijkyz
  	real(kind=prec) :: ribkx !,rijkx
	  real(kind=prec) :: rajky !,rijky

	  real(kind=prec)    :: sg(6) !,   cnorm=R_ZERO !debug
	  real(kind=prec)    :: dx,dy,dz,cxm1,cx,cym1,cy,czm1,cz
	  complex(kind=prec) :: zh(2,2),ciuw,sum
	  ! integer ikk ! debug

    ciuw = ISIGN * CMPLX(0.0, 1.0, 8)*omega*MU_0 !没有用到此变量

    l = mGrid%Nx
	  m = mGrid%Ny
	  n = mGrid%Nz
	  nzn=n
	  n=n+1

    nhx=l*(m-1)*nzn
    nhy=(l-1)*m*nzn
    nhz=(l-1)*(m-1)*nzn
    np=nhx+nhy+nhz

!--------------------------------------------------
!   First, Hx Form
!--------------------------------------------------
    jj=0
    do i=1,l
    
      dx=mGrid%dx(i)
      
      do k=2,n
      
        dz=mGrid%delZ(k)
        czm1=mGrid%dz(k-1)
        
        do j=2,m
          
		      dy=mGrid%delY(j)
		      cym1=mGrid%dy(j-1)
		      cy=mGrid%dy(j)
		      		  

!         y element           
          ribcyx=rho_C(4)%v(i,j-1,k-1) !xy
		      ribcyz=rho_C(6)%v(i,j-1,k-1) !yz
		      rijcyx=rho_C(4)%v(i,j,k-1)   !xy
		      rijcyz=rho_C(6)%v(i,j,k-1)   !yz
!         z element
          ribczx=rho_C(5)%v(i,j-1,k-1) !xz
		      ribczy=rho_C(6)%v(i,j-1,k-1) !yz
		      rijczx=rho_C(5)%v(i,j,k-1)   !xz
		      rijczy=rho_C(6)%v(i,j,k-1)   !yz
! 
          rijcy=(rho_C(2)%v(i,j,k-1) * cy       &
                +rho_C(2)%v(i,j-1,k-1) * cym1)  &
                /(cym1+cy)

          if (k.ne.n) then

		        cz=mGrid%dz(k)

!           y element
		        ribkyx=rho_C(4)%v(i,j-1,k) !xy
		        ribkyz=rho_C(6)%v(i,j-1,k) !yz
		        rijkyx=rho_C(4)%v(i,j,k)   !xy
		        rijkyz=rho_C(6)%v(i,j,k)   !yz
!           z element
	          ribkzx=rho_C(5)%v(i,j-1,k) !xz
		        ribkzy=rho_C(6)%v(i,j-1,k) !yz
		        rijkzx=rho_C(5)%v(i,j,k)   !xz
		        rijkzy=rho_C(6)%v(i,j,k)   !yz
!
            rijky=(rho_C(2)%v(i,j,k) * cy       &
                  +rho_C(2)%v(i,j-1,k) * cym1)  &
                  /(cym1+cy)
!
            ribkz=(rho_C(3)%v(i,j-1,k) * cz       &
                  +rho_C(3)%v(i,j-1,k-1) * czm1)  &
                  /(czm1+cz)
            rijkz=(rho_C(3)%v(i,j,k) * cz       &
                  +rho_C(3)%v(i,j,k-1) * czm1)  &
                  /(czm1+cz)

          end if
        
          jj=jj+1
		  
		      sum=C_ZERO	
	
!-------------Hx--------------

! Hxibc   C1&B1
          if((j.eq.2).or.(k.eq.2)) then
            sum=sum-(ribcyz*dx/2)*inE%x(i,j-1,k-1)
          end if
! Hxijc   C2&B2
          if(k.eq.2) then
            sum=sum+(rijcy*dx*dy/czm1)*inE%x(i,j,k-1)
          end if
! Hxiec   C3&B3
          if((j.eq.m).or.(k.eq.2)) then
            sum=sum+(rijcyz*dx/2)*inE%x(i,j+1,k-1)
          end if
! Hxibk   C4
          if((j.eq.2).and.(k.ne.n)) then
            sum=sum+(ribkz*dx*dz/cym1)*inE%x(i,j-1,k)
          end if
! Hxiek   C6
          if((j.eq.m).and.(k.ne.n)) then
            sum=sum+(rijkz*dx*dz/cy)*inE%x(i,j+1,k)
          end if
! Hxibf   C7
          if((j.eq.2).and.(k.ne.n)) then
            sum=sum+(ribkyz*dx/2)*inE%x(i,j-1,k+1)
          end if
! Hxief   C9
          if((j.eq.m).and.(k.ne.n)) then
            sum=sum-(rijkyz*dx/2)*inE%x(i,j+1,k+1)
          end if

!-------------Hy--------------

! Hyibc   C10&B5
          if((i.eq.1).or.(k.eq.2)) then
            sum=sum-((ribcyx*dx*cym1/czm1  &
                -ribcyz*cym1-ribczx*dx)/4)*inE%y(i,j-1,k-1)
          end if
! Hydbc   C11&B6
          if((i.eq.l).or.(k.eq.2)) then
            sum=sum-((ribcyx*dx*cym1/czm1  &
                +ribcyz*cym1-ribczx*dx)/4)*inE%y(i+1,j-1,k-1)
          end if
! Hyibk   C12&***B7***
          if(i.eq.1) then
            if(k.ne.n) then
              sum=sum-(((ribczx-ribkzx)*dx+(ribkyz-ribcyz)*cym1  &
                  -(ribkyx/cz+ribcyx/czm1)*dx*cym1)/4  &
                  +ribkz*dz)*inE%y(i,j-1,k)
	          else
              do ia=1,6
                sg(ia)=sigma_C(ia)%v(i,j,k-1)
			        end do
              call zhalf(sg,zh)
              sum=sum-((ribczx*dx-ribcyz*cym1  &
                  -(zh(2,2)+ribcyx/czm1)*dx*cym1)/4)*inE%y(i,j-1,k)
	          end if
          end if
! Hydbk   C13&***B8***
          if(i.eq.l) then
            if(k.ne.n) then
              sum=sum-(((ribczx-ribkzx)*dx-(ribkyz-ribcyz)*cym1  &
                  -(ribkyx/cz+ribcyx/czm1)*dx*cym1)/4  &
                  -ribkz*dz)*inE%y(i+1,j-1,k)
	          else
              do ia=1,6
                sg(ia)=sigma_C(ia)%v(i,j,k-1)
			        end do
              call zhalf(sg,zh)
              sum=sum-((ribczx*dx+ribcyz*cym1  &
                  -(zh(2,2)+ribcyx/czm1)*dx*cym1)/4)*inE%y(i+1,j-1,k)
	          end if	        
          end if
! Hyibf   C14
          if((i.eq.1).and.(k.ne.n)) then
            sum=sum-((ribkyx*dx*cym1/cz  &
                +ribkyz*cym1+ribkzx*dx)/4)*inE%y(i,j-1,k+1)
          end if
! Hydbf   C15
          if((i.eq.l).and.(k.ne.n)) then
            sum=sum-((ribkyx*dx*cym1/cz  &
                -ribkyz*cym1+ribkzx*dx)/4)*inE%y(i+1,j-1,k+1)
          end if
! Hyijc   C16&B9
          if((i.eq.1).or.(k.eq.2)) then
            sum=sum-((rijcyx*dx*cy/czm1  &
                -rijcyz*cy+rijczx*dx)/4)*inE%y(i,j,k-1)
          end if
! Hydjc   C17&B10
          if((i.eq.l).or.(k.eq.2)) then
            sum=sum-((rijcyx*dx*cy/czm1  &
                +rijcyz*cy+rijczx*dx)/4)*inE%y(i+1,j,k-1)
          end if
! Hyijk   C18&***B11***
          if(i.eq.1) then
	          if(k.ne.n) then
              sum=sum-(((rijkzx-rijczx)*dx+(rijkyz-rijcyz)*cy  &
                  -(rijkyx/cz+rijcyx/czm1)*dx*cy)/4  &
				          -rijkz*dz)*inE%y(i,j,k)	        
			      else
              sum=sum+((rijczx*dx+rijcyz*cy  &
                  +(zh(2,2)+rijcyx/czm1)*dx*cy)/4)*inE%y(i,j,k)
			      end if	     
          end if
! Hydjk   C19&***B12***
          if(i.eq.l) then
            if(k.ne.n) then
              sum=sum-(((rijkzx-rijczx)*dx-(rijkyz-rijcyz)*cy  &
                  -(rijkyx/cz+rijcyx/czm1)*dx*cy)/4  &
                  +rijkz*dz)*inE%y(i+1,j,k)
	          else
              sum=sum+((rijczx*dx-rijcyz*cy  &
                  +(zh(2,2)+rijcyx/czm1)*dx*cy)/4)*inE%y(i+1,j,k)
			      end if
          end if
! Hyijf   C20
          if((i.eq.1).and.(k.ne.n)) then
            sum=sum-((rijkyx*dx*cy/cz  &
                +rijkyz*cy-rijkzx*dx)/4)*inE%y(i,j,k+1)
          end if
! Hydjf   C21
          if((i.eq.l).and.(k.ne.n)) then
            sum=sum-((rijkyx*dx*cy/cz  &
                -rijkyz*cy-rijkzx*dx)/4)*inE%y(i+1,j,k+1)
          end if

!-------------Hz--------------

! Hzibc   C22&B13
          if((i.eq.1).or.(j.eq.2)) then
            sum=sum-((ribczx*dx*czm1/cym1  &
                -ribczy*czm1-ribcyx*dx)/4)*inE%z(i,j-1,k-1)
          end if
! Hzdbc   C23&B14
          if((i.eq.l).or.(j.eq.2)) then
            sum=sum-((ribczx*dx*czm1/cym1  &
                +ribczy*czm1-ribcyx*dx)/4)*inE%z(i+1,j-1,k-1)
          end if
! Hzijc   C24&B15
          if(i.eq.1) then
            sum=sum-(((ribcyx-rijcyx)*dx+(rijczy-ribczy)*czm1  &
                -(rijczx/cy+ribczx/cym1)*dx*czm1)/4  &
                +rijcy*dy)*inE%z(i,j,k-1)
          end if
! Hzdjc   C25&B16
          if(i.eq.l) then
            sum=sum-(((ribcyx-rijcyx)*dx-(rijczy-ribczy)*czm1  &
                -(rijczx/cy+ribczx/cym1)*dx*czm1)/4  &
                -rijcy*dy)*inE%z(i+1,j,k-1)
          end if
! Hziec   C26&B17
          if((i.eq.1).or.(j.eq.m)) then
            sum=sum-((rijczx*dx*czm1/cy  &
                +rijczy*czm1+rijcyx*dx)/4)*inE%z(i,j+1,k-1)
          end if
! Hzdec   C27&B18
          if((i.eq.l).or.(j.eq.m)) then
            sum=sum-((rijczx*dx*czm1/cy  &
                -rijczy*czm1+rijcyx*dx)/4)*inE%z(i+1,j+1,k-1)
          end if
! Hzibk   C28
          if((k.ne.n).and.((i.eq.1).or.(j.eq.2))) then
            sum=sum-((ribkzx*dx*cz/cym1  &
                -ribkzy*cz+ribkyx*dx)/4)*inE%z(i,j-1,k)
          end if
! Hzdbk   C29
          if((k.ne.n).and.((i.eq.l).or.(j.eq.2))) then
            sum=sum-((ribkzx*dx*cz/cym1  &
                +ribkzy*cz+ribkyx*dx)/4)*inE%z(i+1,j-1,k)
          end if
! Hzijk   C30
          if((k.ne.n).and.(i.eq.1)) then
            sum=sum-(((rijkyx-ribkyx)*dx+(rijkzy-ribkzy)*cz  &
                -(rijkzx/cy+ribkzx/cym1)*dx*cz)/4  &
                -rijky*dy)*inE%z(i,j,k)
          end if
! Hzdjk   C31
          if((k.ne.n).and.(i.eq.l)) then
            sum=sum-(((rijkyx-ribkyx)*dx-(rijkzy-ribkzy)*cz  &
                -(rijkzx/cy+ribkzx/cym1)*dx*cz)/4  &
                +rijky*dy)*inE%z(i+1,j,k)
          end if
! Hziek   C32
          if((k.ne.n).and.((i.eq.1).or.(j.eq.m))) then
            sum=sum-((rijkzx*dx*cz/cy  &
                +rijkzy*cz-rijkyx*dx)/4)*inE%z(i,j+1,k)
          end if
! Hzdek   C33
          if((k.ne.n).and.((i.eq.l).or.(j.eq.m))) then
            sum=sum-((rijkzx*dx*cz/cy  &
                -rijkzy*cz-rijkyx*dx)/4)*inE%z(i+1,j+1,k)
          end if

!         bvec(jj)=sum
          OutE%x(i,j,k)=sum

!		  write(107,'(i7,2x,2f50.9)') jj,sum !debug
!          cnorm=cnorm+conjg(sum)*sum !debug

		    end do
	    end do
	  end do
    
!--------------------------------------------------
!   Second, Hy Form
!--------------------------------------------------
    do j=1,m
      
      dy=mGrid%dy(j)
      
      do k=2,n
        
        dz=mGrid%delZ(k)
		    czm1=mGrid%dz(k-1)
		    
        do i=2,l
          
          dx=mGrid%delX(i)
		      cxm1=mGrid%dx(i-1)
		      cx=mGrid%dx(i)
		      
!         x element
          rajcxy=rho_C(4)%v(i-1,j,k-1) !xy
		      rajcxz=rho_C(5)%v(i-1,j,k-1) !xz
		      rijcxy=rho_C(4)%v(i,j,k-1)
		      rijcxz=rho_C(5)%v(i,j,k-1)		  
!         z element
          rajczx=rho_C(5)%v(i-1,j,k-1) !xz
		      rajczy=rho_C(6)%v(i-1,j,k-1) !yz
		      rijczx=rho_C(5)%v(i,j,k-1)
		      rijczy=rho_C(6)%v(i,j,k-1)
!
          rijcx=(rho_C(1)%v(i,j,k-1) * cx       &
                +rho_C(1)%v(i-1,j,k-1) * cxm1)  &
                /(cx+cxm1)
	 
          if (k.ne.n) then

		        cz=mGrid%dz(k)

!           x element
		        rajkxy=rho_C(4)%v(i-1,j,k) !xy
		        rajkxz=rho_C(5)%v(i-1,j,k) !xz
		        rijkxy=rho_C(4)%v(i,j,k)
		        rijkxz=rho_C(5)%v(i,j,k)
!           z element
            rajkzx=rho_C(5)%v(i-1,j,k) !xz
			      rajkzy=rho_C(6)%v(i-1,j,k) !yz
            rijkzx=rho_C(5)%v(i,j,k)
			      rijkzy=rho_C(6)%v(i,j,k)
!
            rijkx=(rho_C(1)%v(i,j,k) * cx      &
                +rho_C(1)%v(i-1,j,k) * cxm1)   &
                /(cx+cxm1)
!
            rajkz=(rho_C(3)%v(i-1,j,k) * cz       &
                  +rho_C(3)%v(i-1,j,k-1) * czm1)  &
                  /(czm1+cz)
            rijkz=(rho_C(3)%v(i,j,k) * cz       &
                  +rho_C(3)%v(i,j,k-1) * czm1)  &
                  /(czm1+cz)

          end if
        
          jj=jj+1

		      sum=C_ZERO

!-------------Hy--------------

! Hyajc   C13&B9  
          if((i.eq.2).or.(k.eq.2)) then
            sum=sum-(rajcxz*dy/2)*inE%y(i-1,j,k-1)
          end if
! Hyijc   C14&B10 
          if(k.eq.2) then
            sum=sum+(rijcx*dx*dy/czm1)*inE%y(i,j,k-1)
          end if
! Hydjc   C15&B11 
          if((i.eq.l).or.(k.eq.2)) then
            sum=sum+(rijcxz*dy/2)*inE%y(i+1,j,k-1)
          end if
! Hyajk   C16
          if((i.eq.2).and.(k.ne.n)) then
            sum=sum+(rajkz*dy*dz/cxm1)*inE%y(i-1,j,k)
          end if
! Hydjk   C18
          if((i.eq.l).and.(k.ne.n)) then
            sum=sum+(rijkz*dy*dz/cx)*inE%y(i+1,j,k)
          end if
! Hyajf   C19
          if((i.eq.2).and.(k.ne.n)) then
            sum=sum+(rajkxz*dy/2)*inE%y(i-1,j,k+1)
          end if
! Hydjf   C21
          if((i.eq.l).and.(k.ne.n)) then
            sum=sum-(rijkxz*dy/2)*inE%y(i+1,j,k+1)
          end if  

!-------------Hz--------------

! Hzajc   C22&B13
          if((i.eq.2).or.(j.eq.1)) then
            sum=sum-((rajczy*dy*czm1/cxm1  &
                -rajczx*czm1-rajcxy*dy)/4)*inE%z(i-1,j,k-1)
          end if
! Hzijc   C23&B14
          if(j.eq.1) then
            sum=sum-(((rajcxy-rijcxy)*dy+(rijczx-rajczx)*czm1  &
                -(rijczy/cx+rajczy/cxm1)*dy*czm1)/4  &
                +rijcx*dx)*inE%z(i,j,k-1)
          end if  
! Hzdjc   C24&B15
          if((i.eq.l).or.(j.eq.1)) then
            sum=sum-((rijczy*dy*czm1/cx  &
                +rijczx*czm1+rijcxy*dy)/4)*inE%z(i+1,j,k-1)
          end if 
! Hzaec   C25&B16
          if((i.eq.2).or.(j.eq.m)) then
            sum=sum-((rajczy*dy*czm1/cxm1  &
                +rajczx*czm1-rajcxy*dy)/4)*inE%z(i-1,j+1,k-1)
          end if   
! Hziec   C26&B17
          if(j.eq.m) then
            sum=sum-(((rajcxy-rijcxy)*dy-(rijczx-rajczx)*czm1  &
                -(rijczy/cx+rajczy/cxm1)*dy*czm1)/4  &
                -rijcx*dx)*inE%z(i,j+1,k-1)
          end if
! Hzdec   C27&B18
          if((i.eq.l).or.(j.eq.m)) then
            sum=sum-((rijczy*dy*czm1/cx  &
                -rijczx*czm1+rijcxy*dy)/4)*inE%z(i+1,j+1,k-1)
          end if   
! Hzajk   C28
          if(((i.eq.2).or.(j.eq.1)).and.(k.ne.n)) then
            sum=sum-((rajkzy*dy*cz/cxm1  &
                -rajkzx*cz+rajkxy*dy)/4)*inE%z(i-1,j,k)
          end if    
! Hzijk   C29
          if((j.eq.1).and.(k.ne.n)) then
            sum=sum-(((rijkxy-rajkxy)*dy+(rijkzx-rajkzx)*cz  &
                -(rijkzy/cx+rajkzy/cxm1)*dy*cz)/4  &
                -rijkx*dx)*inE%z(i,j,k)
          end if 
! Hzdjk   C30
          if(((i.eq.l).or.(j.eq.1)).and.(k.ne.n)) then
            sum=sum-((rijkzy*dy*cz/cx  &
                +rijkzx*cz-rijkxy*dy)/4)*inE%z(i+1,j,k)
          end if 
! Hzaek   C31
          if(((i.eq.2).or.(j.eq.m)).and.(k.ne.n)) then
            sum=sum-((rajkzy*dy*cz/cxm1  &
                +rajkzx*cz+rajkxy*dy)/4)*inE%z(i-1,j+1,k)
          end if    
! Hziek   C32
          if((j.eq.m).and.(k.ne.n)) then
            sum=sum-(((rijkxy-rajkxy)*dy-(rijkzx-rajkzx)*cz  &
                -(rijkzy/cx+rajkzy/cxm1)*dy*cz)/4  &
                +rijkx*dx)*inE%z(i,j+1,k)
          end if 
! Hzdek   C33
          if(((i.eq.l).or.(j.eq.m)).and.(k.ne.n)) then
            sum=sum-((rijkzy*dy*cz/cx  &
                -rijkzx*cz-rijkxy*dy)/4)*inE%z(i+1,j+1,k)
          end if

!-------------Hx--------------

! Hxajc   C1&B1
          if((j.eq.1).or.(k.eq.2)) then
            sum=sum-((rajcxy*cxm1*dy/czm1  &
                -rajcxz*cxm1-rajczy*dy)/4)*inE%x(i-1,j,k-1)
          end if
! Hxaec   C2&B2
          if((j.eq.m).or.(k.eq.2)) then
            sum=sum-((rajcxy*cxm1*dy/czm1  &
                +rajcxz*cxm1-rajczy*dy)/4)*inE%x(i-1,j+1,k-1)
          end if  
! Hxajk   C3&***B3***
          if(j.eq.1) then
            if(k.ne.n) then
              sum=sum-(((rajczy-rajkzy)*dy+(rajkxz-rajcxz)*cxm1  &
                  -(rajkxy/cz+rajcxy/czm1)*cxm1*dy)/4  &
                  +rajkz*dz)*inE%x(i-1,j,k)
	          else
              do ia=1,6
                sg(ia)=sigma_C(ia)%v(i,j,k-1)
			        end do
              call zhalf(sg,zh)
              sum=sum-((rajczy*dy-rajcxz*cxm1  &
                  +(zh(1,1)-rajcxy/czm1)*cxm1*dy)/4)*inE%x(i-1,j,k)
	          end if
          end if
! Hxaek   C4&***B4***
          if(j.eq.m) then
            if(k.ne.n) then
              sum=sum-(((rajczy-rajkzy)*dy-(rajkxz-rajcxz)*cxm1  &
                  -(rajkxy/cz+rajcxy/czm1)*cxm1*dy)/4  &
                  -rajkz*dz)*inE%x(i-1,j+1,k)
	          else
              do ia=1,6
                sg(ia)=sigma_C(ia)%v(i,j,k-1)
			        end do
              call zhalf(sg,zh)
              sum=sum-((rajczy*dy+rajcxz*cxm1  &
                  +(zh(1,1)-rajcxy/czm1)*cxm1*dy)/4)*inE%x(i-1,j+1,k)
	          end if
          end if
! Hxajf   C5
          if((j.eq.1).and.(k.ne.n)) then
            sum=sum-((rajkxy*cxm1*dy/cz  &
                +rajkxz*cxm1+rajkzy*dy)/4)*inE%x(i-1,j,k+1)
          end if
! Hxaef   C6
          if((j.eq.m).and.(k.ne.n)) then
            sum=sum-((rajkxy*cxm1*dy/cz  &
                -rajkxz*cxm1+rajkzy*dy)/4)*inE%x(i-1,j+1,k+1)
          end if
! Hxijc   C7&B5
          if((j.eq.1).or.(k.eq.2)) then
            sum=sum-((rijcxy*cx*dy/czm1  &
                -rijcxz*cx+rijczy*dy)/4)*inE%x(i,j,k-1)
          end if
! Hxiec   C8&B6
          if((j.eq.m).or.(k.eq.2)) then
            sum=sum-((rijcxy*cx*dy/czm1  &
                +rijcxz*cx+rijczy*dy)/4)*inE%x(i,j+1,k-1)
          end if
! Hxijk   C9&***B7***
	        if(j.eq.1) then
	          if(k.ne.n) then
              sum=sum-(((rijkzy-rijczy)*dy+(rijkxz-rijcxz)*cx  &
                  -(rijkxy/cz+rijcxy/czm1)*cx*dy)/4  &
                  -rijkz*dz)*inE%x(i,j,k)
	          else
              sum=sum+((rijczy*dy+rijcxz*cx  &
                  -(zh(1,1)-rijcxy/czm1)*cx*dy)/4)*inE%x(i,j,k)
	          end if
          end if
! Hxiek   C10&***B8***
          if(j.eq.m) then
            if(k.ne.n) then
              sum=sum-(((rijkzy-rijczy)*dy-(rijkxz-rijcxz)*cx  &
                  -(rijkxy/cz+rijcxy/czm1)*cx*dy)/4  &
                  +rijkz*dz)*inE%x(i,j+1,k)
	          else
              sum=sum+((rijczy*dy-rijcxz*cx  &
                  -(zh(1,1)-rijcxy/czm1)*cx*dy)/4)*inE%x(i,j+1,k)
	          end if
          end if
! Hxijf   C11
          if((j.eq.1).and.(k.ne.n)) then
            sum=sum-((rijkxy*cx*dy/cz  &
                +rijkxz*cx-rijkzy*dy)/4)*inE%x(i,j,k+1)
          end if
! Hxief   C12
          if((j.eq.m).and.(k.ne.n)) then
            sum=sum-((rijkxy*cx*dy/cz  &
                -rijkxz*cx-rijkzy*dy)/4)*inE%x(i,j+1,k+1)
          end if

!         bvec(jj)=sum
          OutE%y(i,j,k)=sum

!          cnorm=cnorm+conjg(sum)*sum !debug
!		  write(107,'(i7,2x,2f50.9)') jj,sum !debug

		    end do
	    end do
	  end do
	
!--------------------------------------------------
!   Third, Hz Form
!--------------------------------------------------	
    do k=1,n-1
      
      dz=mGrid%dz(k)
      
      do j=2,m
        
        dy=mGrid%delY(j)
        cym1=mGrid%dy(j-1)
		    cy=mGrid%dy(j)
		    
        do i=2,l
          
          dx=mGrid%delX(i)
		      cxm1=mGrid%dx(i-1)
		      cx=mGrid%dx(i)
		      
!         x element
          rabkxy=rho_C(4)%v(i-1,j-1,k)
		      rabkxz=rho_C(5)%v(i-1,j-1,k)
	        rajkxy=rho_C(4)%v(i-1,j,k)
		      rajkxz=rho_C(5)%v(i-1,j,k)
!
		      ribkxy=rho_C(4)%v(i,j-1,k)
		      ribkxz=rho_C(5)%v(i,j-1,k)
		      rijkxy=rho_C(4)%v(i,j,k)
		      rijkxz=rho_C(5)%v(i,j,k)
!
          ribkx=(rho_C(1)%v(i,j-1,k) * cx       &
                +rho_C(1)%v(i-1,j-1,k) * cxm1)  &
                /(cx+cxm1)
		      rijkx=(rho_C(1)%v(i,j,k) * cx       &
                +rho_C(1)%v(i-1,j,k) * cxm1)  &
                /(cx+cxm1)
!         y element
          rabkyx=rho_C(4)%v(i-1,j-1,k)
		      rabkyz=rho_C(6)%v(i-1,j-1,k)
		      rajkyx=rho_C(4)%v(i-1,j,k)
		      rajkyz=rho_C(6)%v(i-1,j,k)
!
		      ribkyx=rho_C(4)%v(i,j-1,k)
		      ribkyz=rho_C(6)%v(i,j-1,k)
		      rijkyx=rho_C(4)%v(i,j,k)
		      rijkyz=rho_C(6)%v(i,j,k)
!
          rajky=(rho_C(2)%v(i-1,j,k) * cy       &
                +rho_C(2)%v(i-1,j-1,k) * cym1)  &
                /(cym1+cy)
		      rijky=(rho_C(2)%v(i,j,k) * cy       &
                +rho_C(2)%v(i,j-1,k) * cym1)  &
                /(cym1+cy)
        
          jj=jj+1

          sum=C_ZERO

!-------------Hz--------------

! Hzabk   C25
          if((i.eq.2).or.(j.eq.2)) then
            sum=sum-(rabkxy*dz/2)*inE%z(i-1,j-1,k)
          end if
! Hzibk   C26
          if(j.eq.2) then
            sum=sum+(ribkx*dx*dz/cym1)*inE%z(i,j-1,k)
          end if
! Hzdbk   C27
          if((i.eq.l).or.(j.eq.2)) then
            sum=sum+(ribkxy*dz/2)*inE%z(i+1,j-1,k)
          end if
! Hzajk   C28
          if(i.eq.2) then
            sum=sum+(rajky*dy*dz/cxm1)*inE%z(i-1,j,k)
          end if
! Hzdjk   C30
          if(i.eq.l) then
            sum=sum+(rijky*dy*dz/cx)*inE%z(i+1,j,k)
          end if 
! Hzaek   C31
          if((i.eq.2).or.(j.eq.m)) then
            sum=sum+(rajkxy*dz/2)*inE%z(i-1,j+1,k)
          end if   
! Hziek   C32
          if(j.eq.m) then
            sum=sum+(rijkx*dx*dz/cy)*inE%z(i,j+1,k)
          end if
! Hzdek   C33
          if((i.eq.l).or.(j.eq.m)) then
            sum=sum-(rijkxy*dz/2)*inE%z(i+1,j+1,k)
          end if 

!-------------Hx--------------

! Hxabk   C1
          if((j.eq.2).or.(k.eq.1)) then
            sum=sum-((rabkxz*cxm1*dz/cym1  &
                -rabkxy*cxm1-rabkyz*dz)/4)*inE%x(i-1,j-1,k)
          end if
! Hxajk   C2
          if(k.eq.1) then
            sum=sum-(((rabkyz-rajkyz)*dz+(rajkxy-rabkxy)*cxm1  &
                -(rajkxz/cy+rabkxz/cym1)*cxm1*dz)/4  &
                +rajky*dy)*inE%x(i-1,j,k)
          end if
! Hxaek   C3
          if((j.eq.m).or.(k.eq.1)) then
            sum=sum-((rajkxz*cxm1*dz/cy  &
                +rajkxy*cxm1+rajkyz*dz)/4)*inE%x(i-1,j+1,k)
          end if
! Hxabf   C4
          if(j.eq.2) then
            sum=sum-((rabkxz*cxm1*dz/cym1  &
                +rabkxy*cxm1-rabkyz*dz)/4)*inE%x(i-1,j-1,k+1)
          end if
! Hxaef   C6
          if(j.eq.m) then
            sum=sum-((rajkxz*cxm1*dz/cy  &
                -rajkxy*cxm1+rajkyz*dz)/4)*inE%x(i-1,j+1,k+1)
          end if
! Hxibk   C7
          if((j.eq.2).or.(k.eq.1)) then
            sum=sum-((ribkxz*cx*dz/cym1  &
                -ribkxy*cx+ribkyz*dz)/4)*inE%x(i,j-1,k)
          end if
! Hxijk   C8
          if(k.eq.1) then
            sum=sum-(((rijkyz-ribkyz)*dz+(rijkxy-ribkxy)*cx  &
                -(rijkxz/cy+ribkxz/cym1)*cx*dz)/4  &
                -rijky*dy)*inE%x(i,j,k)
          end if
! Hxiek   C9
          if((j.eq.m).or.(k.eq.1)) then
            sum=sum-((rijkxz*cx*dz/cy  &
                +rijkxy*cx-rijkyz*dz)/4)*inE%x(i,j+1,k)
          end if
! Hxibf   C10
          if(j.eq.2) then
            sum=sum-((ribkxz*cx*dz/cym1  &
                +ribkxy*cx+ribkyz*dz)/4)*inE%x(i,j-1,k+1)
          end if
! Hxief   C12
          if(j.eq.m) then
            sum=sum-((rijkxz*cx*dz/cy  &
                -rijkxy*cx-rijkyz*dz)/4)*inE%x(i,j+1,k+1)
          end if

!-------------Hy--------------

! Hyabk   C13
          if((i.eq.2).or.(k.eq.1)) then
            sum=sum-((rabkyz*cym1*dz/cxm1  &
                -rabkyx*cym1-rabkxz*dz)/4)*inE%y(i-1,j-1,k)
          end if
! Hyibk   C14
          if(k.eq.1) then
            sum=sum-(((rabkxz-ribkxz)*dz+(ribkyx-rabkyx)*cym1  &
                -(ribkyz/cx+rabkyz/cxm1)*cym1*dz)/4  &
                +ribkx*dx)*inE%y(i,j-1,k)
          end if
! Hydbk   C15
          if((i.eq.l).or.(k.eq.1)) then
            sum=sum-((ribkyz*cym1*dz/cx  &
                +ribkyx*cym1+ribkxz*dz)/4)*inE%y(i+1,j-1,k)
          end if
! Hyabf   C16
          if(i.eq.2) then
            sum=sum-((rabkyz*cym1*dz/cxm1  &
                +rabkyx*cym1-rabkxz*dz)/4)*inE%y(i-1,j-1,k+1)
          end if
! Hydbf   C18
          if(i.eq.l) then
            sum=sum-((ribkyz*cym1*dz/cx  &
                -ribkyx*cym1+ribkxz*dz)/4)*inE%y(i+1,j-1,k+1)
          end if
! Hyajk   C19
          if((i.eq.2).or.(k.eq.1)) then
            sum=sum-((rajkyz*cy*dz/cxm1  &
                -rajkyx*cy+rajkxz*dz)/4)*inE%y(i-1,j,k)
          end if
! Hyijk   C20
          if(k.eq.1) then
            sum=sum-(((rijkxz-rajkxz)*dz+(rijkyx-rajkyx)*cy  &
                -(rijkyz/cx+rajkyz/cxm1)*cy*dz)/4  &
                -rijkx*dx)*inE%y(i,j,k)
          end if
! Hydjk   C21
          if((i.eq.l).or.(k.eq.1)) then
            sum=sum-((rijkyz*cy*dz/cx  &
                +rijkyx*cy-rijkxz*dz)/4)*inE%y(i+1,j,k)
          end if
! Hyajf   C22
          if(i.eq.2) then
            sum=sum-((rajkyz*cy*dz/cxm1  &
                +rajkyx*cy+rajkxz*dz)/4)*inE%y(i-1,j,k+1)
          end if
! Hydjf   C24
          if(i.eq.l) then
            sum=sum-((rijkyz*cy*dz/cx  &
                -rijkyx*cy-rijkxz*dz)/4)*inE%y(i+1,j,k+1)
          end if
		  
!         bvec(jj)=sum
          OutE%z(i,j,k)=sum

!          cnorm=cnorm+conjg(sum)*sum !debug
!		  write(107,'(i7,2x,2f50.9)') jj,sum !debug

		    end do
	    end do
	  end do


!    do k=1,n
!      do j=1,m+1
!        do i=1,l
!           write(107,'(2f30.9)') OutE%x(i,j,k) !debug
!		end do
!	  end do
!	end do
!    write(1010,'(f50.9)') sqrt(cnorm)
!	pause !debug

  end subroutine calcb


  ! impedance of anisotropy half space
  subroutine zhalf(sg,z)
    implicit none
	  real(kind=prec)    :: sg(6)
	  complex(kind=prec) :: z(2,2)
    ! local variables
	  real(kind=prec)    :: a1,a2,bs
	  real(kind=prec)    :: axx,axy,ayy
	  real(kind=prec)    :: da0,da,sa0,db0,doubletiny
	  complex(kind=prec) :: ic,commi,k1,k2
	  complex(kind=prec) :: zrot(2,2)
    !
	  ic = ISIGN * cmplx(0.0d0,1.0d0) 
    commi=(1.0d0+ic)*dsqrt(0.5d0*omega*MU_0)     
	  axx=sg(1)-sg(5)*sg(5)/sg(3) 
	  axy=sg(4)-sg(5)*sg(6)/sg(3)
	  ayy=sg(2)-sg(6)*sg(6)/sg(3)
	  sa0=axx+ayy
	  da0=axx-ayy
	  da=dsqrt(da0*da0+4.d0*axy*axy)                         
	  a1=0.5d0*(sa0+da)                  
	  a2=0.5d0*(sa0-da) 
	  k1=commi*dsqrt(a1)
	  k2=commi*dsqrt(a2)
!
    if(da.lt.tiny(doubletiny))then
	    bs=0.d0
	  else                                                
	    db0=da0/da
	    bs=0.5d0*dacos(db0)
	    if(axy.lt.0.d0) bs=-bs
	  endif
!
    zrot(1,1)=0.d0
    zrot(1,2)=k1/a1
    zrot(2,1)=-k2/a2
    zrot(2,2)=0.d0
!
    call rotz(zrot,-bs,z)

  end subroutine zhalf
  

  subroutine rotz(za,betrad,zb)
! =============================
! Rotates the impedance ZA by BETRAD (in radians) to obtain ZB
!
   implicit none
   real(kind=prec)      :: betrad
   complex(kind=prec)   :: za(2,2),zb(2,2)
!
   real(kind=prec)      :: co2,si2
   complex(kind=prec)   :: sum1,sum2,dif1,dif2
!
   co2=dcos(2.d0*betrad)
   si2=dsin(2.d0*betrad)
!
   sum1=za(1,1)+za(2,2)
   sum2=za(1,2)+za(2,1)
   dif1=za(1,1)-za(2,2)
   dif2=za(1,2)-za(2,1)
!
   zb(1,1)=0.5d0*(sum1+dif1*co2+sum2*si2)
   zb(1,2)=0.5d0*(dif2+sum2*co2-dif1*si2)
   zb(2,1)=0.5d0*(-dif2+sum2*co2-dif1*si2)
   zb(2,2)=0.5d0*(sum1-dif1*co2-sum2*si2)

  end subroutine rotz
  
  
    
  ! ********** 3D forward modelling operator for MT anisotropic problem **********************
  subroutine CurlcurleSetUp()
    implicit none
    integer                     :: status     ! for dynamic memory allocation
    integer                     :: ia,i,j,k,l,m,n,nzn ! dummy variables
    integer                     :: nhx,nhy,nhz,np,np5
	  integer                     :: ii,jj,kk

!   for Hx>>
    real(kind=prec) :: ribcyx,ribcyz,ribkyx,ribkyz,rijcyx,rijcyz,rijkyx,rijkyz !total 8
		real(kind=prec) :: ribczx,ribczy,ribkzx,ribkzy,rijczx,rijczy,rijkzx,rijkzy 
		real(kind=prec) :: rijcy,rijky
		real(kind=prec) :: ribkz,rijkz
!   for Hy>>
		real(kind=prec) :: rajcxy,rajcxz,rajkxy,rajkxz,rijcxy,rijcxz,rijkxy,rijkxz
		real(kind=prec) :: rajczx,rajczy,rajkzx,rajkzy !,rijczx,rijczy,rijkzx,rijkzy
		real(kind=prec) :: rijcx,rijkx
		real(kind=prec) :: rajkz !,rijkz
!   for Hz>>
		real(kind=prec) :: rabkxy,ribkxy !,rajkxy,rajkxz,rabkxz,ribkxz,rijkxy,rijkxz
!   real(kind=prec) :: rabkyx,rabkyz,rajkyx,rajkyz,ribkyx,ribkyz,rijkyx,rijkyz
		real(kind=prec) :: ribkx !,rijkx
		real(kind=prec) :: rajky !,rijky

		real(kind=prec)    :: sg(6)
		real(kind=prec)    :: dx,dy,dz,cxm1,cx,cym1,cy,czm1,cz
		complex(kind=prec) :: zh(2,2),ciuw
		! integer ikk ! debug

    ciuw = ISIGN * CMPLX(0.0d0,omega*MU_0)

    l = mGrid%Nx
		m = mGrid%Ny
		n = mGrid%Nz
		nzn=n
		n=nzn+1

    nhx=l*(m-1)*nzn
    nhy=(l-1)*m*nzn
    nhz=(l-1)*(m-1)*nzn
    np=nhx+nhy+nhz

    np5=l*(m-1)*nzn+(l-1)*m*nzn+3*l*(m-2)*(nzn-1)+l*(m-1)*(nzn-1)+  &
        12*(l-1)*(m-1)*(nzn-1)+9*(l-1)*(m-1)*nzn+                   &
				3*(l-2)*m*(nzn-1)+(l-1)*m*(nzn-1)+5*(l-2)*(m-1)*nzn+        &
        5*(l-1)*(m-2)*nzn+4*(l-1)*(m-2)*(nzn-1)+                    &
        4*(l-2)*(m-1)*(nzn-1)+2*(l-2)*(m-2)*nzn+1

		allocate(spa(np5), STAT=status)
		allocate(ija(np5), STAT=status)
	 
    spa=C_ZERO
		ija=0

!--------------------------------------------------
!   First, Hx Form
!--------------------------------------------------
    jj=0
    ija(1)=np+2
    kk=np+1
    do i=1,l
      
      dx=mGrid%dx(i)
      
      do k=2,n
        
        dz=mGrid%delZ(k)
        czm1=mGrid%dz(k-1)
        
        do j=2,m
          
					dy=mGrid%delY(j)
					cym1=mGrid%dy(j-1)
					cy=mGrid%dy(j)
							  

!         y element           
          ribcyx=rho_C(4)%v(i,j-1,k-1) !xy
					ribcyz=rho_C(6)%v(i,j-1,k-1) !yz
					rijcyx=rho_C(4)%v(i,j,k-1)   !xy
					rijcyz=rho_C(6)%v(i,j,k-1)   !yz
!         z element
          ribczx=rho_C(5)%v(i,j-1,k-1) !xz
					ribczy=rho_C(6)%v(i,j-1,k-1) !yz
					rijczx=rho_C(5)%v(i,j,k-1)   !xz
					rijczy=rho_C(6)%v(i,j,k-1)   !yz
! 
          rijcy=(rho_C(2)%v(i,j,k-1) * cy       &
                +rho_C(2)%v(i,j-1,k-1) * cym1)  &
                /(cym1+cy)

          if (k.ne.n) then

						cz=mGrid%dz(k)

!           y element
						ribkyx=rho_C(4)%v(i,j-1,k) !xy
						ribkyz=rho_C(6)%v(i,j-1,k) !yz
						rijkyx=rho_C(4)%v(i,j,k)   !xy
						rijkyz=rho_C(6)%v(i,j,k)   !yz
!           z element
						ribkzx=rho_C(5)%v(i,j-1,k) !xz
						ribkzy=rho_C(6)%v(i,j-1,k) !yz
						rijkzx=rho_C(5)%v(i,j,k)   !xz
						rijkzy=rho_C(6)%v(i,j,k)   !yz
!
            rijky=(rho_C(2)%v(i,j,k) * cy       &
                  +rho_C(2)%v(i,j-1,k) * cym1)  &
                  /(cym1+cy)
!
            ribkz=(rho_C(3)%v(i,j-1,k) * cz       &
                  +rho_C(3)%v(i,j-1,k-1) * czm1)  &
                  /(czm1+cz)
            rijkz=(rho_C(3)%v(i,j,k) * cz       &
                  +rho_C(3)%v(i,j,k-1) * czm1)  &
                  /(czm1+cz)

          end if

        
          jj=jj+1
!*****************************CsDone*****************************************
!------
! Hxijk   C5&***B4***
          if(k.ne.n) then
            spa(jj)=ciuw*dx*dy*dz+(rijcyz-rijkyz+ribkyz-ribcyz)*dx/2   &
                    +(rijky/cz+rijcy/czm1)*dx*dy+(rijkz/cy+ribkz/cym1)*dx*dz                    
					else
!++++
            do ia=1,6
              sg(ia)=sigma_C(ia)%v(i,j,k-1)
						end do
            call zhalf(sg,zh)		                  
            spa(jj)=ciuw*dx*dy*dz+rijcy*dx*dy/czm1-zh(2,1)*dx*dy   &
                    +(rijcyz-ribcyz)*dx/2
!			write(*,*) jj,-zh(2,1) ! debug
          end if
! Hxiek   C6
          if((j.ne.m).and.(k.ne.n)) then
            kk=kk+1
            ii=(m-1)*(n-1)*(i-1)+(m-1)*(k-2)+j
            spa(kk)=-rijkz*dx*dz/cy
            ija(kk)=ii
          end if
!------
! Hxibf   C7
          if((j.ne.2).and.(k.ne.n)) then
            kk=kk+1
            ii=(m-1)*(n-1)*(i-1)+(m-1)*(k-1)+j-2
            spa(kk)=-ribkyz*dx/2
            ija(kk)=ii
          end if
! Hxijf   C8
          if(k.ne.n) then
            kk=kk+1
            ii=(m-1)*(n-1)*(i-1)+(m-1)*(k-1)+j-1
            spa(kk)=-rijky*dx*dy/cz
            ija(kk)=ii
          end if
! Hxief   C9
          if((j.ne.m).and.(k.ne.n)) then
            kk=kk+1
            ii=(m-1)*(n-1)*(i-1)+(m-1)*(k-1)+j
            spa(kk)=rijkyz*dx/2
            ija(kk)=ii
          end if           
!*********************************CsDone*************************************
!------
! Hyibc   C10&B7
          if((i.ne.1).and.(k.ne.2)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-2)+(l-1)*(k-3)+i-1
            spa(kk)=(ribcyx*dx*cym1/czm1-ribcyz*cym1-ribczx*dx)/4
            ija(kk)=ii
          end if
! Hydbc   C11&B8
          if((i.ne.l).and.(k.ne.2)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-2)+(l-1)*(k-3)+i
            spa(kk)=(ribcyx*dx*cym1/czm1+ribcyz*cym1-ribczx*dx)/4
            ija(kk)=ii
          end if
!------
! Hyibk   C12&***B7***
          if(i.ne.1) then
						if(k.ne.n) then
              kk=kk+1
              ii=nhx+(l-1)*(n-1)*(j-2)+(l-1)*(k-2)+i-1
              spa(kk)=((ribczx-ribkzx)*dx+(ribkyz-ribcyz)*cym1  &
                      -(ribkyx/cz+ribcyx/czm1)*dx*cym1)/4+ribkz*dz
              ija(kk)=ii
						else
              kk=kk+1
              ii=nhx+(l-1)*(n-1)*(j-2)+(l-1)*(k-2)+i-1
              spa(kk)=(ribczx*dx-ribcyz*cym1  &
                     -(zh(2,2)+ribcyx/czm1)*dx*cym1)/4
              ija(kk)=ii
            end if
          end if
! Hydbk   C13&***B8***
          if(i.ne.l) then
            if(k.ne.n) then
              kk=kk+1
              ii=nhx+(l-1)*(n-1)*(j-2)+(l-1)*(k-2)+i
              spa(kk)=((ribczx-ribkzx)*dx-(ribkyz-ribcyz)*cym1  &
                      -(ribkyx/cz+ribcyx/czm1)*dx*cym1)/4-ribkz*dz
              ija(kk)=ii
						else
              kk=kk+1
              ii=nhx+(l-1)*(n-1)*(j-2)+(l-1)*(k-2)+i
              spa(kk)=(ribczx*dx+ribcyz*cym1  &
                      -(zh(2,2)+ribcyx/czm1)*dx*cym1)/4
              ija(kk)=ii
            end if
          end if
!------
! Hyibf   C14
          if((i.ne.1).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-2)+(l-1)*(k-1)+i-1
            spa(kk)=(ribkyx*dx*cym1/cz+ribkyz*cym1+ribkzx*dx)/4
            ija(kk)=ii
          end if
! Hydbf   C15
          if((i.ne.l).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-2)+(l-1)*(k-1)+i
            spa(kk)=(ribkyx*dx*cym1/cz-ribkyz*cym1+ribkzx*dx)/4
            ija(kk)=ii
          end if
!------
! Hyijc   C16&B9
          if((i.ne.1).and.(k.ne.2)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-3)+i-1
            spa(kk)=(rijcyx*dx*cy/czm1-rijcyz*cy+rijczx*dx)/4
            ija(kk)=ii
          end if            
! Hydjc   C17&B10
          if((i.ne.l).and.(k.ne.2)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-3)+i
            spa(kk)=(rijcyx*dx*cy/czm1+rijcyz*cy+rijczx*dx)/4
            ija(kk)=ii
          end if 
!------
! Hyijk   C18&***B11***
          if(i.ne.1) then
						if(k.ne.n) then
              kk=kk+1
              ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-2)+i-1
              spa(kk)=((rijkzx-rijczx)*dx+(rijkyz-rijcyz)*cy  &
                      -(rijkyx/cz+rijcyx/czm1)*dx*cy)/4-rijkz*dz
              ija(kk)=ii
						else
              kk=kk+1
              ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-2)+i-1
              spa(kk)=-(rijczx*dx+rijcyz*cy  &
                      +(zh(2,2)+rijcyx/czm1)*dx*cy)/4
              ija(kk)=ii
						end if
          end if
! Hydjk   C19&***B12***
          if(i.ne.l) then
            if(k.ne.n) then
              kk=kk+1
              ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-2)+i
              spa(kk)=((rijkzx-rijczx)*dx-(rijkyz-rijcyz)*cy  &
                      -(rijkyx/cz+rijcyx/czm1)*dx*cy)/4+rijkz*dz
              ija(kk)=ii
						else
              kk=kk+1
              ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-2)+i
              spa(kk)=-(rijczx*dx-rijcyz*cy  &
                      +(zh(2,2)+rijcyx/czm1)*dx*cy)/4
              ija(kk)=ii
						end if
          end if
!------
! Hyijf   C20		  
          if((i.ne.1).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-1)+i-1
            spa(kk)=(rijkyx*dx*cy/cz+rijkyz*cy-rijkzx*dx)/4
            ija(kk)=ii
          end if
! Hydjf   C21
          if((i.ne.l).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-1)+i
            spa(kk)=(rijkyx*dx*cy/cz-rijkyz*cy-rijkzx*dx)/4
            ija(kk)=ii
          end if		  
!*********************************CsDone*************************************
!------
! Hzibc   C22&B13
          if((i.ne.1).and.(j.ne.2)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-3)+i-1
            spa(kk)=(ribczx*dx*czm1/cym1-ribczy*czm1-ribcyx*dx)/4
            ija(kk)=ii
          end if
! Hzdbc   C23&B14
          if((i.ne.l).and.(j.ne.2)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-3)+i
            spa(kk)=(ribczx*dx*czm1/cym1+ribczy*czm1-ribcyx*dx)/4
            ija(kk)=ii
          end if
!------		  		  
! Hzijc   C24&B15
          if(i.ne.1) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-2)+i-1
            spa(kk)=((ribcyx-rijcyx)*dx+(rijczy-ribczy)*czm1  &
                    -(rijczx/cy+ribczx/cym1)*dx*czm1)/4+rijcy*dy
            ija(kk)=ii
          end if
! Hzdjc   C25&B16
          if(i.ne.l) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-2)+i
            spa(kk)=((ribcyx-rijcyx)*dx-(rijczy-ribczy)*czm1  &
                   -(rijczx/cy+ribczx/cym1)*dx*czm1)/4-rijcy*dy
            ija(kk)=ii
          end if
!------
! Hziec   C26&B17
          if((i.ne.1).and.(j.ne.m)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-1)+i-1
            spa(kk)=(rijczx*dx*czm1/cy+rijczy*czm1+rijcyx*dx)/4
            ija(kk)=ii
          end if
! Hzdec   C27&B18
          if((i.ne.l).and.(j.ne.m)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-1)+i
            spa(kk)=(rijczx*dx*czm1/cy-rijczy*czm1+rijcyx*dx)/4
            ija(kk)=ii
          end if
!------
! Hzibk   C28
          if((k.ne.n).and.(i.ne.1).and.(j.ne.2)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-3)+i-1
            spa(kk)=(ribkzx*dx*cz/cym1-ribkzy*cz+ribkyx*dx)/4
            ija(kk)=ii
          end if
! Hzdbk   C29
          if((k.ne.n).and.(i.ne.l).and.(j.ne.2)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-3)+i
            spa(kk)=(ribkzx*dx*cz/cym1+ribkzy*cz+ribkyx*dx)/4
            ija(kk)=ii
          end if
!------
! Hzijk   C30
          if((k.ne.n).and.(i.ne.1)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-2)+i-1
            spa(kk)=((rijkyx-ribkyx)*dx+(rijkzy-ribkzy)*cz  &
                    -(rijkzx/cy+ribkzx/cym1)*dx*cz)/4-rijky*dy
            ija(kk)=ii
          end if
! Hzdjk   C31
          if((k.ne.n).and.(i.ne.l)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-2)+i
            spa(kk)=((rijkyx-ribkyx)*dx-(rijkzy-ribkzy)*cz  &
                    -(rijkzx/cy+ribkzx/cym1)*dx*cz)/4+rijky*dy
            ija(kk)=ii
          end if
!------
! Hziek   C32
          if((k.ne.n).and.(i.ne.1).and.(j.ne.m)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i-1
            spa(kk)=(rijkzx*dx*cz/cy+rijkzy*cz-rijkyx*dx)/4
            ija(kk)=ii
          end if
! Hzdek   C33
          if((k.ne.n).and.(i.ne.l).and.(j.ne.m)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i
            spa(kk)=(rijkzx*dx*cz/cy-rijkzy*cz-rijkyx*dx)/4
            ija(kk)=ii
          end if

          ija(jj+1)=kk+1

				end do
			end do
		end do

!   debug
!   do i=1,kk
!      write(901,*) spa(i)
!      write(902,*) ija(i)
!   end do
!   pause
!    ikk=kk ! debug

!--------------------------------------------------
!   Second, Hy Form
!--------------------------------------------------
    do j=1,m
      
      dy=mGrid%dy(j)
      
      do k=2,n
        
        dz=mGrid%delZ(k)
        czm1=mGrid%dz(k-1)
        
        do i=2,l
        
          dx=mGrid%delX(i)
          cxm1=mGrid%dx(i-1)
					cx=mGrid%dx(i)
					
!         x element
          rajcxy=rho_C(4)%v(i-1,j,k-1) !xy
					rajcxz=rho_C(5)%v(i-1,j,k-1) !xz
					rijcxy=rho_C(4)%v(i,j,k-1)
					rijcxz=rho_C(5)%v(i,j,k-1)		  
!         z element
          rajczx=rho_C(5)%v(i-1,j,k-1) !xz
					rajczy=rho_C(6)%v(i-1,j,k-1) !yz
					rijczx=rho_C(5)%v(i,j,k-1)
					rijczy=rho_C(6)%v(i,j,k-1)
!
          rijcx=(rho_C(1)%v(i,j,k-1) * cx       &
                +rho_C(1)%v(i-1,j,k-1) * cxm1)  &
                /(cx+cxm1)
	 
          if (k.ne.n) then

						cz=mGrid%dz(k)

!           x element
						rajkxy=rho_C(4)%v(i-1,j,k) !xy
						rajkxz=rho_C(5)%v(i-1,j,k) !xz
						rijkxy=rho_C(4)%v(i,j,k)
						rijkxz=rho_C(5)%v(i,j,k)
!           z element
            rajkzx=rho_C(5)%v(i-1,j,k) !xz
						rajkzy=rho_C(6)%v(i-1,j,k) !yz
            rijkzx=rho_C(5)%v(i,j,k)
						rijkzy=rho_C(6)%v(i,j,k)
!
            rijkx=(rho_C(1)%v(i,j,k) * cx      &
                +rho_C(1)%v(i-1,j,k) * cxm1)   &
                /(cx+cxm1)
!
            rajkz=(rho_C(3)%v(i-1,j,k) * cz       &
                  +rho_C(3)%v(i-1,j,k-1) * czm1)  &
                  /(czm1+cz)
            rijkz=(rho_C(3)%v(i,j,k) * cz       &
                  +rho_C(3)%v(i,j,k-1) * czm1)  &
                  /(czm1+cz)

          end if
        
          jj=jj+1	 		  

!**********************************CsDone************************************
!------
! Hyijk   C17&***B12***
          if(k.ne.n) then
            spa(jj)=ciuw*dx*dy*dz+(rajkxz-rijkxz+rijcxz-rajcxz)*dy/2  &
                    +(rijkx/cz+rijcx/czm1)*dx*dy+(rijkz/cx+rajkz/cxm1)*dy*dz
          else
!++++
            do ia=1,6
              sg(ia)=sigma_C(ia)%v(i,j,k-1)
						end do
						call zhalf(sg,zh)
            spa(jj)=ciuw*dx*dy*dz+rijcx*dx*dy/czm1+zh(1,2)*dx*dy  &
                    +(rijcxz-rajcxz)*dy/2

          end if      
! Hydjk   C18
          if((i.ne.l).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-2)+i
            spa(kk)=-rijkz*dy*dz/cx
            ija(kk)=ii
          end if
!------
! Hyajf   C19
          if((i.ne.2).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-1)+i-2
            spa(kk)=-rajkxz*dy/2
            ija(kk)=ii
          end if 
! Hyijf   C20
          if(k.ne.n) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-1)+i-1
            spa(kk)=-rijkx*dx*dy/cz
            ija(kk)=ii
          end if
! Hydjf   C21
          if((i.ne.l).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-1)+i
            spa(kk)=rijkxz*dy/2
            ija(kk)=ii
          end if  
!*********************************CsDone*************************************  
!------ 
! Hzajc   C22&B13
          if((i.ne.2).and.(j.ne.1)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-2)+i-2
            spa(kk)=(rajczy*dy*czm1/cxm1-rajczx*czm1-rajcxy*dy)/4
            ija(kk)=ii
          end if
! Hzijc   C23&B14
          if(j.ne.1) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-2)+i-1
            spa(kk)=((rajcxy-rijcxy)*dy+(rijczx-rajczx)*czm1  &
                    -(rijczy/cx+rajczy/cxm1)*dy*czm1)/4+rijcx*dx
            ija(kk)=ii
          end if 
! Hzdjc   C24&B15
          if((i.ne.l).and.(j.ne.1)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-2)+i
            spa(kk)=(rijczy*dy*czm1/cx+rijczx*czm1+rijcxy*dy)/4
            ija(kk)=ii
          end if
!------
! Hzaec   C25&B16
          if((i.ne.2).and.(j.ne.m)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-1)+i-2
            spa(kk)=(rajczy*dy*czm1/cxm1+rajczx*czm1-rajcxy*dy)/4
            ija(kk)=ii
          end if       
! Hziec   C26&B17
          if(j.ne.m) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-1)+i-1
            spa(kk)=((rajcxy-rijcxy)*dy-(rijczx-rajczx)*czm1  &
                    -(rijczy/cx+rajczy/cxm1)*dy*czm1)/4-rijcx*dx
            ija(kk)=ii
          end if 
! Hzdec   C27&B18
          if((i.ne.l).and.(j.ne.m)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-2)+(l-1)*(j-1)+i
            spa(kk)=(rijczy*dy*czm1/cx-rijczx*czm1+rijcxy*dy)/4
            ija(kk)=ii
          end if 
!------
! Hzajk   C28
          if((i.ne.2).and.(j.ne.1).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-2)+i-2
            spa(kk)=(rajkzy*dy*cz/cxm1-rajkzx*cz+rajkxy*dy)/4
            ija(kk)=ii
          end if   
! Hzijk   C29
          if((j.ne.1).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-2)+i-1
            spa(kk)=((rijkxy-rajkxy)*dy+(rijkzx-rajkzx)*cz  &
                    -(rijkzy/cx+rajkzy/cxm1)*dy*cz)/4-rijkx*dx
            ija(kk)=ii
          end if 
! Hzdjk   C30
          if((i.ne.l).and.(j.ne.1).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-2)+i
            spa(kk)=(rijkzy*dy*cz/cx+rijkzx*cz-rijkxy*dy)/4
            ija(kk)=ii
          end if
!------  
! Hzaek   C31
          if((i.ne.2).and.(j.ne.m).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i-2
            spa(kk)=(rajkzy*dy*cz/cxm1+rajkzx*cz+rajkxy*dy)/4
            ija(kk)=ii
          end if      
! Hziek   C32
          if((j.ne.m).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i-1
            spa(kk)=((rijkxy-rajkxy)*dy-(rijkzx-rajkzx)*cz  &
                    -(rijkzy/cx+rajkzy/cxm1)*dy*cz)/4+rijkx*dx
            ija(kk)=ii
          end if     
! Hzdek   C33
          if((i.ne.l).and.(j.ne.m).and.(k.ne.n)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i
            spa(kk)=(rijkzy*dy*cz/cx-rijkzx*cz-rijkxy*dy)/4
            ija(kk)=ii
          end if 

          ija(jj+1)=kk+1

				end do
			end do
		end do

!   debug
!   do i=ikk,kk
!      write(901,*) spa(i)
!      write(902,*) ija(i)
!   end do
!   pause
!    ikk=kk ! debug

!--------------------------------------------------
!   Third, Hz Form
!--------------------------------------------------
    do k=1,n-1
      
      dz=mGrid%dz(k)
      
      do j=2,m
        
        dy=mGrid%delY(j)
        cym1=mGrid%dy(j-1)
			  cy=mGrid%dy(j)
			  
        do i=2,l
          
          dx=mGrid%delX(i)
					cxm1=mGrid%dx(i-1)
					cx=mGrid%dx(i)
					
!         x element
          rabkxy=rho_C(4)%v(i-1,j-1,k)
					rajkxy=rho_C(4)%v(i-1,j,k)
!         y element
					ribkxy=rho_C(4)%v(i,j-1,k)
					rijkxy=rho_C(4)%v(i,j,k)

          ribkx=(rho_C(1)%v(i,j-1,k) * cx       &
                +rho_C(1)%v(i-1,j-1,k) * cxm1)  &
                /(cx+cxm1)
					rijkx=(rho_C(1)%v(i,j,k) * cx       &
                +rho_C(1)%v(i-1,j,k) * cxm1)  &
                /(cx+cxm1)

          rajky=(rho_C(2)%v(i-1,j,k) * cy       &
                +rho_C(2)%v(i-1,j-1,k) * cym1)  &
                /(cym1+cy)
					rijky=(rho_C(2)%v(i,j,k) * cy       &
                +rho_C(2)%v(i,j-1,k) * cym1)  &
                /(cym1+cy)
        
          jj=jj+1

!*********************************CsDone*************************************
!------
! Hzijk   C29
          spa(jj)=ciuw*dx*dy*dz+(rajkxy-rijkxy+ribkxy-rabkxy)*dz/2  &
                  +(rijkx/cy+ribkx/cym1)*dx*dz+(rijky/cx+rajky/cxm1)*dy*dz        
! Hzdjk   C30
          if(i.ne.l) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-2)+i
            spa(kk)=-rijky*dy*dz/cx
            ija(kk)=ii
          end if  
!------
! Hzaek   C31
          if((i.ne.2).and.(j.ne.m)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i-2
            spa(kk)=-rajkxy*dz/2
            ija(kk)=ii
          end if     
! Hziek   C32
          if(j.ne.m) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i-1
            spa(kk)=-rijkx*dx*dz/cy
            ija(kk)=ii
          end if     
! Hzdek   C33
          if((i.ne.l).and.(j.ne.m)) then
            kk=kk+1
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i
            spa(kk)=rijkxy*dz/2
            ija(kk)=ii
          end if 
             
          ija(jj+1)=kk+1  

				end do
			end do
		end do

!   debug
!   do i=1,kk
!      write(901,*) i,spa(i)
!      write(902,*) i,ija(i)
!   end do
!   pause

  end subroutine CurlcurleSetUp    ! CurlcurleSetUp


  !************ preconditioning operator *******************
  subroutine Mackie_ILU()
    implicit none
    integer                     :: status     ! for dynamic memory allocation
    integer                     :: is,i,j,k,l,m,n,nzn ! dummy variables
    integer                     :: nhx,nhy,nhz,np,np2,np5
		integer                     :: ii,jj,nnz

!   for Hx>>
    real(kind=prec) :: ribcyx,ribcyz,ribkyx,ribkyz,rijcyx,rijcyz,rijkyx,rijkyz !total 8
		real(kind=prec) :: ribczx,ribczy,ribkzx,ribkzy,rijczx,rijczy,rijkzx,rijkzy 
		real(kind=prec) :: rijcy,rijky
		real(kind=prec) :: ribkz,rijkz
!   for Hy>>
		real(kind=prec) :: rajcxy,rajcxz,rajkxy,rajkxz,rijcxy,rijcxz,rijkxy,rijkxz
		real(kind=prec) :: rajczx,rajczy,rajkzx,rajkzy !,rijczx,rijczy,rijkzx,rijkzy
		real(kind=prec) :: rijcx,rijkx
		real(kind=prec) :: rajkz !,rijkz
!   for Hz>>
		real(kind=prec) :: rabkxy,ribkxy !,rajkxy,rajkxz,rabkxz,ribkxz,rijkxy,rijkxz
!   real(kind=prec) :: rabkyx,rabkyz,rajkyx,rajkyz,ribkyx,ribkyz,rijkyx,rijkyz
		real(kind=prec) :: ribkx !,rijkx
		real(kind=prec) :: rajky !,rijky

		real(kind=prec)    :: sg(6)
		real(kind=prec)    :: dx,dy,dz,cxm1,cx,cym1,cy,czm1,cz
		complex(kind=prec) :: zh(2,2),ciuw,fact
! 
		integer,allocatable :: ilu(:)
    integer,allocatable :: levs(:),ia(:),ja(:),jw(:)
    complex(kind=prec),allocatable :: wka(:),spm(:)
    integer lfil,ierr,iwk

    ciuw = ISIGN * CMPLX(0.0d0, omega*MU_0)

    l = mGrid%Nx
		m = mGrid%Ny
		n = mGrid%Nz
		nzn=n
		n=n+1

    nhx=l*(m-1)*nzn
    nhy=(l-1)*m*nzn
    nhz=(l-1)*(m-1)*nzn
    np=nhx+nhy+nhz

    np2=(l+1)*(m+1)*n*3

    np5=l*(m-1)*nzn+(l-1)*m*nzn+3*l*(m-2)*(nzn-1)+l*(m-1)*(nzn-1)+  &
        12*(l-1)*(m-1)*(nzn-1)+9*(l-1)*(m-1)*nzn+                   &
				3*(l-2)*m*(nzn-1)+(l-1)*m*(nzn-1)+5*(l-2)*(m-1)*nzn+        &
        5*(l-1)*(m-2)*nzn+4*(l-1)*(m-2)*(nzn-1)+                    &
        4*(l-2)*(m-1)*(nzn-1)+2*(l-2)*(m-2)*nzn+1

    allocate(ilu(np5),levs(np5*10),ia(np2+1),ja(np5),jw(3*np2)) !integer
    allocate(wka(np2),spm(np5))	!complex
    allocate(alu(np5*10),jlu(np5*10),ju(np2+1)) ! common variables
		 
		ilu=0
		levs=0
		ia=0
		ja=0
		jw=0
		wka=C_ZERO
		spm=C_ZERO

		jlu=0
		ju=0
		alu=C_ZERO

!--------------------------------------------------
!   First, Mxx
!--------------------------------------------------
    jj=0
    nnz=0
    
    do i=1,l
      
      dx=mGrid%dx(i)
      
      do k=2,n
        
        dz=mGrid%delZ(k)
        czm1=mGrid%dz(k-1)
        
        do j=2,m
        
          dy=mGrid%delY(j)
          cym1=mGrid%dy(j-1)
          cy=mGrid%dy(j)
          
!         y element           
          ribcyx=rho_C(4)%v(i,j-1,k-1) !xy
					ribcyz=rho_C(6)%v(i,j-1,k-1) !yz
					rijcyx=rho_C(4)%v(i,j,k-1)   !xy
					rijcyz=rho_C(6)%v(i,j,k-1)   !yz
!         z element
          ribczx=rho_C(5)%v(i,j-1,k-1) !xz
					ribczy=rho_C(6)%v(i,j-1,k-1) !yz
					rijczx=rho_C(5)%v(i,j,k-1)   !xz
					rijczy=rho_C(6)%v(i,j,k-1)   !yz
! 
          rijcy=(rho_C(2)%v(i,j,k-1) * cy       &
                +rho_C(2)%v(i,j-1,k-1) * cym1)  &
                /(cym1+cy)

          if (k.ne.n) then

						cz=mGrid%dz(k)

!           y element
						ribkyx=rho_C(4)%v(i,j-1,k) !xy
						ribkyz=rho_C(6)%v(i,j-1,k) !yz
						rijkyx=rho_C(4)%v(i,j,k)   !xy
						rijkyz=rho_C(6)%v(i,j,k)   !yz
!           z element
						ribkzx=rho_C(5)%v(i,j-1,k) !xz
						ribkzy=rho_C(6)%v(i,j-1,k) !yz
						rijkzx=rho_C(5)%v(i,j,k)   !xz
						rijkzy=rho_C(6)%v(i,j,k)   !yz
!
            rijky=(rho_C(2)%v(i,j,k) * cy       &
                  +rho_C(2)%v(i,j-1,k) * cym1)  &
                  /(cym1+cy)
!
            ribkz=(rho_C(3)%v(i,j-1,k) * cz       &
                  +rho_C(3)%v(i,j-1,k-1) * czm1)  &
                  /(czm1+cz)
            rijkz=(rho_C(3)%v(i,j,k) * cz       &
                  +rho_C(3)%v(i,j,k-1) * czm1)  &
                  /(czm1+cz)

          end if

        
          jj=jj+1
!*****************************CsDone*****************************************
!------
! Hxijk   C5&***B4***
          if(k.ne.n) then
            fact=ciuw*dx*dy*dz+(rijcyz-rijkyz+ribkyz-ribcyz)*dx/2   &
                    +(rijky/cz+rijcy/czm1)*dx*dy+(rijkz/cy+ribkz/cym1)*dx*dz
					else
!++++
            do is=1,6
              sg(is)=sigma_C(is)%v(i,j,k-1)
						end do
            call zhalf(sg,zh)		                  
            fact=ciuw*dx*dy*dz+rijcy*dx*dy/czm1-zh(2,1)*dx*dy   &
                    +(rijcyz-ribcyz)*dx/2
          end if
          nnz=nnz+1
          alu(nnz)=fact
          jlu(nnz)=jj
          ilu(nnz)=jj
! Hxiek   C6
          if((j.ne.m).and.(k.ne.n)) then
            ii=(m-1)*(n-1)*(i-1)+(m-1)*(k-2)+j
            fact=-rijkz*dx*dz/cy
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if
!------
! Hxibf   C7
          if((j.ne.2).and.(k.ne.n)) then
            ii=(m-1)*(n-1)*(i-1)+(m-1)*(k-1)+j-2
            fact=-ribkyz*dx/2
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if
! Hxijf   C8
          if(k.ne.n) then
            ii=(m-1)*(n-1)*(i-1)+(m-1)*(k-1)+j-1
            fact=-rijky*dx*dy/cz
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if
! Hxief   C9
          if((j.ne.m).and.(k.ne.n)) then
            ii=(m-1)*(n-1)*(i-1)+(m-1)*(k-1)+j
            fact=rijkyz*dx/2
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if        

				end do
			end do
		end do

!--------------------------------------------------
!   Second, Myy
!--------------------------------------------------
    do j=1,m
      
      dy=mGrid%dy(j)
      
      do k=2,n
        
        dz=mGrid%delZ(k)
        czm1=mGrid%dz(k-1)
        
        do i=2,l
          
          dx=mGrid%delX(i)
          cxm1=mGrid%dx(i-1)
          cx=mGrid%dx(i)
		  
!         x element
          rajcxy=rho_C(4)%v(i-1,j,k-1) !xy
					rajcxz=rho_C(5)%v(i-1,j,k-1) !xz
					rijcxy=rho_C(4)%v(i,j,k-1)
					rijcxz=rho_C(5)%v(i,j,k-1)		  
!         z element
          rajczx=rho_C(5)%v(i-1,j,k-1) !xz
					rajczy=rho_C(6)%v(i-1,j,k-1) !yz
					rijczx=rho_C(5)%v(i,j,k-1)
					rijczy=rho_C(6)%v(i,j,k-1)
!
          rijcx=(rho_C(1)%v(i,j,k-1) * cx       &
                +rho_C(1)%v(i-1,j,k-1) * cxm1)  &
                /(cx+cxm1)
	 
          if (k.ne.n) then

						cz=mGrid%dz(k)

!           x element
						rajkxy=rho_C(4)%v(i-1,j,k) !xy
						rajkxz=rho_C(5)%v(i-1,j,k) !xz
						rijkxy=rho_C(4)%v(i,j,k)
						rijkxz=rho_C(5)%v(i,j,k)
!           z element
            rajkzx=rho_C(5)%v(i-1,j,k) !xz
						rajkzy=rho_C(6)%v(i-1,j,k) !yz
            rijkzx=rho_C(5)%v(i,j,k)
						rijkzy=rho_C(6)%v(i,j,k)
!
            rijkx=(rho_C(1)%v(i,j,k) * cx      &
                +rho_C(1)%v(i-1,j,k) * cxm1)   &
                /(cx+cxm1)
!
            rajkz=(rho_C(3)%v(i-1,j,k) * cz       &
                  +rho_C(3)%v(i-1,j,k-1) * czm1)  &
                  /(czm1+cz)
            rijkz=(rho_C(3)%v(i,j,k) * cz       &
                  +rho_C(3)%v(i,j,k-1) * czm1)  &
                  /(czm1+cz)

          end if
        
          jj=jj+1	 		  

!**********************************CsDone************************************
!------
! Hyijk   C17&***B12***
          if(k.ne.n) then
            fact=ciuw*dx*dy*dz+(rajkxz-rijkxz+rijcxz-rajcxz)*dy/2  &
                    +(rijkx/cz+rijcx/czm1)*dx*dy+(rijkz/cx+rajkz/cxm1)*dy*dz
          else
!++++
            do is=1,6
              sg(is)=sigma_C(is)%v(i,j,k-1)
						end do
						call zhalf(sg,zh)
            fact=ciuw*dx*dy*dz+rijcx*dx*dy/czm1+zh(1,2)*dx*dy  &
                    +(rijcxz-rajcxz)*dy/2

          end if  
          nnz=nnz+1
          alu(nnz)=fact
          jlu(nnz)=jj
          ilu(nnz)=jj	    
! Hydjk   C18
          if((i.ne.l).and.(k.ne.n)) then
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-2)+i
            fact=-rijkz*dy*dz/cx
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if
!------
! Hyajf   C19
          if((i.ne.2).and.(k.ne.n)) then
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-1)+i-2
            fact=-rajkxz*dy/2
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if 
! Hyijf   C20
          if(k.ne.n) then
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-1)+i-1
            fact=-rijkx*dx*dy/cz
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if
! Hydjf   C21
          if((i.ne.l).and.(k.ne.n)) then
            ii=nhx+(l-1)*(n-1)*(j-1)+(l-1)*(k-1)+i
            fact=rijkxz*dy/2
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if

				end do
			end do
		end do

!--------------------------------------------------
!   Third, Mzz
!--------------------------------------------------
    do k=1,n-1
      
      dz=mGrid%dz(k)
      
      do j=2,m
        
        dy=mGrid%delY(j)
        cym1=mGrid%dy(j-1)
				cy=mGrid%dy(j)
				
        do i=2,l
          
          dx=mGrid%delX(i)
          cxm1=mGrid%dx(i-1)
          cx=mGrid%dx(i)
		  
!         x element
          rabkxy=rho_C(4)%v(i-1,j-1,k)
					rajkxy=rho_C(4)%v(i-1,j,k)
!         y element
					ribkxy=rho_C(4)%v(i,j-1,k)
					rijkxy=rho_C(4)%v(i,j,k)

          ribkx=(rho_C(1)%v(i,j-1,k) * cx       &
                +rho_C(1)%v(i-1,j-1,k) * cxm1)  &
                /(cx+cxm1)
					rijkx=(rho_C(1)%v(i,j,k) * cx       &
                +rho_C(1)%v(i-1,j,k) * cxm1)  &
                /(cx+cxm1)

          rajky=(rho_C(2)%v(i-1,j,k) * cy       &
                +rho_C(2)%v(i-1,j-1,k) * cym1)  &
                /(cym1+cy)
					rijky=(rho_C(2)%v(i,j,k) * cy       &
                +rho_C(2)%v(i,j-1,k) * cym1)  &
                /(cym1+cy)
        
          jj=jj+1

!*********************************CsDone*************************************
!------
! Hzijk   C29
          fact=ciuw*dx*dy*dz+(rajkxy-rijkxy+ribkxy-rabkxy)*dz/2  &
                  +(rijkx/cy+ribkx/cym1)*dx*dz+(rijky/cx+rajky/cxm1)*dy*dz        
          nnz=nnz+1
          alu(nnz)=fact
          jlu(nnz)=jj
          ilu(nnz)=jj
! Hzdjk   C30
          if(i.ne.l) then
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-2)+i
            fact=-rijky*dy*dz/cx
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if  
!------
! Hzaek   C31
          if((i.ne.2).and.(j.ne.m)) then
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i-2
            fact=-rajkxy*dz/2
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if     
! Hziek   C32
          if(j.ne.m) then
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i-1
            fact=-rijkx*dx*dz/cy
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if     
! Hzdek   C33
          if((i.ne.l).and.(j.ne.m)) then
            ii=nhx+nhy+(l-1)*(m-1)*(k-1)+(l-1)*(j-1)+i
            fact=rijkxy*dz/2
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=ii
            ilu(nnz)=jj
            nnz=nnz+1
            alu(nnz)=fact
            jlu(nnz)=jj
            ilu(nnz)=ii
          end if

				end do
			end do
		end do

! debug
!  do i=1,nnz
!     write(101,*) alu(i)
!	 write(102,*) jlu(i)
!  end do
!  pause

    !  
    call coocsr(np,nnz,alu,ilu,jlu,spm,ja,ia)
    !
		lfil=16
    iwk=np5*10

		call iluk(np,spm,ja,ia,lfil,alu,jlu,ju,levs,iwk,wka,jw,ierr)
    print*,'ierr ',ierr

! debug
!  do i=1,nnz
!     write(101,*) alu(i)
!	 write(102,*) jlu(i)
!  end do
!  pause

    deallocate(ilu,levs,ia,ja,jw)
    deallocate(wka,spm)

  end subroutine Mackie_ILU




!    ===================================================================
      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
!----------------------------------------------------------------------- 
      implicit none
	    integer i,j,k,k0,iad,nrow,nnz
      complex(kind=prec) :: a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
!-----------------------------------------------------------------------
!  Coordinate     to   Compressed Sparse Row 
!----------------------------------------------------------------------- 
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry:
!--------- 
! nrow	= dimension of the matrix 
! nnz	= number of nonzero elements in matrix
! a,
! ir, 
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
! 	  the elements, ir(k) = its row number and jc(k) = its column 
!	  number. The order of the elements is arbitrary. 
!
! on return:
!----------- 
! ir 	is destroyed
!
! ao, jao, iao = matrix in general sparse matrix format with ao 
! 	continung the real values, jao containing the column indices, 
!	and iao being the pointer to the begi0nning of the row, 
!	in arrays ao, jao.
!
! Notes:
!------ This routine is NOT in place.  See coicsr
!
!------------------------------------------------------------------------
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
! determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
! starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
! go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
! shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
!------------- end of coocsr ------------------------------------------- 
!----------------------------------------------------------------------- 
      end subroutine coocsr


!***********************************************************************
      subroutine iluk(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ierr)
      implicit none 
      integer n
      complex(kind=prec) :: a(*),alu(*),w(n)
      integer ja(*),ia(n+1),jlu(*),ju(n),levs(*),jw(3*n),lfil,iwk,ierr
!----------------------------------------------------------------------* 
!     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) *
!----------------------------------------------------------------------*
!
! on entry:
!========== 
! n       = integer. The row dimension of the matrix A. The matrix 
!
! a,ja,ia = matrix stored in Compressed Sparse Row format.              
!
! lfil    = integer. The fill-in parameter. Each element whose
!           leve-of-fill exceeds lfil during the ILU process is dropped.
!           lfil must be .ge. 0 
!
! tol     = real*8. Sets the threshold for dropping small terms in the
!           factorization. See below for details on dropping strategy.
!  
! iwk     = integer. The minimum length of arrays alu, jlu, and levs.
!
! On return:
!===========
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
!
! levs    = integer (work) array of size iwk -- which contains the 
!           levels of each element in alu, jlu.
!
! ierr    = integer. Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is .gt.  n.)
!           ierr  = -2   --> The matrix L overflows the array al.
!           ierr  = -3   --> The matrix U overflows the array alu.
!           ierr  = -4   --> Illegal value for lfil.
!           ierr  = -5   --> zero row encountered in A or U.
!
! work arrays:
!=============
! jw      = integer work array of length 3*n.
! w       = real work array of length n 
!
! Notes/known bugs: This is not implemented efficiently storage-wise.
!       For example: Only the part of the array levs(*) associated with
!       the U-matrix is needed in the routine.. So some storage can 
!       be saved if needed. The levels of fills in the LU matrix are
!       output for information only -- they are not needed by LU-solve. 
!        
!----------------------------------------------------------------------
! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
! jw(n+1:2n)  stores the nonzero indicator. 
! 
! Notes:
! ------
! All the diagonal elements of the input matrix must be  nonzero.
!
!----------------------------------------------------------------------* 
!     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,n2,jlev, min 
      complex(kind=prec) :: t, s, fact 
      if (lfil .lt. 0) goto 998
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)
!     and pointer array.
!-----------------------------------------------------------------------
      n2 = n+n 
      ju0 = n+2
      jlu(1) = ju0
!
!     initialize nonzero indicator array + levs array -- 
!
      do 1 j=1,2*n 
         jw(j)  = 0
 1    continue
!-----------------------------------------------------------------------
!     beginning of main loop.
!-----------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
!     
!     unpack L-part and U-part of row of A in arrays w 
!     
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
!
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
            if (cdabs(t) .eq. 0.0) goto 170 
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n2+lenl) = 0 
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
               jw(n2+ii) = 0 
            else
               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k
               w(jpos) = t
               jw(n2+jpos) = 0 
               jw(n+k) = jpos
            endif
 170     continue
!
         jj = 0
!
!     eliminate previous rows
!     
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.
!-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
!     
!     determine smallest column index
!     
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
!
         if (k .ne. jj) then
!     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
!     exchange in jw(n+  (pointers/ nonzero indicator).
            jw(n+jrow) = jj
            jw(n+j) = k
!     exchange in jw(n2+  (levels) 
            j = jw(n2+jj) 
            jw(n2+jj)  = jw(n2+k) 
            jw(n2+k) = j
!     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
!
!     zero out element in row by resetting jw(n+jrow) to zero.
!     
         jw(n+jrow) = 0
!     
!     get the multiplier for row to be eliminated (jrow) + its level
!     
         fact = w(jj)*alu(jrow)
         jlev = jw(n2+jj) 
         if (jlev .gt. lfil) goto 150
!
!     combine current row and row jrow
!
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
!     
!     dealing with upper part.
!     
               if (jpos .eq. 0) then
!
!     this is a fill-in element
!     
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
                  jw(n2+i) = jlev+levs(k)+1 
               else
!
!     this is not a fill-in element 
!
                  w(jpos) = w(jpos) - s
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
               endif
            else
!     
!     dealing with lower part.
!     
               if (jpos .eq. 0) then
!
!     this is a fill-in element
!
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
                  jw(n2+lenl) = jlev+levs(k)+1 
               else
!
!     this is not a fill-in element 
!
                  w(jpos) = w(jpos) - s
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
               endif
            endif
 203     continue
         w(jj) = fact
         jw(jj)  = jrow
         goto 150 
 160     continue 
!     
!     reset double-pointer to zero (U-part) 
!     
         do 308 k=1, lenu
            jw(n+jw(ii+k-1)) = 0
 308     continue
!
!     update l-matrix
!         
         do 204 k=1, lenl 
!------------------------------------------------------------------------
            if (ju0 .gt. iwk) goto 996
!------------------------------------------------------------------------
            if (jw(n2+k) .le. lfil) then
               alu(ju0) =  w(k)
               jlu(ju0) =  jw(k)
               ju0 = ju0+1
            endif
 204     continue
!     
!     save pointer to beginning of row ii of U
!     
         ju(ii) = ju0
!
!     update u-matrix
!
         do 302 k=ii+1,ii+lenu-1 
            if (jw(n2+k) .le. lfil) then
               jlu(ju0) = jw(k)
               alu(ju0) = w(k)
               levs(ju0) = jw(n2+k) 
               ju0 = ju0+1
            endif
 302     continue

         if (w(ii) .eq. 0.0) goto 999 
!     
         alu(ii) = cmplx(1.0d0,0.d0)/ w(ii) 
!     
!     update pointer to beginning of next row of U.
!     
         jlu(ii+1) = ju0
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
 500  continue
      ierr = 0
      return
!
!     incomprehensible error. Matrix must be wrong.
!     
 995  ierr = -1
      return
!     
!     insufficient storage in L.
!     
 996  ierr = -2
      return
!     
!     insufficient storage in U.
!     
 997  ierr = -3
      return
!     
!     illegal lfil entered.
!     
 998  ierr = -4
      return
!     
!     zero row encountered in A or U. 
!     
 999  ierr = -5
      return
!----------------end-of-iluk--------------------------------------------
!-----------------------------------------------------------------------
      end subroutine iluk
      

   !**********************************************************************
   ! SetUp routines do calculations (maybe once; possibly more than once)
   ! DivCorrSetup must be called once for each conductivity distribuition
   !  (i.e., before first forward run; after any change to conductivity)
  subroutine DivCorrSetUp()
    implicit none

    integer :: l,m,n,i,j,k,ii,jj,kk
	  integer :: np,np6,md
	  real (kind=prec) :: dxi,dxa,dyj,dyb,dzk,dzc,dx,dy,dz,dx1,dy1,dz1

    l = mGrid%Nx
	  m = mGrid%Ny
	  n = mGrid%Nz

    np=(l-1)*(m-1)*(n-1)

    np6=(l-1)*(m-1)*(n-1)+(l-2)*(m-1)*(n-1)+  &
        (l-1)*(m-2)*(n-1)+(l-1)*(m-1)*(n-2)+1

	  allocate(sppot(np6),sppot1(np6),ijapot(np6))

    sppot=R_ZERO
	  sppot1=R_ZERO
	  ijapot=0

    jj=0
    ijapot(1)=np+2
    kk=np+1
    
	do k=2,n
      do j=2,m
        do i=2,l
          dx=mGrid%delXinv(i)
          dy=mGrid%delYinv(j)
          dz=mGrid%delZinv(k)
          dx1=1.0/dx
          dy1=1.0/dy
          dz1=1.0/dz
          dxi=mGrid%dx(i)
          dxa=mGrid%dx(i-1)
          dyj=mGrid%dy(j)
          dyb=mGrid%dy(j-1)
          dzc=mGrid%dz(k-1)
          dzk=mGrid%dz(k)

          jj=jj+1

!
! Phi(i,j,k)
!
          sppot(jj)=dy1*dz1/dxi+dy1*dz1/dxa+dx1*dz1/dyj+  &
                    dx1*dz1/dyb+dx1*dy1/dzk+dx1*dy1/dzc 
!
! Phi(i+1,j,k)
!
          if(i.ne.l) then
            kk=kk+1
            ii=(l-1)*(m-1)*(k-2)+(l-1)*(j-2)+i
            sppot(kk)=-dy1*dz1/dxi
            ijapot(kk)=ii
          end if
!
! Phi(i,j+1,k)
!
          if(j.ne.m) then
            kk=kk+1
            ii=(l-1)*(m-1)*(k-2)+(l-1)*(j-1)+i-1
            sppot(kk)=-dx1*dz1/dyj
            ijapot(kk)=ii
          end if
!
! Phi(i,j,k+1)
!
          if(k.ne.n) then
            kk=kk+1
            ii=(l-1)*(m-1)*(k-1)+(l-1)*(j-2)+i-1
            sppot(kk)=-dx1*dy1/dzk
            ijapot(kk)=ii
          end if

          ijapot(jj+1)=kk+1

        end do
	  end do
	end do

    do i=1,kk
      sppot1(i)=sppot(i)
!	  write(107,*) i,sppot1(i) !debug
!	  write(107,*) i,ijapot(i) !debug
    end do
!	pause !debug

! Incomplete Cholesky Decomposition
    do k=1,np
      sppot1(k)=dsqrt(sppot1(k))
      do j=ijapot(k), ijapot(k+1)-1
        sppot1(j)=sppot1(j)/sppot1(k)
      end do
 
      md=ijapot(k+1)-ijapot(k)
      ii=ijapot(k)
      if(md.ge.1) sppot1(ijapot(ii))=sppot1(ijapot(ii))-sppot1(ii)**2
      if(md.ge.2) sppot1(ijapot(ii+1))=sppot1(ijapot(ii+1))-sppot1(ii+1)**2
      if(md.eq.3) sppot1(ijapot(ii+2))=sppot1(ijapot(ii+2))-sppot1(ii+2)**2
 
    end do

!    write(109,*) sppot1 !debug
!	pause !debug

  end subroutine DivCorrSetUp


  ! creat operaters for forward solver and preconditioner as well as
  ! those for magnetic field divergence correction  
  subroutine modelOperatorSetup()
    
    call CurlcurleSetUp()
    call Mackie_ILU()
    call DivCorrSetUp()
  
  end subroutine modelOperatorSetup

  ! clean up 3D operators
  subroutine modelOperatorCleanup()
    
    deallocate(spa,ija)
    deallocate(alu,jlu,ju)
    deallocate(sppot,sppot1,ijapot)
  
  end subroutine modelOperatorCleanup


  !*********************
  subroutine MInv_multB(inE, outE)
    implicit none
		type (cvector), intent(in)     ::  inE
		type (cvector), intent(inout)  ::  outE
	! local variables
    integer                     :: i,k,l,m,n,nzn ! dummy variables
    integer                     :: nhx,nhy,nhz,np
	  logical                     :: flag
	  complex (kind=prec), pointer, dimension(:)    :: x,y

    l = inE%Nx
	  m = inE%Ny
	  n = inE%Nz
	  nzn=n
	  n=n+1

    nhx=l*(m-1)*nzn
    nhy=(l-1)*m*nzn
    nhz=(l-1)*(m-1)*nzn
    np=nhx+nhy+nhz
	
	  allocate(x(np),y(np))

    x=C_ZERO
	  y=C_ZERO

    if (.not.inE%allocated) then
      write(0,*) 'inE in MultA_N not allocated yet'
      stop
    end if

    if (.not.outE%allocated) then
      write(0,*) 'outE in MultA_N not allocated yet'
      stop
    end if

    ! Check whether the bounds are the same
    if ((inE%nx == outE%nx).and.&
         (inE%ny == outE%ny).and.&
         (inE%nz == outE%nz)) then

       if ((inE%gridType == outE%gridType)) then

          flag=.true.
	      call copy_vs(inE,flag,y)

 	      call lusol(np,y,x)

	      flag=.false.
	      call copy_vs(outE,flag,x)

       else
          write (0, *) 'MultA_N: not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-complex vectors for MultA_N are not of same size'
    end if

		deallocate(x,y)

  end subroutine MInv_multB

   
  complex*16 function cddot(n,x)
    implicit none
		integer n,i
		complex*16 x(*)
		cddot=C_ZERO
		do i=1,n
			cddot=cddot+conjg(x(i))*x(i)
		end do
		cddot=cdsqrt(cddot)

  end function cddot


  !**************************************
  subroutine lusol(n, y, x)
    implicit none
		integer n
    complex (kind=prec) :: x(n), y(n)
	
!-----------------------------------------------------------------------
!
! This routine solves the system (LU) x = y, 
! given an LU decomposition of a matrix stored in (alu, jlu, ju) 
! modified sparse row format 
!
!-----------------------------------------------------------------------
! on entry:
! n   = dimension of system 
! y   = the right-hand-side vector
! alu, jlu, ju 
!     = the LU matrix as provided from the ILU routines. 
!
! on return
! x   = solution of LU x = y.     
!-----------------------------------------------------------------------
! 
! Note: routine is in place: call lusol (n, x, x, alu, jlu, ju) 
!       will solve the system with rhs x and overwrite the result on x . 
!
!-----------------------------------------------------------------------
! local variables
!
        integer i,k
!
! forward solve
!
        do 40 i = 1, n
           x(i) = y(i)
           do 41 k=jlu(i),ju(i)-1
              x(i) = x(i) - alu(k)* x(jlu(k))
 41        continue
 40     continue
!
!     backward solve.
!
	do 90 i = n, 1, -1
	   do 91 k=ju(i),jlu(i+1)-1
              x(i) = x(i) - alu(k)*x(jlu(k))
 91	   continue
           x(i) = alu(i)*x(i)
 90     continue
!
  	return
!----------------end of lusol ------------------------------------------
  end subroutine lusol



  ! *************************************
  subroutine MultA_N(inE, adjt, outE)

    implicit none
    type (cvector), intent (in)              :: inE
    logical, intent (in)                     :: adjt
    type (cvector), intent (inout)           :: outE

    if (.not.inE%allocated) then
      write(0,*) 'inE in MultA_N_ANI not allocated yet'
      stop
    end if

    if (.not.outE%allocated) then
      write(0,*) 'outE in MultA_N_ANI not allocated yet'
      stop
    end if

    ! Check whether the bounds are the same
    if ((inE%nx == outE%nx).and.&
         (inE%ny == outE%ny).and.&
         (inE%nz == outE%nz)) then

       if ((inE%gridType == outE%gridType)) then

          Call atimes(inE, adjt, outE)

       else
          write (0, *) 'MultA_N: not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-complex vectors for MultA_N are not of same size'
    end if

  end subroutine MultA_N
  
  
  subroutine atimes(inE, adjt, outE)
    implicit none
    type (cvector), intent (in)              :: inE
    logical, intent (in)                     :: adjt
    type (cvector), intent (inout)           :: outE
		! local variables
    integer                     :: i,k,l,m,n,nzn ! dummy variables
    integer                     :: nhx,nhy,nhz,np
		logical                     :: flag
		complex (kind=prec), pointer, dimension(:)    :: r,x

    l = inE%Nx
		m = inE%Ny
		n = inE%Nz
		nzn=n
		n=n+1

    nhx=l*(m-1)*nzn
    nhy=(l-1)*m*nzn
    nhz=(l-1)*(m-1)*nzn
    np=nhx+nhy+nhz
	
		allocate(r(np),x(np))

    r=C_ZERO
		x=C_ZERO

		flag=.true.
		call copy_vs(inE,flag,x)	
	
    ! 预留接口
    do i=1, np
      r(i)=r(i)+spa(i)*x(i)
      do k=ija(i),ija(i+1)-1
        r(i)=r(i)+spa(k)*x(ija(k))
        r(ija(k))=r(ija(k))+spa(k)*x(i)
			end do
		end do

    flag=.false.
		call copy_vs(OutE,flag,r)
	
		deallocate(r,x)	
	  
  end subroutine atimes

  
  ! **************************************
  subroutine copy_vs(Ev,flag,Es)
    implicit none
	  logical,intent (in)      :: flag
    type (cvector)           :: Ev
    complex (kind=prec)      :: Es(*)
    ! local variables
	  integer ir,i,j,k,l,m,n

		l = mGrid%Nx
		m = mGrid%Ny
		n = mGrid%Nz + 1

		if(flag) then
			ir=1
			! Hx
      do i=1,l
        do k=2,n
          do j=2,m
            Es(ir)=Ev%x(i,j,k)
						ir=ir+1
					end do
				end do
			end do
			! Hy
			do j=1,m
				do k=2,n
					do i=2,l
            Es(ir)=Ev%y(i,j,k)
						ir=ir+1
					end do
				end do
			end do
			! Hz
      do k=1,n-1
        do j=2,m
          do i=2,l
            Es(ir)=Ev%z(i,j,k)
						ir=ir+1
					end do
				end do
			end do
		!
		else
			ir=1
			! Hx
      do i=1,l
        do k=2,n
          do j=2,m
            Ev%x(i,j,k)=Es(ir)
						ir=ir+1
					end do
				end do
			end do
			! Hy
			do j=1,m
        do k=2,n
          do i=2,l
            Ev%y(i,j,k)=Es(ir)
						ir=ir+1
					end do
				end do
			end do
			! Hz
      do k=1,n-1
        do j=2,m
          do i=2,l
            Ev%z(i,j,k)=Es(ir)
						ir=ir+1
					end do
				end do
			end do
		end if

  end subroutine copy_vs


  !************************************
  subroutine DivA(inPhi, outPhi)
    implicit none
    type (rscalar), intent(in)                :: inPhi
    type (rscalar), intent(inout)             :: outPhi
    !
		integer i,k,l,m,n,np
		logical flag
    real (kind=prec),allocatable :: r(:),x(:)

		l = mGrid%Nx
		m = mGrid%Ny
		n = mGrid%Nz
		np=(l-1)*(m-1)*(n-1)
    
		allocate(r(np),x(np))

		r=R_ZERO
		x=R_ZERO
    
		flag=.true.
		call copy_Div(inPhi,flag,x)

    do i=1,np
       r(i)=r(i)+sppot(i)*x(i)
       do k=ijapot(i),ijapot(i+1)-1
         r(i)=r(i)+sppot(k)*x(ijapot(k))
         r(ijapot(k))=r(ijapot(k))+sppot(k)*x(i)
			end do
		end do

		flag=.false.
		call copy_Div(outPhi,flag,r)

		deallocate(r,x)
  end subroutine DivA


  !********************
  subroutine copy_Div(Ev,flag,Es)
    implicit none
	  type (rscalar) :: Ev
	  logical flag
	  real (kind=prec) :: Es(*)
    !
	  integer i,j,k,l,m,n,ir

		l = mGrid%Nx
		m = mGrid%Ny
		n = mGrid%Nz

		if(flag) then
			ir=1
      do k=2,n
        do j=2,m
          do i=2,l
            Es(ir)=Ev%v(i,j,k)
						ir=ir+1
          end do
				end do
			end do
    else
			ir=1
      do k=2,n
        do j=2,m
          do i=2,l
            Ev%v(i,j,k)=Es(ir)
						ir=ir+1
          end do
				end do
			end do
		end if

  end subroutine copy_Div


  !*******************************
  
  subroutine DivMInvB(inE,outE)
    implicit none
	  type (rscalar), intent(in) :: inE
    type (rscalar), intent(inout) :: outE
	  ! local variables
    integer                     :: l,m,n,np
	  logical                     :: flag
	  real (kind=prec), pointer, dimension(:)    :: x,b

    l = mGrid%Nx
  	m = mGrid%Ny
  	n = mGrid%Nz

	  np=(l-1)*(m-1)*(n-1)
	
	  allocate(x(np),b(np))

    x=R_ZERO
	  b=R_ZERO

    if (.not.inE%allocated) then
      write(0,*) 'inE in MultA_N not allocated yet'
      stop
    end if

    if (.not.outE%allocated) then
      write(0,*) 'outE in MultA_N not allocated yet'
      stop
    end if

    ! Check whether the bounds are the same
    if ((inE%nx == outE%nx).and.&
         (inE%ny == outE%ny).and.&
         (inE%nz == outE%nz)) then

       if ((inE%gridType == outE%gridType)) then

          flag=.true.
	      call copy_Div(inE,flag,b)

 	      call asolver(np,b,x)

	      flag=.false.
	      call copy_Div(outE,flag,x)

       else
          write (0, *) 'MultA_N: not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-complex vectors for MultA_N are not of same size'
    end if

	  deallocate(x,b)

  end subroutine DivMInvB



! ======================================================
    subroutine asolver(n,b,x)
      implicit none
! ======================================================
!
!    Solve Ax=b, A preconditioner, determined by incomplete
!    Cholesky decomposition, stored in sppot1 and ijapot arrays.
!
      integer i,k,m,n
      real(kind=prec) :: b(*),x(*),sum
 
      do i=1, n
        x(i)=R_ZERO
      end do
 
!------forward substitution, solve L(t)*y = b, storing y in x
 
      x(1)=b(1)/sppot1(1)
      do i=1, n-1
        do m=ijapot(i), ijapot(i+1)-1
          x(ijapot(m))=x(ijapot(m))+sppot1(m)*x(i)
        end do
        x(i+1)=(b(i+1)-x(i+1))/sppot1(i+1)
      end do
 
! -----backward substitution, solve L*x = y -----------------
 
      x(n)=x(n)/sppot1(n)
      do i=n-1, 1, -1
        sum=R_ZERO
        do k=ijapot(i), ijapot(i+1)-1
          sum=sum+sppot1(k)*x(ijapot(k))
        end do
        x(i)=(x(i)-sum)/sppot1(i)
      end do

      return
    end subroutine asolver
    

  !**********************************************************************
  ! Purpose is to compute div inE (input magnetic field)
  ! where sigma is the edge conductivity. Thus, in practice, this computes
  ! divergence of currents.
  ! NOTE that conductivity in air is modified to SIGMA_AIR for this
  ! subroutine.
  ! This is done as a separate specialized routine to avoid carrying
  ! around multiple edge conductivities
  subroutine DivC(inE, outScr, outSci)

    implicit none
    type (cvector), intent(in)		         :: inE
    type (rscalar), intent(inout)		 :: outScr, outSci
	  type (cscalar) :: outSc
    integer                                      :: ix, iy, iz

    IF(.not.inE%allocated) THEN
 	  WRITE(0,*) 'inE not allocated in DivC'
 	  STOP
    ENDIF

    IF(.not.outScr%allocated) THEN
 	  WRITE(0,*) 'outScr not allocated in DivC'
 	  STOP
    ENDIF

    IF(.not.outSci%allocated) THEN
 	  WRITE(0,*) 'outScr not allocated in DivC'
 	  STOP
    ENDIF

    if (outScr%gridType == CORNER) then

       ! Check whether all the inputs/ outputs involved are even of the same
       ! size
       if ((inE%nx == outScr%nx).and.&
            (inE%ny == outScr%ny).and.&
            (inE%nz == outScr%nz)) then

          call create_cscalar(inE%grid, outSc, CORNER)
          ! computation done only for internal nodes
          do ix = 2, outSc%nx
             do iy = 2, outSc%ny
                do iz = 2,outSc%grid%nz
                   outSc%v(ix, iy, iz) = &
                        (inE%x(ix,iy,iz)-inE%x(ix - 1,iy,iz)) * &
                        inE%grid%delXinv(ix)    &
                        + (inE%y(ix,iy,iz)-inE%y(ix,iy - 1,iz)) * &
                        inE%grid%delYinv(iy)    &
                        + (inE%z(ix,iy,iz)-inE%z(ix,iy,iz - 1)) * &
                        inE%grid%delZinv(iz)
				           outScr%v(ix, iy, iz)=dreal(outSc%v(ix, iy, iz))
				           outSci%v(ix, iy, iz)=dimag(outSc%v(ix, iy, iz))
                enddo   ! iz
             enddo      ! iy
          enddo         ! ix
          
          call deall_cscalar(outSc) ! release memory, outSc is a local variable of the derived type rscalar

       else
          write(0, *) 'Error: DivC: scalars not same size'
       end if

    else
       write(0, *) 'Error: DivC: output scalar not compatible use'
    end if

  end subroutine DivC ! DivC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getResTensor(RhoE)
    implicit none
		type(rscalar), intent(inOut)   :: RhoE(6)
		! local variables
		integer ia

     do ia=1,6

       if(RhoE(ia)%allocated) then
          if((RhoE(ia)%Ny .ne. rho_C(ia)%Ny).or.  &
					(RhoE(ia)%Nx .ne. rho_C(ia)%Nx) .or.   &
					(RhoE(ia)%Nz .ne. rho_C(ia)%Nz)) then
             call deall_rscalar(RhoE(ia))
             call create_rscalar(rho_C(ia)%grid,RhoE(ia),CENTER)
          endif
       else
          call create_rscalar(rho_C(ia)%grid,RhoE(ia),CENTER)
       endif

       call copy_rscalar(RhoE(ia), rho_C(ia))

		end do

  end subroutine getResTensor
         
end  module modelOperator3D
