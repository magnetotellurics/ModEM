! *****************************************************************************
!  This module is a "wrapped" version of W Siripunvaraporn
!  2D forward modeling code ... quick and dirty way to
!  get a working 2D code.  This module is for TE mode only.
!  This module provides storage for all coefficient and work
!   arrays, along with routines for seting these up and using
!   for solution of equations.  The public routines in this module
!   are intended to be used for all steps required for forward
!   modeling.

!  Philosphy: keep as many routines untouched as possible,
!   mostly rewriting only those routines that interface with the
!   outside world.  The old WS routines and variable names
!   are mosly all private to this module.  (Note: most routines
!   still require some modification, but only with regard to
!   parameter lists and declarations)

module fwdtemod
   use math_constants
   use wsfwd2d
   use SolnSpace
   use ModelSpace

   implicit none

   !  routines that are public
   public	:: 	Fwd2DsetupTE, Fwd2DsolveTE, Fwd2DdeallTE, &
			UpdateCondTE,UpdateFreqTE,SetBoundTE
   private	::	multEarrayByArea, IntVecToEarray,EarrayToIntVec

   !  variables declared in module header are available
   !  to all module routines, and are saved between calls
   save

   !  integer parameters are calculated and set by initialization
   !   routine
   !  integers describing grid size (for internal use)
   integer, private	:: Nza,Nz,Ny
   !  real arrays describing grid spacing (for internal use)
   real (kind=prec), allocatable, dimension(:),private	:: Dz, Dy, Cz, Cy
   !  real and complex arrays used for solving the equations (internal use)
   real (kind=prec), allocatable, dimension(:,:),private:: ATE,CCon
   real (kind=prec), allocatable, dimension(:),private	:: BTE
   complex (kind=prec), allocatable, dimension(:,:),private :: AII
   complex (kind=prec), allocatable, dimension(:),private:: EXI
   integer, allocatable, dimension(:),private	:: ipiv
   logical		:: Initialized  = .false.

   Contains

   ! *****************************************************************************

      Subroutine FWD2DsetupTE(grid,Sigma,IER)
      !  routine allocates and initializes everything
      type (grid_t), intent(in)	    :: grid
      type (modelParam_t), intent(in)	:: Sigma
      !  allocates saved module arrays
      !  IER is only output ...
      !     = 0 if everyting works, -1 otherwise
      !         could add different error codes for different
      !         errors ... not much is checked here!
      integer, intent(out)	:: IER

      !  local variables
      integer	::  iy,iz

      if(.not.Initialized) then
         ! initialize ...
         !  first set array sizes using WS names
         Ny = grid%Ny
         Nz = grid%Nz
         Nza = grid%Nza

         ! allocate arrays for use within module
         !   NOTE: the array size parameters are all set by a call
         !        to SetWSparams(Ny,Nz)
         allocate(Dz(NZ0MX))
         allocate(Dy(NY0MX))
         allocate(Cz(NZ0MX))
         allocate(Cy(NY0MX))
         allocate(CCon(NZ0MX,NY0MX))
         allocate(ATE(MMIMX,4))
         allocate(BTE(MMBMX))
         allocate(AII(NZ3MX,MMIMX))
         allocate(EXI(MMIMX))
         allocate(ipiv(MMIMX))

         do iy = 1,Ny
            Dy(iy) = grid%Dy(iy)
         enddo
         do iz = 1,Nz
            Dz(iz) = grid%Dz(iz)
         enddo
         ! compute block center differences (actually 2xDistance!)
         call DistanceBetweenBlocks(Nz,Dz,Cz)
	 call DistanceBetweenBlocks(Ny,Dy,Cy)

         ! set Initialization flag
         Initialized = .true.
         IER = 0

         !  using input conductivity, set Crho (solver representation
         !  of conductivity) and solver coefficients ATE, BTE
         call UpdateCondTE(Sigma)

      else
         ! return error: already initialized
         IER = -1
      endif
      end subroutine FWD2DsetupTE

!**********************************************************************

      Subroutine Fwd2DsolveTE(b,Esol,IER)
      ! this is the solver driver ... derived from Fwd2DTE_LU (WS)
      !  call this after: FWD2DTEsetup and UpdateFreq to set period

      !  This works more or less as in WS code: only solves the system
      !   for boundary conditions specified in HXB, and does not allow for
      !   forcing.  Need to modify to allow specification of source terms.
      !   EXB can be generated in the proper format by a call to setBound2D_TE

      !   Oct 21, 2005: new version that takes a more general rhs data
      !   structure of type rhs2d, allowing for more general forcing +
      !   solution of adjoint problem ... not yet coded to return proper BC
      !   for adjoint solution!

      type(rhsVector_t), intent(in)	:: b
      !  The solution is returned in array Esol ...
      !  could make a derived data type to carry soln and info
      !   about grid size
      complex(kind=prec), intent(out)	:: Esol(Ny+1,Nz+1)
      integer, intent(out)	::IER
      character*1		:: NTC

      ! local variables
      integer mmi,mmb,kl,ku,iz,iy

      !  define solution option for call to lapack back-substitution routine
      select case(b%adj)
         case('FWD')
            NTC = 'N'
         case('TRN')
            NTC = 'T'
         case('ADJ')
            NTC = 'C'
         case default
            print *,'this solution case for ADJ not coded'
            IER = -10
            return
      end select

      mmi = (Ny-1)*(Nz-1)
      mmb = 2*Ny + 2*Nz
      kl = Nz-1
      ku = Nz-1

      if ((b%nonzero_bc).and.(b%adj.eq.'FWD')) then
         call MulAibWithXb(Nz,Ny,BTE,b%bc,EXI)
      else
         EXI = C_ZERO
      endif

      if (b%nonzero_source) then
         if(b%adj.eq.'FWD') then
            ! pre-multiply source by symmetrization weights
            call multEarrayByArea(b%source%v,Esol)
            !  add b%source to EXI
            call addEarrayToIntVec(Esol,EXI)
         else
            ! just add b%source to EXI
            call addEarrayToIntVec(b%source%v,EXI)
	 endif
      endif

      !  same call for forward, adjoint (or transposed) problems
      call ZGBTRS(NTC,mmi,kl,ku,1,AII,NZ3MX,ipiv,EXI,MMIMX,IER)

      if (IER.NE.0) then
        write(6,*) '!!! ERROR WHILE SOLVING FWD TE  !!!'
        return
      endif

      !   copy interior nodes into solution vector
      Esol = C_ZERO
      call IntVecToEarray(EXI,Esol)

      if(b%adj .eq. 'FWD') then
         if(b%nonzero_bc) then
            !  copy boundary nodes into solution vector
            !  left side of domain
            ! NOTE: not coded for adjt BC!
            do iz = 1,Nz+1
               Esol(1,iz) = b%bc(iz)
            enddo
            ! Top of domain
            do iy = 2,Ny
            Esol(iy,1) = b%bc(iy+Nz)
            enddo
            ! bottom of domain
            do iy = 2,Ny
               Esol(iy,Nz+1) = b%bc(iy+Nz+Ny-1)
            enddo
            ! right side of domain
	    do iz = 1,Nz+1
               Esol(Ny+1,iz) = b%bc(Nz+2*Ny-1+iz)
            enddo
	 else
	    Esol(1,:) = 0.
	    Esol(:,1) = 0.
	    Esol(:,Nz+1) = 0.
	    Esol(Ny+1,:) = 0.
	 endif
      else
         ! for adjoint case, post-multiply by area weights
         call multEarrayByArea(Esol,Esol)
         !  coding of output BC for adjoint case not done yet
      endif

      end subroutine Fwd2DsolveTE

!**********************************************************************

      Subroutine UpdateCondTE(Sigma)

!      updates resisistivity (internal solver representation)
!        and then computes ATE, BTE (which depend on resistivity)

      type (modelParam_t), intent(in)	:: Sigma

      !  local variables
      real(kind=prec)	::	Cond2D(Ny,Nz)
      integer	:: iy, iz

         call ModelParamToCell(Sigma,Ny,Nz,Cond2D)
         ! Copy inputs into local grid variables
         !  NOTE: order of indices changed for compatability
         !    with 2D TE routines
         do iy = 1,Ny
            do iz = 1,Nz
               CCon(iz,iy) = Cond2D(iy,iz)
            enddo
         enddo

         ! initialize ATE, BTE   (NOTE: all input and output
         !         arguments are module variables)
         call SetupA_TE(Nza,Nz,Ny,Dz,Dy,Cz,Cy,CCon,ATE,BTE)

      end subroutine UpdateCondTE

!**********************************************************************

      Subroutine UpdateFreqTE(per)

!      updates frequency (period) dependence, modifying AII
      real(kind=prec), intent(in)	:: per

      call FormAII(per,Nz,Ny,ATE,AII,ipiv)

      end subroutine UpdateFreqTE

!**********************************************************************

      Subroutine SetBoundTE(per,EXB)

!     wrapper for WSfwdMod routine SetBound2D_TE
      real(kind=prec),intent(in)		:: per
      complex(kind=prec), intent(inout)	:: EXB(MMBMX)

      call SetBound2D_TE(per,Nza,Nz,Ny,Dz,Dy,CCon,EXI,EXB)

      end subroutine SetBoundTE

!**********************************************************************

      Subroutine Fwd2DdeallTE()

      ! Deallocates arrays used within module
      ! and sets Initialized = .False.
         if (Initialized) then
         deallocate(Dz)
         deallocate(Dy)
         deallocate(Cz)
         deallocate(Cy)
         deallocate(CCon)
         deallocate(ATE)
         deallocate(BTE)
         deallocate(AII)
         deallocate(EXI)
         deallocate(ipiv)
         Initialized = .False.
		 end if

      end subroutine Fwd2DdeallTE

!**********************************************************************

     subroutine EarrayToIntVec(E,V)
        !  copies complex 2D array E (interior nodes only)
        !  into forcing or solution vector
        complex(kind=prec), intent(in)       :: E(Ny+1,Nz+1)
        complex(kind=prec), intent(out)       :: V(:)
        integer :: iy,iz,is

        is = 1
        do iy = 2,Ny
           do iz = 2,Nz
              V(is) = E(iy,iz)
              is = is + 1
           enddo
        enddo
     end subroutine EarrayToIntVec

!**********************************************************************

     subroutine IntVecToEarray(V,E)
        !  copy complex 2D array E (interior nodes only)
        !  into forcing or solution vector
        complex(kind=prec), intent(in)       :: V(:)
        complex(kind=prec), intent(out)       :: E(Ny+1,Nz+1)
        integer :: iy,iz,is

        is = 1
        do iy = 2,Ny
           do iz = 2,Nz
              E(iy,iz) = V(is)
              is = is + 1
           enddo
        enddo
     end subroutine IntVecToEarray

!**********************************************************************

     subroutine addEarrayToIntVec(E,V)
        !  adds complex 2D array E (interior nodes only)
        !  into forcing or solution vector
        complex(kind=prec), intent(in)       :: E(Ny+1,Nz+1)
        complex(kind=prec), intent(inout)       :: V(:)
        integer :: iy,iz,is

        is = 1
        do iy = 2,Ny
           do iz = 2,Nz
              V(is) = V(is) + E(iy,iz)
              is = is + 1
           enddo
        enddo
     end subroutine addEarrayToIntVec

!**********************************************************************

     subroutine multEarrayByArea(Ein,Eout)
        !  multiplies comlex 2D array Ein (interior nodes only)
        !  by area weights Cy*Cz
        !  Input may overwrite output

        complex(kind=prec), intent(in)       :: Ein(Ny+1,Nz+1)
        complex(kind=prec), intent(inout)       :: Eout(Ny+1,Nz+1)
        integer :: iy,iz

        do iy = 2,Ny
           do iz = 2,Nz
              Eout(iy,iz) = Cy(iy)*Cz(iz)*Ein(iy,iz)
           enddo
        enddo
     end subroutine multEarrayByArea

end module
