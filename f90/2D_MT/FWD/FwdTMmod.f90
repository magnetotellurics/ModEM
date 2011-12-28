! *****************************************************************************
!  This module is a "wrapped" version of W Siripunvaraporn
!  2D forward modeling code ... quick and dirty way to
!  get a working 2D code.  This module is for TM mode only.
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

module fwdtmmod
   use math_constants
   use wsfwd2d
   use SolnSpace
   use ModelSpace

   implicit none

   !  routines that are public
   public		:: Fwd2DsetupTM, Fwd2DsolveTM, Fwd2DdeallTM, &
			UpdateCondTM,UpdateFreqTM,SetBoundTM
   private		:: multHarrayByArea,addHarrayToIntVec, &
			HarrayToIntVec
   !  variables declared in module header are available
   !  to all module routines, and are saved between calls
   save

   !  integer parameters are calculated and set by initialization
   !   routine
   !  integers describing grid size (for internal use)
   integer, private	:: Nzb,Nz,Ny
   !  real arrays describing grid spacing (for internal use)
   real (kind=prec), allocatable, dimension(:),private	:: Dzb, Dy, Czb, Cy
   !  real and complex arrays used for solving the equations (internal use)
   real (kind=prec), allocatable, dimension(:,:),private	:: ATM,CRho
   real (kind=prec), allocatable, dimension(:),private	:: BTM
   complex (kind=prec), allocatable, dimension(:,:),private :: AII
   complex (kind=prec), allocatable, dimension(:),private	:: HXI
   integer, allocatable, dimension(:),private	:: ipiv
   logical		:: Initialized  = .false.

   Contains


   ! *****************************************************************************
      Subroutine FWD2DsetupTM(grid,Sigma,IER)
      !  routine allocates and initializes everything
      !  Inputs
      type (grid_t), intent(in)     :: grid
      type (modelParam_t), intent(in) :: Sigma
      !  outputs ... just allocates saved module arrays
      !         IER = 0 if everyting works, -1 otherwise
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
         Nzb =  Nz - grid%Nza

         ! allocate arrays for use within module
         allocate(Dzb(NZ0MX))
         allocate(Dy(NY0MX))
         allocate(Czb(NZ0MX))
         allocate(Cy(NY0MX))
         allocate(Crho(NZ0MX,NY0MX))
         allocate(ATM(MMIMX,4))
         allocate(BTM(MMBMX))
         allocate(AII(NZ3MX,MMIMX))
         allocate(HXI(MMIMX))
         allocate(ipiv(MMIMX))

         do iy = 1,Ny
            Dy(iy) = grid%Dy(iy)
         enddo
         do iz = 1,Nzb
            Dzb(iz) = grid%Dz(iz+grid%Nza)
         enddo
         ! compute block center differences (actually 2xDistance!)
         Call DistanceBetweenBlocks(Nzb,Dzb,Czb)
	 Call DistanceBetweenBlocks(Ny,Dy,Cy)

         ! set Initialization flag
         Initialized = .true.
         IER = 0

         !  using input conductivity, set Crho (solver representation
         !  of conductivity) and solver coefficients ATM, BTM
         call UpdateCondTM(Sigma)

      else
         ! return error: already initialized
         IER = -1
      endif
      end subroutine FWD2DsetupTM

!**********************************************************************
      Subroutine Fwd2DsolveTM(b,Hsol,IER)
	!!!!!!!   DOCUMENTATION KNOWN TO BE INCORRECT ... MODIFIED
	!!!!!!!   EXTENSIVELY BY WS WITHOUT MODIFICATION OF DOCUMENTATION!!!!!
      !Subroutine Fwd2DsolveTM(HXB,Hsol,IER)
      ! this is the solver driver ... derived from fwdtmmod_LU (WS)
      !  call this after: fwdtmmodsetup and UpdateFreq to set period

      !  This works more or less as in WS code: only solves the system
      !   for boundary conditions specified in HXB, and does not allow for
      !   forcing.  Need to modify to allow specification of source terms.
      !   HXB can be generated in the proper format by a call to setBound2D_TM
      type (rhsVector_t), intent(in)   :: b

      !  The solution is returned in array Hsol ...
      !  could make a derived data type to carry soln and info
      !   about grid size
      complex(kind=prec), intent(out)	:: Hsol(Ny+1,Nzb+1)
      integer, intent(out)	::IER
      character*1    :: NTC

      ! local variables
      integer mmi,mmb,kl,ku,iz,iy

      ! define solution option for call to lapack back-substitution routine
      select case (b%adj)
         case('FWD')
           NTC = 'N'
         case('TRN')
           NTC = 'T'
         case('ADJ')
           NTC = 'C'
         case default
           print*,'this solution case for ADJ not coded'
           IER = -10
           return
      end select

      mmi = (Ny-1)*(Nzb-1)
      mmb = 2*Ny + 2*Nzb
      kl = Nzb-1
      ku = Nzb-1

      if ((b%nonzero_bc).and.(b%adj.eq.'FWD')) then
        Call MulAibWithXb(Nzb,Ny,BTM,b%bc,HXI)
      else
	HXI = C_ZERO
      endif

      if (b%nonzero_source) then
         if(b%adj.eq.'FWD') then
            ! pre-multiply source by symmetrization weights
            call multHarrayByArea(b%source%v,Hsol)

            !  add b%source to HXI
            call addHarrayToIntVec(Hsol,HXI)
         else
            ! just add b%source to HXI
            call addHarrayToIntVec(b%source%v,HXI)
         endif
      endif

      ! same call for forward, adjoint (or transposed) problems
      Call ZGBTRS(NTC,mmi,kl,ku,1,AII,NZ3MX,ipiv,HXI,MMIMX,IER)

      if (IER.NE.0) then
        write(6,*) '!!! ERROR WHILE SOLVING FWD TM  !!!'
        return
      endif

      !   copy interior nodes into solution vector
      Hsol = C_ZERO
      call IntVecToHarray(HXI,Hsol)

      if (b%adj.eq. 'FWD') then
	if (b%nonzero_bc) then
          !  copy boundary nodes into solution vector
          ! left
          do iz = 1,Nzb+1
             Hsol(1,iz) = b%bc(iz)
          enddo
          ! Earth surface
          do iy = 2,Ny
             Hsol(iy,1) = b%bc(iy+Nzb)
          enddo
          ! bottom
          do iy = 2,Ny
             Hsol(iy,Nzb+1) = b%bc(iy+Nzb+Ny-1)
          enddo
          ! right
          do iz = 1,Nzb+1
             Hsol(Ny+1,iz) = b%bc(Nzb+2*Ny-1+iz)
          enddo
        else
	  Hsol(1,:) = 0.
	  Hsol(:,1) = 0.
	  Hsol(:,Nzb+1) = 0.
	  Hsol(Ny+1,:)  = 0.
        endif
      else
        ! for adjoint case, post-multiply by area weights
        call multHarrayByArea(Hsol,Hsol)
	! coding of output BC for adjoint case not done yet
      endif

      end subroutine Fwd2DsolveTM


!**********************************************************************
      Subroutine UpdateCondTM(Sigma)
!      updates resisistivity (internal solver representation)
!        and then computes ATM, BTM (which depend on resistivity)

      type (modelParam_t), intent(in)	:: Sigma

      !  local variables
      integer				:: iy, iz

      real (kind=prec)  :: Cond2D(Ny,Nz)

         call ModelParamToCell(Sigma,Ny,Nz,Cond2D)
         ! Copy inputs into local grid variables
         !  NOTE: order of indices changed for compatability
         !    with 2D TM routines
         do iy = 1,Ny
            do iz = 1,Nzb
               Crho(iz,iy) = ONE/Cond2D(iy,iz+Nz-Nzb)
            enddo
         enddo

         ! initialize ATM, BTM   (NOTE: all input and output
         !         arguments are module variables)
         Call SetupA_TM(Nzb,Ny,Dzb,Dy,Czb,Cy,CRho,ATM,BTM)


      end subroutine UpdateCondTM


!**********************************************************************
      Subroutine UpdateFreqTM(per)
!      updates frequency (period) dependence, modifying AII
      real(kind=prec), intent(in)	:: per

      Call FormAII(per,Nzb,Ny,ATM,AII,ipiv)

      end subroutine UpdateFreqTM


!**********************************************************************
      Subroutine SetBoundTM(per,HXB)
!     wrapper for WSfwdMod routine SetBound2D_TM
      real(kind=prec),intent(in)		:: per
      complex(kind=prec), intent(inout)	:: HXB(MMBMX)

      Call SetBound2D_TM(per,Nzb,Ny,Dzb,Dy,CRho,HXI,HXB)

      end subroutine SetBoundTM


!**********************************************************************
      Subroutine Fwd2DdeallTM()

      ! Deallocates arrays used within module
      ! and sets Initialized = .False.
        if (Initialized) then
         deallocate(Dzb)
         deallocate(Dy)
         deallocate(Czb)
         deallocate(Cy)
         deallocate(Crho)
         deallocate(ATM)
         deallocate(BTM)
         deallocate(AII)
         deallocate(HXI)
         deallocate(ipiv)
         Initialized = .False.
		end if
 
      end subroutine Fwd2DdeallTM

!**********************************************************************

      subroutine HarrayToIntVec(H,V)
        !  copies complex 2D array E (interior nodes only)
        !  into forcing or solution vector
        complex(kind=prec), intent(in)       :: H(Ny+1,Nzb+1)
        complex(kind=prec), intent(out)      :: V(:)
        integer :: iy,iz,is

        is = 1
        do iy = 2,Ny
           do iz = 2,Nzb
              V(is) = H(iy,iz)
              is = is + 1
           enddo
        enddo

      end subroutine HarrayToIntVec

!**********************************************************************

      subroutine IntVecToHarray(V,H)
        !  copy complex 2D array H (interior nodes only)
        !  into forcing or solution vector
        complex(kind=prec), intent(out)       :: H(Ny+1,Nzb+1)
        complex(kind=prec), intent(in)      :: V(:)
        integer :: iy,iz,is

        is = 1
        do iy = 2,Ny
           do iz = 2,Nzb
              H(iy,iz) = V(is)
              is = is + 1
           enddo
        enddo

      end subroutine IntVecToHarray

!**********************************************************************

      subroutine addHarrayToIntVec(H,V)
         !  adds complex 2D array E (interior nodes only)
         !  into forcing or solution vector
         complex(kind=prec), intent(in)       :: H(Ny+1,Nzb+1)
         complex(kind=prec), intent(inout)    :: V(:)
         integer :: iy,iz,is

         is = 1
         do iy = 2,Ny
            do iz = 2,Nzb
               V(is) = V(is) + H(iy,iz)
               is = is + 1
            enddo
         enddo
      end subroutine addHarrayToIntVec

!**********************************************************************

     subroutine multHarrayByArea(Hin,Hout)
     !  multiplies comlex 2D array Ein (interior nodes only)
     !  by area weights Cy*Cz
     !  Input may overwrite output

     complex(kind=prec), intent(in)       :: Hin(Ny+1,Nz+1)
     complex(kind=prec), intent(inout)    :: Hout(Ny+1,Nz+1)
     integer :: iy,iz

     do iy = 2,Ny
       do iz = 2,Nzb
         Hout(iy,iz) = Cy(iy)*Czb(iz)*Hin(iy,iz)
       enddo
     enddo

     end subroutine multHarrayByArea

end module
