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
   use wsfwd2dmod
!   use sunperf 
   implicit none

   !  routines that are public
   public		:: Fwd2DsetupTM, Fwd2DsolveTM, Fwd2DdeallTM, &
			UpdateCondTM,UpdateFreqTM,SetBoundTM
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
      Subroutine FWD2DsetupTM(TMgrid,Cond2D,IER)
      !  routine allocates and initializes everything
      !  Inputs
      type (grid2d_t), intent(in)   :: TMgrid
      real (kind=prec), intent(in) :: Cond2D(Ny,Nz)
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
         Ny = TMgrid%Ny
         Nz = TMgrid%Nz
         Nzb =  Nz - TMgrid%Nza

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
            Dy(iy) = TMgrid%Dy(iy)
         enddo
         do iz = 1,Nzb
            Dzb(iz) = TMgrid%Dz(iz+TMgrid%Nza)
         enddo
         ! compute block center differences (actually 2xDistance!)
         Call DistanceBetweenBlocks(Nzb,Dzb,Czb)
	 Call DistanceBetweenBlocks(Ny,Dy,Cy)
      
         ! set Initialization flag
         Initialized = .true.
         IER = 0

         !  using input conductivity, set Crho (solver representation
         !  of conductivity) and solver coefficients ATM, BTM
         call UpdateCondTM(Cond2D)

      else
         ! return error: already initialized
         IER = -1
      endif
      end subroutine FWD2DsetupTM
      
      
!**********************************************************************
      Subroutine Fwd2DsolveTM(HXB,Hsol,IER)
      ! this is the solver driver ... derived from fwdtmmod_LU (WS)
      !  call this after: fwdtmmodsetup and UpdateFreq to set period

      !  This works more or less as in WS code: only solves the system
      !   for boundary conditions specified in HXB, and does not allow for
      !   forcing.  Need to modify to allow specification of source terms.
      !   HXB can be generated in the proper format by a call to setBound2D_TM
      complex(kind=prec), intent(in)	:: HXB(MMBMX)

      !  The solution is returned in array Hsol ...
      !  could make a derived data type to carry soln and info
      !   about grid size
      complex(kind=prec), intent(out)	:: Hsol(Ny+1,Nzb+1)
      integer, intent(out)	::IER
 
      ! local variables
      integer mmi,mmb,kl,ku,is,iz,iy

      mmi = (Ny-1)*(Nzb-1)
      mmb = 2*Ny + 2*Nzb
      kl = Nzb-1
      ku = Nzb-1
      
      Call MulAibWithXb(Nzb,Ny,BTM,HXB,HXI)
      Call ZGBTRS('N',mmi,kl,ku,1,AII,NZ3MX,ipiv,HXI,MMIMX,IER)
      if (IER.NE.0) then
        write(6,*) '!!! ERROR WHILE SOLVING FWD TM  !!!'
        return
      endif

      !   copy interior nodes into solution vector
      is = 1
      do iy = 2,Ny
         do iz = 2,Nzb
            Hsol(iy,iz) = HXI(is)
            is = is + 1
         enddo
      enddo         
      !  copy boundary nodes into solution vector
      ! left
      do iz = 1,Nzb+1
         Hsol(1,iz) = HXB(iz)
      enddo
      ! Earth surface
      do iy = 2,Ny
         Hsol(iy,1) = HXB(iy+Nzb)
      enddo
      ! bottom
      do iy = 2,Ny
         Hsol(iy,Nzb+1) = HXB(iy+Nzb+Ny-1)
      enddo
      ! right
      do iz = 1,Nzb+1
         Hsol(Ny+1,iz) = HXB(Nzb+2*Ny-1+iz)
      enddo


      end subroutine Fwd2DsolveTM


!**********************************************************************
      Subroutine UpdateCondTM(Cond2D)
!      updates resisistivity (internal solver representation)
!        and then computes ATM, BTM (which depend on resistivity)

      real (kind=prec), intent(in) :: Cond2D(Ny,Nz)
      integer	:: iy, iz

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
      end subroutine Fwd2DdeallTM


end module
