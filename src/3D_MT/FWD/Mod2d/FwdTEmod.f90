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
   use wsfwd2dmod
!   use sunperf 
   implicit none

   !  routines that are public
   public	:: Fwd2DsetupTE, Fwd2DsolveTE, Fwd2DdeallTE, &
			UpdateCondTE,UpdateFreqTE,SetBoundTE
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
   real (kind=prec), allocatable, dimension(:,:),private	:: ATE,CCon
   real (kind=prec), allocatable, dimension(:),private	:: BTE
   complex (kind=prec), allocatable, dimension(:,:),private :: AII
   complex (kind=prec), allocatable, dimension(:),private	:: EXI
   integer, allocatable, dimension(:),private	:: ipiv
   logical		:: Initialized  = .false.

   Contains
   
   ! *****************************************************************************
      Subroutine FWD2DsetupTE(TEgrid,IER)
      !  routine allocates and initializes module grid variables ... 
      !    no longer sets conductivity or initializes any coefficients
      !
      !  Inputs
      type (grid2d_t), intent(in)   :: TEgrid
      !  outputs ... just allocates saved module arrays 
      !         IER = 0 if everyting works, -1 otherwise
      !         could add different error codes for different
      !         errors ... not much is checked here!
      integer, intent(out)	:: IER

      !  local variables
      integer	::  iy,iz
!Bug fixing:
! In case of Nx > Ny the current version creates a deallocation error in UpdateCondTE. 
! This is because the Initializition is done ONLY once using Ny=Nx and dy=dx (from the first mode).
! However, for the 2nd mode Ny=Ny and dy=dy
! Fixing this problem is done by moving the defintion of Ny,Nz,Nza,Dy,Dz ouside the if statement

      if(.not.Initialized) then

         ! allocate arrays for use within module
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
         ! set Initialization flag
         Initialized = .true.
         IER = 0

      else
         ! return error: already initialized
         IER = -1
      endif
      
         ! initialize ... 
         !  first set array sizes using WS names
         Ny = TEgrid%Ny
         Nz = TEgrid%Nz
         Nza = TEgrid%Nza
             
          do iy = 1,Ny
            Dy(iy) = TEgrid%Dy(iy)
         enddo
         do iz = 1,Nz
            Dz(iz) = TEgrid%Dz(iz)
         enddo
         ! compute block center differences (actually 2xDistance!)
         Call DistanceBetweenBlocks(Nz,Dz,Cz)
	     Call DistanceBetweenBlocks(Ny,Dy,Cy)     
      
      
      
      end subroutine FWD2DsetupTE


!**********************************************************************
      Subroutine Fwd2DsolveTE(EXB,Esol,IER)
      ! this is the solver driver ... derived from Fwd2DTE_LU (WS)
      !  call this after: FWD2DTEsetup and UpdateFreq to set period

      !  This works more or less as in WS code: only solves the system
      !   for boundary conditions specified in HXB, and does not allow for
      !   forcing.  Need to modify to allow specification of source terms.
      !   EXB can be generated in the proper format by a call to setBound2D_TE
      complex(kind=prec), intent(in)	:: EXB(MMBMX)

      !  The solution is returned in array Esol ...
      !  could make a derived data type to carry soln and info
      !   about grid size
      complex(kind=prec), intent(out)	:: Esol(Ny+1,Nz+1)
      integer, intent(out)	::IER
 
      ! local variables
      integer mmi,mmb,kl,ku,is,iz,iy

      mmi = (Ny-1)*(Nz-1)
      mmb = 2*Ny + 2*Nz
      kl = Nz-1
      ku = Nz-1
      
      Call MulAibWithXb(Nz,Ny,BTE,EXB,EXI)
      Call ZGBTRS('N',mmi,kl,ku,1,AII,NZ3MX,ipiv,EXI,MMIMX,IER)
      if (IER.NE.0) then
        write(6,*) '!!! ERROR WHILE SOLVING FWD TE  !!!'
        return
      endif

      !   copy interior nodes into solution vector
      is = 1
      do iy = 2,Ny
         do iz = 2,Nz
            Esol(iy,iz) = EXI(is)
            is = is + 1
         enddo
      enddo         
      !  copy boundary nodes into solution vector
      ! left
      do iz = 1,Nz+1
         Esol(1,iz) = EXB(iz)
      enddo
      ! Earth surface
      do iy = 2,Ny
         Esol(iy,1) = EXB(iy+Nz)
      enddo
      ! bottom
      do iy = 2,Ny
         Esol(iy,Nz+1) = EXB(iy+Nz+Ny-1)
      enddo
      ! right
      do iz = 1,Nz+1
         Esol(Ny+1,iz) = EXB(Nz+2*Ny-1+iz)
      enddo


      end subroutine Fwd2DsolveTE


!**********************************************************************
      Subroutine UpdateCondTE(Cond2D)

!      updates resisistivity (internal solver representation)
!        and then computes ATE, BTE (which depend on resistivity)

      real (kind=prec), intent(in) :: Cond2D(Ny,Nz)
      integer	:: iy, iz

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
         Call SetupA_TE(Nza,Nz,Ny,Dz,Dy,Cz,Cy,CCon,ATE,BTE)

      end subroutine UpdateCondTE


!**********************************************************************
      Subroutine UpdateFreqTE(per)

!      updates frequency (period) dependence, modifying AII
      real(kind=prec), intent(in)	:: per 

      Call FormAII(per,Nz,Ny,ATE,AII,ipiv)

      end subroutine UpdateFreqTE


!**********************************************************************
      Subroutine SetBoundTE(per,EXB)

!     wrapper for WSfwdMod routine SetBound2D_TE 
      real(kind=prec),intent(in)		:: per
      complex(kind=prec), intent(inout)	:: EXB(MMBMX)

      Call SetBound2D_TE(per,Nza,Nz,Ny,Dz,Dy,CCon,EXI,EXB)

      end subroutine SetBoundTE
      
      
!**********************************************************************
      Subroutine Fwd2DdeallTE()

      ! Deallocates arrays used within module
      ! and sets Initialized = .False.

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
      end subroutine Fwd2DdeallTE


end module
