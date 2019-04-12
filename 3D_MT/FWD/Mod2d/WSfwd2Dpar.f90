! *****************************************************************************
!  This module contains integer parameters formerly defined
!  in an include file; this module is used by routines which
!  need to know sizes of passed (or locally allocated) arrays
!  Allows dynamic allocation without significant changes to
!  existing codes

!  I have also added definitions of numerical constants
!    used in some routines

module wsfwd2dpar
   use math_constants
   implicit none
   save
   !  integer parameters are calculated and set by initialization
   !   routine
   integer		:: NZ3MX, MMIMX, NZ0MX, NY0MX, MMBMX
   integer		:: NZ1MX, NZ2MX

   real(kind=prec),parameter       :: Mue =  MU_0

   !  air conductivity is fixed ... but this can be reset by any
   !  routine that uses this module.
   real(kind=prec)                :: CondAir = SIGMA_AIR

   contains
     subroutine SetWSparams(Ny,Nz)

     integer, intent(in) :: Ny,Nz

     NZ0MX = Nz+1
     NY0MX = Ny+1
     MMIMX = (NY0MX-1)*(NZ0MX-1)
     NZ3MX = 3*NZ0MX+4
     NZ2MX = NZ0MX+2
     NZ1MX = NZ0MX+1
     MMBMX = 2*NY0MX+2*NZ0MX
     end subroutine


end module wsfwd2dpar
