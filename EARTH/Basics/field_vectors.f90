! *****************************************************************************
module field_vectors
  ! Module containing the definitions of the field vectors. This module is only
  ! going to be used by the main program, and the memory allocation/deallocation
  ! subroutines

  use griddef
  implicit none


  save

  real(8),allocatable,dimension(:)			  :: divr,divi	! divH correction
  real(8),allocatable,dimension(:)			  :: divsr,divsi  !adj divH correction

  ! Storing the upper and lower boundary values for the field components.
  ! They are initialized as part of Hx,Hy,Hz in init_ringcurrent, and
  ! never changed from their original values. However, if Hx,Hy,Hz are *not*
  ! global, the boundary conditions need to be explicitly saved.
  !type (sparsevecc)							  :: bvH

  ! This variable contains the magnetic fields Hx,Hy,Hz
  ! They are used for the intermediate calculations only.
  complex(8),allocatable,dimension(:)		  :: hx,hy,hz
  complex(8),allocatable,dimension(:)		  :: sx,sy,sz !adj divH correction


  ! This block contains fields defined on edges, and pre-multiplied
  ! by the edge lengths. The fact that vectorh is global is significant,
  ! since this has the effect of the SAVE attribute in subroutines;
  ! each time vectorh is recomputed, its' previous value is passed on
  ! as the initial value. This feature is preserved from the former
  ! versions of the forward solver. It greatly speeds up the convergence.
  !complex(8),allocatable,dimension(:)		  :: vectorh,vectorb


end module field_vectors
