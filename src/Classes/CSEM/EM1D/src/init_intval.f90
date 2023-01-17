!---------------------------------------------------------
!> EM1D subroutine init_intval
!
!> allocate vectors for Hankel integral values
!>   the number of integrals required depends on the computation job (fwd model / sens.)
!>   and the anisotropy
!
!> Rita Streich 2011
!---------------------------------------------------------
subroutine init_intval(refl_var,bgdat)

  implicit none

  !external variables
  type(refl_struct) :: refl_var   !all variables that have to be remembered while computing 1D fields
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields

  !internal variables
  integer(kind=int32) :: ierr        !error index


  if(bgdat%dowhat .EQ. fwdmodel) then !fields only
    allocate(refl_var%intvalre(refl_var%nrad,nintHED),refl_var%intvalim(refl_var%nrad,nintHED), &
             refl_var%spl_derivre(refl_var%nrad,nintHED),refl_var%spl_derivim(refl_var%nrad,nintHED), stat=ierr)
  else !sensitivities (and possibly fields)
    if(bgdat%aniso .EQ. iso) then !isotropic
      allocate(refl_var%intvalre(refl_var%nrad,nintHEDd),refl_var%intvalim(refl_var%nrad,nintHEDd), &
               refl_var%spl_derivre(refl_var%nrad,nintHEDd),refl_var%spl_derivim(refl_var%nrad,nintHEDd), stat=ierr)
    else !VTI-anisotropic
      allocate(refl_var%intvalre(refl_var%nrad,nintHEDdvti),refl_var%intvalim(refl_var%nrad,nintHEDdvti), &
               refl_var%spl_derivre(refl_var%nrad,nintHEDdvti),refl_var%spl_derivim(refl_var%nrad,nintHEDdvti), stat=ierr)
    endif
  endif
  if(ierr.NE.0) call alloc_error(pid,'init_intval','refl_var%intvalre etc.',ierr)


endsubroutine init_intval
