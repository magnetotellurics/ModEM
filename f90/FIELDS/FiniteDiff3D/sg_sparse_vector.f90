! *****************************************************************************
module sg_sparse_vector
  ! Sparse complex vector operations. At present, sparse vectors are used for
  ! representation of data functionals. Only complex vectors are supported, and
  ! only those routines needed for modeling and initial work on inversion have
  ! been developed.
  ! Belongs to SG_Basics class: staggered cartesian grid, data
  ! types defined on this grid, and operations defined on these data types. Not
  ! specific to EM problem, no dependency on outside (from other classes)
  !  modules.

  use math_constants
  use griddef
  use sg_vector
  implicit none

  INTERFACE create
     module procedure create_sparsevecc
  END INTERFACE

  INTERFACE deall
     module procedure deall_sparsevecc
  END INTERFACE

  INTERFACE reall
     module procedure reall_sparsevecc
  END INTERFACE

  INTERFACE scMult
     module procedure scMult_sparsevecc
  END INTERFACE

  INTERFACE linComb
     module procedure linComb_sparsevecc
  END INTERFACE

  INTERFACE add
     module procedure add_scvector
  END INTERFACE

  INTERFACE dotProd
     module procedure dotProd_scvector_f
     !module procedure dotProd_csvector_f
  END INTERFACE

  INTERFACE dotProd_noConj
     module procedure dotProd_noConj_scvector_f
     !module procedure dotProd_noConj_csvector_f
  END INTERFACE

  INTERFACE conjg
     module procedure conjg_sparsevecc_f
  END INTERFACE

  INTERFACE newValue
     module procedure newValueC_sparsevecc
     module procedure newValueR_sparsevecc
     module procedure copyValue_csvector
  END INTERFACE


  interface assignment (=)
     module PROCEDURE copy_sparsevecc
     module PROCEDURE copy_csvector
  end interface

  public			:: sparsevecc
  public			:: create_sparsevecc, deall_sparsevecc
  public            :: read_sparsevecc, write_sparsevecc
  public            :: scMult_sparsevecc
  public			:: newValueC_sparsevecc, newValueR_sparsevecc
  public			:: copyValue_csvector, conjg_sparsevecc_f
  public			:: copy_sparsevecc, linComb_sparsevecc
  public			:: dotProd_scvector_f, dotProd_csvector_f
  public			:: add_scvector

!**************************************************************************
  type :: sparsevecc

     ! complex vector defined on edge/ face nodes;
     ! store the intention of the use in a character string defined
     ! as in GridDef as a parameter: EDGE or FACE
     character (len=80)	                             	:: gridType=''
     ! nCoeff is number of non-zero nodes
     integer 						:: nCoeff  = 0
     ! xyz = 1,2,3 refers to x, y or z components,
     ! i,j,k are arrays of indices that defines grid location
     integer , pointer, dimension(:) 		:: i,j,k,xyz
     ! c is complex array of coefficients
     complex (kind=prec), pointer, dimension(:) 	:: c
     ! has sparse vector been allocated?
     logical					:: allocated = .false.
     ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
     ! (probably will not be needed in the future when compilers will support
     ! ISO/IEC 15581 - the "allocatable array extension")
     logical					:: temporary = .false.
     ! pointer to the parent grid not needed or set: we can
     ! make full use of the sparse vector without it...
     ! - should be passed explicitly if it is ever required!
     !type (grid_t), pointer                 	:: grid

  end type sparsevecc


Contains

  ! **********************************************************************
  ! delete/ deallocate the sparse vector
  subroutine deall_sparsevecc(oldLC)

    implicit none
    type (sparsevecc)                :: oldLC
    integer                          :: status

    if(oldLC%allocated) then
       deallocate(oldLC%i,STAT=status)
       deallocate(oldLC%j, STAT=status)
       deallocate(oldLC%k, STAT=status)
       deallocate(oldLC%xyz, STAT=status)
       deallocate(oldLC%c, STAT=status)
       oldLC%gridType = ''
       oldLC%allocated = .false.
    end if

       !=================================================================
       !==================================================== Added by Oat
       !=================================================================
       oldLC%nCoeff = 0

  end subroutine deall_sparsevecc


  ! **********************************************************************
  ! create an object of type sparsevecc of length nCoeff
  ! pointer to the grid not needed (and in some cases
  ! when sparse vectors are set up, we don't even have it).
  subroutine create_sparsevecc(nCoeff,newLC, gridType)

    implicit none
    integer, intent(in) 					:: nCoeff
    type (sparsevecc), intent(inout) 		:: newLC
    character (len=80), intent(in)     		:: gridType
    integer					:: status

    ! the old baggage is out of the door
    if(newLC%allocated) then
       deallocate(newLC%i, STAT = status)
       deallocate(newLC%j, STAT = status)
       deallocate(newLC%k, STAT = status)
       deallocate(newLC%xyz, STAT = status)
       deallocate(newLC%c, STAT = status)
       newLC%gridType = ''
       newLC%allocated = .false.
    endif

    newLC%allocated = .true.
    allocate(newLC%i(nCoeff),STAT=status)
    newLC%allocated = newLC%allocated .and. (status .eq. 0)
    allocate(newLC%j(nCoeff),STAT=status)
    newLC%allocated = newLC%allocated .and. (status .eq. 0)
    allocate(newLC%k(nCoeff),STAT=status)
    newLC%allocated = newLC%allocated .and. (status .eq. 0)
    allocate(newLC%xyz(nCoeff),STAT=status)
    newLC%allocated = newLC%allocated .and. (status .eq. 0)
    allocate(newLC%c(nCoeff),STAT=status)
    newLC%allocated = newLC%allocated .and. (status .eq. 0)

    newLC%nCoeff = nCoeff
    newLC%i = 0
    newLC%j = 0
    newLC%k = 0
    newLC%xyz = 0
    newLC%c = C_ZERO
    newLC%gridType = gridType

  end subroutine create_sparsevecc

  ! **********************************************************************
  ! Reallocates an object of type sparsevecc. The object has to already be
  ! allocated. If allocated and shorter than nCoeff, more memory is
  ! allocated at the end and the contents are preserved.
  ! If allocated and longer than nCoeff, truncates to the first nCoeff values.
  ! This is useful when we need to store the information somewhere, but do
  ! not yet know the final length of the vector. Once it is fully read and
  ! the number of coefficients is known, use this routine to truncate
  ! to the correct length preserving all the values already stored.
  subroutine reall_sparsevecc(nCoeff,newLC)

    implicit none
    integer, intent(in) 						:: nCoeff
    type (sparsevecc), intent(inout) 			:: newLC
    type (sparsevecc)                           :: tempLC
    integer					                    :: n, status

    ! the old baggage is out of the door
    if(.not. newLC%allocated) then
    	write(0, *) 'The input sparsevecc has to be allocated in reall_sparsevecc'
    	return
    end if

	tempLC = newLC

	if (tempLC%nCoeff .eq. nCoeff) then
		! do nothing
	else
		call deall_sparsevecc(newLC)
    	newLC%allocated = .true.
    	allocate(newLC%i(nCoeff),STAT=status)
    	newLC%allocated = newLC%allocated .and. (status .eq. 0)
    	allocate(newLC%j(nCoeff),STAT=status)
    	newLC%allocated = newLC%allocated .and. (status .eq. 0)
    	allocate(newLC%k(nCoeff),STAT=status)
    	newLC%allocated = newLC%allocated .and. (status .eq. 0)
    	allocate(newLC%xyz(nCoeff),STAT=status)
    	newLC%allocated = newLC%allocated .and. (status .eq. 0)
    	allocate(newLC%c(nCoeff),STAT=status)
    	newLC%allocated = newLC%allocated .and. (status .eq. 0)
    	newLC%gridType = tempLC%gridType
    	newLC%nCoeff = nCoeff
	end if

	if (tempLC%nCoeff > nCoeff) then
		! new vector will be shorter
		do n = 1,nCoeff
    		newLC%i(n) = tempLC%i(n)
    		newLC%j(n) = tempLC%j(n)
    		newLC%k(n) = tempLC%k(n)
    		newLC%xyz(n) = tempLC%xyz(n)
    		newLC%c(n) = tempLC%c(n)
    	end do
	else if (tempLC%nCoeff < nCoeff) then
		! new vector will be longer; copy the old values
		do n = 1,tempLC%nCoeff
    		newLC%i(n) = tempLC%i(n)
    		newLC%j(n) = tempLC%j(n)
    		newLC%k(n) = tempLC%k(n)
    		newLC%xyz(n) = tempLC%xyz(n)
    		newLC%c(n) = tempLC%c(n)
    	end do
    	! ... then pad with zeroes
    	do n = tempLC%nCoeff+1,nCoeff
    		newLC%i(n) = 0
    		newLC%j(n) = 0
    		newLC%k(n) = 0
    		newLC%xyz(n) = 0
    		newLC%c(n) = C_ZERO
    	end do
	end if

    call deall_sparsevecc(tempLC)

  end subroutine reall_sparsevecc

  !******************************************************************************
  ! write_sparsevecc writes a sparse vector  in a simple ASCII format; vector has
  ! to exist and be allocated before calling this routine, and the file unit
  ! must already be available for writing.
  subroutine write_sparsevecc(fid, SV)

      integer,        intent(in)        :: fid
      type (sparsevecc), intent(in)     :: SV

      !  local variables
      integer                           :: ii, istat

      if(.not. SV%allocated) then
         write(0, *) 'sparse vector must be allocated before call to write_sparsevecc'
         return
      endif

      write(fid,'(i12,a10)',iostat=istat) SV%nCoeff,trim(SV%gridType)

      do ii = 1,SV%nCoeff
         write(fid,'(4i5,2es14.6)',iostat=istat) SV%i(ii),SV%j(ii),SV%k(ii),SV%xyz(ii),real(SV%c(ii)),aimag(SV%c(ii))
      end do

  end subroutine write_sparsevecc

  !**********************************************************************************
  ! read_sparsevecc reads a sparse vector in a simple ASCII format; vector must match
  ! the input grid; file unit must already be available for reading.
  subroutine read_sparsevecc(fid, SV, grid)

      integer,        intent(in)        :: fid
      type (sparsevecc), intent(inout)  :: SV
      type (grid_t), intent(in), optional :: grid

      !  local variables
      integer                           :: nCoeff
      character(80)                     :: gridType
      integer                           :: ii, istat
      real (kind(SV%c))                  :: cr, ci

      read(fid,*,iostat=istat) nCoeff,gridType

      if(.not. SV%allocated) then
        call create_sparsevecc(nCoeff,SV,gridType)
      else
        call reall_sparsevecc(nCoeff,SV)
        SV%gridType = gridType
      endif

      do ii = 1,nCoeff
        ! makes the acceptable formatting more flexible
        read(fid,*,iostat=istat) SV%i(ii),SV%j(ii),SV%k(ii),SV%xyz(ii),cr,ci
        SV%c(ii) = cmplx(cr,ci)
      end do

      ! if grid is available, make sure the two are consistent
      if (present(grid)) then
        if ((maxval(SV%i) > grid%nx) .or. (maxval(SV%j) > grid%ny) .or. (maxval(SV%k) > grid%nz)) then
            write(0, *) 'sparse vector size does not match the grid in read_sparsevecc'
            return
        endif
      endif

  end subroutine read_sparsevecc

  ! **********************************************************************
  ! this will copy a sparse complex vector from SV1 to SV2 ...
  ! interface to =
  ! basically like copy commands for vectors, scalars, BC
  ! check for size consistency (nCoeff), reallocate output if needed
  ! note that before allocation nCoeff = 0
  ! Remember, SV2 = SV1
  subroutine copy_sparsevecc(SV2,SV1)

    implicit none
    type (sparsevecc), intent(in)		:: SV1
    type (sparsevecc), intent(inout)		:: SV2
    integer	                     		:: status

    ! check to see if RHS (SV1) is active (allocated)
    if(.not. SV1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_sparsevecc'
       return
    end if

    ! allocate output if needed, otherwise check for consistency
    if(.not. SV2%allocated) then
       call create_sparsevecc(SV1%nCoeff, SV2, SV1%gridType)
    elseif(SV1%nCoeff .ne. SV2%nCoeff) then
        ! internal memory allocation is strongly discouraged. But this
        ! is an exception
        call deall_sparsevecc(SV2)
        ! ... now allocate for correct number of components
        call create_sparsevecc(SV1%nCoeff, SV2, SV1%gridType)
    end if

    ! happen to have the same specs
    if (SV1%gridType == SV2%gridType) then

             ! just copy the components
             SV2%i = SV1%i
             SV2%j = SV1%j
             SV2%k = SV1%k
             SV2%xyz = SV1%xyz
             SV2%c = SV1%c

	else
             write (0, *) 'not compatible usage for copy_sparsevecc'

    end if

    if(SV1%temporary) then
    	call deall_sparsevecc(SV1)
    end if

  end subroutine copy_sparsevecc

  ! **********************************************************************
  ! linear combination of two sparse vectors, output as a sparse vector
  ! allocates (or reallocates) output sparse vector Loc3
  subroutine linComb_sparsevecc(Lic1,ic1,Lic2,ic2,Loc3)

    implicit none
    type (sparsevecc), intent(in)		:: Lic1,Lic2
    type (sparsevecc), intent(inout)		:: Loc3
    complex (kind=prec), intent(in)		:: ic1,ic2

    integer					:: n,m,nm
    integer					:: nCoeffSum
    integer, allocatable, dimension(:)  	:: Lic1oldIndex
    integer	                     		:: status

    allocate(Lic1oldIndex(Lic2%nCoeff), STAT = status)

    ! it all depends on how many nodes are common
    if(Loc3%allocated) then
       deallocate(Loc3%i, STAT = status)
       deallocate(Loc3%j, STAT = status)
       deallocate(Loc3%k, STAT = status)
       deallocate(Loc3%xyz, STAT = status)
       deallocate(Loc3%c, STAT = status)
       Loc3%gridType = ''
       Loc3%allocated = .false.
    endif

    if (Lic1%gridType == Lic2%gridType) then

       ! count common indices
       nCoeffSum = Lic1%nCoeff+Lic2%nCoeff
       do m = 1,Lic2%nCoeff
          Lic1oldIndex(m) = 0
          do n = 1,Lic1%nCoeff

             if((Lic1%xyz(n).eq.Lic2%xyz(m)).and.(Lic1%i(n).eq.Lic2%i(m)).and.  &
                  (Lic1%j(n).eq.Lic2%j(m)).and.(Lic1%k(n).eq.Lic2%k(m))) then
                nCoeffSum = nCoeffSum-1
                Lic1oldIndex(m) = n
                exit
             endif

          enddo
       enddo

    else
       write (0, *) 'not compatible usage for LinCompSparseVecC'
    end if

    Call create_sparsevecc(nCoeffsum, Loc3, Lic1%gridType)
    nm = Lic1%nCoeff
    Loc3%i(1:nm) = Lic1%i
    Loc3%j(1:nm) = Lic1%j
    Loc3%k(1:nm) = Lic1%k
    Loc3%xyz(1:nm) = Lic1%xyz
    Loc3%c(1:nm) = ic1*Lic1%c

    do m = 1,Lic2%nCoeff
       ! if none of them are common, just concatenate
       if(Lic1oldIndex(m).eq.0) then
          nm = nm+1
          Loc3%i(nm) = Lic2%i(m)
          Loc3%j(nm) = Lic2%j(m)
          Loc3%k(nm) = Lic2%k(m)
          Loc3%xyz(nm) = Lic2%xyz(m)
          Loc3%c(nm) = ic2*Lic2%c(m)
       else
          Loc3%c(Lic1oldIndex(m)) =  Loc3%c(Lic1oldIndex(m))+ic2*Lic2%c(m)
       endif
    enddo

    deallocate(Lic1oldIndex, STAT = status)

  end subroutine linComb_sparsevecc

  ! **********************************************************************
  ! multiply a sparse complex vector SV1 by a complex value to get SV2.
  ! Sometimes computing a full linear combination by linComb_sparsevecc
  ! is too much work for this simple need... SV2 = ic * SV1
  ! output SV2 can overwrite input SV1
  ! allocates for output if necessary
  subroutine scMult_sparsevecc(cs,SV1,SV2)

    implicit none
    complex(kind=prec), intent(in)	    :: cs
    type (sparsevecc), intent(in)		:: SV1
    type (sparsevecc), intent(inout)    :: SV2

    !  make sure SV2 is allocated and of the correct size
    if(SV2%allocated) then
       if(SV2%gridType .ne. SV1%gridType) then
          write(0,*) 'not compatible usage for scMult_sparsevecc'
          return
       elseif(SV2%nCoeff .ne. SV1%nCoeff) then
          call deall_sparsevecc(SV2)
          call create_sparsevecc(SV1%nCoeff,SV2,SV1%gridType)
       endif
    else
       call create_sparsevecc(SV1%nCoeff,SV2,SV1%gridType)
    endif

    SV2%i = SV1%i
    SV2%j = SV1%j
    SV2%k = SV1%k
    SV2%xyz = SV1%xyz
    SV2%c = cs * SV1%c

  end subroutine scMult_sparsevecc

  ! **********************************************************************
  ! compute complex dot product between a sparse vector SV and a vector of
  ! type cvector ... result in c
  function dotProd_scvector_f(SV,V) result(c)

    implicit none
    type (sparsevecc), intent(in)		:: SV
    type (cvector), intent(in)			:: V
    complex(kind=prec)				:: c
    integer					:: i
    integer					:: xi, yi, zi

    c = C_ZERO

    if((.not.SV%allocated).or.(.not.V%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_scvector_f'
       return
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for dotProd_scvector_f'
       return
    endif

    ! sum over  non-zero terms in sparse vector (conjugate sparse)
    ! (need to check xyz the component)
    ! Remember, xyz = 1,2,3 refers to x, y or z components
    do i = 1,SV%nCoeff

       ! generic test for both edge and face (all the components)
       if ((SV%i(i).le.V%grid%nx+1).or.(SV%j(i).le.V%grid%ny+1).or.&
            (SV%k(i).le.V%grid%nz+1)) then

          ! dealing with x-components
          if (SV%xyz(i) == 1) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + conjg(SV%c(i)) * V%x(xi, yi, zi)

             ! dealing with y-component
          else if (SV%xyz(i) == 2) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + conjg(SV%c(i)) * V%y(xi, yi, zi)

             ! dealing with z-component
          else if (SV%xyz(i) == 3) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + conjg(SV%c(i)) * V%z(xi, yi, zi)
          end if

       else
          write(0,*) 'IJK out of bounds for dotProd_scvector_f'
          return
       endif

    enddo

  end function dotProd_scvector_f

  ! **********************************************************************
  ! compute complex dot product between a sparse vector SV and a vector of
  ! type cvector ... result in c
  function dotProd_csvector_f(V,SV) result(c)

    implicit none
    type (sparsevecc), intent(in)		:: SV
    type (cvector), intent(in)			:: V
    complex(kind=prec)				:: c
    integer					:: i
    integer					:: xi, yi, zi

    c = C_ZERO

    if((.not.SV%allocated).or.(.not.V%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_csvector_f'
       return
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for dotProd_csvector_f'
       return
    endif

    ! sum over  non-zero terms in sparse vector (conjugate full)
    ! (need to check xyz the component)
    ! Remember, xyz = 1,2,3 refers to x, y or z components
    do i = 1,SV%nCoeff

       ! generic test for both edge and face (all the components)
       if ((SV%i(i).le.V%grid%nx+1).or.(SV%j(i).le.V%grid%ny+1).or.&
            (SV%k(i).le.V%grid%nz+1)) then

          ! dealing with x-components
          if (SV%xyz(i) == 1) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * conjg(V%x(xi, yi, zi))

             ! dealing with y-component
          else if (SV%xyz(i) == 2) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * conjg(V%y(xi, yi, zi))

             ! dealing with z-component
          else if (SV%xyz(i) == 3) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * conjg(V%z(xi, yi, zi))
          end if

       else
          write(0,*) 'IJK out of bounds for dotProd_csvector_f'
          return
       endif

    enddo

  end function dotProd_csvector_f

  ! **********************************************************************
  ! compute complex dot product between a sparse vector SV and a vector of
  ! type cvector ... result in c
  !   FOR THIS VERSION FIRST VECTOR IS NOT CONJUGATED
  function dotProd_noConj_scvector_f(SV,V) result(c)

    implicit none
    type (sparsevecc), intent(in)		:: SV
    type (cvector), intent(in)			:: V
    complex(kind=prec)				:: c
    integer					:: i
    integer					:: xi, yi, zi

    c = C_ZERO

    if((.not.SV%allocated).or.(.not.V%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_scvector_f'
       return
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for dotProd_scvector_f'
       return
    endif

    ! sum over  non-zero terms in sparse vector
    ! (need to check xyz the component)
    ! Remember, xyz = 1,2,3 refers to x, y or z components
    do i = 1,SV%nCoeff

       ! generic test for both edge and face (all the components)
       if ((SV%i(i).le.V%grid%nx+1).or.(SV%j(i).le.V%grid%ny+1).or.&
            (SV%k(i).le.V%grid%nz+1)) then

          ! dealing with x-components
          if (SV%xyz(i) == 1) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%x(xi, yi, zi)

             ! dealing with y-component
          else if (SV%xyz(i) == 2) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%y(xi, yi, zi)

             ! dealing with z-component
          else if (SV%xyz(i) == 3) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%z(xi, yi, zi)
          end if

       else
          write(0,*) 'IJK out of bounds for dotProd_scvector_f'
          return
       endif

    enddo

  end function dotProd_noConj_scvector_f


  ! **********************************************************************
  ! compute complex dot product between a sparse vector SV and a vector of
  ! type cvector ... result in c
  !   FOR THIS VERSION FIRST VECTOR IS NOT CONJUGATED
  function dotProd_noConj_csvector_f(V,SV) result(c)

    implicit none
    type (sparsevecc), intent(in)		:: SV
    type (cvector), intent(in)			:: V
    complex(kind=prec)				:: c
    integer					:: i
    integer					:: xi, yi, zi

    c = C_ZERO

    if((.not.SV%allocated).or.(.not.V%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_csvector_f'
       return
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for dotProd_csvector_f'
       return
    endif

    ! sum over  non-zero terms in sparse vector
    ! (need to check xyz the component)
    ! Remember, xyz = 1,2,3 refers to x, y or z components
    do i = 1,SV%nCoeff

       ! generic test for both edge and face (all the components)
       if ((SV%i(i).le.V%grid%nx+1).or.(SV%j(i).le.V%grid%ny+1).or.&
            (SV%k(i).le.V%grid%nz+1)) then

          ! dealing with x-components
          if (SV%xyz(i) == 1) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%x(xi, yi, zi)

             ! dealing with y-component
          else if (SV%xyz(i) == 2) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%y(xi, yi, zi)

             ! dealing with z-component
          else if (SV%xyz(i) == 3) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%z(xi, yi, zi)
          end if

       else
          write(0,*) 'IJK out of bounds for dotProd_csvector_f'
          return
       endif

    enddo

  end function dotProd_noConj_csvector_f


  ! **********************************************************************
  ! compute sum cs*SV + V where SV is a complex sparse vector. V is a full
  ! complex vector of type cvector, and cs is a complex scalar; result
  ! overwrites V. Can also be used to construct a complex full vector
  ! from a sparse complex vector if cs = (1.0, 0.0) and V = 0 (initially)
  ! or copy from a sparse vector to a full vector
  subroutine add_scvector(cs,SV,V)

    implicit none
    type (sparsevecc), intent(in)	:: SV
    type (cvector), intent(inout)	:: V
    complex(kind=prec), intent(in)		:: cs
    integer				:: i
    integer				:: xi, yi, zi

    if((.not.SV%allocated).or.(.not.V%allocated)) then
       write(0,*) 'RHS not allocated yet for add_scvector'
       return
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for add_scvector'
       return
    endif

    ! loop over non-zero terms in sparse vector, adding to
    ! corresponding terms in full vector
    ! (need to check component xyz ...)
    do i = 1,SV%nCoeff

       ! generic test for both edge and face (all the components)
       if ((SV%i(i).le.V%grid%nx+1).or.(SV%j(i).le.V%grid%ny+1).or.&
            (SV%k(i).le.V%grid%nz+1)) then

          ! dealing with x-components
          if (SV%xyz(i) == 1) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             V%x(xi, yi, zi) = cs*SV%c(i) + V%x(xi, yi, zi)

             ! dealing with y-component
          else if (SV%xyz(i) == 2) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             V%y(xi, yi, zi) = cs*SV%c(i) + V%y(xi, yi, zi)

             ! dealing with z-component
          else if (SV%xyz(i) == 3) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             V%z(xi, yi, zi) = cs*SV%c(i) + V%z(xi, yi, zi)
          end if

       else
          write(0,*) 'IJK out of bounds for add_scvector'
          return
       endif

    enddo

  end subroutine add_scvector

  ! **********************************************************************
  ! this will conjugate a sparse complex vector SV1 (result SV2)
  ! A.K.
  function conjg_sparsevecc_f(SV1) result (SV2)

    type (sparsevecc), intent(in)		:: SV1
    type (sparsevecc)           		:: SV2

    ! check to see if SV1 is active (allocated)
    if(.not.SV1%allocated) then
       write(0,*) 'Input sparse vector not allocated yet for conjg_sparsevecc_f'
       return
    endif

    !  make sure SV2 is allocated and of the correct size
    if(SV2%allocated) then
       if ((SV2%gridType .ne. SV1%gridType) .or. (SV2%nCoeff .ne. SV1%nCoeff)) then
          call deall_sparsevecc(SV2)
          call create_sparsevecc(SV1%nCoeff,SV2,SV1%gridType)
       endif
    else
       call create_sparsevecc(SV1%nCoeff,SV2,SV1%gridType)
    endif

    SV2%i = SV1%i
    SV2%j = SV1%j
    SV2%k = SV1%k
    SV2%xyz = SV1%xyz
    SV2%c = conjg(SV1%c)

    SV2%temporary = .true.

  end function conjg_sparsevecc_f

  ! **************************************************************************
  ! subroutine to fill a single entry in the sparse vector with a value
  ! from full vector; sparsevecc has to be pre-allocated
  ! A.K.
  subroutine newValueC_sparsevecc(SV,index,c,i,j,k,xyz)

	type (sparsevecc), intent(inout)                :: SV
	complex(kind=prec), intent(in)                  :: c
	integer, intent(in)								:: i,j,k,xyz,index
	integer											:: istat

	if (.not.SV%allocated) then
	  write (0, *) 'Sparse vector in newValueC_sparsevecc is not allocated yet'
	  return
	end if

	if (index.gt.SV%nCoeff) then !ubound(SV%c)
	  write (0, *) 'The chosen index in newValueC_sparsevecc is not allocated yet'
	  return
	end if

	SV%i(index) = i
	SV%j(index) = j
	SV%k(index) = k
	SV%xyz(index) = xyz
	SV%c(index) = c

  end subroutine newValueC_sparsevecc


  ! **************************************************************************
  ! A subroutine to fill a single entry in the sparse vector with a value
  ! from full vector; sparsevecc has to be pre-allocated
  ! A.K.
  subroutine newValueR_sparsevecc(SV,index,r,i,j,k,xyz)

	type (sparsevecc), intent(inout)                :: SV
	real(kind=prec), intent(in)		                :: r
	integer, intent(in)								:: i,j,k,xyz,index
	integer											:: istat

	if (.not.SV%allocated) then
	  write (0, *) 'Sparse vector in newValueR_sparsevecc is not allocated yet'
	  return
	end if

	if (index.gt.SV%nCoeff) then
	  write (0, *) 'The chosen index in newValueR_sparsevecc is not allocated yet'
	  return
	end if

	SV%i(index) = i
	SV%j(index) = j
	SV%k(index) = k
	SV%xyz(index) = xyz
	SV%c(index) = dcmplx(r,0.0d0)

  end subroutine newValueR_sparsevecc

  ! **************************************************************************
  ! A subroutine to fill a single entry in the sparse vector with a value
  ! from full vector; sparsevecc has to be pre-allocated
  ! A.K.
  subroutine copyValue_csvector(SV,index,V,i,j,k,xyz)

	type (sparsevecc), intent(inout)                :: SV
	type (cvector), intent(in)                      :: V
	integer, intent(in)				:: i,j,k,xyz,index
	integer						:: istat

	if (.not.SV%allocated) then
	  write (0, *) 'Sparse vector in copyValue_csvector is not allocated yet'
	  return
	end if

	if (index.gt.SV%nCoeff) then
	  write (0, *) 'The chosen index in copyValue_csvector is not allocated yet'
	  return
	end if

	SV%i(index) = i
	SV%j(index) = j
	SV%k(index) = k
	SV%xyz(index) = xyz
	 if (xyz == 1) then
		SV%c(index) = V%x(i,j,k)
	 else if (xyz == 2) then
		SV%c(index) = V%y(i,j,k)
	 else if (xyz == 3) then
		SV%c(index) = V%z(i,j,k)
	 end if

  end subroutine copyValue_csvector

  ! **********************************************************************
  ! copy from a full vector to a sparse vector, this routine has quite a
  ! limited functionality as it assumes one knows the total number of
  ! non-zero vectors in full vector description
  subroutine copy_csvector(SV,V)

    implicit none
    type (sparsevecc), intent(inout)	:: SV
    type (cvector), target, intent(in)	:: V
    integer				:: i
    integer				:: xi, yi, zi

    if(.not.V%allocated) then
       write(0,*) 'RHS not allocated yet for copy_csvector'
       return
    endif

    if(.not.SV%allocated) then
       write(0,*) 'LHS not allocated yet for copy_csvector'
       return
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for copy_csvector'
       return
    endif

    if (V%gridType == EDGE) then
    i = 0
    ! for x - component
    do xi = 1, V%nx
       do yi = 1, V%ny+1
	  do zi = 1, V%nz+1
	     if (V%x(xi, yi, zi) /= 0.0) then
                i = i + 1
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   return
                end if
                SV%xyz(i) = 1
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%x(xi, yi, zi)
             end if
          end do
       end do
    end do
    ! for y - component
    do xi = 1, V%nx+1
       do yi = 1, V%ny
	  do zi = 1, V%nz+1
	     if (V%y(xi, yi, zi) /= 0.0) then
                i = i + 1
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   return
                end if
                SV%xyz(i) = 2
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%y(xi, yi, zi)
             end if
          end do
       end do
    end do
    ! for z - component
    do xi = 1, V%nx+1
       do yi = 1, V%ny+1
	  do zi = 1, V%nz
	     if (V%z(xi, yi, zi) /= 0.0) then
                i = i + 1
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   return
                end if
                SV%xyz(i) = 3
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%z(xi, yi, zi)
             end if
          end do
       end do
    end do

 else if (V%gridType == FACE) then
    ! for x - component
    do xi = 1, V%nx+1
       do yi = 1, V%ny
	  do zi = 1, V%nz
	     if (V%x(xi, yi, zi) /= 0.0) then
                i = i + 1
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   return
                end if
                SV%xyz(i) = 1
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%x(xi, yi, zi)
             end if
          end do
       end do
    end do
    ! for y - component
    do xi = 1, V%nx
       do yi = 1, V%ny+1
	  do zi = 1, V%nz
	     if (V%y(xi, yi, zi) /= 0.0) then
                i = i + 1
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   return
                end if
                SV%xyz(i) = 2
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%y(xi, yi, zi)
             end if
          end do
       end do
    end do
    ! for z - component
    do xi = 1, V%nx
       do yi = 1, V%ny
	  do zi = 1, V%nz+1
	     if (V%z(xi, yi, zi) /= 0.0) then
                i = i + 1
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   return
                end if
                SV%xyz(i) = 3
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%z(xi, yi, zi)
             end if
          end do
       end do
    end do

 else

    write (0, *) 'Vector (full) use not proper in copy_csvector'

 end if

end subroutine copy_csvector

end module sg_sparse_vector
