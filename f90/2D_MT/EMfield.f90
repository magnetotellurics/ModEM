module EMfield
!   module for basic objects used to represent solutions to
!   the discretized equations, including basic create/destroy
!   routines, vector space arithmetic, dot products
!   Also contains routines for sparse vectors, allowing for
!   similar vector space operations (including operations on
!   full and sparse vectors together)

!  Sparse vectors are used for representation of data functionals,
!  and only complex data types are supported.
! Only those routines needed for modeling and initial work on 2D
!   inversion have been developed here (less than for 3D).
!  NOTE that for 2D we only deal with scalar fields, defined either
!   on nodes or cells ... these are distinguished through
!   structure field gridType, which can be either NODE or CELL
!  For full-sparse operations full vectors are stored as type
!   Vec2D ... gridTypes should match!

use math_constants
use utilities
use griddef
implicit none

 type :: cvector
    !!   generic 2D array for storing solution vectors, etc

    logical           :: allocated = .false.
    ! actual size of array
    integer   :: N1=0
    integer   :: N2=0
    complex(kind=prec), pointer, dimension(:,:)      :: v
    type(grid_t), pointer		:: grid
    !   corners, cells, sides ... full ...interior, whatever
    !   supported types at present:
    !         CELL, NODE, CELL_EARTH, NODE_EARTH, EDGE_EARTH
    character*80			:: gridType = ''
  end type cvector

  type :: sparsevecc

     ! complex sparse vector for 2D EM -- represents scalar field
     character (len=80)			:: gridType=''
     ! nCoeff is number of non-zero nodes
     integer				:: nCoeff  = 0
     ! j,k are arrays of indices that defines grid location
     integer , pointer, dimension(:)	:: j,k
     complex (kind=prec), pointer, dimension(:)      :: c
     ! has sparse vector been allocated?
     logical				:: allocated = .false.
     ! pointer to the parent grid
     type (grid_t), pointer		:: grid

  end type sparsevecc

  !  full storage vector routines
  public	:: create_cvector,deall_cvector
  public	:: dotProd_cvector,zero_cvector,copy_cvector

  !  sparse storage vector routines
  public	:: create_sparsevecc, deall_sparsevecc
  public	:: linComb_sparsevecc
  public	:: copy_sparsevecc, scMult_sparsevecc

  !  combined full/sparse storage vector routines
  public	:: dotProd_scvector
  public	:: add_scvector

  interface assignment (=)
     MODULE PROCEDURE copy_cvector
     MODULE PROCEDURE copy_sparsevecc
  end interface

contains
   !  should probably add : pointwise multiplication,
   !    linear combinations ...

   !************************************************************************
     !  create_cvector allocates and initializes arrays for
     !   Earth-cell conductivity structure;
     !   Pass grid of type grid_t to set array sizes
     subroutine create_cvector(grid,gridType,vec)

       implicit none
       type (grid_t), intent(in), target	:: grid
       character*80				:: gridType
       type (cvector), intent(inout)		:: vec
       !  local variables
       integer 					::Nz,Ny,Nza,NzEarth
       integer                  ::istat

       Nz = grid%Nz
       Nza = grid%Nza
       NzEarth = Nz-Nza
       Ny = grid%Ny

       select case (gridType)
          case (NODE)
            vec%N1 = Ny+1
            vec%N2 = Nz+1
            allocate(vec%v(Ny+1,Nz+1),STAT=istat)
          case (NODE_EARTH)
            vec%N1 = Ny+1
            vec%N2 = NzEarth+1
            allocate(vec%v(Ny+1,NzEarth+1),STAT=istat)
          case (CELL)
            vec%N1 = Ny
            vec%N2 = Nz
            allocate(vec%v(Ny,Nz),STAT=istat)
          case (CELL_EARTH)
            vec%N1 = Ny
            vec%N2 = NzEarth
            allocate(vec%v(Ny,NzEarth),STAT=istat)
          case (EDGE_EARTH)
            ! allocates for all edges, including boundaries
            vec%N1 = Ny+1
            vec%N2 = NzEarth+1
            allocate(vec%v(Ny+1,NzEarth+1),STAT=istat)
          case default
             write(0,*) 'This gridType not coded'
             return
       end select
       vec%v = C_ZERO
       vec%gridType = gridType
       vec%allocated = .true.
       vec%grid => grid
     end subroutine create_cvector

    !************************************************************
     subroutine deall_cvector(vec)
       implicit none
       type (cvector), intent(inout)   :: vec
       integer istat
       if (vec%allocated) then
       if (associated(vec%v)) deallocate(vec%v,STAT=istat)
       vec%allocated = .false.
       vec%gridType = ''
       nullify(vec%grid)
      end if
     end subroutine deall_cvector

   !**********************************************************************
   function dotProd_cvector(e1,e2,conj_Case) result(c)

   !   dot product of two model space (conductivity parameter)
   !    objects; implementation for
     type(cvector), intent(in)		:: e1,e2
     logical, intent(in)		:: conj_case
     complex(kind=prec)	:: c

     ! local variables
     integer				:: j,k

     if((e1%N1 .ne. e2%N1).or. (e1%N2 .ne. e2%N2)) then
        call errStop('size of m1, m2 incompatable in EarthCondDotProd')
     endif

     c = C_ZERO
     if(conj_case) then
        do k = 1,e1%N2
           do j = 1,e1%N1
              c = c + conjg(e1%v(j,k))*e2%v(j,k)
           enddo
        enddo
     else
        do k = 1,e1%N2
           do j = 1,e1%N1
              c = c + e1%v(j,k)*e2%v(j,k)
           enddo
        enddo
     endif

   end function dotProd_cvector

   !**********************************************************************
   subroutine zero_cvector(e)

     type(cvector), intent(inout)	:: e

     e%v = C_ZERO

   end subroutine zero_cvector

   !**********************************************************************
   subroutine copy_cvector(eOut,eIn)

     ! copies eIn to eOut

     type(cvector), intent(in)		:: eIn
     type(cvector), intent(inout)	:: eOut

     ! make sure eOut is allocated, and of same size; if not same
     !   size, deallocate and allocate as correct size; otherwise
     !   do nothing (use m as input ... allowing m to overwrite m1)
     if(eOut%allocated) then
        if((eIn%N1 .ne. eOut%N1).or. (eIn%N2 .ne. eOut%N2)) then
           call deall_cvector(eOut)
           call create_cvector(eIn%grid,eIn%gridType,eOut)
        endif
     else
        call create_cvector(eIn%grid,eIn%gridTYpe,eOut)
     endif
     eOut%v = eIn%v

   end subroutine copy_cvector

!   Sparse vector routines
! **********************************************************************
  ! delete/ deallocate the sparse vector
  subroutine deall_sparsevecc(oldLC)

    implicit none
    type (sparsevecc), intent(inout)		:: oldLC
    integer					:: status

    if(oldLC%allocated) then
       deallocate(oldLC%j, STAT=status)
       deallocate(oldLC%k, STAT=status)
       deallocate(oldLC%c, STAT=status)
       nullify(oldLC%j,oldLC%k,oldLC%c)
       if(associated(oldLC%grid)) then
          nullify(oldLC%grid)
       endif
       oldLC%ncoeff = 0
       oldLC%gridType = ''
       oldLC%allocated = .false.
    end if

  end subroutine deall_sparsevecc

  ! **********************************************************************
  ! create an object of type sparsevecc of length nCoeff
  subroutine create_sparsevecc(grid,gridType,nCoeff,newLC)

    implicit none
    type(grid_t), intent(in), target	:: grid
    character (len=80), intent(in)     	:: gridType
    integer, intent(in) 		:: nCoeff
    type (sparsevecc), intent(inout) 	:: newLC
    integer				:: status

    if(newLC%allocated) then
       call deall_sparsevecc(newLC)
    endif

    newLC%allocated = .true.
    allocate(newLC%j(nCoeff),STAT=status)
    newLC%allocated = newLC%allocated .and. (status .eq. 0)
    allocate(newLC%k(nCoeff),STAT=status)
    newLC%allocated = newLC%allocated .and. (status .eq. 0)
    allocate(newLC%c(nCoeff),STAT=status)
    newLC%allocated = newLC%allocated .and. (status .eq. 0)

    newLC%nCoeff = nCoeff
    newLC%j = 0
    newLC%k = 0
    newLC%c = C_ZERO
    newLC%gridType = gridType
    newLC%grid => grid

  end subroutine create_sparsevecc
  ! **********************************************************************
 ! linear combination of two sparse vectors, output as a sparse vector
  ! allocates (or reallocates) output sparse vector Loc3
  subroutine linComb_sparsevecc(Lic1,c1,Lic2,c2,Loc3)

    implicit none
    type (sparsevecc), target, intent(in)     :: Lic1,Lic2
    type (sparsevecc), intent(inout)          :: Loc3
    complex (kind=prec), intent(in)	:: c1,c2

    integer                                     :: n,m,nm
    integer                                     :: nCoeffSum
    integer, allocatable, dimension(:)          :: Lic1oldIndex
    integer                                     :: status

    allocate(Lic1oldIndex(Lic2%nCoeff), STAT = status)

    if(Loc3%allocated) then
       call deall_sparsevecc(Loc3)
    endif

    if (Lic1%gridType == Lic2%gridType) then
       ! count common indices
       nCoeffSum = Lic1%nCoeff+Lic2%nCoeff
       do m = 1,Lic2%nCoeff
          Lic1oldIndex(m) = 0
          do n = 1,Lic1%nCoeff
             if((Lic1%j(n).eq.Lic2%j(m)).and.(Lic1%k(n).eq.Lic2%k(m)))then
                nCoeffSum = nCoeffSum-1
                Lic1oldIndex(m) = n
                exit
             endif
          enddo
       enddo
    else
       write (0, *) 'gridType incompatability in linComb_sparsevecc'
       stop
    end if

    call create_sparsevecc(Lic1%grid,Lic1%gridType,nCoeffsum,Loc3)
    nm = Lic1%nCoeff
    Loc3%j(1:nm) = Lic1%j
    Loc3%k(1:nm) = Lic1%k
    Loc3%c(1:nm) = c1*Lic1%c

    do m = 1,Lic2%nCoeff
       ! if no indices are common, just concatenate
       if(Lic1oldIndex(m).eq.0) then
          nm = nm+1
          Loc3%j(nm) = Lic2%j(m)
          Loc3%k(nm) = Lic2%k(m)
          Loc3%c(nm) = c2*Lic2%c(m)
       else
          Loc3%c(Lic1oldIndex(m)) =  Loc3%c(Lic1oldIndex(m))+c2*Lic2%c(m)
      endif
    enddo
  end subroutine linComb_sparsevecc

  ! **********************************************************************
  ! multiplication of a sparse vector by a scalar, output as a sparse vector
  ! allocates (or reallocates) output sparse vector Loc
  subroutine scMult_sparsevecc(cs,Lic,Loc)

    implicit none
    complex (kind=prec), intent(in)	    :: cs
    type (sparsevecc), target, intent(in)     :: Lic
    type (sparsevecc), intent(inout)          :: Loc

    integer                                     :: n,m,nm
    integer                                     :: nCoeffSum
    integer                                     :: status

    if(Loc%allocated) then
       call deall_sparsevecc(Loc)
    endif

    call create_sparsevecc(Lic%grid,Lic%gridType,Lic%nCoeff,Loc)
    nm = Lic%nCoeff
    Loc%j(1:nm) = Lic%j
    Loc%k(1:nm) = Lic%k
    Loc%c(1:nm) = cs*Lic%c

  end subroutine scMult_sparsevecc

  ! **********************************************************************
  ! copy a sparse 2D complex vector from SV1 to SV2 ...
  ! 		can interface to =
  ! If necessary deallocate and reinitialize SV2
  subroutine copy_sparsevecc(SV2,SV1)
    implicit none
    !  INPUT is second argument SV1
    type (sparsevecc), target, intent(in)       :: SV1
    !  OUTPUT is second argument SV1
    type (sparsevecc), intent(inout)            :: SV2

    if(SV2%allocated) then
       call deall_sparsevecc(SV2)
    endif

    if(SV1%allocated) then
       call create_sparsevecc(SV1%grid,SV1%gridType,SV1%nCoeff,SV2)
       SV2%j = SV1%j
       SV2%k = SV1%k
       SV2%c = SV1%c
    end if
  end subroutine copy_sparsevecc

  ! **********************************************************************
  ! compute complex dot product between a sparse vector SV and a full
  !   storage vector of type cvector
  function dotProd_scvector(SV,FV,Conj_Case) result(c)

    implicit none
    type (sparsevecc), intent(in)		:: SV  ! sparse vector
    type (cvector), intent(in)			:: FV  ! full vector
    logical, intent(in)				:: conj_Case ! = .true.
				! for standard Hermitian inner product
    complex(kind=prec)			:: c
    integer					:: i
    integer					:: yi, zi

    if(FV%gridType == SV%gridType) then
       c = C_ZERO
       ! sum over  non-zero terms in sparse vector (conjugate sparse)
       do i = 1,SV%nCoeff
          yi = SV%j(i)
          zi = SV%k(i)
          if ((yi.gt.0).and.(yi.le.FV%N1).and. &
		(zi.gt.0).and.(zi.le.FV%N2)) then
             if(Conj_Case) then
                c = c + conjg(SV%c(i)) * FV%v(yi,zi)
             else
                c = c + SV%c(i) * FV%v(yi,zi)
             endif
          else
             write(0,*) 'indices for sparse-full dot product out of bounds '
             stop
          endif
       enddo
    else
       write (0, *) 'gridType incompatibility in dotProd_scvector'
       stop
    endif

  end function dotProd_scvector
  ! **********************************************************************
  ! compute sum cs*SV + V where SV is a complex sparse vector. V is
  ! a complex 2D array of type cvector, and cs is a complex scalar; result
  ! overwrites V. Can also be used to construct a complex full vector
  ! from a sparse complex vector if cs = (1.0, 0.0) and V = 0 (initially)
  ! or copy from a sparse vector to a full vector
  subroutine add_scvector(cs,SV,FV)

    implicit none
    complex(kind=prec), intent(in)	:: cs
    type (sparsevecc), intent(in)		:: SV  ! sparse vector
    type (cvector), intent(inout)		:: FV  ! full vector
    integer					:: i
    integer					:: yi, zi

    if(FV%gridType == SV%gridType) then
       ! loop over non-zero terms in sparse vector, adding to
       ! corresponding terms in full vector
       do i = 1,SV%nCoeff
          yi = SV%j(i)
          zi = SV%k(i)
          if ((yi.gt.0).and.(yi.le.FV%N1).and. &
		(zi.gt.0).and.(zi.le.FV%N2)) then
             FV%v(yi,zi) = cs*SV%c(i) + FV%v(yi,zi)
          else
             write(0,*) 'indices for sparse-full sum out of bounds'
             stop
          endif
       enddo
    else
       write (0, *) 'gridType incompatability in SparseFullDotProd2D'
       stop
    endif

  end subroutine add_scvector

end module EMfield
