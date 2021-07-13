module SolnSpace
!   higher level module to define solution vector and RHS vector objects
!   for the forward solver; plus basic methods, linear algebra, dot products
!
! Defines: solnVector, sparseVector, rhsVector
! Uses: EMfield

use EMfieldInterp
use transmitters

implicit none

interface assignment (=)
   MODULE PROCEDURE copy_rhsVector
   MODULE PROCEDURE copy_solnVrhsV
   MODULE PROCEDURE copy_solnVrhsVMTX
   MODULE PROCEDURE copy_solnVector
   MODULE PROCEDURE copy_solnVectorMTX
   MODULE PROCEDURE copy_sparseVector
end interface

interface create
   MODULE PROCEDURE create_rhsVector
   MODULE PROCEDURE create_rhsVectorMTX
   MODULE PROCEDURE create_solnVector
   MODULE PROCEDURE create_solnVectorMTX
   MODULE PROCEDURE create_sparseVector
end interface

interface deall
   MODULE PROCEDURE deall_rhsVector
   MODULE PROCEDURE deall_rhsVectorMTX
   MODULE PROCEDURE deall_solnVector
   MODULE PROCEDURE deall_solnVectorMTX
   MODULE PROCEDURE deall_sparseVector
end interface

interface dotProd
   MODULE PROCEDURE dotProd_solnVector
   MODULE PROCEDURE dotProd_rhsVsolnV
   MODULE PROCEDURE dotProd_sparseVsolnV
end interface

 type :: solnVector_t
    !!   Generic solution type, same name must be used to allow
    !!   use of higher level inversion modules on different problems.
    !!
    !!   However, the specific implementation in this module
    !!    will depend on the specific problem.  This version is for 2D MT.
    !!   the basic solution object
    type(cvector)			:: vec

    !!  Mode (TE or TM) and omega/period are stored in the transmitter
    !!  dictionary. We avoid multiple potential problems by not duplicating
    !!  this information in EM solution (even when we have it available).
    !!  Then, if we need access to mode or frequency, we refer to the dictionary.
    type(grid_t), pointer		:: grid
    integer 				:: tx = 0
    logical                 :: allocated = .false.
  end type solnVector_t

  type :: solnVectorMTX_t
    !! Generic solution type for storing solutions from multiple transmitters
    integer			:: nTx = 0
    type(solnVector_t), pointer		:: solns(:)
    logical			:: allocated = .false.
    logical         :: temporary = .false.
  end type solnVectorMTX_t

  type :: sparseVector_t
    !!   Generic solution type, same name must be used to allow
    !!   use of higher level inversion modules on different problems.
    type(sparsevecc)			:: L
  end type sparseVector_t


  type :: rhsVector_t
     !!   right hand side for solving both TE and TM mode equations,
     !!   forward or adjoint problems (mode & omega stored in txDict)
     character*3			:: adj = ''
     logical           			:: nonzero_source = .false.
     logical				:: nonzero_bc = .false.
     logical				:: allocated = .false.
     type(cvector)			:: source
     complex(kind=prec), pointer, dimension(:)	:: bc
     type(grid_t), pointer		:: grid
     integer                :: tx = 0
  end type rhsVector_t

  type :: rhsVectorMTX_t
    !! Generic solution type for storing RHS's for multiple transmitters
    integer         :: nTx = 0
    type(rhsVector_t), pointer     :: combs(:)
    logical         :: allocated = .false.
  end type rhsVectorMTX_t

contains

!**********************************************************************
!           Basic solnVector methods
!**********************************************************************
     subroutine create_solnVector(grid,iTx,e)
       ! no need to pass the mode to create_solnVector:
       ! we now get it from the transmitter dict.
       ! now the interface is generic,
       ! and could be used as such

       implicit none
       type(grid_t), intent(in), target	:: grid
       integer, intent(in)              :: iTx
       type (solnVector_t), intent(inout)	:: e

       ! local
       character(2)      :: mode
       character(80)     :: gridType

       mode = txDict(iTx)%mode

       if(mode .eq. 'TE') then
          gridType = NODE
       else if(mode .eq. 'TM') then
          gridType = NODE_EARTH
       else
          call errStop('Unknown mode in create_solnVector')
       endif

       call create_cvector(grid,gridType,e%vec)
       e%tx = iTx
       e%grid => grid
       e%allocated = .true.

     end subroutine create_solnVector

     !************************************************************
     subroutine deall_solnVector(e)
       implicit none
       type (solnVector_t), intent(inout)   :: e

       call deall_cvector(e%vec)
       if(associated(e%grid)) then
           nullify(e%grid)
       endif
       e%allocated = .false.

     end subroutine deall_solnVector

     !************************************************************
     subroutine copy_solnVector(eOut,eIn)

       implicit none
       type (solnVector_t), intent(in)	:: eIn
       type (solnVector_t), intent(inout)	:: eOut

       !  should have some error checking for eIn ...
       call copy_cvector(eOut%vec,eIn%vec)

       eOut%tx = eIn%tx
       eOut%grid => eIn%grid
       eOut%allocated = eIn%allocated

     end subroutine copy_solnVector

     !**********************************************************************
     function dotProd_solnVector(e1,e2,conj_Case) result(c)
     !   dot product of two solution space objects
       type(solnVector_t), intent(in)		:: e1,e2
       logical, intent(in)		:: conj_case
       complex(kind=prec)	:: c

       c = dotProd_cvector(e1%vec,e2%vec,conj_case)

     end function dotProd_solnVector

     !**********************************************************************
     subroutine zero_solnVector(e)
     !  zeros a solution space object

       type(solnVector_t), intent(inout)	:: e

       e%vec%v = C_ZERO

     end subroutine zero_solnVector

   !**********************************************************************
   subroutine random_solnVector(e,eps)

      type(solnVector_t), intent(inout)     :: e
      real(kind=prec), intent(in), optional    :: eps

      !  local variables
      integer                           :: j, istat
      real(kind=prec), dimension(:,:), allocatable :: vr,vi

      if (.not. e%allocated) then
        call errStop('input solution vector not allocated in random_solnVector')
      end if

      allocate(vr(e%vec%N1,e%vec%N2),vi(e%vec%N1,e%vec%N2),STAT=istat)
      call random_number(vr)
      call random_number(vi)

        if (present(eps)) then
            e%vec%v = cmplx(vr,vi) * eps
        else
            e%vec%v = cmplx(vr,vi) * 0.05
        end if

      deallocate(vr,vi,STAT=istat)

   end subroutine random_solnVector

!**********************************************************************
!           Basic solnVectorMTX methods
!**********************************************************************

   subroutine create_solnVectorMTX(nTx,eAll)

      integer, intent(in)               :: nTx
      type(solnVectorMTX_t), intent(inout)  :: eAll

      ! local variables
      integer      :: istat

      if (eAll%allocated) then
         if (eAll%nTx == nTx) then
            ! do nothing
            return
         else
            call deall_solnVectorMTX(eAll)
         end if
      end if

      eAll%nTx = nTx
      allocate(eAll%solns(nTx), STAT=istat)
      eAll%allocated = .true.

   end subroutine create_solnVectorMTX

   !**********************************************************************
   subroutine deall_solnVectorMTX(eAll)

      type(solnVectorMTX_t), intent(inout)     :: eAll

      !  local variables
      integer                           :: j, istat

      if (eAll%allocated) then
	    do j = 1,eAll%nTx
	  	  call deall_solnVector(eAll%solns(j))
	    end do
	  end if

      if (associated(eAll%solns)) deallocate(eAll%solns, STAT=istat)
      eAll%nTx = 0
      eAll%allocated = .false.

   end subroutine deall_solnVectorMTX

   !**********************************************************************
   subroutine random_solnVectorMTX(eAll,eps)

      type(solnVectorMTX_t), intent(inout)     :: eAll
      real(kind=prec), intent(in), optional    :: eps

      !  local variables
      integer                           :: j, istat

      if (.not. eAll%allocated) then
        call errStop('input solution vector not allocated in random_solnVector')
      end if

      do j = 1,eAll%nTx
        if (present(eps)) then
            call random_solnVector(eAll%solns(j),eps)
        else
            call random_solnVector(eAll%solns(j),0.05*ONE)
        end if
      end do

   end subroutine random_solnVectorMTX

   !**********************************************************************
   subroutine copy_solnVectorMTX(eOut,eIn)

      type(solnVectorMTX_t), intent(inout)  :: eOut
      type(solnVectorMTX_t), intent(in)     :: eIn

      !  local variables
      integer                   :: j, istat

      if (eOut%allocated) then
        call deall_solnVectorMTX(eOut)
      end if

      call create(eIn%nTx,eOut)
      do j = 1,eIn%nTx
        if (eIn%solns(j)%allocated) then
            eOut%solns(j) = eIn%solns(j)
        else
            eOut%solns(j)%allocated = .false.
        end if
      end do
      eOut%allocated = .true.

      !if (eIn%temporary) then
      !  call deall_solnVectorMTX(eIn)
      !end if

   end subroutine copy_solnVectorMTX

!**********************************************************************
!           Basic sparseVector methods
!**********************************************************************

     subroutine create_sparseVector(grid,iTx,LC,nCoeff)

       type(grid_t), intent(in)	:: grid
       integer, intent(in)		:: iTx
       type(sparseVector_t)			:: LC
       integer, optional, intent(in)	:: nCoeff

       ! local
       character(2)      :: mode
       character(80)     :: gridType
       integer			 :: nc

       if (present(nCoeff)) then
          nc = nCoeff
       else
          nc = 0 ! will reallocate memory later in the program
       end if

       mode = txDict(iTx)%mode

       if(mode .eq. 'TE') then
          gridType = NODE
       else if(mode .eq. 'TM') then
          gridType = NODE_EARTH
       else
          call errStop('Unknown mode in create_sparseVector')
       endif

       call create_sparsevecc(grid,gridType,nc,LC%L)

     end subroutine create_sparseVector

     !************************************************************
     subroutine deall_sparseVector(LC)

       type (sparseVector_t), intent(inout)   	:: LC

       call deall_sparsevecc(LC%L)

     end subroutine deall_sparseVector

   !************************************************************
   subroutine copy_sparseVector(eOut,eIn)

       !  2D  version
       implicit none
       type (sparseVector_t), intent(in)    :: eIn
       type (sparseVector_t), intent(inout) :: eOut

       ! local variables
       integer              :: k

       call copy_sparsevecc(eOut%L,eIn%L)

   end subroutine copy_sparseVector

     !**********************************************************************
     subroutine linComb_sparseVector(Lin1,c1,Lin2,c2,Lout)
       ! linear combination of two sparseVector objects

       type (sparseVector_t), intent(in)		:: Lin1,Lin2
       complex (kind=prec), intent(in)	:: c1,c2
       type (sparseVector_t), intent(inout)		:: Lout

       call linComb_sparsevecc(Lin1%L,c1,Lin2%L,c2,Lout%L)

     end subroutine linComb_sparseVector

!**********************************************************************
!           combined solnVector/sparseVector methods
!**********************************************************************
     function dotProd_sparseVsolnV(SV,FV,Conj_Case) result(c)

       type (sparseVector_t), intent(in)             :: SV  ! sparse vector
       type (solnVector_t), intent(in)               :: FV  ! full vector
       logical, intent(in)                     :: conj_Case ! = .true.
       complex(kind=prec)		:: c

       c = dotProd_scvector(SV%L,FV%vec,Conj_Case)

     end function dotProd_sparseVsolnV

     !**********************************************************************
     subroutine add_sparseVsolnV(cs,SV,FV)

        complex(kind=prec), intent(in)	:: cs
        type (sparseVector_t), intent(in)	:: SV  ! sparse vector
        type (solnVector_t), intent(inout)	:: FV  ! full vector

       call add_scvector(cs,SV%L,FV%vec)

     end subroutine add_sparseVsolnV

!**********************************************************************
!           rhsVector methods
!**********************************************************************
     !  allocates and initializes arrays for the "rhs" structure
     !   set pointer to grid + mode/source fields of rhs before calling
     subroutine create_rhsVector(grid,iTx,b)

       type(grid_t), intent(in),target	:: grid
       integer, intent(in)              :: iTx
       type (rhsVector_t), intent(inout)   	:: b

       !  local variables
       character(2)     :: mode
       character(80)    :: gridType
       integer ::       Nz,Ny,Nza,Nzi,nBC

       Nz = grid%Nz
       Ny = grid%Ny
       Nza = grid%Nza

       mode = txDict(iTx)%mode

       if(mode .eq. 'TE') then
          gridType = NODE
          Nzi = Nz
       else if(mode .eq. 'TM') then
          gridType = NODE_EARTH
          Nzi = Nz-Nza
       else
          call errStop('Unknown mode in create_rhsVector')
       endif

       nBC = 2*(Ny+1)+2*(Nzi+1)

       if (b%nonzero_BC) then
          ! Initialize array for boundary condition data
          allocate(b%BC(nBC))
       endif

       if (b%nonzero_source) then
           ! Initialize full array for storage of source
           call create_cvector(grid,gridType,b%source)
       endif

       b%tx = iTx
       b%grid => grid
       b%allocated = .true.

     end subroutine create_rhsVector

    !************************************************************
     subroutine deall_rhsVector(b)

       type (rhsVector_t), intent(inout)   :: b

       if (.not.(b%allocated)) then
         return
       endif

       if (b%nonzero_BC) then
          deallocate(b%BC)
          nullify(b%BC)
       endif

       if (b%nonzero_source) then
          call deall_cvector(b%source)
       endif
       b%allocated = .false.

     end subroutine deall_rhsVector

     !**********************************************************************
     subroutine add_sparseVrhsV(cs,SV,FV)

        complex(kind=prec), intent(in)	:: cs
        type (sparseVector_t), intent(in)		:: SV  ! sparse vector
        type (rhsVector_t), intent(inout)		:: FV  ! full vector

       call add_scvector(cs,SV%L,FV%source)
       FV%nonzero_source = .true.

     end subroutine add_sparseVrhsV

     !************************************************************
     subroutine copy_rhsVector(bOut,bIn)

       !  2D  version
       implicit none
       type (rhsVector_t), intent(in)   :: bIn
       type (rhsVector_t), intent(inout)    :: bOut

       ! local variables
       integer              :: k

       if (.not. bIn%allocated) then
         call errStop('input EM RHS not allocated yet in copy_rhsVector')
       endif

       bOut%adj = bIn%adj
       bOut%nonzero_Source = bIn%nonzero_Source
       bOut%nonzero_BC = bIn%nonzero_BC
       call create_rhsVector(bIn%grid,bIn%tx,bOut)

       if (bOut%nonzero_Source) then
        ! source cvector already created
        call copy_cvector(bOut%source,bIn%source)
       end if

       if (bOut%nonzero_BC) then
        ! BC vector already created
        bOut%bc = bIn%bc
       end if

       bOut%allocated = bIn%allocated

       !if (bIn%temporary) then
       !   call deall_rhsVector(bIn)
       !endif

     end subroutine copy_rhsVector

     !**********************************************************************
     subroutine zero_rhsVector(e)
     !  zeros a solution space object

       type(rhsVector_t), intent(inout)	:: e

       if(e%nonzero_source .and. e%allocated) then
          e%source%v = C_ZERO
       else
          if(e%nonzero_bc .and. e%allocated) then
             e%bc = C_ZERO
          endif
       endif
     end subroutine zero_rhsVector

   !**********************************************************************
   subroutine random_rhsVector(e,eps)

      type(rhsVector_t), intent(inout)     :: e
      real(kind=prec), intent(in), optional    :: eps

      !  local variables
      integer                           :: j, istat
      real(kind=prec), dimension(:,:), allocatable :: vr,vi
      real(kind=prec), dimension(:), allocatable   :: br,bi

      if (e%nonzero_source .and. e%allocated) then
          allocate(vr(e%source%N1,e%source%N2),vi(e%source%N1,e%source%N2),STAT=istat)
          call random_number(vr)
          call random_number(vi)

          if (present(eps)) then
              e%source%v = cmplx(vr,vi) * eps
          else
              e%source%v = cmplx(vr,vi) * 0.05
          end if

          deallocate(vr,vi,STAT=istat)
      end if

      if(e%nonzero_bc .and. e%allocated) then
          allocate(br(size(e%bc)),bi(size(e%bc)),STAT=istat)
          call random_number(br)
          call random_number(bi)
          if (present(eps)) then
              e%bc = cmplx(br,bi) * eps
          else
              e%bc = cmplx(br,bi) * 0.05
          end if

          deallocate(br,bi,STAT=istat)
      endif

   end subroutine random_rhsVector

     !**********************************************************************

     subroutine sparse2full_rhsVector(b)

     !   Converts an rhsVector object from sparse to full representation,
     !   by consolidating bc into the full cvector source.
     !   The order of BCs in WSfwd2D is: left, surface, bottom, right.
     ! CURRENTLY THIS IS NOT WORKING SO STEST CANNOT BE RUN. FIX THIS LATER.
     ! AK 24 May 2018

       type (rhsVector_t), intent(inout)             :: b

       !  local variables
       character(2)     :: mode
       character(80)    :: gridType
       integer ::       Nz,Ny,Nza,Nzi,nBC,k1,k2,istat

       Nz = b%grid%Nz
       Ny = b%grid%Ny
       Nza = b%grid%Nza

       mode = txDict(b%tx)%mode

       if(mode .eq. 'TE') then
          gridType = NODE
          Nzi = Nz
       else if(mode .eq. 'TM') then
          gridType = NODE_EARTH
          Nzi = Nz-Nza
       else
          call errStop('Unknown mode in sparse2full_rhsVector')
       endif

       nBC = 2*(Ny+1)+2*(Nzi+1)

       if (.not. b%nonzero_source) then
           ! Initialize full array for storage of source
           call create_cvector(b%grid,gridType,b%source)
       endif

       if (b%nonzero_BC) then
           ! left
           k1 = 1
           k2 = b%source%N2
           b%source%v(1,1:b%source%N2) = b%BC(k1:k2)
           ! surface
           k1 = b%source%N2+1
           k2 = b%source%N2+b%source%N1
           b%source%v(1:b%source%N1,1) = b%BC(k1:k2)
           ! bottom
           k1 = b%source%N2+b%source%N1+1
           k1 = b%source%N2+2*b%source%N1
           b%source%v(1:b%source%N1,b%source%N2) = b%BC(k1:k2)
           ! right
           k1 = b%source%N2+2*b%source%N1+1
           k1 = 2*b%source%N1+2*b%source%N2
           b%source%v(b%source%N1,1:b%source%N2) = b%BC(k1:k2)
       endif

       deallocate(b%BC, STAT=istat)

       b%nonzero_Source = .true.
       b%nonzero_BC = .false.


     end subroutine sparse2full_rhsVector
!**********************************************************************
!           combined solnVector/rhsVector methods
!**********************************************************************
     !**********************************************************************
     subroutine copy_solnVrhsV(b,e)
     !  implements b = e

       type(rhsVector_t), intent(inout) :: b
       type(solnVector_t), intent(in)   :: e

       if (.not. e%allocated) then
         call errStop('input EM soln not allocated yet in copy_solnVrhsV')
       endif

       call create_rhsVector(e%grid,e%tx,b)

       b%nonzero_source = .true.
       b%source = e%vec
       b%nonzero_BC = .false.
       b%allocated = .true.

     end subroutine copy_solnVrhsV

     !**********************************************************************
     function dotProd_rhsVsolnV(comb,FV,Conj_Case) result(c)
       ! computes a dot product between RHS vector and solution vector;
       ! does not compute anything from the boundary conditions!!!
       ! this might have to be changed.

       type (rhsVector_t), intent(in)             :: comb  ! rhs
       type (solnVector_t), intent(in)            :: FV  ! full vector
       logical, intent(in)                     :: conj_Case ! = .true.
       complex(kind=prec)       :: c

       if((.not. comb%allocated) .or. (.not. FV%allocated)) then
            call errStop('RHS and solution vectors have to be allocated in dotProd_rhsVsolnV')
       endif

       if(comb%nonzero_source) then
            c = dotProd_cvector(comb%source,FV%vec,Conj_Case)
       else
            c = C_ZERO
       endif

     end function dotProd_rhsVsolnV

!**********************************************************************
!           Basic rhsVectorMTX methods
!**********************************************************************

   subroutine create_rhsVectorMTX(nTx,bAll)

      integer, intent(in)               :: nTx
      type(rhsVectorMTX_t), intent(inout)  :: bAll

      !  local variables
      integer                           :: istat

      if (bAll%allocated) then
         if (bAll%nTx == nTx) then
            ! do nothing
            return
         else
            call deall_rhsVectorMTX(bAll)
         end if
      end if

      bAll%nTx = nTx
      allocate(bAll%combs(nTx), STAT=istat)
      bAll%allocated = .true.

   end subroutine create_rhsVectorMTX

   !**********************************************************************
   subroutine deall_rhsVectorMTX(bAll)

      type(rhsVectorMTX_t), intent(inout)     :: bAll

      !  local variables
      integer                           :: j, istat

      if (bAll%allocated) then
        do j = 1,bAll%nTx
          call deall_rhsVector(bAll%combs(j))
        end do
      end if

      if (associated(bAll%combs)) deallocate(bAll%combs, STAT=istat)
      bAll%nTx = 0
      bAll%allocated = .false.

   end subroutine deall_rhsVectorMTX

   !**********************************************************************
   subroutine random_rhsVectorMTX(bAll,eps)

      type(rhsVectorMTX_t), intent(inout)      :: bAll
      real(kind=prec), intent(in), optional    :: eps

      !  local variables
      integer                           :: j, istat

      if (.not. bAll%allocated) then
        call errStop('input RHS vector not allocated in random_rhsVectorMTX')
      end if

      do j = 1,bAll%nTx
        if (present(eps)) then
            call random_rhsVector(bAll%combs(j),eps)
        else
            call random_rhsVector(bAll%combs(j),0.05*ONE)
        end if
      end do

   end subroutine random_rhsVectorMTX

   !**********************************************************************
   subroutine copy_solnVrhsVMTX(eOut,eIn)

      type(rhsVectorMTX_t), intent(inout)  :: eOut
      type(solnVectorMTX_t), intent(in)     :: eIn

      !  local variables
      integer                   :: j, istat

      if (eOut%allocated) then
        call deall_rhsVectorMTX(eOut)
      end if

      call create(eIn%nTx,eOut)
      do j = 1,eIn%nTx
        eOut%combs(j)%nonzero_source = .true.
        eOut%combs(j) = eIn%solns(j)
      end do
      eOut%allocated = .true.

      !if (eIn%temporary) then
      !  call deall_solnVectorMTX(eIn)
      !end if

   end subroutine copy_solnVrhsVMTX

end module SolnSpace
