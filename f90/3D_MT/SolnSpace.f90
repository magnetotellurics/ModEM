module SolnSpace
!   higher level module to define EM solution and RHS objects
!   plus basic methods, linear algebra, dot products
!   3D MT version
!
! Defines: solnVector, sparseVector, rhsVector
! Uses: EMfield

use math_constants
use utilities
use sg_vector
use sg_boundary
use sg_sparse_vector
use transmitters

implicit none

interface assignment (=)
   !MODULE PROCEDURE copy_rhsVector - doesn't exist yet
   MODULE PROCEDURE copy_solnVrhsV
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
    !!    will depend on the specific problem.  This version is for 3D MT,
    !!  and collects solutions for both "polarizations" in a single
    !!   structure.  Thus the impedance operator acts on an object of this type.

    !! e is array of electric field solution vectors for nPol source
    !! polarizations; 2 for 3D MT - one electrical field solution for each
    !! allows a variable number of polarizations to accommodate other uses
    !! e.g. active source applications
    integer					:: nPol = 2
    integer                 :: Pol_index(2)
    type(cvector), pointer  :: pol(:)

    !! tx points to information in the transmitter dictionary about the source
    !!   used to compute the solution, e.g. omega/period;
    !!   do not duplicate it here to avoid potential problems
    integer 			:: tx = 0

    !! grid is a pointer to numerical discretization stored in SensMatrix
    type(grid_t), pointer	:: grid

    !! allocated when the solnVector was created but not yet deallocated
    logical			:: allocated = .false.

    !! avoid memory leaks: set this to true for function outputs only
    logical			:: temporary = .false.

  end type solnVector_t

  type :: solnVectorMTX_t
    !! Generic solution type for storing solutions from multiple transmitters
    integer						:: nTx = 0
    type(solnVector_t), pointer		:: solns(:)
    logical						:: allocated = .false.
    logical						:: temporary = .false.
  end type solnVectorMTX_t

  type :: sparseVector_t
    !!   Generic solution type, same name must be used to allow
    !!   use of higher level inversion modules on different problems.
    !!   Here we need nPol sparse vectors, one for each polarization
    !!    to represent linear data functionals on objects of type solnVector
    !!  this doesn't really need the grid... not used anywhere at the moment
    !!  but then the generic create interface shouldn't use it either
    !!  unless sparse vectors require the grid pointer for something.
    integer						:: nPol = 2
    integer						:: nCoeff = 0
    type(sparsevecc), pointer	:: L(:)
    logical						:: allocated = .false.
    logical						:: temporary = .false.
    type(grid_t), pointer       :: grid
    integer						:: tx = 0
  end type sparseVector_t

  type :: RHS_t
     ! merges internal sources and boundary conditions into a single
     ! data structure; this is a compact representation of the right
     !  hand side for the induction equations for ONE mode

     character*3		:: adj = ''
     character*2		:: polName ! Ex or Ey
     logical                    :: nonzero_BC     = .false.
     logical                    :: nonzero_Source = .false.
     logical                    :: sparse_Source  = .false.
     logical			:: allocated      = .false.
    logical					:: temporary = .false.
     type (cvector) 		:: s
     type (sparsevecc) 		:: sSparse
     type (cboundary) 		:: bc
     type(grid_t), pointer	:: grid
     integer                :: tx = 0
  end type RHS_t

  type :: rhsVector_t
     ! rhs data structure for multiple polarizations, the abstract
     !  full rhs used for the abstract full solnVector
     ! pointer to the grid needed for cleaner routines in this module.

     integer				:: nPol = 2
     type (RHS_t), pointer	:: b(:)
     logical                    :: nonzero_BC     = .false.
     logical                    :: nonzero_Source = .false.
     logical                    :: sparse_Source  = .false.
     logical				:: allocated = .false.
     logical				:: temporary = .false.
     type(grid_t), pointer	:: grid
     integer				:: tx = 0
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

     !  generic routine for creating the solnVector type for 3D problems:
     !  number of polarization obtained from the transmitter dictionary

       implicit none
       type(grid_t), intent(in), target	    :: grid
       integer, intent(in)                  :: iTx
       type (solnVector_t), intent(inout)		:: e

       ! local variables
       integer				:: k,istat,iPol

       if (e%allocated) then
          if (associated(e%grid, target=grid) .and. (e%tx == iTx)) then
             ! do nothing
             return
          else
             call deall_solnVector(e)
          end if
       end if

       e%nPol = txDict(iTx)%nPol
       
       do iPol=1,e%nPol
        e%Pol_index(iPol)=iPol
       end do
       
       allocate(e%pol(e%nPol), STAT=istat)
       do k = 1,e%nPol
          call create_cvector(grid,e%pol(k),EDGE)
       enddo
       e%tx = iTx
       e%grid => grid

	   e%allocated = .true.

     end subroutine create_solnVector

     !************************************************************
     subroutine deall_solnVector(e)

       !  3D  version
       implicit none
       type (solnVector_t), intent(inout)   :: e

       ! local variables
       integer				:: k, istat

       if (associated(e%pol)) then
          do k = 1,e%nPol
             call deall_cvector(e%pol(k))
          enddo
          deallocate(e%pol, STAT=istat)
       endif

       if(associated(e%grid)) then
           nullify(e%grid)
       endif

       e%allocated = .false.

     end subroutine deall_solnVector

     !************************************************************
     subroutine copy_solnVector(eOut,eIn)

       !  3D  version
       implicit none
       type (solnVector_t), intent(in)	:: eIn
       type (solnVector_t), intent(inout)	:: eOut

       ! local variables
       integer				:: k

       if (.not. eIn%allocated) then
         call errStop('input EM soln not allocated yet in copy_solnVector')
       endif

       call create_solnVector(eIn%grid,eIn%tx,eOut)

       do k = 1,eIn%nPol
          call copy_cvector(eOut%pol(k),eIn%pol(k))
       enddo

       eOut%allocated = eIn%allocated

       !if (eIn%temporary) then
       !   call deall_solnVector(eIn)
       !endif

     end subroutine copy_solnVector

     !**********************************************************************
     subroutine zero_solnVector(e)
     !  zeros a solution space object

       type(solnVector_t), intent(inout)	:: e

       ! local variables
       integer				:: k

       do k = 1,e%nPol
          call zero_cvector(e%pol(k))
       enddo

     end subroutine zero_solnVector

      ! **********************************************************************
      ! * Creates a random perturbation in the EM soln - used for testing
      subroutine random_solnVector(e,eps)

        implicit none
        type (solnVector_t), intent(inout)               :: e
        real (kind=prec), intent(in), optional           :: eps
        ! local
        integer     :: j,k

        if (.not. e%allocated) then
          call errStop('EM solution not allocated in random_solnVector')
        elseif (present(eps)) then
          do k = 1,e%nPol
            call random_cvector(e%pol(k),eps)
          end do
        else
          do k = 1,e%nPol
            call random_cvector(e%pol(k),0.05*ONE)
          end do
        end if

      end subroutine random_solnVector

     !**********************************************************************
     function dotProd_solnVector(FV1,FV2,Conj_Case) result(c)
       ! computes a dot product between solution vectors

       type (solnVector_t), intent(in)         :: FV1, FV2  ! full vectors
       logical, intent(in)                     :: conj_Case ! = .true.
       complex(kind=prec)       :: c

       !  local variables
       complex(kind=prec)   :: temp
       integer              :: k

       if((.not. FV1%allocated) .or. (.not. FV2%allocated)) then
            call errStop('solution vectors have to be allocated in dotProd_solnVector')
       elseif(FV1%tx .ne. FV2%tx) then
            call errStop('different transmitters on input to dotProd_solnVector')
       endif

       c = C_ZERO

       do k = 1,FV1%nPol
           if(Conj_Case) then
           temp = dotProd_cvector_f(FV1%pol(k),FV2%pol(k))
           else
           temp = dotProd_noConj_cvector_f(FV1%pol(k),FV2%pol(k))
          endif
          c = c + temp
       enddo

     end function dotProd_solnVector

!**********************************************************************
!           Basic solnVectorMTX methods
!**********************************************************************

   subroutine create_solnVectorMTX(nTx,eAll)

      integer, intent(in)               :: nTx
      type(solnVectorMTX_t), intent(inout)  :: eAll

      !  local variables
      integer                           :: istat

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

   !************************************************************
   subroutine copy_solnVectorMTX(eOut,eIn)

       !  3D  version
       implicit none
       type (solnVectorMTX_t), intent(in)	:: eIn
       type (solnVectorMTX_t), intent(inout)	:: eOut

       ! local variables
       integer				:: j

       if (.not. eIn%allocated) then
         call errStop('input multi-transmitter EM soln not allocated yet in copy_solnVectorMTX')
       endif

       call create_solnVectorMTX(eIn%nTx,eOut)

       do j = 1,eIn%nTx
         if (eIn%solns(j)%allocated) then
            eOut%solns(j) = eIn%solns(j)
         else
            eOut%solns(j)%allocated = .false.
         end if
       end do

       eOut%allocated = eIn%allocated

       !if (eIn%temporary) then
       !   call deall_solnVectorMTX(eIn)
       !endif

   end subroutine copy_solnVectorMTX

!**********************************************************************
!           Basic sparseVector methods
!**********************************************************************
!    don't really use linear combinations ... so not implemented

    subroutine create_sparseVector(grid,iTx,LC,nCoeff)

      !  generic routine for creating the sparseVector type for 3D problems:
      !  number of polarization obtained from the transmitter dictionary

       implicit none
       type(grid_t), intent(in), target	    :: grid
       integer, intent(in)                  :: iTx
       type (sparseVector_t), intent(inout)		:: LC
       integer, intent(in), optional		:: nCoeff

       ! local variables
       integer				:: nc,k,istat

       if (present(nCoeff)) then
          nc = nCoeff
       else
          nc = 0 ! will reallocate memory later in the program
       end if

       if (LC%allocated) then
          ! it is safest and quite efficient to reallocate sparse vectors
          call deall_sparseVector(LC)
       end if

       LC%nPol = txDict(iTx)%nPol
       allocate(LC%L(LC%nPol), STAT=istat)
       do k = 1,LC%nPol
          call create_sparsevecc(nc,LC%L(k),EDGE)
       enddo

       LC%nCoeff = nc
       LC%grid => grid
       LC%tx = iTx
	   LC%allocated = .true.

    end subroutine create_sparseVector

    !************************************************************
    subroutine deall_sparseVector(LC)

      ! 3D version
      type (sparseVector_t), intent(inout)   	:: LC

      ! local variables
      integer				:: k,istat

      if (associated(LC%L)) then
         do k = 1,LC%nPol
            call deall_sparsevecc(LC%L(k))
         enddo
         deallocate(LC%L, STAT=istat)
         nullify(LC%grid)
      endif

   end subroutine deall_sparseVector

   !************************************************************
   subroutine copy_sparseVector(eOut,eIn)

       !  3D  version
       implicit none
       type (sparseVector_t), intent(in)	:: eIn
       type (sparseVector_t), intent(inout)	:: eOut

       ! local variables
       integer				:: k

       if (.not. eIn%allocated) then
         call errStop('input EM sparse vector not allocated yet in copy_sparseVector')
       endif

       call create_sparseVector(eIn%grid,eIn%tx,eOut,eIn%nCoeff)

       do k = 1,eIn%nPol
          call copy_sparsevecc(eOut%L(k),eIn%L(k))
       enddo

       eOut%allocated = eIn%allocated

       !if (eIn%temporary) then
       !   call deall_sparseVector(eIn)
       !endif

   end subroutine copy_sparseVector

!**********************************************************************
!           combined solnVector/sparseVector methods
!**********************************************************************
!  subroutine add_sparseVsolnV(cs,SV,FV)  does not appear to be needed

   function dotProd_sparseVsolnV(SV,FV,Conj_Case) result(c)

       type (sparseVector_t), intent(in)             :: SV  ! sparse vector
       type (solnVector_t), intent(in)               :: FV  ! full vector
       logical, intent(in)                     :: conj_Case ! = .true.
       complex(kind=prec)		:: c
       integer					:: k

       c = C_ZERO
       if(conj_case) then
          do k = 1,SV%nPol
             c = c + dotProd_scvector_f(SV%L(k),FV%pol(k))
          enddo
       else
          do k = 1,SV%nPol
             c = c + dotProd_noConj_scvector_f(SV%L(k),FV%pol(k))
          enddo
       endif

   end function dotProd_sparseVsolnV


!**********************************************************************
!           RHS methods
!**********************************************************************
     !  allocates and initializes arrays for the "rhs" structure
     !   set pointer to grid before calling
     !   NO gridType needed for 3DMT

     subroutine create_RHS(grid,iTx,b)
     !   3D version  ...
     !     does not create sparse vectors if sparsesource = .true.
     !       (this would in any event require knowing number of
     !         non-zero coefficients to allow for)
		 ! NOTE: Does not do anything if b%allocated is .true.

       type(grid_t), intent(in),target	:: grid
       integer, intent(in)              :: iTx
       type (RHS_t), intent(inout)   	:: b

			 if (b%allocated) then
					! do nothing - exit the create subroutine
					return
			 endif

       if (b%nonzero_BC) then
          ! create boundary condition data structures for each polarization
          !  NOTE: get rid of "used for" in BC type def
          call create_cboundary(grid,b%bc)
          b%allocated = .true.
       endif

       if (b%nonzero_source) then
          if(b%sparse_source) then
            !   can't create sparse vector without knowing how
            !    many components it has; nothing to do until we
            !    actually use (e.g., add another sparse vector)
          else
              ! Initialize full array for storage of source
              call create_cvector(grid, b%s, EDGE)
              b%allocated = .true.
          endif
       endif

       b%tx = iTx
       b%grid => grid

     end subroutine create_RHS

    !************************************************************
     subroutine deall_RHS(b)

       type (RHS_t), intent(inout)   :: b

       if (b%nonzero_BC) then
          call deall_cboundary(b%bc)

       endif

       if (b%nonzero_source) then
          if(b%sparse_source) then
             call deall_sparsevecc(b%sSparse)
          else
             call deall_cvector(b%s)
          endif
       endif

       nullify(b%grid)
       b%allocated = .false.

     end subroutine deall_RHS

     !**********************************************************************


     !**********************************************************************
     subroutine zero_RHS(b)
     !  zeros a RHS object

       type(RHS_t), intent(inout)	:: b

       if(b%nonzero_source .and. b%allocated) then
          if(b%sparse_source) then
             !  if sparse vector is zeroed, all components are
             !    deleted ...
             call deall_sparsevecc(b%sSparse)
          else
             call zero_cvector(b%s)
          endif
       else if(b%nonzero_bc .and. b%allocated) then
          call zero_cboundary(b%bc)
       else
          if(.not.b%allocated) then
             call errStop('Input not yet allocated in zero_RHS')
          endif
       endif
     end subroutine zero_RHS


!**********************************************************************
!           rhsVector methods
!**********************************************************************

     subroutine create_rhsVector(grid,iTx,b)
     !  allocates and initializes arrays for the "rhs" structure
     !   set pointer to grid before calling
     !   3D version  ...
     !     does not create sparse vectors if sparsesource = .true.
     !       (this would in any event require knowing number of
     !         non-zero coefficients to allow for)

       type(grid_t), intent(in),target	:: grid
       integer, intent(in)              :: iTx
       type (rhsVector_t), intent(inout)   	:: b

       integer				:: k,istat

       if (b%allocated) then
          if (associated(b%grid, target=grid) .and. (b%tx == iTx)) then
             ! do nothing
             return
          else
             call deall_rhsVector(b)
          end if
       end if

       b%nPol = txDict(iTx)%nPol
       allocate(b%b(b%nPol), STAT=istat)
       do k = 1,b%nPol
          b%b(k)%nonzero_bc = b%nonzero_bc
          b%b(k)%nonzero_source = b%nonzero_source
          b%b(k)%sparse_source = b%sparse_source
          call create_RHS(grid,iTx,b%b(k))
       enddo

       b%tx = iTx
       b%allocated = .true.

     end subroutine create_rhsVector

    !************************************************************
     subroutine deall_rhsVector(b)

       type (rhsVector_t), intent(inout)   :: b

       integer			:: k,istat

       if (associated(b%b)) then
          do k = 1,b%nPol
             call deall_RHS(b%b(k))
          enddo
          deallocate(b%b, STAT=istat)
       endif

       b%allocated = .false.

     end subroutine deall_rhsVector

!     !************************************************************
!     ! need copy_RHS for this... too complicated, will write this
!     ! when it is needed
!     subroutine copy_rhsVector(bOut,bIn)
!
!       !  3D  version
!       implicit none
!       type (rhsVector_t), intent(in)	:: bIn
!       type (rhsVector_t), intent(inout)	:: bOut
!
!       ! local variables
!       integer				:: k
!
!       if (.not. bIn%allocated) then
!         call errStop('input EM RHS not allocated yet in copy_rhsVector')
!       endif
!
!       call create_rhsVector(bIn%grid,bIn%tx,bOut)
!
!       do k = 1,bIn%nPol
!          call copy_RHS(bOut%b(k),bIn%b(k))
!       enddo
!
!       bOut%allocated = bIn%allocated
!
!       if (bIn%temporary) then
!          call deall_rhsVector(bIn)
!       endif
!
!     end subroutine copy_rhsVector

     !**********************************************************************
     subroutine zero_rhsVector(b)
     !  zeros a rhsVector object

       type(rhsVector_t), intent(inout)	:: b

       integer			:: k

       do k = 1,b%nPol
          call zero_RHS(b%b(k))
       enddo
     end subroutine zero_rhsVector

  ! **********************************************************************
  ! * Creates a random perturbation in the EM RHS - used for testing
  subroutine random_rhsVector(b,eps)

    implicit none
    type (rhsVector_t), intent(inout)                :: b
    real (kind=prec), intent(in), optional           :: eps
    ! local
    integer     :: k

       do k = 1,b%nPol
        if (.not. (b%b(k)%allocated .and. b%b(k)%nonzero_source)) then
          call errStop('EM RHS not allocated in random_rhsVector')
        elseif (present(eps)) then
          call random_cvector(b%b(k)%s,eps)
        else
          call random_cvector(b%b(k)%s,0.05*ONE)
        end if
       enddo

  end subroutine random_rhsVector

     !**********************************************************************
     subroutine copy_solnVrhsV(b,e)
     !  implements b = e

       type(rhsVector_t), intent(inout) :: b
       type(solnVector_t), intent(in)   :: e

       integer          :: k

       if (.not. e%allocated) then
         call errStop('input EM soln not allocated yet in copy_solnVrhsV')
       endif

       call create_rhsVector(e%grid,e%tx,b)

       do k = 1,b%nPol
          b%b(k)%nonzero_source = .true.
          b%b(k)%s = e%pol(k)
          b%b(k)%sparse_source = .false.
          b%b(k)%nonzero_BC = .false.
          b%b(k)%allocated = .true.
       enddo
     end subroutine copy_solnVrhsV


     !**********************************************************************

     subroutine add_sparseVrhsV(cs,SV,comb)

     !   Forms the linear combination cs*SV+comb where:
     !     cs is complex
     !     SV is an sparseVector object
     !     comb is an rhsVector object
     !   Result is returned in comb
     !
     !   In this implementation, an rhsVector object

       complex(kind=prec), intent(in)  :: cs
       type (sparseVector_t), intent(in)             :: SV  ! sparse vector
       type (rhsVector_t), intent(inout)             :: comb  ! full vector

     !  local variables
       type(sparsevecc)			:: temp
       integer				:: k

       do k = 1,comb%nPol
          comb%b(k)%nonzero_source = .true.
          if(comb%b(k)%sparse_Source) then
             ! use sparse vector storage for output
             if(comb%b(k)%sSparse%allocated) then
                !  add to contentes of comb sparse vector and cs*SV
                call linComb_sparsevecc(comb%b(k)%sSparse,C_ONE, &
			SV%L(k),cs,temp)
             else
                call scMult_sparsevecc(cs,SV%L(k),temp)
             endif
             call copy_sparsevecc(temp,comb%b(k)%sSparse)
          else
             !  might want to check for allocation of comb%b(k)%s first
             call add_scvector(cs,SV%L(k),comb%b(k)%s)
          endif
       enddo

       call deall_sparsevecc(temp)

     end subroutine add_sparseVrhsV

!**********************************************************************
!           combined solnVector/rhsVector methods
!**********************************************************************
     function dotProd_rhsVsolnV(comb,FV,Conj_Case) result(c)
       ! computes a dot product between RHS vector and solution vector;
       ! does not compute anything from the boundary conditions!!!
       ! this might have to be changed.

       type (rhsVector_t), intent(in)             :: comb  ! rhs
       type (solnVector_t), intent(in)            :: FV  ! full vector
       logical, intent(in)                     :: conj_Case ! = .true.
       complex(kind=prec)       :: c

       !  local variables
       complex(kind=prec)   :: temp
       integer              :: k

       if((.not. comb%allocated) .or. (.not. FV%allocated)) then
            call errStop('RHS and solution vectors have to be allocated in dotProd_rhsVsolnV')
       elseif(comb%tx .ne. FV%tx) then
            call errStop('different transmitters on input to dotProd_rhsVsolnV')
       endif

       c = C_ZERO

       do k = 1,comb%nPol
          if(comb%b(k)%nonzero_source) then
             if(comb%b(k)%sparse_source) then
                if(Conj_Case) then
                temp = dotProd_scvector_f(comb%b(k)%sSparse,FV%pol(k))
                else
                temp = dotProd_noConj_scvector_f(comb%b(k)%sSparse,FV%pol(k))
                endif
             else
                if(Conj_Case) then
                temp = dotProd_cvector_f(comb%b(k)%s,FV%pol(k))
                else
                temp = dotProd_noConj_cvector_f(comb%b(k)%s,FV%pol(k))
                endif
             endif
          else
             temp = C_ZERO
          endif
          c = c + temp
       enddo

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

    type(rhsVectorMTX_t), intent(inout)     :: bAll
    real (kind=prec), intent(in), optional  :: eps

    !  local variables
    integer                           :: j

    do j = 1,bAll%nTx
      if (present(eps)) then
        call random_rhsVector(bAll%combs(j),eps)
      else
        call random_rhsVector(bAll%combs(j),0.05*ONE)
      end if
    end do

   end subroutine random_rhsVectorMTX

end module SolnSpace
