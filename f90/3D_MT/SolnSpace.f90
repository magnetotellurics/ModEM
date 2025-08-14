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
   MODULE PROCEDURE copy_rhsVector
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
	!! NM: to make it possibile to have both MT and CSEM in one eAll solution, the nPol must be general.
	!! Since nPol is defined from the Tx dictionary, there is NO need to have an abstract nPol.
	!! AK: added Pol_name for more flexible read/write. 2018-04-27
    integer		          :: nPol  
    integer ,pointer      :: Pol_index(:) 	
    character(20), pointer  :: Pol_name(:)
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
    !! Since nPol is defined from the Tx dictionary, there is NO need to have an abstract nPol.
    integer						:: nPol  
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
     ! this logic is not obvious, BUT: this can contain both source
     ! and BC; however, the source either has a sparse or full
     ! representation, can't have both. So, if nonzero_Source,
     ! check sparse_Source, otherwise, source is zero.

     character*3		:: adj = ''
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
     ! this logic is not obvious, BUT: this can contain both source
     ! and BC; however, the source either has a sparse or full
     ! representation, can't have both. So, if nonzero_Source,
     ! check sparse_Source, otherwise, source is zero.

     integer				:: nPol
     character(20), pointer	:: Pol_name(:) ! Ex or Ey
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
	   allocate(e%Pol_index(e%nPol), STAT=istat)
       allocate(e%Pol_name(e%nPol), STAT=istat)
       
       do iPol=1,e%nPol
        e%Pol_index(iPol)=iPol
       end do
       
       ! set up the mode names based on transmitter type;
       ! for now, only set up for MT. Used for I/O.
       if (trim(txDict(iTx)%tx_type) .eq. 'MT') then
        if (e%nPol == 2) then
            e%Pol_name(1) = 'Ex'
            e%Pol_name(2) = 'Ey'
        else
         call errStop('problem creating MT modes in create_solnVector')
        end if
       end if

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
          deallocate(e%Pol_index, STAT=istat)
          deallocate(e%Pol_name, STAT=istat)
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
     function construct_esoln_fname(prefix, iTx, polarization) result(fname)

       implicit none

       character(len=*), intent(in) :: prefix
       integer, intent(in) :: iTx
       integer, intent(in) :: polarization
       character(len=512)  :: fname

       write(fname, '(A, A, I4.4, A, I2.2, A)') trim(prefix), '.iTx.', iTX, 'pol', polarization, '.cvec'

     end function construct_esoln_fname

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

   ! **********************************************************************
   ! * Creates a random perturbation in the EM soln - used for testing
   subroutine random_solnVectorMTX(eAll,eps)

        implicit none
        type (solnVectorMTX_t), intent(inout)            :: eAll
        real (kind=prec), intent(in), optional           :: eps
        ! local
        integer     :: j,k

        if (.not. eAll%allocated) then
          call errStop('EM solution not allocated in random_solnVectorMTX')
        elseif (present(eps)) then
          do j = 1,eAll%nTx
            call random_solnVector(eAll%solns(j),eps)
          end do
        else
          do j = 1,eAll%nTx
            call random_solnVector(eAll%solns(j),0.05*ONE)
          end do
        end if

   end subroutine random_solnVectorMTX

   !**********************************************************************
   subroutine getBC_solnVectorMTX(eAll,bAll)

      type(solnVectorMTX_t), intent(in)    :: eAll
      type(rhsVectorMTX_t), intent(inout)  :: bAll

      !  local variables
      integer                           :: j,istat

      if (.not. eAll%allocated) then
          call errStop('EM solution not allocated in getBC_solnVectorMTX')
      end if

      call create_rhsVectorMTX(eAll%nTx,bAll)
      do j = 1,eAll%nTx
          call getBC_solnVector(eAll%solns(j),bAll%combs(j))
      end do

   end subroutine getBC_solnVectorMTX

   !**********************************************************************
   subroutine count_solnVectorMTX(eAll,nTx_nPol)

      type(solnVectorMTX_t), intent(in)    :: eAll
      integer, intent(inout)               :: nTx_nPol

      !  local variables
      integer                           :: j,istat

      if (.not. eAll%allocated) then
          call errStop('EM solution not allocated in count_solnVectorMTX')
      end if

      nTx_nPol = 0
      do j = 1,eAll%nTx
          nTx_nPol = nTx_nPol + eAll%solns(j)%nPol
      end do

   end subroutine count_solnVectorMTX

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

   ! **********************************************************************
   ! * Creates a random perturbation in the EM sparse soln - used for testing
   subroutine random_sparseVector(e,eps)

        implicit none
        type (sparseVector_t), intent(inout)             :: e
        real (kind=prec), intent(in), optional           :: eps
        ! local
        integer     :: j,k

        if (.not. e%allocated) then
          call errStop('EM sparse vector not allocated in random_sparseVector')
        elseif (present(eps)) then
          do k = 1,e%nPol
            call random_sparsevecc(e%L(k),eps)
          end do
        else
          do k = 1,e%nPol
            call random_sparsevecc(e%L(k),0.05*ONE)
          end do
        end if

   end subroutine random_sparseVector

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
					!return        !NM: In case of MT and CSEM we need to create b for b%nonzero_source
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

     !************************************************************
     subroutine copy_RHS(bOut,bIn)

       !  3D  version
       implicit none
       type (RHS_t), intent(in)   :: bIn
       type (RHS_t), intent(inout)    :: bOut

       ! local variables
       integer              :: k

       if (.not. bIn%allocated) then
         call errStop('input EM RHS not allocated yet in copy_RHS')
       endif

       bOut%nonzero_Source = bIn%nonzero_Source
       bOut%sparse_Source = bIn%sparse_Source
       bOut%nonzero_BC = bIn%nonzero_BC
       call create_RHS(bIn%grid,bIn%tx,bOut)

       if (bIn%nonzero_BC) then
          call copy_cboundary(bOut%bc,bIn%bc)
       endif

       if (bIn%nonzero_Source) then
          if(bIn%sparse_Source) then
            call copy_sparsevecc(bOut%sSparse,bIn%sSparse)
          else
            call copy_cvector(bOut%s,bIn%s)
          endif
       endif

       bOut%adj = bIn%adj
       bOut%allocated = bIn%allocated

       !if (bIn%temporary) then
       !   call deall_RHS(bIn)
       !endif

     end subroutine copy_RHS

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
       allocate(b%Pol_name(b%nPol), STAT=istat)

       ! set up the mode names based on transmitter type;
       ! for now, only set up for MT. Used for I/O.
       if (trim(txDict(iTx)%tx_type) .eq. 'MT') then
        if (b%nPol == 2) then
            b%Pol_name(1) = 'Ex'
            b%Pol_name(2) = 'Ey'
        else
         call errStop('problem creating MT modes in create_rhsVector')
        end if
       end if

       do k = 1,b%nPol
          b%b(k)%nonzero_bc = b%nonzero_bc
          b%b(k)%nonzero_source = b%nonzero_source
          b%b(k)%sparse_source = b%sparse_source
          call create_RHS(grid,iTx,b%b(k))
       enddo

       b%tx = iTx
       b%grid => grid

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
          deallocate(b%Pol_name, STAT=istat)
          deallocate(b%b, STAT=istat)
       endif

       b%allocated = .false.

     end subroutine deall_rhsVector

     !************************************************************
     ! need copy_RHS for this... done by AK 2018-04-25
     subroutine copy_rhsVector(bOut,bIn)

       !  3D  version
       implicit none
       type (rhsVector_t), intent(in)	:: bIn
       type (rhsVector_t), intent(inout)	:: bOut

       ! local variables
       integer				:: k

       if (.not. bIn%allocated) then
         call errStop('input EM RHS not allocated yet in copy_rhsVector')
       endif

       bOut%nonzero_Source = bIn%nonzero_Source
       bOut%sparse_Source = bIn%sparse_Source
       bOut%nonzero_BC = bIn%nonzero_BC
       call create_rhsVector(bIn%grid,bIn%tx,bOut)

       do k = 1,bIn%nPol
          call copy_RHS(bOut%b(k),bIn%b(k))
       enddo

       bOut%allocated = bIn%allocated

       !if (bIn%temporary) then
       !   call deall_rhsVector(bIn)
       !endif

     end subroutine copy_rhsVector

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
     ! * Creates a random perturbation in the EM RHS; AK 2018-04-25
     subroutine random_rhsVector(b,eps)

         implicit none
         type (rhsVector_t), intent(inout)                :: b
         real (kind=prec), intent(in), optional           :: eps
         ! local
         integer     :: k

         ! allow for multiple types of information if needed
         do k = 1,b%nPol
             if (.not. b%b(k)%allocated) then
                 call errStop('EM RHS not allocated in random_rhsVector')
             else

                 if (b%b(k)%nonzero_Source) then
                     if (b%b(k)%sparse_Source) then
                         ! sparse vector nonzero source
                         if (present(eps)) then
                             call random_sparsevecc(b%b(k)%sSparse,eps)
                         else
                             call random_sparsevecc(b%b(k)%sSparse,0.05*ONE)
                         end if
                     else
                         ! full vector nonzero source
                         if (present(eps)) then
                             call random_cvector(b%b(k)%s,eps)
                         else
                             call random_cvector(b%b(k)%s,0.05*ONE)
                         end if
                     end if
                 end if

                 ! nonzero boundary conditions
                 if (b%b(k)%nonzero_BC) then
                     if (present(eps)) then
                         call random_cboundary(b%b(k)%bc,eps)
                     else
                         call random_cboundary(b%b(k)%bc,0.05*ONE)
                     end if
                 end if
             end if
         enddo

     end subroutine random_rhsVector

     !**********************************************************************
     subroutine copy_solnVrhsV(b,e)
     !  implements b = e; copies to a full source vector

       type(rhsVector_t), intent(inout)     :: b
       type(solnVector_t), intent(in)       :: e

       integer          :: k

       if (.not. e%allocated) then
         call errStop('input EM soln not allocated yet in copy_solnVrhsV')
       endif

       b%nonzero_Source = .true.
       b%sparse_Source = .false.
       b%nonzero_BC = .false.
       call create_rhsVector(e%grid,e%tx,b)

       do k = 1,b%nPol
            b%b(k)%s = e%pol(k)
            b%b(k)%allocated = .true.
       enddo

     end subroutine copy_solnVrhsV

     !**********************************************************************
     subroutine getBC_solnVector(e,b)
     !  extracts the BC from e to save in b

       type(rhsVector_t), intent(inout)     :: b
       type(solnVector_t), intent(in)       :: e

       integer          :: k

       if (.not. e%allocated) then
         call errStop('input EM soln not allocated yet in getBC_solnVrhsV')
       endif

       b%nonzero_Source = .false.
       b%sparse_Source = .false.
       b%nonzero_BC = .true.
       call create_rhsVector(e%grid,e%tx,b)

       do k = 1,b%nPol
            call getBC(e%pol(k),b%b(k)%bc)
            b%b(k)%allocated = .true.
       enddo

     end subroutine getBC_solnVector

     !**********************************************************************

     subroutine sparse2full_rhsVector(b)

     !   Converts an rhsVector object from sparse to full representation,
     !   by consolidating sSparse and bc into a full cvector s

       type (rhsVector_t), intent(inout)             :: b

     !  local variables
       type(sparsevecc)         :: temp
       integer              :: k

       do k = 1,b%nPol
          if (.not. b%b(k)%s%allocated) then
            call create_cvector(b%b(k)%grid,b%b(k)%s,EDGE)
          end if
          if(b%b(k)%nonzero_Source .and. b%b(k)%sparse_Source) then
             ! add sparse vector to full
             call add_scvector(C_ONE,b%b(k)%sSparse,b%b(k)%s)
          end if
          if(b%b(k)%nonzero_BC) then
             ! convert BC to sparse vector and add to full
             call copy_bsparsevecc(b%b(k)%bc,temp)
             call add_scvector(C_ONE,temp,b%b(k)%s)
             call deall_sparsevecc(temp)
          end if
          b%b(k)%nonzero_Source = .true.
          b%b(k)%sparse_Source = .false.
          b%b(k)%nonzero_BC = .false.
       enddo

       b%nonzero_Source = .true.
       b%sparse_Source = .false.
       b%nonzero_BC = .false.


     end subroutine sparse2full_rhsVector

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

       comb%nonzero_source = .true.

       do k = 1,comb%nPol
          comb%b(k)%nonzero_source = comb%nonzero_source
          comb%b(k)%sparse_source = comb%sparse_source
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

   !**********************************************************************
   subroutine count_rhsVectorMTX(bAll,nTx_nPol)

      type(rhsVectorMTX_t), intent(in)    :: bAll
      integer, intent(inout)               :: nTx_nPol

      !  local variables
      integer                           :: j,istat

      if (.not. bAll%allocated) then
          call errStop('RHS not allocated in count_rhsVectorMTX')
      end if

      nTx_nPol = 0
      do j = 1,bAll%nTx
          nTx_nPol = nTx_nPol + bAll%combs(j)%nPol
      end do

   end subroutine count_rhsVectorMTX

end module SolnSpace
