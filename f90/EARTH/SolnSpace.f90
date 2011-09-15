module SolnSpace
!   higher level module to define EM solution and RHS objects
!   plus basic methods, linear algebra, dot products
!   EARTH version
!
! Defines: solnVector, sparseVector, rhsVector
! Uses: EMfield

use math_constants
use file_units
use utilities
use sg_scalar
use sg_vector
use sg_boundary
use sg_sparse_vector
use transmitters
use receivers

implicit none

interface assignment (=)
   !MODULE PROCEDURE copy_rhsVector - doesn't exist yet
   MODULE PROCEDURE copy_solnVrhsV
   MODULE PROCEDURE copy_solnVector
   MODULE PROCEDURE copy_solnVectorMTX
   !MODULE PROCEDURE copy_sparseVector - doesn't exist yet
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
    type(cvector)           :: vec
    type(grid_t), pointer   :: grid
    integer                 :: tx = 0
    integer                 :: errflag = 0
    logical                 :: allocated = .false.
    logical                 :: temporary = .false. ! true for function outputs only
  end type solnVector_t

  type :: solnVectorMTX_t
    !! Generic solution type for storing solutions from multiple transmitters
    integer         :: nTx = 0
    type(solnVector_t), pointer     :: solns(:)
    logical         :: allocated = .false.
    logical         :: temporary = .false. ! true for function outputs only
  end type solnVectorMTX_t


  type :: sparseVector_t
    !!   Generic solution type, same name must be used to allow
    !!   use of higher level inversion modules on different problems.
    type(sparsevecc)            :: L
    logical                     :: temporary = .false. ! true for function outputs only
  end type sparseVector_t


  type :: rhsVector_t
     !!   right hand side for solving Maxwell's equations
     character*3            :: adj = ''
     logical                :: nonzero_source = .false.
     logical                :: nonzero_bc = .false.
     logical                :: allocated = .false.
    logical                 :: temporary = .false. ! true for function outputs only
     type(cvector)          :: source
     complex(kind=prec), pointer, dimension(:)  :: bc
     type(grid_t), pointer      :: grid
     integer                    :: tx = 0
  end type rhsVector_t


  type :: rhsVectorMTX_t
    !! Generic solution type for storing RHS's for multiple transmitters
    integer         :: nTx = 0
    type(rhsVector_t), pointer     :: combs(:)
    logical         :: allocated = .false.
    logical         :: temporary = .false. ! true for function outputs only
  end type rhsVectorMTX_t


  ! ***************************************************************************
  ! * full solution information for a single 'slice' at a fixed radius (in km)
  type :: solution_t

    character(80)                               :: type  ! CELL/NODE
    type (receiver_t), pointer, dimension(:,:)  :: o
    complex(8), pointer, dimension(:,:)         :: x,y,z    !nx,ny

  end type solution_t

Contains

!**********************************************************************
!           Basic solnVectorMTX methods
!**********************************************************************

   subroutine create_solnVectorMTX(nTx,eAll)

      integer, intent(in)               :: nTx
      type(solnVectorMTX_t), intent(inout)  :: eAll

      ! local variables
      integer      :: istat

      call deall_solnVectorMTX(eAll)

      eAll%nTx = nTx
      allocate(eAll%solns(nTx), STAT=istat)
      eAll%allocated = .true.

   end subroutine create_solnVectorMTX

   !**********************************************************************
   subroutine deall_solnVectorMTX(eAll)

      type(solnVectorMTX_t)     :: eAll

      !  local variables
      integer                   :: j, istat

      if (eAll%allocated) then
        do j = 1,eAll%nTx
          call deall_solnVector(eAll%solns(j))
        end do
      end if

      if (associated(eAll%solns)) deallocate(eAll%solns, STAT=istat)
      eAll%nTx = 0
      eAll%allocated = .false.

   end subroutine deall_solnVectorMTX

  ! **********************************************************************
  ! * Creates a random perturbation in the EM soln - used for testing
  subroutine random_solnVectorMTX(eAll,eps)

    implicit none
    type (solnVectorMTX_t), intent(inout)            :: eAll
    real (kind=prec), intent(in), optional           :: eps
    ! local
    integer     :: j

    if (.not. eAll%allocated) then
      call errStop('EM solution not allocated in random_solnVectorMTX')
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

  !****************************************************************************
  ! linComb_solnVectorMTX computes linear combination of two field solutions
  ! stored for multiple transmitters; subroutine, not a function
  ! both input vectors must have the same dimension
  subroutine linComb_solnVectorMTX(c1, E1, c2, E2, E3)

    implicit none
    !   input vectors
    type (solnVectorMTX_t), intent(in)     :: E1, E2
    !  input complex scalars
    complex (kind=prec), intent(in)        :: c1, c2
    type (solnVectorMTX_t), intent(inout)  :: E3
    ! local
    integer                                :: j

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
    	call errStop('inputs not allocated yet for linComb_solnVectorMTX')
    elseif (E1%nTx .ne. E2%nTx) then
    	call errStop('inputs of different sizes for linComb_solnVectorMTX')
    end if

    ! check to see if LHS (E3) is active (allocated)
    if((.not.E3%allocated) .or. (E3%nTx .ne. E1%nTx)) then
    	call create_solnVectorMTX(E1%nTx,E3)
    end if

    do j=1,E1%nTx
    	if (E1%solns(j)%tx == E2%solns(j)%tx) then
    		! initialize...
        	E3%solns(j) = E1%solns(j)
        	! form linear combination
			call linComb_solnVector(c1, E1%solns(j), c2, E2%solns(j), E3%solns(j))
		else
			call errStop('inputs for different transmitters in linComb_solnVectorMTX')
		end if
    end do

  end subroutine linComb_solnVectorMTX ! linComb_solnVectorMTX

  !****************************************************************************
  ! write_solnVectorMTX writes an ASCII data file containing the full solnVectorMTX
  subroutine write_solnVectorMTX(fname, E)

    implicit none
    !   input vectors
    character(*), intent(in)			   :: fname
    type (solnVectorMTX_t), intent(in)         :: E
    ! local
    integer                                :: j,ios

    if(.not.E%allocated) then
    	call errStop('input not allocated yet for write_solnVectorMTX')
    end if

	open(ioWRITE,file=fname,status='unknown',form='formatted',iostat=ios)
	write(ioWRITE,'(a36,i3,a8)') "# Full EM field solution output for ",E%nTx,"periods."
	do j = 1,E%nTx
		write(ioWRITE,'(i3)') E%solns(j)%tx
		call write_cvector(ioWRITE,E%solns(j)%vec)
	end do
	close(ioWRITE)

  end subroutine write_solnVectorMTX ! write_solnVectorMTX

  !****************************************************************************
  ! read_solnVectorMTX reads an ASCII data file containing the full solnVectorMTX
  subroutine read_solnVectorMTX(fname, E, grid)

    implicit none
    !   input vectors
    character(*), intent(in)			   :: fname
    type (solnVectorMTX_t), intent(inout)  :: E
    type (grid_t), intent(in), target	   :: grid
    ! local
    integer                                :: j,iTx,nTx,ios,istat
    character(100)						   :: comment

    if(E%allocated) then
    	call deall_solnVectorMTX(E)
    end if

	open(ioREAD,file=fname,status='unknown',form='formatted',iostat=ios)
	read(ioREAD,'(a35)',iostat=istat,advance='no') comment
	read(ioREAD,*,iostat=istat) nTx
	print *, node_info,'Number of transmitters: ',nTx
    call create_solnVectorMTX(nTx,E)
	do j = 1,nTx
		read(ioREAD,'(i3)',iostat=istat) iTx
        call create_solnVector(grid,iTx,E%solns(j))
		call read_cvector(ioREAD,E%solns(j)%vec,grid)
		E%solns(j)%errflag = 0
		E%solns(j)%grid => grid
	end do
	close(ioREAD)

  end subroutine read_solnVectorMTX ! read_solnVectorMTX

!**********************************************************************
!           Basic solnVector methods
!**********************************************************************
     subroutine create_solnVector(grid,iTx,e)
       ! the interface is generic, and could be used as such
       ! ... note gridType = EDGE for all solution vectors

       implicit none
       type(grid_t), intent(in), target     :: grid
       integer, intent(in)                  :: iTx
       type (solnVector_t), intent(inout)   :: e
       ! local
       character(80)  :: gridType

       call create_cvector(grid,e%vec,EDGE)
       e%tx = iTx
       e%errflag = 0
       e%grid => grid
       e%allocated = .true.

     end subroutine create_solnVector

     !************************************************************
     subroutine deall_solnVector(e)

       type (solnVector_t)   :: e

       call deall_cvector(e%vec)
       if(associated(e%grid)) then
           nullify(e%grid)
       endif
       e%allocated = .false.

     end subroutine deall_solnVector

  !****************************************************************************
  ! linComb_solnVector computes linear combination of two field solutions
  ! stored for multiple transmitters; subroutine, not a function
  ! both input vectors must have the same dimension
  subroutine linComb_solnVector(c1, E1, c2, E2, E3)

    implicit none
    !   input vectors
    type (solnVector_t), intent(in)         :: E1, E2
    !  input complex scalars
    complex (kind=prec), intent(in)        :: c1, c2
    type (solnVector_t), intent(inout)      :: E3
    ! local
    integer                                :: j

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
        call errStop('inputs not allocated yet for linComb_solnVector')
    elseif (E1%tx .ne. E2%tx) then
        call errStop('different transmitters on input to linComb_solnVector')
    end if

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
        call create_solnVector(E1%grid,E1%tx,E3)
    end if

    ! initialize...
    E3%vec = E1%vec
    ! form linear combination
    call linComb_cvector(c1, E1%vec, c2, E2%vec, E3%vec)

    E3%tx = E1%tx
    E3%errflag = max(E1%errflag, E2%errflag)
    E3%grid => E1%grid
    E3%allocated = .true.

  end subroutine linComb_solnVector

     !************************************************************
     subroutine copy_solnVector(eOut,eIn)

       implicit none
       type (solnVector_t), intent(in)  :: eIn
       type (solnVector_t), intent(inout)   :: eOut

       !  should have some error checking for eIn ...
       call copy_cvector(eOut%vec,eIn%vec)

       eOut%tx = eIn%tx
       eOut%errflag = eIn%errflag
       eOut%grid => eIn%grid
       eOut%allocated = eIn%allocated

       if(eIn%temporary) then
            call deall_solnVector(eIn)
       endif

     end subroutine copy_solnVector

     !**********************************************************************
     function dotProd_solnVector(e1,e2,conj_Case) result(c)
     !   dot product of two solution space objects
       type(solnVector_t), intent(in)       :: e1,e2
       logical, intent(in)      :: conj_case
       complex(kind=prec)   :: c

        if(conj_case) then
            c = dotProd_cvector_f(E1%vec,E2%vec)
        else
            c = dotProd_noConj_cvector_f(E1%vec,E2%vec)
        endif

     end function dotProd_solnVector

     !**********************************************************************
     subroutine zero_solnVector(e)
     !  zeros a solution space object

       type(solnVector_t), intent(inout)    :: e

       call zero_cvector(e%vec)

     end subroutine zero_solnVector

  ! **********************************************************************
  ! * Creates a random perturbation in the EM soln - used for testing
  subroutine random_solnVector(e,eps)

    implicit none
    type (solnVector_t), intent(inout)               :: e
    real (kind=prec), intent(in), optional           :: eps
    ! local
    integer     :: j

    if (.not. e%allocated) then
      call errStop('EM solution not allocated in random_solnVector')
    elseif (present(eps)) then
      call random_cvector(e%vec,eps)
    else
      call random_cvector(e%vec,0.05*ONE)
    end if

  end subroutine random_solnVector

  !****************************************************************************
  ! write_solnVector writes an ASCII data file containing one solnVector
  subroutine write_solnVector(fname, E)

    implicit none
    !   input vectors
    character(*), intent(in)               :: fname
    type (solnVector_t), intent(in)         :: E
    ! local
    character(80)                          :: code
    integer                                :: j,ios
    character(100)                         :: fn_output

    if(.not.E%allocated) then
        call errStop('input not allocated yet for write_solnVector')
    end if

    code = freqList%info(E%tx)%code !write (code,'(i3.3)') E%tx
    fn_output = trim(fname)//'_'//trim(code)//'.field'
    write(*,*) 'Writing to file: ',trim(fn_output)
    open(ioWRITE,file=fn_output,status='unknown',form='formatted',iostat=ios)
    write(ioWRITE,'(a45,f9.3,a6)') "# Full EM field solution output for period ",   &
                                        freqList%info(E%tx)%period,' days.'
    write(ioWRITE,'(i3)') E%tx
    call write_cvector(ioWRITE,E%vec)
    close(ioWRITE)

  end subroutine write_solnVector ! write_solnVector

  !****************************************************************************
  ! read_solnVector reads an ASCII data file containing one solnVector
  subroutine read_solnVector(fname, grid, iTx, E)

    implicit none
    !   input vectors
    character(*), intent(in)               :: fname
    type (solnVector_t), intent(inout)  :: E
    integer, intent(in)                 :: iTx
    type (grid_t), intent(in), target      :: grid
    ! local
    integer                                :: j,ios,istat
    character(80)                          :: code
    character(100)                         :: comment, fn_input

    if(E%allocated) then
        call deall_solnVector(E)
    end if
    call create_solnVector(grid,iTx,E)

    code = freqList%info(iTx)%code !write (code,'(i3.3)') iTx
    fn_input = trim(fname)//'_'//trim(code)//'.field'
    write(*,*) node_info,'Reading from file: ',trim(fn_input)
    open(ioREAD,file=fn_input,status='unknown',form='formatted',iostat=ios)
    read(ioREAD,'(a35)',iostat=istat) comment
    read(ioREAD,'(i3)',iostat=istat) j
    if (j .ne. iTx) then
        write(0,*) node_info,'Warning: transmitter ',iTx,' is read from file ',j
    end if
    call read_cvector(ioREAD,E%vec,grid)
    E%errflag = 0
    E%grid => grid
    close(ioREAD)

  end subroutine read_solnVector ! read_solnVector

!**********************************************************************
!           Basic sparseVector methods
!**********************************************************************

     subroutine create_sparseVector(grid,iTx,LC,nCoeff)
       ! ... note gridType = EDGE for all solution vectors

       type(grid_t), intent(in) :: grid
       integer, intent(in)      :: iTx
       type(sparseVector_t)         :: LC
       integer, optional, intent(in)    :: nCoeff

       ! local
       integer           :: nc
       character(80)     :: gridType

       if (present(nCoeff)) then
          nc = nCoeff
       else
          nc = 0 ! will reallocate memory later in the program
       end if

       call create_sparsevecc(nc,LC%L,EDGE)

     end subroutine create_sparseVector

     !************************************************************
     subroutine deall_sparseVector(LC)

       type (sparseVector_t), intent(inout)     :: LC

       call deall_sparsevecc(LC%L)

     end subroutine deall_sparseVector

     !**********************************************************************
     subroutine linComb_sparseVector(Lin1,c1,Lin2,c2,Lout)
       ! linear combination of two sparseVector objects

       type (sparseVector_t), intent(in)        :: Lin1,Lin2
       complex (kind=prec), intent(in)  :: c1,c2
       type (sparseVector_t), intent(inout)     :: Lout

       call linComb_sparsevecc(Lin1%L,c1,Lin2%L,c2,Lout%L)

     end subroutine linComb_sparseVector

!**********************************************************************
!           combined solnVector/sparseVector methods
!**********************************************************************
     function dotProd_sparseVsolnV(SV,FV,Conj_Case) result(c)

       type (sparseVector_t), intent(in)             :: SV  ! sparse vector
       type (solnVector_t), intent(in)               :: FV  ! full vector
       logical, intent(in)                     :: conj_Case ! = .true.
       complex(kind=prec)       :: c

        if(conj_case) then
            c = dotProd_scvector_f(SV%L,FV%vec)
        else
            c = dotProd_noConj_scvector_f(SV%L,FV%vec)
        endif
       !c = dotProd_scvector(SV%L,FV%vec,Conj_Case)

     end function dotProd_sparseVsolnV

     !**********************************************************************
     subroutine add_sparseVsolnV(cs,SV,FV)

        complex(kind=prec), intent(in)  :: cs
        type (sparseVector_t), intent(in)   :: SV  ! sparse vector
        type (solnVector_t), intent(inout)  :: FV  ! full vector

       call add_scvector(cs,SV%L,FV%vec)

     end subroutine add_sparseVsolnV

!**********************************************************************
!           rhsVector methods
!**********************************************************************
     !  allocates and initializes arrays for the "rhs" structure
     !   set pointer to grid + mode/source fields of rhs before calling
     subroutine create_rhsVector(grid,iTx,b)

       type(grid_t), intent(in),target  :: grid
       integer, intent(in)              :: iTx
       type (rhsVector_t), intent(inout)    :: b

       !  local variables
       character(2)     :: mode
       character(80)    :: gridType
       integer          :: Nx,Ny,nBC

       Nx = grid%Nx
       Ny = grid%Ny
       gridType = EDGE

       nBC = 2*Nx*Ny

       if (b%nonzero_BC) then
          ! Initialize array for boundary condition data
          allocate(b%BC(nBC))
       endif

       if (b%nonzero_source) then
           ! Initialize full array for storage of source
           call create_cvector(grid,b%source,gridType)
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

        complex(kind=prec), intent(in)  :: cs
        type (sparseVector_t), intent(in)       :: SV  ! sparse vector
        type (rhsVector_t), intent(inout)       :: FV  ! full vector

       call add_scvector(cs,SV%L,FV%source)
       FV%nonzero_source = .true.
       FV%allocated = .true.

     end subroutine add_sparseVrhsV

     !**********************************************************************
     subroutine zero_rhsVector(b)
     !  zeros a solution space object

       type(rhsVector_t), intent(inout) :: b

       if(b%nonzero_source .and. b%allocated) then
          call zero_cvector(b%source)
       else
          if(b%nonzero_bc .and. b%allocated) then
             b%bc = C_ZERO
          endif
       endif
     end subroutine zero_rhsVector

  ! **********************************************************************
  ! * Creates a random perturbation in the EM RHS - used for testing
  subroutine random_rhsVector(b,eps)

    implicit none
    type (rhsVector_t), intent(inout)                :: b
    real (kind=prec), intent(in), optional           :: eps
    ! local
    integer     :: j

    if (.not. (b%allocated .and. b%nonzero_source)) then
      call errStop('EM RHS not allocated in random_rhsVector')
    elseif (present(eps)) then
      call random_cvector(b%source,eps)
    else
      call random_cvector(b%source,0.05*ONE)
    end if

  end subroutine random_rhsVector

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

       if(e%temporary) then
            call deall_solnVector(e)
       endif

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
            if(conj_case) then
            c = dotProd_cvector_f(comb%source,FV%vec)
            else
            c = dotProd_noConj_cvector_f(comb%source,FV%vec)
            endif
            !c = dotProd_cvector(comb%source,FV%vec,Conj_Case)
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
