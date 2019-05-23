module SymmetryTest
!  higher level unility module that tests the implementation of the
!  Jacobian computations for each of the critical operators in
!  J = L S^{-1} P + Q.
!
!   write(0,*) 'Usage: Test the adjoint implementation for each of the critical'
!   write(0,*) '       operators in J = L S^{-1} P + Q'
!   write(0,*)
!   write(0,*) ' -A  J rFile_Model rFile_Data [wFile_Model wFile_Data]'
!   write(0,*) '  Tests the equality d^T J m = m^T J^T d for any model and data.'
!   write(0,*) '  Optionally, outputs J m and J^T d.'
!   write(0,*)
!   write(0,*) ' -A  L rFile_EMsoln rFile_Data [wFile_EMrhs wFile_Data]'
!   write(0,*) '  Tests the equality d^T L e = e^T L^T d for any EMsoln and data.'
!   write(0,*) '  Optionally, outputs L e and L^T d.'
!   write(0,*)
!   write(0,*) ' -A  S rFile_EMrhs [wFile_EMsoln]'
!   write(0,*) '  Tests the equality b^T S^{-1} b = b^T (S^{-1})^T b for any EMrhs.'
!   write(0,*) '  For simplicity, use one EMrhs for forward and transpose solvers.'
!   write(0,*) '  Optionally, outputs e = S^{-1} b.'
!   write(0,*)
!   write(0,*) ' -A  P rFile_Model rFile_EMsoln [wFile_Model wFile_EMrhs]'
!   write(0,*) '  Tests the equality e^T P m = m^T P^T e for any EMsoln and data.'
!   write(0,*) '  Optionally, outputs P m and P^T e.'
!   write(0,*)
!   write(0,*) ' -A  Q rFile_Model rFile_Data [wFile_Model wFile_Data]'
!   write(0,*) '  Tests the equality d^T Q m = m^T Q^T d for any model and data.'
!   write(0,*) '  Optionally, outputs Q m and Q^T d.'
!   write(0,*)
!   write(0,*) 'Finally, generates random 5% perturbations, if implemented:'
!   write(0,*) ' -A  m rFile_Model wFile_Model [delta]'
!   write(0,*) ' -A  d rFile_Data wFile_Data [delta]'
!   write(0,*) ' -A  e rFile_EMsoln wFile_EMsoln [delta]'
!   write(0,*) ' -A  b rFile_EMrhs wFile_EMrhs [delta]'
!
!    write(*,*) ' -S  rFile_Model rFile_dModel rFile_Data wFile_Data [wFile_Sens]'
!    write(*,*) '  Multiplies by the full Jacobian, row by row, to get d = J m.'
!    write(*,*) '  Compare to the output of [MULT_BY_J] to test [COMPUTE_J]'

!  Before any of these routines may be called, the transmitter (txDict),
!  data type (typeDict) and receiver (rxDict) dictionaries must be
!  created and initialized. These are heavily used by Level II routines
!  in DataFunc, SolverSens and ForwardSolver, inherited by this module;
!  "pointers" to dictionary entries are attached to data vector d.

  use SensComp

  implicit none

  public    :: Jtest, Ltest, Stest, Ptest, Qtest

Contains

  !**********************************************************************
   subroutine Jtest(m0,m,d)

   !  Tests the equality d^T J m = m^T J^T d; outputs J m and J^T d.
   !
   !   m0 is background model parameter
   type(modelParam_t), intent(in)           :: m0
   !   m is the input model parameter perturbation;
   !    on output, this is J^T d
   type(modelParam_t), intent(inout)        :: m
   !   d is the input perturbation in the data vector;
   !    on output, this is J m
   type(dataVectorMTX_t), intent(inout)     :: d

   !  local variables
   type(modelParam_t)      :: JTd
   type(dataVectorMTX_t)   :: dPred,Jm
   type(solnVectorMTX_t)   :: eAll
   real(kind=prec)         :: r1,r2

   write(*,*) 'Symmetry test for operators J and J^T'

   ! initialize
   dPred = d
   Jm = d

   ! compute background solution
   call fwdPred(m0,dPred,eAll)

   ! compute J m
   call Jmult(m,m0,Jm,eAll)

   ! compute J^T d
   call JmultT(m0,d,JTd,eAll)

   ! compute dot product #1: d^T J m
   r1 = dotProd(d,Jm)
   write(*,*) 'dot product #1: ',r1,'[d^T J m]'

   ! compute dot product #2: m^T J^T d
   r2 = dotProd(m,JTd)
   write(*,*) 'dot product #2: ',r2,'[m^T J^T d]'

   ! output the results
   m = JTd
   d = Jm

   ! clean up
   call deall(dPred)
   call deall(eAll)
   call deall(JTd)
   call deall(Jm)

  end subroutine Jtest

  !**********************************************************************
   subroutine Ltest(m0,e,d,ePred)

   !  Tests the equality d^T L e = e^T L^T d; outputs L e and L^T d.
   !
   !   m0 is background model parameter
   type(modelParam_t), intent(in)           :: m0
   !   e is the input solution vector perturbation;
   !    on output, this is L^T d
   type(solnVectorMTX_t), intent(inout)     :: e
   !   d is the input perturbation in the data vector;
   !    on output, this is L e
   type(dataVectorMTX_t), intent(inout)     :: d
   !   eAll is the background solution vector; optional
   type(solnVectorMTX_t), optional, intent(in) :: ePred

   !  local variables
   type(rhsVector_t)       :: LTd
   type(dataVectorMTX_t)   :: dPred,Le
   type(solnVectorMTX_t)   :: eAll
   real(kind=prec)         :: r1
   complex(kind=prec)      :: c2
   integer                 :: j
   logical                 :: Conj_case

   write(*,*) 'Symmetry test for operators L and L^T'

   ! initialize
   Le = d
   Conj_case = .false.

   ! compute background solution
   if (.not. present(ePred)) then
    dPred = d
    call fwdPred(m0,dPred,eAll)
   else
    eAll = ePred
   endif

   do j = 1,d%nTx

       ! compute L e
       call Lmult(eAll%solns(j),m0,e%solns(j),Le%d(j))

       ! compute L^T d
       LTd = eAll%solns(j)
       call LmultT(eAll%solns(j),m0,d%d(j),LTd)

       ! compute dot product #1: d^T L e
       r1 = dotProd(d%d(j),Le%d(j))
       write(*,'(a15,i4,a17)',advance='no') '... frequency #',j,' dot product #1: '
       write(*,*) r1,'[d^T L e]'

       ! compute dot product #2: e^T L^T d
       c2 = dotProd(LTd,e%solns(j),Conj_case)
       write(*,'(a15,i4,a17)',advance='no') '... frequency #',j,' dot product #2: '
       write(*,*) c2,'[e^T L^T d]'

       ! convert RHS to solution vectors for output
       !e%solns(j) = LTd

   enddo  ! tx

   ! clean up
   call deall(dPred)
   call deall(eAll)
   call deall(LTd)
   call deall(Le)

  end subroutine Ltest

  !**********************************************************************
   subroutine Stest(m0,b,eAll)

   !  Tests the equality b^T S^{-1} b = b^T (S^{-1})^T b; outputs e = S^{-1} b.
   !
   !   m0 is background model parameter
   type(modelParam_t), intent(in)           :: m0
   !   b is the input forcing
   type(rhsVectorMTX_t), intent(inout)      :: b
   !   e is the output solution S^{-1} b
   type(solnVectorMTX_t), intent(inout)     :: eAll

   !  local variables
   type(solnVector_t)      :: eFwd,eAdj
   type(rhsVector_t)       :: comb
   real(kind=prec)         :: r1,r2
   integer                 :: j,iTx
   logical                 :: Conj_case

   write(*,*) 'Symmetry test for operators S^{-1} and (S^{-1})^T'

   ! initialize
   call create(b%nTx,eAll)
   Conj_case = .false.

   do j = 1,b%nTx

      !  initialize the temporary data vectors
      iTx = b%combs(j)%tx

      !  manage any necessary initilization for this transmitter
      call initSolver(iTx,m0,b%combs(j)%grid,eFwd)

      ! solve forward problem with source in comb, and save for output
      call sensSolve(iTx,FWD,eFwd,b%combs(j))
      eAll%solns(j) = eFwd

      !  now, initialize for sensitivity ... do nothing to comb, though!
      call initSolver(iTx,m0,b%combs(j)%grid,eAdj)

      ! solve transpose problem with source in comb
      call sensSolve(iTx,ADJ,eAdj,b%combs(j))

      ! convert RHS to full vector for dot product
      comb = b%combs(j)
      call sparse2full_rhsVector(comb)

      ! compute dot product #1: b^T S^{-1} b
      r1 = dotProd(comb,eFwd,Conj_case)
      write(*,'(a15,i4,a17)',advance='no') '... frequency #',j,' dot product #1: '
      write(*,*) r1,'[b^T S^{-1} b]'

      ! compute dot product #2: b^T (S^{-1})^T b
      r2 = dotProd(comb,eAdj,Conj_case)
      write(*,'(a15,i4,a17)',advance='no') '... frequency #',j,' dot product #2: '
      write(*,*) r2,'[b^T (S^{-1})^T b]'

      ! clean up
      call deall_rhsVector(comb)

   enddo  ! tx

   ! clean up
   call exitSolver()
   call deall(eFwd)
   call deall(eAdj)
   call deall(comb)

  end subroutine Stest


  !**********************************************************************
   subroutine Ptest(m0,dTemplate,m,e,ePred)

   !  Tests the equality e^T P m = m^T P^T e; outputs P m and P^T e.
   !
   !   m0 is background model parameter
   type(modelParam_t), intent(in)           :: m0
   !   dTemplate is the template data vector
   type(dataVectorMTX_t), intent(in)        :: dTemplate
   !   m is the input model parameter perturbation;
   !    on output, this is P^T d
   type(modelParam_t), intent(inout)        :: m
   !   e is the input perturbation in the solution vector;
   !    on output, this is P m
   type(solnVectorMTX_t), intent(inout)     :: e
   !   eAll is the background solution vector; optional
   type(solnVectorMTX_t), optional, intent(in) :: ePred

   !  local variables
   type(modelParam_t)      :: PTe,mTemp
   type(rhsVector_t)       :: Pm
   type(dataVectorMTX_t)   :: dPred
   type(solnVectorMTX_t)   :: eAll
   real(kind=prec)         :: r1,r2
   integer                 :: j
   logical                 :: Conj_case

   write(*,*) 'Symmetry test for operators P and P^T'

   ! initialize
   mTemp = m
   call zero(mTemp)
   Conj_case = .false.

   ! compute background solution
   if (.not. present(ePred)) then
    dPred = dTemplate
    call fwdPred(m0,dPred,eAll)
   else
    eAll = ePred
   endif

   do j = 1,e%nTx

       ! compute P m
       Pm = eAll%solns(j)
       call Pmult(eAll%solns(j),m0,m,Pm)

       ! compute P^T d
       PTe = m0
       call PmultT(eAll%solns(j),m0,e%solns(j),PTe)

       ! compute dot product #1: e^T P m
       r1 = dotProd(Pm,e%solns(j),Conj_case)
       write(*,'(a15,i4,a17)',advance='no') '... frequency #',j,' dot product #1: '
       write(*,*) r1,'[e^T P m]'

       ! compute dot product #2: m^T P^T e
       r2 = dotProd(m,PTe)
       write(*,'(a15,i4,a17)',advance='no') '... frequency #',j,' dot product #2: '
       write(*,*) r2,'[m^T P^T e]'

       ! sum up models for output
       call scMultAdd(ONE,PTe,mTemp)

       ! output the results
       !e%solns(j) = Pm

   enddo  ! tx

   ! output the results
   m = mTemp


   ! clean up
   call deall(dPred)
   call deall(eAll)
   call deall(PTe)
   call deall(Pm)
   call deall(mTemp)

  end subroutine Ptest


  !**********************************************************************
   subroutine Qtest(m0,m,d)

   !  Tests the equality d^T Q m = m^T Q^T d; outputs Q m and Q^T d.
   !
   !   m0 is background model parameter
   type(modelParam_t), intent(in)           :: m0
   !   m is the input model parameter perturbation;
   !    on output, this is Q^T d
   type(modelParam_t), intent(inout)        :: m
   !   d is the input perturbation in the data vector;
   !    on output, this is Q m
   type(dataVectorMTX_t), intent(inout)     :: d

   !  local variables
   type(modelParam_t)      :: QTd,mTemp
   type(dataVectorMTX_t)   :: dPred,Qm
   type(solnVectorMTX_t)   :: eAll
   real(kind=prec)         :: r1,r2
   integer                 :: j

   write(*,*) 'Symmetry test for operators Q and Q^T'

   ! initialize
   dPred = d
   Qm = d
   mTemp = m
   call zero(mTemp)

   ! compute background solution
   call fwdPred(m0,dPred,eAll)

   do j = 1,d%nTx

	   ! compute Q m
	   call Qmult(eAll%solns(j),m0,m,Qm%d(j))

	   ! compute Q^T d
	   call QmultT(eAll%solns(j),m0,d%d(j),QTd)

	   ! compute dot product #1: d^T Q m
	   r1 = dotProd(d%d(j),Qm%d(j))
       write(*,'(a15,i4,a17)',advance='no') '... frequency #',j,' dot product #1: '
	   write(*,*) r1,'[d^T Q m]'

	   ! compute dot product #2: m^T Q^T d
	   r2 = dotProd(m,QTd)
       write(*,'(a15,i4,a17)',advance='no') '... frequency #',j,' dot product #2: '
	   write(*,*) r2,'[m^T Q^T d]'

	   ! sum up models for output
	   call scMultAdd(ONE,QTd,mTemp)

   enddo  ! tx

   ! output the results
   m = mTemp
   d = Qm

   ! clean up
   call deall(dPred)
   call deall(eAll)
   call deall(QTd)
   call deall(Qm)
   call deall(mTemp)

  end subroutine Qtest


  !**********************************************************************

end module SymmetryTest
