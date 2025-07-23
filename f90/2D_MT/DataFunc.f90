! *****************************************************************************
module datafunc
  ! 2D MT data functionals ... now for TE + TM impedance
  !
  ! This module contains
  !  (1) routines for evaluation of impedances, and ultimately other
  !       interpretation parameters
  !  (2) routines to compute data functionals for linearized
  !       impedances,  and ultimately other interpretation paramters
  !   The idea:
  !     -> first the dictionaries txDict, typeDict, rxDict are initialized
  !         by calling appropriate initialization/setup routines
  !     -> data are stored in structures (defined in module DataSpace)
  !        which contain indices into transmitter and receiver dictionaries
  !        in addition to actual data values.  These indices are used by
  !        the data functional computation routines to compute predicted data.
  !
  !  This module is specific to 2D MT; similar modules will need to be written
  !     to implement data functionals for other problems

  use EMfieldInterp     !  basic interpolation routines for 2D TE and TM
                        !    solution grids; allow computation of both E and B
                        !    at an arbitrary point in either grid
  use ModelSpace, only:	rhoC => ModelParamToOneCell	!  model parameterization
					!  dependent function needed
					!  for TM mode data
					!  functionals
  use SolnSpace
  use transmitters
  use receivers
  use dataTypes

  implicit none

  !   Names of these routines must be as here, as these are called by
  !    top-level inversion routines
  public			:: dataResp, Lrows, Qrows


Contains

!******************************************************************************
  subroutine dataResp(ef,Sigma,iDT,iRX,Resp)

  !   2D MT impedance; gets mode from solution
  ! given electric field solution  (stored as type cvector)
  !   + indices into receiver dictionary
  ! compute complex scalar imepdance
  !   This now creates the required sparse vectors
  !      (either for TE or TM, as appropriate) for impedance evaluation

  type (solnVector_t), intent(in)			:: ef
  ! model parameter used to computed ef
  type (modelParam_t), intent(in)   :: Sigma
  ! indicies into data type and receiver dictionaries
  integer, intent(in)				:: iDt, iRX
  !  real data returned; note that this will always be real (even when
  !   data are complex ... then real and imag parts come in pairs)
  real(kind=prec), intent(inout)    :: Resp(2)
  !  Z will always be an array in the calling
  !    routine (as in top-level inversion routines) and so needs to
  !    be treated as an array here, even if there is only one element.
  !   As an example: to add tippers to TE mode, dimension on Z
  !     will have to be changed to 2 (and of Resp to 4)!
  complex(kind=prec)	:: Z(1),tempZ(1)

  !  local variables
  type(sparsevecc)		:: Lb,Le,Lbz,Lby
  complex(kind=prec)	:: B,E,Bz,By
  real(kind=prec)	:: omega, x(2)
  logical			:: Conj_Case = .false.
  character(2)			:: mode
  character*80			:: msg
  integer               :: yz  !define which magnetic componenets is required (By=1, Bz=2)
  !  get mode, frequency for transmitter used to compute solution ef
  mode = txDict(ef%tx)%mode
  omega =  txDict(ef%tx)%omega
  ! get location from receiver dictionary
  x = rxDict(iRX)%x
  if( mode.eq.'TE') then
     ! electric field Ex
     call NodeInterpSetup2D(ef%grid,x,mode,Le)
     ! magnetic field By
     yz=1
     call BfromESetUp_TE(ef%grid,x,omega,yz,Lb)
     
     yz=1
     call BfromESetUp_TE(ef%grid,x,omega,yz,Lby)
     ! magnetic field Bz (required if TP is needed)
     yz=2
     call BfromESetUp_TE(ef%grid,x,omega,yz,Lbz)
  elseif(mode .eq.'TM') then
     ! magnetic field Bx
     call NodeInterpSetup2D(ef%grid,x,mode,Lb)
     ! electric field Ey
     call EfromBSetUp_TM(ef%grid,x,omega,Sigma,Le)
  else
     call errStop('option not available in dataResp')
  endif
  
 if(iDt .eq. Tzy_Impedance) then
        ! error checking
       if (.not. Lbz%allocated) then
        write(0,*) 'Sparse vector Lb not allocated in dataResp for tx=',ef%tx,' rx=',iRX
       endif
       ! magnetic field By
       By = dotProd_scvector(Lby,ef%vec,Conj_case)     
       ! magnetic field Bz
       Bz = dotProd_scvector(Lbz,ef%vec,Conj_case) 
       Z(1) = Bz/By 
 else      
      ! error checking
      if (.not. Le%allocated) then
        write(0,*) 'Sparse vector Le not allocated in dataResp for tx=',ef%tx,' rx=',iRX
      endif
      if (.not. Lb%allocated) then
        write(0,*) 'Sparse vector Lb not allocated in dataResp for tx=',ef%tx,' rx=',iRX
      endif
      !  Using sparse vector representations of data functionals,
      !          compute impedance
      E = dotProd_scvector(Le,ef%vec,Conj_case)
      ! magnetic field
      B = dotProd_scvector(Lb,ef%vec,Conj_case)

      ! impedance is trivial for 2D!
      Z(1) = E/B
      ! Now, check if we are working with Rho and Phase, and overwrite Z if needed
      if (iDt .eq. Rho_Phs_TM ) then
          Z(1)  = dcmplx((abs(Z(1))**2*MU_0/omega), atan2(ISIGN*dimag(Z(1)),real(Z(1)))*R2D+180.0d0)
      end if
 end if
 
  ! convert to real output
  Resp(1) = real(Z(1))
  Resp(2) = imag(Z(1))

  ! clean up
   if(iDt .eq. Tzy_Impedance) then
         call deall_sparsevecc(Lby)
         call deall_sparsevecc(Lbz)
   else    
         call deall_sparsevecc(Lb)
         call deall_sparsevecc(Le)
   end if
   

  end subroutine dataResp

!****************************************************************************
  subroutine Lrows(e0,Sigma0,iDT,iRX,Lz)
  !  given input background electric field solution,
  !  index into receiver dictionary for a single site (iRX)
  !  compute sparse complex vector giving coefficients
  !  of linearized impedance functional (complex representation)
  !  For TM mode solution also returns sparse vector Q (model
  !     paramter space) for derivative of data functional with
  !     respect to model paramters; Q is not referenced for TE data
  !   NOTE: Lz has to be declared as an array for
  !     consistency with calling program (in general the
  !     number nFunc of complex data functionals that will
  !     be returned by data functional routines could be > 1)

  !  electric field solutions are stored as type solnVector
  type (solnVector_t), intent(in)		:: e0
  ! model parameter used to computed e0
  type (modelParam_t), intent(in)   :: Sigma0
  ! indicies into data type and receiver dictionaries
  integer, intent(in)			:: iRX,iDT
  !  Lz, Qz will always be arrays in the calling
  !    routine (as in top-level inversion routines) and so need to
  !    be treated as arrays here, even if there is only one element.
  !   As an example: to add tippers to TE mode, dimension on LZ will
  !    have to be changed to 2.
  type(sparseVector_t), intent(inout)		:: Lz(1)

  !  local variables
  complex (kind=prec)		:: B,E,c_E,c_B,By,Bz,c_By,c_Bz,Z
  type(sparsevecc)			:: Le,Lb,Lby,Lbz
  real (kind=prec)		:: x(2), omega
  character(2)				:: mode
  logical				:: Conj_case = .false.
  integer                :: yz

  !  get mode, frequency for transmitter used to compute solution ef
  mode = txDict(e0%tx)%mode
  omega =  txDict(e0%tx)%omega
  ! get location from reciever dictionary
  x = rxDict(iRX)%x

  !   evaluate E, B, for background solution
if(mode .eq.'TE') then
     ! electric field
     call NodeInterpSetup2D(e0%grid,x,mode,Le)
     ! magnetic field
     yz=1
     call BfromESetUp_TE(e0%grid,x,omega,yz,Lb)
     
     ! magnetic field By (its repetation to the previous Interpolation, but, saved in Lby. Leave it like this for now)
     yz=1
     call BfromESetUp_TE(e0%grid,x,omega,yz,Lby)
     ! magnetic field Bz (required if TP is needed)
     yz=2
     call BfromESetUp_TE(e0%grid,x,omega,yz,Lbz)
elseif(mode .eq.'TM') then
     ! magnetic field
     call NodeInterpSetup2D(e0%grid,x,mode,Lb)
     ! electric field
     call EfromBSetUp_TM(e0%grid,x,omega,Sigma0,Le,e0%vec)
  else
     call errStop('option not available in Lrows')
  endif
  
if(iDt .eq. Tzy_Impedance) then
      Bz = dotProd_scvector(Lbz,e0%vec,Conj_case)  
      By = dotProd_scvector(Lby,e0%vec,Conj_case)


      c_By = C_ONE/By
      c_Bz = -Bz/(By*By)
      call linComb_sparsevecc(Lbz,c_by,Lby,c_Bz,Lz(1)%L)
      call deall_sparsevecc(Lby)
      call deall_sparsevecc(Lbz)  
else      
      !  Compute electric, magnetic field for background soln.
      E = dotProd_scvector(Le,e0%vec,Conj_case)
      B = dotProd_scvector(Lb,e0%vec,Conj_case)

      !  compute sparse vector representations of linearized
      !    impedance functional; coefficients depend on E, B
      c_E = C_ONE/B
      c_B = -E/(B*B)
      !  Nominally, Lz = c_E*Le + c_B*Lb
      ! However, note that Lz is of type sparseVector (here just a wrapped
      !    version of a sparsevecc object)
      call linComb_sparsevecc(Le,c_E,Lb,c_B,Lz(1)%L)
     
      ! Now, check if we are working with Rho and Phase: Not ready yet.
      if (iDt .eq. Rho_Phs_TM) then
          
      end if    
      ! clean up
      call deall_sparsevecc(Le)
      call deall_sparsevecc(Lb)
end if

  end subroutine Lrows
!
!****************************************************************************
  subroutine Qrows(e0,Sigma0,iDT,iRX,zeroValued,Qreal,Qimag)
  !  given input background solution vector (e0) and model parameter (Sigma0)
  !  and indices into data type and receiver dictionaries
  !  compute derivative of data functional with respect to model parameters
  !  for all components of the data type ...
  !  For TM mode solution returns sparse vector Q (model
  !     parameter space) for derivative of data functional with
  !     respect to model parameters; Q is not referenced for TE data

  type (solnVector_t), intent(in)          :: e0
  type (modelParam_t), intent(in)      :: Sigma0
  integer, intent(in)                        :: iDT, iRX
  logical, intent(out)                      :: zeroValued
  !   NOTE: Qreal and Qimag have to be declared as arrays for
  !     consistency with calling program (in general the
  !     number nFunc of complex data functionals that will
  !     be returned by data functional routines could be > 1)
  !   NOTE: Qreal and Qimag both exist regardless of whether the data
  !     are real or complex, since Q itself is complex
  type(modelParam_t), intent(inout)     :: Qreal(:), Qimag(:)

  ! local variables
  complex (kind=prec)   :: B,c_E
  type(sparsevecc)      :: Lb,Le
  real (kind=prec)      :: x(2), omega
  character(2)           :: mode
  logical               :: Conj_case = .false.
  logical               :: isComplex = .true.
  integer               :: nFunc = 1
  type(sparseVector_t)  :: Qz(1)
  integer       :: istat,ncomp,iFunc

  ! check for incorrect usage
  ncomp = typeDict(iDT)%nComp
  if((ncomp .ne. 2*nFunc) .or. (typeDict(iDT)%isComplex .neqv. isComplex)) then
    call errStop('incorrect usage of type dictionary in Qrows')
  endif
  if((size(Qreal) .ne. nFunc) .or. (size(Qimag) .ne. nFunc)) then
    call errStop('incorrect output size in Qrows')
  endif

  !  get mode, frequency from transmitter dictionary
  mode = txDict(e0%tx)%mode
  omega =  txDict(e0%tx)%omega
  ! get location from receiver dictionary
  x = rxDict(iRX)%x

  if(mode .eq. 'TE') then
    ! for efficiency, just set the logical to zero and exit
    zeroValued = .true.
	!do iFunc = 1, nFunc
	!  Qreal(iFunc) = Sigma0
	!  call zero(Qreal(iFunc))
    !  Qimag(iFunc) = Sigma0
    !  call zero(Qimag(iFunc))
	!enddo
  elseif(mode.eq. 'TM' ) then
     ! set the logical to false - vectors are non-zero.
     zeroValued = .false.
     ! compute magnetic field for background solution
     call NodeInterpSetup2D(e0%grid,x,mode,Lb)
     ! electric field
     call EfromBSetUp_TM(e0%grid,x,omega,Sigma0,Le,e0%vec,Qz(1)%L)
     B = dotProd_scvector(Lb,e0%vec,Conj_case)
     c_E = C_ONE/B
     !  also need to multiply parameter space sparse vector
     !    (derivative of Le coefficients wrt cell resistivities)
     !   by 1/B  (actually 1/B == 1 !)
     Qz(1)%L%C = c_E*Qz(1)%L%C
     !  map from the grid to the model parameter
     call SparseCellToModelParam(Qz(1)%L,Sigma0,Qreal(1),Qimag(1))
     ! clean up
     call deall_sparsevecc(Le)
     call deall_sparsevecc(Lb)
  else
     call errStop('option not available in Qrows')
  endif

  end subroutine Qrows

end module dataFunc
