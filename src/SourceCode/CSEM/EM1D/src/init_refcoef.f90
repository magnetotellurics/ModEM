!---------------------------------------------------------
! EM1D subroutine init_refcoef
!
! allocate matrices and vectors for reflection coefficients
!   can be optimized for memory by checking in advance if we really need both TE and TM
!
! Rita Streich 2010
!---------------------------------------------------------
subroutine init_refcoef()

  implicit none

  !no external variables

  !internal variables
  integer(kind=int32)   :: ierr  !error index


!-----------------------------------------------------------------
!  prepare interface reflection and transmission coefficient precomputation
!-----------------------------------------------------------------

  allocate(rupallTE(nlay-1,0:filtlen),rupallTM(nlay-1,0:filtlen),rdnallTE(nlay-1,0:filtlen),rdnallTM(nlay-1,0:filtlen), &
           tupallTE(nlay-1,0:filtlen),tupallTM(nlay-1,0:filtlen),tdnallTE(nlay-1,0:filtlen),tdnallTM(nlay-1,0:filtlen), &
           pvertall1(nlay,0:filtlen),pvertall2(nlay,0:filtlen), stat=ierr)
  if(ierr.ne.0) call alloc_error(pid,'init_refcoef','interface reflection and transmission coeff. vectors',ierr)


  ! for each interface --> nr of layers - 1, will be used for both TE and TM mode
!!$  allocate(rup(nlay-1),rdn(nlay-1), tup(nlay-1),tdn(nlay-1), stat=ierr)
!!$  if(ierr.ne.0) call alloc_error(pid,'init_refcoef','interface reflection and transmission coeff. vectors',ierr)


endsubroutine init_refcoef



!---------------------------------------------------------
! EM1D subroutine prepare_refcoef
!
! precompute interface reflection and transmission coefficients
!   for a given max. radius
!
! Rita Streich 2010-2011, Alexander Grayver 2014
!---------------------------------------------------------
subroutine prepare_refcoef(refl_var,r,dipoletype,aniso)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  real(kind=real64),intent(in)    :: r             !base radius
  integer(kind=int32),intent(in)  :: dipoletype    !source type (HED, VED etc, input zero for wire)
  integer(kind=int32),intent(in)  :: aniso         !anisotropy index


  !internal variables
  real(kind=real64)             :: kap0,kappa,kappasq,phsq  !horiz. wavenumber, freq-normalized horiz. wavenumber squared
  integer(kind=int32)           :: fl,i      !filter length
  real(kind=real64)             :: ABSCIS  !constant used for defining wavenumber spacing
  real(kind=real64)             :: E,ER    !constants used for defining wavenumber spacing
  integer(kind=int32)           :: ilay        !layer counter
  integer(kind=int32)           :: ikap        !wavenumber counter
  complex(kind=real128) :: a,b,c,d,f
  DATA E/1.10517091807564762D0/,ER/.904837418035959573D0/ !factor between wavenumbers
  DATA ABSCIS/0.7059431685223780D0/        !factor for determining central wavenumber, taken from Hankel routine


if(refl_var%refcoef_changed) then

 kap0 = ABSCIS / r
 kappasq = kap0**2
 phsq = kappasq / omegasq
 isiso: if(aniso .eq. iso) then
  !get vertical wavenumbers first - use the same wavenumbers as inside Hankel routine
  pvertall1(:,298) = sqrt(epsmuv(:) - phsq)

  kappa = kap0
  do ikap=299,filtlen
    kappa = kappa * E
    kappasq = kappa**2
    phsq = kappasq / omegasq
    pvertall1(:,ikap) = sqrt(epsmuv(:) - phsq)
  enddo
  kappa = kap0
  do ikap=297,1,-1
    kappa = kappa * ER
    kappasq = kappa**2
    phsq = kappasq / omegasq
    pvertall1(:,ikap) = sqrt(epsmuv(:) - phsq)
  enddo
  
  !loop to prevent stack overflow
  do ilay=1,nlay
    pvertall2(ilay,:) = pvertall1(ilay,:)
  enddo

 else !vti

  !get vertical wavenumbers first - use the same wavenumbers as inside Hankel routine
  pvertall1(:,298) = sqrt(epsmuh(:) - phsq)
  pvertall2(:,298) = sqrt(epsmuh(:) - epsmuratio(:)*phsq)

  kappa = kap0
  do ikap=299,filtlen
    kappa = kappa * E
    kappasq = kappa**2
    phsq = kappasq / omegasq
    pvertall1(:,ikap) = sqrt(epsmuh(:) - phsq)
    pvertall2(:,ikap) = sqrt(epsmuh(:) - epsmuratio(:)*phsq)
  enddo
  kappa = kap0
  do ikap=297,1,-1
    kappa = kappa * ER
    kappasq = kappa**2
    phsq = kappasq / omegasq
    pvertall1(:,ikap) = sqrt(epsmuh(:) - phsq)
    pvertall2(:,ikap) = sqrt(epsmuh(:) - epsmuratio(:)*phsq)
  enddo

  !fix against acausal propagation (OK?????)
  do ilay=1,nlay
    pvertall2(ilay,:) = cmplx(real(pvertall2(ilay,:),kind=real64),abs(aimag(pvertall2(ilay,:))),kind=real64)
  enddo

 endif isiso

  fl = filtlen

  !source type VED: TM mode only
  select case (dipoletype)
  case (ved)

    do ilay = 1,nlay-1
      rupallTM(ilay,1:fl) = (1._real64 - (epsh(ilay)*pvertall2(ilay+1,1:fl)) / (epsh(ilay+1)*pvertall2(ilay,1:fl))) / &
                            (1._real64 + (epsh(ilay)*pvertall2(ilay+1,1:fl)) / (epsh(ilay+1)*pvertall2(ilay,1:fl)))
	!WARNING: minus is not in formula 121 of Løseth and Ursin, but setting a minus here
	!  makes field values reasonable for receivers in different layer than source
	!  (checked on 2-layer model: integrals A1TE and A1TM have to nearly cancel out)
      if(pvertall2(ilay,298).eq.pvertall2(ilay+1,298)) then
        tupallTM(ilay,1:fl) = 1._real64
      else
        tupallTM(ilay,1:fl) = -2._real64*sqrt(epsh(ilay)*epsh(ilay+1)*pvertall2(ilay,1:fl)*pvertall2(ilay+1,1:fl)) / &
                      (epsh(ilay+1)*pvertall2(ilay,1:fl) + epsh(ilay)*pvertall2(ilay+1,1:fl))
      endif

      rdnallTM(ilay,1:fl) = -rupallTM(ilay,1:fl)
      tdnallTM(ilay,1:fl) = tupallTM(ilay,1:fl)
    enddo


  !source type VMD: TE mode only
  case (vmd)

    do ilay = 1,nlay-1
      rupallTE(ilay,1:fl) = (1._real64 - pvertall1(ilay,1:fl)/pvertall1(ilay+1,1:fl)) / (1 + pvertall1(ilay,1:fl)/pvertall1(ilay+1,1:fl))
      if(pvertall1(ilay,298).eq.pvertall1(ilay+1,298)) then
        tupallTE(ilay,1:fl) = 1._real64
      else
        tupallTE(ilay,1:fl) = 2._real64*sqrt(pvertall1(ilay,1:fl)*pvertall1(ilay+1,1:fl)) / (pvertall1(ilay+1,1:fl) + pvertall1(ilay,1:fl))
      endif
      
      rdnallTE(ilay,1:fl) = -rupallTE(ilay,1:fl)
      tdnallTE(ilay,1:fl) = tupallTE(ilay,1:fl)
    enddo


  !other source types: both TE and TM mode
  case default

    do ilay = 1,nlay-1
      rupallTM(ilay,1:fl) = -(1._real64 - (epsh(ilay+1)*pvertall2(ilay,1:fl)) / (epsh(ilay)*pvertall2(ilay+1,1:fl)))/ &
                            (1._real64 + (epsh(ilay+1)*pvertall2(ilay,1:fl)) / (epsh(ilay)*pvertall2(ilay+1,1:fl)))
        !rupallTM(ilay,1:fl) = (epsh(ilay+1)*pvertall2(ilay,1:fl)-epsh(ilay)*pvertall2(ilay+1,1:fl))/(epsh(ilay+1)*pvertall2(ilay,1:fl)+epsh(ilay)*pvertall2(ilay+1,1:fl))
      if(pvertall2(ilay,298).eq.pvertall2(ilay+1,298)) then
        tupallTM(ilay,1:fl) = 1._real64
      else
        tupallTM(ilay,1:fl) = -2._real64*sqrt(epsh(ilay)*epsh(ilay+1)*pvertall2(ilay,1:fl)*pvertall2(ilay+1,1:fl)) / &
                      (epsh(ilay+1)*pvertall2(ilay,1:fl) + epsh(ilay)*pvertall2(ilay+1,1:fl))
      endif

      rupallTE(ilay,1:fl) = (1._real64 - pvertall1(ilay,1:fl)/pvertall1(ilay+1,1:fl)) / (1._real64 + pvertall1(ilay,1:fl)/pvertall1(ilay+1,1:fl))
      if(pvertall1(ilay,298).eq.pvertall1(ilay+1,298)) then
        tupallTE(ilay,1:fl) = 1._real64
      else
        tupallTE(ilay,1:fl) = 2._real64*sqrt(pvertall1(ilay,1:fl)*pvertall1(ilay+1,1:fl)) / (pvertall1(ilay+1,1:fl) + pvertall1(ilay,1:fl))
      endif

      rdnallTM(ilay,1:fl) = -rupallTM(ilay,1:fl)
      tdnallTM(ilay,1:fl) = tupallTM(ilay,1:fl)
      rdnallTE(ilay,1:fl) = -rupallTE(ilay,1:fl)
      tdnallTE(ilay,1:fl) = tupallTE(ilay,1:fl)
    enddo

  end select

  refl_var%refcoef_changed = .TRUE.

endif

endsubroutine prepare_refcoef

