!----------------------------------------------------------------------
!  Performs automatic calculation of Bessel transform to specified
!  relative and absolute error
!
!  Argument list:
!
!  BESR,BESI - real and imaginary parts returned by BESAUT
!  ORDER - order of the BEssel function
!  NL - lower limit for Gauss order to start computation
!  NU - upper limit for Gauss order
!     NU,NL01,...7 selects 3,7,15,31,63,127, and 255 point Gauss
!     quadrature between the zero crossings of the Bessel function
!  R - argument of the Bessel function
!  FUNCT - external FUNCTION to evaluate kernel function at
!          a given wavenumber
!          The call is of the form
!              Y = FUNCT(X)
!              YR = real(Y); YI = aimag(Y)
!          where X is supplied by BESQUD and YR and YI are the returned
!          real and imaginary parts. All other arguments to FUNCT must be 
!          passed in common blocks.
!  RERR,AERR - relative and absolute error for termination
!       BESAUT terminates when increasing the Gauss order does not
!       change the result by more than RERR or when the absolute error
!       is less than AERR or when a Gauss order of NU is reached.
!  NPCS - number of pieces into which each partial integraiton is divided.
!       Ordinarily set to one. For very small values of R where
!       the kernel function is appreciable only over the first few
!       loops of the Bessel function, NPCS may be increased to achieve
!       reasonable accuracy.
!  NEW - if NEW=1, the integrations are computed and saved at each Gauss
!       order. If NEW=2, previously computed integrands are used. Note
!       that ORDER, R, and NPCS must not be changed when setting NEW=2.
!  IERR - error parameter
!       IERR=0 -- normal return
!       IERR=1 -- result not accurate to RERR due to too low a Gauss
!       order or convergence not achoeved in BESTRN
!
!  Subroutines reqiured:
!      BESTRN,BESQUD,JBESS,PADECF,CF,ZEROJ,DOT
!
!  A. Chave IGRR/UCSD
!
!  converted to F90: RIta Streich 2008
!
!-------------------------------------------------------------------------
SUBROUTINE BESAUT(BESR,BESI,ORDER,NL,NU,R,FUNCT,RERR,AERR,NPCS,NEW,NSUM,XSUM,IERR)

  implicit none

!  COMMOM /TEST/ contains run statistics
!  NG=highest Gauss order used
!  NF=total function evaluations for all calls to BESTRN
!  NI=total number of partial integrands on last call
!  to BESTRN only
  COMMON /TEST/NG,NF,NI
  integer(kind=int32)    :: NG,NF,NI

  !external variables
  real(kind=real64)      :: BESR,BESI
  real(kind=real64)      :: ORDER
  integer(kind=int32)    :: NL,NU
  real(kind=real64)      :: R
  real(kind=real64)      :: RERR,AERR
  integer(kind=int32)    :: NPCS
  integer(kind=int32)    :: NEW
  integer(kind=int32)    :: NSUM
  real(kind=real64),dimension(NSUM) :: XSUM
  integer(kind=int32)    :: IERR

  complex(kind=real64),external    :: FUNCT

  !internal variables
  integer(kind=int32)              :: NW
  real(kind=real64)                :: OLDR,OLDI
  integer(kind=int32)              :: N


     NF=0
     IF (NL.GT.NU) THEN
       BESR=0.
       BESI=0.
       IERR=1
       RETURN
       ENDIF
     NW=MAX(NEW,1)
     CALL BESTRN(BESR,BESI,ORDER,NL,R,FUNCT,.1*RERR,.1*AERR,NPCS,XSUM,NSUM,NW,IERR)
     IF ((IERR.NE.0).AND.(NL.EQ.7)) THEN
       NG=NL
       RETURN
     ELSE
       OLDR=BESR
       OLDI=BESI
       DO 10 N=NL+1,NU
         CALL BESTRN(BESR,BESI,ORDER,N,R,FUNCT,.1*RERR,.1*AERR,NPCS,XSUM,NSUM,2,IERR)
         IF ((IERR.NE.0).AND.(N.EQ.7)) THEN
           !BESR=OLDR
           !BESI=OLDI
           NG=N
           RETURN
         ELSEIF ((ABS(BESR-OLDR).LE.RERR*ABS(BESR)+AERR).AND. &
                 (ABS(BESI-OLDI).LE.RERR*ABS(BESI)+AERR)) THEN
           NG=N
           RETURN
         ELSE
           OLDR=BESR
           OLDI=BESI
         ENDIF
10       CONTINUE
     ENDIF
     NG=7
     IERR=1
     RETURN

endsubroutine BESAUT


!----------------------------------------------------------------------
!
!  Computes Bessel transform of specified order defined as
!  Integral(FUNCT(X)*J-SUB-ORDER(X*R)*DX) from X=0 to infinity
!  Computation is achieved by integration between the asymptotic
!  zero crossings of the Bessel function using Gauss quadrature.
!  The resulting series of partial integrands is summed by calculating
!  the Pade approximants to speed up convergence.
!
!  Argument list:
!
!  BESR,BESI - real and imaginary parts returned by BESTRN
!  ORDER - order of the Bessel function
!  NG - number of Gauss points to use in the quadrature routine.
!       NG=1 through 7 selects 3,7,15,31,63,127,and 255 terms.
!  R - argument of the Bessel function
!  FUNCT - external routine to evaluate kernel function at argument
!      X. The call is of the form
!         Y = FUNCT(X)
!         YR = real(Y); YI = aimag(Y)
!      where YR and YI are the returned real and imaginary parts
!  RERR,AERR - specified relative and absolute error for the
!      calculation. The integration
!      terminates when an additional term does not change the
!      result by more than RERR*RESULT+AERR
!  NPCS - number of pieces into which each partial integrand is divided,
!      ordinarily set to one. For very small values of range where
!      the kernel function is appreciable only over the first few
!      loops of the Bessel function, NPCS may be increased to achieve
!      reasonable accuracy. Note that NPCS affects only the Pade
!      sum portion of the integration, over X(NSUM) to infinity.
!  XSUM - vector of values of the kernel argument of FUNCT for which
!      explicit calculation of the integral is desired, so that the
!      integral over 0 to XSUM(NSUM) is added to the integral over
!      XSUM(NSUM) to Infinity with the Pade method invoked only for
!      the latter. This allows the Pade summation method to be
!      overridden and some types of singularities to be handled.
!  NSUM - number of values in XSUM, may be zero.
!  NEW - determines method of kernal calculation
!      NEW=0 means calculate but do not save integrands
!      NEW=1 means calculate kernel by calling FUNCT and save kernel
!            times Bessel function
!      NEW=2 means use saved kernels times Bessel functions in
!            COMMON /BESINT/. Note tjat ORDER,R,NPCS,XSUM, and
!            NSUM may not be changed when setting NEW=2.
!  IERR - error parameter
!      0 normal return, integral convergent
!      1 means no convergence after NSTOP terms in the Pade sum
!
!  Subroutines required:
!    BESQUD,PADECF,CF,ZEROJ,DOT.JBESS
!
!  A. CHAVE IGPP/UCSD
!
!  NTERM is maximum number of Bessel function loops stored if NEW.NE.0
!  NSTOP is maximum number of Pade terms
!
!  Converted to F90: Rita Streich 2008
!
!---------------------------------------------------------------------------
SUBROUTINE BESTRN(BESR,BESI,ORDER,NG,R,FUNCT,RERR,AERR,NPCS,XSUM,NSUM,NEW,IERR)

  implicit none

!  COMMON /TEST/ contains run statistics
!  NGAUSS=highest Gauss order used
!  NF=total function evaluations for all calls to BESTRN
!  NI=total number of partial integrands used in PADECF on last call
!     to BESTRN
  COMMON /TEST/NGAUSS,NF,NI
  integer(kind=int32)      :: NGAUSS,NF,NI

!  IF NEW.gt.0, COMMON /BESINT/ contains the kernel function arguments
!    and values KARG(I,N) and KERN(I,N) for N=1 to NPS and I=1 to
!    2**NK(N)+1)-1.
!  KERN(I,N),KERN(I+1,N) contain the real and imaginary parts of the
!    kernel function for the corresponding argument value (I,N).
!  Each value of N corresponds to one entry in XSUM or to one interval
!    between zero crossings of the Bessel function divided by NPCS.
!  Note that up to NTERM values are stored. THe integration proceeds
!    for N.GT.NTERM as if NEW=0. NP is an index passed internally by
!    BESTRN.
  COMMON /BESINT/NK,NP,NPS,KARG,KERN

  !external variables
  real(kind=real64)        :: BESR,BESI
  real(kind=real64)        :: ORDER
  integer(kind=int32)      :: NG
  real(kind=real64)        :: R
  real(kind=real64)        :: RERR,AERR
  integer(kind=int32)      :: NPCS
  integer(kind=int32)      :: NSUM
  real(kind=real64),dimension(NSUM)  :: XSUM
  integer(kind=int32)      :: NEW
  integer(kind=int32)      :: IERR

  complex(kind=real64),external      :: funct

  !internal variables
  integer(kind=int32),parameter      :: NTERM=nmaxsave, NSTOP=nmaxsum
  integer(kind=int32),dimension(NTERM)    :: NK
  integer(kind=int32)                     :: NP,NPS !NPS: Number of values Saved
  real(kind=real64),dimension(255,NTERM)  :: KARG
  real(kind=real64),dimension(510,NTERM)  :: KERN
  real(kind=real64)                  :: LASTR,LASTI
  real(kind=real64),dimension(NSTOP) :: SR,SI
  integer(kind=int32)                :: NPO  !Old Number of saved values
  integer(kind=int32)                :: I,N
  integer(kind=int32)                :: NW,NPB,L
  real(kind=real64)                  :: A,B  !temp integration limits
  real(kind=real64)                  :: XINC,AA,BB
  real(kind=real64)                  :: SUMR,SUMI,XSUMR,XSUMI
  real(kind=real64)                  :: TERMR,TERMI,TR,TI,LR,LI


     IF (NEW.EQ.2) THEN
       NPO=NPS
     ELSE
       DO 5 I=1,NTERM
         NK(I)=0
5      CONTINUE
       NPS=0
       NPO=NTERM
     ENDIF
!  Check for trivial case
     IF ((ORDER.NE.0.0).AND.(R.EQ.0.)) THEN
       BESR=0.
       BESI=0.
       IERR=0
       RETURN
     ENDIF
     NI=0
     NW=NEW
     NP=1
     NPB=1
     L=1
     B=0.0
     SUMR=0.0
     SUMI=0.0
     XSUMR=0.0
     XSUMI=0.0

     IF (NSUM.GT.0) THEN
!  compute Bessel transform explicitly on (0,XSUM(NSUM))
       LASTR=0.0
       LASTI=0.0
       DO 10 N=1,NSUM
         IF ((NW.EQ.2).AND.(NP.GT.NPO)) NW=1
         IF (NP.GT.NTERM) NW=0
         A=B
         B=XSUM(N)
         CALL BESQUD(A,B,TERMR,TERMI,NG,NW,ORDER,R,FUNCT)

         XSUMR=XSUMR+TERMR
         XSUMI=XSUMI+TERMI
         IF ((ABS(XSUMR-LASTR).LE.RERR*ABS(XSUMR)+AERR).AND. &
             (ABS(XSUMI-LASTI).LE.RERR*ABS(XSUMI)+AERR)) THEN
           BESR=XSUMR
           BESI=XSUMI
           IERR=0
           NPS=MAX(NP,NPS)
           !RETURN
         ELSE
           NP=NP+1
           LASTR=XSUMR
           LASTI=XSUMI
         ENDIF
10       CONTINUE

!  Find first zero crossing of Bessel function beyond XSUM(NSUM)
15     CONTINUE
       IF (ZEROJ(NPB,ORDER).GT.XSUM(NSUM)*R) GOTO 20
       NPB=NPB+1
       GOTO 15
     ENDIF

!  Entry point for Pade summation of partial integrants
20   CONTINUE
     LASTR=0.0
     LASTI=0.0
     A=B
     B=ZEROJ(NPB,ORDER)/R

!  Calculate terms and sum with PADECF, quitting when convergence is
!  obtained
     DO 30 N=1,NSTOP
       IF (NPCS.EQ.1) THEN
         IF ((NW.EQ.2).AND.(NP.GT.NPO)) NW=1
         IF (NP.GT.NTERM) NW=0
         CALL BESQUD(A,B,TERMR,TERMI,NG,NW,ORDER,R,FUNCT)
         NP=NP+1

       ELSE
         TERMR=0.
         TERMI=0.
         LR=0.
         LI=0.
         XINC=(B-A)/dble(NPCS)
         AA=A
         BB=A+XINC
         DO 25 I=1,NPCS
           IF ((NW.EQ.2).AND.(NP.GT.NPO)) NW=1
           IF (NP.GT.NTERM) NW=0
           CALL BESQUD(AA,BB,TR,TI,NG,NW,ORDER,R,FUNCT)
           TERMR=TERMR+TR
           TERMI=TERMI+TI
           IF ((ABS(TERMR-LR).LE.RERR*ABS(TERMR)+AERR).AND. &
           (ABS(TERMI-LI).LE.RERR*ABS(TERMI)+AERR)) THEN
             GOTO 26
           ELSE
             LR=TERMR
             LI=TERMI
             AA=BB
             BB=BB+XINC
             NP=NP+1
           ENDIF
25       CONTINUE

       ENDIF

26     CONTINUE

       NI=NI+1
       SR(L)=TERMR
       SI(L)=TERMI
       CALL PADECF(SUMR,SUMI,SR,SI,L)
       IF ((ABS(SUMR-LASTR).LE.RERR*ABS(SUMR)+AERR).AND. &
           (ABS(SUMI-LASTI).LE.RERR*ABS(SUMI)+AERR)) THEN
         BESR=XSUMR+SUMR
         BESI=XSUMI+SUMI
         IERR=0
         NPS=MAX(NP-1,NPS)
         RETURN
       ELSE
         LASTR=SUMR
         LASTI=SUMI
         NPB=NPB+1
         A=B
         B=ZEROJ(NPB,ORDER)/R
         L=L+1
       ENDIF
30   CONTINUE
     BESR=XSUMR+SUMR
     BESI=XSUMI+SUMI
     IERR=1
     NPS=MAX(NP,NPS)
     RETURN

endsubroutine BESTRN



!------------------------------------------------------------------
!  This subroutine calculates the integral of F(X)*J-SUB-N(X*R)
!  over the interval A to B at a specified Gauss order.
!  The result is obtained using a sequence of 1,3,7,15,31,63,
!  127, and 255 point interlacing Gauss formulae so that no integrand
!  evaluations are wasted. The kernel functions may be saved so that
!  Bessel transforms of similar kernels are computed without new
!  evaluation of the kernel.
!  Details on the formulae are given in "The optimum addition of points
!  to quadrature formulae" by T.N.L Patterson, Maths. Comp.
!  22, 847-856 (1968). Gauss weights are taken from
!  Comm. A.C.M. 16,694-699 (1973).
!
!  Argument list:
!
!  A - lower limit of integration
!  B - upper limit of integration
!  BESR,BESI - returned integral value real and imaginary parts
!  NG - number of points in the Gauss formula. NG=1,...7
!       selects 3,7,15,31,63,127, and 255 point quadrature.
!  NEW - selects method of kernel evaluation
!      NEW=0 calculates kernels by calling F - nothing saved
!      NEW=1 calculates kernels by calling F and saves kernel times
!            Bessel function in COMMON /BESINT/
!      NEW=2 uses saved kernel times Bessel functions in
!            COMMON /BESINT/
!  ORDER - order of the Bessel function
!  R - argument of the Bessel function
!  F - F(X) is the external integrand subroutine
!
!  A. Chave IGPP/UCSD
!
!  Converted to F90: Rita Streich 2008
!------------------------------------------------------------------
SUBROUTINE BESQUD(A,B,BESR,BESI,NG,NEW,ORDER,R,F)

  implicit none

!  COMMON /TEST/ contains run statistics
!  NGAUSS=highest Gauss order used
!  NF=total function evaluations for all calls to BESTRN
!  NI=total number of partial integrands on last call
!     to BESTRN only
  COMMON /TEST/NGAUSS,NF,NI
  COMMON /BESINT/NK,NP,NPS,KARG,KERN

  !external variables
  real(kind=real64)         :: A,B
  real(kind=real64)         :: BESR,BESI
  integer(kind=int32)       :: NG
  integer(kind=int32)       :: NEW
  real(kind=real64)         :: ORDER
  real(kind=real64)         :: R

  complex(kind=real64),external :: F

  !internal variables
   ! maximum number of Bessel function loops that can be saved
  integer(kind=int32)              :: NGAUSS,NF,NI
  integer(kind=int32),parameter    :: NTERM=nmaxsave
  integer(kind=int32),dimension(NTERM)     :: NK
  integer(kind=int32)                      :: NP,NPS
  real(kind=real64),dimension(255,NTERM)   :: KARG
  real(kind=real64),dimension(510,NTERM)   :: KERN
  real(kind=real64),dimension(254)  :: WT
  real(kind=real64),dimension(127)  :: WA
  integer(kind=int32),dimension(7)  :: NWT,NWA
  real(kind=real64),dimension(254)  :: FUNCT
  real(kind=real64),dimension(64)   :: FR1,FI1,FR2,FI2,BES1,BES2

  real(kind=real64)                 :: SUM,DIFF,SUMP,SUMM
  real(kind=real64)                 :: FZEROR,FZEROI
  complex(kind=real64)              :: FVAL
  real(kind=real64)                 :: BESF
  integer(kind=int32)               :: LA,LK
  integer(kind=int32)               :: N,KN,K,NA,J,I,NW
  real(kind=real64)                 :: X
  real(kind=real64)                 :: ACUMR,ACUMI


!  starting adresses in arrays WT and WA at face order
     DATA NWT/1,3,7,15,31,63,127/, NWA/1,2,4,8,16,32,64/
!  Array WT contains Gauss quadrature weights stored as needed at
!  each order. Array WA contains the corresponding argument weights.
     DATA (WT(I),I=1,20)/ &
     0.55555555555555555556E+00_real64,0.88888888888888888889E+00_real64, &
     0.26848808986833344073E+00_real64,0.10465622602646726519E+00_real64, &
     0.40139741477596222291E+00_real64,0.45091653865847414235E+00_real64, &
     0.13441525524378422036E+00_real64,0.51603282997079739697E-01_real64, &
     0.20062852937698902103E+00_real64,0.17001719629940260339E-01_real64, &
     0.92927195315124537686E-01_real64,0.17151190913639138079E+00_real64, &
     0.21915685840158749640E+00_real64,0.22551049979820668739E+00_real64, &
     0.67207754295990703540E-01_real64,0.25807598096176653565E-01_real64, &
     0.10031427861179557877E+00_real64,0.84345657393211062463E-02_real64, &
     0.46462893261757986541E-01_real64,0.85755920049990351154E-01_real64/
     DATA (WT(I),I=21,40)/ &
     0.10957842105592463824E+00_real64,0.25447807915618744154E-02_real64, &
     0.16446049854387810934E-01_real64,0.35957103307129322097E-01_real64, &
     0.56979509494123357412E-01_real64,0.76879620499003531043E-01_real64, &
     0.93627109981264473617E-01_real64,0.10566989358023480974E+00_real64, &
     0.11195687302095345688E+00_real64,0.11275525672076869161E+00_real64, &
     0.33603877148207730542E-01_real64,0.12903800100351265626E-01_real64, &
     0.50157139305899537414E-01_real64,0.42176304415588548391E-02_real64, &
     0.23231446639910269443E-01_real64,0.42877960025007734493E-01_real64, &
     0.54789210527962865032E-01_real64,0.12651565562300680114E-02_real64, &
     0.82230079572359296693E-02_real64,0.17978551568128270333E-01_real64/
     DATA (WT(I),I=41,60)/ &
     0.28489754745833548613E-01_real64,0.38439810249455532039E-01_real64, &
     0.46813554990628012403E-01_real64,0.52834946790116519862E-01_real64, &
     0.55978436510476319408E-01_real64,0.36322148184553065969E-03_real64, &
     0.25790497946856882724E-02_real64,0.61155068221172463397E-02_real64, &
     0.10498246909621321898E-01_real64,0.15406750466559497802E-01_real64, &
     0.20594233915912711149E-01_real64,0.25869679327214746911E-01_real64, &
     0.31073551111687964880E-01_real64,0.36064432780782572640E-01_real64, &
     0.40715510116944318934E-01_real64,0.44914531653632197414E-01_real64, &
     0.48564330406673198716E-01_real64,0.51583253952048458777E-01_real64, &
     0.53905499335266063927E-01_real64,0.55481404356559363988E-01_real64/
     DATA (WT(I),I=61,80)/ &
     0.56277699831254301273E-01_real64,0.56377628360384717388E-01_real64, &
     0.16801938574103865271E-01_real64,0.64519000501757369228E-02_real64, &
     0.25078569652949768707E-01_real64,0.21088152457266328793E-02_real64, &
     0.11615723319955134727E-01_real64,0.21438980012503867246E-01_real64, &
     0.27394605263981432516E-01_real64,0.63260731936263354422E-03_real64, &
     0.41115039786546930472E-02_real64,0.89892757840641357233E-02_real64, &
     0.14244877372916774306E-01_real64,0.19219905124727766019E-01_real64, &
     0.23406777495314006201E-01_real64,0.26417473395058259931E-01_real64, &
     0.27989218255238159704E-01_real64,0.18073956444538835782E-03_real64, &
     0.12895240826104173921E-02_real64,0.30577534101755311361E-02_real64/
     DATA (WT(I),I=81,100)/ &
     0.52491234548088591251E-02_real64,0.77033752332797418482E-02_real64, &
     0.10297116957956355524E-01_real64,0.12934839663607373455E-01_real64, &
     0.15536775555843982440E-01_real64,0.18032216390391286320E-01_real64, &
     0.20357755058472159467E-01_real64,0.22457265826816098707E-01_real64, &
     0.24282165203336599358E-01_real64,0.25791626976024229388E-01_real64, &
     0.26952749667633031963E-01_real64,0.27740702178279681994E-01_real64, &
     0.28138849915627150636E-01_real64,0.50536095207862517625E-04_real64, &
     0.37774664632698466027E-03_real64,0.93836984854238150079E-03_real64, &
     0.16811428654214699063E-02_real64,0.25687649437940203731E-02_real64, &
     0.35728927835172996494E-02_real64,0.46710503721143217474E-02_real64/
     DATA (WT(I),I=101,120)/ &
     0.58434498758356395076E-02_real64,0.70724899954335554680E-02_real64, &
     0.83428387539681577056E-02_real64,0.96411777297025366953E-02_real64, &
     0.10955733387837901648E-01_real64,0.12275830560082770087E-01_real64, &
     0.13591571009765546790E-01_real64,0.14893641664815182035E-01_real64, &
     0.16173218729577719942E-01_real64,0.17421930159464173747E-01_real64, &
     0.18631848256138790186E-01_real64,0.19795495048097499488E-01_real64, &
     0.20905851445812023852E-01_real64,0.21956366305317824939E-01_real64, &
     0.22940964229387748761E-01_real64,0.23854052106038540080E-01_real64, &
     0.24690524744487676909E-01_real64,0.25445769965464765813E-01_real64, &
     0.26115673376706097680E-01_real64,0.26696622927450359906E-01_real64/
     DATA (WT(I),I=121,140)/ &
     0.27185513229624791819E-01_real64,0.27579749566481873035E-01_real64, &
     0.27877251476613701609E-01_real64,0.28076455793817246607E-01_real64, &
     0.28176319033016602131E-01_real64,0.28188814180192358694E-01_real64, &
     0.84009692870519326354E-02_real64,0.32259500250878684614E-02_real64, &
     0.12539284826474884353E-01_real64,0.10544076228633167722E-02_real64, &
     0.58078616599775673635E-02_real64,0.10719490006251933623E-01_real64, &
     0.13697302631990716258E-01_real64,0.31630366082226447689E-03_real64, &
     0.20557519893273465236E-02_real64,0.44946378920320678616E-02_real64, &
     0.71224386864583871532E-02_real64,0.96099525623638830097E-02_real64, &
     0.11703388747657003101E-01_real64,0.13208736697529129966E-01_real64/
     DATA (WT(I),I=141,160)/ &
     0.13994609127619079852E-01_real64,0.90372734658751149261E-04_real64, &
     0.64476204130572477933E-03_real64,0.15288767050877655684E-02_real64, &
     0.26245617274044295626E-02_real64,0.38516876166398709241E-02_real64, &
     0.51485584789781777618E-02_real64,0.64674198318036867274E-02_real64, &
     0.77683877779219912200E-02_real64,0.90161081951956431600E-02_real64, &
     0.10178877529236079733E-01_real64,0.11228632913408049354E-01_real64, &
     0.12141082601668299679E-01_real64,0.12895813488012114694E-01_real64, &
     0.13476374833816515982E-01_real64,0.13870351089139840997E-01_real64, &
     0.14069424957813575318E-01_real64,0.25157870384280661489E-04_real64, &
     0.18887326450650491366E-03_real64,0.46918492424785040975E-03_real64/
     DATA (WT(I),I=161,180)/ &
     0.84057143271072246365E-03_real64,0.12843824718970101768E-02_real64, &
     0.17864463917586498247E-02_real64,0.23355251860571608737E-02_real64, &
     0.29217249379178197538E-02_real64,0.35362449977167777340E-02_real64, &
     0.41714193769840788528E-02_real64,0.48205888648512683476E-02_real64, &
     0.54778666939189508240E-02_real64,0.61379152800413850435E-02_real64, &
     0.67957855048827733948E-02_real64,0.74468208324075910174E-02_real64, &
     0.80866093647888599710E-02_real64,0.87109650797320868736E-02_real64, &
     0.93159241280693950932E-02_real64,0.98977475240487497440E-02_real64, &
     0.10452925722906011926E-01_real64,0.10978183152658912470E-01_real64, &
     0.11470482114693874380E-01_real64,0.11927026053019270040E-01_real64/
     DATA (WT(I),I=181,200)/ &
     0.12345262372243838455E-01_real64,0.12722884982732382906E-01_real64, &
     0.13057836688353048840E-01_real64,0.13348311463725179953E-01_real64, &
     0.13592756614812395910E-01_real64,0.13789874783240936517E-01_real64, &
     0.13938625738306850804E-01_real64,0.14038227896908623303E-01_real64, &
     0.14088159516508301065E-01_real64,0.69379364324108267170E-05_real64, &
     0.53275293669780613125E-04_real64,0.13575491094922871973E-03_real64, &
     0.24921240048299729402E-03_real64,0.38974528447328229322E-03_real64, &
     0.55429531493037471492E-03_real64,0.74028280424450333046E-03_real64, &
     0.94536151685852538246E-03_real64,0.11674841174299594077E-02_real64, &
     0.14049079956551446427E-02_real64,0.16561127281544526052E-02_real64/
     DATA (WT(I),I=201,220)/ &
     0.19197129710138724125E-02_real64,0.21944069253638388388E-02_real64, &
     0.24789582266575679307E-02_real64,0.27721957645934509940E-02_real64, &
     0.30730184347025783234E-02_real64,0.33803979910869203823E-02_real64, &
     0.36933779170256508183E-02_real64,0.40110687240750233989E-02_real64, &
     0.43326409680929828545E-02_real64,0.46573172997568547773E-02_real64, &
     0.49843645647655386012E-02_real64,0.53130866051870565663E-02_real64, &
     0.56428181013844441585E-02_real64,0.59729195655081658049E-02_real64, &
     0.63027734490857587172E-02_real64,0.66317812429018878941E-02_real64, &
     0.69593614093904229394E-02_real64,0.72849479805538070639E-02_real64, &
     0.76079896657190565832E-02_real64,0.79279493342948491103E-02_real64/
     DATA (WT(I),I=221,240)/ &
     0.82443037630328680306E-02_real64,0.85565435613076896192E-02_real64, &
     0.88641732094824942641E-02_real64,0.91667111635607884067E-02_real64, &
     0.94636899938300652943E-02_real64,0.97546565363174114611E-02_real64, &
     0.10039172044056840798E-01_real64,0.10316812330947621682E-01_real64, &
     0.10587167904885197931E-01_real64,0.10849844089337314099E-01_real64, &
     0.11104461134006926537E-01_real64,0.11350654315980596602E-01_real64, &
     0.11588074033043952568E-01_real64,0.11816385890830235763E-01_real64, &
     0.12035270785279562630E-01_real64,0.12244424981611985899E-01_real64, &
     0.12443560190714035263E-01_real64,0.12632403643542078765E-01_real64, &
     0.12810698163877361967E-01_real64,0.12978202239537399286E-01_real64/
     DATA (WT(I),I=241,254)/ &
     0.13134690091960152836E-01_real64,0.13279951743930530650E-01_real64, &
     0.13413793085110098513E-01_real64,0.13536035934956213614E-01_real64, &
     0.13646518102571291428E-01_real64,0.13745093443001896632E-01_real64, &
     0.13831631909506428676E-01_real64,0.13906019601325461264E-01_real64, &
     0.13968158806516938516E-01_real64,0.14017968039456608810E-01_real64, &
     0.14055382072649964277E-01_real64,0.14080351962553661325E-01_real64, &
     0.14092845069160408355E-01_real64,0.14094407090096179347E-01_real64/

     DATA (WA(I),I=1,20)/ &
     0.77459666924148337704E+00_real64,0.96049126870802028342E+00_real64, &
     0.43424374934680255800E+00_real64,0.99383196321275502221E+00_real64, &
     0.88845923287225699889E+00_real64,0.62110294673722640294E+00_real64, &
     0.22338668642896688163E+00_real64,0.99909812496766759766E+00_real64, &
     0.98153114955374010687E+00_real64,0.92965485742974005667E+00_real64, &
     0.83672593816886873550E+00_real64,0.70249620649152707861E+00_real64, &
     0.53131974364437562397E+00_real64,0.33113539325797683309E+00_real64, &
     0.11248894313318662575E+00_real64,0.99987288812035761194E+00_real64, &
     0.99720625937222195908E+00_real64,0.98868475754742947994E+00_real64, &
     0.97218287474858179658E+00_real64,0.94634285837340290515E+00_real64/
     DATA (WA(I),I=21,40)/ &
     0.91037115695700429250E+00_real64,0.86390793819369047715E+00_real64, &
     0.80694053195021761186E+00_real64,0.73975604435269475868E+00_real64, &
     0.66290966002478059546E+00_real64,0.57719571005204581484E+00_real64, &
     0.48361802694584102756E+00_real64,0.38335932419873034692E+00_real64, &
     0.27774982202182431507E+00_real64,0.16823525155220746498E+00_real64, &
     0.56344313046592789972E-01_real64,0.99998243035489159858E+00_real64, &
     0.99959879967191068325E+00_real64,0.99831663531840739253E+00_real64, &
     0.99572410469840718851E+00_real64,0.99149572117810613240E+00_real64, &
     0.98537149959852037111E+00_real64,0.97714151463970571416E+00_real64, &
     0.96663785155841656709E+00_real64,0.95373000642576113641E+00_real64/
     DATA (WA(I),I=41,60)/ &
     0.93832039777959288365E+00_real64,0.92034002547001242073E+00_real64, &
     0.89974489977694003664E+00_real64,0.87651341448470526974E+00_real64, &
     0.85064449476835027976E+00_real64,0.82215625436498040737E+00_real64, &
     0.79108493379984836143E+00_real64,0.75748396638051363793E+00_real64, &
     0.72142308537009891548E+00_real64,0.68298743109107922809E+00_real64, &
     0.64227664250975951377E+00_real64,0.59940393024224289297E+00_real64, &
     0.55449513263193254887E+00_real64,0.50768775753371660215E+00_real64, &
     0.45913001198983233287E+00_real64,0.40897982122988867241E+00_real64, &
     0.35740383783153215238E+00_real64,0.30457644155671404334E+00_real64, &
     0.25067873030348317661E+00_real64,0.19589750271110015392E+00_real64/
     DATA (WA(I),I=61,80)/ &
     0.14042423315256017459E+00_real64,0.84454040083710883710E-01_real64, &
     0.28184648949745694339E-01_real64,0.99999759637974846462E+00_real64, &
     0.99994399620705437576E+00_real64,0.99976049092443204733E+00_real64, &
     0.99938033802502358193E+00_real64,0.99874561446809511470E+00_real64, &
     0.99780535449595727456E+00_real64,0.99651414591489027385E+00_real64, &
     0.99483150280062100052E+00_real64,0.99272134428278861533E+00_real64, &
     0.99015137040077015918E+00_real64,0.98709252795403406719E+00_real64, &
     0.98351865757863272876E+00_real64,0.97940628167086268381E+00_real64, &
     0.97473445975240266776E+00_real64,0.96948465950245923177E+00_real64, &
     0.96364062156981213252E+00_real64,0.95718821610986096274E+00_real64/
     DATA (WA(I),I=81,100)/ &
     0.95011529752129487656E+00_real64,0.94241156519108305981E+00_real64, &
     0.93406843615772578800E+00_real64,0.92507893290707565236E+00_real64, &
     0.91543758715576504064E+00_real64,0.90514035881326159519E+00_real64, &
     0.89418456833555902286E+00_real64,0.88256884024734190684E+00_real64, &
     0.87029305554811390585E+00_real64,0.85735831088623215653E+00_real64, &
     0.84376688267270860104E+00_real64,0.82952219463740140018E+00_real64, &
     0.81462878765513741344E+00_real64,0.79909229096084140180E+00_real64, &
     0.78291939411828301639E+00_real64,0.76611781930376009072E+00_real64, &
     0.74869629361693660282E+00_real64,0.73066452124218126133E+00_real64, &
     0.71203315536225203459E+00_real64,0.69281376977911470289E+00_real64/
     DATA (WA(I),I=101,120)/ &
     0.67301883023041847920E+00_real64,0.65266166541001749610E+00_real64, &
     0.63175643771119423041E+00_real64,0.61031811371518640016E+00_real64, &
     0.58836243444766254143E+00_real64,0.56590588542365442262E+00_real64, &
     0.54296566649831149049E+00_real64,0.51955966153745702199E+00_real64, &
     0.49570640791876146017E+00_real64,0.47142506587165887693E+00_real64, &
     0.44673538766202847374E+00_real64,0.42165768662616330006E+00_real64, &
     0.39621280605761593918E+00_real64,0.37042208795007823014E+00_real64, &
     0.34430734159943802278E+00_real64,0.31789081206847668318E+00_real64, &
     0.29119514851824668196E+00_real64,0.26424337241092676194E+00_real64, &
     0.23705884558982972721E+00_real64,0.20966523824318119477E+00_real64/
     DATA (WA(I),I=121,127)/ &
     0.18208649675925219825E+00_real64,0.15434681148137810869E+00_real64, &
     0.12647058437230196685E+00_real64,0.98482396598119202090E-01_real64, &
     0.70406976042855179063E-01_real64,0.42269164765363603212E-01_real64, &
     0.14093886410782462614E-01_real64/

!  check for trivial case
     IF (A.GE.B) THEN
       BESR=0.
       BESI=0.
       RETURN
     ENDIF
     IF (NEW.EQ.2) GOTO 200
!  scale factors
     SUM=(B+A)/2.
     DIFF=(B-A)/2.
!  one point Gauss
     FVAL = F(SUM)
     FZEROR = real(FVAL)
     FZEROI = aimag(FVAL)
     NF=NF+1
     BESF=JBESS(SUM*R,ORDER)
     FZEROR=BESF*FZEROR
     FZEROI=BESF*FZEROI
     IF (NEW.EQ.1) THEN
       KARG(1,NP)=SUM
       KERN(1,NP)=FZEROR
       KERN(2,NP)=FZEROI
       LA=2
       LK=3
     ENDIF
     N=1
     KN=1
5    CONTINUE
!  step through Gauss orders
     DO 40 K=KN,NG
!  compute new function values
       NA=NWA(K)
       DO 10 J=1,NWA(K)
         X=WA(NA)*DIFF
         NA=NA+1
         SUMP=SUM+X
         SUMM=SUM-X
         FVAL = F(SUMP)
         FR1(J) = real(FVAL)
         FI1(J) = aimag(FVAL)
         FVAL = F(SUMM)
         FR2(J) = real(FVAL)
         FI2(J) = aimag(FVAL)
         NF=NF+2
         BES1(J)=JBESS(SUMP*R,ORDER)
         BES2(J)=JBESS(SUMM*R,ORDER)
         IF (NEW.GE.1) THEN
           KARG(LA,NP)=SUMP
           KARG(LA+1,NP)=SUMM
           LA=LA+2
         ENDIF
10     CONTINUE

!  compute products of kernels and Bessel functions
       DO 20 J=1,NWA(K)
         FR1(J)=BES1(J)*FR1(J)
         FI1(J)=BES1(J)*FI1(J)
         FR2(J)=BES2(J)*FR2(J)
         FI2(J)=BES2(J)*FI2(J)
         FUNCT(N)=FR1(J)+FR2(J)
         FUNCT(N+1)=FI1(J)+FI2(J)
         N=N+2
20     CONTINUE
       IF (NEW.GE.1) THEN
         DO 30 J=1,NWA(K)
           KERN(LK,NP)=FR1(J)
           KERN(LK+1,NP)=FI1(J)
           KERN(LK+2,NP)=FR2(J)
           KERN(LK+3,NP)=FI2(J)
           LK=LK+4
30       CONTINUE
       ENDIF
40   CONTINUE

!  compute dot product of weights with integrand values
     NW=NWT(NG)
!  dot should be replaced with SDOT from SCILIB on CRAY-1 if D.P. is
!    not used
     ACUMR=DOT(NW,WT(NW:2*NW-1),1,FUNCT(1:253),2)
     ACUMI=DOT(NW,WT(NW:2*NW-1),1,FUNCT(2:254),2)
     BESR=(ACUMR+WT(2*NW)*FZEROR)*DIFF
     BESI=(ACUMI+WT(2*NW)*FZEROI)*DIFF
     IF (NP.LE.NTERM) NK(NP)=NG
     RETURN
200  CONTINUE
!  construct funct from saved kernels
     DIFF=(B-A)/2.
     FZEROR=KERN(1,NP)
     FZEROI=KERN(2,NP)
     K=MIN(NK(NP),NG)
     NW=NWT(K)
     LK=3
     DO 210 N=1,2*NW,2
       FUNCT(N)=KERN(LK,NP)+KERN(LK+2,NP)
       FUNCT(N+1)=KERN(LK+1,NP)+KERN(LK+3,NP)
       LK=LK+4
210     CONTINUE
     IF (NK(NP).GE.NG) THEN
!  no additional orders required - compute dot product of weitghs with
!    integrand values and return
       ACUMR=DOT(NW,WT(NW:2*NW-1),1,FUNCT(1:253),2)
       ACUMI=DOT(NW,WT(NW:2*NW-1),1,FUNCT(2:254),2)
       BESR=(ACUMR+WT(2*NW)*FZEROR)*DIFF
       BESI=(ACUMI+WT(2*NW)*FZEROI)*DIFF
       RETURN
     ELSE
!  compute additional orders before tkaing dot product
       SUM=(B+A)/2.
       KN=K+1
       LA=LK/2+1
       GOTO 5
     ENDIF

endsubroutine BESQUD



!----------------------------------------------------------------------
!  Computes SUM(S(I)),I=1,...N by computation of Pade approximant
!  using continued fraction expansion. Function is designed to be
!  called sequentially as N is incremented from 1 to its final value.
!  The Nth continued fraction coefficient is calculated and
!  stored and the Nth convergent returned. It is up to the user to
!  stop the calculation when the desired accuracy is achieved.
!  Algorithm from Hänggi et al., Z. Naturforsch. 33A, 402-417 (1977).
!  In their notation, vectors CFCOR,CFCOI are lower case D,vectors DR,
!  DI are upper case D, vectors XR,XI are X, and vectors SR,SI are S.
!
!  A. Chave IGPP/UCSD
!
!  Converted to F90: Rita Streich 2008
!
!----------------------------------------------------------------------
SUBROUTINE PADECF(SUMR,SUMI,SR,SI,N)

  implicit none

  !external variables
  real(kind=real64)                 :: SUMR,SUMI
  real(kind=real64),dimension(1:)   :: SR,SI
  integer(kind=int32)               :: N

  !internal variables
  integer(kind=int32),parameter     :: NC=nmaxsum
  real(kind=real64),dimension(NC)   :: XR,XI
  real(kind=real64),dimension(0:NC) :: DR,DI
  real(kind=real64),dimension(NC)   :: CFCOR,CFCOI
  integer(kind=int32)               :: I,K
  integer(kind=int32)               :: L
  real(kind=real64)                 :: DENOM
  real(kind=real64)                 :: T1,T2

  DATA DR(0)/-1./,DI(0)/0./
  SAVE XR,XI,DR,DI,CFCOR,CFCOI

     IF (N.LT.3) THEN
!  initialize for recursive calculations
       IF (N.EQ.1) THEN
         DO 10 I=1,NC
           XR(I)=0.
           XI(I)=0.
10       CONTINUE
       ENDIF
       DR(N)=SR(N)
       DI(N)=SI(N)
       DENOM=DR(N-1)**2+DI(N-1)**2
       CFCOR(N)=-(DR(N-1)*DR(N)+DI(N-1)*DI(N))/DENOM
       CFCOI(N)=-(DR(N-1)*DI(N)-DR(N)*DI(N-1))/DENOM
       CALL CF(SUMR,SUMI,CFCOR,CFCOI,N)
       RETURN
     ELSE
       L=2*INT((N-1)/2)
!  update X vectors for recursive calculation of coefficients
       DO 20 K=L,4,-2
         XR(K)=XR(K-1)+CFCOR(N-1)*XR(K-2)-CFCOI(N-1)*XI(K-2)
         XI(K)=XI(K-1)+CFCOR(N-1)*XI(K-2)+CFCOI(N-1)*XR(K-2)
20     CONTINUE
       XR(2)=XR(1)+CFCOR(N-1)
       XI(2)=XI(1)+CFCOI(N-1)
!  interchange odd and even parts
       DO 30 K=1,MAX(1,L-1),2
         T1=XR(K)
         T2=XI(K)
         XR(K)=XR(K+1)
         XI(K)=XI(K+1)
         XR(K+1)=T1
         XI(K+1)=T2
30     CONTINUE
!  compute first coefficients
       DR(N)=SR(N)
       DI(N)=SI(N)
       DO 40 K=1,MAX(1,L/2)
         DR(N)=DR(N)+SR(N-K)*XR(2*K-1)-SI(N-K)*XI(2*K-1)
         DI(N)=DI(N)+SI(N-K)*XR(2*K-1)+SR(N-K)*XI(2*K-1)
40     CONTINUE
!  compute new CF coefficient
       DENOM=DR(N-1)**2+DI(N-1)**2
       CFCOR(N)=-(DR(N)*DR(N-1)+DI(N)*DI(N-1))/DENOM
       CFCOI(N)=-(DR(N-1)*DI(N)-DR(N)*DI(N-1))/DENOM
!  evaluate continued fraction
       CALL CF(SUMR,SUMI,CFCOR,CFCOI,N)
       RETURN
     ENDIF

endsubroutine PADECF



!----------------------------------------------------------------------
!  Evaluates a complex continued fraction by recursive division
!  starting at the bottom, as used by PADECF
!  RESR,RESI are real and imaginary parts returned
!  CFCOR,CFCOI are real and imaginary vectors of continued fraction
!  coefficients
!----------------------------------------------------------------------
SUBROUTINE CF(RESR,RESI,CFCOR,CFCOI,N)

  implicit none

  !external variables
  real(kind=real64)                :: RESR,RESI
  real(kind=real64),dimension(1:)  :: CFCOR,CFCOI
  integer(kind=int32)              :: N

  !internal variables
  integer(kind=int32),parameter   :: NC=nmaxsum, NCM=NC-1
  real(kind=real64),dimension(NC) :: ONE
  integer(kind=int32)             :: K
  real(kind=real64)               :: DENOM
  real(kind=real64)               :: RESRO

  DATA ONE/0.,NCM*1./

     RESR=ONE(N)+CFCOR(N)
     RESI=CFCOI(N)
     DO 10 K=N-1,1,-1
       DENOM=RESR**2+RESI**2
       RESRO=RESR
       RESR=ONE(K)+(RESR*CFCOR(K)+RESI*CFCOI(K))/DENOM
       RESI=(RESRO*CFCOI(K)-RESI*CFCOR(K))/DENOM
10   CONTINUE
     RETURN

endsubroutine CF



!----------------------------------------------------------------------
!  Compute dot product of two double precision vectors with nonunit
!  increment allowed. Replacement for BLAS subroutine SDOT.
!----------------------------------------------------------------------
real(kind=real64) FUNCTION DOT(N,X1,INC1,X2,INC2)

  implicit none

  !external variables
  integer(kind=int32)             :: N
  real(kind=real64),dimension(1:) :: X1,X2
  integer(kind=int32)             :: INC1,INC2

  !internal variables
  integer(kind=int32)             :: I,K


     IF (INC2.GT.0) THEN
       K=1
     ELSE
       K=N*ABS(INC2)
     ENDIF
     DOT=0.0
     IF (INC1.GT.0) THEN
       DO 10 I=1,N,INC1
         DOT=DOT+X1(I)*X2(K)
         K=K+INC2
10       CONTINUE
     ELSE
       DO 20 I=N,1,INC1
         DOT=DOT+X1(I)*X2(K)
         K=K+INC2
20       CONTINUE
     ENDIF
     RETURN

endfunction DOT


