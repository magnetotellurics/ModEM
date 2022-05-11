! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."


real(kind=real64) function bessj0(x) result (bj0)

!    =======================================================
!    Purpose: Compute Bessel functions J0(x)
!    Input :  x   --- Argument of J0(x)
!    Output:  BJ0 --- J0(x)
!    =======================================================

  REAL(kind=real64), INTENT(IN)   :: x

  REAL(kind=real64)    :: r, t1, x2
  INTEGER(kind=int32)  :: k

  !internal variables
  real(kind=real64)   :: y,z
  !polynomial coeff. for numerator and denominator for P0(x)
  real(kind=real64)   :: pn0,pn1,pn2,pn3,pn4,pn5,pn6,pn7,pn8,pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8
  !polynomial coeff. for numerator and denominator for Q0(x)
  real(kind=real64)   :: qn0,qn1,qn2,qn3,qn4,qn5,qn6,qn7,qn8,qd0,qd1,qd2,qd3,qd4,qd5,qd6,qd7,qd8,qd9

  real(kind=real64)   :: p0num,p0denom,q0num,q0denom


  data pn0,pn1,pn2,pn3,pn4,pn5,pn6,pn7,pn8 &
       /0.9999999999999999999999995647e+0_real64, 0.5638253933310769952531889297e+1_real64, &
        0.1124846237418285394887270013e+2_real64, 0.1009280644639441488899111404e+2_real64, &
        0.4290591487686900980651458361e+1_real64, 0.8374209971661497198619102718e+0_real64, &
        0.6702347074465611456598882534e-1_real64, 0.1696260729396856143084502774e-2_real64, &
        0.6463970103128382090713889584e-5_real64/
  data pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8 &
       /0.9999999999999999999999999999e+0_real64, 0.5639352566123269952531467562e+1_real64, &
        0.1125463057106955935416066535e+2_real64, 0.1010501892629524191262518048e+2_real64, &
        0.4301396985171094350444425443e+1_real64, 0.8418926086780046799127094223e+0_real64, &
        0.6784915305473610998681570734e-1_real64, 0.1754416614608056207958880988e-2_real64, &
        0.7482977995134121064747276923e-5_real64/

  data qn0,qn1,qn2,qn3,qn4,qn5,qn6,qn7,qn8 &
       /-0.1562499999999999999999995808e-1_real64, -0.1111285583113679178917024959e+0_real64, &
        -0.2877685516355036842789761274e+0_real64, -0.3477683453166454475665803194e+0_real64, &
        -0.2093031978191084473537206358e+0_real64, -0.6209520943730206312601003832e-1_real64, &
        -0.8434508346572023650653353729e-2_real64, -0.4414848186188819989871882393e-3_real64, &
        -0.5768946278415631134804064871e-5_real64/
  data qd0,qd1,qd2,qd3,qd4,qd5,qd6,qd7,qd8,qd9 &
       /0.9999999999999999999999999999e+0_real64, 0.7121383005365046745065850254e+1_real64, &
        0.1848194194302368046679068851e+2_real64, 0.224232752243598371299407153e+2_real64, &
        0.1359286169255959339963319677e+2_real64, 0.408948926810120478008094478e+1_real64, &
        0.5722140925672174525430730669e+0_real64, 0.3219814230905924725810683346e-1_real64, &
        0.5299687475496044642364124073e-3_real64, 0.9423249021001925212258428217e-6_real64/



  x2 = x * x
  IF (x == 0.0_real64) THEN
    bj0 = 1.0_real64
    RETURN
  END IF
  IF (x <= 10.0_real64) THEN
    bj0 = 1.0_real64
    r = 1.0_real64
    DO  k = 1, 35
      r = -0.25_real64 * r * x2 / (k*k)
      bj0 = bj0 + r
      IF (ABS(r) < ABS(bj0)*relerrbes) EXIT
    END DO
  ELSE

    z=8./x  
    y=z**2  
    t1 = x - 0.25_real64 * dpi

    p0num = pn0+y*(pn1+y*(pn2+y*(pn3+y*(pn4+y*(pn5+y*(pn6+y*(pn7+y*pn8)))))))
    p0denom = pd0+y*(pd1+y*(pd2+y*(pd3+y*(pd4+y*(pd5+y*(pd6+y*(pd7+y*pd8)))))))

    q0num = qn0+y*(qn1+y*(qn2+y*(qn3+y*(qn4+y*(qn5+y*(qn6+y*(qn7+y*qn8)))))))
    q0denom = qd0+y*(qd1+y*(qd2+y*(qd3+y*(qd4+y*(qd5+y*(qd6+y*(qd7+y*(qd8+y*qd9))))))))

    bj0=sqrt(rp2/x)*(cos(t1)*(p0num/p0denom)-z*sin(t1)*(q0num/q0denom))  

  END IF

endfunction bessj0
 


real(kind=real64) function bessj1(x) result(bj1)

!    =======================================================
!    Purpose: Compute Bessel function J1(x)
!    Input :  x   --- Argument of J1(x)
!    Output:  BJ1 --- J1(x)
!    =======================================================

  REAL(kind=real64), INTENT(IN)   :: x

  REAL(kind=real64)    :: r, t2, x2
  INTEGER(kind=int32)  :: k


  !internal variables
  real(kind=real64)   :: y,z
  !polynomial coeff. for numerator and denominator for P0(x)
  real(kind=real64)   :: pn0,pn1,pn2,pn3,pn4,pn5,pn6,pn7,pn8,pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8
  !polynomial coeff. for numerator and denominator for Q0(x)
  real(kind=real64)   :: qn0,qn1,qn2,qn3,qn4,qn5,qn6,qn7,qn8,qd0,qd1,qd2,qd3,qd4,qd5,qd6,qd7,qd8,qd9

  real(kind=real64)   :: p0num,p0denom,q0num,q0denom

  data pn0,pn1,pn2,pn3,pn4,pn5,pn6,pn7,pn8 &
       /0.1000000000000000000000000489e+1_real64, 0.5581663300347182292169450071e+1_real64, &
        0.1100186625131173123750501118e+2_real64, 0.9727139359130463694593683431e+1_real64, &
        0.4060011483142278994462590992e+1_real64, 0.7742832212665311906917358099e+0_real64, &
        0.602161775281109875209824863e-1_real64, 0.1482350677236405118074646993e-2_real64, &
        0.6094215148131061431667573909e-5_real64/
  data pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8 &
       /0.9999999999999999999999999999e+0_real64, 0.5579832245659682292169922224e+1_real64, &
        0.109916844773161728897277104e+2_real64, 0.9707206835125961446797916892e+1_real64, &
        0.4042610016540342097334497865e+1_real64, 0.7671965204303836019508430169e+0_real64, &
        0.5893258668794493100786371406e-1_real64, 0.139399364498125685240422253e-2_real64, &
        0.4585597769784750669754696825e-5_real64/

  data qn0,qn1,qn2,qn3,qn4,qn5,qn6,qn7,qn8 &
       /0.4687499999999999999999995275e-1_real64, 0.3302394516691663879252493748e+0_real64, &
        0.8456888491208195767613862428e+0_real64, 0.1008551084218946085420665147e+1_real64, &
        0.5973407972399900690521296181e+0_real64, 0.1737697433393258207540273097e+0_real64, &
        0.230386281481956857389361074e-1_real64, 0.1171224207976250587945594946e-2_real64, &
        0.1486418220337492918307904804e-4_real64/
  data qd0,qd1,qd2,qd3,qd4,qd5,qd6,qd7,qd8,qd9 &
       /0.9999999999999999999999999999e+0_real64, 0.7049380763213049609070823421e+1_real64, &
        0.1807129960468949760845562209e+2_real64, 0.2159171174362827330505421695e+2_real64, &
        0.1283239297740546866114600499e+2_real64, 0.3758349275324260869598403931e+1_real64, &
        0.5055985453754739528620657666e+0_real64, 0.2665604326323907148063400439e-1_real64, &
        0.3821140353404633025596424652e-3_real64, 0.3206696590241261037875154062e-6_real64/



  x2 = x * x
  IF (x == 0.0_real64) THEN
    bj1 = 0.0_real64
    RETURN
  END IF
  IF (x <= 9.0_real64) THEN
    bj1 = 1.0_real64
    r = 1.0_real64
    DO  k = 1, 35
      r = -0.25_real64 * r * x2 / (k*(k+1))
      bj1 = bj1 + r
      IF (ABS(r) < ABS(bj1)*relerrbes) EXIT
    END DO
    bj1 = 0.5_real64 * x * bj1
  ELSE

    z=8./x  
    y=z**2  
    t2 = x - 0.75_real64 * dpi

    p0num = pn0+y*(pn1+y*(pn2+y*(pn3+y*(pn4+y*(pn5+y*(pn6+y*(pn7+y*pn8)))))))
    p0denom = pd0+y*(pd1+y*(pd2+y*(pd3+y*(pd4+y*(pd5+y*(pd6+y*(pd7+y*pd8)))))))

    q0num = qn0+y*(qn1+y*(qn2+y*(qn3+y*(qn4+y*(qn5+y*(qn6+y*(qn7+y*qn8)))))))
    q0denom = qd0+y*(qd1+y*(qd2+y*(qd3+y*(qd4+y*(qd5+y*(qd6+y*(qd7+y*(qd8+y*qd9))))))))

    bj1=sqrt(rp2/x)*(cos(t2)*(p0num/p0denom)-z*sin(t2)*(q0num/q0denom))


  END IF

endfunction bessj1



!----------------------------------------------------------------------
!  Computes zero of Bessel function of the first kind from
!  McMahon's asymptotic expansion
!  NZERO - number of the zero
!  ORDER - order of the Bessel function (0 or 1)
!----------------------------------------------------------------------
real(kind=real64) FUNCTION ZEROJ(NZERO,ORDER)

  implicit none

  !external variables
  integer(kind=int32)      :: NZERO
  real(kind=real64)        :: ORDER

  !internal variables
  real(kind=real64),parameter  :: ZT1=-1.D0/8.D0, ZT2=124.D0/1536.D0,ZT3=-120928.D0/491520.D0, &
                       ZT4=401743168.D0/220200960.D0
  real(kind=real64),parameter  :: OT1=3.D0/8.D0,OT2=-36.D0/1536.D0,OT3=113184.D0/491520.D0, &
                       OT4=-1951209.D0/220200960.D0
  real(kind=real64)            :: BETA

     IF (ORDER.EQ.0.0) THEN
       BETA=(real(NZERO,kind=real64)-.25_real64) * dpi
       ZEROJ=BETA-ZT1/BETA-ZT2/BETA**3-ZT3/BETA**5-ZT4/BETA**7
     ELSEIF (ORDER.EQ.1.0) THEN
       BETA=(real(NZERO,kind=real64)+.25_real64) * dpi
       ZEROJ=BETA-OT1/BETA-OT2/BETA**3-OT3/BETA**5-OT4/BETA**7
     ENDIF

endfunction ZEROJ


!----------------------------------------------------------------------
!  Computes Bessel function of order ORDER and argument X by calling
!  SLATEC routines DBESJ0 and DBESJ1
!  May be replaced with IMSL routines MMBSJ0,MMBSJ1
!  RS: here replaced by Numerical Recipes routines bessj0 and bessj1
!----------------------------------------------------------------------
real(kind=real64) FUNCTION JBESS(X,ORDER)

  implicit none

  !external variables
  real(kind=real64)    :: X,ORDER

  !internal variables

     IF (ORDER.EQ.0.0) THEN
       JBESS=bessj0(X)
     ELSEIF (ORDER.EQ.1.0) THEN
       JBESS=bessj1(X)
     ENDIF
     RETURN

endfunction JBESS

