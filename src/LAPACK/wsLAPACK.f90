!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!module wsLAPACK
!  module containing routines from LAPACK & BLAS
!   use wsBLAS
!   contains

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  ====================================================================

      SUBROUTINE DPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
!
!  -- LAPACK driver routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPPSV computes the solution to a real system of linear equations
!     A * X = B,
!  where A is an N-by-N symmetric positive definite matrix stored in
!  packed format and X and B are N-by-NRHS matrices.
!
!  The Cholesky decomposition is used to factor A as
!     A = U**T* U,  if UPLO = 'U', or
!     A = L * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is a lower triangular
!  matrix.  The factored form of A is then used to solve the system of
!  equations A * X = B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!          See below for further details.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T, in the same storage
!          format as A.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i of A is not
!                positive definite, so the factorization could not be
!                completed, and the solution has not been computed.
!
!  Further Details
!  ===============
!
!  The packed storage scheme is illustrated by the following example
!  when N = 4, UPLO = 'U':
!
!  Two-dimensional storage of the symmetric matrix A:
!
!     a11 a12 a13 a14
!         a22 a23 a24
!             a33 a34     (aij = conjg(aji))
!                 a44
!
!  Packed storage of the upper triangle of A:
!
!  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
!
!  =====================================================================
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DPPTRF, DPPTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPPSV ', -INFO )
         RETURN
      END IF
!
!     Compute the Cholesky factorization A = U'*U or A = L*L'.
!
      CALL DPPTRF( UPLO, N, AP, INFO )
      IF( INFO.EQ.0 ) THEN
!
!        Solve the system A*X = B, overwriting B with X.
!
         CALL DPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
!
      END IF
      RETURN
!
!     End of DPPSV
!
      END SUBROUTINE DPPSV

!  ====================================================================

      SUBROUTINE DPPTRF( UPLO, N, AP, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * )
!     ..
!
!  Purpose
!  =======
!
!  DPPTRF computes the Cholesky factorization of a real symmetric
!  positive definite matrix A stored in packed format.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!          See below for further details.
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**T*U or A = L*L**T, in the same
!          storage format as A.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ======= =======
!
!  The packed storage scheme is illustrated by the following example
!  when N = 4, UPLO = 'U':
!
!  Two-dimensional storage of the symmetric matrix A:
!
!     a11 a12 a13 a14
!         a22 a23 a24
!             a33 a34     (aij = aji)
!                 a44
!
!  Packed storage of the upper triangle of A:
!
!  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JC, JJ
      DOUBLE PRECISION   AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSPR, DTPSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPPTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         JJ = 0
         DO 10 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
!
!           Compute elements 1:J-1 of column J.
!
            IF( J.GT.1 ) &
               CALL DTPSV( 'Upper', 'Transpose', 'Non-unit', J-1, AP, &
                           AP( JC ), 1 )
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = AP( JJ ) - DDOT( J-1, AP( JC ), 1, AP( JC ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               AP( JJ ) = AJJ
               GO TO 30
            END IF
            AP( JJ ) = SQRT( AJJ )
   10    CONTINUE
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         JJ = 1
         DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = AP( JJ )
            IF( AJJ.LE.ZERO ) THEN
               AP( JJ ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            AP( JJ ) = AJJ
!
!           Compute elements J+1:N of column J and update the trailing
!           submatrix.
!
            IF( J.LT.N ) THEN
               CALL DSCAL( N-J, ONE / AJJ, AP( JJ+1 ), 1 )
               CALL DSPR( 'Lower', N-J, -ONE, AP( JJ+1 ), 1, &
                          AP( JJ+N-J+1 ) )
               JJ = JJ + N - J + 1
            END IF
   20    CONTINUE
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = J
!
   40 CONTINUE
      RETURN
!
!     End of DPPTRF
!
      END SUBROUTINE DPPTRF

!  ====================================================================

      SUBROUTINE DPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPPTRS solves a system of linear equations A*X = B with a symmetric
!  positive definite matrix A in packed storage using the Cholesky
!  factorization A = U**T*U or A = L*L**T computed by DPPTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The triangular factor U or L from the Cholesky factorization
!          A = U**T*U or A = L*L**T, packed columnwise in a linear
!          array.  The j-th column of U or L is stored in the array AP
!          as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DTPSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPPTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Solve A*X = B where A = U'*U.
!
         DO 10 I = 1, NRHS
!
!           Solve U'*X = B, overwriting B with X.
!
            CALL DTPSV( 'Upper', 'Transpose', 'Non-unit', N, AP, &
                        B( 1, I ), 1 )
!
!           Solve U*X = B, overwriting B with X.
!
            CALL DTPSV( 'Upper', 'No transpose', 'Non-unit', N, AP, &
                        B( 1, I ), 1 )
   10    CONTINUE
      ELSE
!
!        Solve A*X = B where A = L*L'.
!
         DO 20 I = 1, NRHS
!
!           Solve L*Y = B, overwriting B with X.
!
            CALL DTPSV( 'Lower', 'No transpose', 'Non-unit', N, AP, &
                        B( 1, I ), 1 )
!
!           Solve L'*X = Y, overwriting B with X.
!
            CALL DTPSV( 'Lower', 'Transpose', 'Non-unit', N, AP, &
                        B( 1, I ), 1 )
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of DPPTRS
!
      END SUBROUTINE DPPTRS

!  ====================================================================

      SUBROUTINE DPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
!
!  -- LAPACK driver routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPBSV computes the solution to a real system of linear equations
!     A * X = B,
!  where A is an N-by-N symmetric positive definite band matrix and X
!  and B are N-by-NRHS matrices.
!
!  The Cholesky decomposition is used to factor A as
!     A = U**T * U,  if UPLO = 'U', or
!     A = L * L**T,  if UPLO = 'L',
!  where U is an upper triangular band matrix, and L is a lower
!  triangular band matrix, with the same number of superdiagonals or
!  subdiagonals as A.  The factored form of A is then used to solve the
!  system of equations A * X = B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD).
!          See below for further details.
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**T*U or A = L*L**T of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i of A is not
!                positive definite, so the factorization could not be
!                completed, and the solution has not been computed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked * are not used by the routine.
!
!  =====================================================================
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DPBTRF, DPBTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBSV ', -INFO )
         RETURN
      END IF
!
!     Compute the Cholesky factorization A = U'*U or A = L*L'.
!
      CALL DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
      IF( INFO.EQ.0 ) THEN
!
!        Solve the system A*X = B, overwriting B with X.
!
         CALL DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
!
      END IF
      RETURN
!
!     End of DPBSV
!
      END SUBROUTINE DPBSV

!  ====================================================================

      SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPBTF2 computes the Cholesky factorization of a real symmetric
!  positive definite band matrix A.
!
!  The factorization has the form
!     A = U' * U ,  if UPLO = 'U', or
!     A = L  * L',  if UPLO = 'L',
!  where U is an upper triangular matrix, U' is the transpose of U, and
!  L is lower triangular.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of super-diagonals of the matrix A if UPLO = 'U',
!          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U'*U or A = L*L' of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked * are not used by the routine.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, KLD, KN
      DOUBLE PRECISION   AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSYR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      KLD = MAX( 1, LDAB-1 )
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         DO 10 J = 1, N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = AB( KD+1, J )
            IF( AJJ.LE.ZERO ) &
               GO TO 30
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
!
!           Compute elements J+1:J+KN of row J and update the
!           trailing submatrix within the band.
!
            KN = MIN( KD, N-J )
            IF( KN.GT.0 ) THEN
               CALL DSCAL( KN, ONE / AJJ, AB( KD, J+1 ), KLD )
               CALL DSYR( 'Upper', KN, -ONE, AB( KD, J+1 ), KLD, &
                          AB( KD+1, J+1 ), KLD )
            END IF
   10    CONTINUE
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = AB( 1, J )
            IF( AJJ.LE.ZERO ) &
               GO TO 30
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
!
!           Compute elements J+1:J+KN of column J and update the
!           trailing submatrix within the band.
!
            KN = MIN( KD, N-J )
            IF( KN.GT.0 ) THEN
               CALL DSCAL( KN, ONE / AJJ, AB( 2, J ), 1 )
               CALL DSYR( 'Lower', KN, -ONE, AB( 2, J ), 1, &
                          AB( 1, J+1 ), KLD )
            END IF
   20    CONTINUE
      END IF
      RETURN
!
   30 CONTINUE
      INFO = J
      RETURN
!
!     End of DPBTF2
!
      END SUBROUTINE DPBTF2

!  ====================================================================

      SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPBTRF computes the Cholesky factorization of a real symmetric
!  positive definite band matrix A.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**T*U or A = L*L**T of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked * are not used by the routine.
!
!  Contributed by
!  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 32, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, IB, II, J, JJ, NB
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   WORK( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DPBTF2, DPOTF2, DSYRK, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( ( .NOT.LSAME( UPLO, 'U' ) ) .AND. &
          ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'DPBTRF', UPLO, N, KD, -1, -1 )
!
!     The block size must not exceed the semi-bandwidth KD, and must not
!     exceed the limit set by the size of the local array WORK.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KD ) THEN
!
!        Use unblocked code
!
         CALL DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
      ELSE
!
!        Use blocked code
!
         IF( LSAME( UPLO, 'U' ) ) THEN
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the upper triangle of the matrix in band
!           storage.
!
!           Zero the upper triangle of the work array.
!
            DO 20 J = 1, NB
               DO 10 I = 1, J - 1
                  WORK( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
!
!           Process the band matrix one diagonal block at a time.
!
            DO 70 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL DPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11   A12   A13
!                          A22   A23
!                                A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A12, A22 and
!                 A23 are empty if IB = KD. The upper triangle of A13
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A12
!
                     CALL DTRSM( 'Left', 'Upper', 'Transpose', &
                                 'Non-unit', IB, I2, ONE, AB( KD+1, I ), &
                                 LDAB-1, AB( KD+1-IB, I+IB ), LDAB-1 )
!
!                    Update A22
!
                     CALL DSYRK( 'Upper', 'Transpose', I2, IB, -ONE, &
                                 AB( KD+1-IB, I+IB ), LDAB-1, ONE, &
                                 AB( KD+1, I+IB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Copy the lower triangle of A13 into the work array.
!
                     DO 40 JJ = 1, I3
                        DO 30 II = JJ, IB
                           WORK( II, JJ ) = AB( II-JJ+1, JJ+I+KD-1 )
   30                   CONTINUE
   40                CONTINUE
!
!                    Update A13 (in the work array).
!
                     CALL DTRSM( 'Left', 'Upper', 'Transpose', &
                                 'Non-unit', IB, I3, ONE, AB( KD+1, I ), &
                                 LDAB-1, WORK, LDWORK )
!
!                    Update A23
!
                     IF( I2.GT.0 ) &
                        CALL DGEMM( 'Transpose', 'No Transpose', I2, I3, &
                                    IB, -ONE, AB( KD+1-IB, I+IB ), &
                                    LDAB-1, WORK, LDWORK, ONE, &
                                    AB( 1+IB, I+KD ), LDAB-1 )
!
!                    Update A33
!
                     CALL DSYRK( 'Upper', 'Transpose', I3, IB, -ONE, &
                                 WORK, LDWORK, ONE, AB( KD+1, I+KD ), &
                                 LDAB-1 )
!
!                    Copy the lower triangle of A13 back into place.
!
                     DO 60 JJ = 1, I3
                        DO 50 II = JJ, IB
                           AB( II-JJ+1, JJ+I+KD-1 ) = WORK( II, JJ )
   50                   CONTINUE
   60                CONTINUE
                  END IF
               END IF
   70       CONTINUE
         ELSE
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the lower triangle of the matrix in band
!           storage.
!
!           Zero the lower triangle of the work array.
!
            DO 90 J = 1, NB
               DO 80 I = J + 1, NB
                  WORK( I, J ) = ZERO
   80          CONTINUE
   90       CONTINUE
!
!           Process the band matrix one diagonal block at a time.
!
            DO 140 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL DPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11
!                    A21   A22
!                    A31   A32   A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A21, A22 and
!                 A32 are empty if IB = KD. The lower triangle of A31
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A21
!
                     CALL DTRSM( 'Right', 'Lower', 'Transpose', &
                                 'Non-unit', I2, IB, ONE, AB( 1, I ), &
                                 LDAB-1, AB( 1+IB, I ), LDAB-1 )
!
!                    Update A22
!
                     CALL DSYRK( 'Lower', 'No Transpose', I2, IB, -ONE, &
                                 AB( 1+IB, I ), LDAB-1, ONE, &
                                 AB( 1, I+IB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Copy the upper triangle of A31 into the work array.
!
                     DO 110 JJ = 1, IB
                        DO 100 II = 1, MIN( JJ, I3 )
                           WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
  100                   CONTINUE
  110                CONTINUE
!
!                    Update A31 (in the work array).
!
                     CALL DTRSM( 'Right', 'Lower', 'Transpose', &
                                 'Non-unit', I3, IB, ONE, AB( 1, I ), &
                                 LDAB-1, WORK, LDWORK )
!
!                    Update A32
!
                     IF( I2.GT.0 ) &
                        CALL DGEMM( 'No transpose', 'Transpose', I3, I2, &
                                    IB, -ONE, WORK, LDWORK, &
                                    AB( 1+IB, I ), LDAB-1, ONE, &
                                    AB( 1+KD-IB, I+IB ), LDAB-1 )
!
!                    Update A33
!
                     CALL DSYRK( 'Lower', 'No Transpose', I3, IB, -ONE, &
                                 WORK, LDWORK, ONE, AB( 1, I+KD ), &
                                 LDAB-1 )
!
!                    Copy the upper triangle of A31 back into place.
!
                     DO 130 JJ = 1, IB
                        DO 120 II = 1, MIN( JJ, I3 )
                           AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
  120                   CONTINUE
  130                CONTINUE
                  END IF
               END IF
  140       CONTINUE
         END IF
      END IF
      RETURN
!
  150 CONTINUE
      RETURN
!
!     End of DPBTRF
!
      END SUBROUTINE DPBTRF

!  ====================================================================

      SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPBTRS solves a system of linear equations A*X = B with a symmetric
!  positive definite band matrix A using the Cholesky factorization
!  A = U**T*U or A = L*L**T computed by DPBTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangular factor stored in AB;
!          = 'L':  Lower triangular factor stored in AB.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          The triangular factor U or L from the Cholesky factorization
!          A = U**T*U or A = L*L**T of the band matrix A, stored in the
!          first KD+1 rows of the array.  The j-th column of U or L is
!          stored in the j-th column of the array AB as follows:
!          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DTBSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Solve A*X = B where A = U'*U.
!
         DO 10 J = 1, NRHS
!
!           Solve U'*X = B, overwriting B with X.
!
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
!
!           Solve U*X = B, overwriting B with X.
!
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
   10    CONTINUE
      ELSE
!
!        Solve A*X = B where A = L*L'.
!
         DO 20 J = 1, NRHS
!
!           Solve L*X = B, overwriting B with X.
!
            CALL DTBSV( 'Lower', 'No transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
!
!           Solve L'*X = B, overwriting B with X.
!
            CALL DTBSV( 'Lower', 'Transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of DPBTRS
!
      END SUBROUTINE DPBTRS

!  ====================================================================

      SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DPOTF2 computes the Cholesky factorization of a real symmetric
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U' * U ,  if UPLO = 'U', or
!     A = L  * L',  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n by n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U'*U  or A = L*L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         DO 10 J = 1, N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!
!           Compute elements J+1:N of row J.
!
            IF( J.LT.N ) THEN
               CALL DGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ), &
                           LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = A( J, J ) - DDOT( J-1, A( J, 1 ), LDA, A( J, 1 ), &
                  LDA )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!
!           Compute elements J+1:N of column J.
!
            IF( J.LT.N ) THEN
               CALL DGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ), &
                           LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = J
!
   40 CONTINUE
      RETURN
!
!     End of DPOTF2
!
      END SUBROUTINE DPOTF2

!  ====================================================================

!  =====================================================================

      SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  DGEQRF computes a QR factorization of a real M-by-N matrix A:
!  A = Q * R.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of min(m,n) elementary reflectors (see Further
!          Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is
!          the optimal blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!  and tau in TAU(i).
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEQR2, DLARFB, DLARFT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size.
!
      NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', M, N, -1, &
                       -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code initially
!
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
!
!           Compute the QR factorization of the current block
!           A(i:m,i:i+ib-1)
!
            CALL DGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, &
                         IINFO )
            IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                            A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H' to A(i:m,i+ib:n) from the left
!
               CALL DLARFB( 'Left', 'Transpose', 'Forward', &
                            'Columnwise', M-I+1, N-I-IB+1, IB, &
                            A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                            LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
!
!     Use unblocked code to factor the last or only block.
!
      IF( I.LE.K ) &
         CALL DGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, &
                      IINFO )
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DGEQRF
!
      END SUBROUTINE DGEQRF

!  =====================================================================

      SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEQR2 computes a QR factorization of a real m by n matrix A:
!  A = Q * R.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(m,n) by n upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!  and tau in TAU(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQR2', -INFO )
         RETURN
      END IF
!
      K = MIN( M, N )
!
      DO 10 I = 1, K
!
!        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
         CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, &
                      TAU( I ) )
         IF( I.LT.N ) THEN
!
!           Apply H(i) to A(i:m,i+1:n) from the left
!
            AII = A( I, I )
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                        A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
!
!     End of DGEQR2
!
      END SUBROUTINE DGEQR2

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
!     ..
!
!  Purpose
!  =======
!
!  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.
!
!  Arguments
!  =========
!
!  X       (input) DOUBLE PRECISION
!  Y       (input) DOUBLE PRECISION
!          X and Y specify the values x and y.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
!
!     End of DLAPY2
!
      END FUNCTION DLAPY2

!  =====================================================================

      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) DOUBLE PRECISION array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) DOUBLE PRECISION
!          The value tau in the representation of H.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C
!
         IF( TAU.NE.ZERO ) THEN
!
!           w := C' * v
!
            CALL DGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO, &
                        WORK, 1 )
!
!           C := C - v * w'
!
            CALL DGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( TAU.NE.ZERO ) THEN
!
!           w := C * v
!
            CALL DGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV, &
                        ZERO, WORK, 1 )
!
!           C := C - w * v'
!
            CALL DGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
!
!     End of DLARF
!
      END SUBROUTINE DLARF

!  =====================================================================

      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                         T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ), &
                         WORK( LDWORK, * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFB applies a real block reflector H or its transpose H' to a
!  real m by n matrix C, from either the left or the right.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply H or H' from the Left
!          = 'R': apply H or H' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H' (Transpose)
!
!  DIRECT  (input) CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  K       (input) INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).
!
!  V       (input) DOUBLE PRECISION array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
!  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 ) &
         RETURN
!
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
!
      IF( LSAME( STOREV, 'C' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C1'
!
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2
!
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                              ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,  &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                              -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                              C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                           ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                              ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                              C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                           ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
!
         ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C2'
!
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1
!
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1 * W'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                              -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                           ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                           K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                           ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W
!
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
!
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C1'
!
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                           ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2'
!
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                              C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
                              WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2' * W'
!
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                              V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                              C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C1
!
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                           ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                              ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                              C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
!
            END IF
!
         ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C2'
!
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                           ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1'
!
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                              C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1' * W'
!
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                              V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C2
!
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                           ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
!
            END IF
!
         END IF
      END IF
!
      RETURN
!
!     End of DLARFB
!
      END SUBROUTINE DLARFB

!  =====================================================================

      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H' * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v' ) ,
!                      ( v )
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) DOUBLE PRECISION
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) DOUBLE PRECISION array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) DOUBLE PRECISION
!          The value tau.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = DNRM2( N-1, X, INCX )
!
      IF( XNORM.EQ.ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        general case
!
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN ) &
               GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!
!           If ALPHA is subnormal, it may lose relative accuracy
!
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            ALPHA = BETA
         END IF
      END IF
!
      RETURN
!
!     End of DLARFG
!
      END SUBROUTINE DLARFG

!  =====================================================================

      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFT forms the triangular factor T of a real block reflector H
!  of order n, which is defined as a product of k elementary reflectors.
!
!  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
!  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and
!
!     H  =  I - V * T * V'
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and
!
!     H  =  I - V' * T * V
!
!  Arguments
!  =========
!
!  DIRECT  (input) CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
!  N       (input) INTEGER
!          The order of the block reflector H. N >= 0.
!
!  K       (input) INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.
!
!  V       (input/output) DOUBLE PRECISION array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
!  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  Further Details
!  ===============
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   VII
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
!
!              general case
!
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
!
!                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
!
                  CALL DGEMV( 'Transpose', N-I+1, I-1, -TAU( I ), &
                              V( I, 1 ), LDV, V( I, I ), 1, ZERO, &
                              T( 1, I ), 1 )
               ELSE
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
!
                  CALL DGEMV( 'No transpose', I-1, N-I+1, -TAU( I ), &
                              V( 1, I ), LDV, V( I, I ), LDV, ZERO, &
                              T( 1, I ), 1 )
               END IF
               V( I, I ) = VII
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
                           LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
!
!              general case
!
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
!
                     CALL DGEMV( 'Transpose', N-K+I, K-I, -TAU( I ), &
                                 V( 1, I+1 ), LDV, V( 1, I ), 1, ZERO, &
                                 T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
!
                     CALL DGEMV( 'No transpose', K-I, N-K+I, -TAU( I ), &
                                 V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, &
                                 T( I+1, I ), 1 )
                     V( I, N-K+I ) = VII
                  END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                              T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
!
!     End of DLARFT
!
      END SUBROUTINE DLARFT

!  =====================================================================

      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN, &
                         RND, SFMIN, SMALL, T
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAMC2
!     ..
!     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN, &
                         EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
!
      DLAMCH = RMACH
      RETURN
!
!     End of DLAMCH
!
      END FUNCTION DLAMCH
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!     ..
!
!  Purpose
!  =======
!
!  DLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         B = 1
         C = DLAMC3( A, B )
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) ) &
            LRND = .FALSE.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         LT = 0
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
!+       END WHILE
!
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
!
!     End of DLAMC1
!
      END SUBROUTINE DLAMC1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
!     ..
!
!  Purpose
!  =======
!
!  DLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) DOUBLE PRECISION
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) DOUBLE PRECISION
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, &
                         NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
                         SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX, &
                         LRMIN, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
!
!        Start to find EPS.
!
         B = LBETA
         A = B**( -LT )
         LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS ) &
            B = LEPS
!
         LEPS = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
!+       END WHILE
!
         IF( A.LT.LEPS ) &
            LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
!
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND. &
                  ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
!!!
! Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
!!!
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine DLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         IEEE = IEEE .OR. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
!
!        Finally, call DLAMC5 to compute EMAX and RMAX.
!
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
      RETURN
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-', &
            '  EMIN = ', I8, / &
            ' If, after inspection, the value EMIN looks', &
            ' acceptable please comment out ', &
            / ' the IF block as marked within the code of routine', &
            ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of DLAMC2
!
      END SUBROUTINE DLAMC2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
!     ..
!
!  Purpose
!  =======
!
!  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) DOUBLE PRECISION
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      DLAMC3 = A + B
!
      RETURN
!
!     End of DLAMC3
!
      END FUNCTION DLAMC3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE DLAMC4( EMIN, START, BASE )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
!     ..
!
!  Purpose
!  =======
!
!  DLAMC4 is a service routine for DLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) DOUBLE PRECISION
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Executable Statements ..
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. &
!            ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND. &
          ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of DLAMC4
!
      END SUBROUTINE DLAMC4
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
!     ..
!
!  Purpose
!  =======
!
!  DLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         EMAX = EMAX - 1
      END IF
!
      IF( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         EMAX = EMAX - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE ) &
            OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE ) &
         Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
!
      RMAX = Y
      RETURN
!
!     End of DLAMC5
!
      END SUBROUTINE DLAMC5

!  ====================================================================

      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), &
                         WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  DORMQR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, &
                         MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size.  NB may be at most NBMAX, where NBMAX
!     is used to define the local array T.
!
      NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K, &
           -1 ) )
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K, &
                    -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                      IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. &
             ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), &
                         LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H' is applied to C(i:m,1:n)
!
               MI = M - I + 1
               IC = I
            ELSE
!
!              H or H' is applied to C(1:m,i:n)
!
               NI = N - I + 1
               JC = I
            END IF
!
!           Apply H or H'
!
            CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
                         IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, &
                         WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = IWS
      RETURN
!
!     End of DORMQR
!
      END SUBROUTINE DORMQR

!================================================================

      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORM2R overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q' (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
         RETURN
!
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) &
           THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
         ELSE
!
!           H(i) is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
         END IF
!
!        Apply H(i)
!
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ), &
                     LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of DORM2R
!
      END SUBROUTINE DORM2R
!  ====================================================================

!<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      SUBROUTINE ZGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGBTRF computes an LU factorization of a complex m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ) )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, &
                         JU, K2, KM, KV, NB, NW
      COMPLEX*16         TEMP
!     ..
!     .. Local Arrays ..
      COMPLEX*16         WORK13( LDWORK, NBMAX ), &
                         WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      INTEGER            ILAENV, IZAMAX
      EXTERNAL           ILAENV, IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZGBTF2, ZGEMM, ZGERU, ZLASWP, &
                         ZSCAL, ZSWAP, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'ZGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
!
!        Use unblocked code
!
         CALL ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
!
!        Zero the subdiagonal elements of the work array WORK31
!
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
         JU = 1
!
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
            DO 80 JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
               KM = MIN( KL, M-JJ )
               JP = IZAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
                     IF( JP+JJ-1.LT.J+KL ) THEN
!
                        CALL ZSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1, &
                                    AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                        CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                    WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL ZSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1, &
                                    AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
!
!                 Compute multipliers
!
                  CALL ZSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), &
                              1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ ) &
                     CALL ZGERU( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1, &
                                 AB( KV, JJ+1 ), LDAB-1, &
                                 AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
                  IF( INFO.EQ.0 ) &
                     INFO = JJ
               END IF
!
!              Copy current column of A31 into the work array WORK31
!
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 ) &
                  CALL ZCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1, &
                              WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
!
!              Use ZLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
               CALL ZLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, &
                            IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
!
!              Update the relevant part of the trailing submatrix
!
               IF( J2.GT.0 ) THEN
!
!                 Update A12
!
                  CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                              JB, J2, ONE, AB( KV+1, J ), LDAB-1, &
                              AB( KV+1-JB, J+JB ), LDAB-1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A22
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I2, J2, &
                                 JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
                                 AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
                                 AB( KV+1, J+JB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A32
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I3, J2, &
                                 JB, -ONE, WORK31, LDWORK, &
                                 AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
                                 AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
!
               IF( J3.GT.0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
!
!                 Update A13 in the work array
!
                  CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                              JB, J3, ONE, AB( KV+1, J ), LDAB-1, &
                              WORK13, LDWORK )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A23
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I2, J3, &
                                 JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
                                 WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), &
                                 LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A33
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I3, J3, &
                                 JB, -ONE, WORK31, LDWORK, WORK13, &
                                 LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
!
!                 Copy the lower triangle of A13 back into place
!
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
!
!              Adjust the pivot indices.
!
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
                  IF( JP+JJ-1.LT.J+KL ) THEN
!
!                    The interchange does not affect A31
!
                     CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                 AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
!
!                    The interchange does affect A31
!
                     CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                 WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
!
!              Copy the current column of A31 back into place
!
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 ) &
                  CALL ZCOPY( NW, WORK31( 1, JJ-J+1 ), 1, &
                              AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
!
      RETURN
!
!     End of ZGBTRF
!
      END SUBROUTINE ZGBTRF

!<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      SUBROUTINE ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGBTF2 computes an LU factorization of a complex m-by-n band matrix
!  A using partial pivoting with row interchanges.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U, because of fill-in resulting from the row
!  interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
!     ..
!     .. External Functions ..
      INTEGER            IZAMAX
      EXTERNAL           IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
         RETURN
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
      JU = 1
!
      DO 40 J = 1, MIN( M, N )
!
!        Set fill-in elements in column J+KV to zero.
!
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
         KM = MIN( KL, M-J )
         JP = IZAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
!
!           Apply interchange to columns J to JU.
!
            IF( JP.NE.1 ) &
               CALL ZSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1, &
                           AB( KV+1, J ), LDAB-1 )
            IF( KM.GT.0 ) THEN
!
!              Compute multipliers.
!
               CALL ZSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
!
!              Update trailing submatrix within the band.
!
               IF( JU.GT.J ) &
                  CALL ZGERU( KM, JU-J, -ONE, AB( KV+2, J ), 1, &
                              AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), &
                              LDAB-1 )
            END IF
         ELSE
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
            IF( INFO.EQ.0 ) &
               INFO = J
         END IF
   40 CONTINUE
      RETURN
!
!     End of ZGBTF2
!
      END SUBROUTINE ZGBTF2

!<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IP, IX
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZSWAP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.EQ.0 ) &
         RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I ) &
               CALL ZSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I ) &
               CALL ZSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I ) &
               CALL ZSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLASWP
!
      END SUBROUTINE ZLASWP

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, &
                       N4 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. Executable Statements ..
!
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
             ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                   ( IC.GE.162 .AND. IC.LE.169 ) ) &
                  SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) &
         RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      GO TO ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                  C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
             C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
             C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE 
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
!     End of ILAENV
!
      END FUNCTION ILAENV

!<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
            'an illegal value' )
!
!     End of XERBLA
!
      END SUBROUTINE XERBLA

!<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      SUBROUTINE ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, &
                         INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGBTRS solves a system of linear equations
!     A * X = B,  A**T * X = B,  or  A**H * X = B
!  with a general band matrix A using the LU factorization computed
!  by ZGBTRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B     (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) COMPLEX*16 array, dimension (LDAB,N)
!          Details of the LU factorization of the band matrix A, as
!          computed by ZGBTRF.  U is stored as an upper triangular band
!          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!          the multipliers used during the factorization are stored in
!          rows KL+KU+2 to 2*KL+KU+1.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= N, row i of the matrix was
!          interchanged with row IPIV(i).
!
!  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMV, ZGERU, ZLACGV, ZSWAP, ZTBSV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
          LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGBTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
         RETURN
!
      KD = KU + KL + 1
      LNOTI = KL.GT.0
!
      IF( NOTRAN ) THEN
!
!        Solve  A*X = B.
!
!        Solve L*X = B, overwriting B with X.
!
!        L is represented as a product of permutations and unit lower
!        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!        where each transformation L(i) is a rank-one modification of
!        the identity matrix.
!
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J ) &
                  CALL ZSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL ZGERU( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), &
                           LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
!
         DO 20 I = 1, NRHS
!
!           Solve U*X = B, overwriting B with X.
!
            CALL ZTBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU, &
                        AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
!
      ELSE IF( LSAME( TRANS, 'T' ) ) THEN
!
!        Solve A**T * X = B.
!
         DO 30 I = 1, NRHS
!
!           Solve U**T * X = B, overwriting B with X.
!
            CALL ZTBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, &
                        LDAB, B( 1, I ), 1 )
   30    CONTINUE
!
!        Solve L**T * X = B, overwriting B with X.
!
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL ZGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ), &
                           LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J ) &
                  CALL ZSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
!
      ELSE
!
!        Solve A**H * X = B.
!
         DO 50 I = 1, NRHS
!
!           Solve U**H * X = B, overwriting B with X.
!
            CALL ZTBSV( 'Upper', 'Conjugate transpose', 'Non-unit', N, &
                        KL+KU, AB, LDAB, B( 1, I ), 1 )
   50    CONTINUE
!
!        Solve L**H * X = B, overwriting B with X.
!
         IF( LNOTI ) THEN
            DO 60 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL ZLACGV( NRHS, B( J, 1 ), LDB )
               CALL ZGEMV( 'Conjugate transpose', LM, NRHS, -ONE, &
                           B( J+1, 1 ), LDB, AB( KD+1, J ), 1, ONE, &
                           B( J, 1 ), LDB )
               CALL ZLACGV( NRHS, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J ) &
                  CALL ZSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   60       CONTINUE
         END IF
      END IF
      RETURN
!
!     End of ZGBTRS
!
      END SUBROUTINE ZGBTRS

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      SUBROUTINE ZLACGV( N, X, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         X( * )
!     ..
!
!  Purpose
!  =======
!
!  ZLACGV conjugates a complex vector of length N.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The length of the vector X.  N >= 0.
!
!  X       (input/output) COMPLEX*16 array, dimension
!                         (1+(N-1)*abs(INCX))
!          On entry, the vector of length N to be conjugated.
!          On exit, X is overwritten with conjg(X).
!
!  INCX    (input) INTEGER
!          The spacing between successive elements of X.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IOFF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = DCONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 ) &
            IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = DCONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
!
!     End of ZLACGV
!
      END SUBROUTINE ZLACGV
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      LOGICAL          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME ) &
         RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR. &
             INTA.GE.145 .AND. INTA.LE.153 .OR. &
             INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR. &
             INTB.GE.145 .AND. INTB.LE.153 .OR. &
             INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END FUNCTION LSAME
!******************************************************************************
!  Given a matrix A, with logical dimensions M by N and physical dimensions MP by NP, this
!  routine computes its singular value decomposition, A = U.W.V'. The matrix U replaces 
!  A on output. The diagonal matrix of singular values W is output as a vector W. The matrix 
!  V (not the Transpose V') is output as V. M must be greater or equal to N; if it is smaller,
!  then A should be filled up to square with zero rows.
!******************************************************************************
!  Verified: working on 6 Feb 1998
!
! NOTES: 
! 1. When used to solve Linear Equations with n equations and n unknowns,
!     mp=np=m=n.
! 2. When finding inverses, n = soln and b = IdentityMatrix(n)
!
! Modifications to Original Program:
! 1. Equalties/Inequalities are replaced by Differences and compared with EPS
! 2. The Matrix U in U.W.V' is actually stored in "a"
!
! DEBUG:
! 1.  "IMPLICIT NONE" and "REAL(DBP)"
! 2.  Parameters might need to be larger
!******************************************************************************
SUBROUTINE SVDCMP(a,m,n,mp,np,w,v)
   implicit none
   INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
   
   INTEGER :: m,n,mp,np, i,j,l,k, its, nm, jj
   REAL(dbp) :: g, sscale, anorm, s, f,h,x,z,y,c

   REAL(dbp), parameter :: EPS = 3.0d-15

   INTEGER, PARAMETER :: nmax = 1000        !Maximum anticipated value of N
   REAL(dbp) :: A(mp,np), w(np), v(np,np), rv1(nmax)


!   PRINT *, 'Enter the tolerance or precision required'
!   READ *, eps
   !eps=10E-12

   if (m.lt.n) then
   	PRINT *, 'You must augment A with extra zero rows'
      call exit(10)
   ENDIF 

            !Householder Reduction to bidiagonal form
!(see Forsythe,Malcolm,Moler, "Computer Methods for Mathematical Computations"
   g=0.0d0
   sscale = 0.0d0
   anorm = 0.0d0
   do i = 1,n
      l = i + 1
      rv1(i) = sscale*g
      g = 0.0d0
      s = 0.0d0
      sscale = 0.0d0
      if (i.le.m) then
         do k = i,m
            sscale = sscale + dABS(a(k,i))
         end do       ! k loop
!         if (sscale.ne.0.0d0) then
			if ( dabs(sscale-0.0d0).gt.EPS ) then
            do k = i,m
               a(k,i) = a(k,i) / sscale
               s = s + a(k,i)*a(k,i)
            end do    ! k loop
            f = a(i,i)
            g = - SIGN(SQRT(s),f)
            h = f*g - s
            a(i,i) = f - g
            if (i.ne.n) then
               do j = l,n
                  s = 0.0d0
                  do k = i,m
                     s = s + a(k,i)*a(k,j)
                  end do      ! k loop
                  f = s / h
                  do k = i, m 
                     a(k,j) = a(k,j) + f*a(k,i)
                  end do   ! k loop
               end do      ! j loop
            end if
            do k = i, m 
               a(k,i) = sscale * a(k,i)
            end do         ! k loop
         end if
      end if

      w(i) = sscale * g
      g = 0.0d0
      s = 0.0d0
      sscale = 0.0d0
      if ((i.le.m).AND.(i.ne.n)) then
         do k = l, n
            sscale = sscale + dABS(a(i,k))
         end do         ! k loop
!         if (sscale.ne.0.0d0) then
			if ( dabs(sscale-0.0d0).gt.EPS ) then
            do k = l, n
               a(i,k) = a(i,k) /sscale
               s = s + a(i,k) * a(i,k)
            end do      ! k loop 
            f = a(i,l) 
            g = - SIGN(SQRT(s),f)
            h = f * g - s
            a(i,l) = f - g
            do k = l, n
               rv1(k) = a(i,k) / h
            end do      ! k loop
            if (i.ne.m) then
               do j = l, m 
                  s = 0.0d0
                  do k = l, n 
                     s = s + a(j,k)*a(i,k)
                  end do   ! k loop
                  do k = l, n 
                     a(j,k) = a(j,k) + s*rv1(k)
                  end do   ! k loop
               end do      ! j loop
            end if
				do k = l, n
               a(i,k) = sscale * a(i,k)
           	end do
         end if
      end if
      anorm = MAX(anorm, (dABS(w(i)) + dABS(rv1(i))))
   end do

! Accumulation of right-hand Transformations
   do i = n, 1, -1 
      if (i.lt.n) then
!         if (g.ne.0.0d0) then
			if ( dabs(g-0.0d0).gt.EPS ) then
            do j = l, n       ! Double division to avoid possible overflow
               v(j,i) = (a(i,j) / a(i,l)) / g
            end do      ! j loop
            do j = l, n
               s = 0.0d0
               do k = l, n
                  s = s + a(i,k)*v(k,j)
               end do   ! k loop
               do k = l, n
                  v(k,j) = v(k,j) + s * v(k,i)
               end do   ! k loop
           	end do      ! j loop
         end if
         do j = l, n 
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
         end do
      end if
      v(i,i) = 1.0d0
      g = rv1(i)
      l = i
   end do

! Accumulation of left-hand Transformations
   do i = n, 1, -1
      l = 1 + i
      g = w(i)
      if (i.lt.n) then
         do j = l, n
            a(i,j) = 0.0d0
         end do
      end if
!      if (g.ne.0.0d0) then
      if ( dabs(g-0.0d0).gt.EPS ) then
         g = 1.0d0 / g
         if (i.ne.n) then
            do j = l,n 
               s = 0.0d0
               do k = l, m
                  s = s + a(k,i)*a(k,j)
               end do   ! k loop
               f = (s/a(i,i)) * g
               do k = i, m 
                  a(k,j) = a(k,j) + f * a(k,i)
               end do   ! k loop
            end do      ! j loop
         end if
         do j = i, m 
            a(j,i) = a(j,i) * g
         end do         ! j loop
      else
         do j = i, m
            a(j,i) = 0.0d0
         end do         ! j loop
      end if
      a(i,i) = a(i,i) + 1.0d0
   end do               ! i loop

! Diagonalization of the bidigonal form
   do k = n, 1, -1                  !Loop over singular values
      do its = 1,30                 !Loop over allowed iterations
         do l = k, 1, -1            !Test for splitting
            nm = l - 1              ! Note that rv1(1) is always zero
!           if ( (dABS(rv1(l))+anorm) .eq. anorm ) GO TO 2
!          	if ( (dABS(w(nm))+anorm) .eq. anorm ) GO TO 1
            if ( dabs((dABS(rv1(l))+anorm) - anorm).lt.eps ) GO TO 2
          	if ( dabs((dABS(w(nm))+anorm) - anorm).lt.eps ) GO TO 1
         end do      !  l loop

1        c = 0.0d0                  ! Cancellation of rv1(l), if l>1 :
         s = 1.0d0
         do i = l, k
            f = s * rv1(i)
!            if ( (dABS(f)+anorm) .ne. anorm ) then
            if ( dabs( (dABS(f)+anorm) - anorm) .GT. eps ) then

               g = w(i)
               h = SQRT(f*f + g*g)
               w(i) = h
               h = 1.0d0 / h
               c = g * h
               s = -f * h
               do j = 1, m
                  y = a(j,nm)
                  z = a(j,i)
                  a(j,nm) = (y*c) + (z*s)
                  a(j,i) = -(y*s) + (z*c)
               end do   ! j loop
            end if
         end do         ! i loop
2        z = w(k) 
         if (l .eq. k) then         ! convergence
				if (z .lt. 0.0d0) then  ! Singular value is made non-negative
               w(k) = -z
               do j = 1,n
                  v(j,k) = -v(j,k)
               end do         ! j loop
	    		end if
            GO TO 3
         end if
         if (its.eq.30) then
         	PRINT*, 'No Convergence in 30 iterations'
            call exit(10)
         ENDIF
         x = w(l)          ! Shift from bottom 2-by-2 minor
         nm = k - 1
         y = w(nm)
         g = rv1(nm)
         h = rv1(k)
         f = ( (y-z)*(y+z) + (g-h)*(g+h) ) / ( 2.0d0*h*y)
         g = SQRT(f*f + 1.0d0)
         f = ( (x-z)*(x+z) + h*((y/(f+SIGN(g,f))) - h) ) / x

! Next   QR Transformation
         c = 1.0d0
         s = 1.0d0
         do j = l, nm
          	i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
            z = SQRT(f*f + h*h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c) + (g*s)
            g = -(x*s) + (g*c)
            h = y*s
            y = y*c
            do jj = 1, n
               x = v(jj,j)
               z = v(jj,i)
               v(jj,j) = (x*c) + (z*s)
               v(jj,i) = -(x*s) + (z*c)
           	end do
            z = SQRT(f*f + h*h)
            w(j) = z
!            if (z.ne.0.0d0) then
            if (  dabs(z-0.0d0).gt.eps  ) then
               z = 1.0d0 / z
               c = f*z
               s = h*z
            end if
            f = (g*c) + (y*s)
            x = -(g*s) + (y*c)
            do jj = 1, m
               y = a(jj,j)
               z = a(jj,i)
               a(jj,j) = (y*c) + (z*s)
               a(jj,i) = -(y*s) + (z*c)
            end do
         end do         ! j loop
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
      end do            ! its loop
3  continue
   end do               ! k loop

   return
END SUBROUTINE svdcmp


subroutine SVDSORT(U, W, V,m,n)
    !' Sort U, V, W  by Singular Values in decending order
    ! '   [meaning W(0) will be greatest singular value]
     INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)   
    real(dbp):: U(m,n),w(m),V(m,n)
    Integer :: i,j,k,m,n   
    real(dbp):: S      


    do i = 1 , N - 1
        K = i                    !Find next highest singular value index k
        S = W(K)
        do J = i + 1, N
            If (W(J) .gt. S) Then
                K = J
                S = W(K)
            End If
        end do
        If (K .ne. i) Then
            W(K) = W(i)          !Swap W(k), W(i)
            W(i) = S
            do J = 1, N       !Swap V(Row i), V(Row k)
                S = V(J, i)
                V(J, i) = V(J, K)
                V(J, K) = S
            end do
            do J = 1, M       !Swap U(Row i), U(Row k)
                S = U(J, i)
                U(J, i) = U(J, K)
                U(J, K) = S
            end do
        End If
        end do
end subroutine SVDSORT        

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!end module
