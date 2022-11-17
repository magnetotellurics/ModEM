! *****************************************************************************
!  This module contains routines from rebocc that are used
!  for setting up model operators used for solving 1d TE
!  and TM problems.  Minor modifications from the original
!  code allow dynamic memory allocation. 
module WSfwd1Dmod
use wsfwd2dpar
use wsutils
!use wsLAPACK
!use sunperf

   Contains
   
   ! *****************************************************************************
      SUBROUTINE Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)

      INTEGER Nzb
      REAL(kind=prec)  per,Dzb(*),r1d(*)
      COMPLEX(kind=prec) X1D(*)
      
      INTEGER iz,jj,nzb1,nzb2
      REAL(kind=prec)  Omega,Omue,skd
      REAL(kind=prec)  dz1(NZ1MX),cz1(NZ2MX)
      REAL(kind=prec)  au1d(NZ0MX),ad1d(NZ0MX)
      COMPLEX(kind=prec) ac1d(NZ0MX),xb1d(2),xb(NZ0MX)

      nzb1 = Nzb + 1
      nzb2 = Nzb + 2
      Omega = (TWO*PI)/Per
      Omue  = Omega*Mue
 
      r1d(nzb1) = r1d(Nzb)
      skd  = DSQRT(TWO*r1d(nzb1)/Omue)
      CALL CopyVectorR8(1,Nzb,Dzb,1,Nzb,dz1)
      dz1(nzb1) = skd
      CALL DistanceBetweenBlocks(nzb1,dz1,cz1)

!     Assign operators
      CALL ConstantVectorR8(au1d,nzb1,R_ZERO)
      CALL ConstantVectorR8(ad1d,nzb1,R_ZERO)
      CALL ConstantVectorC16(ac1d,nzb1,R_ZERO)
      jj = 1
      DO iz = 2,nzb1
        au1d(jj) = TWO*r1d(iz-1)/dz1(iz-1)
        ad1d(jj) = TWO*r1d(iz)/dz1(iz)
        ac1d(jj) = DCMPLX(R_ZERO,Omue)*cz1(iz) - au1d(jj) - ad1d(jj)
        jj = jj + 1
      ENDDO

!     Boundary condition for 1D
      xb1d(1) = ONE
      xb1d(2) = R_ZERO

!     solve Aii*Xi = -Aib*Xb
      CALL ConstantVectorC16(xb,nzb,R_ZERO)
      xb(1) = -au1d(1)*xb1d(1)

      CALL Solve1D(nzb,ad1d,ac1d,xb)
      X1D(1) = xb1d(1)
      CALL CopyVectorC16(1,nzb,xb,2,nzb1,X1D)

100   FORMAT(7e11.3)

      Return
      END Subroutine ! Fwd1d_TM


! *****************************************************************************
      SUBROUTINE Fwd1D_TE(per,Nza,Nz,Dz,s1,X1D)

      INTEGER Nza,Nz
      REAL(kind=prec)  per,Dz(*),s1(*)
      COMPLEX(kind=prec) X1D(*)
      
      INTEGER iz,jj,nz1,nz2
      REAL(kind=prec)  Omega,Omue,skd,s1d(NZ1MX)
      REAL(kind=prec)  dz1(NZ1MX),cz1(NZ2MX)
      REAL(kind=prec)  au1d(NZ0MX),ad1d(NZ0MX),ss
      COMPLEX(kind=prec) ac1d(NZ0MX),xb1d(2),xb(NZ0MX)

      nz1 = Nz + 1
      nz2 = Nz + 2
      Omega = (TWO*PI)/Per
      Omue  = Omega*Mue

      DO iz = 1,Nza
        s1d(iz) = CondAir
      ENDDO
      DO iz = Nza+1,Nz
        s1d(iz) = s1(iz-Nza)
      ENDDO
      s1d(nz1) = s1d(Nz)

      skd  = DSQRT(TWO/(s1d(nz1)*Omue))
      DO iz = 1,Nz
        dz1(iz) = Dz(iz)
      ENDDO ! iz

!     CALL CopyVectorR8(1,Nz,Dz,1,Nz,dz1)
      dz1(nz1) = skd
      CALL DistanceBetweenBlocks(nz1,dz1,cz1)

!     Assign operators
      CALL ConstantVectorR8(au1d,nz1,R_ZERO)
      CALL ConstantVectorR8(ad1d,nz1,R_ZERO)
      CALL ConstantVectorC16(ac1d,nz1,R_ZERO)
      jj = 1
      DO iz = 2,nz1
        au1d(jj) = TWO/dz1(iz-1)
        ad1d(jj) = TWO/dz1(iz)
        ss       = (s1d(iz)*dz1(iz)+s1d(iz-1)*dz1(iz-1))
        ac1d(jj) = DCMPLX(R_ZERO,Omue*ss) - au1d(jj) - ad1d(jj)
        jj = jj + 1
      ENDDO

!     Boundary condition for 1D
      xb1d(1) = ONE
      xb1d(2) = R_ZERO

!     solve Aii*Xi = -Aib*Xb
      CALL ConstantVectorC16(xb,Nz,R_ZERO)
      xb(1) = -au1d(1)*xb1d(1)

      CALL Solve1D(Nz,ad1d,ac1d,xb)
      X1D(1) = xb1d(1)
      CALL CopyVectorC16(1,nz,xb,2,nz1,X1D)


100   FORMAT(7e11.3)

      Return
      END Subroutine ! Fwd1d_TE


! *****************************************************************************
      SUBROUTINE Solve1D(nz0,ad1d,ac1d,xb)
 
      INTEGER nz0
      REAL(kind=prec)  ad1d(*)
      COMPLEX(kind=prec) ac1d(*),xb(*)

      INTEGER jj,iPiv(NZ0MX),info
      COMPLEX(kind=prec) a1d(4,NZ0MX)
      character*1 cN
     
      info = 0 
      cN = 'N'
!     make matrix in LAPACK's format
      CALL ConstantMatrixC16(a1d,4,NZ0MX,4,nz0,R_ZERO)

      DO jj = 2,nz0
        a1d(2,jj) = ad1d(jj-1)
      ENDDO
      DO jj = 1,nz0
        a1d(3,jj) = ac1d(jj)
      ENDDO
      DO jj = 1,nz0-1
        a1d(4,jj) = ad1d(jj)
      ENDDO
      CALL ZGBTRF(nz0,nz0,1,1,a1d,4,iPiv,info)
      CALL ZGBTRS(cN,nz0,1,1,1,a1d,4,iPiv,xb,nz0,info)

      END Subroutine ! Solve1D
      
      
end module
