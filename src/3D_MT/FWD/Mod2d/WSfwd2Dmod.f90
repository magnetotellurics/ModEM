! *****************************************************************************
!  This module contains routines from rebocc that are used
!  for setting up model operators used for solving 2d TE
!  and TM problems.  Minor modifications from the original
!  code allow dynamic memory allocation.  Fixed array sizes
!  defined in "parameter.h" in original code are now obtained
!  by using module  wsfwd2dparams.mod .  Integer array size
!  parameters in this module are set by a call to FWD2DTMsetup
!  (FWD2DTEsetup) in FwdTMmod (FwdTEmod).  

!   In a few routines
!  where the parameters were used to dimension local arrays,
!  allocate statements have also need to be added.  

!   Minimal documentation, consistent with original code!

module WSfwd2Dmod
   use wsfwd1dmod
   !use wsLAPACK
   !use sunperf
   implicit none
   
   !  grid2d_t is derived data type used to store basic grid geometry info
   type ::  grid2d_t
      integer	:: Nz,Ny,Nza
      real(kind=prec), pointer, dimension(:) :: Dy,Dz
   end type

   contains
   
   
   ! *****************************************************************************
      SUBROUTINE FormAII(per,Nz0,Ny,AA,AII,iPiv)
 
      INTEGER Nz0,Ny,iPiv(*)
      REAL(kind=prec)  per,AA(MMIMX,4)
      COMPLEX(kind=prec) AII(NZ3MX,MMIMX)

      INTEGER mmi,jj,nz3,kl,ku,kc,info
      REAL(kind=prec)  Omega,Omue

      Omega = (TWO*PI)/per
      Omue  = Omega*Mue

!     Forming matrix Aii in LAPACK's format
      nz3 = 3*(Nz0-1) + 1
      mmi = (Ny-1)*(Nz0-1)
      CALL ConstantMatrixC16(AII,NZ3MX,MMIMX,nz3,mmi,R_ZERO)

      kl = Nz0-1
      ku = Nz0-1
      kc = kl + ku + 1

!     diagonal term
      DO jj = 1,mmi
        AII(kc,jj) = DCMPLX(AA(jj,1),Omue*AA(jj,4))
      ENDDO

!     first lower and upper diagonal terms
      DO jj = 1,mmi-1
        AII(kc-1,jj+1) = AA(jj,2)
        AII(kc+1,jj)   = AA(jj,2)
      ENDDO

!     second lower and upper diagonal terms
      DO jj = 1,mmi-(Nz0-1)
        AII(kc-(Nz0-1),jj+(Nz0-1)) = AA(jj,3)
        AII(kc+(Nz0-1),jj)         = AA(jj,3)
      ENDDO
  
!     LU Decompose 
      CALL ZGBTRF(mmi,mmi,kl,ku,AII,NZ3MX,iPiv,info)
      IF (info.NE.0) THEN
        WRITE(6,*)   &
       '!!! ATTENTION, ERROR WHILE DECOMPOSING MATRIX !!!'
        WRITE(6,*) '!!! Please, check your model and restart !!!'
        WRITE(6,*)   &
       '!!! If problem persists, contact egbert@coas.oregonstate.edu !!'
        STOP
      ENDIF

100   FORMAT(7E11.3)
      
      RETURN
      END  Subroutine ! FormAII()
      
      
      ! *****************************************************************************
      SUBROUTINE MulAibWithXb(Nz0,Ny,Aib,Xb,AXb)

      INTEGER Nz0,Ny
      REAL(kind=prec)  Aib(*)
      COMPLEX(kind=prec) Xb(*),AXb(*)

      INTEGER mmi,mmb,iy,iz,jj,izz

      mmi = (Ny-1)*(Nz0-1)
      mmb = 2*Ny + 2*Nz0

      CALL ConstantVectorC16(AXb,mmi,R_ZERO)

!     left side
      iy = 2
      DO iz = 2,Nz0
        jj = (iy-2)*(Nz0-1) + iz-1
        AXb(jj) = - Aib(iz)*Xb(iz)
      ENDDO

!     right side
      iy = Ny
      izz = 2*Ny + Nz0 + 1
      DO iz = 2,Nz0
        jj = (iy-2)*(Nz0-1) + iz-1
        AXb(jj) = - Aib(izz)*Xb(izz)
        izz = izz + 1
      ENDDO

!     top side
      izz = Nz0 + 2
      iz  = 2
      DO iy = 2,Ny
        jj = (iy-2)*(Nz0-1) + iz-1
        AXb(jj) = AXb(jj) - Aib(izz)*Xb(izz)
        izz = izz + 1
      ENDDO

!     bottom side
      izz = Nz0 + Ny + 1
      iz  = Nz0
      DO iy = 2,Ny
        jj = (iy-2)*(Nz0-1) + iz-1
        AXb(jj) = AXb(jj) - Aib(izz)*Xb(izz)
        izz = izz + 1
      ENDDO

      END Subroutine! MulAibWithXb
!
      SUBROUTINE SetBound2D_TM(per,Nzb,Ny,Dzb,Dy,CRho,HXI0,HXB)

      INTEGER Nzb,Ny
      REAL(kind=prec)  per,Dzb(*),Dy(*),CRho(NZ0MX,NY0MX)
      COMPLEX(kind=prec) HXI0(*),HXB(*)

      INTEGER nzb1,ny1,mmi,mmb,np1,np2,iz,iy,ii,jj
      REAL(kind=prec)  r1d(NZ1MX)
      COMPLEX(kind=prec) X1D(NZ2MX)
!     is this extra complication needed?
!      real (kind=8), pointer, dimension(:)	::  r1d
!      complex (kind=8), pointer, dimension(:)	::  x1d
!      allocate(r1d(NZ1MX))
!      allocate(x1d(NZ2MX))

      nzb1 = Nzb + 1
      ny1  = Ny  + 1

      mmb = 2*Nzb + 2*Ny
      mmi = (Nzb-1)*(Ny-1)
   
      CALL ConstantVectorC16(HXB,mmb,R_ZERO)
      CALL ConstantVectorC16(HXI0,mmi,R_ZERO)

      np1 = NZ0MX
      np2 = NY0MX

!    left fields
      CALL TransMatrixToVectorR8(CRho,np1,np2,1,1,Nzb,1,r1d,1,Nzb)
      CALL Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)
      CALL CopyVectorC16(1,Nzb+1,X1D,1,Nzb+1,HXB)

!     surface fields
      DO ii = Nzb+2,Nzb+Ny
        HXB(ii) = ONE
      ENDDO


!     Interior and bottom fields
      DO iy = 2,Ny
        DO iz = 1,Nzb
          r1d(iz) = (CRho(iz,iy-1)*Dy(iy-1)+CRho(iz,iy)*Dy(iy))/ &
                   (Dy(iy-1) + Dy(iy))
        ENDDO
        CALL Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)
        HXB(Nzb+Ny-1+iy) = X1D(Nzb+1)
        DO iz = 2,Nzb
          jj = (iy-2)*(Nzb-1) + iz-1
          HXI0(jj) = X1D(iz)
        ENDDO
      ENDDO 

!    right fields
      CALL TransMatrixToVectorR8(CRho,np1,np2,1,Ny,Nzb,Ny,r1d,1,Nzb)
      CALL Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)
      CALL CopyVectorC16(1,Nzb+1,X1D,Nzb+2*Ny,2*Nzb+2*Ny,HXB)
      

100   FORMAT(i5,2e11.3)

      RETURN
      END Subroutine! SetBound2D_TM


! *****************************************************************************
      SUBROUTINE SetBound2D_TE(per,Nza,Nz,Ny,Dz,Dy,Cond2d,EXI0,EXB)
!     This routine modified to take Cond2d as input, instead of CRho
!        to minimize changes to WS code, just copy (part of) Cond2d into CCon 

      INTEGER Nza,Nz,Ny
      REAL(kind=prec)  per,Dz(*),Dy(*)
      COMPLEX(kind=prec) EXI0(*),EXB(*)

      INTEGER nz1,ny1,mmi,mmb,np1,np2,iz,iy,ii,jj
      REAL(kind=prec)  s1d(NZ1MX),CCon(NZ0MX,NY0MX),Cond2D(NZ0MX,NY0MX)
      COMPLEX(kind=prec) X1D(NZ2MX)


      nz1 = Nz + 1
      ny1 = Ny + 1

      mmb = 2*Nz + 2*Ny
      mmi = (Nz-1)*(Ny-1)

      DO iy = 1,Ny
        DO iz = 1,Nz-Nza
          CCon(iz,iy) = Cond2D(iz+Nza,iy)
        ENDDO ! iz
      ENDDO ! iy


      CALL ConstantVectorC16(EXB,mmb,R_ZERO)
      CALL ConstantVectorC16(EXI0,mmi,R_ZERO)

      np1 = NZ0MX
      np2 = NY0MX

!     left fields
      CALL TransMatrixToVectorR8(CCon,np1,np2,1,1,Nz-Nza,1,  &
                                s1d,1,Nz-Nza)
      CALL Fwd1D_TE(per,Nza,Nz,Dz,s1d,X1D)
      CALL CopyVectorC16(1,Nz+1,X1D,1,Nz+1,EXB)


!     surface fields
      DO ii = Nz+2,Nz+Ny
        EXB(ii) = ONE
      ENDDO


!     Interior and bottom fields
      DO iy = 2,Ny
        DO iz = 1,Nz-Nza
          s1d(iz) = (CCon(iz,iy-1)*Dy(iy-1)+CCon(iz,iy)*Dy(iy))/  &
                   (Dy(iy-1) + Dy(iy))
        ENDDO
        CALL Fwd1D_TE(per,Nza,Nz,Dz,s1d,X1D)
        EXB(Nz+Ny-1+iy) = X1D(Nz+1)
        DO iz = 2,Nz
          jj = (iy-2)*(Nz-1) + iz-1
          EXI0(jj) = X1D(iz)
        ENDDO
      ENDDO 

!     right fields
      CALL TransMatrixToVectorR8(CCon,np1,np2,1,Ny,Nz-Nza,Ny,  &
                                s1d,1,Nz-Nza)
      CALL Fwd1D_TE(per,Nza,Nz,Dz,s1d,X1D)
      CALL CopyVectorC16(1,Nz+1,X1D,Nz+2*Ny,2*Nz+2*Ny,EXB)

100   FORMAT(i5,2e11.3)

      RETURN
      END Subroutine! SetBound2D_TE
!

      SUBROUTINE SetupA_TM(Nzb,Ny,Dzb,Dy,Czb,Cy,CRho,ATM,BTM)

      INTEGER Nzb,Ny
      real(kind=prec)  Dzb(*),Dy(*),Czb(*),Cy(*),CRho(NZ0MX,NY0MX)
      real(kind=prec)  ATM(MMIMX,4),BTM(MMBMX)
      
!    ATM(:,1) : Real diagonal term
!    ATM(:,2) : Real first upper strip
!     ATM(:,3) : Real second upper strip
!     ATM(:,4) : Imaginary diagonal term
!     BTM      : Boundary

      INTEGER jj,iz,iy,nblx,iblx
      real(kind=prec)  r00,r10,r01,r11
      real(kind=prec)  ar,al,ad,au

      nblx = Ny*Nzb
      iblx = 2*Ny + 2*Nzb

      CALL ConstantMatrixR8(ATM,MMIMX,4,nblx,4,R_ZERO)
      CALL ConstantVectorR8(BTM,iblx,R_ZERO)

      jj = 1
      DO iy = 2,Ny
        DO iz = 2,Nzb
          r00 = CRho(iz,iy)
          r10 = CRho(iz-1,iy)
          r01 = CRho(iz,iy-1)
          r11 = CRho(iz-1,iy-1)

          ar  = TWO*(Dzb(iz)*r00 + Dzb(iz-1)*r10)/Dy(iy)
          al  = TWO*(Dzb(iz)*r01 + Dzb(iz-1)*r11)/Dy(iy-1)
          ad  = TWO*(Dy(iy)*r00  + Dy(iy-1)*r01)/Dzb(iz)
          au  = TWO*(Dy(iy)*r10  + Dy(iy-1)*r11)/Dzb(iz-1)
          ATM(jj,1) = - ar - al - ad - au
          ATM(jj,4) = Cy(iY)*Czb(iz)
          ATM(jj,2) = ad
          IF (iz.EQ.Nzb) THEN
            ATM(jj,2) = R_ZERO
          ENDIF
          IF (iy.LT.NY) THEN 
            ATM(jj,3) = ar
          ENDIF
          jj = jj + 1
        ENDDO ! iz
      ENDDO ! iy

!     boundary fields

!     left side fields
      iy = 2
      DO iz = 2,Nzb
        r01 = CRho(iz,iy-1)
        r11 = CRho(iz-1,iy-1)
        BTM(iz) = TWO*(Dzb(iz)*r01 + Dzb(iz-1)*r11)/Dy(iy-1)
      ENDDO ! iz

!     surface fields
      iz = 2
      DO iy = 2,Ny
        r10 = CRho(iz-1,iy)
        r11 = CRho(iz-1,iy-1)
        BTM(Nzb+iy) = TWO*(Dy(iy)*r10  + Dy(iy-1)*r11)/Dzb(iz-1)
      ENDDO ! iy

!     Bottom fields
      iz = Nzb
      DO iy = 2,Ny
        r00 = CRho(iz,iy)
        r01 = CRho(iz,iy-1)
        BTM(Nzb+Ny+iy-1) = TWO*(Dy(iy)*r00  + Dy(iy-1)*r01)/Dzb(iz)
      ENDDO ! iy

!     Right fields
      iy = Ny
      DO iz = 2,Nzb
        r00 = CRho(iz,iy)
        r10 = CRho(iz-1,iy)
        BTM(Nzb+2*Ny-1+iz) = TWO*(Dzb(iz)*r00 + Dzb(iz-1)*r10)/Dy(iy)
      ENDDO ! iz


      RETURN
      END Subroutine! SetupA_TM()


! *****************************************************************************
      SUBROUTINE SetupA_TE(Nza,Nz,Ny,Dz,Dy,Cz,Cy,Cond2D,ATE,BTE)
!     This routine modified to take Cond2D as input, instead of CRho
!        to minimize changes to WS code, just copy (part of) Cond2d into CCon

      INTEGER Nza,Nz,Ny
      REAL(kind=prec)  Dz(*),Dy(*),Cz(*),Cy(*)
      REAL(kind=prec)  ATE(MMIMX,4),BTE(MMBMX)
      
!     ATE(:,1) : Real diagonal term
!     ATE(:,2) : Real second strip
!     ATE(:,3) : Real thrid strip
!     ATE(:,4) : Imaginary diagonal term
!     BTE      : Boundary

      INTEGER jj,iz,iy,nblx,iblx
      REAL(kind=prec)  s00,s10,s01,s11,s2d
      REAL(kind=prec)  ar,al,ad,au,CCon(NZ0MX,NY0MX),COnd2D(NZ0MX,NY0MX)

      nblx = Ny*Nz
      iblx = 2*Ny + 2*Nz

      CALL ConstantMatrixR8(CCon,NZ0MX,NY0MX,Nza,Ny,CondAir)

      DO iy = 1,Ny
        DO iz = Nza+1,Nz
          CCon(iz,iy) = Cond2D(iz,iy)
        ENDDO ! iz
      ENDDO ! iy

      CALL ConstantMatrixR8(ATE,MMIMX,4,nblx,4,R_ZERO)
      CALL ConstantVectorR8(BTE,iblx,R_ZERO)

      jj = 1
      DO iy = 2,Ny
        DO iz = 2,Nz
          s00 = CCon(iz,iy)
          s10 = CCon(iz-1,iy)
          s01 = CCon(iz,iy-1)
          s11 = CCon(iz-1,iy-1)
          s2d = (s00*Dz(iz)*Dy(iy)   + s01*Dz(iz)*Dy(iy-1) + &
                s10*Dz(iz-1)*Dy(iy) + s11*Dz(iz-1)*Dy(iy-1))/ &
               (Dz(iz)*Dy(iy)   + Dz(iz)*Dy(iy-1) +  &
               Dz(iz-1)*Dy(iy) + Dz(iz-1)*Dy(iy-1)) 

          ar  = TWO*Cz(iz)/Dy(iy)
          al  = TWO*Cz(iz)/Dy(iy-1)
          ad  = TWO*Cy(iy)/Dz(iz)
          au  = TWO*Cy(iy)/Dz(iz-1)
          ATE(jj,1) = - ar - al - ad - au
          ATE(jj,4) = Cy(iy)*Cz(iz)*s2d

          ATE(jj,2) = ad
          IF (iz.EQ.Nz) THEN
            ATE(jj,2) = R_ZERO
          ENDIF
          IF (iy.LT.NY) THEN 
            ATE(jj,3) = ar
          ENDIF
          jj = jj + 1
        ENDDO ! iz
      ENDDO ! iy

!     boundary fields

!     left side fields
      iy = 2
      DO iz = 2,Nz
        BTE(iz) = TWO*Cz(iz)/Dy(iy-1)
      ENDDO ! iz

!     surface fields
      iz = 2
      DO iy = 2,Ny
        BTE(Nz+iy) = TWO*Cy(iy)/Dz(iz-1)
      ENDDO ! iy

!     Bottom fields
      iz = Nz
      DO iy = 2,Ny
        BTE(Nz+Ny+iy-1) = TWO*Cy(iy)/Dz(iz)
      ENDDO ! iy

!     Right fields
      iy = Ny
      DO iz = 2,Nz
        BTE(Nz+2*Ny-1+iz) = TWO*Cz(iz)/Dy(iy)
      ENDDO ! iz

      RETURN
      END subroutine! SetupA_TE()


end module
