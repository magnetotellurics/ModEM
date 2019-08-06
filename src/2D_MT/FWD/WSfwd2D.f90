module wsfwd2d

!  This module contains the basic forward code routines
!   used for 2D rebocc code of Weerachai Siripunvaraporn
!  Integer parameters formerly defined through include files in
!   the original code are static module variables, which
!   are set by calling a setup routine, and are then used by routines which
!  need to know sizes of passed (or locally allocated) arrays.
!   (Formerly this was in WSfwd2Dpar)
!  This enables dynamic allocation (external to this module)
!    without significant changes to the existing 2D modeling codes
!  This module now contains utility routines used by the original
!    forward modeling codes (WSutils), and the contents of WSfwd2Dmod
!    and WSfwd1Dmod

use math_constants
!use wsLAPACK

implicit none
   save
   !  integer parameters are calculated and set by initialization
   !   routine
   integer              :: NZ3MX, MMIMX, NZ0MX, NY0MX, MMBMX
   integer              :: NZ1MX, NZ2MX
   integer              :: NzSave, NzaSave, NySave

   real(kind=prec),parameter       :: Mue =  MU_0

   !  air conductivity is fixed ... but this can be reset by any
   !  routine that uses this module.
   real(kind=prec)                :: CondAir = SIGMA_AIR

   contains

!**********************************************************************
     subroutine SetWSparams(Ny,Nz,Nza)
     ! setup routine: call this to set grid size

     integer, intent(in) :: Ny,Nz,Nza

     NzSave = Nz
     NzaSave = Nza
     NySave = Ny

     NZ0MX = Nz+1
     NY0MX = Ny+1
     MMIMX = (NY0MX-1)*(NZ0MX-1)
     NZ3MX = 3*NZ0MX+4
     NZ2MX = NZ0MX+2
     NZ1MX = NZ0MX+1
     MMBMX = 2*NY0MX+2*NZ0MX
     end subroutine

!>>>>  Routines from WSUtils:
!**********************************************************************
      SUBROUTINE CopyVectorR8(x_1,x_2,vx,y_1,y_2,vy)
      INTEGER x_1,x_2,y_1,y_2
      real(kind=prec)  vx(*),vy(*)

      INTEGER ix,iy

      IF (y_2-y_1.NE.x_2-x_1) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR COPY VECTOR !!!'
        STOP
      ENDIF

      ix = x_1
      DO iy = y_1,y_2
        vy(iy) = vx(ix)
        ix     = ix + 1
      ENDDO

      RETURN
      END subroutine ! CopyVectorR8

!**********************************************************************
      SUBROUTINE CopyVectorC16(x_1,x_2,vx,y_1,y_2,vy)
      INTEGER x_1,x_2,y_1,y_2
      complex(kind=prec)  vx(*),vy(*)

      INTEGER ix,iy

      IF (y_2-y_1.NE.x_2-x_1) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR COPY VECTOR !!!'
        STOP
      ENDIF

      ix = x_1
      DO iy = y_1,y_2
        vy(iy) = vx(ix)
        ix     = ix + 1
      ENDDO

      RETURN
      END  subroutine! CopyVectorC16

!*********************************************************************

      SUBROUTINE CopyMatrixR8(x00,x01,x10,x11,nx1,nx2,mx,  &
                             y00,y01,y10,y11,ny1,ny2,my)
      INTEGER x00,x01,x10,x11,y00,y01,y10,y11,nx1,nx2,ny1,ny2
      real(kind=prec)  mx(nx1,nx2),my(ny1,ny2)

      INTEGER ix1,ix2,iy1,iy2

      IF ((y01-y00.NE.x01-x00).or.(y11-y10.NE.x11-x10)) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR COPY MATRIX !!!'
        STOP
      ENDIF

      ix1 = x00
      DO iy1 = y00,y01
        ix2 = x10
        DO iy2 = y10,y11
          my(iy1,iy2) = mx(ix1,ix2)
          ix2 = ix2 + 1
        ENDDO
        ix1 = ix1 + 1
      ENDDO

      RETURN
      END subroutine ! CopyMatrixR8

!******************************************************************
      SUBROUTINE TransMatrixToVectorR8(mx,np1,np2,s1,s2,e1,e2,vx,sa,ea)
     INTEGER np1,np2,s1,s2,e1,e2,sa,ea
      real(kind=prec)  mx(np1,np2),vx(*)

      INTEGER i_1,i_2,ia

      IF ((e2-s2+1)*(e1-s1+1).NE.(ea-sa+1)) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR COPY MATRIX !!!'
        STOP
      ENDIF

      ia = sa
      DO i_1 = s1,e1
        DO i_2 = s2,e2
          vx(ia) = mx(i_1,i_2)
          ia = ia + 1
        ENDDO
      ENDDO

      END subroutine !

!***************************************************************
      SUBROUTINE ConstantMatrixR8(mx,np1,np2,n1,n2,const_val)
      INTEGER np1,np2,n1,n2
      real(kind=prec)  mx(np1,np2),const_val

      INTEGER i_1,i_2

      DO i_1 = 1,n1
        DO i_2 = 1,n2
          mx(i_1,i_2) = const_val
        ENDDO
      ENDDO

      RETURN
      END subroutine ! ConstantMatrixR8

      SUBROUTINE ConstantMatrixC16(mx,np1,np2,n1,n2,const_val)
      INTEGER np1,np2,n1,n2
      real(kind=prec)  const_val
      complex(kind=prec) mx(np1,np2)

      INTEGER i_1,i_2

      DO i_1 = 1,n1
        DO i_2 = 1,n2
          mx(i_1,i_2) = const_val
        ENDDO
      ENDDO

      RETURN
      END subroutine ! ConstantMatrixC16

!*******************************************************************
      SUBROUTINE ConstantVectorR8(vx,n,const_val)
      INTEGER n
      real(kind=prec)  vx(*),const_val

      INTEGER i

      DO i = 1,n
        vx(i) = const_val
      ENDDO

      RETURN
      END Subroutine

!******************************************************************
      SUBROUTINE ConstantVectorC16(vx,n,const_val)
      INTEGER n
      real(kind=prec)  const_val
      complex(kind=prec)  vx(*)

      INTEGER i

      DO i = 1,n

        vx(i) = const_val
      ENDDO

      RETURN
      END subroutine ! ConstantVectorC16


!******************************************************************
      SUBROUTINE CumulativeDistance(nx,dx,xdis)
      real(kind=prec)  dx(*),xdis(*)
      INTEGER nx,ix

      xdis(1) = 0.0D+0
      DO  ix = 2,nx+1
       xdis(ix) = xdis(ix-1) + dx(ix-1)
      ENDDO

      RETURN
      End Subroutine

!*******************************************************************
      SUBROUTINE DistanceBetweenBlocks(nx,dx,cx)
      real(kind=prec)  dx(*),cx(*),D2
      INTEGER nx,ix

      !!!!!  IS THIS A BUG ?????  !!!!!!
      D2 = 1.0D+0

      DO  ix = 2,nx
       cx(ix) = dx(ix) + dx(ix-1)
      ENDDO
      cx(1)    = D2*dx(1)
      cx(nx+1) = D2*dx(nx)

      RETURN
      END subroutine ! CumulativeDistance

!<<<<<<<<< End WSutils

!>>>>>>>>  Routines from WSfwd1Dmod
   ! *****************************************************************************
      SUBROUTINE Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)

      INTEGER Nzb
      real(kind=prec)  per,Dzb(*),r1d(*)
      complex(kind=prec) X1D(*)

      INTEGER iz,jj,nzb1,nzb2
      real(kind=prec)  Omega,Omue,skd
      real(kind=prec)  dz1(NZ1MX),cz1(NZ2MX)
      real(kind=prec)  au1d(NZ0MX),ad1d(NZ0MX)
      complex(kind=prec) ac1d(NZ0MX),xb1d(2),xb(NZ0MX)

      nzb1 = Nzb + 1
      nzb2 = Nzb + 2
      Omega = (TWO*PI)/Per
      Omue  = Omega*Mue

      r1d(nzb1) = r1d(Nzb)
      skd  = SQRT(TWO*r1d(nzb1)/Omue)
      call CopyVectorR8(1,Nzb,Dzb,1,Nzb,dz1)
      dz1(nzb1) = skd
      call DistanceBetweenBlocks(nzb1,dz1,cz1)

!     Assign operators
      call ConstantVectorR8(au1d,nzb1,R_ZERO)
      call ConstantVectorR8(ad1d,nzb1,R_ZERO)
      call ConstantVectorC16(ac1d,nzb1,R_ZERO)
      jj = 1
      DO iz = 2,nzb1
        au1d(jj) = TWO*r1d(iz-1)/dz1(iz-1)
        ad1d(jj) = TWO*r1d(iz)/dz1(iz)
        ac1d(jj) = CMPLX(R_ZERO,Omue,kind=prec)*cz1(iz) &
			 - au1d(jj) - ad1d(jj)
        jj = jj + 1
      ENDDO

!     Boundary condition for 1D
      xb1d(1) = ONE
      xb1d(2) = R_ZERO

!     solve Aii*Xi = -Aib*Xb
      call ConstantVectorC16(xb,nzb,R_ZERO)
      xb(1) = -au1d(1)*xb1d(1)

      call Solve1D(nzb,ad1d,ac1d,xb)
      X1D(1) = xb1d(1)
      call CopyVectorC16(1,nzb,xb,2,nzb1,X1D)

100   FORMAT(7e11.3)

      Return
      END Subroutine ! Fwd1d_TM


! *****************************************************************************
      SUBROUTINE Fwd1D_TE(per,Nza,Nz,Dz,s1,X1D)

      INTEGER Nza,Nz
      real(kind=prec)  per,Dz(*),s1(*)
      complex(kind=prec) X1D(*)

      INTEGER iz,jj,nz1,nz2
      real(kind=prec)  Omega,Omue,skd,s1d(NZ1MX)
      real(kind=prec)  dz1(NZ1MX),cz1(NZ2MX)
      real(kind=prec)  au1d(NZ0MX),ad1d(NZ0MX),ss
      complex(kind=prec) ac1d(NZ0MX),xb1d(2),xb(NZ0MX)

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

      skd  = SQRT(TWO/(s1d(nz1)*Omue))
      DO iz = 1,Nz
        dz1(iz) = Dz(iz)
      ENDDO ! iz

!     call CopyVectorR8(1,Nz,Dz,1,Nz,dz1)
      dz1(nz1) = skd
      call DistanceBetweenBlocks(nz1,dz1,cz1)

!     Assign operators
      call ConstantVectorR8(au1d,nz1,R_ZERO)
      call ConstantVectorR8(ad1d,nz1,R_ZERO)
      call ConstantVectorC16(ac1d,nz1,R_ZERO)
      jj = 1
      DO iz = 2,nz1
        au1d(jj) = TWO/dz1(iz-1)
        ad1d(jj) = TWO/dz1(iz)
        ss       = (s1d(iz)*dz1(iz)+s1d(iz-1)*dz1(iz-1))
        ac1d(jj) = CMPLX(R_ZERO,Omue*ss,kind=prec) &
			 - au1d(jj) - ad1d(jj)
        jj = jj + 1
      ENDDO

!     Boundary condition for 1D
      xb1d(1) = ONE
      xb1d(2) = R_ZERO

!     solve Aii*Xi = -Aib*Xb
      call ConstantVectorC16(xb,Nz,R_ZERO)
      xb(1) = -au1d(1)*xb1d(1)

      call Solve1D(Nz,ad1d,ac1d,xb)
      X1D(1) = xb1d(1)
      call CopyVectorC16(1,nz,xb,2,nz1,X1D)


100   FORMAT(7e11.3)

      Return
      END Subroutine ! Fwd1d_TE


! *****************************************************************************
      SUBROUTINE Solve1D(nz0,ad1d,ac1d,xb)

      INTEGER nz0
      real(kind=prec)  ad1d(*)
      complex(kind=prec) ac1d(*),xb(*)

      INTEGER jj,iPiv(NZ0MX),info, k
      complex(kind=prec) a1d(4,NZ0MX)
      character*1 cN

      cN = 'N'
!     make matrix in LAPACK's format
      call ConstantMatrixC16(a1d,4,NZ0MX,4,nz0,R_ZERO)

      DO jj = 2,nz0
        a1d(2,jj) = ad1d(jj-1)
      ENDDO
      DO jj = 1,nz0
        a1d(3,jj) = ac1d(jj)
      ENDDO
      DO jj = 1,nz0-1
        a1d(4,jj) = ad1d(jj)
      ENDDO

      call ZGBTRF(nz0,nz0,1,1,a1d,4,iPiv,info)

      call ZGBTRS(cN,nz0,1,1,1,a1d,4,iPiv,xb,nz0,info)

      END Subroutine ! Solve1D

!<<<<<<<<<<  End of routines from WSfwd1Dmod

!>>>>>>>>>> Routines from WSfwd2Dmod
   ! *****************************************************************************
      SUBROUTINE FormAII(per,Nz0,Ny,AA,AII,iPiv)

      INTEGER Nz0,Ny,iPiv(*)
      real(kind=prec)  per,AA(MMIMX,4)
      complex(kind=prec) AII(NZ3MX,MMIMX)

      INTEGER mmi,jj,nz3,kl,ku,kc,info
      real(kind=prec)  Omega,Omue

      Omega = (TWO*PI)/per
      Omue  = Omega*Mue

!     Forming matrix Aii in LAPACK's format
      nz3 = 3*(Nz0-1) + 1
      mmi = (Ny-1)*(Nz0-1)
      call ConstantMatrixC16(AII,NZ3MX,MMIMX,nz3,mmi,R_ZERO)

      kl = Nz0-1
      ku = Nz0-1
      kc = kl + ku + 1

!     diagonal term
      DO jj = 1,mmi
        AII(kc,jj) = CMPLX(AA(jj,1),Omue*AA(jj,4),kind=prec)
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
      call ZGBTRF(mmi,mmi,kl,ku,AII,NZ3MX,iPiv,info)
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
      real(kind=prec)  Aib(*)
      complex(kind=prec) Xb(*),AXb(*)
      INTEGER mmi,mmb,iy,iz,jj,izz

      mmi = (Ny-1)*(Nz0-1)
      mmb = 2*Ny + 2*Nz0

      call ConstantVectorC16(AXb,mmi,R_ZERO)

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
      real(kind=prec)  per,Dzb(*),Dy(*),CRho(NZ0MX,NY0MX)
      complex(kind=prec) HXI0(*),HXB(*)

      INTEGER nzb1,ny1,mmi,mmb,np1,np2,iz,iy,ii,jj
      real(kind=prec)  r1d(NZ1MX)
      complex(kind=prec) X1D(NZ2MX)
!     is this extra complication needed?
!      real (kind=8), allocatable, dimension(:) ::  r1d
!      complex (kind=8), allocatable, dimension(:)      ::  x1d
!      allocate(r1d(NZ1MX))
!      allocate(x1d(NZ2MX))

      nzb1 = Nzb + 1
      ny1  = Ny  + 1

      mmb = 2*Nzb + 2*Ny
      mmi = (Nzb-1)*(Ny-1)

      call ConstantVectorC16(HXB,mmb,R_ZERO)
      call ConstantVectorC16(HXI0,mmi,R_ZERO)

      np1 = NZ0MX
      np2 = NY0MX

!    left fields
      call TransMatrixToVectorR8(CRho,np1,np2,1,1,Nzb,1,r1d,1,Nzb)
      call Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)
      call CopyVectorC16(1,Nzb+1,X1D,1,Nzb+1,HXB)

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
        call Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)
        HXB(Nzb+Ny-1+iy) = X1D(Nzb+1)
        DO iz = 2,Nzb
          jj = (iy-2)*(Nzb-1) + iz-1
          HXI0(jj) = X1D(iz)
        ENDDO
      ENDDO

!    right fields
      call TransMatrixToVectorR8(CRho,np1,np2,1,Ny,Nzb,Ny,r1d,1,Nzb)
      call Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)
      call CopyVectorC16(1,Nzb+1,X1D,Nzb+2*Ny,2*Nzb+2*Ny,HXB)


100   FORMAT(i5,2e11.3)

      RETURN
      END Subroutine! SetBound2D_TM


! *****************************************************************************
      SUBROUTINE SetBound2D_TE(per,Nza,Nz,Ny,Dz,Dy,Cond2d,EXI0,EXB)
!     This routine modified to take Cond2d as input, instead of CRho
!        to minimize changes to WS code, just copy (part of) Cond2d into CCon

      INTEGER Nza,Nz,Ny
      real(kind=prec)  per,Dz(*),Dy(*)
      complex(kind=prec) EXI0(*),EXB(*)

      INTEGER nz1,ny1,mmi,mmb,np1,np2,iz,iy,ii,jj
      real(kind=prec)  s1d(NZ1MX),CCon(NZ0MX,NY0MX),Cond2D(NZ0MX,NY0MX)
      complex(kind=prec) X1D(NZ2MX)
      nz1 = Nz + 1
      ny1 = Ny + 1

      mmb = 2*Nz + 2*Ny
      mmi = (Nz-1)*(Ny-1)

      DO iy = 1,Ny
        DO iz = 1,Nz-Nza
          CCon(iz,iy) = Cond2D(iz+Nza,iy)
        ENDDO ! iz
      ENDDO ! iy


      call ConstantVectorC16(EXB,mmb,R_ZERO)
      call ConstantVectorC16(EXI0,mmi,R_ZERO)

      np1 = NZ0MX
      np2 = NY0MX

!     left fields
      call TransMatrixToVectorR8(CCon,np1,np2,1,1,Nz-Nza,1,  &
                                s1d,1,Nz-Nza)
      call Fwd1D_TE(per,Nza,Nz,Dz,s1d,X1D)
      call CopyVectorC16(1,Nz+1,X1D,1,Nz+1,EXB)


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
        call Fwd1D_TE(per,Nza,Nz,Dz,s1d,X1D)
        EXB(Nz+Ny-1+iy) = X1D(Nz+1)
        DO iz = 2,Nz
          jj = (iy-2)*(Nz-1) + iz-1
          EXI0(jj) = X1D(iz)
        ENDDO
      ENDDO

!     right fields
      call TransMatrixToVectorR8(CCon,np1,np2,1,Ny,Nz-Nza,Ny,  &
                                s1d,1,Nz-Nza)
      call Fwd1D_TE(per,Nza,Nz,Dz,s1d,X1D)
      call CopyVectorC16(1,Nz+1,X1D,Nz+2*Ny,2*Nz+2*Ny,EXB)

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

      call ConstantMatrixR8(ATM,MMIMX,4,nblx,4,R_ZERO)
      call ConstantVectorR8(BTM,iblx,R_ZERO)

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
      real(kind=prec)  Dz(*),Dy(*),Cz(*),Cy(*)
      real(kind=prec)  ATE(MMIMX,4),BTE(MMBMX)

!     ATE(:,1) : Real diagonal term
!     ATE(:,2) : Real second strip
!     ATE(:,3) : Real thrid strip
!     ATE(:,4) : Imaginary diagonal term
!     BTE      : Boundary

      INTEGER jj,iz,iy,nblx,iblx
      real(kind=prec)  s00,s10,s01,s11,s2d
      real(kind=prec)  ar,al,ad,au,CCon(NZ0MX,NY0MX),COnd2D(NZ0MX,NY0MX)

      nblx = Ny*Nz
      iblx = 2*Ny + 2*Nz

      call ConstantMatrixR8(CCon,NZ0MX,NY0MX,Nza,Ny,CondAir)

      DO iy = 1,Ny
        DO iz = Nza+1,Nz
          CCon(iz,iy) = Cond2D(iz,iy)
        ENddo ! IZ
      ENDDO ! iy

      call ConstantMatrixR8(ATE,MMIMX,4,nblx,4,R_ZERO)
      call ConstantVectorR8(BTE,iblx,R_ZERO)

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

end module wsfwd2d
