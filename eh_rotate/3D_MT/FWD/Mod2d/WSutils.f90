! *****************************************************************************
module WSutils
!  module containing routines from utilfunct.f
  
   use math_constants
   
   contains
  
  ! *****************************************************************************
      SUBROUTINE CopyVectorR8(x_1,x_2,vx,y_1,y_2,vy)
      INTEGER x_1,x_2,y_1,y_2
      REAL(kind=prec)  vx(*),vy(*)

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


! *****************************************************************************
      SUBROUTINE CopyVectorC16(x_1,x_2,vx,y_1,y_2,vy)
      INTEGER x_1,x_2,y_1,y_2
      COMPLEX(kind=prec)  vx(*),vy(*)

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
      

! *****************************************************************************                                                                  
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
      

! *****************************************************************************
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


! *****************************************************************************    
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
      

! *****************************************************************************
      SUBROUTINE ConstantVectorR8(vx,n,const_val)
      INTEGER n
      real(kind=prec)  vx(*),const_val

      INTEGER i
      
      DO i = 1,n
        vx(i) = const_val
      ENDDO

      RETURN
      END Subroutine


! *****************************************************************************
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


! ***************************************************************************** 
      SUBROUTINE CumulativeDistance(nx,dx,xdis)
      real(kind=prec)  dx(*),xdis(*)
      INTEGER nx,ix
      
      xdis(1) = R_ZERO
      DO  ix = 2,nx+1
       xdis(ix) = xdis(ix-1) + dx(ix-1)
      ENDDO 

      RETURN
      End Subroutine


! *****************************************************************************
      SUBROUTINE DistanceBetweenBlocks(nx,dx,cx)
      real(kind=prec)  dx(*),cx(*)
      INTEGER nx,ix
      
      DO  ix = 2,nx
       cx(ix) = dx(ix) + dx(ix-1)
      ENDDO 
      cx(1)    = dx(1)
      cx(nx+1) = dx(nx)

      RETURN
      END subroutine ! CumulativeDistance
      
      
end module
