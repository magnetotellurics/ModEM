module anisConstraint
  ! this module is used for anisotropic constraint, added by Kong, 2018-12-19
  use math_constants
  ! Damping factor for anisotropy constraint 
  type :: DampingFactorANI
     real (kind=prec)	:: xy
     real (kind=prec)	:: yz
     real (kind=prec)	:: xz
  end type DampingFactorANI
   
contains

  subroutine anisCsNorm(nx,ny,nz,lambdaANI,rhox,rhoy,rhoz,aniNorm)
    implicit none
    integer nx,ny,nz
    type(DampingFactorANI),intent(in) :: lambdaANI
    real(kind=prec),intent(inout) :: aniNorm
    real(kind=prec),intent(in) :: rhox(nx,ny,*),rhoy(nx,ny,*),rhoz(nx,ny,*)
    ! locals
    integer ix,iy,iz
    real(kind=prec) :: delrxy,delryz,delrxz,rxy2,ryz2,rxz2
    
    aniNorm = R_ZERO
    rxy2 = R_ZERO
    ryz2 = R_ZERO
    rxz2 = R_ZERO
    do iz = 1,nz
      do iy = 1,ny
        do ix = 1,nx
          delrxy = rhox(ix,iy,iz) - rhoy(ix,iy,iz)
          delryz = rhoy(ix,iy,iz) - rhoz(ix,iy,iz)
          delrxz = rhox(ix,iy,iz) - rhoz(ix,iy,iz)
          rxy2 = rxy2 + delrxy**2
          ryz2 = ryz2 + delryz**2
          rxz2 = rxz2 + delrxz**2 
        enddo
      enddo
    enddo
    
    aniNorm = lambdaANI%xy*rxy2 + lambdaANI%yz*ryz2 + lambdaANI%xz*rxz2
  
  end subroutine anisCsNorm
  
  
  subroutine anisCsVec(nx,ny,nz,lambdaANI,rhox,rhoy,rhoz,BTmx,BTmy,BTmz)
    implicit none
    integer nx,ny,nz
    type(DampingFactorANI),intent(in) :: lambdaANI
    real(kind=prec),intent(in) :: rhox(nx,ny,*),rhoy(nx,ny,*),rhoz(nx,ny,*)
    real(kind=prec),intent(inout) :: BTmx(nx,ny,*),BTmy(nx,ny,*),BTmz(nx,ny,*)
    ! locals
    integer ix,iy,iz
    
    BTmx(1:nx,1:ny,1:nz) = R_ZERO
    BTmy(1:nx,1:ny,1:nz) = R_ZERO
    BTmz(1:nx,1:ny,1:nz) = R_ZERO
    
    do iz = 1,nz
      do iy = 1,ny
        do ix = 1,nx
          BTmx(ix,iy,iz) = (lambdaANI%xy + lambdaANI%xz) * rhox(ix,iy,iz) &
               - lambdaANI%xy * rhoy(ix,iy,iz) - lambdaANI%xz * rhoz(ix,iy,iz)
          BTmy(ix,iy,iz) = (lambdaANI%xy + lambdaANI%yz) * rhoy(ix,iy,iz) &
               - lambdaANI%xy * rhox(ix,iy,iz) - lambdaANI%yz * rhoz(ix,iy,iz)
          BTmz(ix,iy,iz) = (lambdaANI%xz + lambdaANI%yz) * rhoz(ix,iy,iz) &
               - lambdaANI%xz * rhox(ix,iy,iz) - lambdaANI%yz * rhoy(ix,iy,iz)
        enddo
      enddo
    enddo
  
  end subroutine anisCsVec 
  
  
  
end module anisConstraint