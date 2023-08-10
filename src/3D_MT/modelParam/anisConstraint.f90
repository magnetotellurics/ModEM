module anisConstraint
  ! this module is used for anisotropic constraint, added by Kong, 2018-12-19
  use math_constants
  ! Damping factor for anisotropy constraint 
  type :: DampingFactorANI
     real (kind=prec)	:: xy
  end type DampingFactorANI
   
contains

  subroutine anisCsNorm(nx,ny,nz,lambdaANI,rhox,rhoy,aniNorm)
    implicit none
    integer nx,ny,nz
    real(kind=prec) ,intent(in) :: lambdaANI
    real(kind=prec),intent(inout) :: aniNorm
    real(kind=prec),intent(in) :: rhox(nx,ny,*),rhoy(nx,ny,*)
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
          rxy2 = rxy2 + delrxy**2
        enddo
      enddo
    enddo
    
    aniNorm = lambdaANI*rxy2 
  
  end subroutine anisCsNorm
  
  
  subroutine anisCsVec(nx,ny,nz,lambdaANI,rhox,rhoy,BTmx,BTmy)
    implicit none
    integer nx,ny,nz
    real(kind=prec) ,intent(in) :: lambdaANI
    real(kind=prec),intent(in) :: rhox(nx,ny,*),rhoy(nx,ny,*)
    real(kind=prec),intent(inout) :: BTmx(nx,ny,*),BTmy(nx,ny,*)
    ! locals
    integer ix,iy,iz
    
    BTmx(1:nx,1:ny,1:nz) = R_ZERO
    BTmy(1:nx,1:ny,1:nz) = R_ZERO
    
    do iz = 1,nz
      do iy = 1,ny
        do ix = 1,nx
          BTmx(ix,iy,iz) = (lambdaANI)*rhox(ix,iy,iz) - (lambdaANI)*rhoy(ix,iy,iz) 
          BTmy(ix,iy,iz) = (lambdaANI)*rhoy(ix,iy,iz) - (lambdaANI)*rhox(ix,iy,iz)
        enddo
      enddo
    enddo
  
  
  end subroutine anisCsVec 
  
  
  
end module anisConstraint