module cSparseVector3D_SG  
  !
  use Constants
  use cVector3D_SG
  !
  type :: cSparsevector3D_SG_t

    ! complex vector defined on edge/ face nodes;
    ! store the intention of the use in a character string defined
    ! as in GridDef as a parameter: EDGE or FACE
    character(len=80)	                        	:: gridType
    ! nCoeff is number of non-zero nodes
    integer 						:: nCoeff
    ! xyz = 1,2,3 refers to x, y or z components,
    ! i,j,k are arrays of indices that defines grid location
    integer , allocatable, dimension(:) 		:: i,j,k,xyz
    ! c is complex array of coefficients
    complex ( kind=prec ), allocatable, dimension(:) 	:: c
    ! has sparse vector been allocated?
    logical					:: allocated
    ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
    ! (probably will not be needed in the future when compilers will support
    ! ISO/IEC 15581 - the "allocatable array extension")
    logical					:: temporary
    ! pointer to the parent grid not needed or set: we can
    ! make full use of the sparse vector without it...
    ! - should be passed explicitly if it is ever required!
    !type (grid_t), pointer              	:: grid
	!
    end type cSparsevector3D_SG_t
	!
contains
	!
	subroutine full2Sparse( self, cvec )
		!
		type( cSparsevector3D_SG_t ), intent( inout ) :: self
		type( cVector3D_SG_t ), intent( in ) :: cvec
		!
		Integer, allocatable,  dimension(:,:,:)  :: Ix,Jx, Kx,XYZ1
		Integer, allocatable,  dimension(:,:,:)  :: Iy,Jy, Ky,XYZ2
		Integer, allocatable,  dimension(:,:,:)  :: Iz,Jz, Kz,XYZ3
		logical, allocatable,  dimension(:,:,:) ::Mx,My,Mz
		!
		Integer :: i,j,k,Nx, Ny, Nz
		!
		self%gridType=''
		self%nCoeff  = 0
		self%allocated = .false.
		self%temporary = .false.
		!What we need from cvec:
		!cvec%x, cvec%y,cvec%z
		!cvec%gridType

		Ix= cvec%x
		Jx= cvec%x
		Kx= cvec%x

		Iy= cvec%y
		Jy= cvec%y
		Ky= cvec%y

		Iz= cvec%z
		Jz= cvec%z
		Kz= cvec%z

		XYZ1= cvec%x
		XYZ2= cvec%y
		XYZ3= cvec%z



		!X component of the cvec%x		 
		do i=1,size(cvec%x,1)  
		Ix(i,:,:)=i
		end do
		do j=1,size(cvec%x,2)
		Jx(:,j,:)=j
		end do			 
		do k=1,size(cvec%x,3)
		kx(:,:,k)=k
		end do	
		XYZ1=1			 

		!Y component of the cvec%y
		do i=1,size(cvec%y,1)
		Iy(i,:,:)=i
		end do
		do j=1,size(cvec%y,2)
		Jy(:,j,:)=j
		end do			 
		do k=1,size(cvec%y,3)
		ky(:,:,k)=k
		end do	
		XYZ2=2

		!Z component of the cvec%z			 
		do i=1,size(cvec%z,1)
		Iz(i,:,:)=i
		end do
		do j=1,size(cvec%z,2)
		Jz(:,j,:)=j
		end do			 
		do k=1,size(cvec%z,3)
		kz(:,:,k)=k
		end do	
		XYZ3=3
				 

		Mx = cvec%x  /= 0
		My = cvec%y  /= 0
		Mz = cvec%z  /= 0

		! Get indices of Non-Zero coefficients
		self%i=(/ pack(Ix,Mx),pack(Iy,My),pack(Iz,Mz) /)
		self%j=(/ pack(Jx,Mx),pack(Jy,My),pack(Jz,Mz) /)
		self%k=(/ pack(Kx,Mx),pack(Ky,My),pack(Kz,Mz) /)
		! Get Values of Non-Zero coefficients		
		self%c=(/ pack(cvec%x,Mx),pack(cvec%y,My),pack(cvec%z,Mz) /)
		! Get Components
		self%xyz=(/ pack(XYZ1,Mx), pack(XYZ2,My),pack(XYZ3,Mz) /)
		! Get number of Non-Zero coefficients
		self%nCoeff=size(self%c)
		! Get gridType
		self%gridType=cvec%gridType

		Write( *, * ) "i", self%i
		Write( *, * ) "j", self%j
		Write( *, * ) "k", self%k
		Write( *, * ) "xyz", self%xyz
		Write( *, * ) "c", self%c
		Write( *, * ) "nCoeff", self%nCoeff
		!
		deallocate( Ix, Jx, Kx, XYZ1 )
		deallocate( Iy, Jy, Ky, XYZ2 )
		deallocate( Iz, Jz, Kz, XYZ3 )
		deallocate( Mx, My, Mz )
		!
	end subroutine full2Sparse
	!
end module cSparseVector3D_SG  

