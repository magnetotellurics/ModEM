! This module contains addinitional routines between cvector_mg data type
!and sparsevectors. Like dot product etc... which are needed by
!solnspace moduloe for instance
module sg_sparse_vector_mg

use math_constants
use griddef
use sg_vector_mg
use sg_sparse_vector

implicit none


interface dotProd
  module procedure dotProd_noConj_scvector_mg_f
  module procedure dotProd_scvector_mg_f
end interface

public  :: dotProd_noConj_scvector_mg_f, dotProd_scvector_mg_f

Contains

! ***********************************************************************************************

! compute complex dot product between a sparse vector SV and a vector of
! type cvector ... result in c
function dotProd_noConj_scvector_mg_f(sv,v) result(c)

implicit none
type (sparsevecc), intent(in)  :: sv
type (cvector_mg), intent(in)  :: v
complex(kind=prec)  :: c

! local
integer  :: imgrid
complex (kind=prec),allocatable  :: sctemp(:)


  allocate(sctemp(v%mgridSize))

  do imgrid = 1, v%mgridSize
    sctemp(imgrid) = dotProd_noConj_scvector_f(sv,v%cvector_mg(imgrid))
  enddo

  c = sum(sctemp)

  deallocate(sctemp)

end function dotProd_noConj_scvector_mg_f

! **********************************************************************
! compute complex dot product between a sparse vector SV and a vector of
! type cvector_mg ... result in c
function dotProd_scvector_mg_f(sv,v) result(c)

implicit none
type (sparsevecc), intent(in)  :: sv
type (cvector_mg), intent(in)  :: v
complex(kind=prec)  :: c

! local
integer  :: imgrid
complex (kind=prec),allocatable  :: sctemp(:)


  allocate(sctemp(v%mgridSize))

  do imgrid = 1, v%mgridSize
    sctemp(imgrid) = dotProd_scvector_f(sv,v%cvector_mg(imgrid))
  enddo

  c = sum(sctemp)

  deallocate(sctemp)

end function dotProd_scvector_mg_f


end module sg_sparse_vector_mg


