!   This is almost just another one of those "wrappers"  used to
!    create vectors that represent diagonal matrices.  The code
!    for the actual calulations for the standard staggered grid
!    formulation are done in module GridCalc.   Note that GridCalcS
!    (in the semiglobal branch) would replace this for the spherical
!     coordinate version.   The idea here is to mimic ModEMM matlab
!    approach as much as possible.  At some point it might make sense
!      to clean this up, and get rid of some of the layers of wrappers!

module MetricElements_SG

   use gridcalc
   use vectranslate

   implicit none
   save
   real(kind=prec),dimension(:),pointer  :: FaceA
   real(kind=prec),dimension(:),pointer  :: EdgeL
   real(kind=prec),dimension(:),pointer  :: DualFaceA
   real(kind=prec),dimension(:),pointer  :: DualEdgeL
   real(kind=prec),dimension(:),pointer  :: Vnode
   real(kind=prec),dimension(:),pointer  :: Vedge

Contains
   subroutine setFaceArea(grid)
      type (grid_t), intent(in)           :: grid
      type (rvector)                      :: temp
      call FaceArea(grid,temp)
      call getRvector(temp,FaceA)
      call deall_rvector(temp)
   end subroutine
   !*******************************************************************
   subroutine setEdgeLength(grid)
      type (grid_t), intent(in)           :: grid
      type (rvector)                      :: temp
      call EdgeLength(grid,temp)
      call getRvector(temp,EdgeL)
      call deall_rvector(temp)
   end subroutine
   !*******************************************************************
   subroutine setDualFaceArea(grid)
      type (grid_t), intent(in)           :: grid
      type (rvector)                      :: temp
      call EdgeArea(grid,temp)
      call getRvector(temp,DualFaceA)
      call deall_rvector(temp)
   end subroutine
   !*******************************************************************
   subroutine setDualEdgeLength(grid)
      type (grid_t), intent(in)           :: grid
      type (rvector)                      :: temp
      call FaceLength(grid,temp)
      call getRvector(temp,DualEdgeL)
      call deall_rvector(temp)
   end subroutine
   !*******************************************************************
   subroutine setVnode(grid)
      type (grid_t), intent(in)           :: grid     ! input model
      type (rscalar)                      :: temp
      call NodeVolume(grid,temp)
      call getRscalar(temp,Vnode)
      call deall_rscalar(temp)
   end subroutine
   !*******************************************************************
   subroutine setVedge(grid)
      type (grid_t), intent(in)           :: grid     ! input model
      type (rvector)                      :: temp
      call EdgeVolume(grid,temp)
      call getRvector(temp,Vedge)
      call deall_rvector(temp)
   end subroutine
   !*******************************************************************
   subroutine deall_MetricElements()
      if(associated(FaceA)) then
          deallocate(FaceA)
      endif
      if(associated(DualFaceA)) then
          deallocate(DualFaceA)
      endif
      if(associated(EdgeL)) then
          deallocate(EdgeL)
      endif
      if(associated(DualEdgeL)) then
          deallocate(DualEdgeL)
      endif
      if(associated(Vedge)) then
          deallocate(Vedge)
      endif
      if(associated(Vnode)) then
          deallocate(Vnode)
      endif
   end subroutine deall_MetricElements
end module
