!   This is almost just another one of those "wrappers"  used to
!    create vectors that represent diagonal matrices.  The code
!    for the actual calculations for the standard staggered grid
!    formulation are done in module GridCalc.   Note that GridCalcS
!    (in the semiglobal branch) would replace this for the spherical
!     coordinate version.   The idea here is to mimic ModEMM matlab
!    approach as much as possible. [AK] This wrapper seems to make sense
!    for the efficiency of calculations; the 3D storage of GridCalc
!    is good for model and data mappings but not for the forward solver.
!   Note that V_C is only used externally by ModelMap, so Vcell is not
!    needed since we only define arrays that are used by modelOperator3D

module MetricElements_SG

   use gridcalc ! spherical coords module selected at compile time from GridCalcS.f90
   use vectranslate

   implicit none
   save
   real(kind=prec),dimension(:),pointer  :: FaceA
   real(kind=prec),dimension(:),pointer  :: EdgeL
   real(kind=prec),dimension(:),pointer  :: DualFaceA
   real(kind=prec),dimension(:),pointer  :: DualEdgeL
   real(kind=prec),dimension(:),pointer  :: Vnode
   real(kind=prec),dimension(:),pointer  :: Vedge
   real(kind=prec),dimension(:),pointer  :: Vcell

Contains
   subroutine setFaceArea(grid)
      type (grid_t), intent(in)           :: grid
      type (rvector)                      :: temp
      if (.not. S_F%allocated) then
        call FaceArea(grid,temp)
        call getRvector(temp,FaceA)
        call deall_rvector(temp)
      else
        call getRvector(S_F,FaceA)
      end if
   end subroutine
   !*******************************************************************
   subroutine setEdgeLength(grid)
      type (grid_t), intent(in)           :: grid
      type (rvector)                      :: temp
      if (.not. l_E%allocated) then
        call EdgeLength(grid,temp)
        call getRvector(temp,EdgeL)
        call deall_rvector(temp)
      else
        call getRvector(l_E,EdgeL)
      end if
   end subroutine
   !*******************************************************************
   subroutine setDualFaceArea(grid)
      type (grid_t), intent(in)           :: grid
      type (rvector)                      :: temp
      if (.not. S_E%allocated) then
        call DualFaceArea(grid,temp)
        call getRvector(temp,DualFaceA)
        call deall_rvector(temp)
      else
        call getRvector(S_E,DualFaceA)
      end if
   end subroutine
   !*******************************************************************
   subroutine setDualEdgeLength(grid)
      type (grid_t), intent(in)           :: grid
      type (rvector)                      :: temp
      if (.not. l_F%allocated) then
        call DualEdgeLength(grid,temp)
        call getRvector(temp,DualEdgeL)
        call deall_rvector(temp)
      else
        call getRvector(l_F,DualEdgeL)
      end if
   end subroutine
   !*******************************************************************
   subroutine setVnode(grid)
      type (grid_t), intent(in)           :: grid     ! input model
      type (rscalar)                      :: temp
      if (.not. V_N%allocated) then
        call NodeVolume(grid,temp)
        call getRscalar(temp,Vnode)
        call deall_rscalar(temp)
      else
        call getRscalar(V_N,Vnode)
      end if
   end subroutine
   !*******************************************************************
   subroutine setVedge(grid)
      type (grid_t), intent(in)           :: grid     ! input model
      type (rvector)                      :: temp
      if (.not. V_E%allocated) then
        call EdgeVolume(grid,temp)
        call getRvector(temp,Vedge)
        call deall_rvector(temp)
      else
        call getRvector(V_E,Vedge)
      end if
   end subroutine
   !*******************************************************************
   subroutine setVcell(grid)
      type (grid_t), intent(in)           :: grid     ! input model
      type (rscalar)                      :: temp
      if (.not. V_C%allocated) then
        call CellVolume(grid,temp)
        call getRscalar(temp,Vcell)
        call deall_rscalar(temp)
      else
        call getRscalar(V_C,Vcell)
      end if
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
      if(associated(Vcell)) then
          deallocate(Vcell)
      endif
   end subroutine deall_MetricElements
end module
