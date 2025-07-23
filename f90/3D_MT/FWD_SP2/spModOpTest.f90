program spModOpTest
  !   Test program for sp modelOperator3D module
  use modelOperator3D
    !   most everything needed for these tests is inherited through
    !    this module

        implicit none


! this is used to set up the numerical grid in SensMatrix
  type(grid_t), save            :: grid
!  storage for the "background" conductivity parameter
  type(modelParam_t), save              :: sigma0
!  storage for EM solutions
   real(kind=prec)                      ::inOmega

   integer, parameter           :: TOPOLOGY = 1
   integer, parameter           :: METRICelements = 2
   integer, parameter           :: CurlCurl = 3
   integer, parameter           :: InitDivCor = 4
   integer, parameter           :: InitModelData = 5
   integer, parameter           :: CC_Precond = 6
   integer, parameter           :: SetupDivCor = 7
   integer      :: testCase = 6
   integer     :: fid,nInterior


   character*80    cfile

     !  load model file
     cfile = 'TEST/mTrue.cpr'
     call read_modelParam(grid,sigma0,cfile)
     !   The read routine actually used is defined in the interface
     !   block  of ModelSpace   (defaults to read_modelParam_WS)
     ! Finish setting up the grid (if that is not done in the read
     ! subroutine)
     call setup_grid(grid)
     if(testCase.le.4) then
        call copy_grid(mGrid,grid)
        call setCurlTopology(grid)
        call setGradTopology(grid)
        call setFaceArea(mGrid)
        call setDualFaceArea(mGrid)
        call setEdgeLength(mGrid)
        call setDualEdgeLength(mGrid)
        call setVnode(mGrid)
        call setVedge(mGrid)
        call boundaryIndex(EDGE,grid,EDGEb,EDGEi)
        call boundaryIndex(CORNER,grid,NODEb,NODEi)
        nInterior = size(EDGEi)
        allocate(VomegaMuSig(nInterior))
     endif

   selectcase(testCase)
      case(TOPOLOGY)
         fid = 95
         ! for  initial tests make all of spModelOperator public
         call setCurlTopology(grid)
         call setGradTopology(grid)
         cfile = 'TEST/CURLtopology.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,T)
         close(fid)
         cfile = 'TEST/GRADtopology.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,G)
         close(fid)
         cfile = 'TEST/GRADtopologyCSR.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSR_real(fid,G)
         close(fid)
         call boundaryIndex(EDGE,grid,EDGEb,EDGEi)
         call boundaryIndex(CORNER,grid,NODEb,NODEi)
         cfile = 'TEST/BdryIntIndex.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         write(0,*) '# EDGEi,b',size(EDGEi),size(EDGEb)
         write(0,*) '# NODEi,b',size(NODEi),size(NODEb)
         write(fid) EDGEi
         write(fid) EDGEb
         write(fid) NODEi
         write(fid) NODEb
         close(fid)
      case(METRICelements)
         call setFaceArea(grid)
         call setDualFaceArea(grid)
         call setEdgeLength(grid)
         call setDualEdgeLength(grid)
         call setVnode(grid)
         call setVedge(grid)
         cfile = 'TEST/MetricElements.bin'
         fid = 95
         open(unit=fid,file=cfile,form ='unformatted')
         write(fid) size(FaceA)
         write(fid) FaceA
         write(fid) size(DualFaceA)
         write(fid) DualFaceA
         write(fid) size(EdgeL)
         write(fid) EdgeL
         write(fid) size(DualEdgeL)
         write(fid) DualEdgeL
         write(fid) size(Vedge)
         write(fid) Vedge
         write(fid) size(Vnode)
         write(fid) Vnode
         close(fid)
      case(CurlCurl)
         call CurlCurlSetup()
         fid = 95 
         cfile = 'TEST/CCii.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,CCii)
         close(fid)
         cfile = 'TEST/CCib.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,CCib)
         close(fid)
         call CurlCurlCleanup()
         write(0,*) 'CCii%allocated',CCii%allocated
         write(0,*) 'CCib%allocated',CCib%allocated
      case(CC_Precond)
         call ModelDataInit(grid)
         write(0,*) 'ModelDataInit'
         inOmega = 0.1_dp
         !call updateFreqCond(inOmega,sigma0)
         call updateOmegaMuSig(inOmega,sigma0)
         write(0,*) 'updatec%OmegaMuSig'
         call PC_setup()
         write(0,*) 'done with PC_setup'
         fid = 95
         cfile = 'TEST/CC_L.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_cmplx(fid,L)
         close(fid)
         cfile = 'TEST/CC_U.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_cmplx(fid,U)
         close(fid)
         cfile = 'TEST/CC_LH.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_cmplx(fid,LH)
         close(fid)
         cfile = 'TEST/CC_UH.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_cmplx(fid,UH)
         close(fid)

      case(InitDivCor)
         call divCorInit()
         !  this computes Vdiv, converts G to actual gradient
         fid = 95
         cfile = 'TEST/Vdiv.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,Vdiv)
         close(fid)
         cfile = 'TEST/grad.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,G)
         close(fid)
      case(InitModelData)
         call ModelDataInit(grid)
         fid = 95
         cfile = 'TEST/CCii.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,CCii)
         close(fid)
         cfile = 'TEST/CCib.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,CCib)
         close(fid)
         cfile = 'TEST/Vdiv.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,Vdiv)
         close(fid)
         cfile = 'TEST/grad.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,G)
         close(fid)
         write(0,*) '# VomegaMuSig = ',size(VomegaMuSig)
      case(SetupDivCor)
         !   next complete setup, using conductivity; in this
         !    test program have to call ModelDataInit first
         call ModelDataInit(grid)
         inOmega = 0.1_dp
         call updateOmegaMuSig(inOmega,sigma0)
         call divCorSetup()
         cfile = 'TEST/VDs.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,VDs)
         close(fid)
         cfile = 'TEST/VDsG.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,VDsG)
         close(fid)
         cfile = 'TEST/VDsG_L.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,VDsG_L)
         close(fid)
         cfile = 'TEST/VDsG_U.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,VDsG_U)
         close(fid)
   end select
end program
