program spOpTestReal

   use spOpTools

   implicit none

   type(spMatCSR_Real)       :: A_CSR,B_CSR,C_CSR,L,U
   type(spMatCSR_CMPLX)       ::  CmplxCSR,L_cmplx,U_cmplx
   type(spMatIJS_Real)       :: A_IJS,B_IJS,C_IJS
   type(spMatIJS_Cmplx)       :: CmplxIJS
   type(spMatCSR_Real),pointer       :: AA(:) 
   type(spMatCSR_Cmplx),pointer       :: CC(:) 

   complex(kind=prec), allocatable, dimension(:)   :: x,y
   real(kind=prec), allocatable, dimension(:)   :: d,col
   integer,         allocatable, dimension(:)   :: idx
   integer              :: fid,nX,nD,nz,m,n,i
   integer, allocatable, dimension(:)   :: r,c

   !    change value of testCase to test different routines in spOpTools
   integer, parameter           :: CONVERT = 1
   integer, parameter           :: MATxVEC = 2
   integer, parameter           :: MATxMAT = 3
   integer, parameter           :: DIAGxMAT = 4
   integer, parameter           :: MATtrans = 5
   integer, parameter           :: LUtriang = 6
   integer, parameter           :: Diag = 7
   integer, parameter           :: LUTsolve = 8
   integer, parameter           :: SubMat = 9
   integer, parameter           :: CSR_2Cd = 10
   integer, parameter           :: CholInc = 11
   integer, parameter           :: Dilu = 12
   integer, parameter           :: DiluCmplx = 13
   integer, parameter           :: BlkDiag = 14
   integer, parameter           :: ilu0 = 15
   integer, parameter           :: sort = 16
   integer      :: testCase = 16

   character*80    cfile

   cfile = 'TEST/xA.bin'
   fid = 95
   open(unit=fid,file=cfile,form ='unformatted')
   read(fid) nX
   allocate(x(nX))
   write(0,*) 'nX = ',nX
   read(fid) x
   close(fid)
   write(0,*) 'x.bin read'
   write(0,*) x
   cfile = 'TEST/D.bin'
   open(unit=fid,file=cfile,form ='unformatted')
   read(fid) nD
   allocate(d(nD))
   read(fid) d
   close(fid)
   write(0,*) 'D.bin read'
   write(0,*) d

   cfile = 'TEST/A_CSR.bin'
   open(unit=fid,file=cfile,form ='unformatted')
   call read_CSR_real(fid,A_CSR)
   close(fid)
   cfile = 'TEST/B_CSR.bin'
   open(unit=fid,file=cfile,form ='unformatted')
   call read_CSR_real(fid,B_CSR)
   close(fid)
   write(0,*) 'A_CSR.bin read'
   write(0,*) 'nRow = ',B_CSR%nRow
   write(0,*) 'nCol = ',B_CSR%nCol
   write(0,*) 'row = ',B_CSR%row
   write(0,*) 'col = ',B_CSR%col
   write(0,*) 'val = ',B_CSR%val

   selectcase(testCase)
      case(CONVERT)
         !    test of conversion: CSR --> IJS
         cfile = 'TEST/A_IJSout.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         m = A_CSR%nRow
         n = A_CSR%nCol
         nz = A_CSR%row(m+1)-1
         call create_spMatIJS_Real(m,n,nz,A_IJS)
         call CSR2IJS_Real(A_CSR,A_IJS)
         call write_IJS_real(fid,A_IJS)
         close(fid)

         !  convert back to IJS
         cfile = 'TEST/A_CSRout.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call create_spMatCSR_Real(m,n,nz,B_CSR)
         call IJS2CSR_Real(A_IJS,B_CSR)
         call write_CSR_real(fid,B_CSR)
         close(fid)
      case(MATxVEC)
         !   A*x 
         allocate(y(A_CSR%nRow))
         call RMATxCVEC(A_CSR,x,y) 
         cfile = 'TEST/yA.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         write(fid) A_CSR%nRow
         write(fid) y
         close(fid)

         !   B*x 
         cfile = 'TEST/xB.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         read(fid) nX
         deallocate(x)
         allocate(x(nX))
         ! write(0,*) 'nX = ',nX
         read(fid) x
         close(fid)
         deallocate(y)
         allocate(y(B_CSR%nRow))
         call RMATxCVEC(B_CSR,x,y) 
         cfile = 'TEST/yB.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         write(fid) B_CSR%nRow
         write(fid) y
         close(fid)

      case(MATxMAT)
         cfile = 'TEST/C_CSRout.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call RMATxRMAT(A_CSR,B_CSR,C_CSR)
         call write_CSR_real(fid,C_CSR)

         cfile = 'TEST/C_IJSout.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         m = C_CSR%nRow
         n = C_CSR%nCol
         nz = C_CSR%row(m+1)-1
         call create_spMatIJS_Real(m,n,nz,C_IJS)
         call CSR2IJS_Real(C_CSR,C_IJS)
         call write_IJS_real(fid,C_IJS)
         close(fid)

      case(DIAGxMAT)
         cfile = 'TEST/Dl_CSRout.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call DIAGxRMAT(d,A_CSR,C_CSR)
         call write_CSR_real(fid,C_CSR)
         close(fid)

         cfile = 'TEST/Dl_IJSout.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         m = C_CSR%nRow
         n = C_CSR%nCol
         nz = C_CSR%row(m+1)-1
         call create_spMatIJS_Real(m,n,nz,C_IJS)
         write(0,*) 'C_IJS created',m,n,nz 
         call CSR2IJS_Real(C_CSR,C_IJS)
         write(0,*) 'converted' 
         call write_IJS_real(fid,C_IJS)
         close(fid)

         cfile = 'TEST/Dr_CSRout.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call RMATxDIAG(A_CSR,d,C_CSR)
         call write_CSR_real(fid,C_CSR)
         close(fid)

         cfile = 'TEST/Dr_IJSout.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call CSR2IJS_Real(C_CSR,C_IJS)
         call write_IJS_real(fid,C_IJS)
         close(fid)

      case(MATtrans)
         cfile = 'TEST/BT_CSRout.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call RMATtrans(B_CSR,C_CSR)
         call write_CSR_real(fid,C_CSR)
         cfile = 'TEST/BT_IJSout.bin'
         close(fid)
         open(unit=fid,file=cfile,form ='unformatted')
         call write_CSRasIJS_real(fid,C_CSR)

      case(LUtriang)
         cfile = 'TEST/A_L.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call lowerTri_Real(A_CSR,L)
         call write_CSRasIJS_real(fid,L)
         cfile = 'TEST/A_U.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call upperTri_real(A_CSR,U)
         call write_CSRasIJS_real(fid,U)

      case(Diag)
         cfile = 'TEST/diagA.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call diag_Real(A_CSR,d)
         write(fid) A_CSR%nRow
         write(fid) d
         close(fid)

      case(LUTsolve)
         cfile = 'TEST/LTsoln.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call lowerTri_Real(A_CSR,L)
         allocate(y(A_CSR%nRow))
         call LTsolve_Real(L,x,y)
         write(fid) A_CSR%nRow
         write(fid) y
         close(fid)
         cfile = 'TEST/UTsoln.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call upperTri_Real(A_CSR,U)
         call UTsolve_Real(U,y,x)
         write(fid) A_CSR%nRow
         write(fid) x
         close(fid)

      case(SubMat)
         allocate(r(3))
         allocate(c(3))
         r(1) = 1
         r(2) = 3
         r(3) = 5
         c = r
         cfile = 'TEST/Asub.bin'
         open(unit=fid,file=cfile,form ='unformatted')
         call subMatrix_Real(A_CSR,r,c,C_CSR)
         call write_CSRasIJS_real(fid,C_CSR)

     case(CSR_2Cd)
         fid = 95
         cfile = 'TEST/Aplus_id.bin'
         open(unit=fid,file=cfile,form ='unformatted')
        call CSR_R2Cdiag(A_CSR,d,CmplxCSR)
        call write_CSRasIJS_cmplx(fid,CmplxCSR)

     case(CholInc)
        fid = 95
        cfile = 'TEST/Asym.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call read_CSR_real(fid,A_CSR)
        close(fid)
        call CholInc_Real(A_CSR,L)
        fid = 95
        cfile = 'TEST/Lcholinc.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call write_CSRasIJS_real(fid,L)

     case(Dilu)
        fid = 95
        cfile = 'TEST/Asym.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call read_CSR_real(fid,A_CSR)
        close(fid)
        call Dilu_Real(A_CSR,L,U)
        cfile = 'TEST/Ldilu.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call write_CSRasIJS_real(fid,L)
        close(fid)
        cfile = 'TEST/Udilu.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call write_CSRasIJS_real(fid,U)

     case(DiluCmplx)
        fid = 95
        cfile = 'TEST/Asym.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call read_CSR_real(fid,A_CSR)
        close(fid)
        call CSR_R2Cdiag(A_CSR,d,CmplxCSR)
        call Dilu_Cmplx(CmplxCSR,L_cmplx,U_cmplx)
        cfile = 'TEST/LdiluC.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call write_CSRasIJS_cmplx(fid,L_cmplx)
        close(fid)
        cfile = 'TEST/UdiluC.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call write_CSRasIJS_cmplx(fid,U_cmplx)

     case(BlkDiag)
        allocate(AA(2))
        AA(1) = A_CSR
        AA(2) = B_CSR
        call BlkDiag_Real(AA,C_CSR)
        cfile = 'TEST/AB_blk_diagReal.bin'
        fid = 95
        open(unit=fid,file=cfile,form ='unformatted')
        call write_CSRasIJS_real(fid,C_CSR)
        close(fid)
        deallocate(AA)
        allocate(CC(2))
        call CSR_R2Cdiag(A_CSR,d,CC(1))
        call R2C_CSR(B_CSR,CC(2)) 
        call BlkDiag_Cmplx(CC,CmplxCSR)
        cfile = 'TEST/AB_blk_diagCmplx.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call write_CSRasIJS_cmplx(fid,CmplxCSR)
        close(fid)
        deallocate(CC)
     case(ilu0)
        fid = 95
        cfile = 'TEST/Aasym.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call read_CSR_real(fid,A_CSR)
        close(fid)
        !cfile = 'TEST/Csym.bin'
        !open(unit=fid,file=cfile,form ='unformatted')
        !call CSR_R2Cdiag(A_CSR,d,CmplxCSR)
        !call write_CSR_Cmplx(fid,CmplxCSR)
        !close(fid)
        call sort_spMatCSR_real(A_CSR)
        call ilu0_Real(A_CSR,L,U)
        cfile = 'TEST/Lilu0.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call write_CSRasIJS_Real(fid,L)
        close(fid)
        cfile = 'TEST/Uilu0.bin'
        open(unit=fid,file=cfile,form ='unformatted')
        call write_CSRasIJS_Real(fid,U)
     case(sort)
        allocate(col(6))
        allocate(idx(6))
        col = (/ 54.0, 2214.0, 6.0, 102.0, 53.0, 55.0/)
        nD = size(col) 
        write(6,*) int(col)
        idx=(/(i, i=1,nD,1)/)
        call QSort(col,idx)
        write(6,*) int(col)
        deallocate(idx)
        deallocate(col)
   end select

end program
