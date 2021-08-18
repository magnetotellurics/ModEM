module spOpTools
!  some tools that manipulate sparse matrices in CSR storage
   use math_constants
   use utilities

!    Generic matrix types and tools, using CSR storage, 
!        but with fortran numbering conventions (starting from 1)

     implicit none

     type :: spMatCSR_Real
        integer                                  :: nRow=0
        integer                                  :: nCol=0
        integer, pointer, dimension(:)           :: row,col
        real(kind=prec), pointer, dimension(:)   :: val
        logical            :: allocated = .false.
        logical            :: lower = .false.
        logical            :: upper = .false.
     end type
     type :: spMatCSR_Cmplx
        integer                                  :: nRow=0
        integer                                  :: nCol=0
        integer, pointer, dimension(:)           :: row,col
        complex(kind=prec), pointer, dimension(:)   :: val
        logical            :: allocated = .false.
        logical            :: lower = .false.
        logical            :: upper = .false.
     end type
     type :: spMatIJS_Real
        integer                                  :: nRow=0
        integer                                  :: nCol=0
        integer, pointer, dimension(:)           :: I,J
        real(kind=prec), pointer, dimension(:)   :: S
        logical            :: allocated = .false.
        logical            :: lower = .false.
        logical            :: upper = .false.
     end type
     type :: spMatIJS_Cmplx
        integer                                  :: nRow=0
        integer                                  :: nCol=0
        integer, pointer, dimension(:)           :: I,J
        complex(kind=prec), pointer, dimension(:)   :: S
        logical            :: allocated = .false.
        logical            :: lower = .false.
        logical            :: upper = .false.
     end type

  INTERFACE create_spMatCSR
     module procedure create_spMatCSR_Real
     module procedure create_spMatCSR_Cmplx
  END INTERFACE
  INTERFACE create_spMatIJS
     module procedure create_spMatIJS_Real
     module procedure create_spMatIJS_Cmplx
  END INTERFACE
  INTERFACE deall_spMatCSR
     module procedure deall_spMatCSR_Real
     module procedure deall_spMatCSR_Cmplx
  END INTERFACE
  INTERFACE deall_spMatIJS
     module procedure deall_spMatIJS_Real
     module procedure deall_spMatIJS_Cmplx
  END INTERFACE
  INTERFACE CSR2IJS
     module procedure CSR2IJS_Real
     module procedure CSR2IJS_Cmplx
  END INTERFACE
  INTERFACE IJS2CSR
     module procedure IJS2CSR_Real
     module procedure IJS2CSR_Cmplx
  END INTERFACE
  INTERFACE lowerTri
     module procedure lowerTri_Real
     module procedure lowerTri_Cmplx
  END INTERFACE
  INTERFACE upperTri
     module procedure upperTri_Real
     module procedure upperTri_Cmplx
  END INTERFACE
  INTERFACE LTsolve
     module procedure LTsolve_Real
     module procedure LTsolve_Cmplx
  END INTERFACE
  INTERFACE UTsolve
     module procedure UTsolve_Real
     module procedure UTsolve_Cmplx
  END INTERFACE
  INTERFACE SubMatrix
     module procedure SubMatrix_Real
     module procedure SubMatrix_Cmplx
  END INTERFACE
  INTERFACE splitMAT
     module procedure splitRMAT
     module procedure splitCMAT
  END INTERFACE
  INTERFACE RMATxVEC
     module procedure RMATxCVEC
     module procedure RMATxRVEC
  END INTERFACE
  INTERFACE sort_spMatCSR
     module procedure sort_spMatCSR_real
     module procedure sort_spMatCSR_cmplx
  END INTERFACE

Contains
   subroutine create_spMatCSR_Real(m,n,nz,A)
      integer, intent(in)               :: m,n,nz
      !   A will be sparse m x n with nz non-zero elements
      type(spMatCSR_Real),  intent(inout)         :: A

      if(A%allocated) then
         call deall_spMatCSR(A)
      endif

      A%nRow = m
      A%nCol = n
      allocate(A%row(m+1))
      allocate(A%col(nz))
      allocate(A%val(nz))
      A%row(m+1)=nz+1
      A%allocated = .true.
      return
   end subroutine
!*************************************************************
   subroutine create_spMatCSR_Cmplx(m,n,nz,A)
      integer, intent(in)               :: m,n,nz
      !   A will be sparse m x n with nz non-zero elements
      type(spMatCSR_Cmplx),  intent(inout)         :: A

      if(A%allocated) then
         call deall_spMatCSR(A)
      endif

      A%nRow = m
      A%nCol = n
      allocate(A%row(m+1))
      allocate(A%col(nz))
      allocate(A%val(nz))
      A%row(m+1)=nz+1
      A%allocated = .true.
      return
   end subroutine
!***************************************************************
   subroutine create_spMatIJS_Real(m,n,nz,A)
      integer, intent(in)               :: m,n,nz
      !   A will be sparse m x n with nz non-zero elements
      type(spMatIJS_Real),  intent(inout)         :: A

      if(A%allocated) then
         call deall_spMatIJS(A)
      endif

      A%nRow = m
      A%nCol = n
      allocate(A%I(nz))
      allocate(A%J(nz))
      allocate(A%S(nz))
      A%allocated = .true.
      return
   end subroutine
!***************************************************************
   subroutine create_spMatIJS_Cmplx(m,n,nz,A)
      integer, intent(in)               :: m,n,nz
      !   A will be sparse m x n with nz non-zero elements
      type(spMatIJS_Cmplx),  intent(inout)         :: A

      if(A%allocated) then
         call deall_spMatIJS(A)
      endif

      A%nRow = m
      A%nCol = n
      allocate(A%I(nz))
      allocate(A%J(nz))
      allocate(A%S(nz))
      A%allocated = .true.
      return
   end subroutine
!**********************************************************
   subroutine deall_spMatCSR_Real(A)
 
      type(spMatCSR_Real)         :: A
      if(A%allocated) then
         deallocate(A%row)
         deallocate(A%col)
         deallocate(A%val)
         A%allocated = .false.
         A%upper = .false.
         A%lower = .false.
         A%nRow = 0
         A%nCol = 0
         return
      endif
   end subroutine
!**********************************************************
   subroutine deall_spMatCSR_Cmplx(A)
      type(spMatCSR_Cmplx)         :: A
      if(A%allocated) then
         deallocate(A%row)
         deallocate(A%col)
         deallocate(A%val)
         A%allocated = .false.
         A%upper = .false.
         A%lower = .false.
         A%nRow = 0
         A%nCol = 0
         return
      endif
   end subroutine
!**********************************************************
   subroutine deall_spMatIJS_Real(A)
 
      type(spMatIJS_Real)         :: A
      if(A%allocated) then
         deallocate(A%I)
         deallocate(A%J)
         deallocate(A%S)
         A%allocated = .false.
         A%nRow = 0
         A%nCol = 0
         return
      endif
   end subroutine
!**********************************************************
   subroutine deall_spMatIJS_Cmplx(A)
 
      type(spMatIJS_Cmplx)         :: A
      if(A%allocated) then
         deallocate(A%I)
         deallocate(A%J)
         deallocate(A%S)
         A%allocated = .false.
         A%nRow = 0
         A%nCol = 0
         return
      endif
   end subroutine
!***************************************************************
   logical function sameSizeCSR_Real(A,B)
      type(spMatCSR_Real),  intent(in)         :: A
      type(spMatCSR_Real),  intent(in)         :: B
      sameSizeCSR_Real = .false.
      if(A%allocated .and. B%allocated .and.     &
           A%nRow.gt.0 .and. B%nRow.gt. 0) then 
         sameSizeCSR_Real = A%nRow.eq.B%nRow .and.  &
          A%nCol.eq.B%nCol .and.                 &
          A%row(A%nRow+1).eq.B%row(B%nRow+1)
      endif
      return
   end function
!***************************************************************
   logical function sameSizeCSR_Cmplx(A,B)
      type(spMatCSR_Cmplx),  intent(in)         :: A
      type(spMatCSR_Cmplx),  intent(in)         :: B
      sameSizeCSR_Cmplx = .false.
      if(A%allocated .and. B%allocated .and.     &
           A%nRow.gt.0 .and. B%nRow.gt. 0) then 
         sameSizeCSR_Cmplx = A%nRow.eq.B%nRow .and.  &
          A%nCol.eq.B%nCol .and.                 &
          A%row(A%nRow+1).eq.B%row(B%nRow+1)
      endif
      return
   end function
!*************************************************************
   integer function maxColumnsR(A) 
      type(spMatCSR_Real),  intent(in)         :: A
      integer           :: i
      maxColumnsR = 0
      do i = 1,A%nRow
         maxColumnsR = max(maxColumnsR,A%row(i+1)-A%row(i))
      enddo
      return
   end function
!*************************************************************
   integer function maxColumnsC(A) 
      type(spMatCSR_Cmplx),  intent(in)         :: A
      integer           :: i
      maxColumnsC = 0
      do i = 1,A%nRow
         maxColumnsC = max(maxColumnsC,A%row(i+1)-A%row(i))
      enddo
      return
   end function
!*************************************************************
   subroutine CSR2IJS_Real(C,S)
      type(spMatCSR_Real),intent(in)            ::  C
      type(spMatIJS_Real),intent(inout)         ::  S
      integer           :: ij,i,j
      ! for now no error checking
      if(.not.S%allocated) then
         call errStop('CSR2IJS: allocate output matrix before call')
      endif
      ij = 0
      do i=1,C%nRow
         do j = C%row(i),C%row(i+1)-1
            ij = ij + 1
            S%I(ij) = i
            S%J(ij) = C%col(j)
            S%S(ij) = C%val(j)
         enddo
      enddo
      return
   end subroutine
!*************************************************************
   subroutine CSR2IJS_Cmplx(C,S)
      type(spMatCSR_Cmplx),intent(in)            ::  C
      type(spMatIJS_Cmplx),intent(inout)         ::  S
      integer           :: ij,i,j
      ! for now no error checking
      if(.not.S%allocated) then
         call errStop('CSR2IJS: allocate output matrix before call')
      endif
      ij = 0
      do i=1,C%nRow
         do j = C%row(i),C%row(i+1)-1
            ij = ij + 1
            S%I(ij) = i
            S%J(ij) = C%col(j)
            S%S(ij) = C%val(j)
         enddo
      enddo
      return
   end subroutine
!*************************************************************
   subroutine IJS2CSR_Real(S,C)
      type(spMatIJS_Real),intent(in)           ::  S
      type(spMatCSR_Real),intent(inout)        ::  C
      integer           :: i,j,nz
      integer, allocatable, dimension(:)  :: rowT

      ! limited error checking

      allocate(rowT(S%nRow+1))

      if(.not.C%Allocated) then
         call errStop('IJS2CSR: allocate output matrix before call')
      endif

      !   first pass: find numbers of columns in each row of output
      rowT = 0
      nz = size(S%I)
      do i = 1,nz
         rowT(S%I(i)) = rowT(S%I(i))+1
      enddo
      !   set row array in output CSR matrix
      C%row(1) = 1
      do i = 1,C%nRow
         C%row(i+1) = C%row(i)+rowT(i)
      enddo

      !    now fill in columns and values
      rowT = 0
      do i = 1,nz
         j = C%row(S%I(i)) +rowT(S%I(i))
         C%col(j) = S%J(i)
         C%val(j) = S%S(i) 
         rowT(S%I(i)) = rowT(S%I(i))+1
      enddo
      deallocate(rowT)
      return
   end subroutine
!*****************************************************************
   subroutine IJS2CSR_Cmplx(S,C)
      type(spMatIJS_Cmplx),intent(in)           ::  S
      type(spMatCSR_Cmplx),intent(inout)        ::  C
      integer           :: i,j,nz
      integer, allocatable, dimension(:)  :: rowT
      ! for now no error checking
      !    first should sort sparse matrix in IJS format, so that
      !    row numbers are strictly non-decreasing; come back to
      !    this, since spMatIJS converted from spMatCSR will
      !    already be ordered

      allocate(rowT(S%nRow+1))

      if(.not.C%Allocated) then
         call errStop('IJS2CSR: allocate output matrix before call')
      endif

      !   first pass: find numbers of columns in each row of output
      rowT = 0
      nz = size(S%I)
      do i = 1,nz
         rowT(S%I(i)) = rowT(S%I(i))+1
      enddo
      !   set row array in output CSR matrix
      C%row(1) = 1
      do i = 1,C%nRow
         C%row(i+1) = C%row(i)+rowT(i)
      enddo

      !    now fill in columns and values
      rowT = 0
      do i = 1,nz
         j = C%row(S%I(i)) +rowT(S%I(i))
         C%col(j) = S%J(i)
         C%val(j) = S%S(i) 
         rowT(S%I(i)) = rowT(S%I(i))+1
      enddo
      deallocate(rowT)
      return
   end subroutine
!*************************************************************
   subroutine RMATxCVEC(A,x,y)
      ! multiply a complex vector x by a real sparse CSR matrix A 
      type(spMatCSR_Real),  intent(in)                  :: A
      complex(kind=prec), dimension(:), intent(in)      :: x
      complex(kind=prec), dimension(:), intent(inout)   :: y

      integer     :: i,j
      
      ! lets start coding this with little checking -- assume
      ! everything is allocated and correct on entry
      
      if(A%nCol.ne.size(x)) then
         write(0,*) ' Error: A%nCol=',A%nCol,' size(x)=',size(x)
         call errStop('RMATxCVEC: matrix and vector sizes incompatible')
      endif

      do i = 1,A%nRow
         y(i) = C_ZERO
         do j = A%row(i),A%row(i+1)-1 
            y(i) = y(i)+A%val(j)*x(A%col(j))
         enddo
      enddo
      return
   end subroutine
!*************************************************************
   subroutine RMATxRVEC(A,x,y)
      ! multiply a complex vector x by a real sparse CSR matrix A 
      type(spMatCSR_Real),  intent(in)                  :: A
      real(kind=prec), dimension(:), intent(in)         :: x
      real(kind=prec), dimension(:), intent(inout)      :: y

      integer     :: i,j
      
      ! lets start coding this with little checking -- assume
      ! everything is allocated and correct on entry
      
      if(A%nCol.ne.size(x)) then
         call errStop('RMATxRVEC: matrix and vector sizes incompatible')
      endif

      do i = 1,A%nRow
         y(i) = 0.0
         do j = A%row(i),A%row(i+1)-1 
            y(i) = y(i)+A%val(j)*x(A%col(j))
         enddo
      enddo
      return
   end subroutine
!*******************************************************************
   subroutine CMATxCVEC(A,x,y)
      ! multiply a complex vector x by a complex sparse CSR matrix A 
      type(spMatCSR_Cmplx),  intent(in)         :: A
      complex(kind=prec), dimension(:), intent(in)   :: x
      complex(kind=prec), dimension(:), intent(inout)   :: y

      integer     :: i,j
      
      ! lets start coding this with little error checking -- assume
      ! everything is allocated and correct on entry

      if(A%nCol.ne.size(x)) then
         call errStop('CMATxCVEC: matrix and vector sizes incompatible')
      endif

      do i = 1,A%nRow
         y(i) = 0.0
         do j = A%row(i),A%row(i+1)-1 
            y(i) = y(i)+A%val(j)*x(A%col(j))
         enddo
      enddo
      return
   end subroutine
!*******************************************************************
   subroutine RMATxRMAT(A,B,C)
     !  matix-vector multiplication
      type(spMatCSR_Real),  intent(in)         :: A
      type(spMatCSR_Real),  intent(in)         :: B
      type(spMatCSR_Real),  intent(inout)      :: C
      type(spMatCSR_Real)                      :: Ctmp

      integer   :: i,j,k,nColMax,nCol,m,n,nz,nnz,jj,j1,j2,i1,i2,l,nzero
      integer, allocatable, dimension(:)       :: colT
      integer, allocatable, dimension(:)       :: rowT
      logical  new

      allocate(rowT(A%nRow+1))
      nColMax = maxColumnsR(A)*maxColumnsR(B)
      allocate(colT(nColMax))

      if(A%nCol.ne.B%nRow) then
         call errStop('RMATxRMAT: matrix sizes incompatible')
      endif

      !   first pass: find numbers of columns in each row of output
      !            matrix C
      rowT(1) = 1
      do i = 1,A%nRow
         nCol = 0
         do j = A%row(i),A%row(i+1)-1
            jj = A%col(j)
            do k = B%row(jj),B%row(jj+1)-1 
               new = .true.
               do l = 1,nCol
                  new = new.and.(colT(l).ne.B%col(k))
               enddo
               if(new) then
                  nCol = nCol+1
                  colT(nCol) = B%col(k)
               endif
            enddo
         enddo
         rowT(i+1) = rowT(i)+nCol 
      enddo
      !  create output sparse matrix
      !  rowT is now row vector for output C
      m = A%nRow
      n = B%nCol
      nz = rowT(m+1)-1
      if(C%allocated) then
          call deall_spMatCSR(C)
      endif
      call create_spMatCSR(m,n,nz,Ctmp)
      !   second pass: fill in columns and values of output
      !            matrix C

      nzero = 0
      Ctmp%row = rowT
      deallocate(rowT)
      do i = 1,A%nRow
         nCol = 0
         i1 = Ctmp%row(i)
         i2 = Ctmp%row(i+1)-1
         Ctmp%val(i1:i2) = 0.0d0
         do j = A%row(i),A%row(i+1)-1
            jj = A%col(j)
            do k = B%row(jj),B%row(jj+1)-1 
               new = .true.
               do l = 1,nCol
                  if(colT(l).eq.B%col(k))then
                     new = .false.
                     exit
                  endif
               enddo
               if(new) then
                  nCol = nCol+1
                  colT(nCol) = B%col(k)
                  Ctmp%val(i1+nCol-1) = A%val(j)*B%val(k)
                  Ctmp%col(i1+nCol-1) = B%col(k)
                  if(Ctmp%val(i1+nCol-1).eq.0.0) then ! new entry
                      nzero = nzero +1
                  endif
               else
                  Ctmp%val(i1+l-1) = Ctmp%val(i1+l-1) + A%val(j)*B%val(k)
                  if(Ctmp%val(i1+l-1).eq.0.0) then ! no new entry 
                      ! but could cancel out nonetheless
                      nzero = nzero +1
                  endif
               endif
            enddo
         enddo
      enddo
      deallocate(colT)
      ! now try to clear zeros in Btmp
      nz = nz - nzero
      call create_spMatCSR(m,n,nz,C)  
      if(nzero.eq.0) then
          C%row=Ctmp%row
          C%col=Ctmp%col
          C%val=Ctmp%val
      else 
          j1 = 1
          k  = 0
          nz = 1
          do i = 1,m
             nnz = 0
             j2 = Ctmp%row(i+1)-1
             do j = j1,j2
                if(abs(Ctmp%val(j)).gt. 0) then
                   nnz = nnz + 1
                   k = k+1
                   C%col(k) = Ctmp%col(j)
                   C%val(k) = Ctmp%val(j)
                endif
             enddo
             j1 = j2+1
             C%row(i) = nz
             nz = nz + nnz
          enddo
          C%row(m+1) = nz
      endif
      call deall_spMatCSR_Real(Ctmp)
      return
   end subroutine
!*******************************************************************
   subroutine CMATxCMAT(A,B,C)
     !  matix-vector multiplication, complex version
      type(spMatCSR_Cmplx),  intent(in)         :: A
      type(spMatCSR_Cmplx),  intent(in)         :: B
      type(spMatCSR_Cmplx),  intent(inout)      :: C

      integer   :: i,j,k,nColMax,nCol,m,n,nz,jj,i1,i2,l
      integer, allocatable, dimension(:)       :: colT
      integer, allocatable, dimension(:)  :: rowT
      logical new

      if(A%nCol.ne.B%nRow) then
         call errStop('CMATxCMAT: matrix sizes incompatible')
      endif

      allocate(rowT(A%nRow+1))
      nColMax = maxColumnsC(A)*maxColumnsC(B)
      allocate(colT(nColMax))

      !   first pass: find numbers of columns in each row of output
      !            matrix C
      rowT(1) = 1
      do i = 1,A%nRow
         nCol = 0
         do j = A%row(i),A%row(i+1)-1
            jj = A%col(j)
            do k = B%row(jj),B%row(jj+1)-1 
               new = .true.
               do l = 1,nCol
                  new = new.and.(colT(l).ne.B%col(k))
               enddo
               if(new) then
                  nCol = nCol+1
                  colT(nCol) = B%col(k)
               endif
            enddo
         enddo
         rowT(i+1) = rowT(i)+nCol 
      enddo

      !  create output sparse matrix
      !  rowT is now row vector for output C
      m = A%nRow
      n = B%nCol
      nz = rowT(m+1)-1
      call create_spMatCSR(m,n,nz,C)

      !   second pass: fill in columns and values of output
      !            matrix C

      C%row = rowT
      deallocate(rowT)
      do i = 1,A%nRow
         nCol = 0
         i1 = C%row(i)
         i2 = C%row(i+1)-1
         C%val(i1:i2) = 0.0d0
         do j = A%row(i),A%row(i+1)-1
            jj = A%col(j)
            do k = B%row(jj),B%row(jj+1)-1 
               new = .true.
               do l = 1,nCol
                  if(colT(l).eq.B%col(k))then
                     new = .false.
                     exit
                  endif
               enddo
               if(new) then
                  nCol = nCol+1
                  colT(nCol) = B%col(k)
                  C%val(i1+nCol-1) = A%val(j)*B%val(k)
                  C%col(i1+nCol-1) = B%col(k)
               else
                  C%val(i1+l-1) = C%val(i1+l-1) + A%val(j)*B%val(k)
               endif
            enddo
         enddo
      enddo
      deallocate(colT)
      return
   end subroutine
!*******************************************************************
   subroutine DIAGxRMAT(D,A,B)
     !  premultiply sparse matrix A by diagonal matrix D
     !        real version
      type(spMatCSR_Real),  intent(in)           :: A
      real(kind=prec), intent(in), dimension(:)  :: D
      type(spMatCSR_Real),  intent(inout)        :: B
      type(spMatCSR_Real)                        :: Btmp

      integer   :: i,j,j1,j2,k,n,m,nz,nnz,nzero

      if(A%nRow.ne.size(D)) then
         call errStop('DIAGxRMAT: matrix sizes incompatible')
      endif
      m = A%nRow
      n = A%nCol
      nz = A%row(A%nRow+1)-1
      if(.not.sameSizeCSR_Real(A,B)) then
         if(B%allocated) then
            call deall_spMatCSR_Real(B)
         endif
         call create_spMatCSR(m,n,nz,B)  
      endif
      call create_spMatCSR(m,n,nz,Btmp)  
      nzero=0
      Btmp%row(1) = 1
      do i = 1,A%nRow 
         Btmp%row(i+1) = A%row(i+1)
         do j = A%row(i),A%row(i+1)-1
            Btmp%val(j) = D(i)*A%val(j)
            Btmp%col(j) = A%col(j)
            if(Btmp%val(j).eq.0.0) then 
                nzero=nzero+1
            endif
         enddo
      enddo
      ! now try to clear zeros in Btmp
      nz = nz - nzero
      if(nzero.eq.0) then
          B%row=Btmp%row
          B%col=Btmp%col
          B%val=Btmp%val
      else 
          j1 = 1
          k  = 0
          nz = 1
          do i = 1,m
             nnz = 0
             j2 = Btmp%row(i+1)-1
             do j = j1,j2
                if(abs(Btmp%val(j)).gt. 0) then
                   nnz = nnz + 1
                   k = k+1
                   B%col(k) = Btmp%col(j)
                   B%val(k) = Btmp%val(j)
                endif
             enddo
             j1 = j2+1
             B%row(i) = nz
             nz = nz + nnz
          enddo
          B%row(m+1) = nz
      endif
      call deall_spMatCSR_Real(Btmp)
      return
   end subroutine
!*******************************************************************
   subroutine DIAGxCMAT(D,A,B)
     !  premultiply sparse matrix A by diagonal matrix D
     !        complex version
      type(spMatCSR_Cmplx),  intent(in)             :: A
      complex(kind=prec), intent(in), dimension(:)  :: D
      type(spMatCSR_Cmplx),  intent(inout)          :: B

      integer   :: i,j,m,n,nz

      if(A%nRow.ne.size(D)) then
         call errStop('DIAGxCMAT: matrix sizes incompatible')
      endif
      if(.not.sameSizeCSR_Cmplx(A,B)) then
         if(B%allocated) then
            call deall_spMatCSR(B)
         endif
         m = A%nRow
         n = A%nCol
         nz = A%row(A%nRow+1)-1
         call create_spMatCSR(m,n,nz,B)  
      endif

      B%row(1) = 1
      do i = 1,A%nRow 
         B%row(i+1) = A%row(i+1)
         do j = A%row(i),A%row(i+1)-1
            B%val(j) = D(i)*A%val(j)
            B%col(j) = A%col(j)
         enddo
      enddo
      return
   end subroutine
!*******************************************************************
   subroutine RMATxDIAG(A,D,B)
     !  postmultiply sparse matrix A by diagonal matrix D
     !        real version
      type(spMatCSR_Real),  intent(in)           :: A
      real(kind=prec), intent(in), dimension(:)  :: D
      type(spMatCSR_Real),  intent(inout)        :: B
      type(spMatCSR_Real)                        :: Btmp

      integer   :: i,j,j1,j2,k,m,n,nz,nnz,nzero

      if(A%nCol.ne.size(D)) then
         call errStop('RMATxDIAG: matrix sizes incompatible')
      endif
      m = A%nRow
      n = A%nCol
      nz = A%row(A%nRow+1)-1
      if(.not.sameSizeCSR_Real(A,B)) then
         if(B%allocated) then
            call deall_spMatCSR_Real(B)
         endif
         call create_spMatCSR(m,n,nz,B)  
      endif
      call create_spMatCSR(m,n,nz,Btmp)  
      nzero = 0
      Btmp%row(1) = 1
      do i = 1,A%nRow 
         Btmp%row(i+1) = A%row(i+1)
         do j = A%row(i),A%row(i+1)-1
            Btmp%val(j) = A%val(j)*D(A%col(j))
            Btmp%col(j) = A%col(j)
            if(Btmp%val(j).eq.0.0) then ! mark zero elements
                nzero = nzero + 1
            endif
         enddo
      enddo
      ! now try to clear zeros in Btmp
      nz = nz - nzero
      call create_spMatCSR(m,n,nz,B)  
      if(nzero.eq.0) then
          B%row=Btmp%row
          B%col=Btmp%col
          B%val=Btmp%val
      else 
          j1 = 1
          k  = 0
          nz = 1
          do i = 1,m
             nnz = 0
             j2 = Btmp%row(i+1)-1
             do j = j1,j2
                if(abs(Btmp%val(j)).gt. 0) then
                   nnz = nnz + 1
                   k = k+1
                   B%col(k) = Btmp%col(j)
                   B%val(k) = Btmp%val(j)
                endif
             enddo
             j1 = j2+1
             B%row(i) = nz
             nz = nz + nnz
          enddo
          B%row(m+1) = nz
      endif
      call deall_spMatCSR_Real(Btmp)
      return
   end subroutine
!*******************************************************************
   subroutine CMATxDIAG(A,D,B)
     !  postmultiply sparse matrix A by diagonal matrix D
     !        complex version
      type(spMatCSR_Cmplx),  intent(in)          :: A
      real(kind=prec), intent(in), dimension(:)  :: D
      type(spMatCSR_Cmplx),  intent(inout)       :: B

      integer   :: i,j,m,n,nz

      if(A%nCol.ne.size(D)) then
         call errStop('CMATxDIAG: matrix sizes incompatible')
      endif
      if(.not.sameSizeCSR_Cmplx(A,B)) then
         if(B%allocated) then
            call deall_spMatCSR(B)
         endif
         m = A%nRow
         n = A%nCol
         nz = A%row(A%nRow+1)-1
         call create_spMatCSR(m,n,nz,B)  
      endif

      B%row(1) = 1
      do i = 1,A%nRow 
         B%row(i+1) = A%row(i+1)
         do j = A%row(i),A%row(i+1)-1
            B%val(j) = A%val(j)*D(A%col(j))
            B%col(j) = A%col(j)
         enddo
      enddo
      return
   end subroutine
!*******************************************************************
   subroutine CMATtrans(A,Atrans,Conj)
      type(spMatCSR_Cmplx),  intent(in)          :: A
      type(spMatCSR_Cmplx),  intent(inout)       :: Atrans
      logical, intent(in),optional               :: Conj
      type(spMatIJS_Cmplx)                       :: B
      integer                   :: i,nz,temp
      logical                   :: conjugate = .true.

      if(present(Conj)) then 
         conjugate = Conj
      endif
      nz = A%row(A%nRow+1)-1
      call create_spMatCSR(A%nCol,A%nRow,nz,Atrans)
      call create_spMatIJS(A%nRow,A%nCol,nz,B)
      call CSR2IJS(A,B)
      do i = 1,nz
          temp = B%I(i)
          B%I(i) = B%J(i)
          B%J(i) = temp
      enddo
      if(conjugate) then
         B%S = conjg(B%S)
      endif
      temp = B%nRow
      B%nRow = B%nCol
      B%nCol = temp

      call IJS2CSR(B,Atrans)
      call deall_spMATIJS(B)
      if(A%lower) then
         Atrans%upper = .true.
      endif 
      if(A%upper) then
         Atrans%lower = .true.
      endif 
      return
   end subroutine
!*******************************************************************
   subroutine RMATtrans(A,Atrans)
      type(spMatCSR_Real),  intent(in)          :: A
      type(spMatCSR_Real),  intent(inout)       :: Atrans
      type(spMatIJS_Real)                       :: B
      integer                   :: i,nz,nz1,temp

      nz = A%row(A%nRow+1)-1
      call create_spMatCSR(A%nCol,A%nRow,nz,Atrans)

      call create_spMatIJS(A%nRow,A%nCol,nz,B)

      call CSR2IJS(A,B)
      do i = 1,nz
          temp = B%I(i)
          B%I(i) = B%J(i)
          B%J(i) = temp
      enddo
      temp = B%nRow
      B%nRow = B%nCol
      B%nCol = temp

      call IJS2CSR_Real(B,Atrans)
      call deall_spMATIJS(B)

      if(A%lower) then
         Atrans%upper = .true.
      endif 
      if(A%upper) then
         Atrans%lower = .true.
      endif 
      return
   end subroutine
!********************************************************************
   subroutine write_CSR_real(fid,A)
      integer, intent(in)          :: fid
      type(spMatCSR_Real),intent(in) :: A
      write(fid) A%nRow,A%nCol,A%row(A%nRow+1)-1
      write(fid) A%row
      write(fid) A%col
      write(fid) A%val
      return
   end subroutine
!********************************************************************
   subroutine read_CSR_real(fid,A)
      integer, intent(in)          :: fid
      type(spMatCSR_Real),intent(inout) :: A
      integer    :: n,m,nz
      read(fid) n,m,nz
      if(A%allocated) then
         call deall_spMatCSR(A)
      endif
      call create_spMatCSR(n,m,nz,A)
      read(fid) A%row
      read(fid) A%col
      read(fid) A%val
      return
   end subroutine
!********************************************************************
   subroutine write_CSR_Cmplx(fid,A)
      integer, intent(in)          :: fid
      type(spMatCSR_Cmplx),intent(in) :: A
      write(fid) A%nRow,A%nCol,A%row(A%nRow+1)-1
      write(fid) A%row
      write(fid) A%col
      write(fid) A%val
      return
   end subroutine
!********************************************************************
   subroutine read_CSR_Cmplx(fid,A)
      integer, intent(in)          :: fid
      type(spMatCSR_Cmplx),intent(inout) :: A
      integer    :: n,m,nz
      read(fid) n,m,nz
      if(A%allocated) then
         call deall_spMatCSR(A)
      endif
      call create_spMatCSR(n,m,nz,A)
      read(fid) A%row
      read(fid) A%col
      read(fid) A%val
      return
   end subroutine
!********************************************************************
   subroutine write_IJS_real(fid,A)
      integer, intent(in)          :: fid
      type(spMatIJS_Real),intent(in) :: A
      integer    :: nz
 
      nz = size(A%I)
      write(fid) A%nRow,A%nCol,nz
      write(fid) A%I
      write(fid) A%J
      write(fid) A%S
      return
   end subroutine
!********************************************************************
   subroutine write_IJS_Cmplx(fid,A)
      integer, intent(in)          :: fid
      type(spMatIJS_Cmplx),intent(in) :: A
      integer    :: nz
 
      nz = size(A%I)
      write(fid) A%nRow,A%nCol,nz
      write(fid) A%I
      write(fid) A%J
      write(fid) A%S
      return
   end subroutine
!********************************************************************
   subroutine read_IJS_real(fid,A)
      integer, intent(in)          :: fid
      type(spMatIJS_Real),intent(inout) :: A
      integer    :: n,m,nz
      read(fid) n,m,nz
      if(A%allocated) then
         call deall_spMatIJS(A)
      endif
      call create_spMatIJS(n,m,nz,A)
      read(fid) A%I
      read(fid) A%J
      read(fid) A%S
      return
   end subroutine
!********************************************************************
   subroutine read_IJS_Cmplx(fid,A)
      integer, intent(in)          :: fid
      type(spMatIJS_Cmplx),intent(inout) :: A
      integer    :: n,m,nz
      read(fid) n,m,nz
      if(A%allocated) then
         call deall_spMatIJS(A)
      endif
      call create_spMatIJS(n,m,nz,A)
      read(fid) A%I
      read(fid) A%J
      read(fid) A%S
      return
   end subroutine
!*****************************************************************
   subroutine write_CSRasIJS_Real(fid,A)
      integer, intent(in)          :: fid
      type(spMatCSR_Real),intent(in) :: A

      type(spMatIJS_Real)            :: B
      integer    :: n,m,nz

      m = A%nRow
      n = A%nCol
      nz = A%row(m+1)-1
      call create_spMatIJS(m,n,nz,B)
      call CSR2IJS(A,B)
      call write_IJS_real(fid,B)
      close(fid)
      call deall_spMatIJS(B)
      return
   end subroutine
!*****************************************************************
   subroutine write_CSRasIJS_Cmplx(fid,A)
      integer, intent(in)          :: fid
      type(spMatCSR_Cmplx),intent(in) :: A

      type(spMatIJS_Cmplx)            :: B
      integer    :: n,m,nz

      m = A%nRow
      n = A%nCol
      nz = A%row(m+1)-1
      call create_spMatIJS(m,n,nz,B)
      call CSR2IJS(A,B)
      call write_IJS_Cmplx(fid,B)
      close(fid)
      call deall_spMatIJS(B)
      return
   end subroutine
!****************************************************************
   subroutine lowerTri_Real(A,L)
   !   extract lower triangular part of matrix A in CSR storage
      type(spMatCSR_Real),intent(in)        ::  A
      type(spMatCSR_Real),intent(inout)     ::  L
      integer           :: kk,n,m,i,j,nz
      integer, allocatable, dimension(:)  :: rowT

      m = A%nRow
      n = A%nCol
      allocate(rowT(m))
 
      !   first pass: find numbers of columns in each row of output
      rowT = 0
      do i = 1,m
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).le.i) then
               rowT(i) = rowT(i)+1
            endif
         enddo
      enddo
      nz = 0
      do i = 1,m
         nz = nz+rowT(i)
      enddo
      call create_spMatCSR(m,n,nz,L)
      L%lower = .true.

      !   set row array in output CSR matrix
      L%row(1) = 1
      do i = 1,L%nRow
         L%row(i+1) = L%row(i)+rowT(i)
      enddo
      deallocate(rowT)

      !    now fill in columns and values
      do i = 1,m
         kk = 0
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).le.i) then
               L%col(L%row(i)+kk) = A%col(j) 
               L%val(L%row(i)+kk) = A%val(j) 
               kk = kk+1
            endif
         enddo
      enddo
      return
   end subroutine
!*******************************************************************
   subroutine upperTri_Real(A,U)
   !   extract lower triangular part of matrix A in CSR storage
      type(spMatCSR_Real),intent(in)        ::  A
      type(spMatCSR_Real),intent(inout)     ::  U
      integer           :: kk,n,m,i,j,nz
      integer, allocatable, dimension(:)  :: rowT

      m = A%nRow
      n = A%nCol
      allocate(rowT(m))
 
      !   first pass: find numbers of columns in each row of output
      rowT = 0
      do i = 1,m
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).ge.i) then
               rowT(i) = rowT(i)+1
            endif
         enddo
      enddo
      nz = 0
      do i = 1,m
         nz = nz+rowT(i)
      enddo
      call create_spMatCSR(m,n,nz,U)
      U%upper = .true.

      !   set row array in output CSR matrix
      U%row(1) = 1
      do i = 1,U%nRow
         U%row(i+1) = U%row(i)+rowT(i)
      enddo
      deallocate(rowT)

      !    now fill in columns and values
      do i = 1,m
         kk = 0
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).ge.i) then
               U%col(U%row(i)+kk) = A%col(j) 
               U%val(U%row(i)+kk) = A%val(j) 
               kk = kk+1
            endif
         enddo
      enddo
      return
   end subroutine
!****************************************************************
   subroutine lowerTri_Cmplx(A,L)
   !   extract lower triangular part of matrix A in CSR storage
      type(spMatCSR_Cmplx),intent(in)        ::  A
      type(spMatCSR_Cmplx),intent(inout)     ::  L
      integer           :: kk,n,m,i,j,nz
      integer, allocatable, dimension(:)  :: rowT

      m = A%nRow
      n = A%nCol
      allocate(rowT(m))
 
      !   first pass: find numbers of columns in each row of output
      rowT = 0
      do i = 1,m
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).le.i) then
               rowT(i) = rowT(i)+1
            endif
         enddo
      enddo
      nz = 0
      do i = 1,m
         nz = nz+rowT(i)
      enddo
      call create_spMatCSR(m,n,nz,L)
      L%lower = .true.

      !   set row array in output CSR matrix
      L%row(1) = 1
      do i = 1,L%nRow
         L%row(i+1) = L%row(i)+rowT(i)
      enddo
      deallocate(rowT)

      !    now fill in columns and values
      do i = 1,m
         kk = 0
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).le.i) then
               L%col(L%row(i)+kk) = A%col(j) 
               L%val(L%row(i)+kk) = A%val(j) 
               kk = kk+1
            endif
         enddo
      enddo
      return
   end subroutine
!*******************************************************************
   subroutine upperTri_Cmplx(A,U)
   !   extract lower triangular part of matrix A in CSR storage
      type(spMatCSR_Cmplx),intent(in)        ::  A
      type(spMatCSR_Cmplx),intent(inout)     ::  U
      integer           :: kk,n,m,i,j,nz
      integer, allocatable, dimension(:)  :: rowT

      m = A%nRow
      n = A%nCol
      allocate(rowT(m))
 
      !   first pass: find numbers of columns in each row of output
      rowT = 0
      do i = 1,m
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).ge.i) then
               rowT(i) = rowT(i)+1
            endif
         enddo
      enddo
      nz = 0
      do i = 1,m
         nz = nz+rowT(i)
      enddo
      call create_spMatCSR(m,n,nz,U)
      U%upper = .true.

      !   set row array in output CSR matrix
      U%row(1) = 1
      do i = 1,U%nRow
         U%row(i+1) = U%row(i)+rowT(i)
      enddo
      deallocate(rowT)

      !    now fill in columns and values
      do i = 1,m
         kk = 0
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).ge.i) then
               U%col(U%row(i)+kk) = A%col(j) 
               U%val(U%row(i)+kk) = A%val(j) 
               kk = kk+1
            endif
         enddo
      enddo
      return
   end subroutine
!*******************************************************************
   subroutine diag_Real(A,D)
   !   extract diagonal part of matrix A in CSR storage
      type(spMatCSR_Real),intent(in)        ::  A
      real(kind=prec),allocatable,intent(inout)  :: D(:)
      integer    :: n,m,i,j

      m = A%nRow
      n = A%nCol
      if(n.ne.m) then
         call errStop('diag only works for square matrices')
      endif
      if(allocated(D)) then
         deallocate(D)
      endif
      allocate(D(m))
      D = 0
      do i = 1,m
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).eq.i) then
               D(i) = A%val(j)
            endif
         enddo
      enddo
   end subroutine  
!*******************************************************************
   subroutine diag_Cmplx(A,D)
   !   extract diagonal part of matrix A in CSR storage
      type(spMatCSR_Cmplx),intent(in)        ::  A
      complex(kind=prec),allocatable,intent(inout)  :: D(:)
      integer    :: n,m,i,j

      m = A%nRow
      n = A%nCol
      if(n.ne.m) then
         call errStop('diag only works for square matrices')
      endif
      if(allocated(D)) then
         deallocate(D)
      endif
      allocate(D(m))
      D = 0
      do i = 1,m
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).eq.i) then
               D(i) = A%val(j)
            endif
         enddo
      enddo
   end subroutine
!*******************************************************************
   subroutine subMatrix_Real(A,r,c,B)
   !   extract submatrix of A with rows and columns given by integer
   !     arrays r and c 
      type(spMatCSR_Real),intent(in)        ::  A
      type(spMatCSR_Real),intent(inout)     ::  B
      integer, intent(in)                   ::  r(:),c(:)
      integer           :: kk,n,m,i,j,nz,k
      integer, allocatable, dimension(:)  :: rowT,colT

      m = size(r)
      n = size(c)
      allocate(rowT(m+1))
      allocate(colT(A%nCol))
      colT = 0
      do i = 1,n
         colT(c(i)) = i
      enddo

     !   count number of entries in each row
      rowT = 0
      nz = 0
      do i =1,m
         do j = A%row(r(i)),A%row(r(i)+1)-1
            if(colT(A%col(j)).gt.0) then
               rowT(i) = rowT(i)+1
               nz = nz+1
            endif
         enddo
      enddo

      call create_spMatCSR(m,n,nz,B)

      !   set row array in output CSR matrix
      B%row(1) = 1
      do i = 1,B%nRow
         B%row(i+1) = B%row(i)+rowT(i)
      enddo

      do i =1,m
         kk = 0
         do j = A%row(r(i)),A%row(r(i)+1)-1
            if(colT(A%col(j)).gt.0) then
               B%col(B%row(i)+kk) = colT(A%col(j))
               B%val(B%row(i)+kk) = A%val(j)
               kk = kk+1
            endif
         enddo
      enddo
      return
   end subroutine
!*******************************************************************
   subroutine subMatrix_Cmplx(A,r,c,B)
   !   extract submatrix of A with rows and columns given by integer
   !     arrays r and c 
      type(spMatCSR_Cmplx),intent(in)        ::  A
      type(spMatCSR_Cmplx),intent(inout)     ::  B
      integer, intent(in)	::      r(:),c(:)
      integer           :: kk,n,m,i,j,nz,k
      integer, allocatable, dimension(:)  :: rowT,colT

      m = size(r)
      n = size(c)
      allocate(rowT(m+1))
      allocate(colT(A%nCol))
      colT = 0
      do i = 1,n
         colT(c(i)) = i
      enddo

     !   count number of entries in each row
      rowT = 0
      nz = 0
      do i =1,m
         do j = A%row(r(i)),A%row(r(i)+1)-1
            if(colT(A%col(j)).gt.0) then
               rowT(i) = rowT(i)+1
               nz = nz+1
            endif
         enddo
      enddo

      call create_spMatCSR(m,n,nz,B)

      !   set row array in output CSR matrix
      B%row(1) = 1
      do i = 1,B%nRow
         B%row(i+1) = B%row(i)+rowT(i)
      enddo

      do i =1,m
         kk = 0
         do j = A%row(r(i)),A%row(r(i)+1)-1
            if(colT(A%col(j)).gt.0) then
               B%col(B%row(i)+kk) = colT(A%col(j))
               B%val(B%row(i)+kk) = A%val(j)
               kk = kk+1
            endif
         enddo
      enddo
      return
   end subroutine
!*******************************************************************
   subroutine BlkDiag_Real(A,B)
   !  merge an array of sparse matrices in CSR storage
   !       into a single block diagonal matrix B
      type(spMatCSR_Real), pointer, intent(in)        ::  A(:)
      type(spMatCSR_Real),intent(inout)     ::  B
 
      integer     :: nBlks,nRowB,nColB,nzB,i,i1,j1,i2,j2,k1

      nBlks = size(A) 
      nRowB = 0
      nColB = 0
      nzB = 0
      do i = 1,nBlks
         nRowB = nRowB + A(i)%nRow
         nColB = nColB + A(i)%nCol
         nzB = nzB + A(i)%row(A(i)%nRow+1)-1
      enddo
      call create_spMatCSR(nRowB,nColB,nzB,B)
      i1 = 1
      j1 = 1
      k1 = 0
      B%upper = .true.
      B%lower = .true.
      do i = 1,nBlks
         B%upper = B%upper .and. A(i)%upper
         B%lower = B%lower .and. A(i)%lower
         i2 = i1 + A(i)%nRow-1
         j2 = j1 + A(i)%row(A(i)%nRow+1)-2
         B%row(i1:i2) = A(i)%row(1:A(i)%nrow)+j1-1  
         B%col(j1:j2) = A(i)%col + k1
         B%val(j1:j2) = A(i)%val
         i1 = i2 + 1
         j1 = j2 + 1
         k1 = k1 + A(i)%nCol
      enddo
      B%row(i1) = j1
   end subroutine
!*******************************************************************
   subroutine BlkDiag_Cmplx(A,B)
   !  merge an array of sparse matrices in CSR storage
   !       into a single block diagonal matrix B
      type(spMatCSR_Cmplx), pointer, intent(in)        ::  A(:)
      type(spMatCSR_Cmplx),intent(inout)     ::  B
 
      integer     :: nBlks,nRowB,nColB,nzB,i,j,i1,j1,i2,j2,k1

      nBlks = size(A) 
      nRowB = 0
      nColB = 0
      nzB = 0
      do i = 1,nBlks
         nRowB = nRowB + A(i)%nRow
         nColB = nColB + A(i)%nCol
         nzB = nzB + A(i)%row(A(i)%nRow+1)-1
      enddo
      call create_spMatCSR(nRowB,nColB,nzB,B)
      i1 = 1
      j1 = 1
      k1 = 0
      B%upper = .true.
      B%lower = .true.
      do i = 1,nBlks
         B%upper = B%upper .and. A(i)%upper
         B%lower = B%lower .and. A(i)%lower
         i2 = i1 + A(i)%nRow-1 
         j2 = j1 + A(i)%row(A(i)%nRow+1)-2
         B%row(i1:i2) = A(i)%row(1:A(i)%nrow)+j1-1  
         B%col(j1:j2) = A(i)%col + k1
         !AK: COMPILER BUG: ifort (IFORT) 19.0.4.233 20190416 produced a segmentation fault with SP2 only (not with SP!)
         !    at the first encounter of the line "B%val(j1:j2) = A(i)%val" in BlkDiag_Cmplx in module spOpTools.f90.
         !    Something to do with vectorization, since this doesn't occur in debug mode. Rewritten to use a loop
         !    with error checking, for now. This is most likely a temporary fix that could be eliminated later.
         !B%val(j1:j2) = A(i)%val
         do j = 1,size(A(i)%val)
            if (isnan(real(A(i)%val(j)))) then
                write(0,*) 'ERROR: NaN in BlkDiag_Cmplx for i=',i,' and j=',j
            endif
            B%val(j1+j-1) = A(i)%val(j)
         enddo
         i1 = i2 + 1
         j1 = j2 + 1
         k1 = k1 + A(i)%nCol
      enddo
      B%row(i1) = j1
   end subroutine
!********************************************************************
   subroutine R2C_CSR(Ar,Ac)
   !   convert real CSR matrix to complex form
      type(spMatCSR_Real),intent(in)        ::  Ar
      type(spMatCSR_Cmplx),intent(inout)    ::  Ac
      integer     :: nRow, nCol, nz, i
    
      nRow = Ar%nRow
      nCol = Ar%nCol
      nz = Ar%row(nRow+1)-1
      call create_spMatCSR(nRow,nCol,nz,Ac)
      Ac%row = Ar%row
      Ac%col = Ar%col
      do i = 1,nz
         Ac%val(i) = cmplx(Ar%val(i),0)
      enddo
   end subroutine
!********************************************************************
   subroutine splitRMAT(A,i,np,B,isizes)
   ! automaticly split a matrix A into row submatrices 
   ! and take only the corresponding row submatrice of B
   ! for the (i+1)th process in n parallel threads
   ! this is used to prepare PETSc AIJ type matrix
   !
   ! Note: this can be easily changed to be split according to the
   ! number of non-zero elements in each submatrix 
   ! important!
   ! the submatrix B is modified to use zero-based index (as Petsc)
      type(spMatCSR_real),intent(in)        ::  A  ! original matrix
      type(spMatCSR_real),intent(inout)     ::  B  ! submatrix
      integer, intent(in)                   ::  i,np
      integer, intent(in), pointer, dimension(:), optional ::  isizes
      integer                               ::  istart,iend,nrow_l,j,k
      integer                               ::  m,n,nz,nz_l,nsub
      real                                  ::  nrow
      integer, allocatable, dimension(:)    ::  rowT,colT

      allocate(colT(A%nCol))
      colT = (/ (j,j=1,A%nCol) /)
      if (A%nrow .lt. np) then
          write(6,*) 'number of process is larger than number of rows!'
          stop
      elseif (np.eq.1) then
          !write(6,*) 'only one process, returning the original Matrix'
          m = A%nRow
          n = A%nCol
          nz = A%row(m+1)-1
          call create_spMatCSR_Real(m,n,nz,B)
          B%nRow=m
          B%nCol=n
          B%row=A%row-1
          B%col=A%col-1
          B%val=A%val
          return
      end if
      nrow = A%nRow
      !nz_l = floor((A%row(nrow+1)-1)/np)
      if (present(isizes)) then ! split into given sizes
          istart = 1
          do k = 1,i 
              istart = istart+ isizes(k)
          enddo
          iend = istart+isizes(k)-1
      else
          nrow_l = floor(nrow/np)
          istart=i*nrow_l+1
          if (i+1 .eq. np)  then! the last sub-matrix
              iend=nrow
          else
              iend=i*nrow_l+nrow_l
          end if
      endif
      allocate(rowT(iend-istart+1))
      rowT = (/ (j,j=istart,iend) /)
      call subMatrix_Real(A,rowT,colT,B)
      B%row=B%row-1
      B%col=B%col-1
      deallocate(colT)
      deallocate(rowT)
   end subroutine splitRMAT   
!********************************************************************
   subroutine splitCMAT(A,i,np,B,isizes)
   ! automaticly split a matrix A into row submatrices 
   ! and take only the corresponding row submatrice of B
   ! for the (i+1)th process in n parallel threads
   ! this is used to prepare PETSc AIJ type matrix
   !
   ! Note: this can be easily changed to be split according to the
   ! number of none zeros elemtents in each submatrix 
   ! the submatrix B is modified to use zero-based index (as Petsc)
      type(spMatCSR_cmplx),intent(in)     ::  A  ! original matrix
      type(spMatCSR_cmplx),intent(inout)  ::  B  ! submatrix
      integer, intent(in)                   ::  i,np
      integer, intent(in), pointer, dimension(:), optional ::  isizes
      integer                               ::  istart,iend,nrow_l,j,k
      integer                               ::  m,n,nz
      real                                  ::  nrow
      integer, allocatable, dimension(:)    ::  rowT,colT

      allocate(colT(A%nCol))
      colT = (/ (j,j=1,A%nCol) /)
      if (A%nrow .lt. np) then
          write(6,*) 'number of processes is larger than number of rows!'
          stop
      elseif (np.eq.1) then
          !write(6,*) 'only one process, returning the original Matrix'
          m = A%nRow
          n = A%nCol
          nz = A%row(m+1)-1
          call create_spMatCSR_Cmplx(m,n,nz,B)
          B%nRow=m
          B%nCol=n
          B%row=A%row-1
          B%col=A%col-1
          B%val=A%val
          return
      end if
      nrow = A%nRow
      !nz_l = floor((A%row(nrow+1)-1)/np)
      if (present(isizes)) then ! split according to the given sizes
          istart = 1
          do k = 1,i 
              istart = istart+ isizes(k)
          enddo
          iend = istart+isizes(k)-1
      else ! split evenly 
          nrow_l = floor(nrow/np)
          istart=i*nrow_l+1
          if (i+1 .eq. np)  then! the last sub-matrix
              iend=nrow
          else
              iend=i*nrow_l+nrow_l
          end if
      endif
      allocate(rowT(iend-istart+1))
      rowT = (/ (j,j=istart,iend) /)
      call subMatrix_Cmplx(A,rowT,colT,B)
      B%row=B%row-1
      B%col=B%col-1
      deallocate(colT)
      deallocate(rowT)
   end subroutine splitCMAT   
!********************************************************************
   subroutine LTsolve_Cmplx(L,b,x)
      ! solve system Lx = b for complex vector x, lower triangular L
      !   here real or cmplx refers to U; x is always complex
      type(spMatCSR_Cmplx),  intent(in)         :: L
      complex(kind=prec), dimension(:), intent(in)   :: b
      complex(kind=prec), dimension(:), intent(inout)   :: x

      complex(kind=prec)               :: d
      integer     :: i,j

      if(L%nRow .ne.L%nCol) then
         call errStop('LTsolve: sparse matrix must be square')
      endif 
      if(.not.L%lower) then
         call errStop('LTsolve: sparse matrix must be lower triangular')
      endif 
      if(size(x).ne.L%nRow) then
         call errStop('LTsolve: output vector x not of correct size')
      endif 
      do i = 1,L%nRow
         x(i) = b(i)
         do j = L%row(i),L%row(i+1)-1
            if(L%col(j).lt.i) then
               x(i) = x(i)-L%val(j)*x(L%col(j))
            else
               !   in this case L%col(j) = i
               d = L%val(j)
            endif
         enddo
         x(i) = x(i)/d
      enddo
      return 
   end subroutine
!********************************************************************
   subroutine UTsolve_Cmplx(U,b,x)
      ! solve system Ux = b for complex vector x, upper triangular U
      !   here real or cmplx refers to U; x is always complex
      type(spMatCSR_Cmplx),  intent(in)         :: U
      complex(kind=prec), dimension(:), intent(in)   :: b
      complex(kind=prec), dimension(:), intent(inout)   :: x
     
      integer                          :: i,j
      complex(kind=prec)               :: d

      if(U%nRow .ne.U%nCol) then
         call errStop('UTsolve: sparse matrix must be square')
      endif 
      if(.not.U%upper) then
         call errStop('UTsolve: sparse matrix must be upper triangular')
      endif 
      if(size(x).ne.U%nRow) then
         call errStop('UTsolve: output vector x not of correct size')
      endif 
      do i = U%nRow,1,-1
         x(i) = b(i)
         do j = U%row(i),U%row(i+1)-1
            if(U%col(j).gt.i) then
               x(i) = x(i)-U%val(j)*x(U%col(j))
            else
               !   in this case U%col(j) = i
               d = U%val(j)
            endif
         enddo
         x(i) = x(i)/d
      enddo
      return 
   end subroutine
!********************************************************************
   subroutine LTsolve_Real(L,b,x)
      ! solve system Lx = b for complex vector x, lower triangular L
      !   here real or cmplx refers to L; x is always complex
      type(spMatCSR_Real),  intent(in)         :: L
      complex(kind=prec), dimension(:), intent(in)   :: b
      complex(kind=prec), dimension(:), intent(inout)   :: x

      real(kind=prec)                                   :: d
      integer                                           :: i,j

      if(L%nRow .ne.L%nCol) then
         call errStop('LTsolve: sparse matrix must be square')
      endif 
      if(.not.L%lower) then
         call errStop('LTsolve: sparse matrix must be lower triangular')
      endif 
      if(size(x).ne.L%nRow) then
         call errStop('LTsolve: output vector x not of correct size')
      endif 
      do i = 1,L%nRow
         x(i) = b(i)
         do j = L%row(i),L%row(i+1)-1
            if(L%col(j).lt.i) then
               x(i) = x(i)-L%val(j)*x(L%col(j))
            else
               !   in this case L%col(j) = i
               d = L%val(j)
            endif
         enddo
         x(i) = x(i)/d
      enddo
      return 
   end subroutine
!********************************************************************
   subroutine UTsolve_Real(U,b,x)
      ! solve system Ux = b for complex vector x, upper triangular U
      !   here real or cmplx refers to L; x is always complex
      type(spMatCSR_Real),  intent(in)         :: U
      complex(kind=prec), dimension(:), intent(in)   :: b
      complex(kind=prec), dimension(:), intent(inout)   :: x
     
      integer                                           :: i,j
      real(kind=prec)                                   :: d

      if(U%nRow .ne.U%nCol) then
         call errStop('UTsolve: sparse matrix must be square')
      endif 
      if(.not.U%upper) then
         call errStop('UTsolve: sparse matrix must be upper triangular')
      endif 
      if(size(x).ne.U%nRow) then
         call errStop('UTsolve: output vector x not of correct size')
      endif 
      do i = U%nRow,1,-1
         x(i) = b(i)
         do j = U%row(i),U%row(i+1)-1
            if(U%col(j).gt.i) then
               x(i) = x(i)-U%val(j)*x(U%col(j))
            else
               !   in this case U%col(j) = i
               d = U%val(j)
            endif
         enddo
         x(i) = x(i)/d
      enddo
      return 
   end subroutine
!********************************************************************
   subroutine CSR_R2Cdiag(A,d,B)
   !    specialized routine to add imaginary D to diagonal of real sparse
   !        matrix A, saving output to complex sparse matrix B
      type(spMatCSR_Real),  intent(in)             :: A
      real(kind=prec), dimension(:), intent(in)    :: d
      type(spMatCSR_Cmplx),  intent(inout)         :: B
      integer                :: n,m,nz,i,j

      m = A%nRow
      n = A%nCol
      nz = A%row(m+1)-1
      call create_spMatCSR_Cmplx(m,n,nz,B)
      B%row  = A%row
      B%col = A%col
      B%val = A%val
      do i =1,m
         do j =B%row(i),B%row(i+1)-1
            !   not sure ISIGN should be here!
            if(B%col(j) .eq. i) then
               B%val(j) = B%val(j)+ISIGN*CMPLX(0.0,1.0,8)*d(i)
               exit
            endif
         enddo
      enddo
   
   end subroutine
!************************************************************************
   subroutine CholInc_Real(A,L)
   !   incomplete cholesky decomposition of a symmetric matrix A
   !   This assumes the matrix is symmetric, but does not check
   !   ALSO: this should only be used with positive definite matrices
   !         Could act up if matrix is not symmetric pos-def

      type(spMatCSR_Real),intent(in)     ::  A
      type(spMatCSR_Real),intent(inout)  ::  L
      integer    :: ii,mMax,n,m,i,j,nij,nji,k,k1,nij1,nji1,i1,j1
      integer, allocatable, dimension(:) :: ij,ji,ji1

      m = A%nRow
      n = A%nCol
      mMax = maxColumnsR(A)
      allocate(ij(mMax))
      allocate(ji(mMax))
      allocate(ji1(mMax))
      call lowerTri_Real(A,L)
      do i = 1,m ! sweep in rows
         !  divide up columns in row i: diagonal, to left, to right
         nij = 0
         nji = 0
         do j = A%row(i),A%row(i+1)-1 !sweep in columns
            if(i.eq.A%col(j)) then
               !   this index should correspond to positions where
               !   columns and values are stored in L
               ii = j - A%row(i) + L%row(i)
            elseif(i.lt.A%col(j)) then
               !    these are indicies into A matrix storage
               nij = nij+1
               ij(nij) = j
            else
               nji = nji+1
               !   these indicies should correspond to positions where
               !   columns and values are stored in L
               ji(nji) = j - A%row(i) + L%row(i)
            endif
         enddo
         !   now compute elements of column i for L
         !    diagonal ii
         do k = 1,nji
            L%val(ii) = L%val(ii) - L%val(ji(k))*L%val(ji(k))
         enddo
         if(L%val(ii).lt.0) then
             L%val(ii) = -L%val(ii)
         endif
         L%val(ii) = sqrt(L%val(ii)) ! diagonal
         do j = 1,nij
            !   these are rows of L that have elements in column i
            i1 = A%col(ij(j))
            nji1 = 0
            !  divide up elements in this row: diagonal, to left
            do k = L%row(i1),L%row(i1+1)-1
               if(i.eq.L%col(k)) then
                  j1 = k
               else
                  nji1 = nji1+1 
                  ji1(nji1) = k
               endif
            enddo
            do k = 1,nji
               do k1  = 1,nji1
                  if(L%col(ji(k)).eq.L%col(ji1(k1))) then 
                     L%val(j1) = L%val(j1) -     &
                              L%val(ji(k))*L%val(ji1(k1))
                  endif
               enddo
            enddo
            L%val(j1) = L%val(j1)/L%val(ii) 
         enddo
      enddo
      deallocate(ji)
      deallocate(ji1)
      deallocate(ij)
   end  subroutine            
!*****************************************************************************
   subroutine Dilu_Real(A,L,U)
   !    this mimics approach used in ModEM --- D-ILU
   !    NOT ILU-0
   !    THIS ASSUMES THE MATRIX IS SYMMETRIC -- but not necessarily Hermitian
      type(spMatCSR_Real),intent(in)           :: A
      type(spMatCSR_Real),intent(inout)        :: L,U

      real(kind=prec),allocatable, dimension(:)    ::d
      integer             ::n,m,nz,i,j

      call lowerTri(A,L)
      call upperTri(A,U)
      m = A%nRow
      allocate(d(m))
      do i = 1,m
         d(i) = 0.0_dp
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).eq.i) then
               d(i) = d(i) + A%val(j)
            elseif(A%col(j).lt.i) then
               d(i) = d(i) - A%val(j)*A%val(j)*d(A%col(j))
            endif
         enddo
         d(i) = 1.0_dp/d(i)
      enddo
      do i = 1,m
         do j = L%row(i),L%row(i+1)-1
            if(L%col(j).eq.i) then
               L%val(j) = 1
            else
               L%val(j) = L%val(j)*d(L%col(j))
            endif
         enddo
      enddo
      do i = 1,m
         do j = U%row(i),U%row(i+1)-1
            if(U%col(j).eq.i) then
               U%val(j) = 1.0_dp/d(i)
               exit
            endif
         enddo
      enddo
      return
   end subroutine
!*****************************************************************************
   subroutine Dilu_Cmplx(A,L,U)
   !    this mimics approach used in ModEM --- D-ILU(diagonal-ILU)
   !    NOT ILU-0
   !    THIS ASSUMES THE MATRIX IS SYMMETRIC -- but not necessarily Hermitian
   !    so it will not work if using modified system equation (as Randy did)
      type(spMatCSR_Cmplx),intent(in)           :: A
      type(spMatCSR_Cmplx),intent(inout)        :: L,U

      complex(kind=prec),allocatable, dimension(:)    ::d
      integer             ::n,m,nz,i,j

      call lowerTri(A,L)
      call upperTri(A,U)
      m = A%nRow
      allocate(d(m))
      do i = 1,m
         d(i) = CMPLX(0.0,0.0,8)
         do j = A%row(i),A%row(i+1)-1
            if(A%col(j).eq.i) then
               d(i) = d(i) + A%val(j)
            elseif(A%col(j).lt.i) then
               d(i) = d(i) - A%val(j)*A%val(j)*d(A%col(j))
            endif
         enddo
         d(i) = CMPLX(1.0,0.0,8)/d(i)
      enddo
      do i = 1,m
         do j = L%row(i),L%row(i+1)-1
            if(L%col(j).eq.i) then
               L%val(j) = 1
            else
               L%val(j) = L%val(j)*d(L%col(j))
            endif
         enddo
      enddo
      do i = 1,m
         do j = U%row(i),U%row(i+1)-1
            if(U%col(j).eq.i) then
               U%val(j) = 1.0_dp/d(i)
               exit
            endif
         enddo
      enddo
      return
   end subroutine
!*****************************************************************************
   subroutine Dilu_Cmplx_AS(A,L,U)
   !    this mimics approach used in ModEM --- D-ILU(diagonal-ILU)
   !    NOT ILU-0
   !    slightly modified to be used with modified system equations 
      type(spMatCSR_Cmplx),intent(in)           :: A
      type(spMatCSR_Cmplx),intent(inout)        :: L,U
      type(spMatCSR_Cmplx)                      :: UT

      complex(kind=prec),allocatable, dimension(:)    ::d
      integer             ::n,m,nz,i,j

      call lowerTri(A,L)
      call upperTri(A,U)
      call CMATtrans(U,UT)
      m = L%nRow
      allocate(d(m))
      do i = 1,m ! loop through rows
         d(i) = C_ZERO
         do j = L%row(i),L%row(i+1)-1 !loop through columns
            if(L%col(j).eq.i) then ! diagonal
               d(i) = d(i) + L%val(j) 
            elseif(L%col(j).lt.i) then ! take the L and U side
               d(i) = d(i) - L%val(j)*UT%val(j)*d(L%col(j))
            endif
         enddo
         d(i) = C_ONE/d(i)
      enddo
      call deall_spMatCSR(UT)
      do i = 1,m
         do j = L%row(i),L%row(i+1)-1
            if(L%col(j).eq.i) then
               L%val(j) = 1
            else
               L%val(j) = L%val(j)*d(L%col(j))
            endif
         enddo
      enddo
      do i = 1,m
         do j = U%row(i),U%row(i+1)-1
            if(U%col(j).eq.i) then
               U%val(j) = 1.0_dp/d(i)
               exit
            endif
         enddo
      enddo
      return
   end subroutine
!*****************************************************************************
   subroutine ilu0_Cmplx(A,L,U)
      type(spMatCSR_Cmplx),intent(in)           :: A
      type(spMatCSR_Cmplx),intent(inout)        :: L,U
      type(spMatCSR_Cmplx)                      :: Atmp
      complex(kind=prec),allocatable, dimension(:)    ::d
      complex(kind=prec)                        :: piv
      integer                                   :: n,m,nz,i,j,j2,k,p
      ! a simple but not at all intuitive ILU0 routine with CSR sparse matrix 
      !
      ! the algorithm is modified from a much fancier version in Yousef Saad's 
      ! 2003 book <Iterative Methods for Sparse Linear Systems>, Chapter 10.2
      !
      ! apparently the original algorithm seeks to minimize the memory 
      ! consumption by storing the L and U in the original sparse matrix 
      ! structure of A (as ILU0 does not have any fill-ins).
      allocate(d(A%nRow))
      n = A%nRow
      m = A%nCol
      nz = A%row(A%nRow+1)-1
      call create_spMatCSR(n,m,nz,Atmp)
      Atmp%row=A%row
      Atmp%col=A%col
      Atmp%val=A%val
      call sort_spMatCSR(Atmp) ! sort the col indices in A
      d = C_ZERO
      do i=1,Atmp%nRow ! loop through rows

          p=Atmp%row(i) ! mark the first none zero element in current row

          do j=Atmp%row(i),Atmp%row(i+1)-1 ! loop through columns

              if(Atmp%col(j).eq.i) then !diagonal
                  d(i) = Atmp%val(j) ! store previous diagonal elements
                  exit ! exit as we reached the last element in L
              elseif(Atmp%col(j).lt.i) then
                  if (d(Atmp%col(j)).eq.C_ZERO) then
                      write(6,*) 'error: zero pivoting in ILU0 '
                      write(6,*) 'in col: ',                              &
     &                            Atmp%col(Atmp%row(i):Atmp%row(i+1)-1)
                      write(6,*) 'in row: ', Atmp%col(j)
                      stop
                  endif 
                  ! first divide each row in L with diagonal elements
                  Atmp%val(j) = Atmp%val(j)/d(Atmp%col(j))
                  piv = Atmp%val(j)
                  do k = p+1,Atmp%row(i+1)-1 ! then adding up each column in U
                      do j2 = Atmp%row(Atmp%col(j)),Atmp%row(Atmp%col(j)+1)-1
                          if(Atmp%col(j2).eq.Atmp%col(k)) then
                              Atmp%val(k) = Atmp%val(k) - piv* Atmp%val(j2)
                              exit
                          endif
                      enddo
                  enddo
                  p = p + 1 ! 
              endif

          enddo

      enddo
      deallocate(d)
      call lowerTri(Atmp,L)
      do i = 2,L%nRow+1 ! set the diagonal element of L to 1
          ! the diagonal element should be the last in each row 
          ! ONLY IF col is properly sorted
          L%val(L%row(i)-1) = C_ONE
      end do
      call upperTri(Atmp,U)
      call deall_spMatCSR(Atmp)
      return
   end subroutine ilu0_Cmplx
!*****************************************************************************
   subroutine ilu0_Real(A,L,U)
      type(spMatCSR_Real),intent(in)           :: A
      type(spMatCSR_Real),intent(inout)        :: L,U
      type(spMatCSR_Real)                      :: Atmp
      real(kind=prec),allocatable, dimension(:):: d
      real(kind=prec)                          :: piv
      integer                                  :: n,m,nz,i,j,j2,k,p
      ! a simple but not at all intuitive ILU0 routine with CSR sparse matrix 
      !
      ! the algorithm is modified from a much fancier version in Yousef Saad's 
      ! 2003 book <Iterative Methods for Sparse Linear Systems>, Chapter 10.2
      !
      ! apparently the original algorithm seeks to minimize the memory 
      ! consumption by storing the L and U in the original sparse matrix 
      ! structure of A (as ILU0 does not have any fill-ins).
      allocate(d(A%nRow))
      n = A%nRow
      m = A%nCol
      nz = A%row(A%nRow+1)-1
      call create_spMatCSR(n,m,nz,Atmp)
      Atmp%row=A%row
      Atmp%col=A%col
      Atmp%val=A%val
      call sort_spMatCSR(Atmp) ! sort the col indices in A
      d = 0.0
      do i=1,Atmp%nRow ! loop through rows
          p=Atmp%row(i) ! mark the first none zero element in current row
 
          do j=Atmp%row(i),Atmp%row(i+1)-1 ! loop through columns
              if(Atmp%col(j).eq.i) then !diagonal
                  d(i) = Atmp%val(j) ! store previous diagonal elements
                  exit ! exit as we reached the last element in L
              elseif(Atmp%col(j).lt.i) then
                  if (d(Atmp%col(j)).eq.0.0) then
                      write(6,*) 'error: zero pivoting in ILU0 '
                      stop
                  endif 
                  ! first divide each row in L with diagonal elements
                  Atmp%val(j) = Atmp%val(j)/d(Atmp%col(j))
                  piv = Atmp%val(j)
                  do k = p+1,Atmp%row(i+1)-1 
                      ! looping through every none-zero column in current row
                      ! write(6,*) Atmp%col(k),Atmp%val(k)
                      do j2 = Atmp%row(Atmp%col(j)),Atmp%row(Atmp%col(j)+1)-1
                      ! adding up each column of previous rows
                          if(Atmp%col(j2).eq.Atmp%col(k)) then ! column matches
                              Atmp%val(k) = Atmp%val(k) - piv* Atmp%val(j2)
                              exit
                          endif
                      enddo
                  enddo
                  p = p + 1 ! 
              endif

          enddo

      enddo
      deallocate(d)
      call lowerTri(Atmp,L)
      do i = 2,L%nRow+1 ! set the diagonal element of L to 1
          ! the diagonal element should be the last in each row 
          ! ONLY IF col is properly sorted
          L%val(L%row(i)-1) = 1.0
      enddo
      call upperTri(Atmp,U)
      call deall_spMatCSR(Atmp)
      return
   end subroutine ilu0_Real
!*******************************************************************
   subroutine RMATplusRMAT(A,B,C,tol)
     !  matix-matrix sum, real version
      type(spMatCSR_Real),  intent(in)         :: A
      type(spMatCSR_Real),  intent(in)         :: B
      type(spMatCSR_Real),  intent(inout)      :: C

      real(kind=prec),intent(in), optional     :: tol
      type(spMatCSR_Real)      :: Ctemp


      integer   :: i,j,k,m,n,nz,jj,i1,j1,j2,nnz,nzero
      logical               new
      real(kind=prec)    :: test,droptol,temp

      !   test for size consistency
      if((A%nCol.ne.B%nCol).and.(A%nRow.ne.B%nRow)) then
         call errStop('RMATpRMAT: matrix sizes incompatible')
      endif
      !    tolerance for dropping small entries derived as sums
      if(present(tol)) then
         droptol = tol
      else
         droptol = 0.0_dp
      endif

      !    create Ctemp  with space for all possible entries
      m = A%nRow
      n = A%nCol
      nz = A%row(m+1)+B%row(m+1)-2
      call create_spMatCSR(m,n,nz,Ctemp)

      !   copy A into Ctemp, and add B, row by row.   
      !   Some entries of the sum may become zero
      j1 = 1
      nzero = 0
      do i=1,m
         Ctemp%row(i) = A%row(i+1)-A%row(i)
         j2 = j1+Ctemp%row(i)-1
         Ctemp%col(j1:j2) = A%col(A%row(i):A%row(i+1)-1)
         Ctemp%val(j1:j2) = A%val(A%row(i):A%row(i+1)-1)
         !    add elements of B in row i
         do k = B%row(i),B%row(i+1)-1
            new = .true.
            do j = A%row(i),A%row(i+1)-1
               if(B%col(k).eq.A%col(j)) then
                 !   add maatrix elements
                  new =.false.
                  jj = j1+j-A%row(i)
                  temp = B%val(k)+A%val(j)
                  if(A%val(j).ne.0.0) then
                      test = abs(temp)/abs(A%val(j))
                  else
                      test = abs(temp)/abs(B%val(k))
                  endif
                  if(test.gt. droptol) then
                     Ctemp%val(jj) = temp
                  else
                     nzero = nzero + 1
                     Ctemp%val(jj) = 0
                  endif
                  exit
               endif
            enddo

            if(new) then
               j2 = j2+1
               Ctemp%row(i) = Ctemp%row(i)+1
               Ctemp%col(j2) = B%col(k)
               Ctemp%val(j2) = B%val(k)
            endif
         enddo 
         j1 = j2+1
      enddo
      !   Now Ctemp contains all sums, but possibly some zeros where
      !    A and B values have cancelled out.   Also Ctemp%row contains
      !      number of elements (including possible zeros), not limits
      !    of col and val arrays.   So clean up Ctemp, put results in C
      nz = sum(Ctemp%row(1:m)) - nzero
      ! write(0,*) 'nz = ',nz, ' nzero = ',nzero
      call create_spMatCSR(m,n,nz,C)
      j1 = 1
      i1 = 0
      nz = 1
      do i = 1,m
         nnz = 0
         j2 = j1+Ctemp%row(i)-1  
         do j = j1,j2
            if(abs(Ctemp%val(j)).gt. 0) then
               nnz = nnz + 1             
               i1 = i1+1
               C%col(i1) = Ctemp%col(j)
               C%val(i1) = Ctemp%val(j)
            endif
         enddo
         j1 = j2+1
         C%row(i) = nz
         nz = nz + nnz  
      enddo
      C%row(m+1) = nz
      call deall_spMatCSR(Ctemp)
   end subroutine
!*******************************************************************
   subroutine sort_spMatCSR_real(A)
   ! sort the CSR col indices (and correspnding vals) in each row into
   ! ascent order
   ! reorder col indices like col:  4  1  3  2 -> 1  2  3  4
   !                          val:  3  1  1 -1 -> 1 -1  1  2 
      implicit none
      type(spMatCSR_Real),  intent(inout)         :: A
      integer                                     :: i,j1,j2,k
      integer,allocatable,dimension(:)            :: idx
      real(kind=prec),allocatable,dimension(:)    :: col
      do i = 2, A%nRow+1
          j1=A%row(i-1)
          j2=A%row(i)-1
          allocate(col(j2-j1+1))
          allocate(idx(j2-j1+1))
          idx=(/(k, k=j1,j2, 1)/)
          col=real(A%col(j1:j2))
          call QSort(col,idx) 
          A%col(j1:j2)=int(col)
          A%val(j1:j2)=A%val(idx)
          deallocate(col)
          deallocate(idx)
      enddo
   end subroutine
!*******************************************************************
   subroutine sort_spMatCSR_Cmplx(A)
   ! sort the CSR col indices (and correspnding vals) in each row into
   ! ascent order
   ! reorder col indices like col:  4  1  3  2 -> 1  2  3  4
   !                          val:  3  1  1 -1 -> 1 -1  1  2 
      implicit none
      type(spMatCSR_Cmplx),  intent(inout)        :: A
      integer                                     :: i,j1,j2,k
      integer,allocatable,dimension(:)            :: idx
      real(kind=prec),allocatable,dimension(:)    :: col
      do i = 2, A%nRow+1
          j1=A%row(i-1)
          j2=A%row(i)-1
          allocate(col(j2-j1+1))
          allocate(idx(j2-j1+1))
          idx=(/(k, k=j1,j2, 1)/)
          col=real(A%col(j1:j2))
          call QSort(col,idx) 
          !write(6,*) 'current row is: ', i-1
          !write(6,*) A%col(j1:j2)
          A%col(j1:j2)=int(col)
          A%val(j1:j2)=A%val(idx)
          deallocate(col)
          deallocate(idx)
      enddo
   end subroutine
!*******************************************************************
   subroutine RMAT2CMAT(R,C)
   ! a silly subroutine to convert the real CSR sp matrix into complex
   ! WARNING: for now this is only used to form PETSc mat structure
   ! 
   ! this assumes the matrix indices starts from ZERO instead of ONE
   ! need to be modified to use in other circumstances
      implicit none
      type(spMatCSR_Real),   intent(in)           :: R
      type(spMatCSR_Cmplx),  intent(inout)        :: C
      integer                                     :: m,n,nnz
      m = R%nRow
      n = R%nCol
      nnz = R%row(R%nRow+1)
      call create_spMatCSR_Cmplx(m,n,nnz,C)
      C%row = R%row
      C%col = R%col
      C%val = R%val
      return
   end subroutine
!*******************************************************************
end module
