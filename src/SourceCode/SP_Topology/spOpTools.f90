!
!> Operations over sparse matrices in CSR storage
!
module SpOpTools
    !
    use Utilities
    use MetricElements
    !
    !> Generic matrix types and tools, using CSR storage, 
    !> but with fortran numbering conventions(starting from 1)
    type :: spMatCSR_Real
        integer :: nRow=0
        integer :: nCol=0
        integer, allocatable, dimension(:) :: row, col
        real( kind=prec ), allocatable, dimension(:) :: val
        logical :: is_allocated = .FALSE.
        logical :: lower = .FALSE.
        logical :: upper = .FALSE.
    end type
    !
    type :: spMatCSR_Cmplx
        integer :: nRow=0
        integer :: nCol=0
        integer, allocatable, dimension(:) :: row, col
        complex( kind=prec ), allocatable, dimension(:) :: val
        logical :: is_allocated = .FALSE.
        logical :: lower = .FALSE.
        logical :: upper = .FALSE.
    end type
    !
    type :: spMatIJS_Real
        integer :: nRow=0
        integer :: nCol=0
        integer, allocatable, dimension(:) :: I, J
        real( kind=prec ), allocatable, dimension(:) :: S
        logical :: is_allocated = .FALSE.
        logical :: lower = .FALSE.
        logical :: upper = .FALSE.
    end type
    !
    type :: spMatIJS_Cmplx
        integer :: nRow=0
        integer :: nCol=0
        integer, allocatable, dimension(:) :: I, J
        complex( kind=prec ), allocatable, dimension(:) :: S
        logical :: is_allocated = .FALSE.
        logical :: lower = .FALSE.
        logical :: upper = .FALSE.
    end type
    !
    interface create_spMatCSR
        module procedure create_spMatCSR_Real
        module procedure create_spMatCSR_Cmplx
    end interface
    !
    interface create_spMatIJS
        module procedure create_spMatIJS_Real
        module procedure create_spMatIJS_Cmplx
    end interface
    !
    interface deall_spMatCSR
        module procedure deall_spMatCSR_Real
        module procedure deall_spMatCSR_Cmplx
    end interface
    !
    interface deall_spMatIJS
        module procedure deall_spMatIJS_Real
        module procedure deall_spMatIJS_Cmplx
    end interface
    !
    interface CSR2IJS
        module procedure CSR2IJS_Real
        module procedure CSR2IJS_Cmplx
    end interface
    !
    interface IJS2CSR
        module procedure IJS2CSR_Real
        module procedure IJS2CSR_Cmplx
    end interface
    !
    interface lowerTri
        module procedure lowerTri_Real
        module procedure lowerTri_Cmplx
    end interface
    !
    interface upperTri
        module procedure upperTri_Real
        module procedure upperTri_Cmplx
    end interface
    !
    interface LTsolve
        module procedure LTsolve_Real
        module procedure LTsolve_Cmplx
    end interface
    !
    interface UTsolve
        module procedure UTsolve_Real
        module procedure UTsolve_Cmplx
    end interface
    !
    interface SubMatrix
        module procedure SubMatrix_Real
        module procedure SubMatrix_Cmplx
    end interface
    !
    interface splitMAT
        module procedure splitRMAT
        module procedure splitCMAT
    end interface
    !
    interface RMATxVEC
        module procedure RMATxCVEC
        module procedure RMATxRVEC
    end interface
    !
    interface sort_spMatCSR
        module procedure sort_spMatCSR_real
        module procedure sort_spMatCSR_cmplx
    end interface
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine writeIJS_Matrix( matrix_csr, i_unit )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: matrix_csr
        
        integer, intent( in ) :: i_unit
        !
        type( spMatIJS_Real ) :: matrix_ijs
        integer :: i
        !
        call create_spMatIJS_Real( size( matrix_csr%row ), size( matrix_csr%col ), size( matrix_csr%val ), matrix_ijs )
        !
        call CSR2IJS( matrix_csr, matrix_ijs )
        !
        do i = 1, matrix_ijs%nCol
            !
            !write( *, * ) matrix_ijs%nRow, matrix_ijs%nCol, size( matrix_csr%row ), size( matrix_csr%col ), size( matrix_csr%val )
            write( i_unit, * ) matrix_ijs%I(i), matrix_ijs%J(i), matrix_ijs%S(i)
            !
        enddo
        !
    end subroutine writeIJS_Matrix
    !
    !> No subroutine briefing
    !
    subroutine create_spMatCSR_Real( m, n, nz, A )
        implicit none
        !
        integer, intent( in ) :: m, n, nz
        ! A will be sparse m x n with nz non-zero elements
        type( spMatCSR_Real ), intent( inout ) :: A
        !
        integer :: istat
        !
        if( A%is_allocated ) then
            call deall_spMatCSR(A)
        endif
        !
        A%nRow = m
        A%nCol = n
        allocate( A%row(m+1), stat = istat )
        allocate( A%col(nz), stat = istat )
        allocate( A%val(nz), stat = istat )
        A%row(m+1) = nz + 1
        !
        A%is_allocated = .TRUE.
        !
    end subroutine create_spMatCSR_Real
    !
    ! No subroutine briefing
    !
    subroutine create_spMatCSR_Cmplx( m, n, nz, A )
        implicit none
        !
        integer, intent( in ) :: m, n, nz
        !   A will be sparse m x n with nz non-zero elements
        type( spMatCSR_Cmplx ), intent( inout ) :: A
        !
        if( A%is_allocated ) then
            call deall_spMatCSR( A )
        endif
        !
        A%nRow = m
        A%nCol = n
        !
        allocate( A%row(m+1) )
        allocate( A%col(nz) )
        allocate( A%val(nz) )
        !
        A%row(m+1)=nz+1
        A%is_allocated = .TRUE.
        !
    end subroutine create_spMatCSR_Cmplx
    !
    !> No subroutine briefing
    !
    subroutine create_spMatIJS_Real( m, n, nz, A )
        implicit none
        !
        integer, intent( in ) :: m, n, nz
        ! A will be sparse m x n with nz non-zero elements
        type(spMatIJS_Real), intent( inout ) :: A
        !
        if( A%is_allocated ) then
            call deall_spMatIJS( A )
        endif
        !
        A%nRow = m
        A%nCol = n
        !
        allocate( A%I(nz) )
        allocate( A%J(nz) )
        allocate( A%S(nz) )
        !
        A%is_allocated = .TRUE.
        !
    end subroutine create_spMatIJS_Real
    !
    !> No subroutine briefing
    !
    subroutine create_spMatIJS_Cmplx( m, n, nz, A )
        implicit none
        !
        integer, intent( in ) :: m, n, nz
        !   A will be sparse m x n with nz non-zero elements
        type(spMatIJS_Cmplx), intent( inout ) :: A
        !
        if( A%is_allocated ) then
            call deall_spMatIJS( A )
        endif
        !
        A%nRow = m
        A%nCol = n
        !
        allocate( A%I(nz) )
        allocate( A%J(nz) )
        allocate( A%S(nz) )
        !
        A%is_allocated = .TRUE.
        !
    end subroutine create_spMatIJS_Cmplx
    !
    !> No subroutine briefing
    !
    subroutine deall_spMatCSR_Real( A )
        implicit none
        !
        type( spMatCSR_Real ) :: A
        !
        if( A%is_allocated ) then
            !
            deallocate(A%row)
            deallocate(A%col)
            deallocate(A%val)
            A%is_allocated = .FALSE.
            A%upper = .FALSE.
            A%lower = .FALSE.
            A%nRow = 0
            A%nCol = 0
            !
        endif
        !
    end subroutine deall_spMatCSR_Real
    !
    !> No subroutine briefing
    !
    subroutine deall_spMatCSR_Cmplx( A )
        implicit none
        !
        type( spMatCSR_Cmplx ) :: A
        !
        if( A%is_allocated ) then
            !
            deallocate(A%row)
            deallocate(A%col)
            deallocate(A%val)
            A%is_allocated = .FALSE.
            A%upper = .FALSE.
            A%lower = .FALSE.
            A%nRow = 0
            A%nCol = 0
            !
        endif
    end subroutine deall_spMatCSR_Cmplx
    !
    !> No subroutine briefing
    !
    subroutine deall_spMatIJS_Real( A )
        implicit none
        !
        type( spMatIJS_Real ) :: A
        !
        if( A%is_allocated ) then
            !
            deallocate(A%I)
            deallocate(A%J)
            deallocate(A%S)
            A%is_allocated = .FALSE.
            A%nRow = 0
            A%nCol = 0
            !
        endif
        !
    end subroutine deall_spMatIJS_Real
    !
    !> No subroutine briefing
    !
    subroutine deall_spMatIJS_Cmplx( A )
        implicit none
        !
        type( spMatIJS_Cmplx ) :: A
        !
        if( A%is_allocated ) then
            !
            deallocate(A%I)
            deallocate(A%J)
            deallocate(A%S)
            A%is_allocated = .FALSE.
            A%nRow = 0
            A%nCol = 0
            !
        endif
        !
    end subroutine deall_spMatIJS_Cmplx
    !
    !> No function briefing
    !
    function sameSizeCSR_Real( A, B ) result( same_size )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A, B
        !
        logical :: same_size
        !
        same_size = .FALSE.
        !
        if( A%is_allocated .AND. B%is_allocated .AND. A%nRow .GT. 0 .AND. B%nRow .GT. 0 ) then 
            same_size = A%nRow .EQ. B%nRow .AND. A%nCol .EQ. B%nCol .AND. &
            A%row( A%nRow + 1 ) .EQ. B%row( B%nRow + 1 )
        endif
        !
    end function sameSizeCSR_Real
    !
    !> No function briefing
    !
    function sameSizeCSR_Cmplx( A, B ) result( same_size )
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A, B
        !
        logical :: same_size
        !
        same_size = .FALSE.
        !
        if( A%is_allocated .AND. B%is_allocated .AND. A%nRow .GT. 0 .AND. B%nRow .GT. 0 ) then 
            same_size = A%nRow .EQ. B%nRow .AND. A%nCol .EQ. B%nCol .AND. &
            A%row( A%nRow + 1 ) .EQ. B%row( B%nRow + 1 )
        endif
        !
    end function sameSizeCSR_Cmplx
    !
    !> No function briefing
    !
    integer function maxColumnsR(A) 
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        integer :: i
        maxColumnsR = 0
        do i = 1, A%nRow
        maxColumnsR = max(maxColumnsR, A%row(i+1)-A%row(i))
        enddo
        !
    end function maxColumnsR
    !
    !> No function briefing
    !
    integer function maxColumnsC(A) 
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        integer :: i
        maxColumnsC = 0
        do i = 1, A%nRow
        maxColumnsC = max(maxColumnsC, A%row(i+1)-A%row(i))
        enddo
        ! 
    end function maxColumnsC
    !
    !> No subroutine briefing
    !
    subroutine CSR2IJS_Real(C, S)
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: C
        type(spMatIJS_Real), intent( inout ) :: S
        integer :: ij, i, j
        ! for now no error checking
        if(.NOT.S%is_allocated) then
            stop "Error: CSR2IJS_Real > allocate output matrix before call"
        endif
        ij = 0
        do i=1, C%nRow
            do j = C%row(i), C%row(i+1)-1
                ij = ij + 1
                S%I(ij) = i
                S%J(ij) = C%col(j)
                S%S(ij) = C%val(j)
            enddo
        enddo
        !
    end subroutine CSR2IJS_Real
    !
    !> No subroutine briefing
    !
    subroutine CSR2IJS_Cmplx(C, S)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: C
        type(spMatIJS_Cmplx), intent( inout ) :: S
        integer :: ij, i, j
        ! for now no error checking
        if(.NOT.S%is_allocated) then
        stop "Error: CSR2IJS_Cmplx > allocate output matrix before call"
        endif
        ij = 0
        do i=1, C%nRow
        do j = C%row(i), C%row(i+1)-1
        ij = ij + 1
        S%I(ij) = i
        S%J(ij) = C%col(j)
        S%S(ij) = C%val(j)
        enddo
        enddo
        !
    end subroutine CSR2IJS_Cmplx
    !
    !> No subroutine briefing
    !
    subroutine IJS2CSR_Real( S, C )
        implicit none
        !
        type( spMatIJS_Real ), intent( in ) :: S
        type( spMatCSR_Real ), intent( inout ) :: C
        !
        integer :: i, j, nz
        integer, allocatable, dimension(:) :: rowT
        !
        allocate( rowT( S%nRow + 1 ) )
        !
        if(.NOT.C%is_allocated) then
            call errStop( "IJS2CSR_Real > allocate output matrix before call" )
        endif
        !
        !> first pass: find numbers of columns in each row of output
        rowT = 0
        nz = size(S%I)
        do i = 1, nz
            rowT(S%I(i)) = rowT(S%I(i))+1
        enddo
        !
        !> set row array in output CSR matrix
        C%row(1) = 1
        do i = 1, C%nRow
            C%row(i+1) = C%row(i)+rowT(i)
        enddo
        !
        !> now fill in columns and values
        rowT = 0
        do i = 1, nz
            !
            j = C%row( S%I(i) ) +rowT( S%I(i) )
            !
            C%col(j) = S%J(i)
            C%val(j) = S%S(i) 
            !
            rowT( S%I(i) ) = rowT( S%I(i) ) + 1
            !
        enddo
        !
    end subroutine IJS2CSR_Real
    !
    !> For now no error checking
    !> first should sort sparse matrix in IJS format, so that
    !> row numbers are strictly non-decreasing; come back to
    !> this, since spMatIJS converted from spMatCSR will
    !> already be ordered
    !
    subroutine IJS2CSR_Cmplx(S, C)
        implicit none
        !
        type(spMatIJS_Cmplx), intent( in ) :: S
        type( spMatCSR_Cmplx ), intent( inout ) :: C
        integer :: i, j, nz
        integer, allocatable, dimension(:) :: rowT

        allocate(rowT(S%nRow+1))

        if(.NOT.C%is_allocated) then
        stop "Error: IJS2CSR_Cmplx > allocate output matrix before call"
        endif

        !   first pass: find numbers of columns in each row of output
        rowT = 0
        nz = size(S%I)
        do i = 1, nz
        rowT(S%I(i)) = rowT(S%I(i))+1
        enddo
        !   set row array in output CSR matrix
        C%row(1) = 1
        do i = 1, C%nRow
        C%row(i+1) = C%row(i)+rowT(i)
        enddo

        !    now fill in columns and values
        rowT = 0
        do i = 1, nz
        j = C%row(S%I(i)) +rowT(S%I(i))
        C%col(j) = S%J(i)
        C%val(j) = S%S(i) 
        rowT(S%I(i)) = rowT(S%I(i))+1
        enddo
        deallocate(rowT)
        !
    end subroutine IJS2CSR_Cmplx
    !
    ! multiply a complex vector x by a real sparse CSR matrix A 
    !
    subroutine RMATxCVEC( A, x, y )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        complex( kind=prec ), dimension(:), intent( in ) :: x
        complex( kind=prec ), dimension(:), intent( inout ) :: y
        !
        integer :: i, j
        !
        !> lets start coding this with little checking -- assume
        !> everything is allocated and correct on entry
        !
        !write( *, * ) "*A%nCol, *size(x), size(y)", A%nCol, size(x), size(y)
        !
        if( A%nCol .NE. size(x) ) then
            write( *, * ) "Error: RMATxCVEC > matrix and vector sizes incompatible = ", A%nCol, size(x)
            stop
        endif
        !
        do i = 1, A%nRow
            y(i) = C_ZERO
            do j = A%row(i), A%row(i+1)-1 
                y(i) = y(i)+A%val(j)*x(A%col(j))
            enddo
        enddo
        !
    end subroutine RMATxCVEC
    !
    !> Multiply a complex vector x by a real sparse CSR matrix A 
    !
    subroutine RMATxRVEC(A, x, y)
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        real( kind=prec ), dimension(:), intent( in ) :: x
        real( kind=prec ), dimension(:), intent( inout ) :: y

        integer :: i, j

        ! lets start coding this with little checking -- assume
        ! everything is allocated and correct on entry

        if(A%nCol.NE.size(x)) then
        stop "Error: RMATxRVEC > matrix and vector sizes incompatible"
        endif

        do i = 1, A%nRow
        y(i) = 0.0
        do j = A%row(i), A%row(i+1)-1 
        y(i) = y(i)+A%val(j)*x(A%col(j))
        enddo
        enddo
        return
    end subroutine RMATxRVEC
    !
    !> Multiply a complex vector x by a complex sparse CSR matrix A 
    !
    subroutine CMATxCVEC(A, x, y)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        complex( kind=prec ), dimension(:), intent( in ) :: x
        complex( kind=prec ), dimension(:), intent( inout ) :: y

        integer :: i, j

        ! lets start coding this with little error checking -- assume
        ! everything is allocated and correct on entry

        if(A%nCol.NE.size(x)) then
        stop "Error: CMATxCVEC > matrix and vector sizes incompatible"
        endif

        do i = 1, A%nRow
        y(i) = 0.0
        do j = A%row(i), A%row(i+1)-1 
        y(i) = y(i)+A%val(j)*x(A%col(j))
        enddo
        enddo
        return
    end subroutine CMATxCVEC
    !
    !> matix-vector multiplication
    !
    subroutine RMATxRMAT( A, B, C )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A, B
        type( spMatCSR_Real ), intent( out ) :: C
        !
        type( spMatCSR_Real ) :: Ctmp
        integer :: i, j, k, nColMax, nCol, m, n, nz, nnz, jj, j1, j2, i1, i2, l, nzero
        integer, allocatable, dimension(:) :: colT, rowT
        logical :: new
        !
        allocate( rowT( A%nRow + 1 ) )
        nColMax = maxColumnsR(A) * maxColumnsR(B)
        allocate( colT( nColMax ) )
        !
        if( A%nCol .NE. B%nRow ) then
            stop "Error: RMATxRMAT > matrix sizes incompatible"
        endif
        !
        !> first pass: find numbers of columns in each row of output matrix C
        rowT(1) = 1
        !
        do i = 1, A%nRow
            nCol = 0
            do j = A%row(i), A%row(i+1)-1
                jj = A%col(j)
                do k = B%row(jj), B%row(jj+1)-1 
                    !
                    new = .TRUE.
                    do l = 1, nCol
                        new = new.and.(colT(l).NE.B%col(k))
                    enddo
                    !
                    if( new ) then
                        nCol = nCol+1
                        colT(nCol) = B%col(k)
                    endif
                enddo
            enddo
            rowT(i+1) = rowT(i)+nCol
        enddo
        !
        !> create output sparse matrix
        !> rowT is now row vector for output C
        m = A%nRow
        n = B%nCol
        nz = rowT(m+1)-1
        if(C%is_allocated) then
        call deall_spMatCSR(C)
        endif
        call create_spMatCSR(m, n, nz, Ctmp)
        !
        !> second pass: fill in columns and values of output matrix C
        !
        nzero = 0
        Ctmp%row = rowT
        deallocate(rowT)
        !
        do i = 1, A%nRow
            nCol = 0
            i1 = Ctmp%row(i)
            i2 = Ctmp%row(i+1)-1
            Ctmp%val(i1:i2) = 0.0d0
            do j = A%row(i), A%row(i+1)-1
                jj = A%col(j)
                do k = B%row(jj), B%row(jj+1)-1 
                    new = .TRUE.
                    do l = 1, nCol
                        if(colT(l).eq.B%col(k))then
                            new = .FALSE.
                            exit
                        endif
                    enddo
                    !
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
        !
        deallocate( colT )
        ! now try to clear zeros in Btmp
        nz = nz - nzero
        !
        call create_spMatCSR(m, n, nz, C)  
        !
        if( nzero .EQ. 0 ) then
            C%row=Ctmp%row
            C%col=Ctmp%col
            C%val=Ctmp%val
        else
            j1 = 1
            k  = 0
            nz = 1
            do i = 1, m
                nnz = 0
                j2 = Ctmp%row(i+1)-1
                do j = j1, j2
                    if(abs(Ctmp%val(j)).GT. 0) then
                        nnz = nnz + 1
                        k = k+1
                        C%col(k) = Ctmp%col(j)
                        C%val(k) = Ctmp%val(j)
                    endif
                enddo
                !
                j1 = j2+1
                C%row(i) = nz
                nz = nz + nnz
            enddo
            !
            C%row(m+1) = nz
            !
        endif
        !
        call deall_spMatCSR_Real(Ctmp)
        !
    end subroutine RMATxRMAT
    !
    !> Matix-vector multiplication, complex version
    !
    subroutine CMATxCMAT( A, B, C )
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A, B
        type( spMatCSR_Cmplx ), intent( out ) :: C

        integer :: i, j, k, nColMax, nCol, m, n, nz, jj, i1, i2, l
        integer, allocatable, dimension(:) :: colT
        integer, allocatable, dimension(:) :: rowT
        logical new

        if(A%nCol.NE.B%nRow) then
        stop "Error: CMATxCMAT > matrix sizes incompatible"
        endif

        allocate(rowT(A%nRow+1))
        nColMax = maxColumnsC(A)*maxColumnsC(B)
        allocate(colT(nColMax))

        !   first pass: find numbers of columns in each row of output
        !            matrix C
        rowT(1) = 1
        do i = 1, A%nRow
        nCol = 0
        do j = A%row(i), A%row(i+1)-1
        jj = A%col(j)
        do k = B%row(jj), B%row(jj+1)-1 
        new = .TRUE.
        do l = 1, nCol
        new = new.and.(colT(l).NE.B%col(k))
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
        call create_spMatCSR(m, n, nz, C)

        !   second pass: fill in columns and values of output
        !            matrix C

        C%row = rowT
        deallocate(rowT)
        do i = 1, A%nRow
        nCol = 0
        i1 = C%row(i)
        i2 = C%row(i+1)-1
        C%val(i1:i2) = 0.0d0
        do j = A%row(i), A%row(i+1)-1
        jj = A%col(j)
        do k = B%row(jj), B%row(jj+1)-1 
        new = .TRUE.
        do l = 1, nCol
        if(colT(l).eq.B%col(k))then
        new = .FALSE.
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
    end subroutine CMATxCMAT
    !
    !> Premultiply sparse matrix A by diagonal matrix D real version
    !
    subroutine DIAGxRMAT( D, A, B )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        real( kind=prec ), intent( in ), dimension(:) :: D
        type( spMatCSR_Real ), intent( out ) :: B
        !
        type( spMatCSR_Real ) :: Btmp
        !
        integer :: i, j, j1, j2, k, n, m, nz, nnz, nzero
        !
        if( A%nRow .NE. size(D) ) then
            call errStop( "DIAGxRMAT > matrix sizes incompatible" )
        endif
        !
        m = A%nRow
        n = A%nCol
        nz = A%row(A%nRow+1)-1
        if(.NOT.sameSizeCSR_Real(A, B)) then
        if(B%is_allocated) then
        call deall_spMatCSR_Real(B)
        endif
        call create_spMatCSR(m, n, nz, B)  
        endif
        call create_spMatCSR(m, n, nz, Btmp)  
        nzero=0
        Btmp%row(1) = 1
        do i = 1, A%nRow 
        Btmp%row(i+1) = A%row(i+1)
        do j = A%row(i), A%row(i+1)-1
        Btmp%val(j) = D(i)*A%val(j)
        Btmp%col(j) = A%col(j)
        if(Btmp%val(j).eq.0.0) then 
        nzero=nzero+1
        endif
        enddo
        enddo
        ! now try to clear zeros in Btmp
        nz = nz - nzero
        if(nzero.EQ.0) then
        B%row=Btmp%row
        B%col=Btmp%col
        B%val=Btmp%val
        else 
        j1 = 1
        k  = 0
        nz = 1
        do i = 1, m
        nnz = 0
        j2 = Btmp%row(i+1)-1
        do j = j1, j2
        if(abs(Btmp%val(j)).GT. 0) then
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
    end subroutine DIAGxRMAT
    !
    !> Premultiply sparse matrix A by diagonal matrix D complex version
    !
    subroutine DIAGxCMAT( D, A, B )
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        complex( kind=prec ), intent( in ), dimension(:) :: D
        type( spMatCSR_Cmplx ), intent( out ) :: B

        integer :: i, j, m, n, nz

        if(A%nRow.NE.size(D)) then
        stop "Error: DIAGxCMAT > matrix sizes incompatible"
        endif
        if(.NOT.sameSizeCSR_Cmplx(A, B)) then
        if(B%is_allocated) then
        call deall_spMatCSR(B)
        endif
        m = A%nRow
        n = A%nCol
        nz = A%row(A%nRow+1)-1
        call create_spMatCSR(m, n, nz, B)  
        endif

        B%row(1) = 1
        do i = 1, A%nRow 
        B%row(i+1) = A%row(i+1)
        do j = A%row(i), A%row(i+1)-1
        B%val(j) = D(i)*A%val(j)
        B%col(j) = A%col(j)
        enddo
        enddo
        return
    end subroutine DIAGxCMAT
    !
    !  postmultiply sparse matrix A by diagonal matrix D
    !        real version
    !
    subroutine RMATxDIAG( A, D, B )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        real( kind=prec ), intent( in ), dimension(:) :: D
        type( spMatCSR_Real ), intent( out ) :: B
        !
        type( spMatCSR_Real ) :: Btmp
        integer :: i, j, j1, j2, k, m, n, nz, nnz, nzero
        !
        if( A%nCol .NE. size( D ) ) then
            call errStop( "RMATxDIAG > matrix sizes incompatible" )
        endif
        !
        m = A%nRow
        n = A%nCol
        nz = A%row(A%nRow+1)-1
        if(.NOT.sameSizeCSR_Real(A, B)) then
            if(B%is_allocated) then
                call deall_spMatCSR_Real(B)
            endif
            call create_spMatCSR(m, n, nz, B)  
        endif
        !
        call create_spMatCSR(m, n, nz, Btmp)  
        nzero = 0
        Btmp%row(1) = 1
        do i = 1, A%nRow 
            Btmp%row(i+1) = A%row(i+1)
            do j = A%row(i), A%row(i+1)-1
                Btmp%val(j) = A%val(j)*D(A%col(j))
                Btmp%col(j) = A%col(j)
                if(Btmp%val(j).eq.0.0) then ! mark zero elements
                    nzero = nzero + 1
                endif
            enddo
        enddo
        ! now try to clear zeros in Btmp
        nz = nz - nzero
        !
        call create_spMatCSR(m, n, nz, B)
        !
        if( nzero .EQ. 0 ) then
            B%row=Btmp%row
            B%col=Btmp%col
            B%val=Btmp%val
        else 
            j1 = 1
            k  = 0
            nz = 1
            do i = 1, m
                nnz = 0
                j2 = Btmp%row(i+1)-1
                do j = j1, j2
                    if(abs(Btmp%val(j)).GT. 0) then
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
        !
        call deall_spMatCSR_Real(Btmp)
        !
    end subroutine RMATxDIAG
    !
    !> Postmultiply sparse matrix A by diagonal matrix D complex version
    !
    subroutine CMATxDIAG(A, D, B)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        real( kind=prec ), intent( in ), dimension(:) :: D
        type( spMatCSR_Cmplx ), intent( inout ) :: B

        integer :: i, j, m, n, nz

        if(A%nCol.NE.size(D)) then
        stop "Error: CMATxDIAG > matrix sizes incompatible"
        endif
        if(.NOT.sameSizeCSR_Cmplx(A, B)) then
        if(B%is_allocated) then
        call deall_spMatCSR(B)
        endif
        m = A%nRow
        n = A%nCol
        nz = A%row(A%nRow+1)-1
        call create_spMatCSR(m, n, nz, B)  
        endif

        B%row(1) = 1
        do i = 1, A%nRow 
        B%row(i+1) = A%row(i+1)
        do j = A%row(i), A%row(i+1)-1
        B%val(j) = A%val(j)*D(A%col(j))
        B%col(j) = A%col(j)
        enddo
        enddo
        !
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine CMATtrans( A, Atrans, Conj )
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        type( spMatCSR_Cmplx ), intent( out ) :: Atrans
        logical, intent( in ), optional :: Conj
        !
        type( spMatIJS_Cmplx ) :: B
        integer :: i, nz, temp
        logical :: conjugate
        !
        conjugate = .TRUE.
        !
        if(present(Conj)) then 
        conjugate = Conj
        endif
        nz = A%row(A%nRow+1)-1
        call create_spMatCSR(A%nCol, A%nRow, nz, Atrans)
        call create_spMatIJS(A%nRow, A%nCol, nz, B)
        call CSR2IJS(A, B)
        do i = 1, nz
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
        !
        call IJS2CSR(B, Atrans)
        call deall_spMATIJS(B)
        if(A%lower) then
        Atrans%upper = .TRUE.
        endif 
        if(A%upper) then
        Atrans%lower = .TRUE.
        endif 
        !
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine RMATtrans( A, Atrans )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        type( spMatCSR_Real ), intent( out ) :: Atrans
        !
        type( spMatIJS_Real ) :: B
        integer :: i, nz, nz1, temp
        !
        nz = A%row( A%nRow + 1 ) - 1
        !
        call create_spMatCSR( A%nCol, A%nRow, nz, Atrans )
        !
        call create_spMatIJS( A%nRow, A%nCol, nz, B )
        !
        call CSR2IJS( A, B )
        !
        do i = 1, nz
            temp = B%I(i)
            B%I(i) = B%J(i)
            B%J(i) = temp
        enddo
        !
        temp = B%nRow
        B%nRow = B%nCol
        B%nCol = temp
        !
        call IJS2CSR_Real( B, Atrans )
        !
        call deall_spMATIJS( B )
        !
        if(A%lower) then
            Atrans%upper = .TRUE.
        endif 
        !
        if(A%upper) then
            Atrans%lower = .TRUE.
        endif 
        !
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine write_CSR_real(fid, A)
        implicit none
        !
        integer, intent( in ) :: fid
        type( spMatCSR_Real ), intent( in ) :: A
        write(fid) A%nRow, A%nCol, A%row(A%nRow+1)-1
        write(fid) A%row
        write(fid) A%col
        write(fid) A%val
        !
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine read_CSR_real( fid, A )
        implicit none
        !
        integer, intent( in ) :: fid
        type( spMatCSR_Real ), intent( out ) :: A
        integer :: n, m, nz
        read(fid) n, m, nz
        if(A%is_allocated) then
        call deall_spMatCSR(A)
        endif
        call create_spMatCSR(n, m, nz, A)
        read(fid) A%row
        read(fid) A%col
        read(fid) A%val
        !
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine write_CSR_Cmplx(fid, A)
        implicit none
        !
        integer, intent( in ) :: fid
        type( spMatCSR_Cmplx ), intent( in ) :: A
        write(fid) A%nRow, A%nCol, A%row(A%nRow+1)-1
        write(fid) A%row
        write(fid) A%col
        write(fid) A%val
        !
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine read_CSR_Cmplx( fid, A )
        implicit none
        !
        integer, intent( in ) :: fid
        type( spMatCSR_Cmplx ), intent( out ) :: A
        integer :: n, m, nz
        read(fid) n, m, nz
        if(A%is_allocated) then
        call deall_spMatCSR(A)
        endif
        call create_spMatCSR(n, m, nz, A)
        read(fid) A%row
        read(fid) A%col
        read(fid) A%val
        !
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine write_IJS_real(fid, A)
        implicit none
        !
        integer, intent( in ) :: fid
        type(spMatIJS_Real), intent( in ) :: A
        integer :: nz

        nz = size(A%I)
        write(fid) A%nRow, A%nCol, nz
        write(fid) A%I
        write(fid) A%J
        write(fid) A%S
        return
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine write_IJS_Cmplx(fid, A)
        implicit none
        !
        integer, intent( in ) :: fid
        type(spMatIJS_Cmplx), intent( in ) :: A
        integer :: nz

        nz = size(A%I)
        write(fid) A%nRow, A%nCol, nz
        write(fid) A%I
        write(fid) A%J
        write(fid) A%S
        return
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine read_IJS_real(fid, A)
        implicit none
        !
        integer, intent( in ) :: fid
        type(spMatIJS_Real), intent( out ) :: A
        integer :: n, m, nz
        read(fid) n, m, nz
        if(A%is_allocated) then
        call deall_spMatIJS(A)
        endif
        call create_spMatIJS(n, m, nz, A)
        read(fid) A%I
        read(fid) A%J
        read(fid) A%S
        !
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine read_IJS_Cmplx(fid, A)
        implicit none
        !
        integer, intent( in ) :: fid
        type(spMatIJS_Cmplx), intent( out ) :: A
        integer :: n, m, nz
        read(fid) n, m, nz
        if(A%is_allocated) then
        call deall_spMatIJS(A)
        endif
        call create_spMatIJS(n, m, nz, A)
        read(fid) A%I
        read(fid) A%J
        read(fid) A%S
        !
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine write_CSRasIJS_Real(A,fid)
        implicit none
        !
        integer, intent( in ) :: fid
        type( spMatCSR_Real ), intent( in ) :: A

        type(spMatIJS_Real) :: B
        integer :: n, m, nz

        m = A%nRow
        n = A%nCol
        nz = A%row(m+1)-1
        call create_spMatIJS(m, n, nz, B)
        call CSR2IJS(A, B)
        call write_IJS_real(fid, B)
        close(fid)
        call deall_spMatIJS(B)
        !
    end subroutine
    !
    !> No subroutine briefing
    !
    subroutine write_CSRasIJS_Cmplx(A,fid)
        implicit none
        !
        integer, intent( in ) :: fid
        type( spMatCSR_Cmplx ), intent( in ) :: A

        type(spMatIJS_Cmplx) :: B
        integer :: n, m, nz

        m = A%nRow
        n = A%nCol
        nz = A%row(m+1)-1
        call create_spMatIJS(m, n, nz, B)
        call CSR2IJS(A, B)
        call write_IJS_Cmplx(fid, B)
        close(fid)
        call deall_spMatIJS(B)
        return
    end subroutine
    !
    !> Extract lower triangular part of matrix A in CSR storage
    !
    subroutine lowerTri_Real( A, L )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        type( spMatCSR_Real ), intent( out ) :: L
        integer :: kk, n, m, i, j, nz
        integer, allocatable, dimension(:) :: rowT

        m = A%nRow
        n = A%nCol
        allocate(rowT(m))

        !   first pass: find numbers of columns in each row of output
        rowT = 0
        do i = 1, m
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).le.i) then
        rowT(i) = rowT(i)+1
        endif
        enddo
        enddo
        nz = 0
        do i = 1, m
        nz = nz+rowT(i)
        enddo
        call create_spMatCSR(m, n, nz, L)
        L%lower = .TRUE.

        !   set row array in output CSR matrix
        L%row(1) = 1
        do i = 1, L%nRow
        L%row(i+1) = L%row(i)+rowT(i)
        enddo
        deallocate(rowT)

        !    now fill in columns and values
        do i = 1, m
        kk = 0
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).le.i) then
        L%col(L%row(i)+kk) = A%col(j) 
        L%val(L%row(i)+kk) = A%val(j) 
        kk = kk+1
        endif
        enddo
        enddo
        return
    end subroutine
    !
    !> Extract lower triangular part of matrix A in CSR storage
    !
    subroutine upperTri_Real( A, U )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        type( spMatCSR_Real ), intent( out ) :: U
        integer :: kk, n, m, i, j, nz
        integer, allocatable, dimension(:) :: rowT

        m = A%nRow
        n = A%nCol
        allocate(rowT(m))

        !   first pass: find numbers of columns in each row of output
        rowT = 0
        do i = 1, m
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).ge.i) then
        rowT(i) = rowT(i)+1
        endif
        enddo
        enddo
        nz = 0
        do i = 1, m
        nz = nz+rowT(i)
        enddo
        call create_spMatCSR(m, n, nz, U)
        U%upper = .TRUE.

        !   set row array in output CSR matrix
        U%row(1) = 1
        do i = 1, U%nRow
        U%row(i+1) = U%row(i)+rowT(i)
        enddo
        deallocate(rowT)

        !    now fill in columns and values
        do i = 1, m
        kk = 0
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).ge.i) then
        U%col(U%row(i)+kk) = A%col(j) 
        U%val(U%row(i)+kk) = A%val(j) 
        kk = kk+1
        endif
        enddo
        enddo
        return
    end subroutine
    !
    !> Extract lower triangular part of matrix A in CSR storage
    !
    subroutine lowerTri_Cmplx( A, L )
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        type( spMatCSR_Cmplx ), intent( out ) :: L
        integer :: kk, n, m, i, j, nz
        integer, allocatable, dimension(:) :: rowT

        m = A%nRow
        n = A%nCol
        allocate(rowT(m))

        !   first pass: find numbers of columns in each row of output
        rowT = 0
        do i = 1, m
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).le.i) then
        rowT(i) = rowT(i)+1
        endif
        enddo
        enddo
        nz = 0
        do i = 1, m
        nz = nz+rowT(i)
        enddo
        call create_spMatCSR(m, n, nz, L)
        L%lower = .TRUE.

        !   set row array in output CSR matrix
        L%row(1) = 1
        do i = 1, L%nRow
        L%row(i+1) = L%row(i)+rowT(i)
        enddo
        deallocate(rowT)

        !    now fill in columns and values
        do i = 1, m
        kk = 0
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).le.i) then
        L%col(L%row(i)+kk) = A%col(j) 
        L%val(L%row(i)+kk) = A%val(j) 
        kk = kk+1
        endif
        enddo
        enddo
        return
    end subroutine
    !
    !> Extract lower triangular part of matrix A in CSR storage
    !
    subroutine upperTri_Cmplx( A, U )
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        type( spMatCSR_Cmplx ), intent( out ) :: U
        integer :: kk, n, m, i, j, nz
        integer, allocatable, dimension(:) :: rowT

        m = A%nRow
        n = A%nCol
        allocate(rowT(m))

        !   first pass: find numbers of columns in each row of output
        rowT = 0
        do i = 1, m
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).ge.i) then
        rowT(i) = rowT(i)+1
        endif
        enddo
        enddo
        nz = 0
        do i = 1, m
        nz = nz+rowT(i)
        enddo
        call create_spMatCSR(m, n, nz, U)
        U%upper = .TRUE.

        !   set row array in output CSR matrix
        U%row(1) = 1
        do i = 1, U%nRow
        U%row(i+1) = U%row(i)+rowT(i)
        enddo
        deallocate(rowT)

        !    now fill in columns and values
        do i = 1, m
        kk = 0
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).ge.i) then
        U%col(U%row(i)+kk) = A%col(j) 
        U%val(U%row(i)+kk) = A%val(j) 
        kk = kk+1
        endif
        enddo
        enddo
        return
    end subroutine
    !
    !> Extract diagonal part of matrix A in CSR storage
    !
    subroutine diag_Real(A, D)
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        real( kind=prec ), allocatable, intent( inout ) :: D(:)
        integer :: n, m, i, j

        m = A%nRow
        n = A%nCol
        if(n.NE.m) then
        stop "Error: diag_Real > diag only works for square matrices"
        endif
        if(allocated(D)) then
        deallocate(D)
        endif
        allocate(D(m))
        D = 0
        do i = 1, m
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).eq.i) then
        D(i) = A%val(j)
        endif
        enddo
        enddo
    end subroutine diag_Real
    !
    !> Extract diagonal part of matrix A in CSR storage
    !
    subroutine diag_Cmplx(A, D)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        complex( kind=prec ), allocatable, intent( inout ) :: D(:)
        integer :: n, m, i, j

        m = A%nRow
        n = A%nCol
        if(n.NE.m) then
        stop "Error: diag_Cmplx > diag only works for square matrices"
        endif
        if(allocated(D)) then
        deallocate(D)
        endif
        allocate(D(m))
        D = 0
        do i = 1, m
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).eq.i) then
        D(i) = A%val(j)
        endif
        enddo
        enddo
    end subroutine diag_Cmplx
    !
    !> Extract submatrix of A with rows and columns given by integer arrays r and c 
    !
    subroutine subMatrix_Real( A, r, c, B )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        type( spMatCSR_Real ), intent( out ) :: B
        integer, intent( in ) :: r(:), c(:)
        !
        integer :: kk, n, m, i, j, nz, k
        integer, allocatable, dimension(:) :: rowT, colT
        !
        m = size(r)
        n = size(c)
        allocate(rowT(m+1))
        allocate(colT(A%nCol))
        colT = 0
        !
        do i = 1, n
            colT(c(i)) = i
        enddo
        !
        !> count number of entries in each row
        rowT = 0
        nz = 0
        do i =1, m
            do j = A%row( r(i) ), A%row( r(i) + 1 ) - 1
                if( colT( A%col(j) ) .GT. 0 ) then
                    rowT(i) = rowT(i)+1
                    nz = nz+1
                endif
            enddo
        enddo
        !
        call create_spMatCSR( m, n, nz, B )
        !
        !   set row array in output CSR matrix
        B%row(1) = 1
        do i = 1, B%nRow
            B%row(i+1) = B%row(i)+rowT(i)
        enddo
        !
        do i =1, m
            kk = 0
            do j = A%row(r(i)), A%row(r(i)+1)-1
                if(colT(A%col(j)).GT.0) then
                    B%col(B%row(i)+kk) = colT(A%col(j))
                    B%val(B%row(i)+kk) = A%val(j)
                    kk = kk+1
                endif
            enddo
        enddo
        !
    end subroutine subMatrix_Real
    !
    !> Extract submatrix of A with rows and columns given by integer arrays r and c 
    !
    subroutine subMatrix_Cmplx(A, r, c, B)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        integer, intent( in ) :: r(:), c(:)
        type( spMatCSR_Cmplx ), intent( out ) :: B
        !
        integer :: kk, n, m, i, j, nz, k
        integer, allocatable, dimension(:) :: rowT, colT
        !
        m = size(r)
        n = size(c)
        allocate(rowT(m+1))
        allocate(colT(A%nCol))
        colT = 0
        do i = 1, n
        colT(c(i)) = i
        enddo
        !
        !   count number of entries in each row
        rowT = 0
        nz = 0
        do i =1, m
        do j = A%row(r(i)), A%row(r(i)+1)-1
        if(colT(A%col(j)).GT.0) then
        rowT(i) = rowT(i)+1
        nz = nz+1
        endif
        enddo
        enddo
        !
        call create_spMatCSR(m, n, nz, B)
        !
        !   set row array in output CSR matrix
        B%row(1) = 1
        do i = 1, B%nRow
        B%row(i+1) = B%row(i)+rowT(i)
        enddo
        !
        do i =1, m
        kk = 0
        do j = A%row(r(i)), A%row(r(i)+1)-1
        if(colT(A%col(j)).GT.0) then
        B%col(B%row(i)+kk) = colT(A%col(j))
        B%val(B%row(i)+kk) = A%val(j)
        kk = kk+1
        endif
        enddo
        enddo
        !
    end subroutine subMatrix_Cmplx
    !
    !> Merge an array of sparse matrices in CSR storage
    !> into a single block diagonal matrix B
    !
    subroutine BlkDiag_Real(A, B)
        implicit none
        !
        type( spMatCSR_Real ), pointer, intent( in ) :: A(:)
        type( spMatCSR_Real ), intent( out ) :: B

        integer :: nBlks, nRowB, nColB, nzB, i, i1, j1, i2, j2, k1

        nBlks = size(A) 
        nRowB = 0
        nColB = 0
        nzB = 0
        do i = 1, nBlks
        nRowB = nRowB + A(i)%nRow
        nColB = nColB + A(i)%nCol
        nzB = nzB + A(i)%row(A(i)%nRow+1)-1
        enddo
        call create_spMatCSR(nRowB, nColB, nzB, B)
        i1 = 1
        j1 = 1
        k1 = 0
        B%upper = .TRUE.
        B%lower = .TRUE.
        do i = 1, nBlks
        B%upper = B%upper .AND. A(i)%upper
        B%lower = B%lower .AND. A(i)%lower
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
        !
    end subroutine BlkDiag_Real
    !
    !> Merge an array of sparse matrices in CSR storage
    !> into a single block diagonal matrix B
    !
    subroutine BlkDiag_Cmplx(A, B)
        implicit none
        !
        type( spMatCSR_Cmplx ), pointer, intent( in ) :: A(:)
        type( spMatCSR_Cmplx ), intent( out ) :: B

        integer :: nBlks, nRowB, nColB, nzB, i, i1, j1, i2, j2, k1

        nBlks = size(A) 
        nRowB = 0
        nColB = 0
        nzB = 0
        do i = 1, nBlks
        nRowB = nRowB + A(i)%nRow
        nColB = nColB + A(i)%nCol
        nzB = nzB + A(i)%row(A(i)%nRow+1)-1
        enddo
        call create_spMatCSR(nRowB, nColB, nzB, B)
        i1 = 1
        j1 = 1
        k1 = 0
        B%upper = .TRUE.
        B%lower = .TRUE.
        do i = 1, nBlks
        B%upper = B%upper .AND. A(i)%upper
        B%lower = B%lower .AND. A(i)%lower
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
        !
    end subroutine BlkDiag_Cmplx
    !
    !   convert real CSR matrix to complex form
    !
    subroutine R2C_CSR(Ar, Ac)
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: Ar
        type( spMatCSR_Cmplx ), intent( out ) :: Ac
        integer :: nRow, nCol, nz, i

        nRow = Ar%nRow
        nCol = Ar%nCol
        nz = Ar%row(nRow+1)-1
        call create_spMatCSR(nRow, nCol, nz, Ac)
        Ac%row = Ar%row
        Ac%col = Ar%col
        do i = 1, nz
        Ac%val(i) = cmplx(Ar%val(i), 0)
        enddo
    end subroutine R2C_CSR
    !
    !> Automaticly split a matrix A into row submatrices 
    !> and take only the corresponding row submatrice of B
    !> for the(i+1)th process in n parallel threads
    !> this is used to prepare PETSc AIJ type matrix
    !> 
    !> Note: this can be easily changed to be split according to the
    !> number of non-zero elements in each submatrix important!
    !> the submatrix B is modified to use zero-based index(as Petsc)
    !
    subroutine splitRMAT(A, i, np, B, isizes)
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A  ! original matrix
        type( spMatCSR_Real ), intent( out ) :: B  ! submatrix
        integer, intent( in ) :: i, np
        integer, intent( in ), pointer, dimension(:), optional :: isizes
        !
        integer :: istart, iend, nrow_l, j, k
        integer :: m, n, nz, nz_l, nsub
        real :: nrow
        integer, allocatable, dimension(:) :: rowT, colT
        !
        allocate(colT(A%nCol))
        colT =(/(j, j=1, A%nCol) /)
        if(A%nrow .LT. np) then
            stop "Error: splitRMAT > number of process is larger than number of rows!"
        elseif(np.EQ.1) then
            !write( *, * ) "only one process, returning the original Matrix"
            m = A%nRow
            n = A%nCol
            nz = A%row(m+1)-1
            call create_spMatCSR_Real(m, n, nz, B)
            B%nRow=m
            B%nCol=n
            B%row=A%row-1
            B%col=A%col-1
            B%val=A%val
            return
        endif
        nrow = A%nRow
        !nz_l = floor((A%row(nrow+1)-1)/np)
        if(present(isizes)) then ! split into given sizes
        istart = 1
        do k = 1, i 
        istart = istart+ isizes(k)
        enddo
        iend = istart+isizes(k)-1
        else
        nrow_l = floor(nrow/np)
        istart=i*nrow_l+1
        if(i+1 .EQ. np)  then! the last sub-matrix
        iend=nrow
        else
        iend=i*nrow_l+nrow_l
        endif
        endif
        allocate(rowT(iend-istart+1))
        rowT =(/(j, j=istart, iend) /)
        call subMatrix_Real(A, rowT, colT, B)
        B%row=B%row-1
        B%col=B%col-1
        deallocate(colT)
        deallocate(rowT)
        !
    end subroutine splitRMAT
    !
    !> Automaticly split a matrix A into row submatrices 
    !> and take only the corresponding row submatrice of B
    !> for the(i+1)th process in n parallel threads
    !> this is used to prepare PETSc AIJ type matrix
    !> 
    !> Note: this can be easily changed to be split according to the
    !> number of none zeros elemtents in each submatrix 
    !> the submatrix B is modified to use zero-based index(as Petsc)
    !
    subroutine splitCMAT(A, i, np, B, isizes)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A  ! original matrix
        type( spMatCSR_Cmplx ), intent( out ) :: B  ! submatrix
        integer, intent( in ) :: i, np
        integer, intent( in ), pointer, dimension(:), optional :: isizes
        integer :: istart, iend, nrow_l, j, k
        integer :: m, n, nz
        real :: nrow
        integer, allocatable, dimension(:) :: rowT, colT
        !
        allocate(colT(A%nCol))
        colT =(/(j, j=1, A%nCol) /)
        if(A%nrow .LT. np) then
            stop "Error: splitCMAT > number of processes is larger than number of rows!"
        elseif(np.EQ.1) then
            !write( *, * ) "only one process, returning the original Matrix"
            m = A%nRow
            n = A%nCol
            nz = A%row(m+1)-1
            call create_spMatCSR_Cmplx(m, n, nz, B)
            B%nRow=m
            B%nCol=n
            B%row=A%row-1
            B%col=A%col-1
            B%val=A%val
            return
        endif
        nrow = A%nRow
        !nz_l = floor((A%row(nrow+1)-1)/np)
        if(present(isizes)) then ! split according to the given sizes
        istart = 1
        do k = 1, i 
        istart = istart+ isizes(k)
        enddo
        iend = istart+isizes(k)-1
        else ! split evenly 
        nrow_l = floor(nrow/np)
        istart=i*nrow_l+1
        if(i+1 .EQ. np)  then! the last sub-matrix
        iend=nrow
        else
        iend=i*nrow_l+nrow_l
        endif
        endif
        allocate(rowT(iend-istart+1))
        rowT =(/(j, j=istart, iend) /)
        call subMatrix_Cmplx(A, rowT, colT, B)
        B%row=B%row-1
        B%col=B%col-1
        deallocate(colT)
        deallocate(rowT)
        !
    end subroutine splitCMAT
    !
    !> Solve system Lx = b for complex vector x, lower triangular L
    !> here real or cmplx refers to U; x is always complex
    !
    subroutine LTsolve_Cmplx( L, b, x )
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: L
        complex( kind=prec ), dimension(:), intent( in ) :: b
        complex( kind=prec ), dimension(:), intent( inout ) :: x
        !
        complex( kind=prec ) :: d
        integer :: i, j
        !
        if( L%nRow .NE. L%nCol ) then
            stop "Error: LTsolve_Cmplx > sparse matrix must be square"
        endif
        !
        if( .NOT. L%lower ) then
            stop "Error: LTsolve_Cmplx > sparse matrix must be lower triangular"
        endif 
        !
        if( size(x) .NE. L%nRow ) then
            stop "Error: LTsolve_Cmplx > output vector x not of correct size"
        endif
        !
        do i = 1, L%nRow
            x(i) = b(i)
            do j = L%row(i), L%row(i+1)-1
                if(L%col(j).LT.i) then
                    x(i) = x(i)-L%val(j)*x(L%col(j))
                else
                    !   in this case L%col(j) = i
                    d = L%val(j)
                endif
            enddo
            x(i) = x(i)/d
        enddo
        !
    end subroutine LTsolve_Cmplx
    !
    !> Solve system Ux = b for complex vector x, upper triangular U
    !> here real or cmplx refers to U; x is always complex
    !
    subroutine UTsolve_Cmplx(U, b, x)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: U
        complex( kind=prec ), dimension(:), intent( in ) :: b
        complex( kind=prec ), dimension(:), intent( inout ) :: x

        integer :: i, j
        complex( kind=prec ) :: d

        if(U%nRow .NE.U%nCol) then
        stop "Error: UTsolve_Cmplx > sparse matrix must be square"
        endif 
        if(.NOT.U%upper) then
        stop "Error: UTsolve_Cmplx > sparse matrix must be upper triangular"
        endif 
        if(size(x).NE.U%nRow) then
        stop "Error: UTsolve_Cmplx > output vector x not of correct size"
        endif 
        do i = U%nRow, 1, -1
        x(i) = b(i)
        do j = U%row(i), U%row(i+1)-1
        if(U%col(j).GT.i) then
        x(i) = x(i)-U%val(j)*x(U%col(j))
        else
        !   in this case U%col(j) = i
        d = U%val(j)
        endif
        enddo
        x(i) = x(i)/d
        enddo
        return
    end subroutine UTsolve_Cmplx
    !
    !> Solve system Lx = b for complex vector x, lower triangular L
    !> here real or cmplx refers to L; x is always complex
    !
    subroutine LTsolve_Real( L, b, x )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: L
        complex( kind=prec ), dimension(:), intent( in ) :: b
        complex( kind=prec ), dimension(:), intent( inout ) :: x

        real( kind=prec ) :: d
        integer :: i, j

        if(L%nRow .NE.L%nCol) then
        stop "Error: LTsolve_Real > sparse matrix must be square"
        endif 
        if(.NOT.L%lower) then
        stop "Error: LTsolve_Real > sparse matrix must be lower triangular"
        endif 
        if(size(x).NE.L%nRow) then
        stop "Error: LTsolve_Real > output vector x not of correct size"
        endif 
        do i = 1, L%nRow
        x(i) = b(i)
        do j = L%row(i), L%row(i+1)-1
        if(L%col(j).LT.i) then
        x(i) = x(i)-L%val(j)*x(L%col(j))
        else
        !   in this case L%col(j) = i
        d = L%val(j)
        endif
        enddo
        x(i) = x(i)/d
        enddo
        return 
    end subroutine LTsolve_Real
    !
    !> Solve system Ux = b for complex vector x, upper triangular U
    !> ere real or cmplx refers to L; x is always complex
    !
    subroutine UTsolve_Real( U, b, x )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: U
        complex( kind=prec ), dimension(:), intent( in ) :: b
        complex( kind=prec ), dimension(:), intent( inout ) :: x
        !
        integer :: i, j
        real( kind=prec ) :: d
        !
        if( U%nRow .NE. U%nCol ) then
            stop "Error: UTsolve_Real > sparse matrix must be square"
        endif
        !
        if( .NOT. U%upper ) then
            stop "Error: UTsolve_Real > sparse matrix must be upper triangular"
        endif
        !
        if( size(x) .NE. U%nRow ) then
            stop "Error: UTsolve_Real > output vector x not of correct size"
        endif
        !
        do i = U%nRow, 1, -1
            x(i) = b(i)
            do j = U%row(i), U%row(i+1)-1
                if(U%col(j).GT.i) then
                    x(i) = x(i)-U%val(j)*x(U%col(j))
                else
                    !   in this case U%col(j) = i
                    d = U%val(j)
                endif
            enddo
            x(i) = x(i)/d
        enddo
        !
    end subroutine UTsolve_Real
    !
    !> Specialized routine to add imaginary D to diagonal of real sparse
    !> matrix A, saving output to complex sparse matrix B
    !
    subroutine CSR_R2Cdiag( A, d, B )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        real( kind=prec ), dimension(:), intent( in ) :: d
        type( spMatCSR_Cmplx ), intent( out ) :: B
        !
        integer :: n, m, nz, i, j
        !
        m = A%nRow
        n = A%nCol
        nz = A%row(m+1)-1
        call create_spMatCSR_Cmplx( m, n, nz, B )
        B%row  = A%row
        B%col = A%col
        B%val = A%val
        do i =1, m
            do j =B%row(i), B%row(i+1)-1
                !   not sure ISIGN should be here!
                if(B%col(j) .EQ. i) then
                    B%val(j) = B%val(j)+ISIGN*CMPLX(0.0, 1.0, 8)*d(i)
                    exit
                endif
            enddo
        enddo
        !
    end subroutine CSR_R2Cdiag
    !
    !> Incomplete cholesky decomposition of a symmetric matrix A
    !> This assumes the matrix is symmetric, but does not check
    !> ALSO: this should only be used with positive definite matrices
    !> Could act up if matrix is not symmetric pos-def
    !
    subroutine CholInc_Real(A, L)
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        type( spMatCSR_Real ), intent( inout ) :: L
        integer :: ii, mMax, n, m, i, j, nij, nji, k, k1, nij1, nji1, i1, j1
        integer, allocatable, dimension(:) :: ij, ji, ji1
        !
        m = A%nRow
        n = A%nCol
        mMax = maxColumnsR(A)
        allocate(ij(mMax))
        allocate(ji(mMax))
        allocate(ji1(mMax))
        call lowerTri_Real(A, L)
        do i = 1, m ! sweep in rows
        !  divide up columns in row i: diagonal, to left, to right
        nij = 0
        nji = 0
        do j = A%row(i), A%row(i+1)-1 !sweep in columns
        if(i.EQ.A%col(j)) then
        !   this index should correspond to positions where
        !   columns and values are stored in L
        ii = j - A%row(i) + L%row(i)
        elseif(i.LT.A%col(j)) then
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
        do k = 1, nji
        L%val(ii) = L%val(ii) - L%val(ji(k))*L%val(ji(k))
        enddo
        if(L%val(ii).LT.0) then
        L%val(ii) = -L%val(ii)
        endif
        L%val(ii) = sqrt(L%val(ii)) ! diagonal
        do j = 1, nij
        !   these are rows of L that have elements in column i
        i1 = A%col(ij(j))
        nji1 = 0
        !  divide up elements in this row: diagonal, to left
        do k = L%row(i1), L%row(i1+1)-1
        if(i.EQ.L%col(k)) then
        j1 = k
        else
        nji1 = nji1+1 
        ji1(nji1) = k
        endif
        enddo
        do k = 1, nji
        do k1  = 1, nji1
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
    end subroutine CholInc_Real
    !
    !> This mimics approach used in ModEM --- D-ILU
    !> NOT ILU-0
    !> THIS ASSUMES THE MATRIX IS SYMMETRIC -- but not necessarily Hermitian
    !
    subroutine Dilu_Real(A, L, U)
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        type( spMatCSR_Real ), intent( inout ) :: L, U

        real( kind=prec ), allocatable, dimension(:) ::d
        integer ::n, m, nz, i, j

        call lowerTri(A, L)
        call upperTri(A, U)
        m = A%nRow
        allocate(d(m))
        do i = 1, m
        d(i) = 0.0_dp
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).eq.i) then
        d(i) = d(i) + A%val(j)
        elseif(A%col(j).LT.i) then
        d(i) = d(i) - A%val(j)*A%val(j)*d(A%col(j))
        endif
        enddo
        d(i) = 1.0_dp/d(i)
        enddo
        do i = 1, m
        do j = L%row(i), L%row(i+1)-1
        if(L%col(j).eq.i) then
        L%val(j) = 1
        else
        L%val(j) = L%val(j)*d(L%col(j))
        endif
        enddo
        enddo
        do i = 1, m
        do j = U%row(i), U%row(i+1)-1
        if(U%col(j).eq.i) then
        U%val(j) = 1.0_dp/d(i)
        exit
        endif
        enddo
        enddo
        return
    end subroutine Dilu_Real
    !
    !> This mimics approach used in ModEM --- D-ILU(diagonal-ILU)
    !> NOT ILU-0
    !> THIS ASSUMES THE MATRIX IS SYMMETRIC -- but not necessarily Hermitian
    !> so it will not work if using modified system equation(as Randy did)
    !
    subroutine Dilu_Cmplx(A, L, U)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        type( spMatCSR_Cmplx ), intent( inout ) :: L, U

        complex( kind=prec ), allocatable, dimension(:) ::d
        integer ::n, m, nz, i, j

        call lowerTri(A, L)
        call upperTri(A, U)
        m = A%nRow
        allocate(d(m))
        do i = 1, m
        d(i) = CMPLX(0.0, 0.0, 8)
        do j = A%row(i), A%row(i+1)-1
        if(A%col(j).eq.i) then
        d(i) = d(i) + A%val(j)
        elseif(A%col(j).LT.i) then
        d(i) = d(i) - A%val(j)*A%val(j)*d(A%col(j))
        endif
        enddo
        d(i) = CMPLX(1.0, 0.0, 8)/d(i)
        enddo
        do i = 1, m
        do j = L%row(i), L%row(i+1)-1
        if(L%col(j).eq.i) then
        L%val(j) = 1
        else
        L%val(j) = L%val(j)*d(L%col(j))
        endif
        enddo
        enddo
        do i = 1, m
        do j = U%row(i), U%row(i+1)-1
        if(U%col(j).eq.i) then
        U%val(j) = 1.0_dp/d(i)
        exit
        endif
        enddo
        enddo
        return
        !
    end subroutine
    !
    !> This mimics approach used in ModEM --- D-ILU(diagonal-ILU)
    !> NOT ILU-0
    !> slightly modified to be used with modified system equations
    !
    subroutine Dilu_Cmplx_AS(A, L, U)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        type( spMatCSR_Cmplx ), intent( inout ) :: L, U
        type( spMatCSR_Cmplx ) :: UT

        complex( kind=prec ), allocatable, dimension(:) ::d
        integer ::n, m, nz, i, j

        call lowerTri(A, L)
        call upperTri(A, U)
        call CMATtrans(U, UT)
        m = L%nRow
        allocate(d(m))
        do i = 1, m ! loop through rows
        d(i) = C_ZERO
        do j = L%row(i), L%row(i+1)-1 !loop through columns
        if(L%col(j).eq.i) then ! diagonal
        d(i) = d(i) + L%val(j) 
        elseif(L%col(j).LT.i) then ! take the L and U side
        d(i) = d(i) - L%val(j)*UT%val(j)*d(L%col(j))
        endif
        enddo
        d(i) = C_ONE/d(i)
        enddo
        call deall_spMatCSR(UT)
        do i = 1, m
        do j = L%row(i), L%row(i+1)-1
        if(L%col(j).eq.i) then
        L%val(j) = 1
        else
        L%val(j) = L%val(j)*d(L%col(j))
        endif
        enddo
        enddo
        do i = 1, m
        do j = U%row(i), U%row(i+1)-1
        if(U%col(j).eq.i) then
        U%val(j) = 1.0_dp/d(i)
        exit
        endif
        enddo
        enddo
        return
        !
    end subroutine Dilu_Cmplx_AS
    !
    !> Simple but not at all intuitive ILU0 routine with CSR sparse matrix 
    !> 
    !> the algorithm is modified from a much fancier version in Yousef Saad"s 
    !> 2003 book <Iterative Methods for Sparse Linear Systems>, Chapter 10.2
    !> 
    !> apparently the original algorithm seeks to minimize the memory 
    !> consumption by storing the L and U in the original sparse matrix 
    !> structure of A(as ILU0 does not have any fill-ins).
    !
    subroutine ilu0_Cmplx(A, L, U)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( in ) :: A
        type( spMatCSR_Cmplx ), intent( inout ) :: L, U
        type( spMatCSR_Cmplx ) :: Atmp
        complex( kind=prec ), allocatable, dimension(:) ::d
        complex( kind=prec ) :: piv
        integer :: n, m, nz, i, j, j2, k, p
        !
        allocate(d(A%nRow))
        n = A%nRow
        m = A%nCol
        nz = A%row(A%nRow+1)-1
        call create_spMatCSR(n, m, nz, Atmp)
        Atmp%row=A%row
        Atmp%col=A%col
        Atmp%val=A%val
        call sort_spMatCSR(Atmp) ! sort the col indices in A
        d = C_ZERO
        do i=1, Atmp%nRow ! loop through rows
            !
            p=Atmp%row(i) ! mark the first none zero element in current row
            !
            do j=Atmp%row(i), Atmp%row(i+1)-1 ! loop through columns
                !
                if(Atmp%col(j).eq.i) then !diagonal
                    d(i) = Atmp%val(j) ! store previous diagonal elements
                    exit ! exit as we reached the last element in L
                elseif(Atmp%col(j).LT.i) then
                    !
                    if(d(Atmp%col(j)).eq.C_ZERO) then
                        write( *, * ) "Error: ilu0_Cmplx > zero pivoting in ILU0 "
                        write( *, * ) "in col: ", Atmp%col(Atmp%row(i):Atmp%row(i+1)-1)
                        write( *, * ) "in row: ", Atmp%col(j)
                        stop
                    endif 
                    ! first divide each row in L with diagonal elements
                    Atmp%val(j) = Atmp%val(j)/d(Atmp%col(j))
                    piv = Atmp%val(j)
                    do k = p+1, Atmp%row(i+1)-1 ! then adding up each column in U
                        do j2 = Atmp%row(Atmp%col(j)), Atmp%row(Atmp%col(j)+1)-1
                            if(Atmp%col(j2).eq.Atmp%col(k)) then
                                Atmp%val(k) = Atmp%val(k) - piv* Atmp%val(j2)
                                exit
                            endif
                        enddo
                    enddo
                    p = p + 1
                    !
                endif
                !
            enddo
            !
        enddo
        deallocate(d)
        call lowerTri(Atmp, L)
        do i = 2, L%nRow+1 ! set the diagonal element of L to 1
        ! the diagonal element should be the last in each row 
        ! ONLY IF col is properly sorted
        L%val(L%row(i)-1) = C_ONE
        enddo
        call upperTri(Atmp, U)
        call deall_spMatCSR(Atmp)
        return
    end subroutine ilu0_Cmplx
    !
    !> A simple but not at all intuitive ILU0 routine with CSR sparse matrix 
    !
    !> The algorithm is modified from a much fancier version in Yousef Saad"s 
    !> 2003 book <Iterative Methods for Sparse Linear Systems>, Chapter 10.2
    !> 
    !> apparently the original algorithm seeks to minimize the memory 
    !> consumption by storing the L and U in the original sparse matrix 
    !> structure of A(as ILU0 does not have any fill-ins).
    !
    subroutine ilu0_Real(A, L, U)
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A
        type( spMatCSR_Real ), intent( inout ) :: L, U
        type( spMatCSR_Real ) :: Atmp
        real( kind=prec ), allocatable, dimension(:):: d
        real( kind=prec ) :: piv
        integer :: n, m, nz, i, j, j2, k, p

        allocate(d(A%nRow))
        n = A%nRow
        m = A%nCol
        nz = A%row(A%nRow+1)-1
        call create_spMatCSR(n, m, nz, Atmp)
        Atmp%row=A%row
        Atmp%col=A%col
        Atmp%val=A%val
        call sort_spMatCSR(Atmp) ! sort the col indices in A
        d = 0.0
        do i=1, Atmp%nRow ! loop through rows
        p=Atmp%row(i) ! mark the first none zero element in current row

        do j=Atmp%row(i), Atmp%row(i+1)-1 ! loop through columns
        if(Atmp%col(j).eq.i) then !diagonal
        d(i) = Atmp%val(j) ! store previous diagonal elements
        exit ! exit as we reached the last element in L
        elseif(Atmp%col(j).LT.i) then
        if(d(Atmp%col(j)).eq.0.0) then
        stop "Error: ilu0_Real > zero pivoting in ILU0 "
        endif 
        ! first divide each row in L with diagonal elements
        Atmp%val(j) = Atmp%val(j)/d(Atmp%col(j))
        piv = Atmp%val(j)
        do k = p+1, Atmp%row(i+1)-1 
        ! looping through every none-zero column in current row
        ! write( *, * ) Atmp%col(k), Atmp%val(k)
        do j2 = Atmp%row(Atmp%col(j)), Atmp%row(Atmp%col(j)+1)-1
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
        call lowerTri(Atmp, L)
        do i = 2, L%nRow+1 ! set the diagonal element of L to 1
        ! the diagonal element should be the last in each row 
        ! ONLY IF col is properly sorted
        L%val(L%row(i)-1) = 1.0
        enddo
        call upperTri(Atmp, U)
        call deall_spMatCSR(Atmp)
        return
    end subroutine ilu0_Real
    !
    !> Matix-Matrix sum, real version
    !
    subroutine RMATplusRMAT( A, B, C, tol )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: A, B
        type( spMatCSR_Real ), intent( out ) :: C
        real( kind=prec ), intent( in ), optional :: tol
        !
        type( spMatCSR_Real ) :: Ctemp
        integer :: i, j, k, m, n, nz, jj, i1, j1, j2, nnz, nzero
        logical               new
        real( kind=prec ) :: test, droptol, temp
        !
        !   test for size consistency
        if((A%nCol.NE.B%nCol).and.(A%nRow.NE.B%nRow)) then
        stop "Error: RMATplusRMAT > matrix sizes incompatible"
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
        call create_spMatCSR(m, n, nz, Ctemp)

        !   copy A into Ctemp, and add B, row by row.   
        !   Some entries of the sum may become zero
        j1 = 1
        nzero = 0
        do i=1, m
        Ctemp%row(i) = A%row(i+1)-A%row(i)
        j2 = j1+Ctemp%row(i)-1
        Ctemp%col(j1:j2) = A%col(A%row(i):A%row(i+1)-1)
        Ctemp%val(j1:j2) = A%val(A%row(i):A%row(i+1)-1)
        !    add elements of B in row i
        do k = B%row(i), B%row(i+1)-1
        new = .TRUE.
        do j = A%row(i), A%row(i+1)-1
        if(B%col(k).eq.A%col(j)) then
        !   add maatrix elements
        new =.FALSE.
        jj = j1+j-A%row(i)
        temp = B%val(k)+A%val(j)
        if(A%val(j).NE.0.0) then
        test = abs(temp)/abs(A%val(j))
        else
        test = abs(temp)/abs(B%val(k))
        endif
        if(test.GT. droptol) then
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
        !      number of elements(including possible zeros), not limits
        !    of col and val arrays.   So clean up Ctemp, put results in C
        nz = sum(Ctemp%row(1:m)) - nzero
        ! write(0, *) "nz = ", nz, " nzero = ", nzero
        call create_spMatCSR(m, n, nz, C)
        j1 = 1
        i1 = 0
        nz = 1
        do i = 1, m
        nnz = 0
        j2 = j1+Ctemp%row(i)-1  
        do j = j1, j2
        if(abs(Ctemp%val(j)).GT. 0) then
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
        !
        end subroutine RMATplusRMAT
        !
        ! sort the CSR col indices(and correspnding vals) in each row into
        ! ascent order
        ! reorder col indices like col:  4  1  3  2 -> 1  2  3  4
        !                          val:  3  1  1 -1 -> 1 -1  1  2 
        subroutine sort_spMatCSR_real(A)
        implicit none
        !
        type( spMatCSR_Real ), intent( inout ) :: A
        integer :: i, j1, j2, k
        integer, allocatable, dimension(:) :: idx
        real( kind=prec ), allocatable, dimension(:) :: col
        do i = 2, A%nRow+1
        j1=A%row(i-1)
        j2=A%row(i)-1
        allocate(col(j2-j1+1))
        allocate(idx(j2-j1+1))
        idx=(/(k, k=j1, j2, 1)/)
        col=real(A%col(j1:j2))
        call QSort(col, idx) 
        A%col(j1:j2)=int(col)
        A%val(j1:j2)=A%val(idx)
        deallocate(col)
        deallocate(idx)
        enddo
    end subroutine
    !
    !> Sort the CSR col indices(and correspnding vals) in each row into ascent order
    !> reorder col indices like col:  4  1  3  2 -> 1  2  3  4
    !>                          val:  3  1  1 -1 -> 1 -1  1  2
    !
    subroutine sort_spMatCSR_Cmplx(A)
        implicit none
        !
        type( spMatCSR_Cmplx ), intent( inout ) :: A
        integer :: i, j1, j2, k
        integer, allocatable, dimension(:) :: idx
        real( kind=prec ), allocatable, dimension(:) :: col
        do i = 2, A%nRow+1
        j1=A%row(i-1)
        j2=A%row(i)-1
        allocate(col(j2-j1+1))
        allocate(idx(j2-j1+1))
        idx=(/(k, k=j1, j2, 1)/)
        col=real(A%col(j1:j2))
        call QSort(col, idx) 
        !write( *, * ) "current row is: ", i-1
        !write( *, * ) A%col(j1:j2)
        A%col(j1:j2)=int(col)
        A%val(j1:j2)=A%val(idx)
        deallocate(col)
        deallocate(idx)
        enddo
    end subroutine
    !
    !> Silly subroutine to convert the real CSR sp matrix into complex
    !> WARNING: for now this is only used to form PETSc mat structure
    !>
    !> this assumes the matrix indices starts from ZERO instead of ONE
    !> need to be modified to use in other circumstances
    !
    subroutine RMAT2CMAT( R, C )
        implicit none
        !
        type( spMatCSR_Real ), intent( in ) :: R
        type( spMatCSR_Cmplx ), intent( out ) :: C
        !
        integer :: m, n, nnz
        !
        m = R%nRow
        n = R%nCol
        nnz = R%row(R%nRow+1)
        call create_spMatCSR_Cmplx(m, n, nnz, C)
        C%row = R%row
        C%col = R%col
        C%val = R%val
        !
    end subroutine
    !
end module SpOpTools
