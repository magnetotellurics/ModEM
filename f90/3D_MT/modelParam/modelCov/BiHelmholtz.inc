!----------------------------------------------------------------------!
! 3D_MT H2 (Bi-Helmholtz) covariance (CovType=2): procedures.          !
!                                                                      !
!   H2 precision operator:  A = kappa^2 I - Delta                      !
!   Matérn order alpha=2:     C_m = A^{-2},  K = A^{-1}               !
!----------------------------------------------------------------------!

! *******************************************************************
  subroutine read_CmSqrt_BiHelm(cfile)

    ! Reads the H2 (Bi-Helmholtz) covariance configuration from file.
    ! Supports the unified cov file format:
    !   Line 17: Nx Ny NzEarth FileType  (or legacy Nx Ny NzEarth)
    !   Sx, Sy, Sz (placeholders)
    !   Parameter array (CovType-specific)
    !   nrules, exception rules (3 values per rule)
    !   mask blocks
    ! When fileType=1 (AR format) is detected with CovType=2 selected,
    ! reads the mask and uses H2 defaults.

    character(*), intent(in)                         :: cfile
    integer                                          :: fid, j, fileType
    integer                                          :: Nx, Ny, NzEarth, nrules
    integer                                          :: istat
    real (kind=prec), allocatable                    :: tempSx(:), tempSy(:)
    real (kind=prec)                                 :: tempSz
    integer                                          :: tempN

    fid = 30
    open(unit=fid,file=cfile,form='formatted',status='old')

    ! skip 16-line header
    do j = 1, 16
       read(fid,*)
    end do

    ! detect filetype and read grid dimensions from line 17
    call detect_cov_filetype(fid, Nx, Ny, NzEarth, fileType)
    if (Nx <= 0 .or. Ny <= 0 .or. NzEarth <= 0) then
       call errStop('Invalid grid dimensions in covariance file '//cfile)
    end if
    CmSqrtBiHelm%Nx = Nx
    CmSqrtBiHelm%Ny = Ny
    CmSqrtBiHelm%NzEarth = NzEarth

    if (fileType == 2) then
        ! normal H2 read: skip Sx, Sy, Sz (placeholders)
       allocate(tempSx(NzEarth), tempSy(NzEarth), STAT=istat)
       read(fid,*) tempSx
       read(fid,*) tempSy
       read(fid,*) tempSz
       deallocate(tempSx, tempSy, STAT=istat)
       ! read parameter array: N kappa pcg_tol pcg_maxIt
       ! All 4 values on one line; read in a single statement so that
       ! list-directed I/O does not advance past the record boundary.
       read(fid,*,iostat=istat) tempN, CmSqrtBiHelm%kappa, &
            CmSqrtBiHelm%pcg_tol, CmSqrtBiHelm%pcg_maxIt
       if (istat /= 0) then
           call errStop('Error reading H2 (Bi-Helmholtz) parameters in '//cfile)
       end if
       if (CmSqrtBiHelm%kappa <= 0.0_prec) then
           call errStop('H2 kappa must be positive in covariance file '//cfile)
       end if
       if (CmSqrtBiHelm%pcg_tol <= 0.0_prec) then
           call errStop('H2 pcg_tol must be positive in covariance file '//cfile)
       end if
       if (CmSqrtBiHelm%pcg_maxIt <= 0) then
           call errStop('H2 pcg_maxIt must be positive in covariance file '//cfile)
       end if
    else
        ! fileType=1 with CovType=2 selected: read mask, use H2 defaults
        write(0,*) 'Warning: AR cov file read with CovType=2.'
        write(0,*) '  Using H2 defaults for kappa, pcg parameters.'
       allocate(tempSx(NzEarth), tempSy(NzEarth), STAT=istat)
       read(fid,*) tempSx
       read(fid,*) tempSy
       read(fid,*) tempSz
       deallocate(tempSx, tempSy, STAT=istat)
       ! skip N (AR smoothing passes)
       read(fid,*) tempN
       CmSqrtBiHelm%kappa = 0.5_prec
       CmSqrtBiHelm%pcg_tol = 1.0e-7_prec
       CmSqrtBiHelm%pcg_maxIt = 250
    end if

    ! read exception rules (3 values per rule; third value unused for CovType=2)
    read(fid,*,iostat=istat) nrules
    if (istat /= 0 .or. nrules < 0) then
       call errStop('Invalid nrules in covariance file '//cfile)
    end if
    CmSqrtBiHelm%nrules = nrules
    CmSqrtBiHelm%tearMatrix = .false.
    if (nrules > 0) then
       allocate(CmSqrtBiHelm%tearPairs(2, nrules), stat=istat)
       if (istat /= 0) call errStop('Memory allocation failed for tearPairs')
       do j = 1, nrules
          read(fid,*) CmSqrtBiHelm%tearPairs(1,j), CmSqrtBiHelm%tearPairs(2,j), tempSz
          if (CmSqrtBiHelm%tearPairs(1,j) < 0 .or. CmSqrtBiHelm%tearPairs(1,j) > 9 .or. &
              CmSqrtBiHelm%tearPairs(2,j) < 0 .or. CmSqrtBiHelm%tearPairs(2,j) > 9) then
             call errStop('Invalid mask code in tear rule (must be 0-9)')
          end if
          CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%tearPairs(1,j), CmSqrtBiHelm%tearPairs(2,j)) = .true.
          CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%tearPairs(2,j), CmSqrtBiHelm%tearPairs(1,j)) = .true.
       end do
    end if

    ! create and read the mask array
    call read_iscalar(fid, CmSqrtBiHelm%mask)

    close(fid)

  end subroutine read_CmSqrt_BiHelm

! *******************************************************************
  subroutine create_CmSqrt_BiHelm(m, cfile)

    ! Initializes the H2 (Bi-Helmholtz) covariance. If cfile is provided, reads
    ! configuration from file; otherwise uses defaults.

    type (modelParam_t), intent(in)     :: m
    character(*), intent(in), optional  :: cfile
    logical                             :: exists

    ! initialize dimensions
    CmSqrtBiHelm%Nx = m%Nx
    CmSqrtBiHelm%Ny = m%Ny
    CmSqrtBiHelm%NzEarth = m%NzEarth

    ! defaults
    CmSqrtBiHelm%kappa = 0.5_prec
    CmSqrtBiHelm%pcg_tol = 1.0e-7_prec
    CmSqrtBiHelm%pcg_maxIt = 250

    ! no exception rules by default
    CmSqrtBiHelm%tearMatrix = .false.
    CmSqrtBiHelm%nrules = 0
    nullify(CmSqrtBiHelm%tearPairs)

    call create_iscalar(m%grid, CmSqrtBiHelm%mask, CELL_EARTH)
    CmSqrtBiHelm%mask%v = FREE

    if (present(cfile)) then
       inquire(FILE=cfile, EXIST=exists)
       if (exists) then
          call read_CmSqrt_BiHelm(cfile)
       else
          call errStop('Unable to find covariance file '//trim(cfile)//' in create_CmSqrt_BiHelm')
       end if
       if ((CmSqrtBiHelm%Nx /= m%Nx) .or. (CmSqrtBiHelm%Ny /= m%Ny) .or. &
           (CmSqrtBiHelm%NzEarth /= m%NzEarth)) then
          call errStop('Grid dimensions do not match in covariance file '//cfile)
       end if
    end if

    ! assemble the sparse operator A = kappa^2 I - Delta
    call assemble_A_BiHelm(m)

    if (output_level > 3) then
       write(6,*) 'create_CmSqrt_BiHelm: CovType=2, Nx=', CmSqrtBiHelm%Nx, &
            ' Ny=', CmSqrtBiHelm%Ny, ' NzEarth=', CmSqrtBiHelm%NzEarth, &
            ' nCells=', CmSqrtBiHelm%Nx * CmSqrtBiHelm%Ny * CmSqrtBiHelm%NzEarth
       write(6,*) '  kappa (cell units)=', CmSqrtBiHelm%kappa, &
            ' (correlation length ~', 1.0_prec/CmSqrtBiHelm%kappa, ' cells)'
       write(6,*) '  pcg_tol=', CmSqrtBiHelm%pcg_tol, &
            ' pcg_maxIt=', CmSqrtBiHelm%pcg_maxIt
       write(6,*) '  operator nnz=', size(CmSqrtBiHelm%A%col)
    end if

    CmSqrtBiHelm%allocated = .true.

  end subroutine create_CmSqrt_BiHelm

! *******************************************************************
  subroutine deall_CmSqrt_BiHelm()

    integer  :: istat

    call deall_spMatCSR(CmSqrtBiHelm%A)
    call deall_iscalar(CmSqrtBiHelm%mask)
    if (associated(CmSqrtBiHelm%tearPairs)) deallocate(CmSqrtBiHelm%tearPairs)
    CmSqrtBiHelm%nrules = 0
    CmSqrtBiHelm%tearMatrix = .false.
    CmSqrtBiHelm%allocated = .false.

  end subroutine deall_CmSqrt_BiHelm

! *******************************************************************
  subroutine assemble_A_BiHelm(m)

    ! Assembles the sparse CSR operator A = kappa^2 I - Delta
    ! over earth cells using a graph-Laplacian stencil in cell-index
    ! space (unit weights).  This is consistent with the AR(1) smoother
    ! which also operates in cell-index coordinates — correlation between
    ! adjacent cells does not depend on physical grid spacing.
    !
    !   A = kappa^2 * I - L_graph
    !
    ! where L_graph has unit off-diagonal entries (all neighbors at
    ! distance 1 in index space).  With kappa > 0 this guarantees a
    ! well-conditioned PD operator whose eigenvalues are O(1) regardless
    ! of the model's physical extent.
    !
    ! Exception rules (tearMatrix) optionally disable edges between
    ! specific mask-code pairs, creating "tears" in the smoothing.
    !
    ! Counting pass first, then fill pass.

    type (modelParam_t), intent(in) :: m
    integer                         :: Nx, Ny, NzEarth, nCells
    integer                         :: i, j, k, idx, idx_nb
    integer                         :: nnz, nz, istat, p
    real (kind=prec)               :: kappa2
    real (kind=prec)               :: diagA
    integer, allocatable            :: rowT(:)

    Nx = m%Nx
    Ny = m%Ny
    NzEarth = m%NzEarth
    nCells = Nx * Ny * NzEarth
    kappa2 = CmSqrtBiHelm%kappa**2

    if (.not. m%grid%allocated) then
       call errStop('Grid must be allocated before assemble_A_BiHelm')
    end if

    if (CmSqrtBiHelm%A%allocated) then
       call deall_spMatCSR(CmSqrtBiHelm%A)
    end if

    ! Pass 1: count non-zeros per row
    allocate(rowT(nCells+1), stat=istat)
    rowT = 0
    nz = 0
    do k = 1, NzEarth
       do j = 1, Ny
          do i = 1, Nx
             idx = (k-1)*Nx*Ny + (j-1)*Nx + i
             nnz = 1
             if (i > 1) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i-1,j,k))) nnz = nnz + 1
             end if
             if (i < Nx) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i+1,j,k))) nnz = nnz + 1
             end if
             if (j > 1) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j-1,k))) nnz = nnz + 1
             end if
             if (j < Ny) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j+1,k))) nnz = nnz + 1
             end if
             if (k > 1) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j,k-1))) nnz = nnz + 1
             end if
             if (k < NzEarth) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j,k+1))) nnz = nnz + 1
             end if
             rowT(idx) = nnz
             nz = nz + nnz
          end do
       end do
    end do

    call create_spMatCSR(nCells, nCells, nz, CmSqrtBiHelm%A)

    CmSqrtBiHelm%A%row(1) = 1
    do idx = 2, nCells+1
       CmSqrtBiHelm%A%row(idx) = CmSqrtBiHelm%A%row(idx-1) + rowT(idx-1)
    end do
    deallocate(rowT, stat=istat)

    ! Pass 2: fill columns and values
    do k = 1, NzEarth
       do j = 1, Ny
          do i = 1, Nx
             idx = (k-1)*Nx*Ny + (j-1)*Nx + i
             p = CmSqrtBiHelm%A%row(idx)

             diagA = kappa2
             if (i > 1) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i-1,j,k))) diagA = diagA + 1.0_prec
             end if
             if (i < Nx) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i+1,j,k))) diagA = diagA + 1.0_prec
             end if
             if (j > 1) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j-1,k))) diagA = diagA + 1.0_prec
             end if
             if (j < Ny) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j+1,k))) diagA = diagA + 1.0_prec
             end if
             if (k > 1) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j,k-1))) diagA = diagA + 1.0_prec
             end if
             if (k < NzEarth) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j,k+1))) diagA = diagA + 1.0_prec
             end if

             ! diagonal
             CmSqrtBiHelm%A%col(p) = idx
             CmSqrtBiHelm%A%val(p) = diagA
             p = p + 1

             ! X-left
             if (i > 1) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i-1,j,k))) then
                   idx_nb = (k-1)*Nx*Ny + (j-1)*Nx + (i-1)
                   CmSqrtBiHelm%A%col(p) = idx_nb
                   CmSqrtBiHelm%A%val(p) = -1.0_prec
                   p = p + 1
                end if
             end if
             ! X-right
             if (i < Nx) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i+1,j,k))) then
                   idx_nb = (k-1)*Nx*Ny + (j-1)*Nx + (i+1)
                   CmSqrtBiHelm%A%col(p) = idx_nb
                   CmSqrtBiHelm%A%val(p) = -1.0_prec
                   p = p + 1
                end if
             end if
             ! Y-left
             if (j > 1) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j-1,k))) then
                   idx_nb = (k-1)*Nx*Ny + (j-2)*Nx + i
                   CmSqrtBiHelm%A%col(p) = idx_nb
                   CmSqrtBiHelm%A%val(p) = -1.0_prec
                   p = p + 1
                end if
             end if
             ! Y-right
             if (j < Ny) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j+1,k))) then
                   idx_nb = (k-1)*Nx*Ny + j*Nx + i
                   CmSqrtBiHelm%A%col(p) = idx_nb
                   CmSqrtBiHelm%A%val(p) = -1.0_prec
                   p = p + 1
                end if
             end if
             ! Z-left
             if (k > 1) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j,k-1))) then
                   idx_nb = (k-2)*Nx*Ny + (j-1)*Nx + i
                   CmSqrtBiHelm%A%col(p) = idx_nb
                   CmSqrtBiHelm%A%val(p) = -1.0_prec
                   p = p + 1
                end if
             end if
             ! Z-right
             if (k < NzEarth) then
                if (.not. CmSqrtBiHelm%tearMatrix(CmSqrtBiHelm%mask%v(i,j,k), &
                    CmSqrtBiHelm%mask%v(i,j,k+1))) then
                   idx_nb = k*Nx*Ny + (j-1)*Nx + i
                   CmSqrtBiHelm%A%col(p) = idx_nb
                   CmSqrtBiHelm%A%val(p) = -1.0_prec
                   p = p + 1
                end if
             end if

          end do
       end do
    end do

  end subroutine assemble_A_BiHelm

! *******************************************************************
  subroutine apply_A_BiHelm(v_in, v_out)

    ! Sparse matvec: v_out = A * v_in over earth cells.

    real (kind=prec), intent(in)  :: v_in(:,:,:)
    real (kind=prec), intent(out) :: v_out(:,:,:)
    integer                       :: Nx, Ny, NzEarth, nCells, i, j, k, idx
    real (kind=prec), allocatable :: x(:), y(:)

    Nx = size(v_in, 1)
    Ny = size(v_in, 2)
    NzEarth = size(v_in, 3)
    nCells = Nx * Ny * NzEarth

    allocate(x(nCells), y(nCells))

    ! flatten 3D -> 1D
    do k = 1, NzEarth
       do j = 1, Ny
          do i = 1, Nx
             idx = (k-1)*Nx*Ny + (j-1)*Nx + i
             x(idx) = v_in(i,j,k)
          end do
       end do
    end do

    ! CSR matvec via spOpTools
    call RMATxRVEC(CmSqrtBiHelm%A, x, y)

    ! unflatten 1D -> 3D
    do k = 1, NzEarth
       do j = 1, Ny
          do i = 1, Nx
             idx = (k-1)*Nx*Ny + (j-1)*Nx + i
             v_out(i,j,k) = y(idx)
          end do
       end do
    end do

    deallocate(x, y)

  end subroutine apply_A_BiHelm

! *******************************************************************
  subroutine apply_freeze_semantics_BiHelm(v)

    ! Zeros out air (mask==0) and ocean (mask==9) cells.

    real (kind=prec), intent(inout) :: v(:,:,:)
    integer                         :: Nx, Ny, NzEarth, i, j, k

    Nx = size(v, 1)
    Ny = size(v, 2)
    NzEarth = size(v, 3)

    do k = 1, NzEarth
       do j = 1, Ny
          do i = 1, Nx
             if (CmSqrtBiHelm%mask%v(i,j,k) == AIR .or. &
                 CmSqrtBiHelm%mask%v(i,j,k) == OCEAN) then
                v(i,j,k) = 0.0_prec
             end if
          end do
       end do
    end do

  end subroutine apply_freeze_semantics_BiHelm

! *******************************************************************
  subroutine solve_A_BiHelm(rhs, sol, tol, maxIt, nIt, ok)

    ! Solve A * sol = rhs using PCG (Jacobi preconditioner).

    real (kind=prec), intent(in)  :: rhs(:,:,:)
    real (kind=prec), intent(out) :: sol(:,:,:)
    real (kind=prec), intent(in)  :: tol
    integer, intent(in)           :: maxIt
    integer, intent(out)          :: nIt
    logical, intent(out)          :: ok

    integer                       :: Nx, Ny, NzEarth, nCells
    integer                       :: i, j, k, idx, it, istat
    real (kind=prec), allocatable :: r(:), p(:), z(:), q(:), x(:), b(:)
    real (kind=prec), allocatable :: diag(:)
    real (kind=prec)              :: alpha, beta, rsold, rsnew, bnrm, denom
    logical                       :: earlyExit

    Nx = size(rhs, 1)
    Ny = size(rhs, 2)
    NzEarth = size(rhs, 3)
    nCells = Nx * Ny * NzEarth

    allocate(b(nCells), x(nCells), r(nCells), p(nCells), z(nCells), q(nCells), stat=istat)
    if (istat /= 0) call errStop('Memory allocation failed in solve_A_BiHelm')
    allocate(diag(nCells), stat=istat)
    if (istat /= 0) call errStop('Memory allocation failed in solve_A_BiHelm (diag)')

    ! flatten rhs
    do k = 1, NzEarth
       do j = 1, Ny
          do i = 1, Nx
             idx = (k-1)*Nx*Ny + (j-1)*Nx + i
             b(idx) = rhs(i,j,k)
          end do
       end do
    end do

    ! extract diagonal for Jacobi preconditioner
    call diag_Real(CmSqrtBiHelm%A, diag)

    ! initial guess = 0
    x = 0.0_prec

    ! r = b - A * 0 = b
    r = b

    ! solve M * z = r  (Jacobi: z = r ./ diag)
    do idx = 1, nCells
       if (abs(diag(idx)) > R_TINY) then
          z(idx) = r(idx) / diag(idx)
       else
          z(idx) = 0.0_prec
       end if
    end do

    p = z
   rsold = dot_product(r, z)
    bnrm = sqrt(dot_product(b, b))

    if (output_level > 3) then
       write(6,*) 'solve_A_BiHelm: starting PCG, tol=', tol, ' maxIt=', maxIt, ' bnrm=', bnrm
    end if

    if (bnrm < R_TINY) then
       sol = 0.0_prec
       ok = .true.
       nIt = 0
       deallocate(b, x, r, p, z, q, diag)
       return
    end if

    if (rsold /= rsold) then
       sol = 0.0_prec
       ok = .false.
       nIt = 0
       deallocate(b, x, r, p, z, q, diag)
       return
    end if

    ok = .false.
    earlyExit = .false.
    do it = 1, maxIt
       ! q = A * p
       call RMATxRVEC(CmSqrtBiHelm%A, p, q)

        denom = dot_product(p, q)
        if ((abs(denom) < R_TINY) .or. (denom /= denom)) then
           ok = .false.
           earlyExit = .true.
           nIt = it
           exit
        end if

        alpha = rsold / denom
        if (alpha /= alpha) then
           ok = .false.
           earlyExit = .true.
           nIt = it
           exit
        end if

       x = x + alpha * p
       r = r - alpha * q

       ! check convergence
        if (sqrt(dot_product(r, r)) / bnrm < tol) then
           if (output_level > 4) then
              write(6,*) 'solve_A_BiHelm: it=', it, &
                   ' rel_res=', sqrt(dot_product(r, r)) / bnrm, ' (converged)'
           end if
           ok = .true.
           nIt = it
           exit
        end if

        if (output_level > 4) then
           write(6,*) 'solve_A_BiHelm: it=', it, &
                ' rel_res=', sqrt(dot_product(r, r)) / bnrm
        end if

        ! solve M * z = r (Jacobi)
       do idx = 1, nCells
          if (abs(diag(idx)) > R_TINY) then
             z(idx) = r(idx) / diag(idx)
          else
             z(idx) = 0.0_prec
          end if
       end do

        rsnew = dot_product(r, z)
        if (rsnew /= rsnew) then
           ok = .false.
           earlyExit = .true.
           nIt = it
           exit
        end if
        if ((abs(rsold) < R_TINY) .or. (rsold /= rsold)) then
           ok = .false.
           earlyExit = .true.
           nIt = it
           exit
        end if

        beta = rsnew / rsold
        if (beta /= beta) then
           ok = .false.
           earlyExit = .true.
           nIt = it
           exit
        end if

       p = z + beta * p
       rsold = rsnew
    end do

    if (.not. ok .and. .not. earlyExit) then
       nIt = maxIt
    end if

    if (output_level > 3) then
       if (ok) then
          write(6,*) 'solve_A_BiHelm: converged in ', nIt, &
               ' iterations, rel_res=', sqrt(dot_product(r, r)) / bnrm
       else
          write(6,*) 'solve_A_BiHelm: FAILED at it=', nIt, &
               ' rel_res=', sqrt(dot_product(r, r)) / bnrm
       end if
    end if

    ! unflatten solution
    do k = 1, NzEarth
       do j = 1, Ny
          do i = 1, Nx
             idx = (k-1)*Nx*Ny + (j-1)*Nx + i
             sol(i,j,k) = x(idx)
          end do
       end do
    end do

    deallocate(b, x, r, p, z, q, diag)

  end subroutine solve_A_BiHelm

! *******************************************************************
  function multBy_CmSqrt_BiHelm(mhat) result(dm)

    ! K = A^{-1}: solve A * dm = mhat

    type (modelParam_t), intent(in) :: mhat
    type (modelParam_t)             :: dm
    integer                         :: nIt
    logical                         :: ok

    if (.not. CmSqrtBiHelm%allocated) then
       call create_CmSqrt_BiHelm(mhat)
    end if

    dm = mhat

    call solve_A_BiHelm(mhat%cellCond%v, dm%cellCond%v, &
         CmSqrtBiHelm%pcg_tol, CmSqrtBiHelm%pcg_maxIt, nIt, ok)

    if (.not. ok) then
       call errStop('PCG failed to converge in multBy_CmSqrt_BiHelm')
    end if

    call apply_freeze_semantics_BiHelm(dm%cellCond%v)

    dm%temporary = .true.

  end function multBy_CmSqrt_BiHelm

! *******************************************************************
  function multBy_Cm_BiHelm(mhat) result(dm)

    ! C_m = A^{-2}: two solves with A

    type (modelParam_t), intent(in) :: mhat
    type (modelParam_t)             :: dm, tmp
    integer                         :: nIt
    logical                         :: ok

    if (.not. CmSqrtBiHelm%allocated) then
       call create_CmSqrt_BiHelm(mhat)
    end if

    tmp = mhat
    call solve_A_BiHelm(mhat%cellCond%v, tmp%cellCond%v, &
         CmSqrtBiHelm%pcg_tol, CmSqrtBiHelm%pcg_maxIt, nIt, ok)

    if (.not. ok) then
       call errStop('PCG failed to converge in multBy_Cm_BiHelm (first solve)')
    end if

    dm = tmp
    call solve_A_BiHelm(tmp%cellCond%v, dm%cellCond%v, &
         CmSqrtBiHelm%pcg_tol, CmSqrtBiHelm%pcg_maxIt, nIt, ok)

    if (.not. ok) then
       call errStop('PCG failed to converge in multBy_Cm_BiHelm (second solve)')
    end if

    call apply_freeze_semantics_BiHelm(dm%cellCond%v)

    dm%temporary = .true.
    call deall_modelParam(tmp)

  end function multBy_Cm_BiHelm

! *******************************************************************
  function multBy_CmSqrtInv_BiHelm(dm) result(mhat)

    ! K^{-1} = A: one sparse matvec

    type (modelParam_t), intent(in) :: dm
    type (modelParam_t)             :: mhat

    if (.not. CmSqrtBiHelm%allocated) then
       call create_CmSqrt_BiHelm(dm)
    end if

    mhat = dm
    call apply_A_BiHelm(dm%cellCond%v, mhat%cellCond%v)

    call apply_freeze_semantics_BiHelm(mhat%cellCond%v)

    mhat%temporary = .true.

  end function multBy_CmSqrtInv_BiHelm
