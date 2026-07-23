program fwd1d

    use field1d, only: conf1d_t, fwd1d_timer, sourceField1d
    use field1d_sunegbert2012, only: sourceField1d_sunegbert2012
    use modelspace
    use sg_vector
    implicit none

    ! ---------------------------------------------------------------------
    ! Solver selection. FWD1D can run EITHER of two independently-written
    ! 1D layered-sphere EM forward solvers, chosen here at compile time.
    ! Only ONE of the two character strings below is meaningful; typos are
    ! caught at runtime (see the errStop check right after the banner).
    !
    !   KELBERT2014 (default) -- EARTH/FWD/field1d.f90, module `field1d`,
    !       subroutine sourceField1d. The ORIGINAL solver, based on
    !       Kelbert et al. (2014), "Global 3-D electromagnetic forward
    !       modelling: a benchmark study", Geophys. J. Int. Historically
    !       validated/compared against the companion MATLAB codebase
    !       (TSModel.m/tanField.m/radField.m in
    !       GlobalDV/ThinSheet/SIEM/@TSModel/) throughout this project.
    !
    !   SUNEGBERT2012 -- EARTH/FWD/field1d_sunegbert2012.f90, module `field1d_sunegbert2012`,
    !       subroutine sourceField1d_sunegbert2012. A SECOND, INDEPENDENTLY-DERIVED
    !       solver, written directly from Sun & Egbert (2012), "A
    !       thin-sheet model for global electromagnetic induction",
    !       Geophys. J. Int. 189, 343-356, Section 2 ("The homogeneous
    !       layered earth") + Appendix A. Written 2026-07 explicitly as a
    !       cross-check: shares NO radial-solver or field-assembly code
    !       with KELBERT2014 (only the already-validated angular/VSH
    !       routines, legendre_norm/vsharm, are reused -- see field1d.f90's
    !       widened public:: list). Independently validated to ~5
    !       significant figures against (a) a closed-form single-sphere
    !       solution (l=2,m=1) and (b) pythonSolver/spherical_em_
    !       induction.py (l=1,m=-1) -- see testing/TESTING_MANUAL.md.
    !
    ! CONVENTION / FUNCTIONALITY DIFFERENCES YOU MUST KNOW BEFORE COMPARING
    ! OUTPUT BETWEEN THE TWO SOLVERS (full detail/derivations in CLAUDE.md):
    !
    !   - Time convention: BOTH are native e^{-i*omega*t} (identical
    !     kl=sqrt(i*omega*mu0*sigma) formula), so their raw output is
    !     directly comparable with no extra conjugation for time-convention
    !     reasons.
    !
    !   - Normalization: KELBERT2014 normalizes against the incident
    !     reference potential at the surface (an internal "tni" quantity)
    !     and bakes extra R0^2/l(l+1) factors into its own H,E assembly
    !     formulas. SUNEGBERT2012 uses the paper's own literal convention:
    !     coeff(l,m)=1 means the coefficient of r^(l+1) in the toroidal
    !     potential T_l(r) is EXACTLY 1 in the air ("unit external
    !     multipole amplitude"), no extra R0^2/l(l+1) factors anywhere.
    !     The two solvers' raw output therefore differs by an overall
    !     complex constant per degree l -- NOT a bug, just a different
    !     (and, for SUNEGBERT2012, more directly physically interpretable)
    !     choice of normalization.
    !
    !   - KELBERT2014 applies a final conjg() to every assembled H/E
    !     component. This was introduced to compensate for a "+m,-m
    !     conjugate-pairing" trick used to reconstruct negative-m terms
    !     from vsharm's own m=0..lmax-only output. As of 2026-07-24 this
    !     conjg() is suspected -- precisely diagnosed, not yet confirmed
    !     fixed -- to also be applied to UNPAIRED terms (e.g. a pure m=0
    !     source), where no pairing-compensation is actually needed: since
    !     conj(z)/z = exp(-2i*arg(z)) and Tnr(T)/Tnsp(T') generically have
    !     DIFFERENT complex phases from each other, the blanket conjg()
    !     introduces a spurious ~90-degree relative phase between
    !     T(r)-based output (Hr,Ephi) and T'(r)-based output (Htheta,Hphi)
    !     for such terms that is NOT present in SUNEGBERT2012 (whose
    !     negative-m terms are built directly via the identity
    !     Y_l^{-m}=(-1)^m*conjg(Y_l^m), term by term -- no blanket
    !     conjugation of the assembled sum is ever needed). See CLAUDE.md,
    !     "field1d.f90 vs field1d_sunegbert2012.f90 direct comparison (l=1,m=0)",
    !     for the full derivation and the ongoing investigation into
    !     whether/how this connects to the historical P10/Mode2 E-field
    !     sign question.
    !
    !   - Coefficient scaling: NEITHER solver applies any extra scaling to
    !     the raw coeff(l,m) read from the source .prm file -- both are used
    !     AS-IS. (A Period/5 rescaling for SUNEGBERT2012, matching MATLAB's
    !     PrimaryField.m "shc" scaling before TSModel.ShcInc, was tried
    !     2026-07-24 and removed 2026-07-25 to keep the two solvers'
    !     driver-level behaviour identical and simple; see CLAUDE.md,
    !     "Matching MATLAB/SIEM's practical coefficient scaling" for the
    !     historical note. If SUNEGBERT2012's absolute amplitude ever needs
    !     to be matched to the historical MATLAB/KELBERT2014 convention for
    !     a specific comparison, apply any such rescaling explicitly at the
    !     call site, not silently here.)
    !
    !   - Both solvers accept the identical conf1d_t earth structure and
    !     the identical flat coeff(:) array format (per-degree blocks
    !     ordered m=0,+1,-1,+2,-2,...,+l,-l -- see field1d_sunegbert2012.f90's
    !     sourceField1d_sunegbert2012 for the exact index formula).
    ! ---------------------------------------------------------------------
    character(len=20), parameter                :: SOLVER = 'KELBERT2014'  ! 'KELBERT2014' (default) or 'SUNEGBERT2012'

    type(conf1d_t)                              :: earth
    type(grid_t)                                :: grid
    type(modelParam_t)                          :: model,source,source_imag
    type(cvector)                               :: h1d, e1d
    character(len=1), parameter                 :: primary_grid = 'E'  ! 'H' or 'E'
    character(80)                               :: period_file,label
    character(80)                               :: layered_model_file
    character(80)                               :: source_model_file
    character(80)                               :: grid_file
    character(80)                               :: fields_output_file
    character(80)                               :: cfile
    character(20)                               :: solver_tag
    real(8), allocatable, dimension(:)          :: depths,coeff_real,coeff_imag,logrho,T
    real(8)                                     :: days
    character(3)                                :: ich
    complex(8), allocatable, dimension(:)       :: coeff
    integer                                     :: i,icoeff,nL,nper,ncoeff,lmax,Nt,Np,Nr,narg,ios,istat

    write(*,*) 'Copyright (c) 2010-2011 Oregon State University'
    write(*,*) 'College of Earth, Ocean and Atmospheric Sciences'
    write(*,*) 'Matlab code written by Jin Sun, last mod. 24 May 2010'
    write(*,*) 'Recoded in Fortran by Anna Kelbert, 11-13 July 2011'
    write(*,*) 'Data scaling updated by Anna Kelbert, 23-28 Nov 2011'
    write(*,*) 'Added E-field; Anna Kelbert with Claude, 9 June 2026'
    write(*,*) 'Primary E grid option; Anna Kelbert w Claude, 26 June 2026'
    write(*,*) 'Second (SUNEGBERT2012) solver + solver-selection option; Anna Kelbert w Claude, 24 July 2026'
    if (SOLVER == 'KELBERT2014') then
        write(*,*) 'Solver: KELBERT2014 (field1d.f90 -- Kelbert et al., 2014, benchmark study) [DEFAULT]'
    else if (SOLVER == 'SUNEGBERT2012') then
        write(*,*) 'Solver: SUNEGBERT2012 (field1d_sunegbert2012.f90 -- Sun & Egbert, 2012, Section 2)'
    else
        call errStop('Unknown SOLVER flag; valid options are KELBERT2014 (default) and SUNEGBERT2012')
    end if
    write(*,*)

    !  parse command line
    narg = command_argument_count()
    if (narg < 4) then
        write(0,*) 'Usage: ./FWD1D layered_model_file source_model_file grid_file fields_output_file'
        stop
    end if

    !call get_command_argument(1, period_file)
    call get_command_argument(1, layered_model_file)
    call get_command_argument(2, source_model_file)
    call get_command_argument(3, grid_file)
    call get_command_argument(4, fields_output_file)

    ! save periods in days
!    open(ioREAD,file=period_file,status='old',form='formatted',iostat=ios)
!    write(6,*) 'Reading from the periods file ',trim(period_file)
!    read(ioREAD,'(a)') label
!    write(6,*) label
!    read(ioREAD,*) nper
!    allocate(T(nper), STAT=istat)
!    do i = 1,nper
!      read(ioREAD,*) days ! reading period in *days*
!      T(i) = days * (24*3600)
!    end do
!    close(ioREAD)

    ! model file should be 1D layered
    call read_modelParam(model,layered_model_file)
    if (model%nL /= model%nc) then
        write(0,*) 'Error in FWD1D: input model file is not layered 1D'
        stop
    end if
    nL = model%nL
    allocate(depths(nL),logrho(nL), STAT=istat)
    do i = 1,nL
        depths(i) = model%L(i)%depth
    end do
    call getParamValues_modelParam(model,logrho)

    ! source file contains the complex source for multiple periods
    call read_modelParam(source,source_model_file,source_imag)
    allocate(coeff_real(source%nc),coeff_imag(source%nc),coeff(source%nc), STAT=istat)
    call getParamValues_modelParam(source,coeff_real)
    call getParamValues_modelParam(source_imag,coeff_imag)
    coeff = dcmplx(coeff_real,coeff_imag)
    !write(*,*) 'coeff: ', coeff

    ! get the periods from the source file
    nper = source%nL
    allocate(T(nper), STAT=istat)
    do i = 1,nper
      days = source%L(i)%period
      T(i) = days * (24*3600)
    end do
    lmax = getDegree_modelParam(source)
    ncoeff = (lmax + 1)**2 ! number of SH

    ! reading grid file (r is in km decreasing from top to bottom)
    call read_grid(grid,grid_file)

    ! set earth radius and domain top radius (in meters)
    earth%r0  = 6371.0e3
    earth%rmax= 1.0e3 * grid%r(1)

    ! set tolerance on toroidal potential
    earth%tol = 1.e-9

    ! set surface conductance (should be small since we're using 3D thinsheet)
    earth%tau = 1.e-4

    ! save model in 1D configuration structure: layers include the core
    allocate(earth%layer(nL+1),earth%sigma(nL+1), STAT=istat)
    earth%layer(1) = 0.0d0
    earth%layer(2:nL+1) = 1.0e3 * depths(1:nL)
    earth%sigma(1:nL) = 10.0**( -logrho(1:nL) )
    earth%sigma(nL+1) = 10.0**( 5.0 ) ! core conductivity

    write(*,*) 'Tops of model layers: ',earth%layer
    write(*,*) 'Conductivity values:  ',earth%sigma

    ! allocate the output cvectors
    if (primary_grid == 'H') then
       ! magnetic field is primary in the global model; this is the default
       write(*,*) 'Grid staggering: H on primary edges (EDGE), E on primary faces (FACE)'
       call create_cvector(grid, h1d, EDGE)
       call create_cvector(grid, e1d, FACE)
    else if (primary_grid == 'E') then
       ! electric field is primary; used e.g. as input to ModEM regional.
       ! Note: coordinate conventions remain global — convert before use in MT.
       write(*,*) 'Grid staggering: E on primary edges (EDGE), H on primary faces (FACE)'
       call create_cvector(grid, e1d, EDGE)
       call create_cvector(grid, h1d, FACE)
    else
       call errStop('Unknown primary_grid flag; valid options are H (default) and E')
    end if

    ! output filename tag, so the two solvers' output never collides
    if (SOLVER == 'KELBERT2014') then
        solver_tag = 'kelbert2014'
    else
        solver_tag = 'sunegbert2012'
    end if

    icoeff = 0

    do i = 1,nper
        write(ich,'(i2.2)') i

        days = T(i)/(24*3600)
        write(*,*) 'Computing the fields for period ',trim(ich),': ',days,' days'

        if (SOLVER == 'KELBERT2014') then
            ! coeff(l,m) used exactly as read from the source .prm file, no extra scaling.
            call sourceField1d(earth,lmax,coeff(icoeff+1:icoeff+ncoeff),T(i),grid,h1d,e1d)
        else
            ! coeff(l,m) used exactly as read from the source .prm file, no extra scaling
            ! (a Period/5 rescaling was tried here 2026-07-24 and removed 2026-07-25 --
            ! see the SOLVER comment block above and CLAUDE.md).
            call sourceField1d_sunegbert2012(earth,lmax,coeff(icoeff+1:icoeff+ncoeff),T(i),grid,h1d,e1d)
        end if
        icoeff = icoeff + ncoeff

        call reset_time(fwd1d_timer)

        cfile = trim(fields_output_file)//'.'//trim(solver_tag)//'.'//primary_grid//'-grid.T'//trim(ich)//'.hfield'
        write(*,*) 'Writing to H file: ',cfile
        open(ioWRITE,file=cfile,status='unknown',form='formatted',iostat=ios)
        write(ioWRITE,'(a47,f9.3,a6)') "# FWD1D full EM field solution H output for period ",   &
                                            T(i),' secs.'
        write(ioWRITE,'(i3)') 1
        call write_cvector(ioWRITE,h1d)
        close(ioWRITE)

        cfile = trim(fields_output_file)//'.'//trim(solver_tag)//'.'//primary_grid//'-grid.T'//trim(ich)//'.efield'
        write(*,*) 'Writing to E file: ',cfile
        open(ioWRITE,file=cfile,status='unknown',form='formatted',iostat=ios)
        write(ioWRITE,'(a47,f9.3,a6)') "# FWD1D full EM field solution E output for period ",   &
                                            T(i),' secs.'
        write(ioWRITE,'(i3)') 1
        call write_cvector(ioWRITE,e1d)
        close(ioWRITE)

        write(*,*) 'Done writing to file: ',elapsed_time(fwd1d_timer),' secs'
    end do

    deallocate(depths,coeff,coeff_real,coeff_imag,logrho,T, STAT=istat)
    deallocate(earth%layer,earth%sigma, STAT=istat)
    call deall_modelParam(model)
    call deall_modelParam(source)
    call deall_modelParam(source_imag)
    call deall_cvector(h1d)
    call deall_cvector(e1d)
    call deall_grid(grid)
    write(*,*) 'Total time taken: ',saved_time(fwd1d_timer),' secs'

end program fwd1d
