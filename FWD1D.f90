program fwd1d

    use field1d, only: conf1d_t, fwd1d_timer, sourceField1d
    use field1d_s2, only: sourceField1d_s2
    use output_convention
    use modelspace
    use sg_vector
    use math_constants, only: PI, D2R, R2D
    implicit none

    ! ---------------------------------------------------------------------
    ! Solver selection. FWD1D can run EITHER of two independently-written
    ! 1D layered-sphere EM forward solvers, chosen here at compile time.
    ! Only ONE of the two character strings below is meaningful; typos are
    ! caught at runtime (see the errStop check right after the banner).
    !
    !   S1 (default) -- EARTH/FWD/field1d.f90, module `field1d`,
    !       subroutine sourceField1d. Written by Anna Kelbert in July 2011 to
    !       replicate the then-current MATLAB 1D code written by Jin Sun; used
    !       for global 3D inversion with complex sources. Originally written
    !       using the native conventions of Kelbert (2006); the best published
    !       reference for this code came much later: Sun, J., Kelbert, A., &
    !       Egbert, G. D. (2015), "Ionospheric current source modeling and
    !       global geomagnetic induction using ground geomagnetic observatory
    !       data", J. Geophys. Res. Solid Earth, 120(10), 6771-6796. In the
    !       2026-07 revisions that added S2 below, S1's native conventions
    !       were adjusted to match S2's for consistency between the two
    !       (see EARTH/FWD/output_convention.f90's native_convention()) --
    !       do NOT cite Kelbert et al. (2014) for this solver's formulation;
    !       that paper is a benchmark study, not the origin of its equations.
    !       Historically validated/compared against the companion MATLAB
    !       codebase (TSModel.m/tanField.m/radField.m in
    !       GlobalDV/ThinSheet/SIEM/@TSModel/) throughout this project.
    !
    !   S2 -- EARTH/FWD/field1d_s2.f90, module `field1d_s2`,
    !       subroutine sourceField1d_s2. Written by Anna Kelbert and Claude in
    !       July 2026, designed to closely follow Sun & Egbert (2012), "A
    !       thin-sheet model for global electromagnetic induction",
    !       Geophys. J. Int. 189, 343-356, Section 2 ("The homogeneous
    !       layered earth") + Appendix A, directly -- natively written in
    !       that paper's own conventions. Built explicitly as a
    !       cross-check: shares NO radial-solver or field-assembly code
    !       with S1 (only the already-validated angular/VSH
    !       routines, legendre_norm/vsharm, are reused -- see field1d.f90's
    !       widened public:: list). Independently validated to ~5
    !       significant figures against (a) a closed-form single-sphere
    !       solution (l=2,m=1) and (b) pythonSolver/spherical_em_
    !       induction.py (l=1,m=-1) -- see testing/TESTING_MANUAL.md.
    !
    ! CONVENTION / FUNCTIONALITY DIFFERENCES YOU MUST KNOW BEFORE COMPARING
    ! OUTPUT BETWEEN THE TWO SOLVERS:
    !
    !   - Time convention: BOTH are native e^{-i*omega*t} (identical
    !     kl=sqrt(i*omega*mu0*sigma) formula), so their raw output is
    !     directly comparable with no extra conjugation for time-convention
    !     reasons.
    !
    !   - Normalization: S1 normalizes against the incident
    !     reference potential at the surface (an internal "tni" quantity)
    !     and bakes extra R0^2/l(l+1) factors into its own H,E assembly
    !     formulas. S2 uses the paper's own literal convention:
    !     coeff(l,m)=1 means the coefficient of r^(l+1) in the toroidal
    !     potential T_l(r) is EXACTLY 1 in the air ("unit external
    !     multipole amplitude"), no extra R0^2/l(l+1) factors anywhere.
    !     The two solvers' raw output therefore differs by an overall
    !     complex constant per degree l -- NOT a bug, just a different
    !     (and, for S2, more directly physically interpretable)
    !     choice of normalization.
    !
    !   - CHANGED (2026-07-26): S1 used to apply a blanket conjg() to every
    !     assembled H/E component. This was NOT a bug -- it was a deliberate
    !     choice, from when S1 was originally written, to make its raw
    !     output match the e^{+i*omega*t} convention used by the global 3D
    !     inversion code S1 was written to feed into (the rest of S1's
    !     internal machinery, including the +m/-m reconstruction via
    !     vsharm's own m=0..lmax-only output, has always used
    !     e^{-i*omega*t} throughout). It has now been dropped: S1's native
    !     convention is meant to match S2's (and Sun & Egbert (2012)'s)
    !     e^{-i*omega*t} convention directly, with no implicit downstream
    !     conjugation baked in -- any e^{+i*omega*t} consumer (e.g. the
    !     global 3D code) is expected to apply that conjugation itself, the
    !     same way testing/test_vs_modem_1D_impedance.f90 conjugates ModEM's own
    !     e^{+i*omega*t}-labeled reference values rather than the solver's
    !     raw output. Separately (and unrelated to this convention change):
    !     S1's three T(r)-valued formulas (Hr, Etheta, Ephi) were fixed to
    !     carry the same explicit sign as S2/Sun & Egbert (2012) eq (5)-(6);
    !     and a genuine bug in the internal +m/-m paired-term reconstruction
    !     (its own conjg() treated the m=0 and m!=0 paths' shared Tnr/Tnsp
    !     radial factor inconsistently) was found and fixed separately.
    !     Both solvers now agree exactly (see testing/test_s1_vs_s2_l1m0.f90)
    !     and match an independent ModEM benchmark with no compensating
    !     sign anywhere (see testing/test_vs_modem_1D_impedance.f90).
    !
    !   - Coefficient scaling: NEITHER solver applies any extra scaling to
    !     the raw coeff(l,m) read from the source .prm file -- both are used
    !     AS-IS. (A Period/5 rescaling for S2, matching MATLAB's
    !     PrimaryField.m "shc" scaling before TSModel.ShcInc, was tried
    !     2026-07-24 and removed 2026-07-25 to keep the two solvers'
    !     driver-level behaviour identical and simple. If S2's absolute
    !     amplitude ever needs to be matched to the historical MATLAB/S1
    !     convention for a specific comparison, apply any such rescaling
    !     explicitly at the call site, not silently here.)
    !
    !   - Both solvers accept the identical conf1d_t earth structure and
    !     the identical flat coeff(:) array format (per-degree blocks
    !     ordered m=0,+1,-1,+2,-2,...,+l,-l -- see field1d_s2.f90's
    !     sourceField1d_s2 for the exact index formula).
    ! ---------------------------------------------------------------------
    character(len=20), parameter                :: SOLVER = 'S1'  ! 'S1' (default) or 'S2'

    ! Desired OUTPUT bookkeeping convention (orthogonal to SOLVER, which
    ! selects the numerics only) -- see EARTH/FWD/output_convention.f90 for
    ! the full definition of the field-bookkeeping dimensions this selects,
    ! and the named presets available here:
    !   'KELBERT2006', 'SUNEGBERT2012', 'EGBERTKELBERT2012'
    !   'EGBERTKELBERT2012_MODEM' (DEFAULT) -- same as EGBERTKELBERT2012 but
    !       ALSO rescales the output E and H by 1/(i*omega*mu0*G) so their
    !       amplitude AND phase match ModEM's boundary-E-field source
    !       normalization (rather than global1d's toroidal-potential
    !       normalization -- the two differ by the potential-vs-field factor
    !       c=i*omega*mu0*G of Faraday's law, G a real geometric constant set
    !       by the air-layer thickness). Preserves Z=E/H (no physics changed);
    !       matches ModEM to a few-% amplitude / <~5deg phase for the P10/MT
    !       modes. Full derivation: docs/source_normalization.md / .pdf.
    character(len=24), parameter                :: OUTPUT_CONVENTION = 'EGBERTKELBERT2012_MODEM'
    type(output_convention_t)                   :: target_conv, native

    ! Desired OUTPUT FIELD FORMAT (orthogonal to both SOLVER and
    ! OUTPUT_CONVENTION -- this controls only how the already-converted
    ! field data is serialized to disk, not any field-bookkeeping dimension):
    !   'MODEM'  (default) -- one binary file per run, byte-for-byte
    !       compatible with ModEM's own .esoln E-field-solution format (see
    !       ModEM's f90/3D_MT/ioMod/ioAscii.f90 FileWriteInit/EfileWrite),
    !       so it can be read directly with ModEM's own tools (e.g.
    !       testing/read_esoln.py) alongside a real ModEM run for direct
    !       comparison. Always writes E only (H is not part of ModEM's
    !       .esoln format), with nMode=1 (global1d computes one full
    !       3-component field per period, not ModEM's separate per-
    !       polarization solves -- ModeName is written as 'GLOBAL1D', not
    !       'Ex'/'Ey', to avoid implying a false correspondence to a single
    !       ModEM polarization). Requires primary_grid='E' (both solvers'
    !       native default) so e1d is EDGE-staggered, matching ModEM's E.
    !   'GLOBAL' -- the original per-period ASCII format (write_cvector),
    !       H and E both written, human-readable, native (phi,theta,r)
    !       component/axis ordering.
    character(len=10), parameter                :: OUTPUT_FORMAT = 'MODEM'

    type(conf1d_t)                              :: earth
    type(grid_t)                                :: grid
    type(modelParam_t)                          :: model,source,source_imag
    type(cvector)                               :: h1d, e1d
    character(len=1)                            :: primary_grid  ! set below from native%primary_grid ('H' or 'E')
    character(80)                               :: period_file,label
    character(80)                               :: layered_model_file
    character(80)                               :: source_model_file
    character(80)                               :: grid_file
    character(160)                              :: fields_output_file
    ! cfile = fields_output_file // solver_tag // OUTPUT_CONVENTION // primary_grid
    !         // '-grid.T' // period-index // '.hfield'/'.efield' -- character(80)
    !         silently truncated this for any non-trivial fields_output_file once
    !         the solver/convention tags were added to the name (2026-07-26).
    character(240)                              :: cfile
    character(20)                               :: solver_tag
    character(220)                              :: hdr
    character(16)                               :: period_str
    real(8), allocatable, dimension(:)          :: depths,coeff_real,coeff_imag,logrho,T
    real(8)                                     :: days
    character(3)                                :: ich, ichmode
    complex(8), allocatable, dimension(:)       :: coeff
    integer                                     :: i,icoeff,nL,nper,ncoeff,lmax,Nt,Np,Nr,narg,ios,istat
    ! n_modes/n_periods split nper (total periods*modes blocks, see below)
    ! into its two factors; iper/imode are this block's own period/mode
    ! index within the flat i=1..nper loop (periods-outer/modes-inner).
    integer                                     :: n_modes, n_periods, iper, imode
    character(30)                               :: argstr
    logical                                     :: apply_fake_center
    real(8)                                     :: fake_center_lat, fake_center_lon

    ! OUTPUT_FORMAT='MODEM' state: modem_nx/ny relabel global1d's own
    ! (phi-cells, theta-cells) as ModEM's own (latitude-cells, longitude-
    ! cells) -- see write_modem_header's header comment for why this is a
    ! relabeling, not a resize. ioMODEM is opened once before the period
    ! loop and closed once after (unlike the GLOBAL format's per-period
    ! open/close), since ModEM's own file format bundles every period into
    ! ONE file.
    integer                                     :: modem_nx, modem_ny, modem_nz
    character(240)                              :: modem_file
    character(len=20)                           :: modem_version
    real(8)                                     :: omega_i

    write(*,*) 'Copyright (c) 2010-2011 Oregon State University'
    write(*,*) 'College of Earth, Ocean and Atmospheric Sciences'
    write(*,*) 'Matlab code written by Jin Sun, last mod. 24 May 2010'
    write(*,*) 'Recoded in Fortran by Anna Kelbert, 11-13 July 2011'
    write(*,*) 'Data scaling updated by Anna Kelbert, 23-28 Nov 2011'
    write(*,*) 'Copyright (c) 2026 University of Colorado Boulder'
    write(*,*) 'Added E-field; Anna Kelbert with Claude, 9 June 2026'
    write(*,*) 'Primary E grid option; Anna Kelbert w Claude, 26 June 2026'
    write(*,*) 'Second (S2) solver + solver-selection option; Anna Kelbert w Claude, 24 July 2026'
    write(*,*) 'Optional "fake pole" grid recentering; Anna Kelbert w Claude, 25 July 2026'
    write(*,*) 'Grid, sign and harmonics convention preset options; Anna Kelbert w Claude, 27 July 2026'
    if (SOLVER == 'S1') then
        write(*,*) 'Solver: S1 (field1d.f90 -- Kelbert, 2011, orig. Jin Sun MATLAB code; see Sun, Kelbert & Egbert, 2015) [DEFAULT]'
    else if (SOLVER == 'S2') then
        write(*,*) 'Solver: S2 (field1d_s2.f90 -- Sun & Egbert, 2012, Section 2)'
    else
        call errStop('Unknown SOLVER flag; valid options are S1 (default) and S2')
    end if

    ! Resolve the solver's native output convention and the requested target
    ! convention (EARTH/FWD/output_convention.f90) -- primary_grid is set
    ! from the SOLVER's own native staggering (see output_convention.f90's
    ! native_convention()), not hardcoded, so switching SOLVER now also
    ! switches which grid staggering is used.
    native      = native_convention(SOLVER)
    target_conv = get_convention(OUTPUT_CONVENTION)
    primary_grid = native%primary_grid
    write(*,*) 'Output convention requested: ', trim(OUTPUT_CONVENTION)
    write(*,*) '  (solver native: time=',trim(native%time_convention),' norm=',trim(native%harmonic_norm), &
               ' CS=',native%condon_shortley,' theta=',trim(native%theta_convention), &
               ' r=',trim(native%r_convention),' primary=',native%primary_grid,')'
    write(*,*) '  (target:        time=',trim(target_conv%time_convention),' norm=',trim(target_conv%harmonic_norm), &
               ' CS=',target_conv%condon_shortley,' theta=',trim(target_conv%theta_convention), &
               ' r=',trim(target_conv%r_convention),' primary=',target_conv%primary_grid,')'
    write(*,*)

    !  parse command line
    narg = command_argument_count()
    if (narg < 4) then
        write(0,*) 'Usage: ./FWD1D layered_model_file source_model_file grid_file fields_output_file &
                    &[fake_center_lat_deg fake_center_lon_deg]'
        write(0,*) '  The two optional trailing arguments, if given, must be given together: they'
        write(0,*) '  request a "fake pole" recentering of the input grid (see recenter_grid_fake_pole'
        write(0,*) '  below) -- useful for evaluating a regional grid (e.g. USA.0.25x0.25.grd) as if'
        write(0,*) '  it were located somewhere else on the sphere (e.g. at the equator/zero meridian,'
        write(0,*) '  or wherever the field of interest peaks), WITHOUT editing the grid file itself.'
        stop
    end if

    !call get_command_argument(1, period_file)
    call get_command_argument(1, layered_model_file)
    call get_command_argument(2, source_model_file)
    call get_command_argument(3, grid_file)
    call get_command_argument(4, fields_output_file)

    ! optional 5th/6th arguments: "fake pole" grid recentering (both or neither)
    apply_fake_center = .false.
    if (narg == 5) then
        call errStop('FWD1D: fake_center_lat_deg and fake_center_lon_deg must both be given together &
                      &(got only one of the two optional trailing arguments)')
    else if (narg >= 6) then
        apply_fake_center = .true.
        call get_command_argument(5, argstr)
        read(argstr,*) fake_center_lat
        call get_command_argument(6, argstr)
        read(argstr,*) fake_center_lon
    end if

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

    ! get the periods from the source file. nper = TOTAL blocks (periods x
    ! modes, flat, periods-outer/modes-inner -- see ModelSpace.f90's
    ! read_modelParam); n_periods/n_modes are the two separate counts used
    ! below to tag output files and drive the MODEM-format nMode dimension.
    ! n_modes defaults to 1 (single-mode files), reproducing the original
    ! one-block-per-period behavior exactly.
    nper = source%nL
    n_modes = source%nModes
    n_periods = nper / n_modes
    allocate(T(nper), STAT=istat)
    do i = 1,nper
      days = source%L(i)%period
      T(i) = days * (24*3600)
    end do
    lmax = getDegree_modelParam(source)
    ncoeff = (lmax + 1)**2 ! number of SH

    ! reading grid file (r is in km decreasing from top to bottom)
    call read_grid(grid,grid_file)

    ! optionally recenter the grid onto a "fake pole" -- see subroutine header
    ! comment below for exactly what this does and does not do.
    if (apply_fake_center) then
        call recenter_grid_fake_pole(grid, fake_center_lat, fake_center_lon)
    end if

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
    if (SOLVER == 'S1') then
        solver_tag = 's1'
    else
        solver_tag = 's2'
    end if

    if (OUTPUT_FORMAT == 'MODEM' .and. primary_grid /= 'E') then
        call errStop('OUTPUT_FORMAT=MODEM requires primary_grid=E (ModEM''s own E-field solution &
                      &is always EDGE-staggered) -- both solvers'' native convention is E already, &
                      &so this should not happen unless native_convention() has been changed')
    end if

    if (OUTPUT_FORMAT == 'MODEM') then
        modem_file = trim(fields_output_file)//'.'//trim(solver_tag)//'.'//trim(OUTPUT_CONVENTION)// &
                     '.esoln'
        write(*,*) 'Writing ModEM-format E-field solution to: ',trim(modem_file)
        call write_modem_header(modem_file, grid, n_periods, n_modes)
    end if

    icoeff = 0

    do i = 1,nper
        ! block i (flat, periods-outer/modes-inner -- see ModelSpace.f90's
        ! read_modelParam) belongs to period-group iper, mode imode within
        ! it; both are used below to tag output filenames/records so a
        ! multi-mode run's files never collide (n_modes=1 -> imode=1
        ! always, reproducing the original single-mode filenames' period
        ! tag exactly, just with an added, harmless '.mode01').
        iper  = (i-1)/n_modes + 1
        imode = mod(i-1,n_modes) + 1
        write(ich,'(i2.2)') iper
        write(ichmode,'(i2.2)') imode

        days = T(i)/(24*3600)
        write(*,*) 'Computing the fields for period ',trim(ich),' mode ',trim(ichmode),': ',days,' days'

        ! (a) rescale source coeffs into the target harmonic basis BEFORE
        ! solving -- see output_convention.f90's rescale_source_coeffs for
        ! exactly what this does (norm/Condon-Shortley only; no-op when
        ! target_conv already matches vsharm's own fixed FULL/CS-included
        ! basis, e.g. OUTPUT_CONVENTION='SUNEGBERT2012').
        call rescale_source_coeffs(coeff(icoeff+1:icoeff+ncoeff), lmax, target_conv)

        if (SOLVER == 'S1') then
            call sourceField1d(earth,lmax,coeff(icoeff+1:icoeff+ncoeff),T(i),grid,h1d,e1d)
        else
            call sourceField1d_s2(earth,lmax,coeff(icoeff+1:icoeff+ncoeff),T(i),grid,h1d,e1d)
        end if
        icoeff = icoeff + ncoeff

        ! (b) map native solver output -> target convention AFTER solving --
        ! time/theta/r only; no-op in every dimension where native already
        ! matches target_conv (see output_convention.f90's
        ! apply_output_convention).
        call apply_output_convention(h1d, e1d, native, target_conv)

        ! (c) optional ModEM source-normalization rescale (only for target
        ! conventions with modem_normalize=.true., e.g. EGBERTKELBERT2012_MODEM)
        ! -- rescales E and H by 1/(i*omega*mu0*G), G=(3/2)sqrt(3/4pi)*d_air,
        ! from the air-layer thickness (grid top radius minus Earth radius), so
        ! the output matches ModEM's boundary-E-field normalization instead of
        ! global1d's toroidal-potential normalization. Preserves Z. See
        ! output_convention.f90's apply_modem_normalization / the
        ! EGBERTKELBERT2012_MODEM preset / docs/source_normalization.pdf.
        if (target_conv%modem_normalize) then
            call apply_modem_normalization(h1d, e1d, 2.0d0*PI/T(i), grid%r(1)*1.0d3 - earth%r0)
        end if

        call reset_time(fwd1d_timer)

        ! NOTE: once theta and/or r differ between native and target_conv,
        ! apply_output_convention has REVERSED the theta and/or r index of
        ! every component array above -- the (i,j,k) indices in the files
        ! written below now correspond to DIFFERENT physical grid points
        ! than grid_file's own th()/r() arrays would suggest under the
        ! solver's native convention. Downstream tooling must be told which
        ! convention was used (see the header line and filename tag below).
        if (OUTPUT_FORMAT == 'MODEM') then
            omega_i = 2.0d0*PI/T(i)
            call write_modem_efield_block(e1d, iper, imode, omega_i)
            write(*,*) 'Wrote period ',trim(ich),' mode ',trim(ichmode),' to ',trim(modem_file)
        else if (OUTPUT_FORMAT == 'GLOBAL') then
            cfile = trim(fields_output_file)//'.'//trim(solver_tag)//'.'//trim(OUTPUT_CONVENTION)// &
                    '.'//primary_grid//'-grid.T'//trim(ich)//'.mode'//trim(ichmode)//'.hfield'
            write(*,*) 'Writing to H file: ',cfile
            open(ioWRITE,file=cfile,status='unknown',form='formatted',iostat=ios)
            write(period_str,'(f9.3)') T(i)
            hdr = "# FWD1D full EM field solution H output for period "//trim(adjustl(period_str))//" secs, mode "// &
                  trim(ichmode)//". OUTPUT_CONVENTION="//trim(OUTPUT_CONVENTION)//" SOLVER="//trim(SOLVER)// &
                  " time="//trim(target_conv%time_convention)//" norm="//trim(target_conv%harmonic_norm)// &
                  " theta="//trim(target_conv%theta_convention)//" r="//trim(target_conv%r_convention)
            write(ioWRITE,'(a)') trim(hdr)
            write(ioWRITE,'(i3)') 1
            call write_cvector(ioWRITE,h1d)
            close(ioWRITE)

            cfile = trim(fields_output_file)//'.'//trim(solver_tag)//'.'//trim(OUTPUT_CONVENTION)// &
                    '.'//primary_grid//'-grid.T'//trim(ich)//'.mode'//trim(ichmode)//'.efield'
            write(*,*) 'Writing to E file: ',cfile
            open(ioWRITE,file=cfile,status='unknown',form='formatted',iostat=ios)
            write(period_str,'(f9.3)') T(i)
            hdr = "# FWD1D full EM field solution E output for period "//trim(adjustl(period_str))//" secs, mode "// &
                  trim(ichmode)//". OUTPUT_CONVENTION="//trim(OUTPUT_CONVENTION)//" SOLVER="//trim(SOLVER)// &
                  " time="//trim(target_conv%time_convention)//" norm="//trim(target_conv%harmonic_norm)// &
                  " theta="//trim(target_conv%theta_convention)//" r="//trim(target_conv%r_convention)
            write(ioWRITE,'(a)') trim(hdr)
            write(ioWRITE,'(i3)') 1
            call write_cvector(ioWRITE,e1d)
            close(ioWRITE)
        else
            call errStop('Unknown OUTPUT_FORMAT flag; valid options are MODEM (default) and GLOBAL')
        end if

        write(*,*) 'Done writing to file: ',elapsed_time(fwd1d_timer),' secs'
    end do

    if (OUTPUT_FORMAT == 'MODEM') then
        close(ioE)
    end if

    deallocate(depths,coeff,coeff_real,coeff_imag,logrho,T, STAT=istat)
    deallocate(earth%layer,earth%sigma, STAT=istat)
    call deall_modelParam(model)
    call deall_modelParam(source)
    call deall_modelParam(source_imag)
    call deall_cvector(h1d)
    call deall_cvector(e1d)
    call deall_grid(grid)
    write(*,*) 'Total time taken: ',saved_time(fwd1d_timer),' secs'

contains

    subroutine recenter_grid_fake_pole(grid_arg, center_lat_deg, center_lon_deg)
        ! Recenters an existing (theta,phi) grid, IN PLACE, onto a "fake"
        ! center point, by a simple additive shift of every theta node value
        ! and phi node value -- NOT a true 3-D spherical rotation. This is
        ! the standard "fake pole" trick for evaluating a small-to-moderate
        ! regional patch (e.g. USA.0.25x0.25.grd) as if it were located
        ! somewhere else on the sphere: the grid's own internal (theta,phi)
        ! structure (cell widths grid%dt/grid%dp, extent, resolution) is
        ! preserved EXACTLY; only its absolute position on the sphere
        ! changes, by translating every node by the same (dtheta,dphi) --
        ! i.e. treating the patch as if the sphere were locally flat over
        ! its extent. This is an approximation, not an exact rotation (an
        ! exact rotation would also change the *shape* of a rectangular
        ! theta x phi patch as it moves away from the equator/moves in
        ! longitude); it is adequate for patches that are not enormous
        ! relative to the sphere, which is the regime this option is meant
        ! for. Only field1d.f90/field1d_s2.f90's OWN sourceField1d
        ! calls (later in this program) are affected -- the source
        ! coefficients (coeff(:), from source_model_file) are NOT touched;
        ! they remain defined relative to the TRUE global pole/meridian, so
        ! the physical multipole source pattern does not move, only the
        ! observation grid does.
        !
        ! center_lat_deg/center_lon_deg (input, degrees, ordinary geographic
        ! latitude [-90,90] and longitude [-180,180] or [0,360)): the point
        ! that the grid's own geometric center -- the midpoint of its
        ! min/max theta extent and of its min/max phi extent, in its
        ! ORIGINAL, as-read-from-file coordinates -- is moved to. E.g.
        ! center_lat_deg=0, center_lon_deg=-90 places the grid's center at
        ! the equator, 90 degrees West.
        type(grid_t), intent(inout) :: grid_arg
        real(8), intent(in)         :: center_lat_deg, center_lon_deg
        real(8) :: theta_c_old, phi_c_old, theta_c_new, phi_c_new, dtheta, dphi
        integer :: ii

        theta_c_old = (minval(grid_arg%th(1:grid_arg%ny+1)) + maxval(grid_arg%th(1:grid_arg%ny+1))) / 2.0d0
        phi_c_old   = (minval(grid_arg%ph(1:grid_arg%nx+1)) + maxval(grid_arg%ph(1:grid_arg%nx+1))) / 2.0d0

        theta_c_new = (90.0d0 - center_lat_deg) * D2R
        phi_c_new   = modulo(center_lon_deg, 360.0d0) * D2R

        dtheta = theta_c_new - theta_c_old
        dphi   = phi_c_new - phi_c_old

        do ii = 1, grid_arg%ny+1
            grid_arg%th(ii) = grid_arg%th(ii) + dtheta
        end do
        do ii = 1, grid_arg%nx+1
            grid_arg%ph(ii) = modulo(grid_arg%ph(ii) + dphi, 2.0d0*PI)
        end do
        ! grid_arg%dt/grid_arg%dp (cell widths) are unaffected by a uniform node shift.

        write(*,*)
        write(*,*) 'Applying fake-pole grid recentering (see recenter_grid_fake_pole in FWD1D.f90):'
        write(*,'(a,f9.4,a,f9.4,a)') '   original grid center (colat,lon) = (', theta_c_old*R2D, ', ', phi_c_old*R2D, ') deg'
        write(*,'(a,f9.4,a,f9.4,a)') '   requested center (lat,lon)       = (', center_lat_deg, ', ', center_lon_deg, ') deg'
        write(*,'(a,f9.4,a,f9.4,a)') '   shift applied (dtheta,dphi)      = (', dtheta*R2D, ', ', dphi*R2D, ') deg'
        write(*,*)

        if (minval(grid_arg%th(1:grid_arg%ny+1)) <= 0.0d0 .or. maxval(grid_arg%th(1:grid_arg%ny+1)) >= PI) then
            call errStop('recenter_grid_fake_pole: requested center pushes the grid across a pole &
                          &(theta out of the open interval (0,180) deg) -- choose a fake center &
                          &farther from the poles, or a smaller/narrower grid extent')
        end if

    end subroutine recenter_grid_fake_pole

    subroutine write_modem_header(fname, grid_arg, nper_arg, nmode_arg)
        ! Opens ioE (unformatted, sequential) and writes ModEM's own 4-record
        ! .esoln header (see ModEM's f90/3D_MT/ioMod/ioAscii.f90's
        ! FileWriteInit): version, nPer, nMode, nx, ny, nz, nzAir, ox, oy,
        ! oz, rotdeg -- then dx(1:nx), dy(1:ny), dz(1:nz). nper_arg/
        ! nmode_arg come directly from the source .prm file's own "periods
        ! P modes M" header (see ModelSpace.f90's read_modelParam) -- when
        ! the source file doesn't specify modes, nmode_arg=1, exactly
        ! matching ModEM's own single-mode-per-period .esoln files.
        !
        ! AXIS RELABELING (not a resize): ModEM's own grid convention is
        ! x=latitude(north-south)-related, y=longitude(east-west)-related,
        ! z=depth (confirmed via ModEM's f90/3D_MT/FWD/boundary_wsS.f90's
        ! cos(latitude) longitude-spacing correction, only applied to its
        ! "y"/imode=2 slice). global1d's own grid convention is x=phi
        ! (longitude), y=theta (colatitude), z=r -- the OPPOSITE axis
        ! assignment. So ModEM's "nx" (latitude cells) = global1d's
        ! grid%ny (theta cells), and ModEM's "ny" (longitude cells) =
        ! global1d's grid%nx (phi cells) -- a relabeling of which physical
        ! axis is "1" vs "2", not a different grid resolution.
        character(*), intent(in)       :: fname
        type(grid_t), intent(in)       :: grid_arg
        integer, intent(in)            :: nper_arg, nmode_arg
        real(8), allocatable           :: modem_dx(:), modem_dy(:), modem_dz(:)
        integer                        :: ios2

        modem_nx = grid_arg%ny   ! latitude cells  = global1d's theta cells
        modem_ny = grid_arg%nx   ! longitude cells = global1d's phi cells
        modem_nz = grid_arg%nz   ! depth cells (unchanged, same axis)

        allocate(modem_dx(modem_nx), modem_dy(modem_ny), modem_dz(modem_nz))
        ! dx: latitude cell widths, degrees. Per an independent trace of
        ! ModEM's own GridCalcS.f90/boundary_wsS.f90/modelParam_IO_WS.f90,
        ! ModEM's "x"(latitude) index runs SOUTH(1)->NORTH(nx+1), the
        ! OPPOSITE of global1d's native th (COLAT_N2S, north(1)->south).
        ! This is exactly the same reversal apply_output_convention already
        ! applied to the field DATA under EGBERTKELBERT2012's
        ! theta_convention=LAT_S2N ("latitude, South->North") -- so dx must
        ! be reversed relative to grid_arg%dt here too, to stay consistent
        ! with the (already-reversed) field arrays below.
        modem_dx(1:modem_nx) = grid_arg%dt(modem_nx:1:-1) * R2D
        ! dy: longitude cell widths, degrees -- no reversal (ModEM's "y"
        ! index runs WEST->EAST, same direction as global1d's native ph;
        ! there is no "phi_convention" dimension in output_convention.f90,
        ! so this axis is never touched by apply_output_convention either).
        modem_dy = grid_arg%dp(1:modem_ny) * R2D
        ! dz: depth cell widths, METERS -- from global1d's own dr (radius
        ! cell widths, km), same top-to-bottom order. Confirmed both by an
        ! independent trace of ModEM's own dz convention (top-to-bottom,
        ! thick air layers first) AND empirically on the global1d side:
        ! global1d's own field magnitude is ~0 at k=1 (top of the air
        ! column, matching ModEM's negligible-air-conductivity top BC) and
        ! grows toward the surface -- i.e. despite r_convention differing
        ! (R_UP native vs R_DOWN target) and apply_output_convention
        ! reversing the r-index, k=1 remains "top" on both sides here, so
        ! no additional index reversal is needed for dz.
        modem_dz = grid_arg%dr(1:modem_nz) * 1.0d3

        modem_version = 'GLOBAL1D'
        open(ioE, file=fname, status='unknown', form='unformatted', iostat=ios2)
        if (ios2 /= 0) call errStop('write_modem_header: could not open '//trim(fname))

        ! ox,oy,oz,rotdeg: per the same independent trace of ModEM's own
        ! grid_t/setup_grid/modelParam_IO_WS.f90 -- ox = latitude (deg) of
        ! the SOUTH edge (index 1, matching the dx reversal above), oy =
        ! longitude (deg) of the WEST edge (index 1, no reversal), oz = 0.0
        ! (ModEM's own default -- an Earth-surface-referenced offset, not
        ! used in any coordinate transform for a nominal EARTH_R model),
        ! rotdeg = 0.0 (confirmed dead/vestigial in ModEM -- read and
        ! compared for file-consistency checking only, never applied to
        ! any coordinate). NOTE: this derivation could not be checked
        ! against a real ModEM-written spherical .esoln file (none exists
        ! in that repo to compare against) -- cross-check against one if
        ! available. This does NOT affect the E-field DATA itself (shape/
        ! order/sign, verified separately against known ASCII-format
        ! reference values) -- only this file's own descriptive geographic
        ! header fields.
        write(ioE) modem_version, int(nper_arg,4), int(nmode_arg,4), int(modem_nx,4), int(modem_ny,4), &
                   int(modem_nz,4), int(grid_arg%nzAir,4), &
                   90.0d0 - grid_arg%th(grid_arg%ny+1)*R2D, grid_arg%ph(1)*R2D, 0.0d0, 0.0d0
        write(ioE) modem_dx
        write(ioE) modem_dy
        write(ioE) modem_dz

        deallocate(modem_dx, modem_dy, modem_dz)
    end subroutine write_modem_header

    subroutine write_modem_efield_block(e1d_arg, ifreq, imode_arg, omega)
        ! Writes one ModEM-format period/mode block: Omega/iFreq/iMode/
        ! ModeName header record + Ex/Ey/Ez arrays (see ModEM's EfileWrite).
        ! iFreq/iMode come from the caller's own period/mode-group indices
        ! (see the periods-outer/modes-inner block ordering established in
        ! ModelSpace.f90's read_modelParam); ModeName is written as
        ! 'GLOBAL1D#' (# = imode_arg) rather than 'Ex'/'Ey' to avoid
        ! implying a false correspondence to a specific ModEM polarization
        ! -- global1d's coeff(:) source for a given mode can be an
        ! arbitrary spherical-harmonic pattern, not necessarily a pure
        ! Ex/Ey-equivalent one.
        !
        ! COMPONENT MAPPING: e1d_arg's own %x/%y/%z are global1d's native
        ! phi/theta/r slots. Under EGBERTKELBERT2012 (or any target
        ! convention whose theta_convention/r_convention have already been
        ! applied via apply_output_convention BEFORE this is called), %y
        ! directly equals the MT "north" component and %z the MT "down"
        ! component (the whole point of that transform -- see
        ! output_convention.f90's NATIVE CONVENTIONS note -- is that it
        ! replaces the old manual Ex=-Etheta/Hz=-Hr relabeling, so no
        ! further sign flip is applied here); %x (phi/east) is unaffected
        ! either way. So: ModEM's Ex (north) <- e1d_arg%y, ModEM's Ey
        ! (east) <- e1d_arg%x, ModEM's Ez (down) <- e1d_arg%z -- each
        ! transposed on its first two axes, since global1d stores every
        ! component (phi,theta,r) and ModEM stores every component
        ! (latitude,longitude,depth) = (theta,phi,r), the SAME relabeling
        ! as write_modem_header's nx/ny swap, applied here to the data
        ! itself component-by-component (each of Ex/Ey/Ez keeps its own
        ! true EDGE extent -- see sg_vector.f90's create_cvector -- no
        ! zero-padding, unlike the GLOBAL ASCII format).
        !
        ! r-AXIS REVERSAL (bug found 2026-07-27, fixed same day): the
        ! theta/r reversal apply_output_convention already applied to
        ! e1d_arg (native R_UP -> target R_DOWN) reverses the r-index of
        ! the FIELD DATA -- e1d_arg's own k=1 is native's k=nz(+1), i.e.
        ! DEEP (near the conductive core, correctly near-zero there), and
        ! k=nz(+1) is native's k=1, i.e. the TOP of the 1000km air column
        ! (correctly the largest value there, for this "external multipole,
        ! amplitude grows with r" source convention -- see CLAUDE.md). But
        ! write_modem_header's own dz array is built from grid_arg%dr
        ! DIRECTLY, UNREVERSED (confirmed correct: it matches a real
        ! ModEM run's own dz to 6 sig figs, i.e. dz index 1 = the THICK
        ! top-of-air layer, same as ModEM). Writing e1d_arg's r-axis
        ! as-is alongside that dz would silently pair "index 1" in the
        ! geometry (top of air) with "index 1" in the field data (bottom
        ! of the core) -- exactly backwards. Reverse the field data's
        ! r-axis here (on top of the existing 1,2-transpose) so both
        ! agree: file index 1 = top of air (for BOTH dz and the field
        ! values), file index nz = bottom/core.
        type(cvector), intent(in)  :: e1d_arg
        integer, intent(in)        :: ifreq, imode_arg
        real(8), intent(in)        :: omega
        character(len=20)          :: mode_name

        write(mode_name,'(a,i0)') 'GLOBAL1D', imode_arg
        write(ioE) omega, int(ifreq,4), int(imode_arg,4), mode_name
        write(ioE) transpose_12_reverse_3(e1d_arg%y)   ! ModEM's Ex (north)
        write(ioE) transpose_12_reverse_3(e1d_arg%x)   ! ModEM's Ey (east)
        write(ioE) transpose_12_reverse_3(e1d_arg%z)   ! ModEM's Ez (down)
    end subroutine write_modem_efield_block

    function transpose_12_reverse_3(A) result(B)
        ! Swaps the first two axes of a rank-3 complex array AND reverses
        ! the third: B(j,i,n3+1-k) = A(i,j,k). The 1,2-swap converts
        ! global1d's (phi,theta,r)-ordered component arrays into ModEM's
        ! (latitude,longitude,depth) = (theta,phi,r)-ordered arrays; the
        ! 3rd-axis reversal undoes apply_output_convention's own r-index
        ! reversal (see write_modem_efield_block's header comment) so the
        ! result's r-axis matches write_modem_header's UNREVERSED dz
        ! (both: index 1 = top of air, index n3 = bottom/core). Explicit
        ! loop, not RESHAPE(...,ORDER=...), to avoid any ambiguity about
        ! RESHAPE's element-order semantics for a correctness-critical,
        ! one-shot conversion.
        complex(8), intent(in)  :: A(:,:,:)
        complex(8), allocatable :: B(:,:,:)
        integer :: n1, n2, n3, i1, i2, i3
        n1 = size(A,1); n2 = size(A,2); n3 = size(A,3)
        allocate(B(n2,n1,n3))
        do i3 = 1,n3
            do i2 = 1,n2
                do i1 = 1,n1
                    B(i2,i1,n3+1-i3) = A(i1,i2,i3)
                end do
            end do
        end do
    end function transpose_12_reverse_3

end program fwd1d
