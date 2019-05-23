! *****************************************************************************
program FWDtest
!  program for testing "wrapped" WS 2D forward modeling code

!     use the interface or "wrapping" modules
        use wsfwd2dpar, only:  SetWSparams
	use fwdtmmod
	use fwdtemod
	implicit none
	
	real (kind=prec), allocatable, dimension(:,:) :: Cond2D
	integer (kind=4) :: Nzb, IER, nPer, i, nHXB
	real (kind=prec), allocatable, dimension(:) :: periods
        !  Solution, boundary conditions for TM test
	complex (kind=prec), allocatable, dimension(:,:) :: Hsol
	complex (kind=prec), allocatable, dimension(:)	:: HXB
        !  Solution, boundary conditions for TE test
	complex (kind=prec), allocatable, dimension(:,:) :: Esol
	complex (kind=prec), allocatable, dimension(:)	:: EXB

        type(grid2d_t) 	:: TMgrid

	character*80 	cfile
	
	! open file for grid/conductivity model
	cfile = 'CondTest'
	open(1,file = cfile, form = 'unformatted')
	! open file for TM test solution
	cfile = 'Hsol.out'
	open(2,file = cfile, form = 'unformatted')
	! open file for TE test solution
	cfile = 'Esol.out'
	open(3,file = cfile, form = 'unformatted')

	!  Read in grid geometry definitions, store in structure TMgrid
	!    first grid dimensions ...
        !  Same grid2d_t structure is used for TM/TE
	read(1) TMgrid%Nz, TMgrid%Nza, TMgrid%Ny, nPer
        Nzb = TMgrid%Nz-TMgrid%Nza
	nHXB = 2*(TMgrid%Ny+1)+2*(TMgrid%Nz+1)

	!    then allocate arrays
        allocate(TMgrid%Dz(TMgrid%Nz))
        allocate(TMgrid%Dy(TMgrid%Ny))
        allocate(Cond2D(TMgrid%Ny,TMgrid%Nz))
        allocate(Hsol(TMgrid%Ny+1,Nzb+1))
        allocate(Esol(TMgrid%Ny+1,TMgrid%Nz+1))
        allocate(periods(nPer))
	allocate(HXB(nHXB))
	allocate(EXB(nHXB))

	!    read in grid spacings
        read(1) TMgrid%Dy
        read(1) TMgrid%Dz
	!    .... and conductivity (not part of TMgrid)
        read(1) Cond2D
        ! Finally, read in list of periods
        read(1) periods
        !  and close input file
        close(1)

        !  set parameters formerly specified through include file
        call setWSparams(TMgrid%Ny,TMgrid%Nz)
 
        !  Testing ... just compute for one period for now 
        !   these are the initial grid setup calls.  Only need to do
        !    once for each mode.  For changes in conductivity, need to call
        !    UpdateCondTM(Cond2D) (or UpdateCondTE) before doing next steps:1
        call FWD2DSetupTM(TMgrid,Cond2D,IER)
        if(IER.lt.0) then
           write(0,*) 'Error Initializing for TM mode soln:IER=',IER
           stop
        endif
        call FWD2DSetupTE(TMgrid,IER)
        if(IER.lt.0) then
           write(0,*) 'Error Initializing for TE mode soln:IER=',IER
           stop
        endif

        ! TM test 
        ! complete setting up TM equations for this period
        !   this routine must be called after each change of
        !     frequency
        call UpdateFreqTM(periods(1))
        ! compute boundary conditions for 2D TM by solving 1D
        !  equations at boundary ... put results in HXB
        call SetBoundTM(periods(1),HXB)
        !  solve 2D TM equations, return solution in Hsol
        call Fwd2DsolveTM(HXB,Hsol,IER)

        ! write solution to file
        write(2) TMgrid%Ny+1,Nzb+1
        write(2) TMgrid%Dy
        write(2) (TMgrid%Dz(i),i=TMgrid%Nza+1,TMgrid%Nz)
        write(2) Hsol

        ! deallocate
        call Fwd2DdeallTM()
        close(2)

        ! TE test 
        ! complete setting up TE equations for this period
        !   this routine must be called after each change of
        !     frequency
        call UpdateFreqTE(periods(1))
        ! compute boundary conditions for 2D TE by solving 1D
        !  equations at boundary ... put results in EXB
        call SetBoundTE(periods(1),EXB)
        !  solve 2D TE equations, return solution in Hsol
        call Fwd2DsolveTE(EXB,Esol,IER)

        ! write solution to file
        write(3) TMgrid%Ny+1,TMgrid%Nz+1
        write(3) TMgrid%Dy
        write(3) TMgrid%Dz
        write(3) Esol

        ! deallocate
        call Fwd2DdeallTE()
        close(3)


end program
