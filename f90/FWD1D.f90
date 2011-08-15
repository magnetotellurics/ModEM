program fwd1d

    use field1d
    use modelspace
    implicit none

    type(conf1d_t)                              :: earth
    type(grid_t)                                :: grid
    type(modelParam_t)                          :: model,source
    type(cvector)                               :: h1d
    character(80)                               :: period_file,label
    character(80)                               :: layered_model_file
    character(80)                               :: source_model_file
    character(80)                               :: grid_file
    character(80)                               :: fields_output_file
    character(80)                               :: cfile
    real(8), allocatable, dimension(:)          :: depths,coeff,logrho,T
    real(8)                                     :: days
    character(3)                                :: ich
    integer                                     :: i,nL,nper,lmax,Nt,Np,Nr,narg,ios,istat

    write(*,*) 'Copyright (c) 2010-2011 Oregon State University'
    write(*,*) 'College of Oceanic and Atmospheric Sciences'
    write(*,*) 'Matlab code written by Jin Sun, last mod. 24 May 2010'
    write(*,*) 'Recoded in Fortran by Anna Kelbert, 11-13 July 2011'
    write(*,*)

    !  parse command line
    narg = command_argument_count()
    if (narg < 5) then
        write(0,*) 'Usage: ./FWD1D period_file layered_model_file source_model_file grid_file fields_output_file'
        stop
    end if

    call get_command_argument(1, period_file)
    call get_command_argument(2, layered_model_file)
    call get_command_argument(3, source_model_file)
    call get_command_argument(4, grid_file)
    call get_command_argument(5, fields_output_file)

    ! save periods in days
    open(ioREAD,file=period_file,status='old',form='formatted',iostat=ios)
    write(6,*) 'Reading from the periods file ',trim(period_file)
    read(ioREAD,'(a)') label
    write(6,*) label
    read(ioREAD,*) nper
    allocate(T(nper), STAT=istat)
    do i = 1,nper
      read(ioREAD,*) days ! reading period in *days*
      T(i) = days * (24*3600)
    end do
    close(ioREAD)

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

    ! source file should only have one layer
    call read_modelParam(source,source_model_file)
    if (source%nL /= 1) then
        write(0,*) 'Error in FWD1D: source file should have exactly one layer'
        stop
    end if
    allocate(coeff(source%nc), STAT=istat)
    call getParamValues_modelParam(source,coeff)
    lmax = getDegree_modelParam(source)

    ! reading grid file (r is in km decreasing from top to bottom)
    call read_grid(grid,grid_file)

    ! set earth radius and domain top radius (in meters)
    earth%r0  = 6371.0e3
    earth%rmax= 1.0e3 * grid%r(1)

    ! set tolerance on toroidal potential
    earth%tol = 1.e-9

    ! set surface conductance (should be small since we're using 3D thinsheet)
    earth%tau = 1.e2

    ! save model in 1D configuration structure: layers include the core
    allocate(earth%layer(nL+1),earth%sigma(nL+1), STAT=istat)
    earth%layer(1) = 0.0d0
    earth%layer(2:nL+1) = 1.0e3 * depths(1:nL)
    earth%sigma(1:nL) = 10.0**( -logrho(1:nL) )
    earth%sigma(nL+1) = 10.0**( 5.0 ) ! core conductivity

    write(*,*) 'Tops of model layers: ',earth%layer
    write(*,*) 'Conductivity values:  ',earth%sigma

    ! allocate the output cvector
    call create_cvector(grid, h1d, EDGE)

    do i = 1,nper
        days = T(i)/(24*3600)
        write(*,*) 'Computing the fields for period ',trim(ich),': ',days,' days'
        call sourceField1d(earth,lmax,coeff,T(i),grid,h1d)

        call reset_time(fwd1d_timer)

        write(ich,'(i3.3)') i
        cfile = trim(fields_output_file)//'_'//trim(ich)//'.field'
        write(*,*) 'Writing to file: ',cfile
        open(ioWRITE,file=cfile,status='unknown',form='formatted',iostat=ios)
        write(ioWRITE,'(a45,f9.3,a6)') "# FWD1D full EM field solution output for period ",   &
                                            days,' days.'
        write(ioWRITE,'(i3)') 1
        call write_cvector(ioWRITE,h1d)
        close(ioWRITE)

        write(*,*) 'Done writing to file: ',elapsed_time(fwd1d_timer),' secs'
    end do

    deallocate(depths,coeff,logrho,T, STAT=istat)
    deallocate(earth%layer,earth%sigma, STAT=istat)
    call deall_modelParam(model)
    call deall_modelParam(source)
    call deall_cvector(h1d)
    call deall_grid(grid)
    write(*,*) 'Total time taken: ',saved_time(fwd1d_timer),' secs'

end program fwd1d
