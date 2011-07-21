program fwd1d

    use field1d
    use modelspace
    use solnspace
    implicit none

    type(conf1d_t)                              :: earth
    type(userdef_control)                       :: userdef
    type(grid_t)                                :: grid
    type(modelParam_t)                          :: model,source
    type(cvector)                               :: h1d
    character(80)                               :: period_file !,label
    character(80)                               :: layered_model_file
    character(80)                               :: source_model_file
    character(80)                               :: grid_file
    character(80)                               :: fields_output_file
    !complex(8), allocatable, dimension(:,:,:)   :: Hp,Ht,Hr
    !real(8), allocatable, dimension(:)          :: R
    real(8), allocatable, dimension(:)          :: depths,coeff,logrho,days
    !integer, parameter                          :: ioREAD=10,ioWRITE=11
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
        return
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
    allocate(days(nper), STAT=istat)
    do i = 1,nper
      read(ioREAD,*) days(i) ! reading period in *days*
    end do
    close(ioREAD)

    ! model file should be 1D layered
    call read_modelParam(model,layered_model_file)
    if (model%nL /= model%nc) then
        write(0,*) 'Error in FWD1D: input model file is not layered 1D'
        return
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
        return
    end if
    allocate(coeff(source%nc), STAT=istat)
    call getParamValues_modelParam(source,coeff)
    lmax = getDegree_modelParam_f(source)

    ! reading grid file (r is in km decreasing from top to bottom)
    call read_grid(grid,grid_file)

    ! set earth radius and domain top radius (in meters)
    earth%r0  = 6371.0e3
    earth%rmax= 1.0e3 * grid%r(1) + 1.0e0

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



    ! number of 1D layers
    !nL = 12
    !depth of layered earth at the top of each layer
    !earth%layer = (/0., 100.e3, 250.e3, 410.e3, 670.e3, 900.e3, 1200.e3, 1600.e3, 2000.e3, 2400.e3, 2900.e3, 3500.e3/);
    !conductivity of each layer
    !earth%sigma = (/-2.2, -2.0, -1.3, -1.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 5.0/);
    !earth%sigma = 10.0**earth%sigma
    !period of the field
    !earth%T = 1*3600.;


    ! period of the field
    !earth%T = days(i)

    !spherical harmonic sources
    !lmax = 1;
    !coef = (/0., 1., 0., 0./);
    !Ns = [1 2];
    !coef = {[1,0,0], [0,0,1,0,0,1,0]};

    !# of grid points in co-latitude (0~180 degrees)
    !Nt = 90;
    !Nt = 36;

    !# of grid points in longitude (0~360 degress)
    !Np = 90;
    !Np = 72;

    !test run
    !Nr = 2;
    !allocate(R(Nr), STAT=istat)
    !R =  (/5371000.,  10833900./);

!    Nr = 11;
!    allocate(R(Nr), STAT=istat)
!    R =  (/5371000.,  10833900.,  16296800.,  21759700.,  27222600., &
!          32685500.,  38148400.,  43611300.,  49074200.,  54537100.,  60000000./);

    !    R = ( 19113.1 19113.0 14563. 10467. 8419. 7395. 6883. &
    !        6627. 6499. 6435. 6403. 6387. 6379. 6375. 6373. 6372. 6371.5 &
    !        6371.0 6370.5 6370. 6368. 6366. 6364. 6362. 6360.5 6358.35 &
    !        6321. 6271. 6221. 6171. 6121. 6071. 6021. &
    !        5991. 5961. 5931. 5901. 5871. 5851. &
    !        5831. 5801. 5771. 5741. 5711. 5701. &
    !        5671. 5621. 5571. 5521. 5471. &
    !        5421. 5371. 5321. 5271. 5221. 5171. 5071. 4971. &
    !        4871. 4771. 4671. 4571. 4471. &
    !        4371. 4271. 4171. 4071. 3971. &
    !        3871. 3771. 3671. 3571. &
    !        3481. 3471. 3371. 3271. 3171. &
    !        3071. 2971. 2921. 2900. 2871.5 2871. );

    call create_cvector(grid, h1d, EDGE)

    !allocate(Hp(Np,Nt+1,Nr+1),Ht(Np+1,Nt,Nr+1),Hr(Np+1,Nt+1,Nr), STAT=istat)

    !allocate(Hp(Np+1,Nt+1,Nr+1),Ht(Np+1,Nt+1,Nr+1),Hr(Np+1,Nt+1,Nr+1), STAT=istat)
    i = 1
    call sourceField1d(earth,lmax,coeff,days(i),grid,h1d)

    write(*,*) 'Writing to file: ',trim(fields_output_file)
    open(ioWRITE,file=fields_output_file,status='unknown',form='formatted',iostat=ios)
    write(ioWRITE,'(a45,f9.3,a6)') "# Full EM field solution output for period ",   &
                                        days(i),' days.'
    write(ioWRITE,'(i3)') 1
    call write_cvector(ioWRITE,h1d)
    close(ioWRITE)

    !deallocate(R, STAT=istat)
    !deallocate(Hp,Ht,Hr, STAT=istat)
    deallocate(depths,coeff,logrho,days, STAT=istat)
    deallocate(earth%layer,earth%sigma, STAT=istat)
    call deall_modelParam(model)
    call deall_modelParam(source)
    call deall_cvector(h1d)
    call deall_grid(grid)

end program fwd1d
