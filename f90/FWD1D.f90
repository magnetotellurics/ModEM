program fwd1d

    use field1d
    implicit none

    type(conf1d_t)                              :: earth
    complex(8), allocatable, dimension(:,:,:)   :: Hp,Ht,Hr
    real(8), allocatable, dimension(:)          :: R
    real(8), dimension(4)                       :: coef
    integer                                     :: nL,lmax,Nt,Np,Nr,istat

    !earth radius and domain top radius
    earth%r0 = 6371.0e3; earth%rmax=7.0e7;

    ! number of 1D layers
    nL = 12

    allocate(earth%layer(nL),earth%sigma(nL), STAT=istat)

    !depth of layered earth at the top of each layer
    earth%layer = (/0., 100.e3, 250.e3, 410.e3, 670.e3, 900.e3, 1200.e3, 1600.e3, 2000.e3, 2400.e3, 2900.e3, 3500.e3/);

    !conductivity of each layer
    earth%sigma = (/-2.2, -2.0, -1.3, -1.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 5.0/);
    earth%sigma = 10.0**earth%sigma

    !surface conductance
    earth%tau = 1.e3;

    !period of the field
    earth%T = 1*3600.;

    !tolarence on toroidal potential
    earth%tol=1.e-9;

    !max degree considered
    earth%Nmax=100;

    !spherical harmonic sources
    lmax = 1;
    coef = (/0., 1., 0., 0./);
    !Ns = [1 2];
    !coef = {[1,0,0], [0,0,1,0,0,1,0]};

    !# of grid points in co-latitude (0~180 degrees)
    !Nt = 90;
    Nt = 36;

    !# of grid points in longitude (0~360 degress)
    !Np = 90;
    Np = 72;

    !test run
    !Nr = 2;
    !allocate(R(Nr), STAT=istat)
    !R =  (/5371000.,  10833900./);

    Nr = 11;
    allocate(R(Nr), STAT=istat)
    R =  (/5371000.,  10833900.,  16296800.,  21759700.,  27222600., &
          32685500.,  38148400.,  43611300.,  49074200.,  54537100.,  60000000./);

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

    allocate(Hp(Np,Nt,Nr),Ht(Np,Nt,Nr),Hr(Np,Nt,Nr), STAT=istat)


    call sourceField1d(earth,lmax,coef,Np,Nt,R,Hp,Ht,Hr)

    deallocate(R, STAT=istat)
    deallocate(Hp,Ht,Hr, STAT=istat)
    deallocate(earth%layer,earth%sigma, STAT=istat)

end program fwd1d
