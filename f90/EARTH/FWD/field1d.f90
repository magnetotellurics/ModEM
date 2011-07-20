! *****************************************************************************
module field1d
	! Implements 1D analytic forward modelling

  use utilities
  use polpak
  use griddef
  use sg_vector
  implicit none

  ! ***************************************************************************
  ! * type conf1d_t contains the configuration for the 1d forward solver

  type :: conf1d_t

      real(8), dimension(:), pointer    :: layer ! tops of layers
      real(8), dimension(:), pointer    :: sigma ! layer conductivities
      real(8)                           :: tau ! near-surface conductance
      real(8)                           :: T ! period in seconds
      real(8)                           :: tol ! tolerance
      real(8)                           :: r0 ! Earth's radius in meters
      real(8)                           :: rmax ! boundary of the domain in meters
      logical                           :: allocated

  end type conf1d_t ! conf1d_t


Contains

subroutine vsharm(lmax,cost,phi,Y,Yt,YpDst)
! vector spherical harmonic components
! Ynm=sqrt((2n+1)/4pi*(n-m)!/(n+m)!)*Pnm(cost)(exp(imphi))

	integer, intent(in)		                                    :: lmax
	real(8), intent(in)		                                    :: cost,phi
	complex(8), dimension(lmax,lmax+1), intent(inout)		    :: Y
	complex(8), dimension(lmax,lmax+1), intent(inout), optional :: Yt,YpDst
	! local
	real(8), dimension(lmax,lmax+1)     :: P_lm
	real(8)     :: tiny
	integer     :: l,m
	complex(8)  :: ep

	if (lmax <= 0) then
	    write(0,*) 'Error in vsharm: this only works for lmax >= 1'
	end if

    !compute Schmidt Seminormalized Associated Legendre Functions, defined as:
    !Q_n^m(x) = P_n(x) for m=0
    !Q_n^m(x) = (-1)^m sqrt ( 2 * (n-m)! / (n+m)! ) P_n^m for m>0
    P_lm(:,:) = 0.0d0
    call legendre_associated(lmax,0,cost,P_lm(:,1))
	do m = 1,lmax
        call legendre_associated(lmax,m,cost,P_lm(:,m+1))
        do l = m,lmax
            P_lm(l,m+1) = P_lm(l,m+1) * ( (-1)**m ) * &
                sqrt ( 2.0d0 * ( d_factorial(l - m) ) / ( d_factorial(l + m) ) )
        end do
    end do

    !now convert this to vector spherical harmonics
    Y(:,:) = cmplx(0.0d0,0.0d0)
    Y(:,1) = P_lm(:,1)
    do l = 1,lmax
        do m = 1,l
            Y(l,m+1) = ((-1)**(m)/sqrt(2.0d0)) * P_lm(l,m+1)
        end do
    end do

    if (.not. present(Yt)) then
        return
    end if

    Yt(:,:) = cmplx(0.0d0,0.0d0)
    do l = 1,lmax

        Yt(l,1)=sqrt((l)*(l+1.0d0))*Y(l,2);
        do m=1,l-1
            Yt(l,m+1)=(sqrt((l-m)*(l+m+1.0d0))*Y(l,m+2) - sqrt((l+m)*(l-m+1.0d0))*Y(l,m))/2
        end do
        Yt(l,l+1)=-sqrt(l/2.0d0)*Y(l,l)

        ep=exp(cmplx(0.0d0,1.0d0)*phi);
        do m=1,l
            Y(l,m+1)=(ep**m)*Y(l,m+1);
            Yt(l,m+1)=(ep**m)*Yt(l,m+1);
        end do

    end do

    if (.not. present(YpDst)) then
        return
    end if

    YpDst(:,:) = cmplx(0.0d0,0.0d0)
    tiny=1e-30
    do l = 1,lmax

        !presumably, YpDst(:,1) = 0
        do m = 1,l-1
            if (abs(abs(cost)-1.0d0) < tiny) then
                YpDst(l,m+1)=sqrt((l-m)*(l+m+1.0d0))*(Y(l,m+2) &
                  +sqrt((l+m)*(l-m+1.0d0))*Y(l,m))/(cmplx(0.0d0,2.0d0)*cost)
            else
                YpDst(l,m+1)=cmplx(0.0d0,1.0d0)*m*Y(l,m+1)/sqrt(1.0d0-cost**2)
            end if
        end do
        if (abs(abs(cost)-1.0d0) < tiny) then
            YpDst(l,l+1)=sqrt(2.0d0*l)*Y(l,l)/(cmplx(0.0d0,2.0d0)*cost)
        else
            YpDst(l,l+1)=cmplx(0.0d0,1.0d0)*l*Y(l,l+1)/sqrt(1.0d0-cost**2)
        end if

    end do

end subroutine

subroutine airprop(lmax,r0,rn0,rnp0,r,rni,rnip,sumup)
!airprop propagates into the air(r0) 
!given the BC at the surface r0
!rn0,rnp0 have to have the dimension of number of sph. harm. degrees

	integer, intent(in)             		        :: lmax
    real(8), intent(in)                             :: r0,r
	complex(8), dimension(lmax), intent(in)         :: rn0,rnp0
	complex(8), dimension(lmax), intent(inout)	    :: rni,rnip
    logical, intent(in)                             :: sumup
    ! local
    complex(8), dimension(lmax)   :: rnr,rnrp
    complex(8)                      :: rr,rrn,rn1,rn1p,rn2,rn2p
    complex(8)                      :: A1,A2
    integer                         :: l,istat

    if (size(rn0) /= (lmax)) then
        write(0,*) 'Error in airprop: rn0 should have dimension of lmax'
    end if
    if (size(rnp0) /= (lmax)) then
        write(0,*) 'Error in airprop: rnp0 should have dimension of lmax'
    end if

    rni(:) = 0
    rnip(:) = 0
    rnr(:) = 0
    rnrp(:) = 0

    rr = r/r0; rrn=1

    do l=1,lmax

        rrn  = rrn*rr
        rn1  = rrn*rr
        rn1p = (l+1)*rn1/r
        rn2  = 1./(rrn)
        rn2p = -l*rn2/r

        A1 = (l*rn0(l)+r0*rnp0(l))/(2*l+1)

        A2 = ((l+1)*(rn0(l)-A1) - r0*(rnp0(l)-A1*(l+1)/r0))/(2*l+1)

        rni(l) = A1*rn1
        rnr(l) = A2*rn2
        rnip(l) = A1*rn1p
        rnrp(l) = A2*rn2p

    end do

    if (sumup) then
        rni=rni+rnr
        rnip=rnip+rnrp
    end if

end subroutine

subroutine rbsls0(lmax,z0,z,type,Rbl,Rblp)
!Computes the ricatti bessel function and
!the derivative scaled by exp(-+abs(imag(z0)))
!and by (2l+-1)!!^+-1
!and by z0^-+l
!to lessen possible overflow and underflow

	integer, intent(in)		                :: lmax
	complex(8), intent(in)		            :: z0,z
	integer, intent(in)		                :: type
	complex(8), dimension(:), intent(inout)	:: Rbl,Rblp
	! local
    complex(8), dimension(:), allocatable   :: Rl
	complex(8)                              :: Sl,Sl0
    complex(8)                              :: cone,ctwo
	integer                                 :: l,istat

	cone = cmplx(0.0d0,1.0d0)
	ctwo = cmplx(0.0d0,2.0d0)

    if ((size(Rbl) .ne. (lmax)) .or. (size(Rblp) .ne. (lmax))) then
        write(0,*) 'Error in rbsls0: output array sizes must be equal to the number of sph. harm. degrees'
    end if
    Rbl(:) = 0
    Rblp(:)= 0

    select case (type)

        case(1)

            allocate(Rl(lmax), STAT=istat)
            call rrbessel(lmax,z,Rl)
            !write(*,*) 'Rl=',Rl
            Sl = 1./ctwo*(exp(cone*z-dimag(z0))-exp(-cone*z-dimag(z0)))
            do l=1,lmax
                Sl0=(2*l+1)/z0*Sl
                Sl =Sl0*Rl(l)
                Rbl(l) =Sl
                Rblp(l)=Sl0-l*Sl/z
            end do
            deallocate(Rl, STAT=istat)

        case(2)

            allocate(Rl(1), STAT=istat)
            Sl=-cone*exp(cone*z+dimag(z0))
            Rl(1)=(1./z)-cone
            do l=1,lmax
                Sl0 = z0/(2*l+1)*Sl
                Sl = Sl0*Rl(1)
                Rbl(l) = Sl
                Rblp(l)= Sl0-l*Sl/z
                Rl(1)=((2*l+1)/z)-1./Rl(1)
            end do
            deallocate(Rl, STAT=istat)
    end select

end subroutine


subroutine rrbessel(lmax,z,Rl)
! ratio of riccati bessel functions 
! between adjascent orders rl=SH_(l)/SH_(l-1)

	integer, intent(in)		                :: lmax
	complex(8), intent(in)		            :: z
	complex(8), dimension(:), intent(inout)	:: Rl
	! local
	complex(8)  :: rl1
	integer     :: l

    if (size(Rl) .ne. (lmax)) then
        write(0,*) 'Error in rrbessel: output array size must be equal to the number of sph. harm. degrees'
    end if

    rl1=1./contFrac(a,lmax,z)
    do l=lmax,1,-1
        Rl(l)=rl1
        rl1=1./((2*l-1)/z-rl1)
    end do

end subroutine

function contFrac(a,lmax,z) result (f)
!continued fraction

    complex(8), external            	    :: a
    integer, intent(in)                     :: lmax
    complex(8), intent(in)                  :: z
    complex(8)                              :: f
    ! local
    complex(8)              :: g,h,r,ak
    real(8)                 :: tiny
    logical                 :: conv
    integer                 :: k

    tiny=1e-30;
    f = a(0,lmax,z)
    if (abs(f) < tiny) then
        f = tiny
    end if
    g = f
    h = 0
    r = 1
    conv=.false.
    k=1
    do
        ak=a(k,lmax,z)
        g=ak+1./g
        if (abs(g) < tiny) then
            g = tiny
        end if
        h=ak+1.*h
        if (abs(h) < tiny) then
            h = tiny
        end if
        h=1./h
        r=g*h
        f=f*r

        if (abs(r-1)<=1.e-14) then
            conv=.true.
            exit
        end if
        k=k+1;
    end do

end function


complex(8) function a(n,lmax,z)

    integer, intent(in)           :: n,lmax
    complex(8), intent(in)        :: z

    a = ((-1)**n)*(2*lmax+1+2*n)/z

end function


subroutine rbslprop(lmax,z0,phn0,phnp0,z,phn,phnp)
!propagate ricatti bessel function from z0 to z
!given its value and derivative at z0

	integer, intent(in)		    :: lmax
	complex(8), intent(in)		:: z0,z
	complex(8), dimension(lmax), intent(in)    :: phn0,phnp0
	complex(8), dimension(lmax), intent(inout)	:: phn,phnp
    ! local
    complex(8), dimension(lmax) :: Sn,Snp,Cn,Cnp
    complex(8), dimension(lmax) :: Sn0,Snp0,Cn0,Cnp0
    complex(8), dimension(lmax) :: phn0h,phnp0h,An,Bn
    complex(8)  :: wr
    integer     :: l,istat

    call rbsls0(lmax,z0,z0,1,Sn0,Snp0)
    call rbsls0(lmax,z0,z0,2,Cn0,Cnp0)
    !write(*,*) 'Sn0 type 1:',Sn0
    !write(*,*) 'Cn0 type 2:',Cn0

    wr = -cmplx(0.,1.);

    do l = 1,lmax
        An(l) = (Cnp0(l)*phn0(l)-Cn0(l)*phnp0(l))*wr;

        phn0h(l)=phn0(l)-An(l)*Sn0(l);
        phnp0h(l)=phnp0(l)-An(l)*Snp0(l);

        Bn(l) = (Sn0(l)*(phnp0h(l))-Snp0(l)*(phn0h(l)))*wr;
    end do

    call rbsls0(lmax,z0,z,1,Sn,Snp)
    call rbsls0(lmax,z0,z,2,Cn,Cnp)
    
    do l = 1,lmax
        phn(l) = Sn(l)*An(l)+Cn(l)*Bn(l)
        phnp(l) = Snp(l)*An(l)+Cnp(l)*Bn(l)
    end do

end subroutine

subroutine sourcePotential(earth,lmax,Rr,Rs,Tnr,Tnsp)
!in Matlab, optionally shift to mid-faces

    type (conf1d_t), intent(in)         :: earth ! configuration structure
    integer, intent(in)                 :: lmax ! maximum sph. harm. degree
    real(8), dimension(:), intent(in)   :: Rr,Rs ! radii for vertical and lateral potentials
    complex(8), dimension(:,:), intent(inout)  :: Tnr,Tnsp ! potential coefficients
    ! local
    integer, dimension(lmax)    :: Ns
    complex(8), dimension(:), allocatable      :: rn0,rnp0,phn0,phnp0
    complex(8), dimension(:), allocatable      :: tnr1,tnsp1,tn,tnp,tlmp,tmp
    real(8), dimension(:), allocatable         :: rl
    complex(8), dimension(:), allocatable      :: kl
    real(8)             :: mu0,pi,omega,rmax
    integer             :: idr,idrmin,idrmax,ids,idsmin,idsmax
    integer             :: i,j,k,l,m,ll,Nl,Nrr,Nrs,istat
    logical             :: sumup

    !Physical model
    mu0 = 1.256637e-6 !(H/m)
    pi = 3.14159265357898

    !No computations are performed for l=0
    if (lmax <= 0) then
        write(0,*) 'Error in sourcePotential: no potentials exist for degree 0 (no magnetic monopoles)'
        return
    end if

    Nl=size(earth%layer) ! number of layers
    if (size(earth%sigma) .ne. Nl) then
      write(0,*) 'Error in sourcePotential: Number of elements in arrays layer and sigma must be equal'
    end if
    if (earth%layer(1) > 0) then
      write(0,*) 'Error in sourcePotential: Depth of 1st layer must be zero'
    end if

    allocate(rl(Nl),kl(Nl),STAT=istat)
    omega = 2*pi/earth%T
    do j = 1,Nl
        rl(j) = earth%r0 - earth%layer(j)
        kl(j) = sqrt(cmplx(0.,1.)*omega*mu0*earth%sigma(j))
    end do
    rmax = earth%rmax
    !-----------------------------------------------------------!

    !Radial response
    Nrs=size(Rs)
    Nrr=size(Rr)

    allocate(tnr1(lmax),tnsp1(lmax),tn(lmax),tnp(lmax),tlmp(lmax),tmp(lmax),STAT=istat)
    allocate(rn0(lmax),rnp0(lmax),phn0(lmax),phnp0(lmax),STAT=istat)

    !within the inner core
    call find_index(Rr,0.d0,rl(Nl),idrmin,idrmax)
    !write(*,*) 'Layer ',Nl,'(Core): ',idrmin,idrmax
    if ((idrmin > 0) .and. (idrmax > 0)) then
        do idr=idrmin,idrmax
            call rbsls0(lmax,kl(Nl)*rl(Nl),kl(Nl)*Rr(idr),1,tnr1,tmp)
            Tnr(idr,:)=tnr1(:)
            !write(*,*) 'core,idr,tnr1: ',Nl,idr,tnr1
        end do
    end if

    call find_index(Rs,0.d0,rl(Nl),idsmin,idsmax)
    !write(*,*) 'Layer ',Nl,'(Core): ',idsmin,idsmax
    if ((idsmin > 0) .and. (idsmax > 0)) then
        do ids=idsmin,idsmax
            call rbsls0(lmax,kl(Nl)*rl(Nl),kl(Nl)*Rs(ids),1,tmp,tnsp1)
            Tnsp(ids,:)=kl(Nl)*tnsp1(:)
            !write(*,*) 'core,ids,tnsp1: ',Nl,ids,tnsp1
        end do
    end if

    call rbsls0(lmax,kl(Nl)*rl(Nl),kl(Nl)*rl(Nl),1,tn,tnp)
    !write(*,*) 'tn: ',tn
    !write(*,*) 'tnp: ',tnp
    do i = 1,lmax
        tnp(i) = kl(Nl)*tnp(i)/tn(i)
    end do

    do i = 1,lmax
        Tnr(:,i)=Tnr(:,i)/tn(i)
        Tnsp(:,i)=Tnsp(:,i)/tn(i)
    end do

    !between layers
    do ll=Nl-1,1,-1

        do i = 1,lmax
            phn0(i) = cmplx(1.0d0,0.0d0)
            phnp0(i) = tnp(i)/kl(ll)
        end do

        call find_index(Rr,rl(ll+1),rl(ll),idrmin,idrmax)
        !write(*,*) 'Layer ',ll,': ',idrmin,idrmax
        if ((idrmin > 0) .and. (idrmax > 0)) then
            do idr=idrmin,idrmax
                call rbslprop(lmax,kl(ll)*rl(ll+1),phn0,phnp0,kl(ll)*Rr(idr),tnr1,tmp)
                Tnr(idr,:)=tnr1(:)
                !write(*,*) 'll,idr,tnr1: ',ll,idr,tnr1
            end do
        end if

        call find_index(Rs,rl(ll+1),rl(ll),idsmin,idsmax)
        !write(*,*) 'Layer ',ll,': ',idsmin,idsmax
        if ((idsmin > 0) .and. (idsmax > 0)) then
            do ids=idsmin,idsmax
                call rbslprop(lmax,kl(ll)*rl(ll+1),phn0,phnp0,kl(ll)*Rs(ids),tmp,tnsp1)
                Tnsp(ids,:)=kl(ll)*tnsp1(:)
                !write(*,*) 'll,ids,tnsp1: ',ll,ids,tnsp1
            end do
        end if

        call rbslprop(lmax,kl(ll)*rl(ll+1),phn0,phnp0,kl(ll)*rl(ll),tn,tnp)

        do i = 1,lmax
            tnp(i) = kl(ll)*tnp(i)/tn(i)
        end do

        do i = 1,lmax
            Tnr(:,i)=Tnr(:,i)/tn(i)
            Tnsp(:,i)=Tnsp(:,i)/tn(i)
        end do

    end do

    !in the air
    do i = 1,lmax
        rn0(i) = cmplx(1.0d0,0.0d0)
        rnp0(i) = tnp(i)-cmplx(0.0d0,1.0d0)*omega*mu0*earth%tau
    end do

    call find_index(Rr,rl(1),rmax,idrmin,idrmax)
    !write(*,*) 'Layer 1 (Air): ',idrmin,idrmax
    if ((idrmin > 0) .and. (idrmax > 0)) then
        do idr=idrmin,idrmax
            sumup = .true.
            call airprop(lmax,rl(1),rn0,rnp0,Rr(idr),tnr1,tmp,sumup)
            Tnr(idr,:)=tnr1(:)
            !write(*,*) 'air,idr,tnr1: ',1,idr,tnr1
        end do
    end if

    call find_index(Rs,rl(1),rmax,idsmin,idsmax)
    !write(*,*) 'Layer 1 (Air): ',idsmin,idsmax
    if ((idsmin > 0) .and. (idsmax > 0)) then
        do ids=idsmin,idsmax
            sumup = .true.
            call airprop(lmax,rl(1),rn0,rnp0,Rs(ids),tmp,tnsp1,sumup)
            Tnsp(ids,:)=tnsp1(:)
            !write(*,*) 'air,ids,tnsp1: ',1,ids,tnsp1
        end do
    end if

    !renormalize against outer boundary
    sumup = .false.
    call airprop(1,rl(1),rn0(1),rnp0(1),earth%rmax,tmp,tlmp,sumup)

    do i = 1,lmax
        Tnr(:,i)=Tnr(:,i)*(-earth%rmax/tlmp(i))
        Tnsp(:,i)=Tnsp(:,i)*(-earth%rmax/tlmp(i))
    end do
    !-----------------------------------------------------------!
    write(*,*) 'Tnr: '
    do j = 1,Nrr
        write(*,*) Tnr(j,:)
    end do
    write(*,*) 'Tnsp: '
    do j = 1,Nrs
        write(*,*) Tnsp(j,:)
    end do
    !-----------------------------------------------------------!

    deallocate(tnr1,tnsp1,tn,tnp,tlmp,tmp,STAT=istat)
    deallocate(rn0,rnp0,phn0,phnp0,STAT=istat)
    deallocate(rl,kl,STAT=istat)

end subroutine


subroutine sourceField1d(earth,lmax,coeff,Np,Nt,R,Hp,Ht,Hr)
!in Matlab, optionally shift to mid-faces

	type (conf1d_t), intent(in)	        :: earth ! configuration structure
	integer, intent(in)		            :: lmax ! maximum sph. harm. degree
	real(8), dimension(:), intent(in)   :: coeff ! vector of sph. harm. coeff.
    integer, intent(in)                 :: Np,Nt ! lateral dimensions of the grid
	real(8), dimension(:), intent(in)	:: R ! radii for the grid
	complex(8), dimension(:,:,:), intent(inout):: Hp,Ht,Hr
	! local
	real(8), dimension(size(R)) :: Rr,Rs
	real(8), dimension(Np+1)    :: ph
	real(8), dimension(Nt+1)    :: th
    integer, dimension(lmax)    :: Ns
    complex(8), dimension(lmax,lmax+1)         :: Yp,Yt,Yr ! indices (l,m+1), m=0,..,lmax
	complex(8), dimension(:,:), allocatable    :: Tnr,Tnsp
    real(8), dimension(:), allocatable         :: coefl
	real(8)				:: dp,dt
    integer             :: idr,idrmin,idrmax,ids,idsmin,idsmax
    integer             :: i,j,k,l,m,Nrr,Nrs,Nd,istat,ncoeff,icoeff
    complex(8)          :: C

    !No computations are performed for l=0 (no magnetic monopoles) so zero coeff is never used
    ncoeff=0
    do l=0,lmax
      ncoeff = ncoeff + (2*l+1)
    end do
    if (size(coeff) .ne. ncoeff) then
        write(0,*) 'Error in sourceField1d: bad sph. harm. coeffs vector (size ',size(coeff),'); must be size ',ncoeff
    end if

    !-----------------------------------------------------------!

    !compute source potentials
    Rs=R(:) ! radii for lateral components
    Nrs=size(Rs)
    Rr=Rs ! radii for vertical components
    Nrr=Nrs
    Nd=lmax ! total number of degrees in sph. harm. expansion

    allocate(Tnr(Nrr,Nd),Tnsp(Nrs,Nd),STAT=istat)
    
    call sourcePotential(earth,lmax,Rr,Rs,Tnr,Tnsp)
    
    Hp=0
    Ht=0
    Hr=0

    !compute field components
    dt=pi/Nt
    dp=2*pi/Np

    do j = 1,Nt+1
        th(j) = dt*(j-1)
    end do
    do i = 1,Np+1
        ph(i) = dp*(i-1)
    end do
    write(*,*) 'Theta: ',th*180/pi
    write(*,*) 'Phi: ',ph*180/pi

    Yp(:,:) = cmplx(0.0d0,0.0d0)
    Yt(:,:) = cmplx(0.0d0,0.0d0)
    Yr(:,:) = cmplx(0.0d0,0.0d0)

    !compute vector spherical harmonics for all degrees and orders at one node
    do j = 1,Nt+1
        do i = 1,Np+1

            !Yp at longitudinal mid-edges
            call vsharm(lmax,cos(th(j)),ph(i)+dp/2,Yr,Yt,Yp)

            !Yt at latitudinal mid-edges
            call vsharm(lmax,cos(th(j)+dt/2),ph(i),Yr,Yt)

            !Yr at vertical mid-edges
            call vsharm(lmax,cos(th(j)),ph(i),Yr)

        end do
    end do

    ! ph component of the field
    do j = 1,Nt+1
        do i = 1,Np
            do k = 1,Nrs

                !ignore l=0 (no magnetic monopoles)
                icoeff = 1

                do l = 1,lmax

                    allocate(coefl(2*l+1), STAT=istat)

                    !ordered m=0,1,-1,2,-2 etc (NOT like in Matlab)
                    coefl = coeff(icoeff+1:icoeff+2*l+1)

                    !m=0 goes first
                    Hp(i,j,k) = Yp(l,1)*coefl(1)*(Tnsp(k,l)/Rs(k))

                    do m = 1,l
                        C = (Yp(l,m+1)*coefl(2*m) + Yp(l,m+1)*coefl(2*m+1))
                        Hp(i,j,k) = Hp(i,j,k) + C*(Tnsp(k,l)/Rs(k))
                    end do

                    icoeff = icoeff+2*l+1

                    deallocate(coefl, STAT=istat)

                end do ! degrees

            end do ! r
        end do ! ph
    end do ! th

    ! th component of the field
    do j = 1,Nt
        do i = 1,Np+1
            do k = 1,Nrs

                !ignore l=0 (no magnetic monopoles)
                icoeff = 1

                do l = 1,lmax

                    allocate(coefl(2*l+1), STAT=istat)

                    !ordered m=0,1,-1,2,-2 etc (NOT like in Matlab)
                    coefl = coeff(icoeff+1:icoeff+2*l+1)

                    !m=0 goes first
                    Ht(i,j,k) = Yt(l,1)*coefl(1)*(Tnsp(k,l)/Rs(k))

                    do m = 1,l
                        C = (Yt(l,m+1)*coefl(2*m) + Yt(l,m+1)*coefl(2*m+1))
                        Ht(i,j,k) = Ht(i,j,k) + C*(Tnsp(k,l)/Rs(k))
                    end do

                    icoeff = icoeff+2*l+1

                    deallocate(coefl, STAT=istat)

                end do ! degrees

            end do ! r
        end do ! ph
    end do ! th

    ! vertical component of the field
    do j = 1,Nt+1
        do i = 1,Np+1
            do k = 1,Nrr

                !ignore l=0 (no magnetic monopoles)
                icoeff = 1

                do l = 1,lmax

                    allocate(coefl(2*l+1), STAT=istat)

                    !ordered m=0,1,-1,2,-2 etc (NOT like in Matlab)
                    coefl = coeff(icoeff+1:icoeff+2*l+1)

                    !m=0 goes first
                    Hr(i,j,k) = -l*(l+1)*Yr(l,1)*coefl(1)*(Tnr(k,l)/Rr(k)**2)

                    do m = 1,l
                        C = - l*(l+1)*(Yr(l,m+1)*coefl(2*m) + Yr(l,m+1)*coefl(2*m+1))
                        Hr(i,j,k) = Hr(i,j,k) + C*(Tnr(k,l)/Rr(k)**2)
                    end do

                    icoeff = icoeff+2*l+1

                    deallocate(coefl, STAT=istat)

                end do ! degrees

            end do ! r
        end do ! ph
    end do ! th

    deallocate(Tnr,Tnsp,STAT=istat)

end subroutine

end module field1d
