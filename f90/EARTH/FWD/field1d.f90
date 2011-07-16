! *****************************************************************************
module field1d
	! Implements 1D analytic forward modelling

  use polpak
  implicit none

Contains

subroutine vsharm(n,cost,phi,Y,Yt,YpDst)
! vector spherical harmonic components
! Ynm=sqrt((2n+1)/4pi*(n-m)!/(n+m)!)*Pnm(cost)(exp(imphi))

	integer, intent(in)		            :: n
	real(8), intent(in)		            :: cost,phi
	real(8), intent(inout)		        :: Y
	real(8), intent(inout), optional    :: Yt,YpDst

    Y=diag([1,(-1).^(1:n)/sqrt(2)])*legendre(n,cost,'sch');

    if (.not. present(Yt)) then
        return
    end if

    Yt = zeros(size(Y));
    if n>0
        Yt(1,:)=sqrt((n)*(n+1))*Y(2,:);
        do m=1,n-1
            Yt(m+1,:)=(sqrt((n-m)*(n+m+1))*Y(m+2,:)...
                       -sqrt((n+m)*(n-m+1))*Y(m,:))/2;
        end do

        Yt(n+1,:)=-sqrt(n/2)*Y(n,:);

        if nargin==3 && ~isempty(phi)
            ep=exp(1i*phi);
            do m=1,n
                Y(m+1,:)=ep.^m.*Y(m+1,:);
                Yt(m+1,:)=ep.^m.*Yt(m+1,:);
            end do
        end
    end

    if (.not. present(YpDst)) then
        return
    end if

    YpDst=zeros(size(Y));
    if (n>0) then
        ids=find(abs(cost)==1);
        do m=1,n-1
            YpDst(m+1,:)=1i*m*Y(m+1,:)./sqrt(1-cost.^2);

            if(.not. isempty(ids)) then
                YpDst(m+1,ids)=sqrt((n-m)*(n+m+1))*(Y(m+2,ids)...
                  +sqrt((n+m)*(n-m+1))*Y(m,ids))./(2i*cost(ids));
            end if
        end do
        YpDst(n+1,1:end)=1i*n*Y(n+1,1:end)./sqrt(1-cost(1:end).^2);
        if (.not. isempty(ids)) then
            YpDst(n+1,ids)=sqrt(2*n)*Y(n,ids)./(2i*cost(ids));
        end if
    end if

end subroutine

subroutine airprop(N,r0,rn0,rnp0,r,rni,rnip,rnr,rnrp,sumup)
!airprop propagates into the air(r0) 
!given the BC at the surface r0
!rn0,rnp0 have to have the dimension of N

	integer, intent(in), dimension(:)		:: N
	real(8), intent(in), dimension(:)		:: r0,rn0,rnp0,r
	real(8), intent(inout), dimension(:)	:: rni,rnip,rnr,rnrp
    logical, intent(in)                     :: sumup

    lr = size(r)
    ln = size(N)
    if (size(rn0) /= ln) then
        write(0,*) 'Error: rn0 is a scalar in airprop'
    end if
    if (size(rnp0) /= ln) then
        write(0,*) 'Error: rnp0 is a scalar in airprop'
    end if

    if (lr>0 .and. ln>0) then

        rni = zeros(lr,ln);
        rnip = zeros(lr,ln);
        rnr = zeros(lr,ln);
        rnrp = zeros(lr,ln);
        rr = r/r0; rrn=1; n=0;

        do j=1,ln

            rrn=rrn.*rr.^(N(j)-n); n=N(j);
            rn1  = rrn.*rr;
            rn1p = (n+1)*rn1./r;
            rn2  = 1./(rrn);
            rn2p = -n*rn2./r;

            A1 = (n*rn0(j)+r0*rnp0(j))/(2*n+1);

            A2 = ((n+1)*(rn0(j)-A1)...
                  -r0*(rnp0(j)-A1*(n+1)/r0))/(2*n+1);

            rni(:,j) = A1*rn1;
            rnr(:,j) = A2*rn2;
            rnip(:,j) = A1*rn1p;
            rnrp(:,j) = A2*rn2p;

        end do
    else
        rni=[];  rnip=[];  rnr=[];  rnrp=[];
    end if

    if (sumup) then
        rni=rni+rnr; rnip=rnip+rnrp;
    end

end subroutine

subroutine rbsls0(N,z0,z,type,Rbn,Rbnp)
!Computes the ricatti bessel function and 
!the derivative scaled by exp(-+abs(imag(z0)))
!and by (2n+-1)!!^+-1
!and by z0^-+n
!to lessen possible overflow and underflow

	integer, intent(in)		:: N
	real(8), intent(in)		:: z0,z
	integer, intent(in)		:: type
	real(8), intent(inout)		:: Rbn,Rbnp

    LN=length(N); Nm=max(N);
    Lz=numel(z);
    Rbn =zeros(Lz,LN);
    Rbnp=zeros(Lz,LN);
    
    select case type
        
        case(1)

            call rrbessel(1:Nm,z,RN);
            Sn = 1/2i*(exp(1i*z-imag(z0))...
                       -exp(-1i*z-imag(z0)));
            z=1./z;
            c=1;
            do n=1,Nm
                Sn0=(2*n+1)/z0*Sn;
                Sn =Sn0.*RN(:,n);
                if (ismember(n,N)) then
                    Rbn(:,c) =Sn;
                    Rbnp(:,c)=Sn0-n*Sn.*z;
                    c=c+1;
                end if
            end do

        case(2)

            Sn=-1i*exp(1i*z+imag(z0));
            z=1./z;
            Rn=z-1i;
            c=1;
            do n=1,Nm
                Sn0 = z0/(2*n+1)*Sn;
                Sn = Sn0.*Rn;
                if (ismember(n,N)) then
                    Rbn(:,c) = Sn;
                    Rbnp(:,c)= Sn0-n*Sn.*z;
                    c=c+1;
                end if
                Rn=(2*n+1)*z-1./Rn;
            end do
    end select

end subroutine

subroutine rrbessel(N,z,RN)
! ratio of riccati bessel functions 
! between adjascent orders rn=SH_(n)/SH_(n-1)

	integer, intent(in)		:: N
	real(8), intent(in)		:: z
	real(8), intent(inout)		:: RN

    RN=zeros(length(z(:)),length(N));

    z=1./z(:);
    rn=1./contFrac(@(n)(-1)^n*(2*max(N)+1+2*n)*z);
    c=length(N);
    do n=max(N),min(N),-1
        if (ismember(n,N)) then
            RN(:,c)=rn; c=c-1;
        end if
        rn=1./((2*n-1).*z-rn);
    end do

end subroutine

function    contFrac(a) result(f)
!continued fraction

    real(8), intent(in)		:: a(:)
    real(8)				:: f
    ! local
    real(8)             :: tiny
    logical             :: conv

    tiny=1e-30;
    f=a(0); f(abs(f)<tiny)=tiny;
    g=f; h=zeros(size(f));
    r=ones(size(f)); idn=find(r);
    conv=.false.
    k=1
    do
        ak=a(k); akidn=ak(idn);
    
        g(idn)=akidn+1./g(idn);  g(abs(g)<tiny)=tiny;
        h(idn)=akidn+1.*h(idn);  h(abs(h)<tiny)=tiny;
        h(idn)=1./h(idn);
        r(idn)=g(idn).*h(idn);
        f(idn)=f(idn).*r(idn);

        idn=find(abs(r-1)>1e-14);
        if (isempty(idn)) then
            conv=.true.
            exit
        end if
        k=k+1;
    end do

end function 


subroutine rbslprop(N,z0,phn0,phnp0,z,phn,phnp)
!propagate ricatti bessel function from z0 to z
!given its value and derivative at z0

	integer, intent(in)		:: N
	real(8), intent(in)		:: z0,phn0,phnp0
	real(8), intent(inout)		:: z,phn,phnp

    call rbsls0(N,z0,[z0;z(:)],1,Sn,Snp);
    call rbsls0(N,z0,[z0;z(:)],2,Cn,Cnp);

    Sn0=Sn(1,:);
    Cn0=Cn(1,:);
    Snp0=Snp(1,:);
    Cnp0=Cnp(1,:);

    wr = -1i;
        An = (Cnp0.*phn0-Cn0.*phnp0).*wr;

        phn0h=phn0-An.*Sn0;
        phnp0h=phnp0-An.*Snp0;

        Bn = (Sn0.*(phnp0h)-Snp0.*(phn0h)).*wr;

    
    phn = Sn(2:end,:)*diag(An)...
         +Cn(2:end,:)*diag(Bn);

    phnp = Snp(2:end,:)*diag(An)...
          +Cnp(2:end,:)*diag(Bn);

end subroutine


subroutine sourceField1d(earth,Ns,coef,Np,Nt,R,Hp,Ht,Hr)
!in Matlab, optionally shift to mid-faces

	type (conf1d_t), intent(in)	:: earth
	integer, intent(in)		:: Ns,Np,Nt
	real(8), intent(in)		:: coef(:),R
	real(8), intent(inout)		:: Hp,Ht,H
	! local
	real(8), dimension(:,:), allocatable    :: Tnr,Tnsp
	real(8)				:: mu0,omega
    integer             :: Nrr,Nrs,Ls,istat

    !Physical model
    mu0 = 1.256637e-6 !(H/m)

    lr=numel(earth%layer)
    if (lr .ne. numel(earth%sigma)) then
      write(0,*) 'Error in sourceField1d: Number of elements in ld and lc must be equal'
    end if
    if (earth%layer(1) .ne. 0) then
      write(0,*) 'Error in sourceField1d: Depth of 1st layer must be zero'
    end if

    rl = earth%r0 - earth%layer
    omega = 2*pi/earth%T
    kl = sqrt(1i*omega*mu0*earth%sigma)
    !-----------------------------------------------------------!

    !Radial response
    Rs=R(:)
    Nrs=size(Rs)
    Rr=Rs
    Nrr=Nrs
    Ls=size(Ns)

    allocate(Tnr(1:Nrr,1:Ls),Tnsp(1:Nrs,1:Ls),STAT=istat)

    !within the inner core
    idr=find(Rr<=rl(lr));
    if (.not. isempty(idr)) then
        call rbsls0(Ns,kl(lr)*rl(lr),kl(lr)*Rr(idr),1,tnrll,tmptnsp);
        Tnr(idr,:)=tnrll;
    end if

    ids=find(Rs<=rl(lr));
    if (.not. isempty(ids)) then
        call rbsls0(Ns,kl(lr)*rl(lr),kl(lr)*Rs(ids),1,tmptnr,tnspll);
        Tnsp(ids,:)=kl(lr)*tnspll;
    end if
    
    [tn,tnp]=rbsls0(Ns,kl(lr)*rl(lr),kl(lr)*rl(lr),1);
    tnp = kl(lr)*tnp./tn;
    
    Tnr=Tnr*diag(1./tn);
    Tnsp=Tnsp*diag(1./tn);

    !between layers
    do ll=lr-1,1,-1

        idr=find(Rr<=rl(ll)&Rr>rl(ll+1));
        if (.not. isempty(idr)) then
            call rbslprop(Ns,kl(ll)*rl(ll+1),1,tnp/kl(ll),kl(ll)*Rr(idr),tnrll,tmptnsp);
            Tnr(idr,:)=tnrll;
        end if
        
        ids=find(Rs<=rl(ll)&Rs>rl(ll+1));
        if (.not. isempty(ids)) then
            call rbslprop(Ns,kl(ll)*rl(ll+1),1,tnp/kl(ll),kl(ll)*Rs(ids),tmptnr,tnspll);
            Tnsp(ids,:)=kl(ll)*tnspll;
        end if

        call rbslprop(Ns,kl(ll)*rl(ll+1),1,tnp/kl(ll),kl(ll)*rl(ll),tn,tnp);

        tnp = kl(ll)*tnp./tn;

        Tnr=Tnr*diag(1./tn);
        Tnsp=Tnsp*diag(1./tn);

    end do

    !in the air
    idr=find(Rr>rl(1));
    if (.not. isempty(idr)) then
        sumup = .true.
        call airprop(Ns,rl(1),1,tnp-1i*omega*mu0*earth%tau,Rr(idr),tnrll,tmptnsp,tmp,tmp,sumup);
        Tnr(idr,:)=tnrll;
    end if

    ids=find(Rs>rl(1));
    if (.not. isempty(ids)) then
        sumup = .true.
        call airprop(Ns,rl(1),1,tnp-1i*omega*mu0*earth%tau,Rs(ids),tmptnr,tnspll,tmp,tmp,sumup);
        Tnsp(ids,:)=tnspll;
    end if
    !renormalize against outter boundary
    sumup = .false.
    call airprop(1,rl(1),1,tnp(1)-1i*omega*mu0*earth%tau,earth%rmax,tmptnr,tlmp,tmp,tmp,sumup);

    Tnr=Tnr*(-earth%rmax/t1mp);
    Tnsp=Tnsp*(-earth%rmax/t1mp);
    !-----------------------------------------------------------!

    !field componets
    theta=linspace(0,pi,Nt+1); theta(end)=[];
    phi=linspace(0,2*pi,Np+1); phi(end)=[];

    dt=pi/Nt;
    dp=pi/Np;

    [p,t]=ndgrid(phi,theta);

    Hp=0;
    Ht=0;
    Hr=0;

    do l=1,Ls

        coefl=coef{l}(:);

        call vsharm(Ns(l),cos(t(:)),p(:)+dp/2,tmp,tmp,Yp);

        Hps=(Yp.'*coefl(1:Ns(l)+1)+Yp(2:end,:)'*coefl(Ns(l)+2:end))*(Tnsp(:,l)./Rs).';

        call vsharm(Ns(l),cos(t(:)+dt/2),p(:),tmp,Yt);

        Hts=(Yt.'*coefl(1:Ns(l)+1)+Yt(2:end,:)'*coefl(Ns(l)+2:end))*(Tnsp(:,l)./Rs).';

        call vsharm(Ns(l),cos(t(:)),p(:),Yr);

        Hrs=-Ns(l)*(Ns(l)+1)*(Yr.'*coefl(1:Ns(l)+1) + Yr(2:end,:)'*coefl(Ns(l)+2:end))*(Tnr(:,l)./(Rr).^2).';

        Hp=Hp+Hps;
        Ht=Ht+Hts;
        Hr=Hr+Hrs;
    end do

    Hp=reshape(conj(Hp),Np,Nt,[]);
    Ht=reshape(conj(Ht),Np,Nt,[]);
    Hr=reshape(conj(Hr),Np,Nt,[]);

    deallocate(Tnr,Tnsp,STAT=istat)

end subroutine

end module field1d
