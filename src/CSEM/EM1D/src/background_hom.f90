!**********************************************************************
!  FD EM subroutine background_hom
!
!  Purpose:  analytically compute background EM field for homogeneous space
!      new "general" routine for all instances of calling background field computations
!
!  Rita Streich 2009-2011
!**********************************************************************
subroutine background_hom_unified(bgdat,src,ifreq,icur)

  implicit none

  !external variables
  type(backgrounddata),intent(inout) :: bgdat        !EM fields at receivers
  type(sorec)          :: src       !pointer to source
  integer(kind=int32),intent(in)   :: ifreq        !index of frequency component
  integer(kind=int32),intent(in)   :: icur         !source current counter

  !internal variables
  integer(kind=int32)  :: nelem     !number of dipole or wire elements in source
  integer(kind=int32)  :: ngrp      !number of "source element groups" (i.e. number of wires for wire source, 1 for dipole source)
  integer(kind=int32)  :: igrp,ielem      !counters for "source element groups" and dipole/wire elements
  real(kind=real64)    :: omega     !frequency*2*pi
  complex(kind=real64) :: j_om_mu   !j*omega*mu0
  complex(kind=real64) :: eta       !electric medium parameter, contains conductivity and permittivity
  complex(kind=real64) :: zeta      !magnetic medium parameter, contains magnetic permeability
  complex(kind=real64) :: gamma     !propagation parameter, contains eta and zeta
  complex(kind=real64) :: etainv    !1/(4*pi*eta)
  complex(kind=real64) :: zetainv   !1/(4*pi*zeta)
  logical              :: elsrc,magsrc  !indicates if there are electric and magnetic sources
  real(kind=real64)    :: lx,ly           !x,y lengths of wire
  real(kind=real64)    :: dlwx,dlwy       !lengths of wire elements
  complex(kind=real64) :: cur       !temp source current
  real(kind=real64)    :: sx,sy,sz  !source position, does not have to be on grid here
  integer(kind=int32)  :: irec      !iReceiver counter
  real(kind=real64)    :: x,y,z     !grid point coordinates, full grid, relative to source!
  real(kind=real64)    :: dx,dy,dz  !dx,dy,dz
  real(kind=real64)    :: r         !radius, source - grid point distance
  real(kind=real64)    :: dv        !size of volume element (used only for iReceiver exactly at source point)
  complex(kind=real64),dimension(3) :: src_j,src_k        !source signatures for a single point
  complex(kind=real64),dimension(3) :: gej,gek,ghj,ghk    !lines of greens function tensors, re-used for lines 1,2,3


  !iReceiver data are collected on process 0, but computation for long wires or many dipole elements can be lengthy
  !--> share the work!

  !some initialization...
  ngrp = size(src%nelem) !number of dipole elements or wires

  !dummy cell volume, only used for the special case if source and iReceiver coincide exactly
  dv = 1._real64

  !frequency
  omega = bgdat%omega
  j_om_mu = dci * omega * dmu0


  !medium parameters: have been searched before, just re-assign them here
  eta = bgdat%sigmabghom
  zeta = cmplx(0._real64,omega*dmu0)
  gamma = sqrt(eta*zeta)
  etainv = 1._real64/(dfourpi*eta)
  zetainv = 1._real64/(dfourpi*zeta)


  !initialize fields
  bgdat%Ex = 0._real64
  bgdat%Ey = 0._real64
  bgdat%Ez = 0._real64
  bgdat%Hx = 0._real64
  bgdat%Hy = 0._real64
  bgdat%Hz = 0._real64


  do igrp=1,ngrp
    nelem = src%nelem(igrp)

    select case (src%type)
    case (dipole)
      elsrc = src%elsrc
      magsrc = src%magsrc

    case (wire)

      !x, y wire lengths
      lx = src%wire(igrp)%endpos(1,2) - src%wire(igrp)%endpos(1,1)
      ly = src%wire(igrp)%endpos(2,2) - src%wire(igrp)%endpos(2,1)

      !x, y lengths of wire elements
      dlwx = lx / real(src%nelem(igrp))
      dlwy = ly / real(src%nelem(igrp))

      !z coordinate of the entire wire
      sz = src%wire(igrp)%endpos(3,1)

      !source current for this wire
      cur = src%cur((icur-1)*ngrp+igrp,ifreq)

      !length of dipole elements times source current
      src_j(1) = dlwx * cur
      src_j(2) = dlwy * cur
      src_j(3) = 0._real64

      elsrc = .true.
      magsrc = .false.
    end select

    do ielem=1,nelem

      select case (src%type)
      case (dipole)

        !source coordinates
        sx = src%pos(1,ielem)
        sy = src%pos(2,ielem)
        sz = src%pos(3,ielem)

        cur = src%cur(ielem,ifreq)

        !source signature
        if (elsrc) then
          !source dipole length times (complex) source current
          src_j(1) = src%ljx(ielem) * cur
          src_j(2) = src%ljy(ielem) * cur
          src_j(3) = src%ljz(ielem) * cur
        endif
        if (magsrc) then
          !components of the source loop area times (complex) source current
          !factor -j*omega*mu0 as in Loeseth and Ursin, 2007
          src_k(1) = - j_om_mu * src%akx(ielem) * cur
          src_k(2) = - j_om_mu * src%aky(ielem) * cur
          src_k(3) = - j_om_mu * src%akz(ielem) * cur
        endif

      case (wire)

        !source x,y coordinates only
        sx = src%wire(igrp)%elempos(ielem,1)
        sy = src%wire(igrp)%elempos(ielem,2)

      end select

      !same coordinates for all field components?
      samecoord: if (bgdat%allcomp_samecoord) then

        !loop over receivers
        do irec = 1,bgdat%nExy

          x = bgdat%Exypos(irec,1) - sx
          y = bgdat%Exypos(irec,2) - sy
          z = bgdat%Exypos(irec,3) - sz
          r = radius0(x,y,z)

          !Ex = Gej11*J1 + Gej12*J2 + Gej13*J3  +  Gek11*K1 + Gek12*K2 + Gek13*K3

          !contribution of electric sources
          if (elsrc) then

            !put in dv in case point is right at source
            !Ex
            gej = greens_ej1(r,x,y,z,etainv,gamma,dv)
            bgdat%Ex(irec) = bgdat%Ex(irec) + sum(gej*src_j)
            !Ey
            gej = greens_ej2(r,x,y,z,etainv,gamma,dv)
            bgdat%Ey(irec) = bgdat%Ey(irec) + sum(gej*src_j)
            !Ez
            gej = greens_ej3(r,x,y,z,etainv,gamma,dv)
            bgdat%Ez(irec) = bgdat%Ez(irec) + sum(gej*src_j)
            !Hx
            ghj = greens_ek1(r,x,y,z,gamma)
            bgdat%Hx(irec) = bgdat%Hx(irec) - sum(ghj*src_j)
            !Hy
            ghj = greens_ek2(r,x,y,z,gamma)
            bgdat%Hy(irec) = bgdat%Hy(irec) - sum(ghj*src_j)
            !Hz
            ghj = greens_ek3(r,x,y,z,gamma)
            bgdat%Hz(irec) = bgdat%Hz(irec) - sum(ghj*src_j)
          endif

          !contribution of magnetic sources
          if (magsrc) then

            !Ex
            gek = greens_ek1(r,x,y,z,gamma)
            bgdat%Ex(irec) = bgdat%Ex(irec) + sum(gek*src_k)
            !Ey
            gek = greens_ek2(r,x,y,z,gamma)
            bgdat%Ey(irec) = bgdat%Ey(irec) + sum(gek*src_k)
            !Ez
            gek = greens_ek3(r,x,y,z,gamma)
            bgdat%Ez(irec) = bgdat%Ez(irec) + sum(gek*src_k)
            !Hx
            ghk = greens_ej1(r,x,y,z,zetainv,gamma,dv)
            bgdat%Hx(irec) = bgdat%Hx(irec) + sum(ghk*src_k)
            !Hy
            ghk = greens_ej2(r,x,y,z,zetainv,gamma,dv)
            bgdat%Hy(irec) = bgdat%Hy(irec) + sum(ghk*src_k)
            !Hz
            ghk = greens_ej3(r,x,y,z,zetainv,gamma,dv)
            bgdat%Hz(irec) = bgdat%Hz(irec) + sum(ghk*src_k)
          endif
        enddo  !iReceiver points

      else !different coordinates for different field components

        !Ex / Ey
        !common coordinates for Ex and Ey
        common_Exy: if (bgdat%nExy.gt.0) then

          !loop over receivers
          do irec = 1,bgdat%nExy

            x = bgdat%Exypos(irec,1) - sx
            y = bgdat%Exypos(irec,2) - sy
            z = bgdat%Exypos(irec,3) - sz
            r = radius0(x,y,z)

            !contribution of electric sources
            if (elsrc) then
              !Ex
              gej = greens_ej1(r,x,y,z,etainv,gamma,dv)
              bgdat%Ex(irec) = bgdat%Ex(irec) + sum(gej*src_j)
              !Ey
              gej = greens_ej2(r,x,y,z,etainv,gamma,dv)
              bgdat%Ey(irec) = bgdat%Ey(irec) + sum(gej*src_j)
            endif

            !contribution of magnetic sources
            if (magsrc) then
              !Ex
              gek = greens_ek1(r,x,y,z,gamma)
              bgdat%Ex(irec) = bgdat%Ex(irec) + sum(gek*src_k)
              !Ey
              gek = greens_ek2(r,x,y,z,gamma)
              bgdat%Ey(irec) = bgdat%Ey(irec) + sum(gek*src_k)
            endif
          enddo  !iReceiver points

       else

          !Ex only
          do irec = 1,bgdat%nEx

            x = bgdat%Expos(irec,1) - sx
            y = bgdat%Expos(irec,2) - sy
            z = bgdat%Expos(irec,3) - sz
            r = radius0(x,y,z)

            !contribution of electric sources
            if (elsrc) then
              gej = greens_ej1(r,x,y,z,etainv,gamma,dv)
              bgdat%Ex(irec) = bgdat%Ex(irec) + sum(gej*src_j)
            endif
            !contribution of magnetic sources
            if (magsrc) then
              gek = greens_ek1(r,x,y,z,gamma)
              bgdat%Ex(irec) = bgdat%Ex(irec) + sum(gek*src_k)
            endif
          enddo  !iReceiver points
         
          !Ey only
          do irec = 1,bgdat%nEy

            x = bgdat%Eypos(irec,1) - sx
            y = bgdat%Eypos(irec,2) - sy
            z = bgdat%Eypos(irec,3) - sz
            r = radius0(x,y,z)

            !contribution of electric sources
            if (elsrc) then
              gej = greens_ej2(r,x,y,z,etainv,gamma,dv)
              bgdat%Ey(irec) = bgdat%Ey(irec) + sum(gej*src_j)
            endif

            !contribution of magnetic sources
            if (magsrc) then
              !Ey
              gek = greens_ek2(r,x,y,z,gamma)
              bgdat%Ey(irec) = bgdat%Ey(irec) + sum(gek*src_k)
            endif
          enddo  !iReceiver points
        endif common_Exy

        !Ez
        do irec = 1,bgdat%nExy

          x = bgdat%Ezpos(irec,1) - sx
          y = bgdat%Ezpos(irec,2) - sy
          z = bgdat%Ezpos(irec,3) - sz
          r = radius0(x,y,z)

          !contribution of electric sources
          if (elsrc) then
            gej = greens_ej3(r,x,y,z,etainv,gamma,dv)
            bgdat%Ez(irec) = bgdat%Ez(irec) + sum(gej*src_j)
          endif
          !contribution of magnetic sources
          if (magsrc) then
            gek = greens_ek3(r,x,y,z,gamma)
            bgdat%Ez(irec) = bgdat%Ez(irec) + sum(gek*src_k)
          endif
        enddo  !iReceiver points


        !Hx / Hy
        !common coordinates for Hx and Hy
        common_Hxy: if (bgdat%nHxy.gt.0) then

          !loop over receivers
          do irec = 1,bgdat%nHxy

            x = bgdat%Hxypos(irec,1) - sx
            y = bgdat%Hxypos(irec,2) - sy
            z = bgdat%Hxypos(irec,3) - sz
            r = radius0(x,y,z)

            !contribution of electric sources
            if (elsrc) then
              !Hx
              ghj = greens_ek1(r,x,y,z,gamma)
              bgdat%Hx(irec) = bgdat%Hx(irec) - sum(ghj*src_j)
              !Hy
              ghj = greens_ek2(r,x,y,z,gamma)
              bgdat%Hy(irec) = bgdat%Hy(irec) - sum(ghj*src_j)
            endif

            !contribution of magnetic sources
            if (magsrc) then
              !Hx
              ghk = greens_ej1(r,x,y,z,zetainv,gamma,dv)
              bgdat%Hx(irec) = bgdat%Hx(irec) + sum(ghk*src_k)
              !Hy
              ghk = greens_ej2(r,x,y,z,zetainv,gamma,dv)
              bgdat%Hy(irec) = bgdat%Hy(irec) + sum(ghk*src_k)
            endif
          enddo  !iReceiver points

       else

          !Hx only
          do irec = 1,bgdat%nHx

            x = bgdat%Hxpos(irec,1) - sx
            y = bgdat%Hxpos(irec,2) - sy
            z = bgdat%Hxpos(irec,3) - sz
            r = radius0(x,y,z)

            !contribution of electric sources
            if (elsrc) then
              ghj = greens_ek1(r,x,y,z,gamma)
              bgdat%Hx(irec) = bgdat%Hx(irec) - sum(ghj*src_j)
            endif
            !contribution of magnetic sources
            if (magsrc) then
              ghk = greens_ej1(r,x,y,z,zetainv,gamma,dv)
              bgdat%Hx(irec) = bgdat%Hx(irec) + sum(ghk*src_k)
            endif
          enddo  !iReceiver points
         
          !Hy only
          do irec = 1,bgdat%nHy

            x = bgdat%Hypos(irec,1) - sx
            y = bgdat%Hypos(irec,2) - sy
            z = bgdat%Hypos(irec,3) - sz
            r = radius0(x,y,z)

            !contribution of electric sources
            if (elsrc) then
              ghj = greens_ek2(r,x,y,z,gamma)
              bgdat%Hy(irec) = bgdat%Hy(irec) - sum(ghj*src_j)
            endif

            !contribution of magnetic sources
            if (magsrc) then
              ghk = greens_ej2(r,x,y,z,zetainv,gamma,dv)
              bgdat%Hy(irec) = bgdat%Hy(irec) + sum(ghk*src_k)
            endif
          enddo  !iReceiver points
        endif common_Hxy

        !Hz
        do irec = 1,bgdat%nHz

          x = bgdat%Hzpos(irec,1) - sx
          y = bgdat%Hzpos(irec,2) - sy
          z = bgdat%Hzpos(irec,3) - sz
          r = radius0(x,y,z)

          !contribution of electric sources
          if (elsrc) then
            ghj = greens_ek3(r,x,y,z,gamma)
            bgdat%Hz(irec) = bgdat%Hz(irec) - sum(ghj*src_j)
          endif
          !contribution of magnetic sources
          if (magsrc) then
            ghk = greens_ej3(r,x,y,z,zetainv,gamma,dv)
            bgdat%Hz(irec) = bgdat%Hz(irec) + sum(ghk*src_k)
          endif
        enddo  !iReceiver points
      endif samecoord

    enddo  !source (i.e. dipole or wire) elements
  enddo !"source element groups"

endsubroutine background_hom_unified

