!**********************************************************************
!>  FD EM subroutines for homogeneous space Green's functions
!
!>  Purpose           :  analytic expressions for homogeneous space greens functions
!>                       --> fields due to infinitesimal dipole sources
!
!>  the following are 3x3 Green's function tensors:
!>  Gej: Greens functions for the electric field due to electric sources
!>  Gek: Greens functions for the electric field due to magnetic sources
!>  Ghj: Greens functions for the magnetic field due to electric sources
!>  Ghk: Greens functions for the magnetic field due to magnetic sources
!
!>  subdivision:
!>    greens_ej1: components of Gej needed to compute E1, i.e. G11,G12,G13
!>    greens_ej2: components of Gej needed to compute E2
!>    greens_ej3: components of Gej needed to compute E3 etc.
!
!>  functions are subdivided since they are called separately because of staggered grid,
!>    E1,E2,E3 are all needed on slightly different coordinates
!
!>  the formulas are published manyfold, but in the form used they are taken from:
!>    Jan van der Kruk, Three-dimensional imaging of multicomponent ground 
!>    penetrating radar data, PhD thesis, TU Delft, the Netherlands, 2001
!>    (equations 4.15a to 4.15d)
!
!>  we have Gek = - Ghj and eta*Gej = zeta*Ghk
!>    --> the same functions are used for computing Gek and Ghj, the minus
!>        for Ghj has to be taken care of outside these functions
!>    --> the same functions are used for computing Gej and Ghk, input
!>        parameters etainv or zetainv have to be specified
!
!>  Rita Streich 2009
!
!**********************************************************************
function greens_ej1(r,x,y,z,etainv,gamma,dv)

  implicit none

  !external variables
  complex(kind=real64),dimension(3) :: greens_ej1
  real(kind=real64) :: r       !distance from source
  real(kind=real64) :: x,y,z   !coordinates relative to source point
  complex(kind=real64) :: etainv  !electric medium parameter, contains conductivity and permittivity
  complex(kind=real64) :: gamma   !propagation parameter, contains eta and zeta
  real(kind=real64) :: dv      !size of volume element

  !internal variables
  complex(kind=real64) :: expo    !exponential term
  real(kind=real64) :: r2,r3,r4,r5  !powers of r
  real(kind=real64) :: x2      !power of x
  real(kind=real64) :: xy,xz
  complex(kind=real64) :: gam2    !power of gamma


  !special case for source point
  !equations follow Lee, S.W., Boersma, J., Law, C.-L. and Deschamps, G.A., 1980:
  !>  Singularity in Green's function and its numerical evaluation, IEEE Transactions on
  !>  Antennas and propagation AP-28(3), 311-317
  !and PhD thesis of J. van der Kruk
  !but not sure if taking the limit for a point source from formula for source volume was correct...
  if(r .EQ. 0._real64) then
    !remember etainv = 1/(4 pi eta)
    greens_ej1(1) = ((gamma**3 - dfourpi/3._real64) * etainv) / dv
    greens_ej1(2) = (( - dfourpi/3._real64) * etainv) / dv !0._real64
    greens_ej1(3) = (( - dfourpi/3._real64) * etainv) / dv !0._real64

  else

    r2 = r*r
    r3 = r2*r
    r4 = r3*r
    r5 = r4*r
    x2 = x*x
    xy = x*y
    xz = x*z
    expo = exp(-gamma*r) * etainv !factor 4*pi is included in etainv
    gam2 = gamma*gamma

    !remember damping coefficients are already in form 1/ex --> multiply them
    greens_ej1(1) = expo * (-gam2/r  + (- gamma/r2 + (gam2*x2-1)/r3 + (3._real64*gamma*x2)/r4 + (3._real64*x2)/r5))
    greens_ej1(2) = expo * (gam2*xy/r3 + 3._real64*gamma*xy/r4 + 3._real64*xy/r5)
    greens_ej1(3) = expo * (gam2*xz/r3 + 3._real64*gamma*xz/r4 + 3._real64*xz/r5)
  endif

  return
end function greens_ej1


function greens_ej2(r,x,y,z,etainv,gamma,dv)

  implicit none

  !external variables
  complex(kind=real64),dimension(3) :: greens_ej2
  real(kind=real64) :: r       !distance from source
  real(kind=real64) :: x,y,z   !coordinates relative to source point
  complex(kind=real64) :: etainv  !electric medium parameter, contains conductivity and permittivity
  complex(kind=real64) :: gamma   !propagation parameter, contains eta and zeta
  real(kind=real64) :: dv      !size of volume element

  !internal variables
  complex(kind=real64) :: expo    !exponential term
  real(kind=real64) :: r2,r3,r4,r5  !powers of r
  real(kind=real64) :: y2      !power of y
  real(kind=real64) :: xy,yz
  complex(kind=real64) :: gam2    !power of gamma


  !special case for source point
  if(r .EQ. 0._real64) then
    !remember etainv = 1/(4 pi eta)
    greens_ej2(1) = (( - dfourpi/3._real64) * etainv) / dv !0._real64
    greens_ej2(2) = ((gamma**3 - dfourpi/3._real64) * etainv) / dv
    greens_ej2(3) = (( - dfourpi/3._real64) * etainv) / dv !0._real64

  else

    r2 = r*r
    r3 = r2*r
    r4 = r3*r
    r5 = r4*r
    y2 = y*y
    xy = x*y
    yz = y*z
    expo = exp(-gamma*r) * etainv !factor 4*pi is included in etainv
    gam2 = gamma*gamma

    greens_ej2(1) = expo * (gam2*xy/r3 + 3._real64*gamma*xy/r4 + 3._real64*xy/r5)
    greens_ej2(2) = expo * (-gam2/r + (- gamma/r2 + (gam2*y2-1)/r3 + (3._real64*gamma*y2)/r4 + (3._real64*y2)/r5))
    greens_ej2(3) = expo * (gam2*yz/r3 + 3._real64*gamma*yz/r4 + 3._real64*yz/r5)
  endif

  return
end function greens_ej2


function greens_ej3(r,x,y,z,etainv,gamma,dv)

  implicit none

  !external variables
  complex(kind=real64),dimension(3) :: greens_ej3
  real(kind=real64) :: r       !distance from source
  real(kind=real64) :: x,y,z   !coordinates relative to source point
  complex(kind=real64) :: etainv  !electric medium parameter, contains conductivity and permittivity
  complex(kind=real64) :: gamma   !propagation parameter, contains eta and zeta
  real(kind=real64) :: dv      !size of volume element

  !internal variables
  complex(kind=real64) :: expo    !exponential term
  real(kind=real64) :: r2,r3,r4,r5  !powers of r
  real(kind=real64) :: z2      !power of y
  real(kind=real64) :: xz,yz
  complex(kind=real64) :: gam2    !power of gamma


  !special case for source point
  if(r .EQ. 0._real64) then
    !remember etainv = 1/(4 pi eta)
    greens_ej3(1) = (( - dfourpi/3._real64) * etainv) / dv !0._real64
    greens_ej3(2) = (( - dfourpi/3._real64) * etainv) / dv !0._real64
    greens_ej3(3) = ((gamma**3 - dfourpi/3._real64) * etainv) / dv

  else

    r2 = r*r
    r3 = r2*r
    r4 = r3*r
    r5 = r4*r
    z2 = z*z
    xz = x*z
    yz = y*z
    expo = exp(-gamma*r) * etainv !factor 4*pi is included in etainv
    gam2 = gamma*gamma

    greens_ej3(1) = expo * (gam2*xz/r3 + 3._real64*gamma*xz/r4 + 3._real64*xz/r5)
    greens_ej3(2) = expo * (gam2*yz/r3 + 3._real64*gamma*yz/r4 + 3._real64*yz/r5)
    greens_ej3(3) = expo * (-gam2/r + (- gamma/r2 + (gam2*z2-1)/r3 + (3._real64*gamma*z2)/r4 + (3._real64*z2)/r5))
  endif

  return
end function greens_ej3


function greens_ek1(r,x,y,z,gamma)

  implicit none

  !external variables
  complex(kind=real64),dimension(3) :: greens_ek1
  real(kind=real64) :: r       !distance from source
  real(kind=real64) :: x,y,z     !coordinates relative to source point
  complex(kind=real64) :: gamma   !propagation parameter, contains eta and zeta

  !internal variables
  complex(kind=real64) :: expo    !exponential term
  real(kind=real64) :: r2,r3   !powers of r


  !special case for source point
  if(r .EQ. 0._real64) then
    greens_ek1 = 0._real64
  else

    r2 = r*r
    r3 = r2*r
    expo = exp(-gamma*r) / dfourpi

    greens_ek1(1) = 0._real64
    greens_ek1(2) = - expo * (gamma*z/r2 + z/r3)
    greens_ek1(3) = expo * (gamma*y/r2 + y/r3)
  endif

  return
end function greens_ek1


function greens_ek2(r,x,y,z,gamma)

  implicit none

  !external variables
  complex(kind=real64),dimension(3) :: greens_ek2
  real(kind=real64) :: r       !distance from source
  real(kind=real64) :: x,y,z     !coordinates relative to source point
  complex(kind=real64) :: gamma   !propagation parameter, contains eta and zeta

  !internal variables
  complex(kind=real64) :: expo    !exponential term
  real(kind=real64) :: r2,r3   !powers of r


  !special case for source point
  if(r .EQ. 0._real64) then
    greens_ek2 = 0._real64
  else

    r2 = r*r
    r3 = r2*r
    expo = exp(-gamma*r) / dfourpi

    greens_ek2(1) = expo * (gamma*z/r2 + z/r3)
    greens_ek2(2) = 0._real64
    greens_ek2(3) = - expo * (gamma*x/r2 + x/r3)
  endif

  return
end function greens_ek2


function greens_ek3(r,x,y,z,gamma)

  implicit none

  !external variables
  complex(kind=real64),dimension(3) :: greens_ek3
  real(kind=real64) :: r       !distance from source
  real(kind=real64) :: x,y,z     !coordinates relative to source point
  complex(kind=real64) :: gamma   !propagation parameter, contains eta and zeta

  !internal variables
  complex(kind=real64) :: expo    !exponential term
  real(kind=real64) :: r2,r3   !powers of r


  !special case for source point
  if(r .EQ. 0._real64) then
    greens_ek3 = 0._real64
  else

    r2 = r*r
    r3 = r2*r
    expo = exp(-gamma*r) / dfourpi

    greens_ek3(1) = - expo * (gamma*y/r2 + y/r3)
    greens_ek3(2) = expo * (gamma*x/r2 + x/r3)
    greens_ek3(3) = 0._real64
  endif

  return
end function greens_ek3

