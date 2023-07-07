!****************************************************************
!
!  FD EM module for computing background field for 1D medium
!         using reflectivity method
!
!  Programmed: Rita Streich 2009
!
!**************************************************************

module refl_mod_new

  use machine_dep_mod
  use util_mod
  use io_common_mod
  use hankel_mod
  use adaptive_int_mod
#ifdef USE_MPI
  use counter_mod
#endif

  implicit none

  private

  !public subroutines
  public :: reflectivity_unified, conjugate_fields

  !public structure
  public :: refl_struct


  !container for all variables that have to be remembered while computing 1D fields
  type refl_struct
    !source description
    integer(kind=int32)                      :: nzsrc      !number of unique source element depths within one source
    integer(kind=int32),dimension(:),pointer :: isrcperz   !indices for source elements at each depth
    real(kind=real64),dimension(:),pointer   :: betasrc    !angles for HED or HMD source elements
    real(kind=real64),dimension(:),pointer   :: zsrc       !unique source (element) depths
    integer(kind=int32),dimension(:),pointer :: nsrcperz   !how many source elements at each source depth
    integer(kind=int32)                      :: isrcstart,isrcend !indices for source range for one depth
    real(kind=real64),dimension(:),pointer   :: xs,ys      !horizontal positions of dipole or wire elements

    !receiver description
    !number of unique receiver depths within one receiver group, separate for each field component
    !combine Exy and Hxy because they require the same integral evaluations
    integer(kind=int32)                      :: nzrecExy,nzrecEz,nzrecHxy,nzrecHz
    !indices for receivers at each depth, separate for different field components
    integer(kind=int32),dimension(:),pointer :: irecperzExy,irecperzEz,irecperzHxy,irecperzHz
    !unique receiver depths, separate for different field components
    real(kind=real64),dimension(:),pointer   :: zrecExy,zrecEz,zrecHxy,zrecHz
    !number of receivers at each receiver depth, separate for different field components
    integer(kind=int32),dimension(:),pointer :: nrecperzExy,nrecperzEz,nrecperzHxy,nrecperzHz

    integer(kind=int32)                      :: irecstart,irecend !indices for source range for one depth
    real(kind=real64),dimension(:),pointer   :: xr,yr      !horizontal positions of receivers

    !integrals
    integer(kind=int32)                      :: nrad       !number of radius values
    real(kind=real64)                        :: rmax,rmin  !max and min (non-zero) radii
    real(kind=real64),dimension(:),pointer       :: radlog    !radii, logarithmically spaced for Hankel transforms
    complex(kind=real64),dimension(:,:),pointer  :: intvaltmp !temp integral values for one zs-zr combination,dimension nr*NREL

    !need to split integral values because spline interpolation only takes reals
    !dimension: nrad * (up to 10) since we need up to 10 integrals, depending on source type
    !re-use for Exy, Ez, Hxy, Hz
    real(kind=real64),dimension(:,:),pointer :: intvalre      !real part of integral values
    real(kind=real64),dimension(:,:),pointer :: intvalim      !imag part of integral values
    real(kind=real64),dimension(:,:),pointer :: spl_derivre   !derivatives of spline function, real part
    real(kind=real64),dimension(:,:),pointer :: spl_derivim   !derivatives of spline function, imag part
    type(receiverdata),dimension(:),pointer  :: EHwire        !for wire sources: separate EM fields for each wire
    !for wire sources: separate EM fields for derivatives for each wire, for isotropic/horizontal conductivity, nlayers x nwires
    type(receiverdata),dimension(:,:),pointer  :: EHwirederiv 
    type(receiverdata),dimension(:,:),pointer  :: EHwirederivv !for wire sources: separate EM fields for derivatives for each wire (VERTICAL)

    logical                                  :: refcoef_changed  !indicates if reflection coeff. need to be recomputed
    integer(kind=int32)                      :: infolevel     !determined the amount of log output from 1D computations
  end type refl_struct

  !number of integrals for different source types
  integer(kind=int32),parameter,public  :: nintHED=10, nintVED=3, nintHMD=10, nintVMD=3, nintWire=6
!!$  integer(kind=int32),parameter,public  :: nintHED_Exy=4, nintHED_Ez=1
  integer(kind=int32),parameter,public  :: nintHEDd=11, nintVEDd=6, nintHMDd=11, nintVMDd=3, nintWired=7 !for isotropic derivatives
  !for vti derivatives
  integer(kind=int32),parameter,public  :: nintHEDdvti=16, nintVEDdvti=9, nintHMDdvti=16, nintVMDdvti=3, nintWiredvti=10
  !will allocate space for max number of integrals ever required - the overhead is negligible
!!$  integer(kind=int32),parameter,public  :: nintHED=16, nintVED=9, nintHMD=16, nintVMD=3, nintWire=10


  integer(kind=int32)               :: nlay     !number of layers
  integer(kind=int32)               :: aniso    !anisotropy index

  real(kind=real64)                 :: omega    !angular frequency
  real(kind=real64)                 :: omegasq  !angular frequency squared
  complex(kind=real64)              :: jomega   !j*omega
  complex(kind=real64)              :: j_om_mu  !j*omega*mu0
  complex(kind=real64),dimension(:),allocatable :: epsv,epsh, epsmuv,epsmuh  !effective permittivity, product of medium properties
  complex(kind=real64),dimension(:),allocatable :: epsmuratio  !epsmuh / epsmuv
  complex(kind=real64),dimension(:),pointer     :: pvert !vertical slowness for all layers
  complex(kind=real64),dimension(:),allocatable :: omsq_epsmuv,omsq_epsmuh  !omega**2 * epsilon * mu0
  real(kind=real64),dimension(:),allocatable    :: branchpt !projections of branch points for all layers to real axis
  integer(kind=int32)               :: ilaysrc,ilayrec !index of the layers where source/receivers are in

  complex(kind=real64),dimension(:),pointer :: rup  !elm. 11 or 22 of upgoing refl. coeff. for all layer boundaries
  complex(kind=real64),dimension(:),pointer :: rdn  !elm. 11 or 22 of downgoing refl. coeff. for all layer boundaries
  complex(kind=real64),dimension(:),pointer :: tup  !elm. 11 or 22 of upgoing transm. coeff. for all layer boundaries
  complex(kind=real64),dimension(:),pointer :: tdn  !elm. 11 or 22 of downgoing transm. coeff. for all layer boundaries
  complex(kind=real64),dimension(:,:),pointer :: rupallTE  !TE upgoing refl. coeff. for all layer boundaries and all wavenumbers
  complex(kind=real64),dimension(:,:),pointer :: rupallTM  !TM upgoing refl. coeff. for all layer boundaries and all wavenumbers
  complex(kind=real64),dimension(:,:),pointer :: rdnallTE  !TE downgoing refl. coeff. for all layer boundaries, all wavenumbers
  complex(kind=real64),dimension(:,:),pointer :: rdnallTM  !TM downgoing refl. coeff. for all layer boundaries, all wavenumbers
  complex(kind=real64),dimension(:,:),pointer :: tupallTE  !TE upgoing transm. coeff. for all layer boundaries, all wavenumbers
  complex(kind=real64),dimension(:,:),pointer :: tupallTM  !TM upgoing transm. coeff. for all layer boundaries, all wavenumbers
  complex(kind=real64),dimension(:,:),pointer :: tdnallTE  !TE downgoing transm. coeff. for all layer boundaries, all wavenumbers
  complex(kind=real64),dimension(:,:),pointer :: tdnallTM  !TM downgoing transm. coeff. for all layer boundaries, all wavenumbers
  complex(kind=real64),dimension(:,:),pointer :: pvertall1,pvertall2     !vertical wavenumbers for all layers


  !dipole source types
  integer(kind=int32),parameter   :: hed = 11
  integer(kind=int32),parameter   :: ved = 12
  integer(kind=int32),parameter   :: hmd = 13
  integer(kind=int32),parameter   :: vmd = 14
  
  real(kind=real64),dimension(:),allocatable    :: dz_below_src, dz_above_src !layer thicknesses below/above sources
  real(kind=real64),dimension(:),allocatable    :: dz_below_rec, dz_above_rec !layer thicknesses below/above receivers
  complex(kind=real64),dimension(:),allocatable :: trans_above_src  !transmission through layers above source
  complex(kind=real64),dimension(:),allocatable :: trans_below_src  !transmission through layers below source
  complex(kind=real64),dimension(:),allocatable :: trans_above_rec  !transmission through layers above receivers
  complex(kind=real64),dimension(:),allocatable :: trans_below_rec  !transmission through layers below receivers

  integer(kind=int32)  :: ilaym          !layer counter for derivatives, take dE/depsilon_m
  real(kind=real64)    :: kappatmp       !temp storage of wavenumber for testing derivatives
  real(kind=real64)                 :: rsplmin  !for radii smaller than this, do not use spline interpolation,
                                                ! but compute integrals for exact radius values		


			
contains
#include "quadpack_dbl.f90"
  include 'reflectivity.f90'
  include 'init_refcoef.f90'
  include 'find_recdepths.f90'
#include "find_radii.f90"
  include 'init_intval.f90'
  include 'conjugate_fields.f90'
  include 'find_srcdepths.f90'
#include "get_dz.f90"
  include 'getilay.f90'
  include 'precompute_intval.f90'
  include 'precomp_intvals_deriv.f90'
  include 'recursion.f90'
  include 'recursion_deriv.f90'
  include 'interp_1val.f90'
  include 'interp_intvals_HED.f90'
  include 'interp_intvals_VED.f90'
  include 'interp_intvals_HMD.f90'
  include 'interp_intvals_VMD.f90'
  include 'interp_intvals_wire.f90'
  include 'wirefields.f90'
  include 'addwirefields.f90'

endmodule refl_mod_new
