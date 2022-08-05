!****************************************************************
!
!  Module containing stuff (somehow) related to I/O and
!    required by more than one modeling code
!
!  Rita Streich 2011
!**************************************************************

module io_common_mod

  use machine_dep_mod
  use util_mod
  use constants_mod
  use spline_mod

  implicit none

  !private

  ! Public subroutines
  public :: receiverdata, backgrounddata, freqdata, wavelet, sorec, wirespec

  type receiverdata
    !final EM data for output
    complex(kind=real64),dimension(:),pointer :: Ex
    complex(kind=real64),dimension(:),pointer :: Ey
    complex(kind=real64),dimension(:),pointer :: Ez
    complex(kind=real64),dimension(:),pointer :: Hx
    complex(kind=real64),dimension(:),pointer :: Hy
    complex(kind=real64),dimension(:),pointer :: Hz
  end type receiverdata


  !general structure for background model, background data and coordinates at which bg data are computed
  !should be "universal" for homogeneous and 1d background and all places where background fields are computed
  type backgrounddata
    real(kind=real64),dimension(:,:),pointer   :: Exypos   !positions where both Ex and Ey are computed 
                                                   !  (Ex and Ey use the same integrals, so we gain efficiency if we combine them)
    real(kind=real64),dimension(:,:),pointer   :: Expos    !positions where only Ex is computed
    real(kind=real64),dimension(:,:),pointer   :: Eypos    !positions where only Ey is computed
    real(kind=real64),dimension(:,:),pointer   :: Ezpos    !positions where Ez is computed
    real(kind=real64),dimension(:,:),pointer   :: Hxypos   !positions where both Hx and Hy are computed 
                                                   !  (Hx and Hy use the same integrals, so we gain efficiency if we combine them)
    real(kind=real64),dimension(:,:),pointer   :: Hxpos    !positions where only Hx is computed
    real(kind=real64),dimension(:,:),pointer   :: Hypos    !positions where only Hy is computed
    real(kind=real64),dimension(:,:),pointer   :: Hzpos    !positions where Hz is computed
    complex(kind=real64),dimension(:),pointer  :: Ex       !Ex data
    complex(kind=real64),dimension(:),pointer  :: Ey       !Ey data
    complex(kind=real64),dimension(:),pointer  :: Ez       !Ez data
    complex(kind=real64),dimension(:),pointer  :: Hx       !Hx data
    complex(kind=real64),dimension(:),pointer  :: Hy       !Hy data
    complex(kind=real64),dimension(:),pointer  :: Hz       !Hz data
    complex(kind=real64),dimension(:,:),pointer  :: dExdm  !Ex sensitivities, isotropic or with respect to HORIZONTAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dEydm  !Ey sensitivities, isotropic or with respect to HORIZONTAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dEzdm  !Ez sensitivities, isotropic or with respect to HORIZONTAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dHxdm  !Hx sensitivities, isotropic or with respect to HORIZONTAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dHydm  !Hy sensitivities, isotropic or with respect to HORIZONTAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dHzdm  !Hz sensitivities, isotropic or with respect to HORIZONTAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dExdmv  !Ex sensitivities, isotropic or with respect to VERTICAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dEydmv  !Ey sensitivities, isotropic or with respect to VERTICAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dEzdmv  !Ez sensitivities, isotropic or with respect to VERTICAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dHxdmv  !Hx sensitivities, isotropic or with respect to VERTICAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dHydmv  !Hy sensitivities, isotropic or with respect to VERTICAL cond.
    complex(kind=real64),dimension(:,:),pointer  :: dHzdmv  !Hz sensitivities, isotropic or with respect to VERTICAL cond.
    integer(kind=int32)                        :: nExy     !number of points at which both Ex and Ey are computed
    integer(kind=int32)                        :: nEx      !number of points at which only Ex is computed
    integer(kind=int32)                        :: nEy      !number of points at which only Ey is computed
    integer(kind=int32)                        :: nEz      !number of points at which Ez is computed
    integer(kind=int32)                        :: nHxy     !number of points at which both Hx and Hy are computed
    integer(kind=int32)                        :: nHx      !number of points at which only Hx is computed
    integer(kind=int32)                        :: nHy      !number of points at which only Hy is computed
    integer(kind=int32)                        :: nHz      !number of points at which Hz is computed
    logical   :: allcomp_samecoord   !flag indicating if coordintes are equal for different field components
    logical   :: allzrec_samecoord   !flag indicating if horizontal coordinates are equal for different depths (true for grid)

    !background model description
    real(kind=real64),dimension(:),pointer   :: sigh    ! background conductivities for each layer, isotropic or horizontal
    real(kind=real64),dimension(:),pointer   :: sigv    ! background conductivity vector, isotropic or vertical
    real(kind=real64),dimension(:),pointer   :: epsrh   ! background permittivity vector, isotropic or horizontal
    real(kind=real64),dimension(:),pointer   :: epsrv   ! background permittivity vector, isotropic or vertical
    integer(kind=int32)                      :: nlay    ! number of layers
    real(kind=real64),dimension(:),pointer   :: zbound  !depths of layer boundaries in m, positive downward
    real(kind=real64)                        :: omega   !angular frequency
    integer(kind=int32)                      :: aniso   !anisotropy index
    real(kind=real64)                        :: rsplmin !minimum radius at which spline interpolation is allowed
    integer(kind=int32)                      :: dowhat  !what to compute: fields or sensitivities or both
    integer(kind=int32)                      :: infolevel  !how much output info to display
    complex(kind=real64)                     :: sigmabghom !complex conductivity of homogeneous background
  end type backgrounddata


  !frequency-dependent specifications:
  ! the frequencies to model and source currents
  type freqdata
    integer(kind=int32)                         :: nfreq      ! number of frequency components
    real(kind=real64),dimension(:),pointer      :: omega      ! angular frequencies (2*pi*f, f in Hertz)
  end type freqdata


  !frequency domain wavelet
  type wavelet
    integer(kind=int32)                     :: nf          !number of frequency components in wavelet
    real(kind=real64)                       :: domega      !(angular) frequency sampling
    real(kind=real64),dimension(:),pointer  :: re,im       !real and imag part of wavelet (split for spline interpolation)
    real(kind=real64),dimension(:),pointer  :: spline_re,spline_im !spline derivatives for real and imag part of wavelet
    real(kind=real64),dimension(:),pointer  :: omega               !angular frequency vector, 2*pi*f, f in Hertz
  end type wavelet


  !source and receiver specifications
  type sorec
    integer(kind=int32)                         :: type   !source type for each source: dipoles, wire, star
    character(len=namlen)                       :: srcname   !name of source, or dummy name of receiver groups
    integer(kind=int32)                         :: nwire  !for wire sources: number of wires composing a source
    integer(kind=int32)                         :: ncur   !number of input currents for "star" sources
    integer(kind=int32),dimension(:),pointer    :: nelem  !number of source/receiver points within each "shot":
                                                                  ! for dipole sources: number of dipole elements
                                                                  ! for wire sources: number of wire elements on each wire
                                                                  ! for receivers: number of receivers
    character(len=namlen),dimension(:),allocatable :: recnames    !names of individual receivers
    real(kind=real64),dimension(:,:),pointer    :: pos    ! dipole source/receiver positions, will be first dim: xyz, 2dn dim: source index
                                                          !  presently first dim: source index, 2dn dim: xyz!!!
    integer(kind=int32),dimension(:,:),pointer  :: ipos   ! positions in grid units, relative to field points
    integer(kind=int32),dimension(:,:),pointer  :: iposstag ! positions in grid units on staggered grids, rel. to field points
    real(kind=real64),dimension(:),pointer      :: ljx    ! x-length of electric dipole sources (former "sor_jx")
    real(kind=real64),dimension(:),pointer      :: ljy    ! y-length of electric dipole sources
    real(kind=real64),dimension(:),pointer      :: ljz    ! z-length of electric dipole sources
    real(kind=real64),dimension(:),pointer      :: akx    ! x-area of magnetic loop sources (projection of loop into yz-plane)
    real(kind=real64),dimension(:),pointer      :: aky    ! y-area of magnetic loop sources (projection of loop into xz-plane)
    real(kind=real64),dimension(:),pointer      :: akz    ! z-area of magnetic loop sources (projection of loop into xy-plane)

    type(wirespec),dimension(:),pointer         :: wire     !specification of wire sources
    real(kind=real64),dimension(:),pointer      :: omegarot !source current rotation frequencies in "star" sources
    real(kind=real64),dimension(:),pointer      :: fi0      !constant phase angle (degrees) for source current in "star" sources
    character(len=namlen),dimension(:),pointer  :: wavnames !names of wavelet files

    logical                                     :: elsrc    ! electric sources input?
    logical                                     :: magsrc   ! magnetic sources input?

    complex(kind=real32),dimension(:,:),pointer :: cur      ! source currents: source elements x frequencies

    !extra parameters for 2.5D modeling - can be ignored for 1D and 3D
    !these parameters have the format of model grids, but are unique for each source --> keep them in sorec structure!
    real(kind=real64)                           :: shiftx   !shift of x coord to move center of source to x=0
    integer(kind=int32)                         :: ndif     !nr of grid nodes where 2D model differs from background
    integer(kind=int32)                         :: ndifhi   !same as ndif, but excluding overlapping lower edges of my domain
    integer(kind=int32),dimension(:),pointer    :: iydif,izdif  !indices of nodes where 2D model differs from background
    integer(kind=int32),dimension(:),pointer    :: idifhi   !indices of nodes within idif not on lower overlapping edge of my domain
    integer(kind=int32)                         :: idif     !counter, used when computing 1D bg fields
    logical,dimension(:),pointer  :: isstagy  !indicates if Ey is staggered - false at higher edges of anomalies for 1D background
    logical,dimension(:),pointer  :: isstagz  !indicates if Ez is staggered - false at bottom edges of anomalies for 1D background
    logical,dimension(:,:),allocatable          :: needfft  !index for where back transform needs to be done

    !background model for each source is stored individually
    real(kind=real64)  :: sigbghom   ! homogeneous background conductivity
    real(kind=real64)  :: epsbghom   ! homogeneous background permittivity
    !!! 1D has to be added too

  end type sorec


  !definition of long wire source
  type wirespec
    real(kind=real64),dimension(3,2)            :: endpos    !wire start and end points for 1 wire only
    real(kind=real64),dimension(:,:),pointer    :: elempos   !positions of wire elements
    real(kind=real64)                           :: dlw       !length of wire elements for 1 wire only
    integer(kind=int32)                         :: nelem     !number of wire elements
  end type wirespec


  !source types
  integer(kind=int32),parameter   :: receiver = 0
  integer(kind=int32),parameter   :: dipole = 1
  integer(kind=int32),parameter   :: wire = 2
  integer(kind=int32),parameter   :: star = 3

  !source type descriptions
  character(len=10),dimension(receiver:star),parameter  :: srctypestr = (/'receiver','dipole  ','wire    ','star    '/)

  !domains for wavelet definition
  integer(kind=int32),parameter   :: timedom = 1
  integer(kind=int32),parameter   :: freqdom = 2

  integer(kind=int32),parameter,public    :: fwdmodel = 1  !job: forward computation of EM fields
  integer(kind=int32),parameter,public    :: deriv = 2     !job: compute derivatives dEx/deps_m
  integer(kind=int32),parameter,public    :: fwd_deriv = 3 !job: forward computation and derivatives dEx/deps_m

  integer(kind=int32),parameter   :: iso = 0 !flag for isotropic medium
  integer(kind=int32),parameter   :: vti = 1 !flag for VTI-anisotropic medium
  integer(kind=int32),parameter   :: orth = 2    ! orthorhombic: can have 3 different diagonal components of conductivity tensor

  !options for amount of output
  integer(kind=int32),parameter   :: output_final=0 !output final results only
  integer(kind=int32),parameter   :: output_more=1  !output final and intermediate results
  integer(kind=int32),parameter   :: output_all=2   !output final and intermediate results and linear system vectors
  integer(kind=int32),parameter   :: withcond=3 !output everything AND compute cond. number of preprocessed matrix - EXPENSIVE!

contains

#include "AvailableUnit.f90"
#include "readsorec.f90"
include 'get_outputlength.f90'
#ifndef NOFFTW
#include "readsorcur.f90"
#include "readwavelet.f90"
#endif
include 'getfreqs.f90'
#include "fileoperations.f90"

endmodule io_common_mod

