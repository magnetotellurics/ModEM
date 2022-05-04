!==============================================================================!
!================================================================ Dipole1D.f90 !
!==============================================================================!
!
!    Copyright 2007-2010
!    Kerry Key
!    Scripps Institution of Oceanography
!    kkey@ucsd.edu
!
!    This file is part of Dipole1D.
!
!    Dipole1D is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Dipole1D is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Dipole1D.  If not, see <http://www.gnu.org/licenses/>.
!
!
! Dipole1D is a subroutine for computing electric and magnetic fields in a layered 1D
! model from an arbitrarily oriented electric dipole source.  Any number of
! layers can be used.  The source and receivers can be located anywhere 
! in the stack of layers, and do not need to reside in the same layer.  Fields
! can be output in the spatial (x,y,z) domain or wavenumber (kx,y,z) domain.
! Fields in the wavenumber domain can only be computed using dipoles aligned 
! with any of x,y,z.
!
! A Lorenz gauged vector potential formulation similar to that in Wait (1982)
! is used. However, I have derived recursions using decaying exponentials rather
! than using tanh functions.
!
! An exp(-i*omega*t) time dependence is used by default but can be changed to exp(+i*omega*t)
! by setting the phaseConvention parameter.
!
! Revision Notes:
!
! Version 7.3  February 10, 2010.   Added Gauss quadrature option for finite dipole integrations.
! Version 7.2  Feb 5, 2010.         David Myer added support for finite dipole
! Version 7.1  November 5, 2009.    Fixed some bugs in PreComputePotCoeffs that affected the curvature terms
!                                   of the spline interpolation coefficients. Thanks to James Gunning for identifying this.
!                                   Also, I finally modified the interp coeff arrays so that only layers with receivers are stored, 
!                                   resulting in a huge savings in unused memory.  The Canonical_RealImag_BxEyEz example 
!                                   inversion now uses only 5MB instead of 85MB.
! Version 7.0  October 20, 2009.    Rearranged the code to encapsulate subroutines into a Fortran module.
! Version 6.7  May 8, 2009.         A few minor changes to reduce compiler temporary array allocations.
! Version 6.6  December 22, 2008.   Fixed a bug in phase lead code. Changed a few intrinsics to use F2003 
!                                   standard.  A few other minor changes.                                   
! Version 6.5  November 22, 2008.   Added phaseConvention parameter so that 
!                                   the phase lag or lead convention can be specified.
! Version 6.4  September 17, 2008.  Fixed an allocation error in PreComputePotCoeffs that sometimes resulted 
!                                   in lots of extra spline interpolations and wasted CPU time.
! Version 6.3  June 17, 2008.       Fixed sign error in jz_ved_kx kernel function.
! Version 6.2  May, 20, 2008.       Rearranged computations and array shapes for derivative comp speedup. 
!                                   Added spline interpolation for derivatives. Added test-drive function to
!                                   determine range of lambda for interpolation.
! Version 6.1  April 18, 2008.      Reformatted source for compatibility with g95 and gfortran free compilers. 
! Version 6.0  March 14, 2008.      Added model response derivatives with respect to conductivity of
!                                   each layer, to use for Jacobians for 1D inversion.  Currently
!                                   these do not support the spline interpolation speedup. 
! Version 5.3  February 29, 2008.   Added optional spline interpolation for major 
!                                   speedup when modeling many layers, many receivers.
! Version 5.2  February 26, 2008.   Restructured filter modules, added new HT and CT filters                              
! Version 5.1  December 14, 2007.   Modified all routines for exp(-iwt) time dependence
! Version 5.0  June 22, 2007.       Option to output fields in (kx,y,z) domain for 2.5D modeling.
!                                   Uses my own sine/cosine transform digital filters created using
!                                   a similar technique to F.N. Kong, Geophys. Prospect., 2007.
! Version 4.0  June 14, 2007.       N-layer version for arbitrary dipole orientation.
! Version 3.0  June 12, 2007.       N-layer version for horizontal dipole.
! Version 2.0  June 1, 2007.        3 layer version for horizontal dipole. 
! Version 1.0  May 25, 2007.        2 layer version for horizontal dipole.  
!
!------------------------------------------------------------------------------!
! Reference:
!------------------------------------------------------------------------------!
!
! Please reference this paper if you publish results using this code, or if 
! you include Dipole1D as a 1D kernel in your own codes (subject to the 
! licensing terms of the GNU General Public License):
!
!    Key, K., 2009, One-dimensional inversion of multi-component, 
!    multi-frequency marine CSEM data: Methodology and synthetic studies for 
!    resolving thin resistive layers: Geophysics, Vol.74,No.2, March-April 2009,
!    doi:  10.1190/1.3058434 
!
!
!------------------------------------------------------------------------------!
! Acknowledgments:
!------------------------------------------------------------------------------!
!
! This work was supported by:
!
!    The Seafloor Electromagnetic Methods Consortium   
!       at Scripps Institution of Oceanography  
!    
! See this URL for a list of SEMC funding sponsors:
!
!    http://marineemlab.ucsd.edu/semc.html
!     
! John Burkardt is thanked for making subroutine LEGENDRE_COMPUTE_DR
! freely available under the LGPL.  
!
!------------------------------------------------------------------------------!
! Using Dipole1D in your own codes:
!------------------------------------------------------------------------------!
!
! Dipole1D is a fortran module with a main subroutine that can be called from 
! external routines.   All variable passing to/from Dipole1D is done in the 
! module dipole1d, located just below here.  You need to initialize all the public
! variables listed there and can then call subroutine comp_dipole1d to obtain
! the electric and magnetic fields. See example in file Call_Dipole1D.f90.
!
!==============================================================================!
!==================================================================== dipole1d !
!==============================================================================!     

module dipole1d

    implicit none
    
    
!---------------------------------------------------------------------
! Public data:
!---------------------------------------------------------------------

!
! These can be accessed outside Dipole1D by placing a "use dipole1d" statement in your subroutine.
!
!
! Transmitter parameters
!
    real(8), public    :: sdm1D             ! (Am)      Source dipole moment 
    real(8), public    :: ftx1D             ! (Hz)      Source frequency 
    real(8), public    :: azimuthTx1D       ! (degrees) Source azimuth from x axis (positive clockwise)
    real(8), public    :: dipTx1D           ! (degrees) Vertical dip angle of source along azimuthTx1D, positive down 
    real(8), public    :: xTx1D,yTx1D,zTx1D ! (m)       Transmitter x, y and z location   
    real(8), public    :: lenTx1D           ! (m)       Dipole length centered on x,y,z above & having azimuth & dip
    integer, public    :: numIntegPts       !           Number of points to use for Gauss quadrature integration   


!
! Model parameters
!
    integer, public                             :: nlay1D    ! Number of layers
    real(8), dimension(:), allocatable, public  :: sig1D     ! (S/m) Layer conductivities 
    real(8), dimension(:), allocatable, public  :: zlay1D    ! (m)   Depth to top of each layer, first layer ignored 
!
! Site locations (where E and B are computed)
!
    integer, public                             :: n1D          ! number of receiver sites
    real(8), dimension(:), allocatable, public  :: x1D, y1D, z1D   
!
! Option to only compute electric fields.  lbcomp = .false. or .true. 
!
    logical, public :: lbcomp  
!
! Fast Hankel Transform Filter method.  Choose the filter coefficients:
!
    character(32), public    :: HTmethod1D   ! Filter options listed in FilterModules.f90

!
! Output domain (either spatial (x,y,z) or (kx,y,z) wavenumber domain for 2.5D modeling primary fields)
!
    character(32), public    :: outputdomain1D     ! Use 'spatial' or 'kx'
    

!
! If using outputdomain1D='kx', then set kx1d, kxmode1D and CTmethod1D:
!
! kxmode1D is used to define the principle axis for the Tx
! 
    real(8), public          :: kx1D     ! If outputdomain1D='kx', then the fields are 
                                         ! computed at at this wavenumber.  
                              
    integer, public          :: kxmode1D ! 1,2,3 for x,y,z dipoles.  Ignored if spatial computation.
!
! (Co)Sine transform using digital filter method. 
!
    character(32), public    :: CTmethod1D  ! Set to a filter method listed in FilterModules.f90

!
! Use Spline interpolation of coefficients for massive speed up when there are many receivers
!
    logical, public          :: lUseSpline1D 

!
! Phase Convention: 'lag' (default) or 'lead'
! Phase Lag  uses the exp(-i w t) convention and phase increases positively with Tx-Rx range. 
! Phase Lead uses the exp(+i w t) convention and phase becomes more negative with Tx-Rx range.
! 
    character(4), public :: phaseConvention
!
!  Electric and magnetic field vectors.  Note that I output vertical current:
!
    complex(8), dimension(:), allocatable, public :: ex1D,ey1D,jz1D,bx1D,by1D,bz1D  ! units are V/m and T

!
! Inversion sensitivity computations:
!

    logical, public :: linversion  ! set =.true. to output the field derivatives with respect to conductivity
!
!  Electric and magnetic field derivatives with respect to layer conductivities
!  This are computed when linversion = .true.
! 
    complex(8),dimension(:,:), allocatable, public :: dexdsig,deydsig,djzdsig,dbxdsig,dbydsig,dbzdsig !dexdsig(i,j) = df_i / dsig_j   
    
!---------------------------------------------------------------------
! Public subroutines callable from outside the module:
!---------------------------------------------------------------------

    public ::  init_defaults_Dipole1D    ! initializes default values for various parameters
    public ::  comp_Dipole1D             ! The routine you call to perform the 1D computations
    
        
!---------------------------------------------------------------------
! Private data:
! 
! You can't use anything listed below here in external routines.
!---------------------------------------------------------------------  

!
! Constants:
!
    complex(8), parameter, private :: ii = (0d0,1d0)
    real(8), parameter, private    :: pi = 3.141592653589793d0
    real(8), parameter, private    :: mu0 = 4d-7*pi
    real(8), parameter, private    :: eps = 8.8541878176d-12
    
    integer, private, parameter    :: finite_integ_method = 2 ! 1 = 'box', 2 = 'gauss'.  1 is mostly for testing, use 2 for accuracy
    
!
! Local model parameters:
!
    real(8),    dimension(:), allocatable, private :: h    ! layer thicknesses
    complex(8), dimension(:), allocatable, private :: csig ! Locally uses complex conductivity to 
                                                           ! handle propagation in highly resistive layers.
    real(8), private                               :: depthTx,heightTx ! depth/height of Tx in its layer
    real(8), private                               :: depthRx,heightRx ! depth/height of Rx in its layer
  
!
! Vector potentials Ay and Az (lambda_z) potentials, their derivates 
! and coefficients
!
    complex(8), dimension(:), allocatable, private :: a,b,c,d, expmgh    
    complex(8), dimension(:), allocatable, private :: gamma,k2,Rm,Rp,Sm,Sp 
  
    real(8), private     :: omega    ! used by HT kernels
    complex(8), private  :: ay, az
    complex(8), private  :: daydz, dazdz, dazdz2, aypdazdz, daydzpdazdz2
    
!
! Local quantities used for derivatives with respect to conductivity
!
    complex(8), dimension(:,:), allocatable, private :: dRmdsig, dRpdsig, dSmdsig, dSpdsig
    complex(8), dimension(:,:), allocatable, private :: dadsig, dbdsig, dcdsig, dddsig    ! coeff_layer by sigma_layer
    complex(8), dimension(:), allocatable, private   :: dAydsig, dAzdsig, dAydsigdz, dAzdsigdz, dAzdsigdz2
    complex(8), dimension(:), allocatable, private   :: dae1dsig, dbe2dsig, dce1dsig, dde2dsig
    
!
! Site parameters:
!
    complex(8), private  :: csigsite  ! stores site conductivity for the Hankel transform kernels
    real(8), private     :: r,z,theta ! distance and angle from x axis to site location, used by HT kernels
    real(8), private     :: dx,dy     ! signed distances from site to transmitter
    integer, private     :: iTxlayer  ! source layer
    integer, private     :: iRxlayer  ! receiver layer
    real(8), private     :: rotang, thetaRx ! Angle from x to rotated x axis,  Angle from x to site
    integer, private     :: isgnsrc, isgndsrcdz
    
    real(8), private     :: sintheta, costheta, sin2theta, cos2theta

!
! Fast Hankel Transform filters:
!
    integer, private                    :: ndhtfc
    real(8), dimension(1601), private   :: base, htj0, htj1 !1601 is largest possible filter length
!
! Fast (Co)Sine Transform filters:
!
    integer, private                  :: ncsfc
    real(8), dimension(1601), private :: basecsfc, cosfc, sinfc
    
!
! VED and HED specifiers
!
    logical, private  :: lved, lhed
!
! Mode for (kx,y,z) computations (only used when outputdomain1D=kx is set on input)
!
    integer, private  :: kxmode

!
! Fields passed from comp_hed_fht and comp_ved_fht to comp_spatial
!    
    complex(8), private  :: exh,eyh,jzh,bxh,byh,bzh
    complex(8), private  :: exv,eyv,jzv,bxv,byv,bzv,ext,bxt  
    
    ! includes these if linversion = .true.
    complex(8), dimension(:), allocatable, private :: dexhdsig,deyhdsig,djzhdsig
    complex(8), dimension(:), allocatable, private :: dbxhdsig,dbyhdsig,dbzhdsig
    complex(8), dimension(:), allocatable, private :: dexvdsig,deyvdsig,djzvdsig
    complex(8), dimension(:), allocatable, private :: dbxvdsig,dbyvdsig,dbzvdsig 
    complex(8), dimension(:), allocatable, private :: dextdsig,dbxtdsig ! working arrays for ved 
    
    ! and these are there for passing in/out from kernel functions:
    complex(8), dimension(:), allocatable, private :: dfdsigj0, dfdsigj1  
        
!
! Index points to use for interpolations:
!
    integer, private    :: nlam_interp   ! number of interpolation index points in lambda
    integer, private    :: iptr          ! a pointer to the last interpolation indices, used for quickly finding next indices
    real(8), dimension(:), allocatable, private :: lam_interp      ! interpolation index points 
    integer, dimension(:), allocatable, private :: iRxLayerInterp  ! Mapping from layer i to interpolation array column, so that
                                                                   ! we only store interp coeffs in layers with receivers...can be
                                                                   ! a huge savings in memory for models with 1000's layers
!
! Potential coefficients and second derivatives at index points for interpolation
! 
!  These are nlayers by nlam_interp in size
!
  ! Potential coefficients for interpolation:
    complex(8), dimension(:,:), allocatable, private :: a_hed_interp, b_hed_interp, c_hed_interp, d_hed_interp  
    complex(8), dimension(:,:), allocatable, private :: a2_hed_interp,b2_hed_interp,c2_hed_interp,d2_hed_interp 
    complex(8), dimension(:,:), allocatable, private :: c_ved_interp, d_ved_interp    
    complex(8), dimension(:,:), allocatable, private :: c2_ved_interp,d2_ved_interp 
    !
    ! Terms for inversion sensitivity:
    ! These are nlam_interp by nlayers by nlayers in size, where the last index is the sensitivity layer    
    complex(8), dimension(:,:,:), allocatable, private :: dadsig_hed_interp, dbdsig_hed_interp
    complex(8), dimension(:,:,:), allocatable, private :: dcdsig_hed_interp, dddsig_hed_interp 
    complex(8), dimension(:,:,:), allocatable, private :: dadsig2_hed_interp,dbdsig2_hed_interp
    complex(8), dimension(:,:,:), allocatable, private :: dcdsig2_hed_interp,dddsig2_hed_interp 
    complex(8), dimension(:,:,:), allocatable, private :: dcdsig_ved_interp, dddsig_ved_interp   
    complex(8), dimension(:,:,:), allocatable, private :: dcdsig2_ved_interp,dddsig2_ved_interp  
    
!--------------------------------------------------------
! Private subroutines only called from within the module:
!--------------------------------------------------------
    private :: comp_spatial, comp_kx
    private :: comp_kx_x, comp_kx_y, comp_kx_z
    private :: comp_hed_fht, comp_ved_fht, rotate_fields
    private :: allocate_dipole1d, deallocate_dipole1d
    private :: PreComputePotCoeffs, deallocate_spline
    private :: initialize_dipole1d, initialize_Tx, setupsite
    private :: ex_hed_x_kx, ey_hed_x_kx, jz_hed_x_kx
    private :: bx_hed_x_kx, by_hed_x_kx, bz_hed_x_kx
    private :: ex_hed_y_kx, ey_hed_y_kx, jz_hed_y_kx    
    private :: bx_hed_y_kx, by_hed_y_kx, bz_hed_y_kx
    private :: ex_ved_kx, ey_ved_kx, jz_ved_kx
    private :: bx_ved_kx, by_ved_kx, bz_ved_kx
    private :: ex_ved_j1, dexdsig_ved_j1
    private :: bx_ved_j1, dbxdsig_ved_j1
    private :: jz_ved_j0, djzdsig_ved_j0
    private :: ex_hed_j0, dexdsig_hed_j0, ex_hed_j1, dexdsig_hed_j1
    private :: ey_hed_j0, deydsig_hed_j0, ey_hed_j1, deydsig_hed_j1
    private :: jz_hed_j1, djzdsig_hed_j1
    private :: bx_hed_j0, dbxdsig_hed_j0
    private :: bx_hed_j1, dbxdsig_hed_j1
    private :: by_hed_j0, dbydsig_hed_j0
    private :: by_hed_j1, dbydsig_hed_j1
    private :: bz_hed_j1, dbzdsig_hed_j1
    private :: hed_coef_nlayer, hed_pot_nlayer, ved_coef_nlayer, ved_pot_nlayer
    private :: spline_i1d, splint_ved, splint_hed
!
!------------------------
! Contained subroutines:
!------------------------
! Thes variables above have global scope within the contained subroutines below.
!
    contains
 
!==============================================================================!
!====================================================== init_defaults_Dipole1D !
!==============================================================================!    
    subroutine init_defaults_Dipole1D 
    
!
! Specify some defaults for parameters required by Dipole1D:
!
    HTmethod1D      = 'kk_ht_201'    ! Use 201 point HT digital filters.
    outputdomain1D  = 'spatial'      ! Assume spatial domain comps
    lbcomp          = .false.        ! This is changed to true if magnetics in data file
    sdm1D           = 1.0            ! (Am), dipole moment. Normalize to unit source moment
    lUseSpline1D    = .true.         ! Use spline interpolation for faster 1D computations
    linversion      = .false.        ! Compute derivatives with respect to sigma(layers)
    phaseConvention = 'lag'          ! The usual default is lag, where phase becomes larger positive values with increasing range.
    lenTx1D         = 0.d0           ! (m) Dipole length 0 = point dipole
    numIntegPts     = 10             ! Number of points to use for Gauss quadrature integration for finite dipole
 
    end subroutine init_defaults_Dipole1D 
    
   
    
!==============================================================================!
!==================================================================== DIPOLE1D !
!==============================================================================!
    subroutine comp_dipole1d
!
! This is the main subroutine that you need to call after initializing all
! the public variables at listed at the top of this module.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
 
    implicit none       
!
! local variables used for finite dipole calculations
!
    real(8) :: finiteStepSize  ! (m) distance between point dipoles summed when lenTx1D /= 0  
    real(8) :: dip, azm        ! converted to radians for convenience
    real(8) :: ctr_xTx1D, ctr_yTx1D, ctr_zTx1D 
    
    real(8), dimension(:), allocatable      ::  weights, xint  ! weights used for integration and local position along dipole
    complex(8), dimension(:), allocatable   :: w_ex1D, w_ey1D, w_jz1D, w_bx1D, w_by1D, w_bz1D
    complex(8), dimension(:,:), allocatable :: w_dexdsig, w_deydsig, w_djzdsig, w_dbxdsig, w_dbydsig, w_dbzdsig
    integer :: nActual1D, nActualnlay1D

!
! local variables for loops, etc:
!
    integer    :: i,j

!
! Allocate local arrays
!
    call allocate_dipole1d 
! 
! Initialize local arrays, check for correct input model, setup digital filters:
! Sets up transmitter specific variables
!   
    call initialize_dipole1d
    
!
! DGM Feb 2010 implement finite length dipole by using linear superposition
!       on a series of point dipoles spaced finiteStepSize meters apart along 
!       the dipole length. 
!         Note that IF the dip is 0, this could be made faster with reciprocity.
!       But for now let's just get it working, shall we?  Dip is usually
!       non-zero in real data anyway.
!   
! KWK Feb 10, 2010  Added Gauss quadrature option to David's finite dipole integration loop
!
    if (lenTx1D /= 0.d0) then
           
        ! For convenience, get dip & azm in radians
        dip = dipTx1D * pi / 180.d0
        azm = azimuthTx1D * pi / 180.d0
        
        ! Save center position of the dipole
        ctr_xTx1D   = xTx1D
        ctr_yTx1D   = yTx1D
        ctr_zTx1D   = zTx1D
        
        !
        ! Compute weights and integration points:
        !
        if (numIntegPts < 1) numIntegPts = 1  ! stupid user, don't do that!
        
        allocate(weights(numIntegPts), xint(numIntegPts))
        
        select case (finite_integ_method)
        
        case (1)  ! brute force mate!
        
            finiteStepSize = lenTx1D / numIntegPts 
            weights = finiteStepSize
            do i = 1,numIntegPts
                xint(i) =  -1.d0 * lenTx1D / 2.d0  + finiteStepSize/2.d0 + (i-1)*finiteStepSize
            enddo
            
        case (2)  ! Gauss quadrature
        
            call legendre_compute_dr (numIntegPts,xint,weights)  ! returns weights on -1,1 interval
            xint = xint*lenTx1D / 2.d0                           ! move from -1,1 to  -1/2 lenTx1D to +1/2 lenTx1D
            weights = weights * lenTx1D / 2.d0  
            
        end select 
 
        
        ! NOTE: Do NOT use n1D and nlay1D to allocate the size of the working
        !       sum variables.  If this is a multi-frequency inversion, then
        !       each freq may not have the same number of sites and computeFwd_CSEM
        !       will change n1D for each call even though ex1D, etc... stay the
        !       maximum size.  Then a seemingly simple thing like ex1D = w_ex1D
        !       will blow chunks (access violation) because these arrays are not
        !       the same size!  It is simpler to make the working vars the same
        !       size than to always remember to code ex1D(1:n1D) = w_ex1D(1:n1D)
        nActual1D = size(ex1D,1)
        if (linversion) then
            nActualnlay1D  = size(dexdsig,2)
        endif
        
        ! Allocate the working variables into which everything will be summed.
        allocate (w_ex1D(nActual1D), w_ey1D(nActual1D), w_jz1D(nActual1D))
        allocate (w_bx1D(nActual1D), w_by1D(nActual1D), w_bz1D(nActual1D))
        w_ex1D  = (0.d0,0.d0)
        w_ey1D  = (0.d0,0.d0)
        w_jz1D  = (0.d0,0.d0)
        w_bx1D  = (0.d0,0.d0)
        w_by1D  = (0.d0,0.d0)
        w_bz1D  = (0.d0,0.d0)
        if (linversion) then
            allocate (w_dexdsig(nActual1D,nActualnlay1D), w_dbxdsig(nActual1D,nActualnlay1D))
            allocate (w_deydsig(nActual1D,nActualnlay1D), w_dbydsig(nActual1D,nActualnlay1D))
            allocate (w_djzdsig(nActual1D,nActualnlay1D), w_dbzdsig(nActual1D,nActualnlay1D))
            w_dexdsig   = (0.d0,0.d0)
            w_deydsig   = (0.d0,0.d0)
            w_djzdsig   = (0.d0,0.d0)
            w_dbxdsig   = (0.d0,0.d0)
            w_dbydsig   = (0.d0,0.d0)
            w_dbzdsig   = (0.d0,0.d0)
        endif
    else
        allocate(xint(1),weights(1))
        xint   = 0.d0
        weights = 1.d0
    endif
    
    !
    ! Loop through xint and sum up weighted response values:
    !
    do j = 1,size(xint,1)
       
        ! Calculate this dipole position based on azimuth, dip, & distance along the dipole
        ! Coord system is NED: x=north, y=east, z=+ve down
        ! Azimuth is degrees clockwise from north of the HEAD
        ! Dip is degrees DOWN from horizontal of the HEAD
        if (lenTx1D /= 0.d0) then
            xTx1D   = ctr_xTx1D + (xint(j) * cos(dip) * cos(azm))
            yTx1D   = ctr_yTx1D + (xint(j) * cos(dip) * sin(azm))
            zTx1D   = ctr_zTx1D + (xint(j) * sin(dip))
        endif    
       
        ! Init Tx position related variables
        call initialize_Tx
        
        !
        !  Pre-compute potential coefficients for later use in spline interpolation:
        !
        if ( n1d==1 ) then ! Don't use spline if only a single receiver 
            lUseSpline1D = .false.
        endif
        if (lUseSpline1D) then     
            call PrecomputePotCoeffs
        endif
!
! Loop over sites and compute fields for hed and ved:
!
        do i = 1,n1D   
!           
! Get site layer and geometric information
!
            call setupsite(i)   
!
! Compute the fields in the spatial or kx domain:
!
            select case (trim(outputdomain1D))
            case ('spatial')
                call comp_spatial(i)
                    
            case  ('kx') 
                call comp_kx(i)
             
            end select !case (trim(outputdomain1D))
    
        enddo ! loop over sites
    
!
! Apply Phase Convention:
!
        if (phaseConvention(1:4) == 'lead') then 
            ! Change from Dipole1D's lag (exp(-iwt)) to lead (exp(+iwt))
            ex1D = conjg(ex1D)
            ey1D = conjg(ey1D)
            jz1D = conjg(jz1D)
            bx1D = conjg(bx1D)
            by1D = conjg(by1D)
            bz1D = conjg(bz1D)
            if  (linversion) then
                dexdsig = conjg(dexdsig)
                deydsig = conjg(deydsig)
                djzdsig = conjg(djzdsig)
                dbxdsig = conjg(dbxdsig)
                dbydsig = conjg(dbydsig)
                dbzdsig = conjg(dbzdsig)
            endif 
        endif
!
! Zero out low values near limits of HT filters:
!
!	   where (abs(ex1d) < 1d-25) ex1d = 0d0
!	   where (abs(ey1d) < 1d-25) ey1d = 0d0
!	   where (abs(jz1d) < 1d-25) jz1d = 0d0
!	   where (abs(bx1d) < 1d-25) bx1d = 0d0
!	   where (abs(by1d) < 1d-25) by1d = 0d0
!	   where (abs(bz1d) < 1d-25) bz1d = 0d0
!	   
!	   if  (linversion) then
!		   where (abs(dexdsig) < 1d-25) dexdsig = 0d0
!		   where (abs(deydsig) < 1d-25) deydsig = 0d0
!		   where (abs(djzdsig) < 1d-25) djzdsig = 0d0
!		   where (abs(dbxdsig) < 1d-25) dbxdsig = 0d0
!		   where (abs(dbydsig) < 1d-25) dbydsig = 0d0
!		   where (abs(dbzdsig) < 1d-25) dbzdsig = 0d0  
!	   endif

        ! If doing a finite dipole, add weighted contribution to the sum
        if (lenTx1D /= 0.d0) then
            w_ex1D = w_ex1D + ex1D*weights(j)
            w_ey1D = w_ey1D + ey1D*weights(j)
            w_jz1D = w_jz1D + jz1D*weights(j)
            w_bx1D = w_bx1D + bx1D*weights(j)
            w_by1D = w_by1D + by1D*weights(j)
            w_bz1D = w_bz1D + bz1D*weights(j)
            if (linversion) then
                w_dexdsig = w_dexdsig + dexdsig*weights(j)
                w_deydsig = w_deydsig + deydsig*weights(j)
                w_djzdsig = w_djzdsig + djzdsig*weights(j)
                w_dbxdsig = w_dbxdsig + dbxdsig*weights(j)
                w_dbydsig = w_dbydsig + dbydsig*weights(j)
                w_dbzdsig = w_dbzdsig + dbzdsig*weights(j)
            endif 
        endif
        
        ! Dump the spline-related variables
        ! Will need to recompute them for each transmitter because layer 
        ! of transmitter may change and max/min range to Rx WILL change.
        call deallocate_spline
        
    enddo  ! j = 1,size(xint,1)       
   
!
! If we did a finite dipole...
!
    if (lenTx1D /= 0.d0) then
        ! Return Tx position to its center value so future calls (for other
        ! frequencies) have the correct center.
        xTx1D   = ctr_xTx1D
        yTx1D   = ctr_yTx1D
        zTx1D   = ctr_zTx1D
        
        ! Scale integral by lenTx1D, then move values to return variables & deallocate the summing variables
        ex1D = w_ex1D / lenTx1D
        ey1D = w_ey1D / lenTx1D
        jz1D = w_jz1D / lenTx1D
        bx1D = w_bx1D / lenTx1D
        by1D = w_by1D / lenTx1D
        bz1D = w_bz1D / lenTx1D
        ! Allocate the working variables into which everything will be summed.
        deallocate (w_ex1D, w_ey1D, w_jz1D, w_bx1D, w_by1D, w_bz1D)
        if (linversion) then
            dexdsig = w_dexdsig /lenTx1D
            deydsig = w_deydsig /lenTx1D
            djzdsig = w_djzdsig /lenTx1D
            dbxdsig = w_dbxdsig /lenTx1D
            dbydsig = w_dbydsig /lenTx1D
            dbzdsig = w_dbzdsig /lenTx1D
            deallocate (w_dexdsig, w_deydsig, w_djzdsig, w_dbxdsig, w_dbydsig, w_dbzdsig)
        endif
    end if
    
    deallocate(xint,weights)
    
!
! Deallocate local arrays
!
    call deallocate_dipole1d    
!
! Goodbye.
!
    end subroutine comp_dipole1d

subroutine legendre_compute_dr ( order, xtab, weight )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_DR computes a Gauss-Legendre rule, Davis-Rabinowitz method.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    ORDER must be greater than 0.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric, and should sum to 2.
!
    implicit none
    
    integer ( kind = 4 ) order
    
    real    ( kind = 8 ) d1
    real    ( kind = 8 ) d2pn
    real    ( kind = 8 ) d3pn
    real    ( kind = 8 ) d4pn
    real    ( kind = 8 ) dp
    real    ( kind = 8 ) dpn
    real    ( kind = 8 ) e1
    real    ( kind = 8 ) fx
    real    ( kind = 8 ) h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) iback
    integer ( kind = 4 ) k
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mp1mi
    integer ( kind = 4 ) ncopy
    integer ( kind = 4 ) nmove
    real    ( kind = 8 ) p
    real    ( kind = 8 ) :: pi = 3.141592653589793D+00
    real    ( kind = 8 ) pk
    real    ( kind = 8 ) pkm1
    real    ( kind = 8 ) pkp1
    real    ( kind = 8 ) t
    real    ( kind = 8 ) u
    real    ( kind = 8 ) v
    real    ( kind = 8 ) x0
    real    ( kind = 8 ) xtab(order)
    real    ( kind = 8 ) xtemp
    real    ( kind = 8 ) weight(order)
    
    if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
    end if
    
    e1 = real ( order * ( order + 1 ), kind = 8 )
    
    m = ( order + 1 ) / 2
    
    do i = 1, m
    
    mp1mi = m + 1 - i
    
    t = real ( 4 * i - 1, kind = 8 ) * pi &
      / real ( 4 * order + 2, kind = 8 )
    
    x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( order, kind = 8 ) ) &
      / real ( 8 * order * order, kind = 8 ) )
    
    pkm1 = 1.0D+00
    pk = x0
    
    do k = 2, order
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do
    
    d1 = real ( order, kind = 8 ) * ( pkm1 - x0 * pk )
    
    dpn = d1 / ( 1.0D+00 - x0 * x0 )
    
    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 * x0 )
    
    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0 * x0 )
    
    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) / &
      ( 1.0D+00 - x0 * x0 )
    
    u = pk / dpn
    v = d2pn / dpn
    !
    !  Initial approximation H:
    !
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn / &
      ( 3.0D+00 * dpn ) ) ) )
    !
    !  Refine H using one step of Newton's method:
    !
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )
    
    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )
    
    h = h - p / dp
    
    xtemp = x0 + h
    
    xtab(mp1mi) = xtemp
    
    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )
    
    weight(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp * xtemp ) / ( fx * fx )
    
    end do
    
    if ( mod ( order, 2 ) == 1 ) then
    xtab(1) = 0.0D+00
    end if
    !
    !  Shift the data up.
    !
    nmove = ( order + 1 ) / 2
    ncopy = order - nmove
    
    do i = 1, nmove
    iback = order + 1 - i
    xtab(iback) = xtab(iback-ncopy)
    weight(iback) = weight(iback-ncopy)
    end do
    !
    !  Reflect values for the negative abscissas.
    !
    do i = 1, order - nmove
    xtab(i) = - xtab(order+1-i)
    weight(i) = weight(order+1-i)
    end do
    
    return
    end subroutine 
!==============================================================================!
!================================================================ comp_spatial !
!==============================================================================!
    subroutine comp_spatial(i)
!
! Subroutine that computes the spatial domain fields for the input 
! model and dipole parameters
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!

    implicit none
    
    integer, intent(in) :: i
    real(8)             :: ang, ss, cc  

!
! Zero out the field quantities:
!
    exh = 0d0; eyh = 0d0; jzh = 0d0
    bxh = 0d0; byh = 0d0; bzh = 0d0
    exv = 0d0; eyv = 0d0; jzv = 0d0
    bxv = 0d0; byv = 0d0; bzv = 0d0
    ext = 0d0; bxt = 0d0
    
    
    if (linversion) then
        dexhdsig = 0d0; deyhdsig = 0d0; djzhdsig = 0d0
        dbxhdsig = 0d0; dbyhdsig = 0d0; dbzhdsig = 0d0
        dexvdsig = 0d0; deyvdsig = 0d0; djzvdsig = 0d0
        dbxvdsig = 0d0; dbyvdsig = 0d0; dbzvdsig = 0d0    
        dextdsig = 0d0; dbxtdsig = 0d0;
    endif
!
! Compute fields for horizontal electric dipole at azimuthTx1D  
!
    if (lhed)   then
        call comp_hed_fht   
    endif
!
! Compute fields for vertical electric dipole:  
!           
    if (lved) then
        call comp_ved_fht 
    endif
!
! Add horizontal and vertical fields using superposition to get correct fields for 
! input dipole azimuthTx1D and dipTx1D
!
    ang = dipTx1D*pi/180.d0
    cc = cos(ang)
    ss = sin(ang)
    ex1D(i) = exh*cc + exv*ss ! save fields for site i in output arrays
    ey1D(i) = eyh*cc + eyv*ss
    jz1D(i) = jzh*cc + jzv*ss
    bx1D(i) = bxh*cc + bxv*ss
    by1D(i) = byh*cc + byv*ss
    bz1D(i) = bzh*cc + bzv*ss   
    
    if (linversion) then
    
        dexdsig(i,:) = dexhdsig*cc + dexvdsig*ss  ! arrays are n1D (#rxs) by nlay1D
        deydsig(i,:) = deyhdsig*cc + deyvdsig*ss
        djzdsig(i,:) = djzhdsig*cc + djzvdsig*ss
        dbxdsig(i,:) = dbxhdsig*cc + dbxvdsig*ss
        dbydsig(i,:) = dbyhdsig*cc + dbyvdsig*ss
        dbzdsig(i,:) = dbzhdsig*cc + dbzvdsig*ss  
    
    endif
    
    end subroutine comp_spatial

!==============================================================================!
!===================================================================== comp_kx !
!==============================================================================!
    subroutine comp_kx(i)
!
! Subroutine that computes the (kx,y,z) domain fields for the input 
! model and dipole parameters at site i
!
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    
    implicit none
    
    integer, intent(in) :: i 
!
! Compute the fields in the (kx,y,z) domain for the particular source direction:
!
     
!
! HED pointing along x
!
    if (kxmode.eq.1) then  
        call comp_kx_x(kx1D,ex1D(i),ey1D(i),jz1D(i),bx1D(i),by1D(i),bz1D(i),lbcomp)
!
! HED pointing along y
!       
    elseif (kxmode.eq.2) then 
        call comp_kx_y(kx1D,ex1D(i),ey1D(i),jz1D(i),bx1D(i),by1D(i),bz1D(i),lbcomp)
!
! VED pointing along z
!   
    elseif (kxmode.eq.3) then 
        call comp_kx_z(kx1D,ex1D(i),ey1D(i),jz1D(i),bx1D(i),by1D(i),bz1D(i),lbcomp)
    else
        write(*,*) 'Error in Dipole1D, kxmode is unknown! kxmode: ',kxmode
        write(*,*) 'Stopping.'
        stop
    endif

    end subroutine comp_kx
    
!==============================================================================!
!=================================================================== comp_kx_x !
!==============================================================================!
    subroutine comp_kx_x(kx,ex,ey,jz,bx,by,bz,lbcomp)   
!
! Compute the (kx,y,z) domain fields for a dipole oriented along the x-axis
!
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    use SinCosFilters
  
    implicit none
    
    complex(8), intent(out) :: ex,ey,jz,bx,by,bz  ! Field components for current site
    real(8), intent(in)     :: kx
    logical,intent(in)      :: lbcomp   
!
! Local variables:
!
    integer    :: j, ysign
    real(8)    :: lam, ky, ruse
    complex(8) :: f  
    
    complex(8), parameter :: coeffcos = 1d0/pi
    complex(8), parameter :: coeffsin = II/pi 

!
! Initialize sums to zero:
!
    ex = 0d0;    ey = 0d0;    jz = 0d0
    bx = 0d0;    by = 0d0;    bz = 0d0

    iptr = 1  ! reset for spline interpolation
    
!
! Check for |dy| < 1
!
   if (dy.lt.0) then
      ysign = -1     !  we only need to use the sign term for odd integrals
   elseif (dy.gt.0) then
      ysign = 1
   elseif (dy.eq.0) then
      ysign = 0 
   endif
   ruse = maxval((/abs(dy),1d0/)) ! always use positive range, adjust sign afterwards on odd integrals
!
! Loop over digital Fourier transform filter coefficients:
!       
    do j = 1,ncsfc

        ky = basecsfc(j)/ruse  
        
        lam = sqrt( kx**2 + ky**2 )     ! kx is input, ky is determined from filter base
!
! N-layer HED potential:
!
        call hed_pot_nlayer(lam)
!
! Compute the field kernels and take product with filter:
!
! For an x directed dipole field symmetries are:
! Even in ky: Ex,By,Ez
! Odd in ky:  Bx,Ey,Bz
!
        ! Ex: 
        call ex_hed_x_kx(kx,ky,f)
        ex = ex + f*cosfc(j)
                
        ! Ey:   
        call ey_hed_x_kx(kx,ky,f)
        ey = ey + f*sinfc(j)

        ! Jz:   
        call jz_hed_x_kx(kx,ky,f)
        jz = jz + f*cosfc(j)
        
        if (lbcomp) then
            ! Bx: 
            call bx_hed_x_kx(kx,ky,f)
            bx = bx + f*sinfc(j)
                    
            ! By:   
            call by_hed_x_kx(kx,ky,f)
            by = by + f*cosfc(j)
    
            ! Bz:   
            call bz_hed_x_kx(kx,ky,f)
            bz = bz + f*sinfc(j)
        endif

    enddo  ! loop over filter coefficients  
    
!
! Apply coefficients and renormalize:
!
    ex = coeffcos*ex/ruse
    ey = coeffsin*ey/ruse*ysign
    jz = coeffcos*jz/ruse
    bx = coeffsin*bx/ruse*ysign
    by = coeffcos*by/ruse 
    bz = coeffsin*bz/ruse*ysign
    
    
    end subroutine comp_kx_x
    
!==============================================================================!
!=================================================================== comp_kx_y !
!==============================================================================!
    subroutine comp_kx_y(kx,ex,ey,jz,bx,by,bz,lbcomp)   
!
! Compute the (kx,y,z) domain fields for a dipole oriented along the y-axis
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    use SinCosFilters
    
    implicit none
    
    complex(8), intent(out) :: ex,ey,jz,bx,by,bz  ! Field components for current site
    real(8), intent(in)     :: kx
    logical,intent(in)      :: lbcomp   
!
! Local variables:
!
    integer    :: j, ysign
    real(8)    :: lam, ky, ruse
    complex(8) :: f

    complex(8), parameter :: coeffcos = 1d0/pi
    complex(8), parameter :: coeffsin = II/pi 

!
! Initialize sums to zero:
!
    ex = 0d0;    ey = 0d0;    jz = 0d0
    bx = 0d0;    by = 0d0;    bz = 0d0
    
    iptr = 1  ! reset for spline interpolation
    
!
! Check for |dy| < 1
!
   if (dy.lt.0) then
      ysign = -1     !  we only need to use the sign term for odd integrals
   elseif (dy.gt.0) then
      ysign = 1
   elseif (dy.eq.0) then
      ysign = 0
   endif
   ruse = maxval((/abs(dy),1.d0/)) ! always use positive range, adjust sign afterwards on odd integrals

!
! Loop over digital Fourier transform filter coefficients:
!       
    do j = 1,ncsfc
    
        ky = basecsfc(j)/ruse  
        
        lam = sqrt( kx**2 + ky**2 )     ! kx is input, ky is determined from filter base
!
! N-layer HED potential:
!
        call hed_pot_nlayer(lam)
!
! Compute the field kernels and take product with filter:
!
! For a y directed dipole, field symmetries are:
! Even in ky: Bx,Ey,Bz
! Odd in ky:  Ex,By,Ez
!
        ! Ex: 
        call ex_hed_y_kx(kx,ky,f)
        ex = ex + f*sinfc(j)
                
        ! Ey:   
        call ey_hed_y_kx(kx,ky,f)
        ey = ey + f*cosfc(j)
        
        ! Jz:   
        call jz_hed_y_kx(kx,ky,f)
        jz = jz + f*sinfc(j)
        
        if (lbcomp) then
            ! Bx: 
            call bx_hed_y_kx(kx,ky,f)
            bx = bx + f*cosfc(j)
                    
            ! By:   
            call by_hed_y_kx(kx,ky,f)
            by = by + f*sinfc(j)
    
            ! Bz:   
            call bz_hed_y_kx(kx,ky,f)
            bz = bz + f*cosfc(j)
        endif

    enddo  ! loop over filter coefficients  
    
!
! Apply coefficients and renormalize:
!
    ex = coeffsin*ex/ruse*ysign
    ey = coeffcos*ey/ruse
    jz = coeffsin*jz/ruse*ysign
    bx = coeffcos*bx/ruse
    by = coeffsin*by/ruse*ysign 
    bz = coeffcos*bz/ruse
    
   
    end subroutine comp_kx_y
!==============================================================================!
!=================================================================== comp_kx_z !
!==============================================================================!
    subroutine comp_kx_z(kx,ex,ey,jz,bx,by,bz,lbcomp)   
!
! Compute the (kx,y,z) domain fields for a dipole oriented along the z-axis
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    use SinCosFilters

    
    implicit none
    
    complex(8), intent(out) :: ex,ey,jz,bx,by,bz  ! Field components for current site
    real(8), intent(in)     :: kx
    logical,intent(in)      :: lbcomp   
!
! Local variables:
!
    integer    :: j, ysign
    real(8)    :: lam, ky, ruse
    complex(8) :: f
    
    complex(8), parameter :: coeffcos = 1d0/pi
    complex(8), parameter :: coeffsin = II/pi

!
! Initialize fields:
!
    ex = 0d0;    ey = 0d0;    jz = 0d0
    bx = 0d0;    by = 0d0;    bz = 0d0
  
    iptr = 1  ! reset for spline interpolation
    
!
! Check for |dy| < 1
!
   if (dy.lt.0) then
      ysign = -1     !  we only need to use the sign term for odd integrals
   elseif (dy.gt.0) then
      ysign = 1
   elseif (dy.eq.0) then
      ysign = 0
   endif   
   ruse = maxval((/abs(dy),1d0/)) ! always use positive range, adjust sign afterwards on odd integrals

!
! Loop over digital Fourier transform filter coefficients:
!       
    do j = 1,ncsfc

        ky = basecsfc(j)/ruse  
        
        lam = sqrt( kx**2 + ky**2 )     ! kx is input, ky is determined from filter base
!
! N-layer HED potential:
!
        call ved_pot_nlayer(lam)
!
! Compute the field kernels and take product with filter:
!
! For an z directed dipole field symmetries are:
! Even in ky: Ex,By,Ez
! Odd in ky:  Bx,Ey,Bz=0
!
        ! Ex: 
        call ex_ved_kx(kx,ky,f)
        ex = ex + f*cosfc(j)
                
        ! Ey:   
        call ey_ved_kx(kx,ky,f)
        ey = ey + f*sinfc(j)

        ! Jz:   
        call jz_ved_kx(kx,ky,f)
        jz = jz + f*cosfc(j)

        if (lbcomp) then
            ! Bx: 
            call bx_ved_kx(kx,ky,f)
            bx = bx + f*sinfc(j)
                    
            ! By:   
            call by_ved_kx(kx,ky,f)
            by = by + f*cosfc(j)
    
            ! Bz:   is zero for ved
        endif
        

    enddo  ! loop over filter coefficients  
    
!
! Apply coefficients and renormalize:
!
    ex = coeffcos*ex/ruse
    ey = coeffsin*ey/ruse*ysign
    jz = coeffcos*jz/ruse
    bx = coeffsin*bx/ruse*ysign
    by = coeffcos*by/ruse 
    bz = 0d0  ! always for vertical Tx in 1D layers
 
 
    end subroutine comp_kx_z
    
!==============================================================================!
!================================================================ comp_hed_fht !
!==============================================================================!
    subroutine comp_hed_fht
!
! Compute the hed fields at the current site using a fast Hankel transform
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    
    implicit none 
    
!
! Local variables:
!
    integer    :: i,j   
    real(8)    :: lam
    complex(8) :: exj0, exj1, eyj0, eyj1, jzj1 
    complex(8) :: bxj0, bxj1, byj0, byj1, bzj1 

    iptr = 1  ! reset for spline interpolation
    
!
! Loop over digital Hankel transform filter coefficients:
!       
    do j = 1,ndhtfc
        
        lam = base(j)/r     
         
      ! N-layer HED potential:

        call hed_pot_nlayer(lam)
        
        ! Ex: 
        call ex_hed_j0(lam,exj0)
        call ex_hed_j1(lam,exj1)   
        exh = exh + exj0*htj0(j)  + exj1*htj1(j)  
                
        ! Ey:   
        call ey_hed_j0(lam,eyj0)
        call ey_hed_j1(lam,eyj1)         
        eyh = eyh + eyj0*htj0(j)  + eyj1*htj1(j)

        ! Jz
        call jz_hed_j1(lam,jzj1)   
        jzh = jzh + jzj1*htj1(j)      
             
        if (lbcomp) then
        
            ! Bx:
            call bx_hed_j0(lam,bxj0)
            call bx_hed_j1(lam,bxj1)   
            bxh = bxh + bxj0*htj0(j)  + bxj1*htj1(j)              
            
            ! By:
            call by_hed_j0(lam,byj0)
            call by_hed_j1(lam,byj1)   
            byh = byh + byj0*htj0(j)  + byj1*htj1(j)         
            
            ! Bz:
            call bz_hed_j1(lam,bzj1)   
            bzh = bzh + bzj1*htj1(j)   

        endif ! lbcomp
!
! Compute derivatives with respect to conductivity:
!
        if (linversion) then 
            
            ! dExdSig:
            call dexdsig_hed_j0(lam)
            call dexdsig_hed_j1(lam)         
            dexhdsig = dexhdsig + dfdsigj0*htj0(j)  + dfdsigj1*htj1(j)
            
            ! dEydSig:
            call deydsig_hed_j0(lam)
            call deydsig_hed_j1(lam)         
            deyhdsig = deyhdsig + dfdsigj0*htj0(j)  + dfdsigj1*htj1(j)
            
            ! dJzdSig:
            call djzdsig_hed_j1(lam)         
            djzhdsig = djzhdsig + dfdsigj1*htj1(j)
            
            if (lbcomp) then    
            
                ! dBxdsig:
                call dbxdsig_hed_j0(lam)
                call dbxdsig_hed_j1(lam)         
                dbxhdsig = dbxhdsig + dfdsigj0*htj0(j)  + dfdsigj1*htj1(j)          
                
                
                ! dBydsig:
                call dbydsig_hed_j0(lam)
                call dbydsig_hed_j1(lam)         
                dbyhdsig = dbyhdsig + dfdsigj0*htj0(j)  + dfdsigj1*htj1(j)              
                
                ! dBzdsig:
                call dbzdsig_hed_j1(lam)         
                dbzhdsig = dbzhdsig + dfdsigj1*htj1(j)
                            
            endif ! lbcomp
        
        endif  ! linversion


    enddo  ! loop over filter coefficients         
  
!  
! Apply coefficients omitted from kernels, and normalize back by r (from defn of digital HT) 
!   
    exh =  sin2theta/(4d0*pi*mu0*csigsite)*exh/r ! need to normalize back by r value
    eyh =  eyh/r 
    jzh = -sintheta/2d0/pi*jzh/r 
    bxh =  bxh/r/2d0/pi
    byh =  sin2theta/4d0/pi*byh/r
    bzh = -costheta/2d0/pi*bzh/r
    
    if (linversion) then 
    
        dexhdsig =  sin2theta/(4d0*pi*mu0*csigsite)*dexhdsig/r ! need to normalize back by r value
        dexhdsig(iRxlayer) = dexhdsig(iRxlayer) - exh/csigsite  ! augment iRxlayer since csigsite was outside kernel function
        deyhdsig =  deyhdsig/r 
        djzhdsig = -sintheta/2d0/pi*djzhdsig /r 
        dbxhdsig =  dbxhdsig/r/2d0/pi
        dbyhdsig =  sin2theta/4d0/pi*dbyhdsig /r
        dbzhdsig = -costheta/2d0/pi*dbzhdsig /r      
        
    endif
    
! Note that x and y components above are with respect to dipole, not cardinal directions.
! x is at -90 degrees from dipole azimuth, y points along dipole (positive in direction of dipole azimuth)
! so now we need to rotate the fields back to the cardinal directions
!   
    if (rotang.ne.0d0) then
        call rotate_fields(exh,eyh,-1.d0*rotang)
        call rotate_fields(bxh,byh,-1.d0*rotang)   
    endif
    
    if (linversion) then 
    
        do i=1,nlay1D
            call rotate_fields(dexhdsig(i),deyhdsig(i),-1.d0*rotang)
            call rotate_fields(dbxhdsig(i),dbyhdsig(i),-1.d0*rotang)   
        enddo
        
    endif
    
    end subroutine comp_hed_fht
    
!==============================================================================!
!================================================================ comp_ved_fht !
!==============================================================================!
    subroutine comp_ved_fht
!
! Compute the ved fields at current site using a fast Hankel transform
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!          
    implicit none
    
!
! Local variables:
!   
    integer    :: j   
    real(8)    :: lam
    complex(8) :: fj0, fj1  
    
    iptr = 1  ! reset for spline interpolation
    
!
! Loop over digital Hankel transform filter coefficients:
!       
    do j = 1,ndhtfc

        lam = base(j)/r     
        
      ! N-layer VED potential:
        call ved_pot_nlayer(lam)

        ! Ex: 
        call ex_ved_j1(lam,fj1)   
        ext = ext + fj1*htj1(j)  
                
        ! Note that Ex and Ey, and Bx and By are related by (co)sines, 
        ! so we only need to compute one component of each

        ! Jz
        call jz_ved_j0(lam,fj0)   
        jzv = jzv + fj0*htj0(j)           
        
        if (lbcomp) then
            ! Bx:
            call bx_ved_j1(lam,fj1)
            bxt = bxt + fj1*htj1(j)               
            
            ! Bz is always zero for a ved    
        
        endif
        
        if (linversion) then
        
            ! dExdsig: 
            call dexdsig_ved_j1(lam)   
            dextdsig = dextdsig + dfdsigj1*htj1(j)  
                    
            ! Note that Ex and Ey, and Bx and By are related by (co)sines, 
            ! so we only need to compute one component of each
    
            ! dJzdsig
            call djzdsig_ved_j0(lam)   
            djzvdsig = djzvdsig + dfdsigj0*htj0(j)           
            
            if (lbcomp) then
                ! dBxdsig:
                call dbxdsig_ved_j1(lam)
                dbxtdsig = dbxtdsig + dfdsigj1*htj1(j)               
                
                ! Bz is always zero for a ved    
            
            endif
        
        endif ! linversion

    enddo  ! loop over filter coefficients          
!  
! Apply coefficients omitted from kernels, and normalize back by r (from defn of digital HT) 
!   
    exv = -dcos(thetaRx)/(2d0*pi*mu0*csigsite)*ext/r 
    eyv = -dsin(thetaRx)/(2d0*pi*mu0*csigsite)*ext/r 
    jzv =  1d0/(2d0*pi)*jzv/r 
    bxv = -dsin(thetaRx)/2d0/pi*bxt/r
    byv =  dcos(thetaRx)/2d0/pi*bxt/r
    
    if (linversion) then
    
        dexvdsig = -dcos(thetaRx)/(2d0*pi*mu0*csigsite)*dextdsig/r 
        dexvdsig(iRxlayer) = dexvdsig(iRxlayer) - exv/csigsite
        
        deyvdsig = -dsin(thetaRx)/(2d0*pi*mu0*csigsite)*dextdsig/r 
        deyvdsig(iRxlayer) = deyvdsig(iRxlayer) - eyv/csigsite
        
        djzvdsig =  1d0/(2d0*pi)*djzvdsig/r 
        dbxvdsig = -dsin(thetaRx)/2d0/pi*dbxtdsig/r
        dbyvdsig =  dcos(thetaRx)/2d0/pi*dbxtdsig/r     
    
    endif
    
    end subroutine comp_ved_fht 
    
!==============================================================================!
!=============================================================== rotate_fields !
!==============================================================================!
    subroutine rotate_fields(ex,ey,rotang)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    complex(8),intent(inout)    :: ex, ey
    real(8),intent(in)          :: rotang
!
! Local variables:
!
    complex(8) :: tmpx, tmpy
    real(8)    :: cc, ss    
!    
! Rotate fields by rotang angle
!
    cc = dcos(rotang)
    ss = dsin(rotang)
    tmpx  =  ex*cc + ey*ss
    tmpy  = -ex*ss + ey*cc
    ex    =  tmpx
    ey    =  tmpy
    
    end subroutine rotate_fields

!==============================================================================!
!========================================================== allocate_dipole1d  !
!==============================================================================!
    subroutine allocate_dipole1d 
!
! Allocates arrays needed by Dipole1D
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
    
    implicit none
    
    integer i
!
! First check that the model depths are correctly defined
!
    do i = 3, nlay1D ! first layer top depth is never used
        if (zlay1D(i).le.zlay1D(i-1)) then
           write(*,*) ''        
           write(*,*) 'Dipole1D Error: layer depths not increasing, layers: ',i,i-1
           write(*,*) 'Stopping.'
           write(*,*) ''
           stop
        endif
    enddo
!
! Allocate storage
!
    allocate ( a(nlay1D),b(nlay1D),c(nlay1D),d(nlay1D),csig(nlay1D),h(nlay1D) )
    allocate ( gamma(nlay1D),k2(nlay1D),expmgh(nlay1D)) 
    allocate ( Rm(nlay1D),Rp(nlay1D),Sm(nlay1D),Sp(nlay1D) )   
!
!  Allocate storage if derivatives also being computed   
!
    if (linversion) then
        allocate ( dRmdsig(nlay1D,nlay1D), dRpdsig(nlay1D,nlay1D) , dSmdsig(nlay1D,nlay1D), dSpdsig(nlay1D,nlay1D) )
        allocate ( dadsig(nlay1D,nlay1D), dbdsig(nlay1D,nlay1D), dcdsig(nlay1D,nlay1D), dddsig(nlay1D,nlay1D) )
        allocate ( dAydsig(nlay1D), dAzdsig(nlay1D)  )
        allocate ( dAydsigdz(nlay1D), dAzdsigdz(nlay1D),dAzdsigdz2(nlay1D)   )   
        allocate ( dae1dsig(nlay1D), dbe2dsig(nlay1D), dce1dsig(nlay1D), dde2dsig(nlay1D) )
        
        allocate ( dexhdsig(nlay1D), deyhdsig(nlay1D), djzhdsig(nlay1D) )
        allocate ( dbxhdsig(nlay1D), dbyhdsig(nlay1D), dbzhdsig(nlay1D) )
        allocate ( dexvdsig(nlay1D), deyvdsig(nlay1D), djzvdsig(nlay1D) )
        allocate ( dbxvdsig(nlay1D), dbyvdsig(nlay1D), dbzvdsig(nlay1D) )
        allocate ( dfdsigj0(nlay1D), dfdsigj1(nlay1D) )
        allocate ( dextdsig(nlay1D), dbxtdsig(nlay1D) )

    endif
        
    
    
    end subroutine allocate_dipole1d  

    
!==============================================================================!
!========================================================= initialize_dipole1d !
!==============================================================================!
    subroutine initialize_dipole1d
!
! Initializes a few things for the 1D computations. Checks to make 
! sure the model is properly defined and also looks for whether
! spatial or (kx,y,z) domain fields requested.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu

    use HankelFilters
    use SinCosFilters   
    
    implicit none
    
    integer i
        
 !
 ! Initialize a few parameters for the potential computations
 ! 
 
    omega         = 2d0*pi*ftx1D
    csig          = sig1D - II*omega*eps ! use complex conductivity locally         
    k2(1:nlay1D)  = -ii*omega*mu0*csig(1:nlay1D) !  using exp(-iwt) time dependence
    
    lhed = .false.
    lved = .false.
!
! Compute layer thicknesses 
!
    h = 1d150 
    do i = 2,nlay1D-1
        h(i) = zlay1D(i+1) - zlay1D(i) 
    enddo

!
! Setup some parameters depending on which output domain requested
!
    select case (trim(outputdomain1D))
!
! Spatial domain solution requested:
!
    case ('spatial')
!
! Initialize HT filter coefficients:
!
        call init_HankelFilters
!
! Set which Hankel Transform filters to use:
!
        select case (trim(HTmethod1D))

        case ('fk_ht_61') 
            ndhtfc = 61
            base(1:61) = fk_ht_base_61
            htj0(1:61) = fk_ht_j0_61
            htj1(1:61) = fk_ht_j1_61
            
        case ('fk_ht_241') 
            ndhtfc = 241
            base(1:241) = fk_ht_base_241
            htj0(1:241) = fk_ht_j0_241
            htj1(1:241) = fk_ht_j1_241       

        case ('kk_ht_101') 
            ndhtfc = 101
            base(1:101) = kk_ht_base_101
            htj0(1:101) = kk_ht_j0_101
            htj1(1:101) = kk_ht_j1_101       

        case ('kk_ht_201') 
            ndhtfc = 201
            base(1:201) = kk_ht_base_201
            htj0(1:201) = kk_ht_j0_201
            htj1(1:201) = kk_ht_j1_201   

        case ('kk_ht_401')   ! Beware
            ndhtfc = 401
            base(1:401) = kk_ht_base_401
            htj0(1:401) = kk_ht_j0_401
            htj1(1:401) = kk_ht_j1_401               
           
                        
        case default ! Use 201 as the default
            write(*,*) ' Default value: Using kk_ht_201 point filters'
            ndhtfc = 201
            base(1:201) = kk_ht_base_201
            htj0(1:201) = kk_ht_j0_201
            htj1(1:201) = kk_ht_j1_201   
          
        end select ! case (trim(HTmethod1D))
        
!
! Set logical flags for computing spatial domain ved and hed
!
        if ( sin(dipTx1D*pi/180.).eq.0. ) then  ! if dip is 0 or 180, then no VED computation needed
            lved = .false.
            lhed = .true.
        elseif ( sin(dipTx1D*pi/180.).eq.1. ) then  ! kwk debug: change this and above to tolerance statements
            lved = .true.
            lhed = .false. 
        else  ! both ved and hed needed
            lved = .true.
            lhed = .true.
        endif
!
! (kx,y,z) domain solution requested:
!
    case ('kx')
!
!  Get mode to compute:
!
        kxmode = kxmode1D  ! assign to local array

!
! kxmode  = 1   HED pointing along x
!         = 2   HED pointing along y
!         = 3   VED pointing along z
!
        if (kxmode /= 3) lhed = .true.
        if (kxmode == 3) lved = .true.

!
! Initialize CT filter coefficients:
!
        call init_SinCosFilters
!
! Set the (co)sine transform method:
!
        select case (trim(CTmethod1D))
        
        case ('kk_ct_81') 
            ncsfc    = 81
            basecsfc(1:81) = kk_ct_base_81
            cosfc(1:81)    = kk_ct_cos_81
            sinfc(1:81)    = kk_ct_sin_81
            
        case ('kk_ct_241') 
            ncsfc    = 241
            basecsfc(1:241) = kk_ct_base_241
            cosfc(1:241)    = kk_ct_cos_241
            sinfc(1:241)    = kk_ct_sin_241      
            
        case ('kk_ct_601') 
            ncsfc    = 601
            basecsfc(1:601) = kk_ct_base_601
            cosfc(1:601)    = kk_ct_cos_601
            sinfc(1:601)    = kk_ct_sin_601    
         
         case default
            ncsfc    = 601
            basecsfc(1:601) = kk_ct_base_601
            cosfc(1:601)    = kk_ct_cos_601
            sinfc(1:601)    = kk_ct_sin_601             
        end select
                 
    case default
        write(*,*) 'Error in Dipole1D, bad outputdomain1D parameter: ', &
         &  outputdomain1D
         stop
         
    end select  ! case (trim(outputdomain1D))
    
    
    end subroutine initialize_dipole1d  

!==============================================================================!
!=============================================================== initialize_Tx !
!==============================================================================!
    subroutine initialize_Tx
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! DGM Feb 2010 - separated from initialize_dipole1d for finite dipole support
    
    implicit none
    
    integer j
!
! Get the layer with the transmitter
!
    iTxlayer = 1 ! default
    do j = 2,nlay1D
       if (zTx1D.gt.zlay1D(j))  then
            iTxlayer = j
       endif
    enddo
!
! Set transmitter depth in iTxlayer
!
    if (iTxlayer.eq.1) then
        depthTx = 1d150 ! large value allows for zeroing source term if in top layer
    else
        depthTx = zTx1D - zlay1D(iTxlayer)   
    endif
!
! Set transmitter height in iTxlayer:
!
    if (iTxlayer.eq.nlay1D) then
        heightTx = 1d150 ! large value allows for zeroing source term if in bottom layer
    else
        heightTx = zlay1D(iTxlayer+1) - zTx1D;   
    endif  
    
    end subroutine initialize_Tx

!==============================================================================!
!=========================================================== deallocate_spline !
!==============================================================================!
    subroutine deallocate_spline
! DGM Feb 2010 - separated from deallocate_dipole1d to support finite dipole
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
    
    implicit none

    if (.not. lUseSpline1D) then
        return
    endif
    
    deallocate( iRxLayerInterp )
    deallocate(  lam_interp )
    deallocate(  a_hed_interp,  b_hed_interp )
    deallocate(  c_hed_interp,  d_hed_interp )
    deallocate( a2_hed_interp, b2_hed_interp )
    deallocate( c2_hed_interp, d2_hed_interp )  
    deallocate(  c_ved_interp,  d_ved_interp )
    deallocate( c2_ved_interp, d2_ved_interp )
    
    if (linversion) then
        deallocate(  dadsig_hed_interp,  dbdsig_hed_interp )
        deallocate(  dcdsig_hed_interp,  dddsig_hed_interp )
        deallocate( dadsig2_hed_interp, dbdsig2_hed_interp )
        deallocate( dcdsig2_hed_interp, dddsig2_hed_interp )
        deallocate(  dcdsig_ved_interp,  dddsig_ved_interp )
        deallocate( dcdsig2_ved_interp, dddsig2_ved_interp )
    endif
    
    end subroutine deallocate_spline

!==============================================================================!
!========================================================= deallocate_dipole1d !
!==============================================================================!
    subroutine deallocate_dipole1d
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
    
    implicit none
    
!
! Deallocate storage
!
    deallocate ( a,b,c,d,csig,h )
    deallocate ( gamma,k2 ) 
    deallocate ( Rm,Rp,Sm,Sp,expmgh  ) 
!
!  Deallocate storage if derivatives also being computed   
!
    if (linversion) then
        deallocate ( dRmdsig, dRpdsig, dSmdsig, dSpdsig )
        deallocate ( dadsig, dbdsig, dcdsig, dddsig, dAydsig, dAzdsig, dAydsigdz, dAzdsigdz, dAzdsigdz2)
        deallocate ( dae1dsig, dbe2dsig, dce1dsig, dde2dsig )
        
        deallocate ( dexhdsig, deyhdsig, djzhdsig )
        deallocate ( dbxhdsig, dbyhdsig, dbzhdsig )
        deallocate ( dexvdsig, deyvdsig, djzvdsig )
        deallocate ( dbxvdsig, dbyvdsig, dbzvdsig )
        deallocate ( dfdsigj0, dfdsigj1 )
        deallocate ( dextdsig, dbxtdsig )
    endif
        
    end subroutine deallocate_dipole1d  

!==============================================================================!
!=================================================================== setupsite !
!==============================================================================!
    subroutine setupsite(i)
!
! Sets up a few geometric parameters for site i:
! location, rotation of dipole, etc.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    
    implicit none
    
    integer, intent(in) :: i
    integer             :: j

!
! Get r and theta to the site in the rotated coordinate system
!
    dx         = x1D(i) - xTx1D
    dy         = y1D(i) - yTx1D
    z          = z1D(i)                ! internal potential routines use r and z
    r          = sqrt( dx**2 + dy**2)  ! Horizontal range to site 
    
    if ((r).lt.1d0) then  
        r = 1d0 ! Move r to 1 so that division by r doesn't explode anywhere
        ! For now this will be a nudge to y=1,x=0 but keep in mind that this will cause small 
        ! error in r=0 fields for inversions...but then again this is where point dipole approximation
        ! is inaccurate, and source-receiver navigation uncertainty will be relatively greatest.
        thetaRx    = 90.d0 ! Angle from x to site fixed to 90 degrees
    else  
        thetaRx    = 180d0/pi*atan2(dy,dx) ! Angle from x to site
    endif
    rotang     = azimuthTx1D - 90.d0     ! Angle from x to rotated x axis
    theta      = thetaRx - rotang        ! Angle from rotated x axis to site
    rotang     = rotang*pi/180.d0        ! convert to radians
    theta      = theta*pi/180.d0         ! convert to radians
    thetaRx    = thetaRx*pi/180.d0


    sintheta   = sin(theta)
    sin2theta  = sin(2.d0*theta)
    costheta   = cos(theta)  
    cos2theta  = cos(2.d0*theta)
    !
! Get layer current site resides in
!
    iRxlayer = 1
    do j= 2,nlay1D
      if (z.gt.zlay1D(j))  then   ! Sites on boundaries use the top layer 
        iRxlayer = j
      endif
    enddo
!
! Set receiver depth in iRxlayer
!
    if (iRxlayer.eq.1) then
        depthRx = 1d150 ! large value allows for zeroing source term if in top layer
    else
        depthRx = z - zlay1D(iRxlayer)   
    endif
!
! Set receiver height in iRxlayer:
!
    if (iRxlayer.eq.nlay1D) then
        heightRx = 1d150 ! large value allows for zeroing source term if in bottom layer
    else
        heightRx = zlay1D(iRxlayer+1) - z;   
    endif  
    
    
!       
! Set sgn operator for source layer primary term and it's derivative 
!
    if (iRxlayer == iTxlayer) then
      isgnsrc = 1                   ! potential has primary source term
      if (z > zTx1D) then
        isgndsrcdz = -1             ! negative source term derivative
      elseif (z < zTx1D) then
        isgndsrcdz = 1              ! positive source term derivative
      elseif (z == zTx1D) then
        isgndsrcdz = 0              ! no source term derivative
      endif
    else
      isgnsrc = 0                 ! potential does not have primary source term
      isgndsrcdz = 0              ! no source term derivative
    endif           

    
! Debug print checks:
    ! write(*,*) 'kxmode: ',kxmode
    ! write(*,*) 'x,y,z,yTz,zTx,layers ',x1D(i),y1D(i),z1D(i),yTx1D,zTx1D,iTxlayer,iRxlayer
    ! write(*,*) 'dx,dy,z,r:',dx,dy,z,r
    ! write(*,*) 'theta,thetaRx,rotang: ',theta*180/pi,thetaRx*180/pi,rotang*180/pi
    ! write(*,*) 'heightTx,depthTx: ',heightTx,depthTx
    ! write(*,*) 'lved,lhed: ',lved,lhed
    ! write(*,*) 'h: ',h
    ! write(*,*) 'isgnsrc, isgndsrcdz:',isgnsrc, isgndsrcdz
        
        
    end subroutine setupsite    
    
!==============================================================================!
!================================================================= ex_hed_x_kx !
!==============================================================================!
    subroutine ex_hed_x_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = II*omega*ay - (kx*kx)/(mu0*csigsite)*(aypdazdz)
    
    end subroutine ex_hed_x_kx  
    
!==============================================================================!
!================================================================= ey_hed_x_kx !
!==============================================================================!
    subroutine ey_hed_x_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = -(kx*ky)/(mu0*csigsite)*(aypdazdz)
    
    end subroutine ey_hed_x_kx      

!==============================================================================!
!================================================================= jz_hed_x_kx !
!==============================================================================!
    subroutine jz_hed_x_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = csigsite*II*omega*az*II*kx + II*kx/mu0*(daydzpdazdz2) ! Using Jz = sig1D*Ez
    ! note that the first term needs the ik_x since az is actually Lambda_z 
    
    end subroutine jz_hed_x_kx      

!==============================================================================!
!================================================================= bx_hed_x_kx !
!==============================================================================!
    subroutine bx_hed_x_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = -kx*ky*az  
    
    end subroutine bx_hed_x_kx  
    
!==============================================================================!
!================================================================= by_hed_x_kx !
!==============================================================================!
    subroutine by_hed_x_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = kx*kx*az + daydz
    
    end subroutine by_hed_x_kx      

!==============================================================================!
!================================================================= bz_hed_x_kx !
!==============================================================================!
    subroutine bz_hed_x_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = -II*ky*ay  
    
    end subroutine bz_hed_x_kx      
    
!==============================================================================!
!================================================================= ex_hed_y_kx !
!==============================================================================!
    subroutine ex_hed_y_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = -(kx*ky)/(mu0*csigsite)*(aypdazdz)
    
    end subroutine ex_hed_y_kx  
    
!==============================================================================!
!================================================================= ey_hed_y_kx !
!==============================================================================!
    subroutine ey_hed_y_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = II*omega*ay-(ky*ky)/(mu0*csigsite)*(aypdazdz)
    
    end subroutine ey_hed_y_kx      

!==============================================================================!
!================================================================= jz_hed_y_kx !
!==============================================================================!
    subroutine jz_hed_y_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = csigsite*II*omega*az*II*ky + II*ky/mu0*(daydzpdazdz2)
    ! note that the first term needs the ik_y since az is actually Lambda_z
    
    end subroutine jz_hed_y_kx      

!==============================================================================!
!================================================================= bx_hed_y_kx !
!==============================================================================!
    subroutine bx_hed_y_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = -ky*ky*az - daydz
    
    end subroutine bx_hed_y_kx  
    
!==============================================================================!
!================================================================= by_hed_y_kx !
!==============================================================================!
    subroutine by_hed_y_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = kx*ky*az
    
    end subroutine by_hed_y_kx      

!==============================================================================!
!================================================================= bz_hed_y_kx !
!==============================================================================!
    subroutine bz_hed_y_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = II*kx*ay  
    
    end subroutine bz_hed_y_kx          

!==============================================================================!
!================================================================= ==ex_ved_kx !
!==============================================================================!
    subroutine ex_ved_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = II*kx/(mu0*csigsite)*dazdz
    
    end subroutine ex_ved_kx    
    
!==============================================================================!
!=================================================================== ey_ved_kx !
!==============================================================================!
    subroutine ey_ved_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = II*ky/(mu0*csigsite)*dazdz

    end subroutine ey_ved_kx        

!==============================================================================!
!=================================================================== jz_ved_kx !
!==============================================================================!
    subroutine jz_ved_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = csigsite*II*omega*az + 1.d0/mu0*dazdz2 ! June 2008, fixed sign error here!
    
    end subroutine jz_ved_kx        

!==============================================================================!
!=================================================================== bx_ved_kx !
!==============================================================================!
    subroutine bx_ved_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = II*ky*az 
    
    end subroutine bx_ved_kx    
    
!==============================================================================!
!=================================================================== by_ved_kx !
!==============================================================================!
    subroutine by_ved_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = -II*kx*az
    
    end subroutine by_ved_kx        

!==============================================================================!
!=================================================================== bz_ved_kx !
!==============================================================================!
    subroutine bz_ved_kx(kx,ky,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: kx,ky
    
    f = 0
    
    end subroutine bz_ved_kx        
    
!==============================================================================!
!=================================================================== ex_ved_j1 !
!==============================================================================!
    subroutine ex_ved_j1(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = dazdz*lam**2
    
    end subroutine ex_ved_j1    
    
!==============================================================================!
!============================================================== dexdsig_ved_j1 !
!==============================================================================!
    subroutine dexdsig_ved_j1(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj1 = dazdsigdz*lam**2
    
    end subroutine dexdsig_ved_j1    

!==============================================================================!
!=================================================================== bx_ved_j1 !
!==============================================================================!
    subroutine bx_ved_j1(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = az*lam**2
    
    end subroutine bx_ved_j1    
 
!==============================================================================!
!============================================================== dbxdsig_ved_j1 !
!==============================================================================!
    subroutine dbxdsig_ved_j1(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj1 = dazdsig*lam**2
    
    end subroutine dbxdsig_ved_j1    
    
!==============================================================================!
!=================================================================== jz_ved_j0 !
!==============================================================================!
    subroutine jz_ved_j0(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = (csigsite*II*omega*az + 1.d0/mu0*dazdz2 )*lam
    
    end subroutine jz_ved_j0        
    
!==============================================================================!
!============================================================== djzdsig_ved_j0 !
!==============================================================================!
    subroutine djzdsig_ved_j0(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj0 = (csigsite*II*omega*dazdsig + 1.d0/mu0*dazdsigdz2 )*lam
    dfdsigj0(iRxlayer) =  dfdsigj0(iRxlayer) + (II*omega*az )*lam
    
    end subroutine djzdsig_ved_j0      

!==============================================================================!
!=================================================================== ex_hed_j0 !
!==============================================================================!
    subroutine ex_hed_j0(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = -(aypdazdz)*lam**3
    
    end subroutine ex_hed_j0
 
!==============================================================================!
!============================================================== dexdsig_hed_j0 !
!==============================================================================!
    subroutine dexdsig_hed_j0(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj0 = -(daydsig + dazdsigdz)*lam**3
    
    end subroutine dexdsig_hed_j0
    
!==============================================================================!
!=================================================================== ex_hed_j1 !
!==============================================================================!
    subroutine ex_hed_j1(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = (aypdazdz)*2d0*lam**2/r
    
    end subroutine ex_hed_j1        
    
!==============================================================================!
!============================================================== dexdsig_hed_j1 !
!==============================================================================!
    subroutine dexdsig_hed_j1(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj1 = (daydsig + dazdsigdz)*2d0*lam**2/r
    
    end subroutine dexdsig_hed_j1           
    
!==============================================================================!
!=================================================================== ey_hed_j0 !
!==============================================================================!
    subroutine ey_hed_j0(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = ii*omega*ay*lam/2d0/pi &
        -1d0/(2d0*pi*mu0*csigsite)*(aypdazdz)*lam**3*(sintheta)**2
    
    end subroutine ey_hed_j0
!==============================================================================!
!============================================================== deydsig_hed_j0 !
!==============================================================================!
    subroutine deydsig_hed_j0(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj0 = ii*omega*dAydsig*lam/2d0/pi &
             -1d0/(2d0*pi*mu0*csigsite)*(dAydsig+dAzdsigdz)*lam**3*(sintheta)**2 
        
    dfdsigj0(iRxlayer) = dfdsigj0(iRxlayer) + 1d0/(2d0*pi*mu0*csigsite**2)*(aypdazdz)*lam**3*(sintheta)**2
       
    end subroutine deydsig_hed_j0
        
!==============================================================================!
!=================================================================== ey_hed_j1 !
!==============================================================================!
    subroutine ey_hed_j1(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = -1d0/(2d0*pi*mu0*csigsite)*(aypdazdz)*lam**2*cos2theta/r
    
    end subroutine ey_hed_j1    
    
!==============================================================================!
!============================================================== deydsig_hed_j1 !
!==============================================================================!
    subroutine deydsig_hed_j1(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    

    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj1 = -1d0/(2d0*pi*mu0*csigsite)*(dAydsig+dAzdsigdz)*lam**2*cos2theta/r !
    dfdsigj1(iRxlayer) =  dfdsigj1(iRxlayer) + 1d0/(2d0*pi*mu0*csigsite**2)*(aypdazdz)*lam**2*cos2theta/r
    
    end subroutine deydsig_hed_j1  
    
!==============================================================================!
!=================================================================== jz_hed_j1 !
!==============================================================================!
    subroutine jz_hed_j1(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = (csigsite*ii*omega*az + 1.d0/mu0 *(daydzpdazdz2) )*lam**2
    
    end subroutine jz_hed_j1    

!==============================================================================!
!============================================================== djzdsig_hed_j1 !
!==============================================================================!
    subroutine djzdsig_hed_j1(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj1 = (csigsite*ii*omega*dazdsig + 1.d0/mu0 *(daydsigdz + dazdsigdz2) )*lam**2
    dfdsigj1(iRxlayer) =  dfdsigj1(iRxlayer) + ii*omega*az*lam**2
    
    end subroutine djzdsig_hed_j1 
        
!==============================================================================!
!=================================================================== bx_hed_j0 !
!==============================================================================!
    subroutine bx_hed_j0(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = -lam*daydz - az*lam**3*sintheta**2
    
    end subroutine bx_hed_j0
        
!==============================================================================!
!============================================================== dbxdsig_hed_j0 !
!==============================================================================!
    subroutine dbxdsig_hed_j0(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj0 = -lam*daydsigdz - dazdsig*lam**3*sintheta**2
    
    end subroutine dbxdsig_hed_j0
     
!==============================================================================!
!=================================================================== bx_hed_j1 !
!==============================================================================!
    subroutine bx_hed_j1(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = -az*lam**2*cos2theta/r
    
    end subroutine bx_hed_j1        

     
!==============================================================================!
!============================================================== dbxdsig_hed_j1 !
!==============================================================================!
    subroutine dbxdsig_hed_j1(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj1 = -dazdsig*lam**2*cos2theta/r
    
    end subroutine dbxdsig_hed_j1        

!==============================================================================!
!=================================================================== by_hed_j0 !
!==============================================================================!
    subroutine by_hed_j0(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = az*lam**3
    
    end subroutine by_hed_j0
 
!==============================================================================!
!============================================================== dbydsig_hed_j0 !
!==============================================================================!
    subroutine dbydsig_hed_j0(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj0 = dazdsig*lam**3
    
    end subroutine dbydsig_hed_j0
       
!==============================================================================!
!=================================================================== by_hed_j1 !
!==============================================================================!
    subroutine by_hed_j1(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = -az*2d0*lam**2/r
    
    end subroutine by_hed_j1    
       
!==============================================================================!
!============================================================== dbydsig_hed_j1 !
!==============================================================================!
    subroutine dbydsig_hed_j1(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj1 = -dazdsig*2d0*lam**2/r
    
    end subroutine dbydsig_hed_j1   

!==============================================================================!
!=================================================================== bz_hed_j1 !
!==============================================================================!
    subroutine bz_hed_j1(lam,f)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    complex(8), intent(out) :: f
    real(8), intent(in)     :: lam
    
    f = ay*lam**2
    
    end subroutine bz_hed_j1    

!==============================================================================!
!============================================================== dbzdsig_hed_j1 !
!==============================================================================!
    subroutine dbzdsig_hed_j1(lam)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!    
    
    implicit none
    real(8), intent(in)     :: lam
    
    dfdsigj1 = daydsig*lam**2
    
    end subroutine dbzdsig_hed_j1    

!==============================================================================!
!============================================================= hed_coef_nlayer !
!==============================================================================!
    subroutine hed_coef_nlayer(lam)
!
! Generates the vector potential coefficients for Ay and Az
! for an N layered model with a horizontal electric dipole transmitter 
! located in any layer.
! Note that the Az coefficients are for lambda_z, where 
! Az = d/dy (lambda_z), so there is a ik_y term left off that
! is accounted for later in the Hankel or (Co)Sine transform formulations.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Version 3.1  May 2008.          Modified df/dsig recursions for speedup.
! Version 3.0  March 11, 2008.    Added recursions for df/dsig kernels for inversion
! Version 2.1  December 14, 2007. Modified to use exp(-iwt) time dependence
! Version 2.0  October 19, 2007.  Using new formula for Lambda_z potential
! Version 1.1  Sept, 2007.        Modifications for speed-up and stability.
! Version 1.0  June 12, 2007.
!
    
    implicit none
!
! Input arguments:
!
    real(8), intent(in)   :: lam    ! Hankel transform parameter
! 
! Local variables:
! 
    integer    :: i

! Forward variables:                
    complex(8) :: gmogp, rjexp, onemrmrp, rmrp, srcp,srcm, rhs
    complex(8) :: sgmosgp, sjexp, smsp, aytop  
    
! Derivative variables:
    complex(8) :: drdsig, dgammadsig, dsdsig, p,dp,rr,dr,q,s,t,u,du,ff,srcpInv,srcmInv, db
    integer    :: isgnRm, isgnRp   
!
! Compute gamma for each layer
!
    gamma(1:nlay1D) = sqrt( lam**2 + k2(1:nlay1D) )       
!
! Compute attenuation factor for each layer
!
    expmgh = exp(-gamma*h)
    
    where (abs(expmgh) <1d-80)  expmgh = 0d0    ! avoid tiny numbers
!
! Initialize some arrays
!
    Rm = 0d0
    Sm = 0d0
    Rp = 0d0
    Sp = 0d0
    a  = 0d0
    b  = 0d0
    c  = 0d0
    d  = 0d0
    
    if (linversion) then
    
        dRmdsig = 0d0 !  note that for bounding layers these terms start out as zero
        dRpdsig = 0d0
        dSmdsig = 0d0
        dSpdsig = 0d0   
         
         dadsig = 0d0
         dbdsig = 0d0
         dcdsig = 0d0
         dddsig = 0d0
        
    endif
    
!
! Compute the recursion coefficients for the Ay and Lambda_z potentials
! for the layers above the transmitter:
!

    do i = 2,iTxlayer
      Rm(i-1) = Rm(i-1)*expmgh(i-1) ! note that I post apply expmgh(i-1)
      Sm(i-1) = Sm(i-1)*expmgh(i-1) ! this so that for the iTxlayer the exponent
      ! is omitted and we then don't need the positive exponent in the formula for a_itxlayer and b_itxlayer
   
      gmogp   = (-ii*omega*mu0*(csig(i)-csig(i-1) )) / (gamma(i)+gamma(i-1))**2  
      sgmosgp =  ( gamma(i)*csig(i-1) - gamma(i-1)*csig(i) ) / ( gamma(i)*csig(i-1) + gamma(i-1)*csig(i) )
      
      rjexp  = Rm(i-1)*expmgh(i-1)   ! then we apply attenuation to bottom of layer
      sjexp  = Sm(i-1)*expmgh(i-1)         
    
      Rm(i)  =  (gmogp + rjexp)   / ( 1.d0 + gmogp*rjexp)   ! and compute the recursion coeff for the current layer,
                                                            ! sans expmgh(i) term
      Sm(i)  =  (sgmosgp + sjexp) / ( 1.d0 + sgmosgp*sjexp)     
        
    
        if (linversion) then
             
            dgammadsig   = -ii*omega*mu0 / (2.d0 *gamma(i))
            drdsig       = 2.d0*gamma(i-1) / ((gamma(i) + gamma(i-1) )**2) *dgammadsig
            
            ! Post-apply expmgh(i-1) to  all recursion coefficients (done to avoid exp(+) term in final formula):
            dRmdsig(1:i-1,i-1) = dRmdsig(1:i-1,i-1)*expmgh(i-1)
    
            ! Initialize dRm(i)/dsig(i)    
            dRmdsig(i,i) = drdsig* (1.d0 - Rm(i)*rjexp) / (1.d0 + gmogp*rjexp) - (1-int(i/nlay1d) )*Rm(i)*h(i)*dgammadsig   
            ! the (int) term deals with Tx in bottom layer
          
            ! Next step is to compute dRm(:)/dsig(i-1) for previous layers
            dgammadsig   = -ii*omega*mu0 / (2.d0 *gamma(i-1))
            drdsig       = -2.d0*gamma(i) / ((gamma(i) + gamma(i-1) )**2) *dgammadsig ! actually dr(i)dsig(i-1)
            
            db           = dRmdsig(i-1,i-1)*expmgh(i-1) - rjexp*h(i-1)*dgammadsig

            dRmdsig(i-1,i)  =   ( drdsig*(1.d0 - Rm(i)*rjexp)  + db*(1.d0 - Rm(i)*gmogp) ) / (1.d0 + gmogp*rjexp)
            
            ! Lastly, apply chain rule for layers 1:i-2 to propagate deriv for layers already initialized on previous iterations:
            dRmdsig(1:i-2,i)    = dRmdsig(1:i-2,i-1) * ( expmgh(i-1)*(1.d0 - Rm(i)*gmogp) )/(1.d0 + gmogp*rjexp) 
    
            ! Now do the same thing for Sm terms:
            dgammadsig      = -ii*omega*mu0 / (2.d0 *gamma(i))
            dsdsig          = (dgammadsig*csig(i-1)*(1.d0 - sgmosgp)  - gamma(i-1)*(1.d0+ sgmosgp) ) /  & 
                            & ( gamma(i)*csig(i-1) + gamma(i-1)*csig(i) )
            
            ! Post-apply expmgh(i-1) to  all recursion coefficients (done to avoid exp(+) term in final formula):
            dSmdsig(1:i-1,i-1) = dSmdsig(1:i-1,i-1)*expmgh(i-1)
            
            ! Initialize dSm(i)/dsig(i):        
            dSmdsig(i,i)    = dsdsig*(1.d0 - Sm(i)*sjexp) / (1.d0 + sgmosgp*sjexp) - (1- int(i/nlay1d) )*Sm(i)*h(i) *dgammadsig 
            ! the (int) term deals with Tx in bottom layer
            
            ! Next step is to compute dSm(:)/dsig(i-1) for previous layers
            dgammadsig  = -ii*omega*mu0 / (2.d0 *gamma(i-1))
            
            dsdsig      = (  gamma(i)*(1.d0 - sgmosgp) - dgammadsig*csig(i)*(1.d0 + sgmosgp) ) /  &
                        & ( gamma(i)*csig(i-1) + gamma(i-1)*csig(i) )   ! actually ds(i)dsig(i-1)
            
            db          = dSmdsig(i-1,i-1)* expmgh(i-1) - sjexp*h(i-1)*dgammadsig
            
            dSmdsig(i-1,i) =    ( dsdsig*(1.d0  - Sm(i)*sjexp)  + db*(1.d0  - Sm(i)*sgmosgp) ) / (1.d0 + sgmosgp*sjexp) 
            
            ! Apply chain rule to propagate deriv for layers already initialized on previous iterations:
            dSmdsig(1:i-2,i) = dSmdsig(1:i-2,i-1) * ( expmgh(i-1)*(1.d0 - Sm(i)*sgmosgp) )/(1.d0 + sgmosgp*sjexp) 
            
            ! At this point both dRmdsig(iTxlayer,1:iTxlayer) and dSmdsig(iTxlayer,1:iTxlayer) 
            ! have an expmgh(iTxlayer) factor omitted
            
        endif ! linversion
        
    enddo
    
!
! Compute the recursion coefficients for the Ay and Lambda_z potentials
! for the layers below the transmitter:
!
  
    do i = nlay1D-1,iTxlayer,-1
       Rp(i+1) = Rp(i+1)*expmgh(i+1)  ! note that I post apply 
       Sp(i+1) = Sp(i+1)*expmgh(i+1)        
      
      ! this so that for the iTxlayer the exponent is omitted and we then don't need the positive exponent in the 
      ! formula for a_i and b_i     
       gmogp   = (-ii*omega*mu0*(csig(i)-csig(i+1) )) / (gamma(i)+gamma(i+1))**2
       sgmosgp = ( gamma(i)*csig(i+1) - gamma(i+1)*csig(i) ) / ( gamma(i)*csig(i+1) + gamma(i+1)*csig(i) )      
          
       rjexp   = Rp(i+1)*expmgh(i+1)  ! since Rp(nlay1D) is zero, no need to worry about height of bottom layer
       sjexp   = Sp(i+1)*expmgh(i+1)
       
       Rp(i)   =  (gmogp + rjexp)   / ( 1.d0 + gmogp*rjexp)
       Sp(i)   =  (sgmosgp + sjexp) / ( 1.d0 + sgmosgp*sjexp)    
   
        if (linversion) then

            dgammadsig = -ii*omega*mu0 / (2.d0 *gamma(i))
            drdsig     = 2.d0*gamma(i+1) / ((gamma(i) + gamma(i+1) )**2) *dgammadsig
            
            ! Post-apply expmgh(i+1) to  all recursion coefficients (done to avoid exp(+) term in final formula):       
            dRpdsig(i+1:nlay1d,i+1) = dRpdsig(i+1:nlay1d,i+1)*expmgh(i+1)       
            
            ! Initialize dRp(i)/dsig(i)    
             dRpdsig(i,i) = drdsig* (1.d0 - Rp(i)*rjexp) / (1.d0 + gmogp*rjexp) -  &
                          & ( 1- int( (nlay1d-i+1)/nlay1d) )*Rp(i)*h(i)*dgammadsig  ! the (int) term deals with Tx in top layer
         
            ! Next step is to compute dRm(:)/dsig(i+1) for previous layers
            dgammadsig  = -ii*omega*mu0 / (2.d0 *gamma(i+1))
            drdsig      = -2.d0*gamma(i) / ((gamma(i) + gamma(i+1) )**2) *dgammadsig ! actually dr(i)dsig(i+1)
            
            db          = dRpdsig(i+1,i+1)*expmgh(i+1) - rjexp*h(i+1)*dgammadsig

            dRpdsig(i+1,i) =    ( drdsig*( 1.d0 - Rp(i)*rjexp)  + db*( 1.d0 - Rp(i)*gmogp) ) / (1.d0 + gmogp*rjexp)
            
            ! Lastly, apply chain rule for layers 1:i-2 to propagate deriv for layers already initialized on previous iterations:
            dRpdsig(i+2:nlay1d,i) = dRpdsig(i+2:nlay1d,i+1) * &
                                           &  ( expmgh(i+1)*(1.d0 - Rp(i)*gmogp) )/(1.d0 + gmogp*rjexp) 
            
            ! Now do the same thing for Sm terms:
            dgammadsig = -ii*omega*mu0 / (2.d0 *gamma(i))
            dsdsig     = (dgammadsig*csig(i+1)*(1.d0 - sgmosgp) - (gamma(i+1)*(1.d0+ sgmosgp)) ) /  &
                       & ( gamma(i)*csig(i+1) + gamma(i+1)*csig(i) )

            ! Post-apply expmgh(i+1) to  all recursion coefficients (done to avoid exp(+) term in final formula):       
            dSpdsig(i+1:nlay1d,i+1)     = dSpdsig(i+1:nlay1d,i+1)*expmgh(i+1)   
            
            ! Initialize dSm(i)/dsig(i):        
            dSpdsig(i,i) = dsdsig* (1.d0 - Sp(i)*sjexp) / (1.d0 + sgmosgp*sjexp) - & 
                         &  ( 1- int( (nlay1d-i+1)/nlay1d) )*Sp(i)*h(i) *dgammadsig ! the (int) term deals with Tx in top layer
            
        ! Next step is to compute dSm(:)/dsig(i+1) for previous layers
            dgammadsig  = -ii*omega*mu0 / (2.d0 *gamma(i+1))
            dsdsig      = ( gamma(i)*(1.d0 - sgmosgp) - dgammadsig*csig(i)*(1.d0 + sgmosgp)  ) / &
                        & ( gamma(i)*csig(i+1) + gamma(i+1)*csig(i) )   ! actually ds(i)dsig(i+1)
                        
            db          = dSpdsig(i+1,i+1)* expmgh(i+1) - sjexp*h(i+1)*dgammadsig
                
            dSpdsig(i+1,i) =    ( dsdsig*( 1.d0 - Sp(i)*sjexp)  + db*( 1.d0 - Sp(i)*sgmosgp) ) / (1.d0 + sgmosgp*sjexp)         
            
            ! Apply chain rule to propagate deriv for layers already initialized on previous iterations:
            dSpdsig(i+2:nlay1d,i) = dSpdsig(i+2:nlay1d,i+1) * ( expmgh(i+1)*(1.d0 - Sp(i)*sgmosgp) ) / & 
                                              (1.d0 + sgmosgp*sjexp) 
        
            ! At this point both dRpdsig(iTxlayer,iTxlayer:nlay1d) and dSpdsig(iTxlayer,iTxlayer:nlay1d) have an 
            ! expmgh(iTxlayer) factor omitted
            
        endif ! linversion
       
    enddo   
  
!
! Compute Ay potential coefficients in the source layer:
!   
    rmrp     = Rm(iTxlayer)*Rp(iTxlayer)
    onemrmrp = 1.d0 - rmrp*expmgh(iTxlayer)*expmgh(iTxlayer)  ! added exp term since I omitted it earlier
    rhs      = mu0*sdm1D/ (2.d0*gamma(iTxlayer))
    srcp     = exp(-gamma(iTxlayer)*heightTx)*rhs   
    srcm     = exp(-gamma(iTxlayer)*depthTx )*rhs 

    a(iTxlayer) = ( rmrp*srcm*expmgh(iTxlayer)  + Rp(iTxlayer)*srcp  ) / onemrmrp  
    b(iTxlayer) = ( rmrp*srcp*expmgh(iTxlayer)  + Rm(iTxlayer)*srcm  ) / onemrmrp  
    
!
! If inversion, compute da/dsig,db/dsig
!
    if (linversion) then
         
         ! Make sure nothing is tiny:
!        where (abs(dRmdsig) <1d-80) dRmdsig = 0d0
!        where (abs(dRpdsig) <1d-80) dRpdsig = 0d0
!        where (abs(dSmdsig) <1d-80) dSmdsig = 0d0
!        where (abs(dSpdsig) <1d-80) dSpdsig = 0d0

  !
  ! Catch to omit recursion derivatives if Tx in upper or lowermost layers
  !
        if (iTxlayer == nlay1D) then
            isgnRm = 0 
            isgnRp = 0          
        elseif (iTxlayer == 1)  then
            isgnRm = 0
            isgnRp = 0  
        else
            isgnRm = 1
            isgnRp = 1
        endif
        
        dgammadsig = -ii*omega*mu0 / (2.d0 *gamma(iTxlayer)) 
        
      ! da/dsig: 
        p  = exp(-gamma(iTxlayer)*heightTx)
        rr = exp(-gamma(iTxlayer)*depthTx )
        dp = -p*heightTx*dgammadsig     
        dr = -rr*depthTx*dgammadsig 
        
        q = Rm(iTxlayer)*expmgh(iTxlayer)
        
        s = Rp(iTxlayer)*expmgh(iTxlayer)
        
        t = Rp(iTxlayer)
        
        u = rhs  
        du = -dgammadsig*rhs/gamma(iTxlayer)  ! term only for iTxlayer deriv
        
        ff = a(iTxlayer)
        
        dadsig(:,iTxlayer) = ( (t*u*rr + ff*s)*dRmdsig(:,iTxlayer)*expmgh(iTxlayer) +  &  
                           &    q*ff*dRpdsig(:,iTxlayer)*expmgh(iTxlayer) + (p+q*rr)*u*dRpdsig(:,iTxlayer) ) / onemrmrp
        
        ! add in terms for deriv in source layer:
        dadsig(iTxlayer,iTxlayer) = dadsig(iTxlayer,iTxlayer) +  &
                                  &  ( t*u*(dp + q*dr) + (p+q*rr)*t*du + isgnRp*(p+q*rr)*u*t*h(iTxlayer)*dgammadsig )/ onemrmrp
        
         
       ! db/dsig:
        p  = exp(-gamma(iTxlayer)*depthTx )
        dp = -p*depthTx*dgammadsig
        rr = exp(-gamma(iTxlayer)*heightTx)
        dr = -rr*heightTx*dgammadsig
        
        q  = Rp(iTxlayer)*expmgh(iTxlayer)  
 
        s  = Rm(iTxlayer)*expmgh(iTxlayer)   
        
        t = Rm(iTxlayer)
        
        ff = b(iTxlayer)    
        
        dbdsig(:,iTxlayer) = ( (t*u*rr + ff*s)*dRpdsig(:,iTxlayer)*expmgh(iTxlayer) + & 
                           &   q*ff*dRmdsig(:,iTxlayer)*expmgh(iTxlayer) +   (p+q*rr)*u*dRmdsig(:,iTxlayer) ) / onemrmrp

        ! add in terms for deriv in source layer:
        dbdsig(iTxlayer,iTxlayer) = dbdsig(iTxlayer,iTxlayer) + &
                                  &  ( t*u*(dp + q*dr) + (p+q*rr)*t*du + isgnRm*(p+q*rr)*u*t*h(iTxlayer)*dgammadsig )/ onemrmrp
         
        !
        ! Compute the source terms for the first step of the upward/downward continuations:
        !   
        
        srcpInv     = -srcp/gamma(iTxlayer)*(gamma(iTxlayer)*heightTx + 1.d0) * dgammadsig
        srcmInv     = -srcm/gamma(iTxlayer)*(gamma(iTxlayer)*depthTx  + 1.d0) * dgammadsig
        
    endif
    
!
! Compute Ay potential coefficients for layers above the source layer:
!
    do i = iTxlayer-1,1,-1
       
        aytop =  (a(i+1)*expmgh(i+1) + b(i+1) + srcm ) 
        a(i)  = aytop / ( 1.d0 + Rm(i)*expmgh(i) )  
        b(i)  = a(i)*Rm(i)
        srcm  = 0d0 ! primary source term only for first layer above source layer
    
    enddo
    
    if (linversion) then    

        do i = iTxlayer-1,1,-1
    
            ! base form for dadsig:
            dadsig(:,i) =  dbdsig(:,i+1) + dadsig(:,i+1)*expmgh(i+1)  - a(i)*dRmdsig(:,i)*expmgh(i)
            
            ! add modified terms for dsig in layer i-1:
            dadsig(i+1,i) = dadsig(i+1,i) - a(i+1)*h(i+1)*expmgh(i+1)*(-ii*omega*mu0 / (2.d0 *gamma(i+1))) 
                
            ! add modified terms for dsig in layer i-1 if i-1 is source layer:
            dadsig(i+1,i)  = dadsig(i+1,i) + srcmInv
            
            ! add terms for dsig in layer i
            dadsig(i,i) = dadsig(i,i) + b(i)*h(i)*expmgh(i)*(-ii*omega*mu0 / (2.d0 *gamma(i)))  
            
            ! lastly divide all terms by this:
            dadsig(:,i) = dadsig(:,i) / (1.d0 + Rm(i)*expmgh(i) )
            
            ! Finally compute dbdsig:
            dbdsig(:,i) = a(i)*dRmdsig(:,i) + Rm(i)*dadsig(:,i)
            
            ! zero out source term for remaining layers
            srcmInv = 0d0           
        
        enddo
        
    endif 
        

!
! Compute Ay potential coefficients for layers below the source layer:
!   
    do i = iTxlayer+1,nlay1D
        
        aytop = (a(i-1) + b(i-1)*expmgh(i-1) + srcp )
        b(i)  = aytop/ (1.d0 + Rp(i)*expmgh(i) )
        a(i)  = b(i)*Rp(i)
        srcp  = 0d0 ! primary source term only for first layer below source layer
    
    enddo
    
    if (linversion) then       
        
        do i = iTxlayer+1,nlay1D    
        
            ! base form for dbdsig:    
            dbdsig(:,i) =  dadsig(:,i-1) + dbdsig(:,i-1)*expmgh(i-1)  - b(i)*dRpdsig(:,i)*expmgh(i)
            
            ! add modified terms for dsig in layer i-1:     
            dbdsig(i-1,i) = dbdsig(i-1,i) - b(i-1)*h(i-1)*expmgh(i-1)*(-ii*omega*mu0 / (2.d0 *gamma(i-1))) 
        
            ! add modified terms for dsig in layer i-1 if i-1 is source layer:          
            dbdsig(i-1,i)  = dbdsig(i-1,i) + srcpInv 
            
            ! add terms for dsig in layer i
            dbdsig(i,i) = dbdsig(i,i) + a(i)*h(i)*expmgh(i)*(-ii*omega*mu0 / (2.d0 *gamma(i)))  
            
            ! lastly divide all terms by this:
            dbdsig(:,i) = dbdsig(:,i) / (1.d0 + Rp(i)*expmgh(i) )
            
            ! Finally compute dadsig:
            dadsig(:,i) = b(i)*dRpdsig(:,i) + Rp(i)*dbdsig(:,i)
            
            ! zero out source term for remaining layers
            srcpInv = 0d0     
     
        enddo
    
    endif


   
!
! Compute \Lambda_z potential coefficients c and d in the source layer:
!
    smsp     = Sm(iTxlayer)*Sp(iTxlayer)
    onemrmrp = 1.d0 - smsp*expmgh(iTxlayer)*expmgh(iTxlayer)  ! added exp term since I omitted it earlier
    rhs      = mu0*sdm1D/ (2.d0*lam**2)
    srcp     = -exp(-gamma(iTxlayer)*heightTx)*rhs   
    srcm     =  exp(-gamma(iTxlayer)*depthTx )*rhs 

    c(iTxlayer) = ( smsp*srcm*expmgh(iTxlayer)  + Sp(iTxlayer)*srcp  ) / onemrmrp  
    d(iTxlayer) = ( smsp*srcp*expmgh(iTxlayer)  + Sm(iTxlayer)*srcm  ) / onemrmrp   
    

!
! If inversion, compute dc/dsig,dd/dsig
!
    if (linversion) then
    
        dgammadsig = -ii*omega*mu0 / (2.d0 *gamma(iTxlayer)) 
    
       ! dc/dsig:
        p  = -exp(-gamma(iTxlayer)*heightTx)
        rr = exp(-gamma(iTxlayer)*depthTx )
        dp = -p*heightTx*dgammadsig 
        dr = -rr*depthTx*dgammadsig 
        
        q  = Sm(iTxlayer)*expmgh(iTxlayer)
        
        s  = Sp(iTxlayer)*expmgh(iTxlayer)
        
        t = Sp(iTxlayer)
        u  = rhs    
        du = 0d0
        
        ff = c(iTxlayer)
        
        dcdsig(:,iTxlayer) = ( (t*u*rr + ff*s)*dSmdsig(:,iTxlayer)*expmgh(iTxlayer) + & 
                           &  q*ff*dSpdsig(:,iTxlayer)*expmgh(iTxlayer) +   (p+q*rr)*u*dSpdsig(:,iTxlayer) ) / onemrmrp
       
       ! add in the heinous extra terms for derivative with respect to Txlayer conductivity  
        dcdsig(iTxlayer,iTxlayer) = dcdsig(iTxlayer,iTxlayer) + ( t*u*(dp + q*dr) + (p+q*rr)*t*du + &
                                  & isgnRp*(p+q*rr)*u*t*h(iTxlayer)*dgammadsig )/ onemrmrp
         
       ! dd/dsig:
        p  = exp(-gamma(iTxlayer)*depthTx )
        rr = -exp(-gamma(iTxlayer)*heightTx)
        dp = -p*depthTx*dgammadsig 
        dr = -rr*heightTx*dgammadsig 
 
        
        q  = Sp(iTxlayer)*expmgh(iTxlayer)
        
        s  = Sm(iTxlayer)*expmgh(iTxlayer)       
        
        t = Sm(iTxlayer)
    
        ff = d(iTxlayer)

        dddsig(:,iTxlayer) = ( (t*u*rr + ff*s)*dSpdsig(:,iTxlayer)*expmgh(iTxlayer) + & 
                           & q*ff*dSmdsig(:,iTxlayer)*expmgh(iTxlayer) +   (p+q*rr)*u*dSmdsig(:,iTxlayer) ) / onemrmrp
        
        ! add in the heinous extra terms for derivative with respect to Txlayer conductivity    
        dddsig(iTxlayer,iTxlayer) = dddsig(iTxlayer,iTxlayer) + ( t*u*(dp + q*dr) + (p+q*rr)*t*du + &
                                  & isgnRm*(p+q*rr)*u*t*h(iTxlayer)*dgammadsig )/ onemrmrp
        
        !
        ! Compute the source terms for the first step of the upward/downward continuations:
        !       
        srcpInv     = -srcp*heightTx*dgammadsig
        srcmInv     = -srcm*depthTx*dgammadsig
         
    endif
    
!
! Compute \Lambda_z potential coefficients for layers above the source layer:
!
    do i = iTxlayer-1,1,-1
        aytop =   c(i+1)*expmgh(i+1) + d(i+1)   + srcm  !   
        c(i)  = aytop / ( 1.d0 + Sm(i)*expmgh(i) )  
        d(i)  = c(i)*Sm(i)
        srcm = 0d0
    enddo
        
    if (linversion) then
    
        do i = iTxlayer-1,1,-1      
        ! base form for dcdsig:
            dcdsig(:,i) =  dddsig(:,i+1) + dcdsig(:,i+1)*expmgh(i+1)  - c(i)*dSmdsig(:,i)*expmgh(i)  
            
            ! add modified terms for dsig in layer i-1:
            dcdsig(i+1,i) = dcdsig(i+1,i) - c(i+1)*h(i+1)*expmgh(i+1)*(-ii*omega*mu0 / (2.d0 *gamma(i+1)))                  
     
            ! add modified terms for dsig in layer i-1 if i-1 is source layer:
            dcdsig(i+1,i)  = dcdsig(i+1,i) + srcmInv
            
            ! add terms for dsig in layer i  
            dcdsig(i,i) = dcdsig(i,i) + d(i)*h(i)*expmgh(i)*(-ii*omega*mu0 / (2.d0 *gamma(i)))  
            
            ! lastly divide all terms by this:
            dcdsig(:,i) = dcdsig(:,i) / (1.d0 + Sm(i)*expmgh(i) )
            
            ! Finally compute dbdsig:
            dddsig(:,i) = c(i)*dSmdsig(:,i) + Sm(i)*dcdsig(:,i)
            
            ! zero out source term for remaining layers
            srcmInv = 0d0         
            
        enddo        
        
    endif

!
! Compute \Lambda_z potential coefficients for layers below the source layer:
!   
    do i = iTxlayer+1,nlay1D
        aytop =  c(i-1) + d(i-1)*expmgh(i-1) + srcp   
        d(i)  = aytop / (1.d0 + Sp(i)*expmgh(i) )
        c(i)  = d(i)*Sp(i)
        srcp = 0d0
     enddo    
     
    if (linversion) then
    
        do i = iTxlayer+1,nlay1D
        ! base form for dddsig:
            dddsig(:,i) =  dcdsig(:,i-1) + dddsig(:,i-1)*expmgh(i-1)  - d(i)*dSpdsig(:,i)*expmgh(i)  
            
            ! add modified terms for dsig in layer i-1:
            dddsig(i-1,i) = dddsig(i-1,i) - d(i-1)*h(i-1)*expmgh(i-1)*(-ii*omega*mu0 / (2.d0 *gamma(i-1))) 
            
            ! add modified terms for dsig in layer i-1 if i-1 is source layer:
            dddsig(i-1,i)  = dddsig(i-1,i) + srcpInv
            
            ! add terms for dsig in layer i  
            dddsig(i,i) = dddsig(i,i) + c(i)*h(i)*expmgh(i)*(-ii*omega*mu0 / (2.d0 *gamma(i))) 
            
            ! lastly divide all terms by this:
            dddsig(:,i) = dddsig(:,i) / (1.d0 + Sp(i)*expmgh(i) )
            
            ! Finally compute dbdsig:
            dcdsig(:,i) = d(i)*dSpdsig(:,i) + Sp(i)*dddsig(:,i)
            
            ! zero out source term for remaining layers
            srcpInv = 0d0                
        enddo     
        
        
!        where (abs(dadsig) < 1d-80) dadsig = 0d0
!        where (abs(dbdsig) < 1d-80) dbdsig = 0d0
!        where (abs(dcdsig) < 1d-80) dcdsig = 0d0
!        where (abs(dddsig) < 1d-80) dddsig = 0d0    


    endif

        
!
! Zero out small values
!
!    where (abs(a) < 1d-80) a = 0d0
!    where (abs(b) < 1d-80) b = 0d0
!    where (abs(c) < 1d-80) c = 0d0
!    where (abs(d) < 1d-80) d = 0d0    

    end subroutine hed_coef_nlayer
    
!==============================================================================!
!============================================================== hed_pot_nlayer !
!==============================================================================!
    subroutine hed_pot_nlayer(lam)  
!
! Generates the vector potential Ay and Az, and their z derivatives
! for an N layered model with a horizontal electric dipole transmitter
! located in any layer.
! Note that the Az coefficients are for lambda_z, where 
! Az = d/dy (lambda_z), so there is an ik_y term left off that
! is accounted for later in the Hankel transform formulation.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Version 4.0  March 11, 2008.    Added df/dsig coefficients for inversion.
! Version 3.0  February 29, 2008. Added optional spline interpolation for speed. 
! Version 2.0  October 19, 2007.  Using new formula for Lambda_z potential
! Version 1.0  June 12, 2007.
!
    
    implicit none
!
! Input arguments:
!
    real(8), intent(in)  :: lam  ! Hankel transform variable
    
!
! Local variables:
!
    complex(8)  :: expcoef1, expcoef2, expsrc1, srcterm, dzsrcterm 
    complex(8)  :: aa,bb,cc,dd,gg  ! local copies of potential coefficients for iRxLayer of interest
    complex(8)  :: ae1,be2,ae1pbe2,ae1mbe2,ce1,de2,lam2,gg2,dgammadsig    
    
!
! Either compute potential coeffs or interpolate using precomputed values
!
    if (lUseSpline1D) then
!
! Get interpolated value of coefficients at current lambda:
!             
        call splint_hed(lam,aa,bb,cc,dd) 
                
!
! Compute gamma for Rx and Tx layers: (also computed in hed_coef)
!
        gamma(iRxlayer) = sqrt( lam**2 + k2(iRxlayer) )   
        gamma(iTxlayer) = sqrt( lam**2 + k2(iTxlayer) )  
        
    else
!
! Generate the coefficients for the current lam value:
!
        call hed_coef_nlayer(lam)
        
        aa = a(iRxlayer)
        bb = b(iRxlayer)
        cc = c(iRxlayer)
        dd = d(iRxlayer)
        
    endif

    gg = gamma(iRxlayer)
    
!
! Now compute the Ay and Az potentials and their derivatives    
!

!
! Layer up/down attenuation terms
!
    expcoef1 = exp(-gg*heightRx) ! upgoing term   
    expcoef2 = exp(-gg*depthRx)  ! downgoing term
!
! Source term if same layer as transmitter
!
    expsrc1  = exp(-gamma(iTxlayer)*dabs(z-zTx1D))
    

!
! Catch tiny numbers:
!
!    if ( abs(aa).lt.1d-150) aa = 0d0
!    if ( abs(bb).lt.1d-150) bb = 0d0
!    if ( abs(cc).lt.1d-150) cc = 0d0
!    if ( abs(dd).lt.1d-150) dd = 0d0  
!    if ( abs(expcoef1).lt.1d-150) expcoef1 = 0d0
!    if ( abs(expcoef2).lt.1d-150) expcoef2 = 0d0
!    if ( abs(expsrc1).lt.1d-150)  expsrc1  = 0d0

    ae1 = aa*expcoef1
    be2 = bb*expcoef2
    
    ae1pbe2 = aa*expcoef1 + bb*expcoef2
    ae1mbe2 = aa*expcoef1 - bb*expcoef2
    
    ce1 = cc*expcoef1
    de2 = dd*expcoef2
    
    lam2 = lam*lam
    gg2  = gg*gg
    
    srcterm  = isgnsrc*mu0*sdm1D/2d0/gamma(iTxlayer)*expsrc1  !isgnsrc is 0|1 
    
    ay       = ae1pbe2 + srcterm 
  
    dzsrcterm  = isgndsrcdz*mu0*sdm1D/2d0*expsrc1        ! isgndsrcdz is -1|0|1 for derivative term
     
    daydz    = gg*(ae1mbe2)  + dzsrcterm
       
    az       = ce1 + de2 - gg*( ae1mbe2)/lam2
    
    dazdz    = gg*(ce1 - de2) - ae1pbe2*gg2/lam2  
    
    dazdz2   = gg2*az
    
    aypdazdz = gg*(ce1 - de2) - ae1pbe2*k2(iRxlayer)/lam2 + srcterm 
    
    daydzpdazdz2 = gg2*(ce1 + de2) - gg*k2(iRxlayer)*ae1mbe2/lam2   + dzsrcterm
    
    csigsite  = csig(iRxlayer)

    if (linversion) then
    
  ! Form the common components:
       
        dgammadsig = -ii*omega*mu0 / (2.d0*gg)  
            
  ! Derivatives for a terms:        
        
        dae1dsig = dadsig(:,iRxlayer)*expcoef1
        ! add on terms for derivative in iRxlayer
        dae1dsig(iRxlayer) = dae1dsig(iRxlayer) - ae1*heightRx*dgammadsig ! note that ae1 is zero in bottom layer, 
                                                                          ! so no problem with heightRx term     
        
   ! Derivatives for b terms:           
        
        dbe2dsig = dbdsig(:,iRxlayer)*expcoef2
        ! add on terms for derivative in iRxlayer
        dbe2dsig(iRxlayer) = dbe2dsig(iRxlayer) - be2*depthRx*dgammadsig  ! note that be2 is zero in bottom layer, 
                                                                          ! so no problem with depthRx term
                

     ! dAydsig:
        
        dAydsig   = dae1dsig + dbe2dsig
        
        ! if iRxlayer = iTxlayer, we need to add source term derivatives:
        srcterm  = -isgnsrc*mu0*sdm1D/2d0*expsrc1*(1.d0/gg2 + dabs(z-zTx1D)/gg)*dgammadsig  !isgnsrc is 0|1 
        dAydsig(iRxlayer) = dAydsig(iRxlayer) + srcterm
        
            
     ! dAydsigdz:
        
        dAydsigdz =  gg*(dae1dsig - dbe2dsig) 
        
        ! Add on terms for derivative in iRxlayer:
        dAydsigdz(iRxlayer) = dAydsigdz(iRxlayer) + dgammadsig*ae1mbe2 
        
        ! if iRxlayer = iTxlayer, we need to add source term derivatives (this is harsh):
        ! note that we need source term for both sides of gg*(dae1dsig - dbe2dsig) 
        srcterm  = -isgndsrcdz*mu0*sdm1D/2d0*expsrc1 *dgammadsig * dabs(z-zTx1D)    
        
        ! isgndsrcdz is -1|0|1 for derivative term
        dAydsigdz(iRxlayer) = dAydsigdz(iRxlayer) + srcterm
        
 
    ! Derivatives for c terms:      
  
        dce1dsig = dcdsig(:,iRxlayer)*expcoef1
        ! add on terms for derivative in iRxlayer
        dce1dsig(iRxlayer) = dce1dsig(iRxlayer) - ce1*heightRx*dgammadsig           
        
    ! Derivatives for d terms:          
        
        dde2dsig = dddsig(:,iRxlayer)*expcoef2
        ! add on terms for derivative in iRxlayer
        dde2dsig(iRxlayer) = dde2dsig(iRxlayer) - de2*depthRx*dgammadsig  
            
 
     ! dAzdsig:
    
        dAzdsig = dce1dsig + dde2dsig - gg/lam2 * (dae1dsig - dbe2dsig )
        ! add on terms for derivative in iRxlayer.  Note that a,b,c,d derivatives have this already, but gg term doesn't
        dAzdsig(iRxlayer) =  dAzdsig(iRxlayer) - dgammadsig/lam2*(ae1mbe2)  
        ! no source terms in Az
        
    ! dAzdsigdz:    
        dAzdsigdz = gg*(dce1dsig - dde2dsig) - gg2/lam2 * (dae1dsig + dbe2dsig )
        ! add on terms for derivative in iRxlayer.  Note that a,b,c,d derivatives have this already, but gg terms don't
        dAzdsigdz(iRxlayer) = dAzdsigdz(iRxlayer) + dgammadsig*(ce1 - de2) - ae1pbe2*(-ii*omega*mu0 / lam2)   
        
    ! Okay,  here we have the final beast to finish off this arduous journey into derivative hell:
    ! dAzdsigdz2:  
        dAzdsigdz2 = gg2*dAzdsig
        dAzdsigdz2(iRxlayer) = dAzdsigdz2(iRxlayer) +  2*gg*dgammadsig*az

    endif
    
    end subroutine hed_pot_nlayer
    
!==============================================================================!
!============================================================= ved_coef_nlayer !
!==============================================================================!
    subroutine ved_coef_nlayer(lam)
!
! Generates the vector potential coefficients for Az
! for an N layered model with a vertical electric dipole transmitter 
! located in any layer.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Version 2.1  May 2008.          Modified df/dsig recursions for speedup.
! Version 2.0  March 12, 2008.    Added df/dsig recursions for inversion.
! Version 1.2  December 14, 2007. Modified for exp(-iwt) time dependence.
! Version 1.1  Sept, 2007.        Modifications for speed-up and stability.
! Version 1.0  June 14, 2007.
!

    
    implicit none
!
! Input arguments:
!
    real(8), intent(in) :: lam    ! Hankel transform variable
! 
! Local variables:
! 
    integer     :: i            
    complex(8)  :: sgmosgp, sjexp, onemrmrp, rmrp, srcp,srcm, rhs
 
   ! Derivative variables:
    complex(8)  :: dgammadsig,dsdsig,p,dp,rr,dr,q,s,t,u,du,ff,srcpInv,srcmInv,db
    integer     :: isgnRm, isgnRp
    
!
! Compute gamma for each layer:
!
    gamma(1:nlay1D) = sqrt( lam**2 + k2(1:nlay1D) )

!
! Compute attenuation factor for each layer
!
    expmgh = exp(-gamma*h)
    
    where (abs(expmgh) <1d-80)  expmgh = 0d0    ! avoid tiny numbers
!
! Initialize some arrays
!
    Sm = 0d0
    Sp = 0d0
    c  = 0d0
    d  = 0d0   
    
    if (linversion) then
    
        dSmdsig = 0d0 !  note that for bounding layers these terms start out as zero
        dSpdsig = 0d0
        
        dcdsig = 0d0
        dddsig = 0d0
        
    endif    
    
!
! Compute the recursion coefficients for the Az potential
! for the layers above the transmitter:
!
    do i = 2,iTxlayer
    
      Sm(i-1) = Sm(i-1)*expmgh(i-1) ! note that I post apply 
      ! this so that for the iTxlayer the exponent is omitted and we then 
      ! don't need the positive exponent in the formula for c_i and d_i
      
      sgmosgp = ( gamma(i)*csig(i-1) - gamma(i-1)*csig(i) ) / ( gamma(i)*csig(i-1) + gamma(i-1)*csig(i) )
      sjexp   = Sm(i-1)*expmgh(i-1) 
      Sm(i)   =  (sgmosgp + sjexp)  / ( 1.d0 + sgmosgp*sjexp)
      
      if (linversion) then
      
            dgammadsig      = -ii*omega*mu0 / (2.d0 *gamma(i))
            dsdsig          = (dgammadsig*csig(i-1)*(1.d0 - sgmosgp)  - &
                            & gamma(i-1)*(1.d0+ sgmosgp) )/( gamma(i)*csig(i-1) + gamma(i-1)*csig(i) )
            
            ! Post-apply expmgh(i-1) to  all recursion coefficients (done to avoid exp(+) term in final formula):
            dSmdsig(1:i-1,i-1) = dSmdsig(1:i-1,i-1)*expmgh(i-1) ! note that rows are now perturbed sigma, 
                                                                ! columns are now the rx location
            
            ! Initialize dSm(i)/dsig(i):        
            dSmdsig(i,i)    = dsdsig*(1.d0 - Sm(i)*sjexp) / (1.d0 + sgmosgp*sjexp) - (1- int(i/nlay1d) )*Sm(i)*h(i) *dgammadsig
            ! the (int) term deals with Tx in bottom layer
            
            ! Next step is to compute dSm(:)/dsig(i-1) for previous layers
            dgammadsig  = -ii*omega*mu0 / (2.d0 *gamma(i-1))
            
            dsdsig      = (  gamma(i)*(1.d0 - sgmosgp) - dgammadsig*csig(i)*(1.d0 + sgmosgp) ) / &
                        &  ( gamma(i)*csig(i-1) + gamma(i-1)*csig(i) )   ! actually ds(i)dsig(i-1)
            
            db          = dSmdsig(i-1,i-1)* expmgh(i-1) - sjexp*h(i-1)*dgammadsig
            
            dSmdsig(i-1,i) = ( dsdsig*(1.d0  - Sm(i)*sjexp)  + db*(1.d0  - Sm(i)*sgmosgp) ) / (1.d0 + sgmosgp*sjexp)    

            ! Apply chain rule to propagate deriv for layers already intialized on previous iterations:
            dSmdsig(1:i-2,i) = dSmdsig(1:i-2,i-1) * ( expmgh(i-1)*(1.d0 - Sm(i)*sgmosgp) )/(1.d0 + sgmosgp*sjexp)
            
            ! At this point  dSmdsig(iTxlayer,1:iTxlayer) has an expmgh(iTxlayer) factor omitted
            
      endif
      
    enddo
!
! Compute the recursion coefficients for the Az potential
! for the layers below the transmitter:
!
    do i = nlay1D-1,iTxlayer,-1
    
        Sp(i+1) = Sp(i+1)*expmgh(i+1)  ! note that I post apply 
        ! this so that for the iTxlayer the exponent is omitted and we then
        ! don't need the positive exponent in the formula for c_i and d_i       
        sgmosgp =  ( gamma(i)*csig(i+1) - gamma(i+1)*csig(i) ) / ( gamma(i)*csig(i+1) + gamma(i+1)*csig(i) )   
        sjexp  = Sp(i+1)*expmgh(i+1) ! since Sp(nlay1D) is zero, no need to worry about height of bottom layer
        Sp(i)  =  (sgmosgp + sjexp) / ( 1.d0 + sgmosgp*sjexp)
        
        if (linversion) then

            dgammadsig = -ii*omega*mu0 / (2.d0 *gamma(i))
            dsdsig     = (dgammadsig*csig(i+1)*(1.d0 - sgmosgp) - (gamma(i+1)*(1.d0+ sgmosgp)) ) / &
                       &  ( gamma(i)*csig(i+1) + gamma(i+1)*csig(i) )

            ! Post-apply expmgh(i+1) to  all recursion coefficients (done to avoid exp(+) term in final formula):       
            dSpdsig(i+1:nlay1d,i+1)     = dSpdsig(i+1:nlay1d,i+1)*expmgh(i+1)   
            
            ! Initialize dSm(i)/dsig(i):        
            dSpdsig(i,i) = dsdsig* (1.d0 - Sp(i)*sjexp) / (1.d0 + sgmosgp*sjexp) - & 
                         &  ( 1- int( (nlay1d-i+1)/nlay1d) )*Sp(i)*h(i) *dgammadsig ! the (int) term deals with Tx in top layer
            
            ! Next step is to compute dSm(:)/dsig(i+1) for previous layers
            dgammadsig  = -ii*omega*mu0 / (2.d0 *gamma(i+1))
            dsdsig      = ( gamma(i)*(1.d0 - sgmosgp) - dgammadsig*csig(i)*(1.d0 + sgmosgp)  ) / &
                        & ( gamma(i)*csig(i+1) + gamma(i+1)*csig(i) )   ! actually ds(i)dsig(i+1)
                        
            db          = dSpdsig(i+1,i+1)* expmgh(i+1) - sjexp*h(i+1)*dgammadsig

            dSpdsig(i+1,i) = ( dsdsig*( 1.d0 - Sp(i)*sjexp)  + db*( 1.d0 - Sp(i)*sgmosgp) ) / (1.d0 + sgmosgp*sjexp)    
             
            ! Apply chain rule to propagate deriv for layers already intialized on previous iterations:
            dSpdsig(i+2:nlay1d,i) = dSpdsig(i+2:nlay1d,i+1) * ( expmgh(i+1)*(1.d0 - Sp(i)*sgmosgp) )/(1.d0 + sgmosgp*sjexp)
            
            ! At this point dSpdsig(iTxlayer,iTxlayer:nlay1d) has an expmgh(iTxlayer) factor omitted
            
        endif
      
    enddo   
  
!
! Compute Az potential coefficients in the source layer:
!   
    rmrp     = Sm(iTxlayer)*Sp(iTxlayer)
    onemrmrp = 1.d0 - rmrp*expmgh(iTxlayer)*expmgh(iTxlayer) ! added exp term since I omitted it earlier
    rhs      = mu0*sdm1D/ (2.d0*gamma(iTxlayer))
    srcp     = exp(-gamma(iTxlayer)*heightTx)*rhs  
    srcm     = exp(-gamma(iTxlayer)*depthTx )*rhs  

    c(iTxlayer) = ( rmrp*srcm*expmgh(iTxlayer)  + Sp(iTxlayer)*srcp  ) / onemrmrp    
    d(iTxlayer) = ( rmrp*srcp*expmgh(iTxlayer)  + Sm(iTxlayer)*srcm  ) / onemrmrp    
    
    if (linversion) then
    
        ! Make sure nothing is tiny:
!        where (abs(dSmdsig) <1d-80) dSmdsig = 0d0
!        where (abs(dSpdsig) <1d-80) dSpdsig = 0d0
   
      !
      ! Catch to omit recursion derivatives if Tx in upper or lowermost layers
      !
        if (iTxlayer == nlay1D) then
            isgnRm = 0 
            isgnRp = 0          
        elseif (iTxlayer == 1)  then
            isgnRm = 0
            isgnRp = 0  
        else
            isgnRm = 1
            isgnRp = 1
        endif
        
        dgammadsig = -ii*omega*mu0 / (2.d0 *gamma(iTxlayer)) 
    
       ! dc/dsig:
        p  = exp(-gamma(iTxlayer)*heightTx)
        rr = exp(-gamma(iTxlayer)*depthTx )
        dp = -p*heightTx*dgammadsig 
        dr = -rr*depthTx*dgammadsig 
        
        q  = Sm(iTxlayer)*expmgh(iTxlayer)
        
        s  = Sp(iTxlayer)*expmgh(iTxlayer)
        
        t = Sp(iTxlayer)

        u  = rhs    
        du = -dgammadsig*rhs/gamma(iTxlayer)  ! term only for iTxlayer deriv
        
        ff = c(iTxlayer)
        
        dcdsig(:,iTxlayer) = ( (t*u*rr + ff*s)*dSmdsig(:,iTxlayer)*expmgh(iTxlayer) + & 
                           & q*ff*dSpdsig(:,iTxlayer)*expmgh(iTxlayer) +   (p+q*rr)*u*dSpdsig(:,iTxlayer) ) / onemrmrp
       
       ! add in the heinous extra terms for derivative with respect to Txlayer conductivity  
        dcdsig(iTxlayer,iTxlayer) = dcdsig(iTxlayer,iTxlayer) + ( t*u*(dp + q*dr) + (p+q*rr)*t*du + & 
                                  & isgnRp*(p+q*rr)*u*t*h(iTxlayer)*dgammadsig )/ onemrmrp
         
       ! dd/dsig:
        p  = exp(-gamma(iTxlayer)*depthTx )
        rr = exp(-gamma(iTxlayer)*heightTx)
        dp = -p*depthTx*dgammadsig 
        dr = -rr*heightTx*dgammadsig 
        
        q  = Sp(iTxlayer)*expmgh(iTxlayer)
        
        s  = Sm(iTxlayer)*expmgh(iTxlayer)       
        
        t = Sm(iTxlayer)
        
        ff = d(iTxlayer)

        dddsig(:,iTxlayer) = ( (t*u*rr + ff*s)*dSpdsig(:,iTxlayer)*expmgh(iTxlayer) + q*ff*dSmdsig(:,iTxlayer)*expmgh(iTxlayer) +  &
                           &  (p+q*rr)*u*dSmdsig(:,iTxlayer) ) / onemrmrp
        
        ! add in the heinous extra terms for derivative with respect to Txlayer conductivity    
        dddsig(iTxlayer,iTxlayer) = dddsig(iTxlayer,iTxlayer) + ( t*u*(dp + q*dr) + (p+q*rr)*t*du + &
                                  & isgnRm*(p+q*rr)*u*t*h(iTxlayer)*dgammadsig )/ onemrmrp
        
        !
        ! Compute the source terms for the first step of the upward/downward continuations:
        !       
        srcpInv     = -srcp/gamma(iTxlayer)*(gamma(iTxlayer)*heightTx + 1.d0) * dgammadsig
        srcmInv     = -srcm/gamma(iTxlayer)*(gamma(iTxlayer)*depthTx  + 1.d0) * dgammadsig
            
    endif
    
!
! Compute Az potential coefficients for layers above the source layer:
!
    do i = iTxlayer-1,1,-1
        rhs = (c(i+1)*expmgh(i+1) + d(i+1) + srcm ) 
        c(i) = rhs / ( 1.d0 + Sm(i)*expmgh(i) )
        d(i) = c(i)*Sm(i)
        srcm = 0d0 ! primary source term only for first layer above source layer
    enddo

    if (linversion) then
    
        do i = iTxlayer-1,1,-1      
        ! base form for dcdsig:
            dcdsig(:,i) =  dddsig(:,i+1) + dcdsig(:,i+1)*expmgh(i+1)  - c(i)*dSmdsig(:,i)*expmgh(i)  
            
            ! add modified terms for dsig in layer i-1:
            dcdsig(i+1,i) = dcdsig(i+1,i) - c(i+1)*h(i+1)*expmgh(i+1)*(-ii*omega*mu0 / (2.d0 *gamma(i+1)))                  
     
            ! add modified terms for dsig in layer i-1 if i-1 is source layer:
            dcdsig(i+1,i)  = dcdsig(i+1,i) + srcmInv
            
            ! add terms for dsig in layer i  
            dcdsig(i,i) = dcdsig(i,i) + d(i)*h(i)*expmgh(i)*(-ii*omega*mu0 / (2.d0 *gamma(i)))  
            
            ! lastly divide all terms by this:
            dcdsig(:,i) = dcdsig(:,i) / (1.d0 + Sm(i)*expmgh(i) )
            
            ! Finally compute dbdsig:
            dddsig(:,i) = c(i)*dSmdsig(:,i) + Sm(i)*dcdsig(:,i)
            
            ! zero out source term for remaining layers
            srcmInv = 0d0         
            
        enddo        
        
    endif


!
! Compute Az potential coefficients for layers below the source layer:
!   
    do i = iTxlayer+1,nlay1D
        ! Store value of Ay potential at top of layer i
        rhs = (c(i-1) + d(i-1)*expmgh(i-1) + srcp ) 
        d(i) = rhs / (1+Sp(i)*expmgh(i) )
        c(i) = d(i)*Sp(i)
        srcp = 0d0 ! primary source term only for first layer below source layer
    enddo
    
    if (linversion) then
    
        do i = iTxlayer+1,nlay1D
        ! base form for dddsig:
            dddsig(:,i) =  dcdsig(:,i-1) + dddsig(:,i-1)*expmgh(i-1)  - d(i)*dSpdsig(:,i)*expmgh(i)  
            
            ! add modified terms for dsig in layer i-1:
            dddsig(i-1,i) = dddsig(i-1,i) - d(i-1)*h(i-1)*expmgh(i-1)*(-ii*omega*mu0 / (2.d0 *gamma(i-1))) 
            
            ! add modified terms for dsig in layer i-1 if i-1 is source layer:
            dddsig(i-1,i)  = dddsig(i-1,i) + srcpInv
            
            ! add terms for dsig in layer i  
            dddsig(i,i) = dddsig(i,i) + c(i)*h(i)*expmgh(i)*(-ii*omega*mu0 / (2.d0 *gamma(i))) 
            
            ! lastly divide all terms by this:
            dddsig(:,i) = dddsig(:,i) / (1.d0 + Sp(i)*expmgh(i) )
            
            ! Finally compute dbdsig:
            dcdsig(:,i) = d(i)*dSpdsig(:,i) + Sp(i)*dddsig(:,i)
            
            ! zero out source term for remaining layers
            srcpInv = 0d0                
        enddo     
        
        
!        where (abs(dadsig) < 1d-80) dadsig = 0d0
!        where (abs(dbdsig) < 1d-80) dbdsig = 0d0
!        where (abs(dcdsig) < 1d-80) dcdsig = 0d0
!        where (abs(dddsig) < 1d-80) dddsig = 0d0    

    endif

!
! Zero out small values
!
!    where (abs(c) < 1d-80) c = 0d0
!    where (abs(d) < 1d-80) d = 0d0    
    
    
    end subroutine ved_coef_nlayer

!==============================================================================!
!============================================================== ved_pot_nlayer !
!==============================================================================!
    subroutine ved_pot_nlayer(lam)  
!
! Generates the vector potential Az, and z derivatives
! for an N layered model with a vertical electric dipole transmitter
! located in any layer.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Version 3.0  March 12, 2008.    Added df/dsig coefficients for inversion.
! Version 2.0  February 29, 2008. Added optional spline interpolation for speed. 
! Version 1.0  June 14, 2007.
!
    
    implicit none
!
! Input arguments:
!
    real(8), intent(in) :: lam  ! Hankel transform variable
!
! Local variables:
!
    complex(8)  :: expcoef1, expcoef2, expsrc1, srcterm
    complex(8)  :: cc,dd,gg  ! local copies of potential coefficients for iRxLayer of interest
    complex(8)  :: ce1,de2,dgammadsig 
    
!
! First generate the coefficients for the current lam value:
!
!
! Either compute potential coeffs or interpolate using precomputed values
!
    if (lUseSpline1D) then
!
! Get interpolated value of coefficients at current lambda:
!             
        call splint_ved(lam,cc,dd) 
    
!
! Compute gamma for Rx and Tx layers: (also computed in hed_coef)
!
        gamma(iRxlayer) = sqrt( lam**2 + k2(iRxlayer) )   
        gamma(iTxlayer) = sqrt( lam**2 + k2(iTxlayer) )  
        
    else
!
! Generate the coefficients for the current lam value:
!
        call ved_coef_nlayer(lam)

        cc = c(iRxlayer)
        dd = d(iRxlayer)
        
    endif

    gg = gamma(iRxlayer)

!
! Now compute the Az potential and z derivatives    
!

!
! Layer up/down attenuation terms:
!
    expcoef1 = exp(-gg*heightRx)  ! upgoing term 
    expcoef2 = exp(-gg*depthRx)  ! downgoing term
!
! Source term if same layer as transmitter
!   
    expsrc1  = exp(-gamma(iTxlayer)*dabs(z-zTx1D))
    
!    if ( abs(cc).lt.1d-150) cc=0d0
!    if ( abs(dd).lt.1d-150) dd=0d0    
!    if ( abs(expcoef1).lt.1d-150) expcoef1=0d0
!    if ( abs(expcoef2).lt.1d-150) expcoef2=0d0
!    if ( abs(expsrc1).lt.1d-150) expsrc1=0d0
    
    srcterm  = isgnsrc*mu0*sdm1D/2d0/gamma(iTxlayer)*expsrc1  !isgnsrc is 0|1 
    
    ce1 = cc*expcoef1
    
    de2 = dd*expcoef2
    
    az       = ce1 + de2 + srcterm 
    
    srcterm  = isgndsrcdz*mu0*sdm1D/2d0*expsrc1        ! isgndsrcdz is -1|0|1 for derivative term   
    
    dazdz    = gg*(ce1 - de2)  + srcterm
    
    srcterm  = isgnsrc*mu0*sdm1D/2d0*expsrc1*(gamma(iTxlayer)*abs(isgndsrcdz)) ! second derivative of source term
    
    dazdz2   = gg*gg*(ce1 + de2 ) + srcterm
  
    csigsite = csig(iRxlayer)
    
    if (linversion) then
    
  ! Form the common components:   

        dgammadsig = -ii*omega*mu0 / (2.d0*gg)  
                
    ! Derivatives for c terms:      
  
        dce1dsig = dcdsig(:,iRxlayer)*expcoef1
        ! add on terms for derivative in iRxlayer
        dce1dsig(iRxlayer) = dce1dsig(iRxlayer) - ce1*heightRx*dgammadsig           
        
    ! Derivatives for d terms:          
        
        dde2dsig = dddsig(:,iRxlayer)*expcoef2
        ! add on terms for derivative in iRxlayer
        dde2dsig(iRxlayer) = dde2dsig(iRxlayer) - de2*depthRx*dgammadsig  
                        
     ! dAzdsig:
    
        dAzdsig = dce1dsig + dde2dsig 
        ! if iRxlayer = iTxlayer, we need to add source term derivatives:
        srcterm  = -isgnsrc*mu0*sdm1D/2d0*expsrc1*(1.d0/(gg*gg) + dabs(z-zTx1D)/gg)*dgammadsig  !isgnsrc is 0|1 
        dAzdsig(iRxlayer) = dAzdsig(iRxlayer) + srcterm
        
     ! dAzdsigdz:

        dAzdsigdz = gg*(dce1dsig - dde2dsig)
        ! add on terms for d(gg)/dsig in iRxlayer
        dAzdsigdz(iRxlayer) = dAzdsigdz(iRxlayer) + dgammadsig*(ce1-de2)
        ! if iRxlayer = iTxlayer, we need to add source term derivatives:
        srcterm = -isgndsrcdz*mu0*sdm1D/2d0*expsrc1 *dgammadsig * dabs(z-zTx1D)    
        dAzdsigdz(iRxlayer) = dAzdsigdz(iRxlayer) + srcterm
        
      ! dAzdsigdz2:

        dAzdsigdz2 = gg*gg*(dce1dsig + dde2dsig)
        ! add on terms for d(gg)/dsig in iRxlayer
        dAzdsigdz2(iRxlayer) = dAzdsigdz2(iRxlayer) + 2*gg*dgammadsig*(ce1+de2)
        ! if iRxlayer = iTxlayer, we need to add source term derivatives:
         ! second derivative of source term
        srcterm = isgnsrc*mu0*sdm1D/2d0*expsrc1*(abs(isgndsrcdz) - dabs(z-zTx1D)*gamma(iTxlayer)*abs(isgndsrcdz) )*dgammadsig   
       

        dAzdsigdz2(iRxlayer) = dAzdsigdz2(iRxlayer) + srcterm
      
    
    endif
        
    end subroutine ved_pot_nlayer
            
!==============================================================================!
!========================================================= PreComputePotCoeffs ! 
!==============================================================================! 
    subroutine PreComputePotCoeffs
!
! Subroutine to precompute potential coefficients across lambda range so that many sites can be
! computed rapidly by spline interpolation in lambda of the precomputed coefficients.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Version 1.2       November 5, 2009        Fixed a bug J. Gunning spotted in array index for 2nd order interpolation coeffs.
!                                           Added code that only saves interpolation coefficients for layers that have receivers.
!                                           Huge memory savings! 12x less for Canonical model with 75 unknown layers.
!
! Version 1.1       September 17,2008       Fixed allocation bug for lLayerHasRx that 
!                                           could result in lots of extra spline interpolations!
!                                           
! Version 1.0       Spring 2008
!
    
    implicit none
!
! Local variables
!
    integer     :: i,j, icnt, nlayersInterp
    real(8)     :: r_min,r_max,min_lam,max_lam, min_base,max_base, dloglam
    
    real(8)     :: lamtest, lam_hed_min, lam_hed_max, lam_ved_min, lam_ved_max, pot_coef_tol

    integer, dimension(:), allocatable :: ilayersInterp

    logical :: linversionIn
    
!-------------------------------
! 1. Define range of lam to use:
!-------------------------------

!
! 1a. Loop over sites and get min/max ranges to Tx
!
    r_min = 1d99
    r_max = 0d0
    
    do i = 1,n1D           
    
        dx         = x1D(i) - xTx1D
        dy         = y1D(i) - yTx1D
        r          = sqrt( dx**2 + dy**2)  ! Horizontal range to site
        if (r<1d0) r = 1d0
        if (r>r_max) r_max = r
        if (r<r_min) r_min = r
        
    enddo ! loop over sites
    

!
! Create mapping from input layer number to interp array layers, so that we only store interp coeffs for layers with Rx's:
!
    allocate ( iRxLayerInterp(nlay1D) )
    iRxLayerInterp = 0
    do i = 1,n1D  ! Loop over receivers
        iRxlayer = 1            ! Get layer current Rx resides in:
        do j= 2,nlay1D
          if (z1D(i).gt.zlay1D(j))  then   ! Sites on boundaries use the top layer 
            iRxlayer = j
          endif
        enddo
        iRxLayerInterp(iRxlayer) = 1
    enddo
    ! now modify iRxLayerInterp to have interp array column index
    icnt = 0
    do i = 1,nlay1D
        if ( iRxLayerInterp(i) /= 0) then
            icnt = icnt + 1
            iRxLayerInterp(i) = icnt
        endif
    enddo
    ! store number of interp layers:
    nlayersInterp = icnt
    
    allocate( ilayersInterp(nlayersInterp) )  ! this small array is used to extract the required layers for interpolation
    icnt = 0
    do i = 1,nlay1D
        if ( iRxLayerInterp(i) /= 0) then
            icnt = icnt  + 1
            ilayersInterp(icnt) = i      
        endif
    enddo 
!
! 1b.  Set min and max lam_interp values
!
! Filters already defined from call to initialize_dipole1d
!
    select case (trim(outputdomain1D))
        
    case ('spatial')
        min_base = base(1)
        max_base = base(ndhtfc)     
        dloglam = ( log10(base(2)) - log10(base(1)) )  ! spacing is same as filter
        
        min_lam = min_base/r_max 
        max_lam = max_base/r_min 
    
       
    case ('kx')  
        min_base = basecsfc(1)
        max_base = basecsfc(ncsfc)
        dloglam = ( log10(basecsfc(2)) - log10(basecsfc(1)) ) ! spacing is same as filter
        

        min_lam = min_base/r_max 
        max_lam = max_base/r_min 
        
    end select
    
        
!
! Now test drive through every decade of lambda range and identify where potential coefficients are of reasonable size:
!
     pot_coef_tol = 1d-30

    linversionIn = linversion
    linversion = .false.
    
    if (lhed)  then
!
! Test for lam_hed_min:
!
        lamtest = min_lam
        lam_hed_min = 0d0
        do while ( (lam_hed_min==0).and.(lamtest <= max_lam) )

            ! compute potential coefficients    
            call hed_coef_nlayer(lamtest)
            
            ! Test all layers with Rxs:
            do i = 1,nlay1d
                if  ( iRxLayerInterp(i) /= 0)  then
                   if ( (abs(a(i)) > pot_coef_tol) .or. (abs(b(i)) > pot_coef_tol) .or. &
                      & (abs(c(i)) > pot_coef_tol) .or. (abs(d(i)) > pot_coef_tol) ) then  
                        lam_hed_min = lamtest   
                    endif
                endif
            enddo
            
            ! Increment lamtest up a decade:
            lamtest = 10**(log10(lamtest)+.25)
        enddo
!
! Test for lam_hed_max:
!       
        lamtest = max_lam
        
        lam_hed_max = 0d0
        do while ((lam_hed_max==0).and.(lamtest >= min_lam))

            ! compute potential coefficients    
            call hed_coef_nlayer(lamtest)
            
            ! Test all layers with Rxs:
            do i = 1,nlay1d
                if  ( iRxLayerInterp(i) /= 0)  then
                   if ( (abs(a(i)) > pot_coef_tol) .or. (abs(b(i)) > pot_coef_tol) .or. &
                      & (abs(c(i)) > pot_coef_tol) .or. (abs(d(i)) > pot_coef_tol) ) then             
                        lam_hed_max = lamtest   
                    endif
                endif
            enddo
            
            ! Increment lamtest down a decade:
            lamtest = 10**(log10(lamtest)-.25)
        enddo   
 
    endif


    if (lved)  then
!
! Test for lam_ved_min:
!
        lamtest = min_lam
        lam_ved_min = 0d0
        do while ( (lam_ved_min==0).and.(lamtest <= max_lam) )

            ! compute potential coefficients    
            call ved_coef_nlayer(lamtest)
            
            ! Test all layers with Rxs:
            do i = 1,nlay1d
                if  ( iRxLayerInterp(i) /= 0) then
                   if ( (abs(c(i)) > pot_coef_tol) .or. (abs(d(i)) > pot_coef_tol) ) then  
                        lam_ved_min = lamtest   
                    endif
                endif
            enddo
            
            ! Increment lamtest up a decade:
            lamtest = 10**(log10(lamtest)+.25)
        enddo
!
! Test for lam_ved_max:
!       
        lamtest = max_lam
        lam_ved_max = 0d0
        do while  ((lam_ved_max==0).and.(lamtest >= min_lam) )
            ! compute potential coefficients    
            call ved_coef_nlayer(lamtest)
            
            ! Test all layers with Rxs:
            do i = 1,nlay1d
                if  ( iRxLayerInterp(i) /= 0)  then
                   if ( (abs(c(i)) > pot_coef_tol) .or. (abs(d(i)) > pot_coef_tol) ) then                 
                        lam_ved_max = lamtest   
                    endif
                endif
            enddo
            
            ! Increment lamtest down a decade:
            lamtest = 10**(log10(lamtest)-.25)
        enddo       

    endif
    
    linversion = linversionIn  ! turns this back to the original setting, had to force this off for lambda test drive efficiency
!
! Update min/max lam values based on the test drive results:
!
    if (lhed.and.lved) then
        min_lam = min(lam_hed_min,lam_ved_min)
        max_lam = min(lam_hed_max,lam_ved_max)  
    elseif (lhed) then
        min_lam = lam_hed_min
        max_lam = lam_hed_max       
    elseif (lved) then
        min_lam = lam_ved_min
        max_lam = lam_ved_max           
    endif

!
! 1c.  Set up lam_interp and allocate coefficient arrays
!
!   
! Compute total number if index points:
!
    if ( ( min_lam == 0) .and. (max_lam == 0) ) then
        nlam_interp = 1 ! kwk debug, fix me
        
    else
        nlam_interp = ceiling( (log10(max_lam)-log10(min_lam) )/dloglam) + 1
    endif

!    
! Allocating both hed and ved arrays although only one might be used (clean up later...)
!
    allocate( lam_interp(nlam_interp) )
    allocate(  a_hed_interp(nlam_interp,nlayersInterp),  b_hed_interp(nlam_interp,nlayersInterp) )
    allocate(  c_hed_interp(nlam_interp,nlayersInterp),  d_hed_interp(nlam_interp,nlayersInterp) )
    allocate( a2_hed_interp(nlam_interp,nlayersInterp), b2_hed_interp(nlam_interp,nlayersInterp) )
    allocate( c2_hed_interp(nlam_interp,nlayersInterp), d2_hed_interp(nlam_interp,nlayersInterp) )    

    allocate(  c_ved_interp(nlam_interp,nlayersInterp),  d_ved_interp(nlam_interp,nlayersInterp) )
    allocate( c2_ved_interp(nlam_interp,nlayersInterp), d2_ved_interp(nlam_interp,nlayersInterp) )
    
    
    if (linversion) then
        allocate(  dadsig_hed_interp(nlam_interp,nlay1D,nlayersInterp),  dbdsig_hed_interp(nlam_interp,nlay1D,nlayersInterp) )
        allocate(  dcdsig_hed_interp(nlam_interp,nlay1D,nlayersInterp),  dddsig_hed_interp(nlam_interp,nlay1D,nlayersInterp) )
        allocate( dadsig2_hed_interp(nlam_interp,nlay1D,nlayersInterp), dbdsig2_hed_interp(nlam_interp,nlay1D,nlayersInterp) )
        allocate( dcdsig2_hed_interp(nlam_interp,nlay1D,nlayersInterp), dddsig2_hed_interp(nlam_interp,nlay1D,nlayersInterp) )    
    
        allocate(  dcdsig_ved_interp(nlam_interp,nlay1D,nlayersInterp),  dddsig_ved_interp(nlam_interp,nlay1D,nlayersInterp) )
        allocate( dcdsig2_ved_interp(nlam_interp,nlay1D,nlayersInterp), dddsig2_ved_interp(nlam_interp,nlay1D,nlayersInterp) )     
        
    endif   

!           
! Create lam_interp values:
!
    lam_interp(1) = min_lam
    
    do i = 2,nlam_interp
        lam_interp(i) = lam_interp(i-1) * 10.d0**dloglam  
    enddo
    

    
!------------------------------------------------------------------------   
! 2.Compute coefficients and get spline coefficients of second derivative
!------------------------------------------------------------------------ 

! 
! 2a. Compute horizontal coefficients if needed
!
    if (lhed)  then
    
        do i = 1,nlam_interp    
        
            ! compute potential coefficients    
            call hed_coef_nlayer(lam_interp(i))

            
            ! Store them into the array 
            a_hed_interp(i,1:nlayersInterp) = a(ilayersInterp)
            b_hed_interp(i,1:nlayersInterp) = b(ilayersInterp)
            c_hed_interp(i,1:nlayersInterp) = c(ilayersInterp)
            d_hed_interp(i,1:nlayersInterp) = d(ilayersInterp)

            if (linversion) then
                dadsig_hed_interp(i,1:nlay1d,1:nlayersInterp) = dadsig(1:nlay1d,ilayersInterp)
                dbdsig_hed_interp(i,1:nlay1d,1:nlayersInterp) = dbdsig(1:nlay1d,ilayersInterp)
                dcdsig_hed_interp(i,1:nlay1d,1:nlayersInterp) = dcdsig(1:nlay1d,ilayersInterp)
                dddsig_hed_interp(i,1:nlay1d,1:nlayersInterp) = dddsig(1:nlay1d,ilayersInterp)
            endif
            
        enddo

        ! Compute the spline interpolation coefficients:
        do i = 1,nlayersInterp 
        
            call spline_i1d(lam_interp,a_hed_interp(:,i),nlam_interp,a2_hed_interp(:,i) )
            call spline_i1d(lam_interp,b_hed_interp(:,i),nlam_interp,b2_hed_interp(:,i) )
            call spline_i1d(lam_interp,c_hed_interp(:,i),nlam_interp,c2_hed_interp(:,i) )
            call spline_i1d(lam_interp,d_hed_interp(:,i),nlam_interp,d2_hed_interp(:,i) )

        !
        ! KWK Note: These take up most of the CPU time for each call to Dipole1D with linversion = true.
        ! This is the place to concentrate future optimization efforts.
        !
            if (linversion) then                
                do j = 1,nlay1d  ! interpolate sensitivity to each layer (i.e., the 2nd index):
                    call spline_i1d(lam_interp,dadsig_hed_interp(:,j,i),nlam_interp,dadsig2_hed_interp(:,j,i) )
                    call spline_i1d(lam_interp,dbdsig_hed_interp(:,j,i),nlam_interp,dbdsig2_hed_interp(:,j,i) )
                    call spline_i1d(lam_interp,dcdsig_hed_interp(:,j,i),nlam_interp,dcdsig2_hed_interp(:,j,i) )
                    call spline_i1d(lam_interp,dddsig_hed_interp(:,j,i),nlam_interp,dddsig2_hed_interp(:,j,i) )

                enddo
            endif

        enddo

    endif
! 
! 2b. Compute vertical coefficients if needed
!       
    if (lved)  then
    
        do i = 1,nlam_interp    
        
            ! compute potential coefficients    
            call ved_coef_nlayer(lam_interp(i))
            
            ! Store them into the array 
            c_ved_interp(i,1:nlayersInterp) = c(ilayersInterp)
            d_ved_interp(i,1:nlayersInterp) = d(ilayersInterp)
            
            if (linversion) then
                dcdsig_ved_interp(i,1:nlay1d,1:nlayersInterp) = dcdsig(1:nlay1d,ilayersInterp)
                dddsig_ved_interp(i,1:nlay1d,1:nlayersInterp) = dddsig(1:nlay1d,ilayersInterp)
            endif
                    
        enddo

        ! Compute the spline interpolating coefficients:
        
        do i = 1,nlayersInterp

            call spline_i1d(lam_interp,c_ved_interp(:,i),nlam_interp,c2_ved_interp(:,i) )
            call spline_i1d(lam_interp,d_ved_interp(:,i),nlam_interp,d2_ved_interp(:,i) )       

            if (linversion) then
                
                do j = 1,nlay1d  !   nlam, nlay, nlayersInterp
                    call spline_i1d(lam_interp,dcdsig_ved_interp(:,j,i),nlam_interp,dcdsig2_ved_interp(:,j,i) )
                    call spline_i1d(lam_interp,dddsig_ved_interp(:,j,i),nlam_interp,dddsig2_ved_interp(:,j,i) )
                enddo
            endif
            
        enddo

    endif
    
!
! Deallocate:
!
    deallocate( ilayersInterp  ) 
     
    end subroutine PreComputePotCoeffs
    
!==============================================================================!
!================================================================== spline_i1d ! 
!==============================================================================!            
     subroutine spline_i1d(x,y,n,y2)
!
! Subroutine to generate second derivative of spline interpolating function
! Modified from Numerical Recipes section 3.3.  Use this before calling splint.
!
! x (real, length n), y (complex, length n) are the precomputed values of 
! the function to be splined.  y2 (complex, length n) is the second derivative 
! of the interpolating function
!
    implicit none
    
    
    integer, intent(in)     :: n
    real(8), intent(in)     :: x(n)
    complex(8), intent(in)  :: y(n)  
    complex(8), intent(out) :: y2(n)    
    
    ! Local:
    
    integer i,k
    real(8) sig,x2(:)
    real(8) p,un,qn,u(:)
    real(8) yreal(:), yimag(:), y2imag(:), y2real(:)
    
    allocatable u,x2, yreal, yimag, y2imag, y2real
    allocate (u(n),x2(n),yreal(n), yimag(n), y2imag(n), y2real(n))
    
    x2 = log10(x) ! KWK use log10 x values for better accuracy

!
! Compute real parts:
!
    y2real = 0d0
    u      = 0d0
    yreal  = dble(y)
       
    do i=2,n-1
        sig = ( x2(i) - x2(i-1) ) / ( x2(i+1) - x2(i-1) )
        p = sig*y2real(i-1) + 2.0d0
        y2real(i) = (sig-1.0d0)/p
        u(i) = (6.d0*( (yreal(i+1) - yreal(i)) / ( x2(i+1) - x2(i) ) - (yreal(i)-yreal(i-1)) &
        &       /( x2(i) - x2(i-1) ) )/( x2(i+1) - x2(i-1) )-sig*u(i-1) )/p
    enddo
 
    qn= 0d0 
    un= 0d0 
    
    y2real(n)=( un - qn*u(n-1) ) / ( qn*y2real(n-1) + 1.d0)
 
    do  k=n-1,1,-1
        y2real(k) = y2real(k)*y2real(k+1) + u(k)
    enddo
    
!
! Compute imaginary parts:
!
    y2imag = 0d0
    u      = 0d0
    yimag  = aimag(y)
       
    do  i=2,n-1
        sig = ( x2(i) - x2(i-1) ) / ( x2(i+1) - x2(i-1) )
        p = sig*y2imag(i-1) + 2.0d0
        y2imag(i) = (sig-1.0d0)/p
        u(i) = (6.d0*( (yimag(i+1) - yimag(i)) / ( x2(i+1) - x2(i) ) - (yimag(i)-yimag(i-1)) &
        &       /( x2(i) - x2(i-1) ) )/( x2(i+1) - x2(i-1) )-sig*u(i-1) )/p
    enddo
 
    qn= 0d0 
    un= 0d0 
    
    y2imag(n)=( un - qn*u(n-1) ) / ( qn*y2imag(n-1) + 1.d0)
 
    do  k=n-1,1,-1
        y2imag(k) = y2imag(k)*y2imag(k+1) + u(k)
    enddo
        
    y2 = cmplx(y2real,y2imag)
    
    deallocate (u,x2,yreal,yimag,y2real,y2imag)

    return
    
    end subroutine spline_i1d
      
!==============================================================================!
!================================================================== splint_hed ! 
!==============================================================================!      
      subroutine splint_hed(x,aa,bb,cc,dd)
!
! Performs spline interpolation using lam_interp and input arrays of values.
! Uses cubic spline interpolation in log10(x)
!
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Version 1.2       November 5, 2009        Updated interp coeff arrays for new efficient memory arrays
!    
    implicit none
    
    complex(8), intent(out) :: aa,bb,cc,dd
    real(8), intent(in)     :: x
!
! Local variables
! 
    integer  klo, khi, istr,i, ilay
    real(8) ai, bi, h2, a3ma, b3mb
    
    aa = 0d0
    bb = 0d0
    cc = 0d0
    dd = 0d0
    
!    if (linversion) then
!        dadsig(iRxlayer,:) = 0d0
!        dbdsig(iRxlayer,:) = 0d0
!        dcdsig(iRxlayer,:) = 0d0
!        dddsig(iRxlayer,:) = 0d0
!    endif

!
! get the index pointer:
!
    istr = iptr
    do i = istr,nlam_interp-1
        if ( (lam_interp(i) <= x) .and. (lam_interp(i+1) >= x) ) then
            iptr = i
            klo = iptr
            khi = iptr+1
      
        !
        ! Compute weights:
        !
            h2 = log10(lam_interp(khi)) - log10(lam_interp(klo))        
            ai  = (log10(lam_interp(khi)) - log10(x))/h2
            bi  = (log10(x) - log10(lam_interp(klo)))/h2
            a3ma = (ai**3-ai)*(h2**2)/6.d0
            b3mb = (bi**3-bi)*(h2**2)/6.d0
        !
        ! Compute spline interpolation:
        !
            ilay = iRxLayerInterp(iRxlayer)
        
            aa =    ai *  a_hed_interp(klo,ilay)  +  bi *  a_hed_interp(khi,ilay) + &
         &       (a3ma * a2_hed_interp(klo,ilay) + b3mb * a2_hed_interp(khi,ilay)) 
            bb =    ai *  b_hed_interp(klo,ilay)  +  bi *  b_hed_interp(khi,ilay) + &
         &       (a3ma * b2_hed_interp(klo,ilay) + b3mb * b2_hed_interp(khi,ilay)) 
            cc =    ai *  c_hed_interp(klo,ilay)  +  bi *  c_hed_interp(khi,ilay) + &
         &       (a3ma * c2_hed_interp(klo,ilay) + b3mb * c2_hed_interp(khi,ilay)) 
            dd =    ai *  d_hed_interp(klo,ilay)  +  bi *  d_hed_interp(khi,ilay) + &
         &       (a3ma * d2_hed_interp(klo,ilay) + b3mb * d2_hed_interp(khi,ilay))              



                        
            if (linversion) then
            
                dadsig(1:nlay1d,iRxlayer) =    ai *  dadsig_hed_interp(klo,1:nlay1d,ilay)  +  & 
                 &                             bi *  dadsig_hed_interp(khi,1:nlay1d,ilay)  +  &
                 &       (a3ma * dadsig2_hed_interp(klo,1:nlay1d,ilay) + b3mb * dadsig2_hed_interp(khi,1:nlay1d,ilay)) 
                 
                dbdsig(1:nlay1d,iRxlayer) =    ai *  dbdsig_hed_interp(klo,1:nlay1d,ilay)  +  &
                 &                             bi *  dbdsig_hed_interp(khi,1:nlay1d,ilay)  +  &
                 &       (a3ma * dbdsig2_hed_interp(klo,1:nlay1d,ilay) + b3mb * dbdsig2_hed_interp(khi,1:nlay1d,ilay)) 
                
                dcdsig(1:nlay1d,iRxlayer) =    ai *  dcdsig_hed_interp(klo,1:nlay1d,ilay)  +  &
                 &                             bi *  dcdsig_hed_interp(khi,1:nlay1d,ilay)  +  &
                 &       (a3ma * dcdsig2_hed_interp(klo,1:nlay1d,ilay) + b3mb * dcdsig2_hed_interp(khi,1:nlay1d,ilay)) 
                    
                dddsig(1:nlay1d,iRxlayer) =    ai *  dddsig_hed_interp(klo,1:nlay1d,ilay)  +  &
                 &                             bi *  dddsig_hed_interp(khi,1:nlay1d,ilay)  +  &
                 &       (a3ma * dddsig2_hed_interp(klo,1:nlay1d,ilay) + b3mb * dddsig2_hed_interp(khi,1:nlay1d,ilay))     
    
                
            endif
            
            exit
            
        endif
        
    enddo
 
      
    end subroutine splint_hed

!==============================================================================!
!================================================================== splint_ved ! 
!==============================================================================!      
    subroutine splint_ved(x,cc,dd)
!
! Performs spline interpolation using lam_interp and input arrays of values.
! Uses cubic spline interpolation in log10(x)
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Version 1.2       November 5, 2009        Updated interp coeff arrays for new efficient memory arrays
!    
    implicit none
    
    complex(8), intent(out) :: cc,dd
    real(8), intent(in)     :: x
!
! Local variables
! 
    integer  klo, khi, istr, i, ilay
    real(8) ai, bi, h2, a3ma, b3mb
 
 
    cc = 0d0
    dd = 0d0
    
!    if (linversion) then
!        dcdsig(1:nlay1d,iRxlayer) = 0d0
!        dddsig(1:nlay1d,iRxlayer) = 0d0
!    endif
!
! get the index pointer:
!
    istr = iptr
    do i = istr,nlam_interp-1
        !
        ! Find bounding points:
        !   
        if ( (lam_interp(i) <= x) .and. (lam_interp(i+1) >= x) ) then
            iptr = i
            klo = iptr
            khi = iptr+1
        !
        ! Compute weights:
        !
            h2 = log10(lam_interp(khi)) - log10(lam_interp(klo))        
            ai  = (log10(lam_interp(khi)) - log10(x))/h2
            bi  = (log10(x) - log10(lam_interp(klo)))/h2
            a3ma = (ai**3-ai)*(h2**2)/6.d0
            b3mb = (bi**3-bi)*(h2**2)/6.d0
        !
        ! Compute spline interpolation:
        !
            ilay = iRxLayerInterp(iRxlayer)
                
            cc =    ai *  c_ved_interp(klo,ilay)  +   bi * c_ved_interp(khi,ilay) + &
         &       (a3ma * c2_ved_interp(klo,ilay) + b3mb * c2_ved_interp(khi,ilay)) 
            dd =    ai *  d_ved_interp(klo,ilay)  +   bi * d_ved_interp(khi,ilay) + &
         &       (a3ma * d2_ved_interp(klo,ilay) + b3mb * d2_ved_interp(khi,ilay))     
                        
                        
            if (linversion) then
                        
        
                dcdsig(1:nlay1d,iRxlayer) =    ai *  dcdsig_ved_interp(klo,1:nlay1d,ilay)  +  & 
                 &                             bi *  dcdsig_ved_interp(khi,1:nlay1d,ilay)  +  &
                 &            (a3ma * dcdsig2_ved_interp(klo,1:nlay1d,ilay)  + b3mb * dcdsig2_ved_interp(khi,1:nlay1d,ilay))
                 
                dddsig(1:nlay1d,iRxlayer) =    ai *  dddsig_ved_interp(klo,1:nlay1d,ilay)  +  &
                 &                             bi *  dddsig_ved_interp(khi,1:nlay1d,ilay)  +  &
                 &            (a3ma * dddsig2_ved_interp(klo,1:nlay1d,ilay) +  b3mb * dddsig2_ved_interp(khi,1:nlay1d,ilay))
    
                    
            endif
            
            exit
            
        endif
    enddo
      

    end subroutine splint_ved
    
    
    
    end module dipole1d
                                             