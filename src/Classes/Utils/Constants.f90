!**
!
! math_constants.f90 defines math constants and other parameters
!
! AUTHORS  Gary Egbert, Anna Kelbert & Naser Meqbel
!          College of Earth, Ocean and Atmospheric Sciences.
!*
module Constants
  implicit none
  
  character (len = 4), parameter  :: NODE   = 'NODE'
  character (len = 4), parameter  :: FACE   = 'FACE'
  character (len = 4), parameter  :: EDGE   = 'EDGE'
  character (len = 6), parameter  :: CENTER = 'CELL'
  character (len = 6), parameter  :: CORNER = 'NODE'
  character (len = 4), parameter  :: CELL = 'CELL'  
  character (len = 10), parameter :: CELL_EARTH = 'CELL EARTH'

  !**
  ! Possible node types:
  !*  
  character(len = 5), parameter :: XFACE = 'XFACE'
  character(len = 5), parameter :: XEDGE = 'XEDGE'
  character(len = 5), parameter :: YFACE = 'YFACE'
  character(len = 5), parameter :: YEDGE = 'YEDGE'
  character(len = 5), parameter :: ZFACE = 'ZFACE'
  character(len = 5), parameter :: ZEDGE = 'ZEDGE'
  
  ! Use this to select single or double precision
  integer, parameter :: SP = selected_real_kind (6, 37)
  integer, parameter :: DP = selected_real_kind (15, 307)
  integer, parameter :: prec = DP

  real (kind = prec), parameter :: PI   = 3.14159265357898_prec
  complex (kind = prec), parameter :: MU_0 = PI*.0000004_prec
  
  ! Important: sign convention used throughout the program
  integer, parameter :: ISIGN = -1
  
  ! Conductivity of the air for computational purposes
  real (kind = prec), parameter :: SIGMA_AIR = 1.0e-10
  
  ! Variable used to decide if a cell is air or ground
  ! in the presence of topography:
  ! Starting from the top of the air layers, program looks
  ! for the first cell in each column with conductivity
  ! exceeding SIGMA_MIN ... top of this cell is taken to be
  ! Earth's surface for this column)
  real (kind = prec), parameter :: SIGMA_MIN = 1.0e-6
  
  ! Useful geophysical constants (thinsheet: use 6358.35 for
  ! Kuvshinov; 6321.0 for Constable)
  real (kind = prec), parameter :: EARTH_R = 6371.0_prec !6378.164
  real (kind = prec), parameter :: CRUST_R = 6358.35_prec
  
  ! Useful conversion constants
  real (kind = prec), parameter :: D2R = PI/180._prec
  real (kind = prec), parameter :: R2D = 180._prec/PI
  real (kind = prec), parameter :: KM2M = 1000.0_prec
  real (kind = prec), parameter :: M2KM = 0.001_prec
  
  ! Smallest possible distance in radians or meters, that
  ! is used to prevent problems with machine errors in if statements.
  real (kind = prec), parameter :: EPS_GRID = 1.0e-4
  
  ! Real and complex precision constants / tolerance
  real (kind = prec), parameter :: LARGE_REAL = 1.0e13
  real (kind = prec), parameter :: TOL4 = 0.0001_dp
  real (kind = prec), parameter :: TOL6 = 0.000001_dp
  real (kind = prec), parameter :: TOL8 = 0.00000001_dp
  real (kind = prec), parameter :: R_TINY = 1.0e-13
  
  real (kind = prec), parameter :: EIGHT = 8.0_prec
  real (kind = prec), parameter :: THREE = 3.0_prec
  real (kind = prec), parameter :: TWO = 2.0_prec
  real (kind = prec), parameter :: ONE = 1.0_prec
  real (kind = prec), parameter :: R_ZERO = 0.0_prec
  real (kind = prec), parameter :: R_ONE = 1._prec    
  real (kind = prec), parameter :: MinusONE = -1.0_prec
  real (kind = prec), parameter :: MinusTWO = -2.0_prec
  
  complex (kind = prec), parameter :: C_ONE = (1.0_prec,0.0_prec)
  complex (kind = prec), parameter :: ONE_I = (0.0_prec,1.0_prec)
  complex (kind = prec), parameter :: C_ZERO = (0.0_prec, 0.0_prec)
  complex (kind = prec), parameter :: C_MinusOne = (-1.0_prec, 0.0_prec)
  
  ! Useful character constants
  character (len = 40), parameter :: DATA_FILE_TITLE = 'Synthetic 3D MT data written by ModEM-OO'
  character (len = 2), parameter :: TE = 'TE'
  character (len = 2), parameter :: TM = 'TM'
  
  character (len = 3), parameter :: FWD = 'FWD'
  character (len = 3), parameter :: TRN = 'TRN'
  character (len = 3), parameter :: ADJ = 'ADJ'
  
  ! Coordinate type options
  character (len = 80), parameter :: CARTESIAN = 'Cartesian'
  character (len = 80), parameter :: SPHERICAL = 'Spherical'
  
  ! Grid type options (global will always use spherical coords)
  character (len = 80), parameter :: REGION = 'Regional'
  character (len = 80), parameter :: SPHERE = 'Global'
  
  ! Job options
  character(*), parameter :: JOB_KEY        = 'JOB'
  character(*), parameter :: JOB_FORWARD    = 'FORWARD'
  character(*), parameter :: JOB_READ_WRITE = 'READ_WRITE'
  
  ! Grid Options
  character(*), parameter :: GRID_KEY                 = 'GRID TYPE'
  character(*), parameter :: GRID_STANDARD_CSG        = 'CARTESIAN'   
  character(*), parameter :: GRID_MULTIRESOLUTION_CSG = 'CARTESIAN_MR' 
  ! Multiresolution grid parameters
  integer     , parameter :: MAX_NGRIDS = 100
  character(*), parameter :: GRID_NGRIDS_KEY     = 'NGRIDS'
  character(*), parameter :: GRID_COARSENESS_KEY = 'CS'
  
  ! Resistivity model file
  character(*), parameter :: RMODEL_KEY = 'RESISTIVITY MODEL'
  
  ! Data file
  character(*), parameter :: DATA_FILE_KEY = 'DATA FILE'
  character(*), parameter :: WDATA_FILE_KEY = 'OUTPUT DATA FILE'
  
  ! Forward models
  character(*), parameter :: FWD_KEY  = 'FORWARDEL'
  character(*), parameter :: FWD_CCGD = 'FWD_CCGD'
  
  ! Finite difference (FD) topology operators options
  character(*), parameter :: FD_TOPOLOGY_KEY = 'FD_TOPOLOGY'
  character(*), parameter :: FD_TOPOLOGY_CSG = 'FD_TOPOLOGY_CSG'
  character(*), parameter :: FD_TOPOLOGY_CMR = 'FD_TOPOLOGY_CMR'
  
  ! Finite differences (FD) metric operators options
  character(*), parameter :: FD_METRIC_KEY = 'FD_METRIC'
  character(*), parameter :: FD_METRIC_CSG = 'FD_METRIC_CSG'
  character(*), parameter :: FD_METRIC_CMR = 'FD_METRIC_CMR'
  
end module Constants
