! *****************************************************************************
! math_constants.f90 defines math constants and other parameters

module math_constants
  implicit none

  ! Use this to select single or double precision
  integer, parameter             :: SP = selected_real_kind(6,37)
  integer, parameter             :: DP = selected_real_kind(15,307)
  integer, parameter			 :: prec = DP

  real (kind=prec), parameter	 :: PI = 3.14159265357898_prec
  real (kind=prec), parameter	 :: MU_0 = PI*.0000004_prec

  ! Important: sign convention used throughout the program
  integer, parameter             :: ISIGN = -1

  ! Conductivity of the air for computational purposes
  real (kind=prec), parameter	 :: SIGMA_AIR = 1.0e-10

  ! Variable used to decide if a cell is air or ground in the presence of topography:
  ! starting from the top of the air layers, program looks for the first cell in each
  ! column with conductivity exceeding SIGMA_MIN ... top of this cell is taken to be
  ! Earth's surface for this column)
  real (kind=prec), parameter    :: SIGMA_MIN = 1.0e-6

  ! Useful geophysical constants (thinsheet: use 6358.35 for Kuvshinov; 6321.0 for Constable)
  real (kind=prec), parameter	 :: EARTH_R = 6371.0_prec !6378.164
  real (kind=prec), parameter	 :: CRUST_R = 6358.35_prec

  ! Useful conversion constants
  real (kind=prec), parameter	 :: D2R = PI/180._prec
  real (kind=prec), parameter	 :: R2D = 180._prec/PI
  real (kind=prec), parameter	 :: KM2M = 1000.0_prec
  real (kind=prec), parameter	 :: M2KM = 0.001_prec


  ! Smallest possible distance in radians or meters, that is used to prevent
  ! problems with machine errors in if statements
  real (kind=prec), parameter    :: EPS_GRID = 1.0e-4

  ! Real and complex precision constants / tolerance
  real (kind=prec), parameter	 :: LARGE_REAL = 1.0e13
  real (kind=prec), parameter    :: TOL4= 0.0001_dp
  real (kind=prec), parameter    :: TOL6= 0.000001_dp
  real (kind=prec), parameter    :: TOL8= 0.00000001_dp
  real (kind=prec), parameter    :: R_TINY= 1.0e-13

  real (kind=prec), parameter	 :: EIGHT = 8.0_prec
  real (kind=prec), parameter	 :: THREE = 3.0_prec
  real (kind=prec), parameter	 :: TWO = 2.0_prec
  real (kind=prec), parameter	 :: ONE = 1.0_prec
  real (kind=prec), parameter	 :: R_ZERO = 0.0_prec
  real (kind=prec), parameter	 :: MinusONE = -1.0_prec
  real (kind=prec), parameter	 :: MinusTWO = -2.0_prec

  complex (kind=prec), parameter :: C_ONE = (1.0_prec,0.0_prec)
  complex (kind=prec), parameter :: C_ZERO = (0.0_prec, 0.0_prec)
  complex (kind=prec), parameter :: C_MinusOne = (-1.0_prec, 0.0_prec)

  ! Useful character constants
  character (len=2), parameter	 :: TE = 'TE'
  character (len=2), parameter	 :: TM = 'TM'

  character (len=3), parameter	 :: FWD = 'FWD'
  character (len=3), parameter	 :: TRN = 'TRN'
  character (len=3), parameter	 :: ADJ = 'ADJ'

  character (len=80), parameter  :: Cartesian = 'Cartesian'
  character (len=80), parameter  :: Spherical = 'Spherical'

end module math_constants !math_constants
