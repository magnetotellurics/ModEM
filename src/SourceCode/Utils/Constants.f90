!
!> math_constants.f90 defines math constants and other parameters
!
!> AUTHORS    Gary Egbert, Anna Kelbert & Naser Meqbel
!>     College of Earth, Ocean and Atmospheric Sciences.
!
module Constants
    implicit none
    !
    !> Program Version
    character( len=5 ), parameter :: VERSION = "1.0.1"
    !
    !>
    integer :: aux_counter = 0
    !
    !> FIELD TYPES
    integer, parameter :: real_t = 1
    integer, parameter :: complex_t = 2
    integer, parameter :: integer_t = 3
    !
    !> INDEXES TYPES
    integer, parameter :: ind_boundary = 1
    integer, parameter :: ind_interior = 2
    !
    !> VECTOR
    character( len=4 ), parameter :: EDGE = "EDGE"
    character( len=4 ), parameter :: FACE = "FACE"
    !
    !> SCALAR
    character( len=4 ), parameter :: NODE = "NODE"
    character( len=4 ), parameter :: CELL = "CELL"
    !
    character( len=4 ), parameter :: CELL_EARTH = "EART"
    !
    !> Possible node types:
    character( len=5 ), parameter :: XFACE = "XFACE"
    character( len=5 ), parameter :: XEDGE = "XEDGE"
    character( len=5 ), parameter :: YFACE = "YFACE"
    character( len=5 ), parameter :: YEDGE = "YEDGE"
    character( len=5 ), parameter :: ZFACE = "ZFACE"
    character( len=5 ), parameter :: ZEDGE = "ZEDGE"
    !
    !> Use this to select single or double precision
    integer, parameter :: SP = selected_real_kind( 6, 37 )
    integer, parameter :: DP = selected_real_kind( 15, 307 )
    integer, parameter :: prec = DP
    !
    real( kind=prec ), parameter :: PI = 3.14159265357898_prec
    complex( kind=prec ), parameter :: mu_0 = PI * .0000004_prec
    !
    !> Important: sign convention used throughout the program
    integer, parameter :: isign = -1
    !
    !> Conductivity of the air for computational purposes
    real( kind=prec ), parameter :: SIGMA_AIR = 1.0d-10
    !
    !> Variable used to decide if a cell is air or ground
    !> in the presence of topography:
    !> Starting from the top of the air layers, program looks
    !> for the first cell in each column with conductivity
    !> exceeding SIGMA_MIN ... top of this cell is taken to be
    !> Earth"s surface for this column)
    real( kind=prec ), parameter :: SIGMA_MIN = 1.0e-6
    !
    !> Useful geophysical constants (thinsheet: use 6358.35 for
    !> Kuvshinov; 6321.0 for Constable)
    real( kind=prec ), parameter :: EARTH_R = 6371.0_prec !6378.164
    real( kind=prec ), parameter :: CRUST_R = 6358.35_prec
    !
    !> Useful conversion constants
    real( kind=prec ), parameter :: D2R = PI/180._prec
    real( kind=prec ), parameter :: R2D = 180._prec/PI
    real( kind=prec ), parameter :: KM2M = 1000.0_prec
    real( kind=prec ), parameter :: M2KM = 0.001_prec
    !
    !> Smallest possible distance in radians or meters, that
    !> is used to prevent problems with machine errors in if statements.
    real( kind=prec ), parameter :: EPS_GRID = 1.0e-4
    !
    !> Real and complex precision constants / tolerance
    real( kind=prec ), parameter :: R_LARGE = 1.0e13
    real( kind=prec ), parameter :: TOL4 = 0.0001_dp
    real( kind=prec ), parameter :: TOL6 = 0.000001_dp
    real( kind=prec ), parameter :: TOL8 = 0.00000001_dp
    real( kind=prec ), parameter :: R_TINY = 1.0e-13
    !
    real( kind=prec ), parameter :: EIGHT = 8.0_prec
    real( kind=prec ), parameter :: THREE = 3.0_prec
    real( kind=prec ), parameter :: TWO = 2.0_prec
    real( kind=prec ), parameter :: ONE = 1.0_prec
    real( kind=prec ), parameter :: R_ZERO = 0.0_prec
    real( kind=prec ), parameter :: R_ONE = 1._prec
    real( kind=prec ), parameter :: MinusONE = -1.0_prec
    real( kind=prec ), parameter :: MinusTWO = -2.0_prec
    !
    complex( kind=prec ), parameter :: C_ONE = (1.0_prec,0.0_prec)
    complex( kind=prec ), parameter :: ONE_I = (0.0_prec,1.0_prec)
    complex( kind=prec ), parameter :: C_ZERO = (0.0_prec, 0.0_prec)
    complex( kind=prec ), parameter :: C_MinusOne = (-1.0_prec, 0.0_prec)
    !
    !> Grid type options (global will always use spherical coords)
    character( len=8 ), parameter :: REGION = "Regional"
    character( len=6 ), parameter :: SPHERE = "Global"
    !
    character( len=7 ), parameter :: FORWARD = "FORWARD"
    character( len=7 ), parameter :: INVERSE = "INVERSE"
    character( len=5 ), parameter :: DERIV = "DERIV"
    !
    character( len=40 ), parameter :: DATA_FILE_TITLE_MT = "Synthetic 3D MT data written by ModEM-OO"
    character( len=42 ), parameter :: DATA_FILE_TITLE_CSEM = "Synthetic 3D CSEM data written by ModEM-OO"
    !
    character( len=4 ), parameter :: LOGE = "LOGE"
    character( len=6 ), parameter :: LINEAR = "LINEAR"
    character( len=6 ), parameter :: LOG_10 = "LOG_10"
    !
end module Constants
