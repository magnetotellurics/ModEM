File Formats
=============

There are four principal file types: Model (providing grid+conductivity; input
or output); Data (input or output); Covariance (required to specify topography
or bathymetry, or to “freeze” conductivity; e.g., in the ocean, or to modify
default smoothing parameters; input); EM solution (electric fields; input or
output). There are also two simple optional control files. Some additional files
associated with the nonstandard (development/testing) options (``-J``, ``-M``,
``-T``, etc.) are not discussed here. The Model, Data, Covariance, and control
files are all flat ASCII files. The EM solution file is very large, and thus
uses a binary format. In typical usage (e.g., nested modelling; see below), the
EM solution files from one run are used as inputs to a second run. Care should
be taken if the files are moved between computers, as binary formats may
possibly be incompatible.  There are matlab functions that can read the binary
EM solution files, at least when all work is done on a fixed computer (see
section on matlab I/O) We summarize formats of the basic ASCII input and output
files next.

.. _data-files:

Data Files
-----------

ModEM data files use the so-called “list format”. This ASCII format (with
essentially each observation a separate item with full metadata, in a list)
results in larger files than necessary, but is very flexible. As data files are
typically not so large (at least compared to other files encountered in 3D MT
inversion) the size is not usually much of an issue.  Note that data files are
required both as inputs and outputs, for both ``-F`` and ``-I`` options.  In the
case of forward modeling on input the data file provides essentially meta data:
the list of periods to be run, as well as a list of the types of responses to
compute (impedance, vertical field TF, etc.), and locations where these
predicted responses will be evaluated.  The list format is specific to ModEM, so
it requires a detailed description.

ModEM “list format” data files contain one or more blocks of data. Each block
corresponds to a single “type” of data (e.g., the full 2 × 2 impedance
tensor), generally for all sites and all periods. Each block consists of an
8-line header followed by the actual data.  As an example, the start of a single
data block is as follows:

.. code-block:: text

    # Cascadia data from XML files, error floor 5%, no topography
    # Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
    > Full_Impedance
    > exp(+i\omega t)
    > [mV/km]/[nT]
    > 0.00
    7
    > 45.276 -119.634
    > 10 109
    1.000000E+00 006-006 0.000 0.000 41000.000 39375.000 0.000 ZXX -1.700416E+03 3.455034E+02 6.167116E+02
    1.000000E+00 006-006 0.000 0.000 41000.000 39375.000 0.000 ZXY 1.056768E+04 -1.259106E+04 6.167116E+02
    1.000000E+00 006-006 0.000 0.000 41000.000 39375.000 0.000 ZYX -1.761286E+04 1.617083E+04 6.167116E+02
    1.000000E+00 006-006 0.000 0.000 41000.000 39375.000 0.000 ZYY 2.497673E+03 -1.202580E+03 6.167116E+02
    1.000000E+01 006-006 0.000 0.000 41000.000 39375.000 0.000 ZXX -9.602744E+02 3.031632E+02 1.981138E+02
    1.000000E+01 006-006 0.000 0.000 41000.000 39375.000 0.000 ZXY 3.743773E+03 -3.790398E+03 1.981138E+02
    1.000000E+01 006-006 0.000 0.000 41000.000 39375.000 0.000 ZYX -6.208288E+03 5.470545E+03 1.981138E+02
    1.000000E+01 006-006 0.000 0.000 41000.000 39375.000 0.000 ZYY 2.279075E+03 -1.333596E+03 1.981138E+02

The actual data are written in a list following the 8-line header. Each data
line (which are split over two lines in the example displayed here, but should
be on a single line in the actual data file) in the list contains one data
component appropriate to the data type, following the format described in the
second line of the header. For example, for the ``Full_Impedance`` data type in
the example above, ``ZXX``, ``ZXY``, ``ZYX`` and ``ZYY`` would be data
components; for each period/site these data components would be listed on a
separate line. Thus a complete set of data components for a data type for even a
single site will typically require multiple lines (4 in the example above);
however, real and imaginary parts of an intrinsically complex data type (such as
impedance components) are given together on a single line, with a common error
bar.

Some (or all) components can be missing at any sites and/or for any periods, and
any ordering of the list elements is allowed. This is an advantage of the
format: it is easy to edit a file and delete or add observations, or to have a
program that does this.  Note however, that if components are missing for all
sites (e.g., you are inverting only off-diagonal impedances for all, or a large
group of sites), it is usually better to use a more restrictive data type to
minimize internal storage and computations (e.g., use ``Off_Diagonal_Impedance``
in place of ``Full_Impedance``). Multiple data blocks are allowed (typically
corresponding to different types), so if a large group of sites had only
off-diagonal impedances, and another group had full impedances, it would be most
efficient to use two blocks, with distinct data types.

The first two lines of the header are descriptive only (i.e., are primarily for
a human reading the file, as indicated by the # at the start of the line. But
note that this symbol is not parsed by the program.)

Line 1: # comment
^^^^^^^^^^^^^^^^^

Comment line, intended to describe your data; will be read and replicated in
output responses. Max 100 chars.

Line 2: # comment
^^^^^^^^^^^^^^^^^

Comment line, intended to describe the contents of each data item in the block; will be
ignored by the code.

The remaining lines give details about the data type. Not all of these are used at present.

Line 3:  data_type
^^^^^^^^^^^^^^^^^^^

Keyword for the data type stored in this data block. Only specific keywords are
supported. For each data type there are a set of additional keywords denoting
specific components for the type. These may be real or complex.

For the 3D MT program these are:

+--------------------------+----------------------------+-----------------+
| Keyword                  | Allowed Components         | Real or Complex |
+==========================+============================+=================+
| Full_Impedance           | ZXX, ZXY, ZYX, ZYY         | Complex         |
+--------------------------+----------------------------+-----------------+
| Off_Diagonal_Impedance   | ZXY, ZYX                   | Complex         |
+--------------------------+----------------------------+-----------------+
| Full_Vertical_Components | TX, TY                     | Complex         |
+--------------------------+----------------------------+-----------------+
| Full_Interstation_TF     | MXX, MXY, MYX, MYY         | Complex         |
+--------------------------+----------------------------+-----------------+
| Off_Diagonal_Rho_Phase   | RHOXY, PHSXY, RHOYX, PHSYX | Real            |
+--------------------------+----------------------------+-----------------+
| Phase_Tensor             | PTXX, PTXY, PTYX, PTYY     | Real            |
+--------------------------+----------------------------+-----------------+

For 2D MT these are:

+--------------------------+----------------------------+-----------------+
| Keyword                  | Allowed Components         | Real or Complex |
+==========================+============================+=================+
| TE_Impedance             | TE                         | Complex         |
+--------------------------+----------------------------+-----------------+
| TM_Impedance             | TM                         | Complex         |
+--------------------------+----------------------------+-----------------+

Line 4: > sign_convention
^^^^^^^^^^^^^^^^^^^^^^^^^

This defines the sign convention assumed for time dependence, which is
effectively assumed in Fourier transform (or other) data processing, which can
be plus or minus. Here this is written as, ``exp(+i omega t)`` for clarity. In
fact, ModEM searches for a minus sign in this line; if not found, plus is
assumed. The data are then converted to the sign convention used internally in
the code (which is determined at compile time by the global variable ``ISIGN``).
Any output data files will follow the convention used in the input file,
although internal calculations will follow the convention defined in ``ISIGN``.

Line 5: > units
^^^^^^^^^^^^^^^

Units for this data type. Internally in the code, the SI units for E/B are used
for MT impedances. However, the user is given a choice for this data type.
Supported options are:

1. SI units for E/B: ``[V/m]/[T]`` (used internally in ModEM)
2. practical units for E/B: ``[mV/km]/[nT]``
3. SI units for E/H: ``[V/m]/[A/m] = Ohm``

The output data will be given in the same units as the input, with the program
converting to the internally used units after input/before output. The units
supplied in the data file must of course match the actual data. Dimensionless
data (such as the vertical magnetic field transfer functions, or inter-station
transfer functions) are identified by empty "[]" units. Apparent resistivity and
phase are loosely viewed as being dimensionless. It is assumed that apparent
resistivity is given in ohm-m; there is no user option on units for

Line 6: > orientation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Orientation of coordinate axis. **THIS IS NOT PRESENTLY IMPLEMENTED IN THE
CODE!!**  At present all impedances are assumed to be rotated into the same
Cartesian coordinate system used by the model grid. This header information may
be useful for other programs that plot data, so it would be good practice to set
the orientation to a meaningful value, even though the inversion code does not
reference this number.

.. warning::

   Coordinate orientation is currently not implemented in the code! This line
   will be ignored.

Line 7: > geographic_origin
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Latitude and longitude for the ``(0, 0, 0)`` point in the data file. Not used in
ModEM, but could be useful for plotting, and to match the origin of the model to
that of the data.  See discussion of coordinate origins in description of model
file formats. Can set to zero if unknown since it does not effect program
function.

Line 8: > nperiods nsites
^^^^^^^^^^^^^^^^^^^^^^^^^

These two integers give the number of periods and sites in the data set. More
properly, these give the *maximum number* of periods and sites. The program will
parse the list of data components, and figure out the actual numbers. Thus it
will never hurt to set these numbers to a larger value than you actually need.
However, if these numbers are smaller than required to read your data block, a
segmentation fault will occur. If unsure, use the command

.. code-block:: bash
    
    ./Mod3DMT -R model_input data_input model_output data_output

to validate your model and data files. This will read your files, store them in
the internal dictionary, data and model structure, write them out again, and
exit. If this does not work, you need to increase the size of the parameters.
(This should be fixed!)

In the example given above, each line in the actual data list contains period,
site code (up to 12 chars), geographic latitude and longitude, ``X/Y/Z``
location in meters, component code (e.g., ``ZXY``), real and imaginary parts,
and the error bar.  Geographic latitude and longitude are not used by the
program, which requires the user to make sure that the correct locations (i.e.,
``X/Y/Z`` in meters are provided, registered to the grid; see further
discussion in :ref:`model-files`.

Of course each specific data type requires a specific format, and consistent
header entries, such as units (described above). List element formats for
currently supported data types are summarized here.Some data file examples are
provided in the `examples/` directory. The ``matlab/ioAscii/`` directory contains
read and write routines in Matlab (at present these reading/writing routines are
restricted to the most common data types:

* Off_Diagonal_Impedance
* Full_Impedance
* Full_Impedance_plus_Full_Vertical_Components

For complex impedances and vertical components the data format is::

    Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error

For apparent resistivity and phase the data format looks like::

    Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Value Error

For inter-station transfer functions, information about the reference site is also required::

    Period(s) Code GG_LatLon XYZ(m) Ref_Code Ref_LatLon Ref_XYZ(m) Component Real Imag Error

For 2D MT, the data format is identical, but the X coordinate is not used or
stored in the program.

.. _model-files:

Model Files
------------

ModEM support two model file formats, based on previously used formats for 3D
resistivity (not conductivity!) models: The first (and the current default) is
based on that used in Weerachai Siripunvaraporn’s WSINV3DMT program. There are a
few differences in the ModEM implementation: (1) ModEM does not support the
option to use a set of resistivity codes (i.e., integers 1-9, used to specify a
small number of fixed resistivity values); (2) we have added the option to store
the resistivity values in natural log or linear formats; (3) an extra line at
the end of the file is used to specify a coordinate origin, as discussed below.
The second format is identical to that used for 3D models in WingLink (referred
to in the ModEM code, perhaps inappropriately, as the “Mackie format”.) These
two formats are similar, but not identical. At present switching between these
model formats requires editing and recompiling the program.

(For completeness here is a very brief summary of the modifications required:
read routines are in ``f90/3D_MT/modelParam/modelParamIO/*.inc``.  All of these
files are included in the compiled ModelSpace module (with compiler directives
``#include``). The names of the read and write routines are
``read_modelParam_WS(mackie)`` and similarly for write. The read/write routines
in use are overloaded on ``read_modelParam`` (and similarly for the write
routine; this is done in ``f90/3D_MT/modelParam/ModelSpace.f90``), and calls to
these generic read/write subroutines are then used in the program to input or
output model parameter objects. Both the "WS" and "Mackie" read/write routines
(and others that one might want to implement) must have identical interfaces.
Then, by changing the names of routine overloaded in the ``ModelSpace``
interface, the expected file format used for model parameter I/O is changed.
Yes, better approaches could be imagined!)

Both of the currently supported model file formats begin with a description of
the tensor-product Cartesian grid: cell dimensions in ``x``, ``y``, and ``z``
are given (in meters, i.e., ``∆xi,i = 1,...,Nx; ∆yi,i = 1,...,Ny ; ∆zi,i =
1,...,Nz;``), followed by a list of linear or log resistivity values. The two
supported formats differ primarily in the ordering of the list of resistivity
values.

A critical issue which is independent of file format, is how coordinates
implicit in the model grid are related to coordinates specified in the separate
data file. The arrays ∆x, ∆y, ∆z define a natural Cartesian coordinate system for
the grid, with origin at the outer corner (of model parameter cell (1,1,1)). We
refer to this here as the “natural grid coordinate system”. The origin of the
coordinate system used in the data file (the “data coordinate system” does not
have to be (and in fact most often is not) the same. The relationship between
the two origins must obviously be defined. This is done in ModEM by giving, in
the model file, the coordinates (in natural grid coordinates) of the origin of
the data coordinate system. A common approach might be to center the grid on a
point chosen near the center of the array of observations. Specifically, this
could correspond to the location given by the arithmetic mean of site
coordinates, or perhaps the coordinates of a centrally located site. In the data
file Cartesian coordinates of all data sites will be specified (in meters)
relative to this point. Then, if the coordinates of the reference point given in
the model parameter file corresponds to the center of the grid, the array of
sites will be centered on the numerical grid. If no origin is specified in the
model file (i.e., the final lines are missing, as would be the case if an actual
WSINV3DMT model input file were used), the horizontal (x,y) origin is taken to
be center of the grid, and the vertical origin the top of the first model layer
(nominally the surface of the Earth, but see :ref:`include-topo`). It is the user’s
responsibility to ensure that a sensible data origin is used, and that data
locations are given consistently relative to this origin. When there are
geographically fixed features (e.g., a coastline, topography) that need to be
represented in the numerical model, it is of course also essential to ensure
that consistent projections and coordinate systems are used for model and data
file setup. Tools for input setup (e.g., Grid3D; N. Meqbel) will commonly take
care of these details. One issue to be aware of is that if you change the grid
discretization, and the origin somehow changes (e.g., relative to a fixed
feature like the ocean) coordinates in the data file might need to be
checked/modified. Another issue to be aware of is the vertical coordinate of the
origin.  For a flat Earth (no topography) it is natural to take the surface to
be z = 0. However, if topography is included then the choice is less clear: 0
might be sea level, or the highest land elevation point in the domain. We
discuss this further in :ref:`include-topo`.

Following are further details on file formats:

.. code-block:: text

    # Written by JavaScript within EM-WebApp <--(1) header line
    59 54 25 0 LOGE <-- (2) grid dimensions, and more
    165357 152127 99214.01 ... <-- (3) grid dimensions: x
    157079 144512 119380 ... <-- grid dimensions: y
    100 120 140 160 ... <-- grid dimensions: z
    <-- (4) resistivity codes (blank line often)
    -1.20397276458951 -1.20397276458951 ... <-- (5) cell resistivities
    ...
    ...
    ...
    -1183956 -1088486 0 <-- (6) coordinate origin (last lines optional)
    0 <-- (7) coordinate rotation (not used)

Explanation/Notes:

(1) Header 
^^^^^^^^^^^

Header is not used by program; useful for user to describe model in some
fashion.  Here, the program that wrote the file is documented.

(2) Grid dimensions and more
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the three integers defining grid size (``Nx``, ``Ny``, ``Nz``),
there is a fourth integer in this line, for consistency with the original
WSINV3DMT model file format.  In the original, this was used to specify the
number of integer resistivity codes, with nCodes = 0 used to indicate that
resistivities are specified as actual values. This is the only option supported
by ModEM, and this integer should be left set to 0.  (Otherwise the program will
stop). The next entry on the line is a character string specifying ``LINEAR``,
``LOGE`` or ``LOG10`` resistivity. In fact the inversion always works in the
natural log domain, even if input and output use the linear option. Note that
output files will match the convention used in the input file.

(3) Grid Dimensions: x, y, z
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The next lines give the grid spacing, in meters. Format is free, and the list for any
of the three arrays can span multiple lines. Use a new line for each of x, y, and z.

(4) Resistivity Codes - (Not used by ModEM)
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The next line is where resistivity codes were to be specified in the WSINV3DMT
files. This line is left blank here, as ModEM does not support the resistivity code
option.

.. note::

    Resistivity codes are not used in ModEM.


(5) Cell Resistivities
^^^^^^^^^^^^^^^^^^^^^^

Now the list of cell resistivities (or integer codes, if this option is being
used) are given. This is free format, and the full list of cell resistivities
will typically span many lines. To fill the array ``rho(Nx:-1:1,1:Ny,1:Nz)`` use
these read statements (note the reversal of order of the index ``ix``):

.. code-block:: text

    DO iz = 1,Nz
    DO iy = 1,Ny
        DO ix = Nx,1,-1
            READ(10,*) rho(ix,iy,iz)
        ENDDO
    ENDDO
    ENDDO

.. note::

   Note the reversal of order of the index ``ix`` in the above pseudocode.

(6) Coordinate Origin (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These three numbers provide the origin (x, y and z), in 'natural grid
coordinates' used for Cartesian data locations.  This line is an addition to the
WSINV3DMT file format, and is optional; default is 0.  In any event, one has to
make sure that the XYZ coordinates of the sites are referred to the same origin
(see :ref:`data-files`).

(7) Coordinate Orientation (Not used/Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The last (also optional) line in the file is for orientation of the model grid.
This can be useful for plotting programs, but is not used by ModEM. It is
important to note that while it is possible to rotate the grid to any
orientation, all transfer function types must be rotated into the grid
coordinate system: x in the impedance is the same as x in the model grid. Again,
it is the users responsibility to ensure consistency between grid and data
rotations.

For 2D MT, the WingLink (Mackie) model format is used, with the same
modification for natural log resistivity. The WingLink 3D format can also be
used for model input/output to ModEM, but this requires slight code
modification, and recompiling, as discussed above.

.. _model-cov-files:

Model Covariance Files
-----------------------

Model regularization is handled through a positive-definite smoothing operator,
that can be viewed as a prior covariance on the model perturbation–i.e., the
deviation of the model parameter (log resistivity) from an assumed model. A
prior (or prejudice) model is required at present in the ModEM regularization
approach. In practice most users will usually want to use something simple, such
as a half space of constant resistivity, but any resistivity model can be used
for the prior. In the simplest usage, no covariance file is required; there are
default values for all parameters. However, for many purposes (including
topography, or a fixed conductivity feature such as an ocean) the covariance
file must be used. The covariance file also allows the user to control smoothing
length scales separately in all directions, and as a function of depth in each
of the vertical layers of the model. An additional parameter defines vertical
smoothing independently. One can also switch off (or reduce) the smoothing
across a boundary. Air and ocean cells are also specified through the model
covariance file.

Example covariance files are provided in ``"examples/3D_MT/``, with file
extension ``.cov``. The file has a 16-line header (Fixed length!) that is
ignored by the code. This is followed by grid dimensions, and parameters that
define the auto-regressive smoothing scheme. These ⍺ between 0 and 1; 0 implies
no smoothing, larger values correspond to longer smoothing (correlation) length
scales. Typical values that we use are 0.1 - 0.3.  The smoothing can be
specified independently for each vertical level, and for x and y directions. A
single parameter defines the vertical smoothing. Note that the smoothing is
applied without reference to the actual variable grid spacing. In addition to
the parameters ⍺, one can specify how many times the smoothing is repeated
(nSmooth). Typically this number is 1-3. Repeating the smoothing implies longer
decoration length scales—these increase length scales by a factor of
sqrt(nSmooth).  Using a smaller ⍺ and more smoothing repetitions will result in
a more "circular" effective smoothing kernel. With nSmooth = 1 the covariance
will be somewhat anisotropic, with longer length scales associated with grid
coordinate axes (i.e., features may tend to be elongated in x and y directions,
relative to a 45 degree angle). Default values for the smoothing parameters
(used if no model covariance file is provided) are :math:`⍺=0.3` and
:math:`nSmooth = 1`. More circular smoothing, with similar length scales would
be achieved with :math:`⍺=0.2` and :math:`nSmooth = 2`.

The covariance file also is used to subdivide the model domain, by defining an
integer index for each cell in the grid. The index can be used to enable
freezing of some regions (so that cells in the region are not changed during the
inversion), or to turn off (or reduce) smoothing across the boundaries between
adjacent regions. To turn off smoothing, one provides "exception rules" which
define how the covariance should smooth across a boundary between regions tagged
with two different indices. All of this information (array of indices for all
grid cells; exception rules) is provided after defining the smoothing
parameters. The most common situation is including the ocean (with conductivity
fixed; index = 9) and inclusion of topography (where cells that are in the air
should be frozen to the very high air resistivity; index 0). It is possible to
freeze other parts of the model domain, for example in hypothesis testing
applications. Thus, to include an ocean or topography in an inversion, it is
necessary to provide a covariance file.

Following is a simple example of a covariance file for a small model: 

.. code-block:: text

    +-----------------------------------------------------------------------------+
    | This file defines model covariance for a recursive autoregression scheme.   |
    | The model space may be divided into distinct areas using integer masks.     |
    | Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between |
    | air, ocean and the rest of the model is turned off automatically. You can   |
    | also define exceptions to override smoothing between any two model areas.   |
    | To turn off smoothing set it to zero. This header is 16 lines long.         |
    | 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |
    | 2. Smoothing in the X direction (NzEarth real values)                       |
    | 3. Smoothing in the Y direction (NzEarth real values)                       |
    | 4. Vertical smoothing (1 real value)                                        |
    | 5. Number of times the smoothing should be applied (1 integer >= 0)         |
    | 6. Number of exceptions (1 integer >= 0)                                    |
    | 7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 3 & 4) |
    | 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|
    +-----------------------------------------------------------------------------+

    21 28 11    <-- Grid dimensions excluding air layers (Nx, Ny, NzEarth)

    0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 <-- Smoothing in the X direction
    0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 <-- Smoothing in the Y direction
    0.1   <-- Vertical smoothing

    2     <-- Number of times smoothing is applied

    1     <-- Number of exception rules
    2 4 0.     <-- Exception rule #1: smoothing across region 2, 4 boundary 
                    this should be a real number; 0. means no smoothing
                    repeat these lines for each exception
    1 11   <-- Start index array specification; see more details below.
            
    9 9 9 9 9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 9 9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 4 4 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 4 4 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 4 4 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 4 4 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 4 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4

The index array is specified through a series of blocks, which specify indices
for a range of vertical layers. First the vertical layer limits are given ``(k1,
k2)`` then an array ``temp(i,j)``. Then ``index(i,j,k1:k2)`` = temp. Here the
array temp (for a single range of layers) is read as:

.. code-block:: text

    DO i = 1,Nx
    DO j = 1,Ny
        read(fid,*) temp(i,j)
    ENDDO
    ENDDO

Blocks are repeated as necessary — a single block is used in the simple example,
indicating that the regions extend over all 11 layers (top-to-bottom) of this
very small grid.

H2 (Bi-Helmholtz) Covariance (CovType=2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ModEM now supports an alternative regularization method based on the 
H2 (Bi-Helmholtz) operator, which provides isotropic smoothing with a 
(Matérn) covariance structure.

To use H2 smoothing, set ``CovType=2`` in the inversion control file
or command line, and provide a covariance file in the new unified format.

The unified covariance file format uses line 17 to specify both grid dimensions
and the CovType:

- **3 values** (``Nx Ny NzEarth``): legacy format, defaults to CovType=1 (AR)
- **4 values** (``Nx Ny NzEarth 2``): new format, CovType=2 (H2)

The file structure for CovType=2 is:

.. code-block:: text

    Line 17:  Nx Ny NzEarth 2          <-- grid dimensions + CovType
    Line 18:  Sx(1:NzEarth)             <-- placeholder (not used by H2)
    Line 19:  Sy(1:NzEarth)             <-- placeholder (not used by H2)
    Line 20:  Sz                         <-- placeholder (not used by H2)
    Line 21:  N kappa pcg_tol pcg_maxIt  <-- parameter array
    Line 22:  nrules                     <-- number of exception rules
    Line 23+: mask1 mask2 smoothing      <-- exception rules (3 values per rule)
    Line 24+: mask blocks                <-- same format as CovType=1

The H2 parameters are:

- **kappa**: dimensionless smoothing parameter; correlation length ~ 1/kappa cells
- **pcg_tol**: PCG solver tolerance (typically 1e-7)
- **pcg_maxIt**: maximum PCG iterations (typically 250)
- **N**: number of smoothing passes (typically 1, not used by H2)

Example CovType=2 covariance file:

.. code-block:: text

    +-----------------------------------------------------------------------------+
    | This file defines model covariance for the H2 (Bi-Helmholtz) smoother.      |
    | ...                                                                         |
    +-----------------------------------------------------------------------------+
    21 28 11 2   <-- Grid dimensions (Nx, Ny, NzEarth) + covtype
    0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 <-- Sx, not used by H2
    0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 <-- Sy, not used by H2
    0.1   <-- Sz, not used by H2

    1 0.5 1e-7 250 <-- Number of smoothings, kappa, PCG tolerance, max PCG iter

    1     <-- Number of exception rules
    2 4 0.     <-- Exception rule #1: smoothing across region 2, 4 boundary 
                    this should be a real number; 0. means no smoothing
                    repeat these lines for each exception
    1 11   <-- Start index array specification
            
    9 9 9 9 9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 9 9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 9 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    9 9 9 4 4 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 4 4 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 4 4 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 4 4 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 4 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 1 1 1 1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4

The Sx, Sy, Sz values are placeholders for future extensions and are not used in 
the current H2 implementation. Exception rules use the same 3-value 
format as CovType=1; the third value (smoothing) is read but not used for
CovType=2 (tears are determined by rule presence).

**Cross-compatibility**: The code is tolerant when file and CovType don't match.
If a CovType=2 file is read with CovType=1 selected (or vice versa), a warning
is printed and the mask is read from the file while defaults are used for the
other parameters.


For 2D MT, Weerachai Siripuvaraporn's model covariance is used. No configuration
file is currently allowed.

Including an Ocean
^^^^^^^^^^^^^^^^^^^

To include an ocean domain, which will have conductivity that is not changed
("frozen") during the inversion run, the covariance file must be used. As noted
in the header to the example covariance file index 9 is used to denote the ocean
domain that will be frozen to the value specified in the prior model. Note that
you also need to specify the conductivity of seawater that you want to use in
the prior model; this would typically be roughly 0.3 ohm-m (but in some cases
might be quite different e.g., when modeling a highly saline body of water
such as the Dead Sea). The same index 9 would be used to freeze a structure in
hypothesis testing. Future versions of ModEM will allow multiple indices to
denote frozen structure, but at present you have to use index 9.

Here is an example of the covariance mask blocks used to specify an ocean: 

.. code-block:: text

    1 1   <-- first layer
            
    9 9 9 9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

    2 2   <-- second layer
            
    9 9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

    3 11   <-- remaining layers
            
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

In this example there are just 2 ocean layers, and the remainder of the model is
indexed by the default mask 1. Assuming a 100 ohm-m prior (for example) the same
pattern would be given in the prior model, with cells that are marked as 9 set
to 0.3 (or ln(0.3)) and those marked with 1 set to 100 (or ln(100)). (But note
that the order of elements in the arrays differs between covariance and model
files!)


.. _include-topo:

Including Topography
^^^^^^^^^^^^^^^^^^^^^

To include topography you need to add extra layers to the top of the model, with
some (typically most) cells in these layers marked as air. As with ocean cells,
the air cells need to be indicated in both the prior model (set to air
resistivity–a large positive number and in the covariance file, marked with
index 0. Here is a simple example of the covariance file for a simple circular
hill (just showing the mask blocks). 

.. note::

    Note: We use resistivity in the model file 

.. code-block:: text

    1 1   <-- top layer -- mostly air +_  top of hill
            
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

    2 2   <-- second layer -- mostly air +_ more hill
            
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

    3 13   <-- remaining layers
            
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1


Here we have essentially added the hill by including two additional layers (in
the 11 Earth layers in the original example model). In this example all of the
Earth cells are marked by the default index 1, and only the air cells in the
additional layers have a different mask index, namely 0. Of course you can use
multiple indices with topography (e.g., including also an ocean).

If you do include topography you need to be careful about the definition of the
vertical coordinates of data sites (specified in the data file), and make sure
that you define the data origin (given in "natural grid coordinates" in the
model file) consistently. The key points are these: (1) z is positive downwards;
(2) the model parameter now explicitly includes the extra layers used to define
topography. Thus, the origin for the "natural grid coordinates" (see
:ref:`model-files`) is the outer corner at the top of the uppermost topography
layer–this corresponds to z=0 in these coordinates. (3) data sites should be
projected onto the surface of the (model) Earth.

One approach would be to define the data origin at sea level (or perhaps at the
first layer that has no air cells). Then in data coordinates all of the layers
above this surface would have negative z coordinate, and any data sites sitting
on this topography should be projected to this level. If you choose the data
origin as the highest point in the model (i.e., coinciding with natural grid
z-coordinate = 0) then all data sites will have non-negative z-coordinates. A
very common error in using topography is to project sites to the wrong
z-coordinate, so that sites are effectively inside the Earth, or in the air. The
inversion obviously will not work correctly if you do this, so if you are using
topography and having problems, double check that you are defining vertical
levels consistently in the model and data files.  Again, using tools for model
and data file setup (such as Grid3D) can help ensure consistency–if these are
used correctly!
