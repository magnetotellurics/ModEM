.. _modem-namelist:

ModEM Namelist - Additional Settings
=====================================

ModEM now has a Fortran Namelist which can be used to supply ModEM with runtime
options. So far, this name list is completely optional; however, it may provide
some features which are convenient and useful.

A Fortran namelist is a file that can provide input for Fortran executables.
The ModEM namelist must be named ``modem.namelist.nl`` and must be present in
the same directory as you run the ModEM executable.

A copy of the ModEM namelist is listed below, but you can also use a ModEM
executable by using the ``-N`` command line argument:

.. code-block:: bash

    $ ./Mod3DMT_SP2 -N

This will create a file named ``modem.namelist.nl`` in the current directory.
If this file already exists, ModEM will *not* create this file.

Each namelist section is optional and can be omitted completely if desired.

.. code-block:: text
   :caption: Example ModEM Namelist file

    &settings
        ! 'output_level' below overrides -V passed in on the command line
        output_level = 'regular'
        ! Covariance Type: 'AR' (default),  'H2' (Bi-Helmholtz) 
        ! Currently performs no actions - use Inv control file or command line arg
        covariance_type = 'AR'
    /
    &io
         data_input_format = 'ASCII'
         data_output_format = 'ASCII'
         model_input_format = 'WS'
         model_output_format = 'WS'
         efield_input_format = 'binary'
         efield_output_format = 'binary'
     /
     &forward
         SFF = .false.
         primary_model = 'file'
         primary_model_file = 'none'
         primary_field = 'file'
         primary_field_file = 'none'
         esoln_output = 'none'
     /
     &spherical_mt
         ! This section is only generated and processed if
         ! ModEM is built with spherical coordinates
         fake_grid = .false.
         fake_center_latlon = 0.0d0, 90.0d0
         fake_data_orientation = .true.
         external_source_file = 'none'

     /


.. |br| raw:: html

   <br />


Namelist Options Reference
-----------------------------

settings Section
~~~~~~~~~~~~~~~~~

.. table :: &settings namelist options 
    :widths: 15 15 7 30

    +-----------------+-----------------------------+-----------+---------------------------------------------+
    | Name            | Values                      | Default   | Description                                 |
    +=================+=============================+===========+=============================================+
    | output_level    | 'debug', 'full', 'regular', | 'regular' | Indicates the desire level of output level  |
    |                 | 'compact', 'result', 'none' |           | to the screen and to files.                 |
    |                 |                             |           | Note: This value will override -v passed on |
    |                 |                             |           | the command line.                           |
    +-----------------+-----------------------------+-----------+---------------------------------------------+
    | covariance_type | 'AR', 'H2'                  | 'AR'      | Covariance type for inversion.              |
    |                 |                             |           | Currently has no effect, see note below.    |
    +-----------------+-----------------------------+-----------+---------------------------------------------+


.. note::

   The `covariance_type` Namelist option currently has no effect. To set it,
   use the 'Covariance backend CovType' option in the inversion control file.
   Additional reference on the covariance type can be found at:
   :ref:`h2-covariance-type`.


forward Section
~~~~~~~~~~~~~~~~


.. table :: &forward namelist options 
    :widths: 15 15 7 30

    +--------------------+-------------------+-----------+---------------------------------------------------------------+
    | Name               | Values            | Default   | Description                                                   |
    +====================+===================+===========+===============================================================+
    | SFF                | .true., .false.   | .false.   | .true. forces the secondary field formulation                 |
    +--------------------+-------------------+-----------+---------------------------------------------------------------+
    | primary_model      | 'file', 'average' | 'file'    | Source of the 1D model.                                       |
    |                    |                   |           | 'file': Read from the specified at primary_model_file.        |
    |                    |                   |           | 'average': 'Avg Description'.                                 |
    |                    |                   |           | Only applicable if 'SFF' is .true.                            |
    +--------------------+-------------------+-----------+---------------------------------------------------------------+
    | primary_model_file | Path to primary   | 'none'    | Path to primary model file to use for SFF.                    |
    |                    | model file        |           | Only applicable if 'SFF' is .true. AND 'primary_model'        |
    |                    |                   |           | is set to 'file'                                              |
    +--------------------+-------------------+-----------+---------------------------------------------------------------+
    | primary_field      | 'file', 'compute' | 'file'    | Source of the primary E-Field: 'file': read from file         |
    |                    |                   |           | specified in primary_field_file. 'compute': compute           |
    |                    |                   |           | the primary field internally.                                 |
    |                    |                   |           | Only applicable if 'SFF' is .ture.                            |
    +--------------------+-------------------+-----------+---------------------------------------------------------------+
    | primary_field_file | Path to primary   | 'none'    | Path to primary E-Field file to use for SFF. Only             |
    |                    | field file        |           | applicable if SFF is true and primary_field_file              |
    |                    |                   |           | is set to 'file'                                              |
    +--------------------+-------------------+-----------+---------------------------------------------------------------+
    | esoln_output       | 'ALL',            | 'ALL'     | EM solution output granularity.                               |
    |                    | 'PER_PERIOD',     |           | 'All': Write all esolns output into one file.                 |
    |                    | 'PER_PERIOD_MODE' |           | 'PER_PERIOD': Write esolns output into a file for each period |
    |                    |                   |           | 'PER_PERIOD_MODE': Write esolns into a file for each period   |
    |                    |                   |           | and mode.                                                     |
    +--------------------+-------------------+-----------+---------------------------------------------------------------+

io Section
~~~~~~~~~~~

.. table :: &io namelist options 
    :widths: 15 15 7 30

    +------------------------+-----------------------------+-----------+----------------------------------------------------+
    | Name                   | Values                      | Default   | Description                                        |
    +========================+=============================+===========+====================================================+
    | data_input_format      | 'ASCII', 'HDF5'             | 'ASCII'   | Input format for data                              |
    +------------------------+-----------------------------+-----------+----------------------------------------------------+
    | data_output_format     | Same as data_input_format   | 'ASCII'   | Output format for data                             |
    +------------------------+-----------------------------+-----------+----------------------------------------------------+
    | model_input_format     | 'WS', 'HDF5', 'binary',     | 'WS'      | Input format for the model.                        |
    |                        |                             |           | 'WS': ModEM version of Weerachai Siripunvaraporn’s.|
    |                        |                             |           | WSINV3DMT model version.                           |
    |                        |                             |           | 'HDF5': Only available if ModEM                    |
    |                        |                             |           | is built with HDF5.                                |
    |                        |                             |           | 'binary': Fortran Binary Format.                   |
    |                        |                             |           | 'RM': Based on Randall Mackie                      |
    |                        |                             |           | model format.                                      |
    +------------------------+-----------------------------+-----------+----------------------------------------------------+
    | model_output_format    | Same as model_input_format  | 'WS'      | Output format for the model.                       |
    |                        |                             |           | Same options as model_input_format                 |
    +------------------------+-----------------------------+-----------+----------------------------------------------------+
    | efield_input_format    | 'binary', 'HDF5', 'ASCII'   | 'binary'  | Input file format for electric field               |
    +------------------------+-----------------------------+-----------+----------------------------------------------------+
    | efield_output_format   | Same as efield_input_format | 'binary'  | Output file format for electric field              |
    +------------------------+-----------------------------+-----------+----------------------------------------------------+
    | jacobian_output_format | 'binary', 'ASCII'           | 'binary'  | Output file format for the Jacobian                |
    +------------------------+-----------------------------+-----------+----------------------------------------------------+

sphereical_mt Section
~~~~~~~~~~~~~~~~~~~~~~

.. note::

   The ``&spherical_mt`` is only generated and processed if ModEM is built with
   spherical coordinates. See :ref:`building_modem_spherical` for more
   information on building ModEM in spherical coordinates.


.. table :: &spherical_mt namelist options 
    :widths: 15 15 7 30

    +------------------------+-----------------------------+---------------+---------------------------------------------+
    | Name                   | Values                      | Default       | Description                                 |
    +========================+=============================+===============+=============================================+
    | fake_grid              | .true., .false.             | .false.       | .false.: Do not apply grid rotation to      |
    |                        |                             |               | equator                                     |
    +------------------------+-----------------------------+---------------+---------------------------------------------+
    | fake_center_latlon     | lat, lon                    | 0.0d0, 90.0d0 | Fake grid center if fake_grid is .true.     |
    +------------------------+-----------------------------+---------------+---------------------------------------------+
    | fake_data_orientation  | .true., .false.             | .true.        | Orient the data to geographic using         |
    |                        |                             |               | fake_center_latlon                          |
    +------------------------+-----------------------------+---------------+---------------------------------------------+


