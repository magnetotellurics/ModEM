Change Log
============

Release 1.2.0 - April 12, 2019
-------------------------------

Flexible Air Layers
^^^^^^^^^^^^^^^^^^^^

See: :ref:`flexible-air-layers`.

Log10 Resistivity Model
^^^^^^^^^^^^^^^^^^^^^^^

See: :ref:`log10-res-model`.

Unit Testing For Developers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See: :ref:`unit-testing`.

Other comments
^^^^^^^^^^^^^^

Relative to the previous stable release, the value of ``LARGE_REAL`` that
defines the error bars for modeled data has changed from 2.0e15 to 1.0e13. This
does not affect any computations but it may be confusing to the user as these
values are output in place of the error bars for modeled MT data. Otherwise,
impedances are exactly the same unless the air layers are changed. With modified
air layers values are noticeably different and probably more accurate... because
now the diagonals are exactly zero over a 1D area. Tippers are less noticeably
affected by the air layers modifications.

In addition to these user facing changes, this code has undergone substantial
development, specifically to allow for flexibly boundary conditions in the
future. We have also cleaned up the handling of non-zero sources in
EMsolve3D.f90, specifically checked the physics and the maths, and added
comments to clarify things.

Known Issues
^^^^^^^^^^^^

**[POSSIBLY NOT AN ISSUE]** As of revision 485 on July 24, 2014 , apparent
resistivities in 3D MT are specified in the file as actual resistivity, but are
converted to the log domain before inversion in the program. There is no way
presently to invert apparent resistivity itself, only logs. Probably OK, but
obtuse.

**[POSSIBLY OBSOLETE]** Also, in the stable version that most users have there is a
bug in error propagation for the apparent resistivities.
