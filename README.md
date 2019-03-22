# ModEM-ON
ModEM-ON is a software package for the joint inversion of 3D electromagnetic
data with emphasis on geophysical applications. It is based on the ModEM code
(http://www.modem-geophysics.com/) from Dr. Gary Egbert.

# Summary
We propose to implement workflows for the integration of multiphysics data using
a 3D code for the joint inversion of CSEM, and seismic data. Such workflows will
be built based on extensions to the software ModEM. The developed workflows will
help unveil complex geological structures, e.g. sub salt regions in sedimentary basins.
We expect this new software to be utilized by EM methods experts, both in Brazil
and abroad, to contribute for the understanding of complex geological formations
associated with salt tectonics in sedimentary basins, today one of the greatest
chalenges for Brazil's oil & gas industry.

# Methodology
ModEM-ON will be targeted on joint and cooperative inversion of electromagnetic
data. It will have facilities for the inclusion of seismic and other geophysical
data as prior information in order to better constrain the inversions. The next
sections show some of the main areas where ModEM-ON project will concentrate its
efforts.

## Forward modeling
1.  Improved parallelization
2.  Finalize works on the cur cur equations in order to get better conditioned systems
3.  Conversion of code from MATLAB to Fortran for new interpolation functions
4.  Conversion of code from MATLAB to Fortran for multi resolution finite difference grids
5.  Explore multigrid and preconditioners for Krylov solvers
6.  Implementation of more efficient boundary conditions

## Parameterization and regularization
1.  Interchangeable parameterizations and regularization strategies


## Flexible joint inversion

## User interface