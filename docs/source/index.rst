.. ModEM documentation master file, created by
   sphinx-quickstart on Fri Nov 28 12:31:29 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:tocdepth: 4

.. image:: ./imgs/ModEM_logo_wide_on_white_simple.png
   :alt: ModEM Logo
   :align: center

ModEM Documentation
===================

Welcome to the ModEM's documentation! ModEM is a modular system for inversion of
electromagnetic geophysical data.  

ModEM was designed as a flexible electromagnetic modeling and inversion system,
written in Fortran 95. Although the code can be (and has been) extended for
inversion of more general types of EM data (e.g., controlled source, DC; see
Meqbel and Ritter, 2015), here we describe the stable, core system developed for
2D and 3D magnetotelluric (MT) problems.

While a primary design goal of the system was to allow simplified extension with
regard to data types, modeling codes, parameterization and regularization, and
inversion search algorithms (see Egbert and Kelbert, 2012; Kelbert et al.
2014). The stable MT program has a command-line interface, which controls
available program options, and specifies required and optional input and output
files.

User's Guide
------------

This documentation contains information about compiling, running and using
ModEM, including information about the file formats, options for forward modeling,
options for inversion, etc.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   obtaining_modem
   getting_started
   building_with_cmake
   file_formats
   forward_and_inversion
   beyond_the_basics
   installing_deps
   release_log
   citations
