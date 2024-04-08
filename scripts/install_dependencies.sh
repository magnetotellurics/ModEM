########################################################
# ModEM dependencies installation guide for Windows 11 #
########################################################
#!/bin/bash
#
# Install a Windows Subsystem for Linux (WSL)
# Ex.: Ubuntu 22.04.6 LTS from Microsoft Store
#
#
#Update|Upgrade ubuntu
    sudo apt-get update
    sudo apt-get upgrade
#
#Install make
    sudo apt install make
#
#Install MPI
    sudo apt-get -y install mpich
#
#Install gfortran
    sudo apt-get install gfortran
#
#Install blas and lapack libs
    sudo apt-get install libblas-dev liblapack-dev
#
#Install ftw3 lib
    sudo apt-get install libfftw3-dev libfftw3-doc
#
#update ubuntu
    sudo apt-get update
#