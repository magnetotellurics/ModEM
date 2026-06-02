FROM ubuntu:latest

RUN apt-get update
RUN apt-get -y install gcc
RUN apt-get -y install g++
RUN apt-get -y install gfortran
RUN apt-get -y install gdb
RUN apt-get -y install valgrind
RUN apt-get -y install libblas-dev
RUN apt-get -y install liblapack-dev
RUN apt-get -y install make
RUN apt-get -y install cmake
RUN apt-get -y install wget
RUN apt-get -y install python3
RUN apt-get -y install python3-pip
RUN apt-get -y install pipx
RUN apt-get -y install vim
RUN apt-get -y install git
RUN apt-get -y install openmpi-bin openmpi-common libopenmpi-dev

# Create a ModEM user and give them sudo access
# so it's easier for people to install things
RUN apt-get -y install sudo
RUN useradd -m -r -g users -s /bin/bash modem
RUN adduser modem sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

ENV LD_LIBRARY_PATH=/usr/local/lib/
# For Python Installs
ENV PATH="$PATH:/root/.local/bin"

COPY . /home/modem/ModEM
RUN chown -R modem /home/modem/ModEM
RUN chgrp -R users /home/modem/ModEM

USER modem
WORKDIR /home/modem/ModEM/f90

# Make MF, SP and SP2 using Configure
RUN ./CONFIG/configure -m mpi -l MF -g release Makefile gfortran && make && mv Mod3DMT Mod3DMT_MF && make clean
RUN ./CONFIG/configure -m mpi -l SP -g release Makefile gfortran && make && mv Mod3DMT Mod3DMT_SP && make clean
RUN ./CONFIG/configure -m mpi -l SP2 -g release Makefile gfortran && make && mv Mod3DMT Mod3DMT_SP2 && make clean

# Make MF, SP and SP2 Using CMake
RUN mkdir -p /home/modem/modem-build
WORKDIR /home/modem/modem-build
RUN cmake ~/ModEM -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DFORWARD_FORMULATION=MF && make
RUN cmake ~/ModEM -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DFORWARD_FORMULATION=SP && make
RUN cmake ~/ModEM -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DFORWARD_FORMULATION=SP2 && make

WORKDIR /home/modem/
