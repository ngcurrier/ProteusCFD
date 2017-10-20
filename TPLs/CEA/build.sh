#!/bin/bash

# unzip all the relevant directories, we need the libraries so we take exec as well
tar -xvf CEA+Fortran.tar.Z
tar -xvf CEAexec-linux.tar.Z

# build CEA exec
gfortran cea2.f -o FCEA2
