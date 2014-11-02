ProteusCFD
==========

Proteus is a computational fluid dynamics (CFD) solver 
with the goal of providing a complete, parallel, 
multi-physics platform for advanced simulation.

Proteus is distributed without any warranty expressed or
implied for any purpose. 

The author, Nicholas Currier, distributes the software from
the work accumulated during the course of his doctoral 
studies. It is provided free for academic and commercial
use with the following provisions:

1) The software is distributed under the GPL v3 license.
Commercial licenses which are not copyleft may be obtained
by contacting the above author and are provided on a one
time basis after contract negotiation.  This means you CANNOT
bundle our software inside another package without providing
your source as well. This is good for the community and for
science! Please respect these terms.

2) If you find this work useful please cite the author's 
dissertation.
Nicholas Currier, "Reacting Plume Inversion on Urban Geometries through 
Gradient Based Design Methodologies", Ph.D. dissertation,
University of Tennessee at Chattanooga, Chattanooga, Tennessee, August 2014.


Dependencies:
  This solver toolchain has dependencies on MPI, HDF5, METIS, and TINYXML
  They are all included here and their location must be correct in 
  make.opts for compilation to succeed.

Compilation instructions:
  In the root directory type - "make TARGET=local"
  This should work for most installs as long as the above libraries are
  in your execution path. Other options are
  
  make TARGET=bluetick (must also switch toolchain to XLC)
  make TARGET=papertape (with infiniband support)

Usage:
  Proteus (like most CFD solvers) requires a volume grid of a geometry. The grid in this
  case must not contain hanging nodes but is fully unstructured in the typical sense.
  Using Proteus with a given mesh involves:
  
  1) Using ./decomp tool to decompose the geometry into multiple parallel partitions.
     This tool currently expect .ugrid (MSU) and .crunch mesh files. Others may be
     added in the future if there is a request and possibly user support (testers)
     for them.
  
  2) Defining <casename>.bc and <casename>.param files (see ProteusTestSuite for examples)
     This sets up boundary conditions and the solver runtime parameters for relevant physics.
  
  3) Running the solver with mpiexec -np <number of processors> ./ucs <casename>
     There is also a run script in the ./tools directory to modify should you need 
     parallel job scripts.
  
  4) Using ./recomp tool to recompose the parallel geometry and solution to .vtk binary files.
     VTK legacy files are the only supported output at this time.
  
  5) Using paraview or visit (or any other VTK capable visualization tool) to view and query results.

