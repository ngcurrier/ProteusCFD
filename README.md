ProteusCFD
==========

<a href="https://scan.coverity.com/projects/ngcurrier-proteuscfd">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/10218/badge.svg"/>
</a>

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


Compilation instructions:
  In the root directory type - "make TARGET=local"
  This should work for most installs as long as the above libraries are
  in your execution path. Other options are
  
  make TARGET=bluetick (must also switch toolchain to XLC)
  make TARGET=papertape (with infiniband support)
  
  Also, you must build the decomposition/recomposition tools to write files Proteus understands.
  "make tools"  will build these for you
  All executable code will be built in the ProteusCFD/bin directory. All executables are meant to 
  be run via command line and will output usage information if you run the executable with no arguments.

Usage:
  Proteus (like most CFD solvers) requires a volume grid of a geometry. The grid in this
  case must not contain hanging nodes but is fully unstructured in the typical sense.
  We recommend Salome or GMSH as an opensource alternative to commercial gridding tools.
  Using Proteus with a given mesh involves:
  
  1) Using ./decomp.x tool to decompose the geometry into multiple parallel partitions.
     This tool currently expect .ugrid (MSU), .su2 (Stanford SU2) and .crunch mesh files. 
     Others may be added in the future if there is a request and possibly user support (testers)
     for them.
  
  2) Defining <casename>.bc and <casename>.param files (see ProteusTestSuite for examples)
     This sets up boundary conditions and the solver runtime parameters for relevant physics.
  
  3) Run the solver with mpiexec -np <number of processors> ./ucs.x <casename>
     There is also a run script in the ./tools directory to modify should you need 
     parallel job scripts.
  
  4) Using ./recomp.x tool to recompose the parallel geometry and solution to .vtk binary files.
     VTK legacy files are the only supported output at this time.
  
  5) Using paraview or visit (or any other VTK capable visualization tool) to view and query results.

Dependencies
============

Proteus CFD has several depencies. An effort has been made to rely only on the smallest subset
of required packages in order to make the compilation of this software simple for new users.
The following packages are required for core functionality:

* MPI (version dependent on local machine infrastructure, we suggest openMPI if you get a choice)
* HDF5 (all output files are stored in this format, included in TPLs directory)
* METIS (mesh partitioning for parallel runs, included in TPLs directory)
* TINYXML (used for I/O in some places)
* GTest (used for unit testing framework, included in TPLs directory)

Dependencies for GUI
====================

Proteus CFD has not traditionally had a GUI. However, development effort is being made to make this a
reality as it is a prerequisite for many users.  This is experimental at this point only.  The following
packages are required:

* pyqt4 (Basis for all widgets and GUI framework)
* vtk with python bindings (Framework for 3D visualization capability) - VTK must be built with the following enabled BUILD_SHARED_LIBS and VTK_WRAP_PYTHON 
* h5py for reading partitioned mesh files
