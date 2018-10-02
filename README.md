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


Compilation targets:
  In the root directory type - "make TARGET=local"
  This should work for most installs as long as the above libraries are
  in your execution path. Other options are
  
  * make TARGET=bluetick (must also switch toolchain to XLC)
  * make TARGET=papertape (with infiniband support)
  
  Also, you must build the decomposition/recomposition tools to write files Proteus understands.
  "make tools"  will build these for you. All executable code will be built in the ProteusCFD/bin directory. All executables are meant to 
  be run via command line and will output usage information if you run the executable with no arguments.

Installation instructions:
  * Install dependencies for OS (list below) to your machine - openMPI is the big one (you need the devel package with headers).
  * Build ucs.x, udecomp.x, urecomp.x and tools required to run
  * make TARGET=local all (on some machines you may have to try this several times to complete all tasks)
      * This also builds the chemical database. You can build it manually if required by.
          * cd chemdata
          * ./chemDB.py
          * cp chemdb.hdf5 ~/.proteusCFD/database/chemdb.hdf5

  * Check that the tests pass after the installation
      * cd bin
      * ./tests.x

Usage:
  Proteus (like most CFD solvers) requires a volume grid of a geometry. The grid in this
  case must not contain hanging nodes but is fully unstructured in the typical sense.
  We recommend Salome or GMSH as an opensource alternative to commercial gridding tools.
  Using Proteus with a given mesh involves:
  
  1) Use udecomp.x <gridfile - e.g. test.ugrid> <number of processors/cores> tool to 
     decompose the geometry into multiple parallel partitions. This tool currently expects .ugrid (MSU), 
     .su2 (Stanford SU2), .msh (GMSH), .cgns (CGNS - i.e. ANSYS Meshing) or .crunch mesh files. Other formats may be added in the future 
     if there is a request and possibly user support (testers) for them. The proteus test suite has several examples. 
     You can check it out by typing at the command line:
     * git clone https://github.com/ngcurrier/ProteusTestSuite/
     * A suggested starting point is the 15DegreeRamp test case. It runs quickly and should get you started easily.
       You can run it by:
         * cd ./ProteusTestSuite
         * Create a symlink to the ProteusCFD tools as below. Alternatively, you can add ProteusCFD/bin to your .bashrc PATH.
             * ln -s ../ProteusCFD/bin/ucs.x ucs.x
             * ln -s ../ProteusCFD/bin/udecomp.x udecomp.x
             * ln -s ../ProteusCFD/bin/urecomp.x urecomp.x
         * udecomp.x 15degramp.crunch 8
         * mpiexec -np 8 ucs.x 15degramp
         * -- Code will run and finish --
         * urecomp.x 15degramp
         * -- Suggested to use paraview to open the 15degramp.vtk file for visualization --         
  
  2) Define <casename>.bc and <casename>.param files (see ProteusTestSuite for examples)
     This sets up boundary conditions and the solver runtime parameters for relevant physics.
     You can get all of the options/flags for the .param file that Proteus currently accepts by typing ucs.x
     with no arguments.  ProteusCFD will write all options to the terminal. Options available for the 
     .bc file can be found in /docs/master.bc.
     * ProteusCFD has the option of defining user based boundary conditions, mesh movement, and initial conditions. Examples
       of these scripts can be found in the code source tree under the pythonScripts folder. Placing them in the runtime
       folder activates these capabilities in turn.
  
  3) Run the solver with mpiexec -np <number of processors> ucs.x <casename>
     There is also a run script in the ./tools directory to modify should you need 
     parallel job scripts.
  
  4) Use urecomp.x <casename -- e.g. test> tool to recompose the parallel geometry and solution to .vtk binary files.
     VTK legacy files are the only supported output at this time.
  
  5) Use paraview or visit (or any other VTK capable visualization tool) to view and query results.

Dependencies included with repository
=====================================

Proteus CFD has several depencies. An effort has been made to rely only on the smallest subset
of required packages in order to make the compilation of this software simple for new users.
The following packages are required for core functionality (and are distributed with the tool):

* HDF5 (all output files are stored in this format, included in TPLs directory)
* METIS (mesh partitioning for parallel runs, included in TPLs directory)
* TINYXML (used for I/O in some places)
* GTest (used for unit testing framework, included in TPLs directory)

Dependencies for OS (you need to install these on the system prior to building)
====================
* Cmake for building Metis
* Make for building ProteusCFD and tools
* GNU/C/C++
* GNU FORTRAN (gfortran)
* MPI (version dependent on local machine infrastructure, we suggest openMPI if you get a choice) - requires dev version headers
* Dev version of python for inclusion of python headers libpython-dev
* H5PY (package python-h5py on Debian) - needed for chemistry database
* Package libz-dev on Debian
* Packages for HDF5 > 1.10: libhdf5-dev AND libhdf5-serial-dev on Debian (you can run script hdf5_fix.sh as root if you are getting linking errors due to serial/parallel HDF5 package installs)
* Package for CGNS support: libcgns-dev on Debian

Dependencies for GUI
====================

Proteus CFD has not traditionally had a GUI. However, development effort is being made to make this a
reality as it is a prerequisite for many users.  This is experimental at this point only.  The following
packages are required:

* pyqt4 (Basis for all widgets and GUI framework)
* vtk with python bindings (Framework for 3D visualization capability) - VTK must be built with the following enabled BUILD_SHARED_LIBS and VTK_WRAP_PYTHON 
* h5py for reading partitioned mesh files
