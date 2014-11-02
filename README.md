Dependencies:
  This solver toolchain has dependencies on HDF5, METIS, and TINYXML
  They are all included here and their location must be correct in 
  make.opts for compilation to succeed.

Compilation instructions:
  In the root directory type - make TARGET=local
  This should work for most installs as long as the above libraries are
  in your execution path. Other options are
  
  make TARGET=bluetick (must also switch toolchain to XLC)
  make TARGET=papertape (with infiniband support)
