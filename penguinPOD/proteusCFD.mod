#%Module1.0#####################################################################
##
## dot modulefile
##
##
set name        ProteusCFD
set version     ProteusCFD-git-Linux64 
set path        $env(HOME)/bin/ProteusCFD/


# for Tcl script use only
set     dotversion      1.0.0


## load dependent modules

## apply settings
prepend-path    PATH                     $path/bin
prepend-path    MANPATH                  $path/share/man
