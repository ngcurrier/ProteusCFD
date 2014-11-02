#ifndef IO_H__
#define IO_H__

#include <iostream>
#include <fstream>
#include <string>
#include "element_lib.h"

namespace STRUCTDYN{

int WriteVTK_Ascii(std::string casename, int* nelem, Element* elems, int nnode, 
		   double* xyz, double* variables, int nvars, std::string* varNames);

}
#endif
