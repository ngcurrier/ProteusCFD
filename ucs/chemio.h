#ifndef CHEM_IO__
#define CHEM_IO__

#include <iostream>
#include <fstream>
#include <string>
#include "general.h"

Int ReadGrimechData(std::string filename, std::string* species, Int* stoichCoeffLeft, Int* stoicCoeffRight, Int nspecies);

#endif
