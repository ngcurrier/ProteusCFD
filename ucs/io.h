#ifndef IO_H__
#define IO_H__

#include "general.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>

//forward declarations
template <class Type> class Mesh;
template <class Type> class PObj;
template <class Type> class SolutionField;

Int ReadUGRID_Ascii(Mesh<Real> &m, std::string filename);
Int ReadCRUNCH_Ascii(Mesh<Real> &m, std::string filename);

Int WriteCRUNCH_Ascii(Mesh<Real> &m, std::string casename);
Int WriteVTK_Ascii(Mesh<Real> &m, std::string casename, std::vector<SolutionField<Real>*>& fields);
Int WriteVTK_Binary(Mesh<Real> &m, std::string casename, std::vector<SolutionField<Real>*>& fields);

void TranslateWinding(Int* nodes, Int translation[6][8], Int num_nodes, Int etype, Int to_other_format);

//filebase is the name minus the extensions... extensions are handled internally
//this file format can be written in parallel and requires a parllel object handle
Int WriteGridXDMF(Mesh<Real> &m, PObj<Real> & p, std::string filebase, std::string meshname="mesh");

#endif
