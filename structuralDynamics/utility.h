#ifndef UTILITY_H__
#define UTILITY_H__

#include <cmath>
#include <fstream>
#include <iostream>
#include <string.h>
#include <typeinfo>
#include "forces.h"

namespace STRUCTDYN{

#define PI 3.14159265

class SForces;

int AnalyticSolution(int dof, int nsteps, double dt, double* m, double* c, 
		     double* k, double* x, double* xd, double* xdd, SForces* forces);

void WriteSolution(int dof, int nsteps, double dt, double* x, double* xd, double* xdd);

void Read1DCase(int* dof_, int nsteps, double dt, double** x_, double** xd_, double** xdd_, 
		double** m_, double** c_, double** k_, SForces* forces);

//solve diagonalized mass matrix, with rhs vector vin, returns in vout
void MSolve(double* m, double* vin, double* vout, int ndof);

}
#endif
