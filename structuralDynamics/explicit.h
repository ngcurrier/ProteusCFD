#ifndef EXPLICIT_H__
#define EXPLICIT_H__

#include "utility.h"
#include "forces.h"
#include "bc.h"
#include "matrix.h"

namespace STRUCTDYN{

//mode == 0 Raulston's
//mode == 1 Trapezoidal rule
int RK2ndOrder(int dof, int nsteps, double dt, double* m, double* c, double* k,
	       double* x, double* xd, double* xdd, SForces* forces, BC* bcs, int mode);

int RK3rdOrder(int dof, int nsteps, double dt, double* m, double* c, double* k,
	       double* x, double* xd, double* xdd, SForces* forces, BC* bcs);

int RK4thOrder(int dof, int nsteps, double dt, double* m, double* c, double* k,
	       double* x, double* xd, double* xdd, SForces* forces, BC* bcs);

//This method requires a diagonal mass and damping matrix or the cost is the
//same as an implicit method, diagonal mass and damping assumed internally
int CentralDiff(int dof, int nsteps, double dt, double* m, double* diagm, double* c, 
		double* k, double* x, double* xd, double* xdd, SForces* forces, BC* bcs);

//This method alleviates the need for a diagonal damping matrix, we just need a
//diagonal mass matrix
int HalfStepCentral(int dof, int nsteps, double dt, double* m, double* diagm, double* c, 
		    double* k, double* x, double* xd, double* xdd, SForces* forces, BC* bcs);

//used for RK methods
void GetAccelFunction(double* F, double t, double* m, double* p, double* c, double* k, 
		      double* x, double* xd, SForces* forces, BC* bcs, int dof);


}
#endif
