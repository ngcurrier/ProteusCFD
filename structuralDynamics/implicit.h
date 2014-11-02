#ifndef IMPLICIT_H__
#define IMPLICIT_H__

#include "bc.h"
#include "forces.h"
#include "utility.h"
#include "matrix.h"

namespace STRUCTDYN{

class BC;
class SForces;

int NewmarkBeta(int dof, int nsteps, double dt, double gamma, double beta, double* m, 
		double* c, double* k, double* x, double* xd, double* xdd, SForces* forces, BC* bcs);


int WilsonTheta(int dof, int nsteps, double dt, double theta, double* m, 
		double* c, double* k, double* x, double* xd, double* xdd, SForces* forces, BC* bcs);

//same routine as above but formulated for netwon iterations, returns delta of variables
int NewmarkBetaDx(int dof, double dt, double gamma, double beta, double* m, 
		  double* c, double* k, double* dx, double* x_n, double* xd_n, double* xdd_n, 
		  double* f, BC* bcs = NULL);

//build keff matrix for NewmarkBeta style integration
void NewmarkBetaKeff(int dof, double dt, double gamma, double beta, double* m, double* c, 
		     double* k, double* keff);

//updates velocity and acceleration given new position and old velocity and acceleration vectors
void NewmarkBetaUpdateVelAcc(int dof, double dt, double gamma, double beta, double* dx, 
			     double* xd_n, double* xdd_n, double* xd_np1, double* xdd_np1);

}
#endif
