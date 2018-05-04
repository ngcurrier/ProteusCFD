#ifndef PORT_FUNCTIONS_H__
#define PORT_FUNCTIONS_H__

#include "portGlobal.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>

#define GRAD_METHOD 2
#define MAX_BUFFER 1024

//NOTE: if function or grad cannot be evaluated at x the change nf to zero
//nf - error return value
//ndv - number of design variables
//x - vector at which f(x) should be computed
//grad - gradient vector
//uiparm - integer parameters 
//urparm - real parameters
#define GRAD_ARGS int* ndv, double* x, int* nf, double* grad, int* uiparm, double* urparm
#define FUNC_ARGS int* ndv, double* x, int* nf, double* f, int* uiparm, double* urparm

void function_(FUNC_ARGS);
void gradient_(GRAD_ARGS);

int GetNptsPointFile(std::string filename);
void ReadPointMoveFile(std::string filename, int npts, int ndv, 
		       int* pts, double* xyz, double* dxyz, int* betaid);
void WritePointMoveFile(std::string filename, int npts, int ndv, 
			int* pts, double* xyz, double* dxyz, int* betaid);

//calls to interface with solver
void ExternComputeObjectiveFunction(int np);
void ExternMoveMesh(int np);
void ExternSaveOriginalMesh(int np);
void ExternSaveOriginalDesignFile();
void ExternCopyBasePartitionFilesDown(int np);
void ExternComputeMeshSensitivities(int np);
void ExternComputeGradient(int np, int method = 2);

//call to execute a command and get the response returned in a string
std::string GetOutputFromCommand(std::string cmd);
//will make a remote execution call and then wait to return until the job completes
void ExecuteRemoteAndWait(std::string cmd);

#endif
