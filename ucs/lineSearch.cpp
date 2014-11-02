#include "lineSearch.h"
#include <cmath>
#include <iostream>
#include <cstring>
#include <iomanip>

void lineSearch(int ndv, double* x, double* bounds, void(*function)(FUNC_ARGS), 
		void(*gradient)(GRAD_ARGS), int* uiparm, double* urparm)
{
  int maxLineSteps = 10;
  int maxIterations = 10;
  int maxFuncEval = 100;
  double tol = 1.0e-12;
  
  //error function return
  int nf;

  int nFuncEval = 0;
  
  double* grad = new double[ndv];
  double* f = new double[ndv];
  double dq = 999.0;
  double dqlocal;

  //stepsize
  double alpha = 0.2;
  double maxStep = 1.0;

  //get starting location
  function(&ndv, x, &nf, f, uiparm, urparm);
  double forig = *f;
  nFuncEval++;

  //search direction
  double* s = new double[ndv];
  double* dx = new double[ndv];
  double* xold = new double[ndv];

  //print header 
  std::cout << "IT\tNF\tF\t\tdeltaF\t\tSTEP" << std::endl;
  std::cout << std::setprecision(5);
  std::cout << std::scientific;
  
  int lineSteps;
  int nIter = 0;
  double fold = *f;
  //enter the outer convergence loop
  do{
    //get search direction
    gradient(&ndv, x, &nf, grad, uiparm, urparm);
    double smax = 0.0;
    for(int i = 0; i < ndv; i++){
      if(fabs(grad[i]) > fabs(smax)){
	smax = fabs(grad[i]);
      }
    }
    double gradmag = 0.0;
    for(int i = 0; i < ndv; i++){
      s[i] = grad[i] / smax;
      gradmag += s[i]*s[i];
    }
    gradmag = sqrt(gradmag);
    if(nIter == 0){
      //try to make our step size uniform based on maxStep and nstep
      alpha = maxStep / gradmag;
    }
    //enter the line search loop
    lineSteps = 0;
    do{
      std::cout << "\t\t Alpha: " << alpha << std::endl;
      memcpy(xold, x, sizeof(double)*ndv);
      //update change based on gradients and step size
      for(int i = 0; i < ndv; i++){
	dx[i] = -s[i]*alpha;
      }
      //check for updates that move past maxstep and rescale if needed
      double mag = 0.0;
      for(int i = 0; i < ndv; i++){
	mag += dx[i]*dx[i];
      }
      mag = sqrt(mag);
      if(fabs(mag) > maxStep){
	double scale = maxStep / fabs(mag);
	for(int i = 0; i < ndv; i++){
	  dx[i] *= scale;
	}
      }
      std::cout << "\t\tUpdates: ";
      for(int i = 0; i < ndv; i++){
	std::cout << dx[i] << " ";
      }
      //update evaluation location
      for(int i = 0; i < ndv; i++){
	x[i] += dx[i];
      }
      //check the bounds
      for(int i = 0; i < ndv; i++){
	if(x[i] < bounds[i*2 + 0]){
	  x[i] = bounds[i*2 + 0];
	}
	if(x[i] > bounds[i*2 + 1]){
	  x[i] = bounds[i*2 + 1];
	}
      }
      function(&ndv, x, &nf, f, uiparm, urparm);
      //if we increased the function value, reduce the step size and try again
      if(fold < *f){
	alpha = 0.7*alpha;
	memcpy(x, xold, sizeof(double)*ndv);
      }
      //if we did well, increase it and try again
      else{
	dq = fabs(fold - *f);
	fold = *f;
	alpha = 1.2*alpha;
      }
      std::cout << "\t\t" << *f << std::endl;
      nFuncEval++;
      lineSteps++;
    }while(lineSteps < maxLineSteps && (dq >= tol));
    nIter++;
    std::cout << nIter << "\t" << nFuncEval << "\t" << fold << "\t" << dq << "\t" << alpha*smax << std::endl;
  }while((nFuncEval < maxFuncEval) && (nIter < maxIterations) && (dq >= tol));
  

  delete [] s;
  delete [] dx;
  delete [] xold;
  delete [] grad;
  delete [] f;
	   
  return;
}
