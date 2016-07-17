#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

//example: demonstration of a zero point function
class BAR{
public:
  BAR(int a){this->a = a;};
  int a;
  double f(double x){return 3*x*x + x - (double)a;};
};

//examples: demonstration of a zero point function
class FOO{
public:
  FOO(){};
  double f(double x, double a, double b, double c){return a*x*x + b*x - c;};
};

//function which takes four doubles as arguments
template <class Type>
class Functor{
  typedef double (Type::*function_type)(double, double, double, double);
public:
  Functor(Type& classRef, function_type function, double constVarA, double constVarB, double constVarC):
    classRef(classRef), function(function), constA(constVarA), constB(constVarB), constC(constVarC)
  { };
  ~Functor(){};
  double wrapper(double x)
 {
    return (classRef.*function)(x, constA, constB, constC);
 };
  
private:
  double constA;
  double constB;
  double constC;
  Type& classRef;
  function_type function;
};


// this function utilizes a Newton's (secant) method with a trust region limiter and
// should be general enough to find roots for most single variable functions.
// You may have to create a dummy functor class to create the correct function prototype to use this
// example of a functor provided in this file --NGC
template <class Type>
double rootFinder(Type& owningClass, double (Type::*fixedPointFunction)(double), double guess, double tol, bool verbose = false)
{
 
  //TUNEABLE PARAMS
  int maxit = 130;
  double h = 1.0e-4;
 
  //do not touch here
  double root = guess;
  double trustRegion = root/5.0;
  double f0 = 9e16;
  double lastGoodF0 = 9e16;

  if(verbose){
    std::cout << "GUESS: " << guess << std::endl;
  }
  
  int it = 0;
  do{
    //compute derivative
    f0 = (owningClass.*fixedPointFunction)(root);
    double f1 = (owningClass.*fixedPointFunction)(root+h);
    double deriv = (f1-f0)/h;
    //silence the output
    if(verbose){
      std::cout << "f: " << f0 << "\tx: " << root << std::endl;
    }
    double step = -f0/deriv;
    double sgn = 1.0;
    if(step < 0.0) sgn = -1.0;
    root += sgn*std::min(fabs(step), fabs(trustRegion));
    //check for improved step, if improved increase trust region
    if(fabs(f0) < fabs(lastGoodF0)){
      lastGoodF0 = f0;
      trustRegion *= 1.5;
    }
    else{
      trustRegion *= 0.7;
    }
    ++it;
  }while(fabs(f0) > tol && it < maxit);
 
  if(it >= maxit) throw std::runtime_error("WARNING: rootFinder() did not converge in maximum number of iterations");
 
  return root;
};
