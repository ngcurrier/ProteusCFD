#ifndef FLUID_STRUCTURE_H__
#define FLUID_STRUCTURE_H__

#include <cmath>
#include "mem_util.h"
#include "bc.h"
#include "implicit.h"
#include "element_lib.h"
#include "structparam.h"
#include "solutionSpaceBase.h"
#include "../ucs/param.h"

#include <string>

namespace STRUCTDYN{

//returns the stiffness(k) and mass(m) (2x2) matrices and rhs for typical section
//airfoil flutter analysis - this is dimensional in input
//when this system is solved the first DOF is h/b (plunge) and the second DOF is alpha (rotation)
class TypicalSection : public SolutionSpaceBase<double> 
{

public:	
  struct SectionMovement{
    double alpha;
    double h;
  };

  Param<double>* param;
  int dof;

  //GIVEN
  double* k;      //k - stiffness matrix (non-dimensional)
  double* m;      //m - mass matrix (non-dimensional)
  double* c;      //c - damping matrix (NOT CURRENTLY UTILIZED, i.e. 0!)
  double* rhs;    //rhs - returns rhs of system (non-dimensional)
  double ma;      //ma - mass of wing (kg)
  double chord;       //c - wing chord (m)
  double kalpha;  //kalpha - torsional stiffness of wing (N/m)
  double kh;      //kh - bending stiffness of wing (N.m/rad)
  double xalpha;  //xalpha - dist. from C.G. to elastic axis (non-dimensional by c/2)
  double icg;     //icg - mass moment of inertia of airfoil about C.G.
  double cl;      //cl - lift coefficient - nondimensional
  double cm;      //cm - moment coefficient at elastic axis (normally 1/4 chord) - nondimensional
  double rhoinf;  //rhoinf - freestream density
  double uinf;    //uinf - freestream velocity

  //DERIVED
  double iea;          //mass moment of inertia about elastic axis (parallel axis thm.)
  double omega_alpha;  //natural frequency rotation
  double b;            //semi-chord

  //SOLUTION
  double* dx;    //update to position which is most recent
  double* x_n;   //position at (n) - previous time level
  double* xd_n;  //velocity at (n) - previous time level
  double* xdd_n; //acceleration at (n) - previous time level
  
  //implicit solution parameters for average acceleration
  double dt;
  double gamma;
  double beta;

  TypicalSection(Param<double>* param, std::string name, TemporalControl<double>& temporalControl);
  ~TypicalSection();

  //print the current solution state
  void Print(std::ostream& sout);
  //forced oscillation function, returns dx, position, velocity, and acceleration
  //given a particular time step
  void ForcedMovement(double* dx, double* x_old, double* x, double* xd, double* xdd);

  //set ICs and solution parameters
  void PreIterate();
  void PreTimeAdvance(){};
  //this is the core of the Newton loop
  void NewtonIterate();
  //update solution storage for next physical timestep, i.e. next Newton iteration
  void PostTimeAdvance();
  void WriteSolution(){};
  void WriteRestartFile();
  void ReadRestartFile();

  void PrintTimers(){};
  void WriteAvailableFields() const{};

private:
  //set the initial conditions
  void IC(double* x, double* xd);
  //build the matrices [m], [k], and [c]
  void BuildSystem();
  //update the RHS of the equation for a given lift/moment coefficient
  void BuildRHS();


};

class Structure : public SolutionSpaceBase<double>
{
public:

  Structure(Param<double>* param, std::string name, TemporalControl<double>& temporalControl);
  ~Structure();

  //this call results in mesh, bcs, forces, etc. being read from
  //file, also calls param->init()
  void Init();

  //dump the needed restart variables
  void WriteRestart(std::ofstream& fout);
  //read the needed restart variables
  void ReadRestart(std::ifstream& fin);
  //compute the updated solution
  void ComputeDx();
  //update solution storage for next physical timestep, i.e. next Newton iteration
  void NextTimestep();

  //build and store the field connectivity to a CFD (wetted) surface
  void BuildFieldConnectivity(double* xyz, int npts);
  //returns the point displacements for the CFD solver to move the wetted surface with
  void GetPointsDisplacements(double* dxyz);
  //return the element id which is closest to the point given
  //also returns a vector which goes from the closest point on the element
  //to the point in question
  int GetClosestElement(double* xyz, double* radVector);
  //apply the force vector given at each of the CFD (wetted) surface points to FEA model
  void ApplyCFDForce(double* force);

  void PrintTimers(){};
  void WriteAvailableFields()const {};

  double* k;      //k - stiffness matrix (dimensional)
  double* m;      //m - mass matrix (dimensional)
  double* c;      //c - damping matrix (NOT CURRENTLY UTILIZED, i.e. 0!)
  double* rhs;    //rhs - returns rhs of system (non-dimensional)

  //structural parameters file, holds mesh, etc.
  Param<double>* param;
  STRUCTDYN::SParam* sparam;
  //holds boundary conditions for structural solver
  STRUCTDYN::BC* bc;

  //SOLUTION
  double* dx;    //update to position which is most recent
  double* x_n;   //position at (n) - previous time level
  double* xd_n;  //velocity at (n) - previous time level
  double* xdd_n; //acceleration at (n) - previous time level
  double* x_np1; //position at (n+1) - most recent time level

  //connectivity data
  bool connect_alloc;
  int nwettedpts;  //number of points on wetted CFD surface
  double* rays;    //rays which go from FEA point to CFD point
  double* feaxyz;  //xyz location of cfd interpolated data
  int* elem;       //element which contains CFD point

  int dof;  //degrees of freedom

  //implicit solution parameters for average acceleration
  double dt;
  double gamma;
  double beta;

  //These are the required member functions for the solutionspace
  void PreIterate(){};
  void PreTimeAdvance(){};
  void NewtonIterate(){};
  void PostTimeAdvance(){};
  void WriteSolution(){};
  void WriteRestartFile(){};
  void ReadRestartFile(){};

private:

};

}
#endif


