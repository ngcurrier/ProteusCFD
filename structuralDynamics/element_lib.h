#ifndef ELEMENT_LIB_H__
#define ELEMENT_LIB_H__

#include <iostream>
#include <cmath>
#include "structparam.h"
#include "utility.h"
#include "geometry.h"
#include "matrix.h"

namespace STRUCTDYN{

class SParam;

class Element
{
public:
  //vector in the local xy plane of the element, used to define global system
  //length is not important
  double vecxy[3];

  //material ID to lookup material properties
  int materialId;

  //list of global node ids
  int* nodes;

  //don't call default element contructor
  //call specific elements instead
  Element();
  virtual ~Element();

  //initialize internal element data
  void Init(int materialId, double* vecxy, int* nodes);

  //compute element mass and element stiffness matrices
  virtual void Compute(double* elm, double* els, bool dynamic, bool gravity, SParam* param)
  {
    std::cerr << "WARING: Compute() not implemented for current element type!!!" << std::endl;
    return;
  };

  //transform element mass and stiffness matrices to global coordinates
  void Transform(double* elm, double* eldiagm, double* els, bool dynamic, SParam* param);

  //place element mass and element stiffness matrices in global mass and stiffness locations
  //if damping is turned on alpha and beta must be specified for proportional damping
  //c = alpha*m + beta*k
  //tNodes is the total # of nodes in the mesh, used for offsets in global matrix
  virtual void Assemble(double* elm, double* eldiagm, double* els, double* gm, 
			double* diagm, double* gs, double* gc, bool damping, 
			bool gravity, double alpha, double beta, SParam* param)
  {
    std::cerr << "WARING: Assemble() not implemented for current element type!!!" << std::endl;
    return;
  };

  //Form diagonal mass matrix via HRZ lumping.
  //     Uses the elements of the consistent mass matrix
  //     that have been scaled so as to preserve the
  //     total element mass.
  //  ...Ref.
  //     Hinton E., Rock T., Zienkiewicz O.C., "A note on mass lumping
  //     and related processes in the finite element method," Earthquake
  //     Engineering and Structural Dynamics, Vol. 4, No. 3, 1976, pp. 245
  //     - 249.
  virtual void LumpHRZ(double* elm, SParam* param)
  {
    std::cerr << "WARING: LumpHRZ() not implemented for current element type!!!" << std::endl;
    return;
  };

  //returns number of dof per node in element
  virtual int GetNodeDOF()
  {
    std::cerr << "WARING: GetNodeDOF() not implemented for current element type!!!" << std::endl;
    return -1;
  };

  //returns number of dof per element
  virtual int GetElemDOF()
  {
    std::cerr << "WARING: GetElemDOF() not implemented for current element type!!!" << std::endl;
    return -1;
  };

  //returns number of nodes in element
  virtual int GetNnodes()
  {
    std::cerr << "WARING: GetNnodes() not implemented for current element type!!!" << std::endl;
    return -1;
  };

  //returns dof for node that corresponds to DX
  virtual int GetDxDOF()
  {
    std::cerr << "WARING: GetDxDOF() not implemented for current element type!!!" << std::endl;
    return -1;
  };
  //returns dof for node that corresponds to DY
  virtual int GetDyDOF()
  {
    std::cerr << "WARING: GetDyDOF() not implemented for current element type!!!" << std::endl;
    return -1;
  };
  //returns dof for node that corresponds to DZ
  virtual int GetDzDOF()
  {
    std::cerr << "WARING: GetDzDOF() not implemented for current element type!!!" << std::endl;
    return -1;
  };
  //returns dof for node that corresponds to RX
  virtual int GetRxDOF()
  {
    std::cerr << "WARING: GetRxDOF() not implemented for current element type!!!" << std::endl;
    return -1;
  };
  //returns dof for node that corresponds to RY
  virtual int GetRyDOF()
  {
    std::cerr << "WARING: GetRyDOF() not implemented for current element type!!!" << std::endl;
    return -1;
  };
  //returns dof for node that corresponds to RZ
  virtual int GetRzDOF()
  {
    std::cerr << "WARING: GetRzDOF() not implemented for current element type!!!" << std::endl;
    return -1;
  };

  //returns an interpolated values for an point which lies on the element
  //this interpolation is naturally different for each element type
  //sol - solution vector for the given mesh
  //xyzPt - point on element
  //param - parameters object to store solution info
  //ndof - returns ndof for the solution at this point
  //values - returns ndof values of the interpolated solution at point
  virtual void InterpolatePoint(double* xyzPt, SParam* param, double* sol, int* ndof, 
				double* values)
  {
    std::cerr << "WARNING: InterpolatePoint not implemented for current element type!!!" 
	      << std::endl;
    return;
  };

  virtual void InterpolateValuesToNodes(double* xyzPt, SParam* param, double* values, 
					int nval, double* rhs)
  {
    std::cerr << "WARNING: InterpolateForcesToNodes not implemented for current element type!!!" 
	      << std::endl;
    return;
  };

private:
};


//
// Each Beam element consists of 2 nodes with 6 DOF each
//
// The DOF for each node are:
// ========================================
// 1 - x translation 
// 2 - y translation 
// 3 - z translation 
// 3 - rotation about x axis
// 4 - rotation about y axis
// 5 - rotation about z axis
//
class Beam : public Element
{
public:
  Beam();

  //compute element mass and element stiffness matrices
  void Compute(double* elm, double* els, bool dynamic, bool gravity, SParam* param);

  //place element mass and element stiffness matrices in global mass and stiffness locations
  //if damping is turned on alpha and beta must be specified for proportional damping
  //c = alpha*m + beta*k
  void Assemble(double* elm, double* eldiagm, double* els, double* gm, double* diagm, 
		double* gs, double* gc, bool damping, bool gravity, double alpha, 
		double beta, SParam* param);

  //compute lumped element matrix using HRZ lumping
  void LumpHRZ(double* elm, SParam* param);


  //returns number of dof per node in element
  int GetNodeDOF()
  {
    return (6);
  };
  //returns number of dof per element
  int GetElemDOF()
  {
    return (12);
  };

  //returns number of nodes in element
  int GetNnodes()
  {
    return (2);
  };

  //returns dof for node that corresponds to DX
  int GetDxDOF()
  {
    return (0);
  };
  //returns dof for node that corresponds to DY
  int GetDyDOF()
  {
    return (1);
  };
  //returns dof for node that corresponds to DZ
  int GetDzDOF()
  {
    return (2);
  };
  //returns dof for node that corresponds to RX
  int GetRxDOF()
  {
    return (3);
  };
  //returns dof for node that corresponds to RY
  int GetRyDOF()
  {
    return (4);
  };
  //returns dof for node that corresponds to RZ
  int GetRzDOF()
  {
    return (5);
  };

  virtual void InterpolatePoint(double* xyzPt, SParam* param, double* sol, int* ndof, 
				double* values);

  virtual void InterpolateValuesToNodes(double* xyzPt, SParam* param, double* values, 
					int nval, double* rhs);
private:

};

}
#endif
