#ifndef BCOBJ_H__
#define BCOBJ_H__

#include "general.h"

template <class Type>
class BCObj
{
 public:

  BCObj();
  ~BCObj();

  
  Int neqn;          //number of eqn vars held in qref
  Int nvars;         //number of total vars held in qref
  Type* Qref;        //reference state to evaluate from - may be dynamic
  Bool QrefFromFile; //0 - if qref points to eqnset qinf
                     //1 - if set from file
                     //avoids a double free memory issue

  Int periodic; //flag to set periodic bc 0 - off, 1 - on

  Type* Qmin;   //minimum state if using periodic bc
  Type* Qmax;   //maximum state if using periodic bc
  
  Type period;  //non-dimensional time over which a full wave cycles

  /***************************************************************/
  //Variables stored in native (w.r.t. file) state for use in 
  //constructing an appropriate reference state
  /***************************************************************/
 
  Type flowDirection[3];

  Int movement;  //0 - static wall, 1 - noSlip wall is sliding, 2 - noSlip wall is rotating
  Type slipDirection[3];
  Type slipSpeed;        //non-dimensional speed of wall slipping
  Type rotationAxis[3];  //direction vector pointing in axial direction
  Type rotationPoint[3]; //point on the axis
  Type omega;    //non-dimensional rotational speed (i.e. rad/s)

  Type velocity;
  Type backPressure;
  Type density;
  Type* massFractions;
  Int nspecies;
  Int factag;        //holds facetag for which the surface can be identified
  Int movingBC;      //flag to identify computational design active or movement bcs

  Int bleedSteps;    //number of steps to close off boundary condition
                     //used to slowly stop flow on viscous BCs

  Type twall;        //used to specify the kind of viscous wall condition we are using
                     //twall > 0.0 is specified temp. ratio, twall < 0.0 is adiabatic

  //function which will modify Qref according to set dynamic boundary
  void UpdateBCObj(Type tstep);  
  //helper functions so we can return the correct type regardless of what is in BCObj
  void GetQref(Real* Qref);
  void GetQref(RCmplx* Qref);

  //returns >0 if a bc is associated with movement from a design
  //specification, i.e. we know the exact displacements of the boundary node
  Int IsMovingBC();
  
 private:
};

#include "bcobj.tcc"

#endif
