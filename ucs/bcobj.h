#ifndef BCOBJ_H__
#define BCOBJ_H__

#include "general.h"
#include "bc_defines.h"
#include <string>
#include <sstream>

template <class Type>
class BCObj
{
 public:

  BCObj();
  ~BCObj();

  //function which will modify Qref according to set dynamic boundary
  void UpdateBCObj(Type tstep);  
  //helper functions so we can return the correct type regardless of what is in BCObj
  void GetQref(Real* Qref);
  void GetQref(RCmplx* Qref);

  //returns >0 if a bc is associated with movement from a design
  //specification, i.e. we know the exact displacements of the boundary node
  bool IsMovingBC();

  void SetBCType(Int setType);                    //sets the numeric type of bc to apply
  void SetName(std::string name){this->name = name;};   //sets the human readable name of the BC
  Int GetBCFactag(){return factag;};              //gets the surface tag

  Int GetBCType(){return bcType;};                //returns the enum type of the BC
  std::string GetName(){return name;};            //returns the name of the boundary condition
  void SetBCFactag(Int factag){this->factag = factag;}; //sets the surface tag
  
  Int neqn;          //number of eqn vars held in qref
  Int nvars;         //number of total vars held in qref
  Type* Qref;        //reference state to evaluate from - may be dynamic
  Bool QrefFromFile; //0 - if qref points to eqnset qinf
                     //1 - if set from file
                     //avoids a double free memory issue

  Type* Qmin;   //minimum state if using periodic bc
  Type* Qmax;   //maximum state if using periodic bc

  /***************************************************************/
  //Variables stored in native (w.r.t. file) state for use in 
  //constructing an appropriate reference state
  /***************************************************************/
 
  Type flowDirection[3];

  Int movement;  //0 - static wall, 1 - noSlip wall is sliding, 2 - noSlip wall is rotating
  Type slipDirection[3];
  Type slipSpeed;        // m/s speed of wall slipping
  Type rotationAxis[3];  //direction vector pointing in axial direction
  Type rotationPoint[3]; //point on the axis
  Type omega;    // rotational speed (i.e. rad/s)

  Type velocity;        // m/s
  Type backPressure;    // pa
  Type density;         // kg/m^3
  Type* massFractions;
  Int nspecies;
  Int movingBC;      //flag to identify computational design active or movement bcs

  Int bleedSteps;    //number of steps to close off boundary condition
                     //used to slowly stop flow on viscous BCs

  Type twall;        // (K) used to specify the kind of viscous wall condition we are using
                     // twall > 0.0 is specified temp. ratio, twall < 0.0 is adiabatic
  Type flux;         // (W/m^2) used for heat conduction object - heat flux

  
 private:

  Int bcType;       //holds an integer which identifies the type of boundary condition which is applied (see bc_defines.h)
  Int factag;       //holds facetag for which the surface can be identified
  std::string name; //holds the name of the bc described by this bcobj
};

#include "bcobj.tcc"

#endif
