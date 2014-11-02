#ifndef SENSORS_H__
#define SENSORS_H__

#include "general.h"
#include "parallel.h"
#include "geometry.h"
#include "interpolation.h"
#include <string>
#include <fstream>
#include <vector>
#include <sstream>

template <class Type>
class Sensors
{
public:
   //reads sensor file and performs initial search for each sensor location 
  //initializes internal data for lookups later
  Sensors(std::string filename, Int neqn, Int nstride, Type* Q, 
	  PObj<Type>* p, Int* ipsp, Int* psp, Type* mxyz, Int nnodes);
  ~Sensors();
  
  //performs initial search through mesh for lookup data on nodes
  //this may have to be called repetitively if the mesh is moving
  void Search();

  //calculates each sensor value and stores results internally
  void Sample();

  //write sensor values at last sampling interval to stdout
  //each sensor has its own file
  Int Write(std::string filename, Int stepTag);

  //writes sensor values at last sampling interval to stdout
  void Print();

  //tolerance for bounding box
  Real tol;

  //flag for memory init
  Bool init;

  //number of sensors
  Int count;

  //number of equations we are sampling values from
  Int neqn;

  //the size of the stride between nodal data
  Int nstride;

  //psp array -- needed if inverse distance weighting is used
  Int* psp;

  //psp crs array -- need if inverse distance weighting is used
  Int* ipsp;

  //xyz locations of all the points in the mesh
  Type* mxyz;
  
  //number of nodes in the mesh for above list
  Int mnnode;

  //the array we are sampling data from
  Type* q;

  //list of xyz location of each sensor, linear array
  Type* xyz;

  //values of sensor at last sampling interval
  Type* values;

  //list of the domains each sensor is in after first search
  Int* domain; 

  //list of local node each sensor is closest to
  Int* nodes;

  //parallel object used for searches
  PObj<Type>* p;

private:
  
  Int ReadFile(std::string filename);

};

//include implementations
#include "sensors.tcc"

#endif
