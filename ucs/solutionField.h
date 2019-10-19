#ifndef SOLUTIONFIELD_H__
#define SOLUTIONFIELD_H__

#include "general.h"
#include "dataInfo.h"
#include "exceptions.h"
#include "h5layer.h"
#include <string>
#ifdef _OPENMP
#include <pthread.h>
#include <omp.h>
#endif

template <class Type> class Mesh;
template <class Type> class PObj;

//states control how the field behaves in time
namespace FIELDS{
  enum statetype{
    STATE_TIME,
    STATE_NONE,
    STATE_NP1,
    STATE_N,
    STATE_NM1, 
    VAR_INTERIOR,
    VAR_EVERYWHERE,
  };
}

template <class Type>
class SolutionField
{
public:
  template <class Type2>
  SolutionField<Type> & operator= (const SolutionField<Type2>& fieldToCopy);

  //standard construtor using a mesh and datainfo, etc.
  SolutionField(Mesh<Type>& mesh, DataInfo<Type> dataInfo, Int stateType, Int varLocation);
  //constructor using and HDF5 dataset and a mesh, useful for recomp
  SolutionField(Mesh<Type>& mesh, hid_t fileId, std::string HDFPath); 

  ~SolutionField();

  //returns a pointer to the data
  Int GetNdof() const;
  Type* GetData(Int state);
  const Type* GetData(Int state) const;

  std::string GetDofName(Int dof) const;
  Bool DofIsVector(Int dof) const;
  Bool DofIsScalar(Int dof) const;

  Type GetMax() const;
  void Normalize();

  //will copy down values as appropriate if the state
  //value is temporal, i.e np1 to n and n to nm1
  void EvolveInTime();
  
  //get the name of the field
  std::string GetName() const;
  //check if the name matches
  Bool IsNamed(std::string testname) const;
  //check if the field has temporal data (i.e. needs to be preserved for restart)
  Bool IsTemporal() const;

  //Reorders based on a new ordering integer list
  void Reorder(Int* newOrder);
  void Fill(Type value);

  void WriteH5(std::string filename, std::string directory, Int states = FIELDS::STATE_NONE, Bool dimensionalize = false);
  void ReadH5(std::string filename, std::string directory, Int states = FIELDS::STATE_NONE);
  
  //get parallel object pointer for the field
  PObj<Type>* GetPObj();

private:
  SolutionField();    //DO NOT USE

  Mesh<Type> & mesh;  //reference to spatial field (mesh) to which data is tied
  std::string name;   //name of field
  Int ndof;           //number of values (data) stored at each spatial location
  DataInfo<Type> dataInfo;  //data info descriptor
  Int mystate;        //FIELDS:: enum which determines whether or not the field is time aware, etc.
  Type * data;        //raw flat array of data
  Int varLoc;         //Enum value which indicates internal, boundaries, or all locations are stored
  Int nentities;      //number of spatial locations at which data is stored

#ifdef _OPENMP
  pthread_mutex_t* mutexes;  //declare global mutexes for each cv
#endif
};


//include implementation
#include "solutionField.tcc"

#endif
