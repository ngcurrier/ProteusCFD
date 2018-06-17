#ifndef SOLUTION_SPACE_BASE_H__
#define SOLUTION_SPACE_BASE_H__

#include "general.h"
#include "timer.h"
#include "temporalControl.h"
#include <vector>

template <class Type>
class ScalarField
{
public:
  ScalarField(std::string name) :
    name(name), scalar(0.0)
  { };
  void SetField(Type value)
  { scalar = value; };
  Type GetField()
  { return scalar; };
  Bool IsNamed(std::string testname)
  {
    if(testname == name){
      return true;
    }
    return false;
  };
  std::string GetName()
  { return name;};
private:
  std::string name;
  Type scalar;
};

template <class Type>
class SolutionSpaceBase
{
public:
  std::string name;
  Int nSteps;
  Int iter;  //keeps track of the current unsteady or steady iteration
  Type time; //if solver is time accurate, this keeps track of the locally dimensioned/non-dimensioned time state
  TemporalControl<Real>& temporalControl;

  SolutionSpaceBase(std::string name, TemporalControl<Real>& temporalControl) :
  name(name), temporalControl(temporalControl), iter(0), time(0.0)
  {};
  virtual ~SolutionSpaceBase()
  {
    for(typename std::vector<ScalarField<Type>*>::iterator it = scalarFields.begin();
	it != scalarFields.end(); ++it){
      delete (*it);
    }
  };
  //work that must be once before timestepping starts
  virtual void PreIterate() = 0;
  //work that must be done once, before each timestep occurs
  virtual void PreTimeAdvance() = 0;
  //work that must happen to advance inside the newton loop
  virtual void NewtonIterate() = 0;
  //work to advance in time, copy down solution levels, post process, etc.
  virtual void PostTimeAdvance() = 0;

  //output functions
  virtual void WriteAvailableFields()const = 0;
  virtual void WriteSolution() = 0;
  virtual void WriteRestartFile() = 0;
  virtual void ReadRestartFile() = 0;
  void AddScalarField(std::string name)
  {
    for(typename std::vector<ScalarField<Type>*>::iterator it = scalarFields.begin(); 
	it != scalarFields.end(); ++it){
      ScalarField<Type>& field = **it;
      if(field.IsNamed(name)){
	std::cerr << "Field of name " + name + " already exists... must have unique name";
	return;
      }
    }
    ScalarField<Type>* temp =  new ScalarField<Type>(name);
    scalarFields.push_back(temp);
  };
  ScalarField<Type>& GetScalarField(std::string name)
  {
    for(typename std::vector<ScalarField<Type>*>::iterator it = scalarFields.begin();
	it != scalarFields.end(); ++it){
      ScalarField<Type>& field = **it;
      if(field.IsNamed(name)){
	return field;
      }
    }
    std::cerr << "WARNING: In SolutionSpaceBase::GetScalarField() field named " << name << " not found";
    std::cerr << "\tReturning bogus reference" << std::endl;
    return **scalarFields.begin();
  };

  //print timers
  virtual void PrintTimers() = 0;

  std::vector<ScalarField<Type>* > scalarFields;
  
  TimerList timers;

private:
  
};

#endif
