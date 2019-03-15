#ifndef TEMPORAL_CONTROL_H__
#define TEMPORAL_CONTROL_H__

#include "general.h"
#include "parameterParser.h"
#include "exceptions.h"
#include <string>
#include <fstream>
#include <vector>

template <class Type>
class TemporalControl
{
public:
  TemporalControl();
  ~TemporalControl();
  Int Read(std::string casename, std::string pathname);
  void Print();

  Int nSteps;
  Int newtonIter;

  Type newtonConvergence; // convergence criterion for inner/newton iterates
  Type dt;

private:
  void ReadTemporalControlSegment(std::ifstream& fin, std::streampos locBegin, std::streampos locEnd);
  Int ParseLine(std::string& line);

  //a list of all of the possible parameters which are checked for
  std::vector<Parameter<Bool> > paramListBool;
  std::vector<Parameter<Int> > paramListInt;
  std::vector<Parameter<Type> > paramListReal;
};

//include implementations
#include "temporalControl.tcc"

#endif
