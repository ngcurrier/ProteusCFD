#ifndef TEMPORAL_CONTROL_H__
#define TEMPORAL_CONTROL_H__

#include "general.h"
#include "parameterParser.h"
#include <string>
#include <fstream>
#include <vector>

template <class Type>
class TemporalControl
{
public:
  TemporalControl();
  ~TemporalControl();
  void ReadTemporalControl(std::ifstream& fin, std::streampos locBegin, std::streampos locEnd);
  Int ParseLine(std::string& line);
  void Print();

  Int nSteps;
  Int newtonIter;
  Type dt;

private:
  //a list of all of the possible parameters which are checked for
  std::vector<Parameter<Bool> > paramListBool;
  std::vector<Parameter<Int> > paramListInt;
  std::vector<Parameter<Type> > paramListReal;
};

//free function which will parse temporal control things from the param file
template <class Type>
Int ReadTemporalControl(TemporalControl<Type>& temporalControl, std::string casename, std::string pathname);

//include implementations
#include "temporalControl.tcc"

#endif
