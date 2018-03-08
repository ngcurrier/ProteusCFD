#ifndef PYTHON_INTERFACE_H__
#define PYTHON_INTERFACE_H__

#ifdef _HAS_PYTHON
#include <Python.h>
#include <numpy/arrayobject.h>
#include <vector>
#include <string>
#include <iostream>

class PythonWrapper
{
 public:
  PythonWrapper(std::string fileRoot, std::string functionName);
  ~PythonWrapper();

  //These are demo stubouts to build from
  std::vector<double> CallDoubleVectorFunction(std::vector<double>& input);
  double CallTwoIntFunction(int a, int b);
  
 protected:

 private:
  PythonWrapper();  //hide default implementation
  PyObject* pName;   //This is the python file name "i.e. the *.py"
  PyObject* pModule; //This is the python module loaded of that file
  PyObject* pFunc;   //This is the active python function handle for calling out
};



#endif // end HAS_PYTHON

#endif // end HEADER DEF
