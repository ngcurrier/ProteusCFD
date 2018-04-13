#ifndef PYTHON_INTERFACE_H__
#define PYTHON_INTERFACE_H__

#ifdef _HAS_PYTHON

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <vector>
#include <string>
#include <iostream>

class PythonWrapper
{
  
 public:
  static int refCounter;
  static int numpyLoaded;
  
  PythonWrapper(std::string path, std::string fileRoot, std::string functionName);
  ~PythonWrapper();

  //Utilities
  void AddToPath(std::string path);
  
  //These are demo stubouts to build from
  std::vector<double> CallDoubleVectorFunction(std::vector<double>& input);
  int CallTwoIntFunction(int a, int b);
  void CallBlank();

  //These are the real solver interfaces
  template <class Type>
  void SetInitialConditions(Type* Qinf, int neqn, int nauxvars, Type* Qreturn, Type* coordsXYZ);
  template <class Type>
  void GetBoundaryVariables(Type* QL, Type* QR, int neqn, int nauxvars, Type* wallXYZ);
  
 protected:

 private:
  PythonWrapper();  //hide default implementation
  PyObject* pName;   //This is the python file name "i.e. the *.py"
  PyObject* pModule; //This is the python module loaded of that file
  PyObject* pFunc;   //This is the active python function handle for calling out
};

#include "pythonInterface.tcc"

#endif // end HAS_PYTHON

#endif // end HEADER DEF
