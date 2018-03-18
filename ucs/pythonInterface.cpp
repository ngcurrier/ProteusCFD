#include "pythonInterface.h"
#include "exceptions.h"
#include <sstream>

#ifdef _HAS_PYTHON

int PythonWrapper::refCounter = 0;
int PythonWrapper::numpyLoaded = 0;

PythonWrapper::PythonWrapper(std::string path, std::string fileRoot, std::string functionName)
{
  Py_Initialize();
  refCounter++;
  //This is a hack b/c numpy never really cleans up if you try to unload it... so we just let it ride
  //It probably leaks memory but upstream won't fix...
  if(!numpyLoaded){
    std::cout << "PythonWrapper:: loading NumPy" << std::endl;
    numpyLoaded++;
    import_array();
  }

  AddToPath(path);
  
  pName = PyString_FromString(fileRoot.c_str());
  pModule = PyImport_Import(pName);

  std::cout << "PythonWrapper:: imported" << std::endl;
  
  if(pModule != NULL){
    pFunc = PyObject_GetAttrString(pModule, functionName.c_str());

    if(pFunc && PyCallable_Check(pFunc)){
      //don't do anything
    }
    else{
      if(PyErr_Occurred()){
	PyErr_Print();
      }
      std::stringstream ss;
      ss << "Cannot load function in PythonWrapper " << functionName << " from file " << fileRoot << std::endl;
      Abort << ss.str();
    }
  }
  else{
    PyErr_Print();
    std::stringstream  ss;
    ss << "PythonWrapper failed to load " << fileRoot << std::endl;
    Abort << ss.str();
  }

  //this is so we can use Numpy arrays
  std::cout << "PythonWrapper:: initialized for " << fileRoot << ":"<< functionName << std::endl;
}

PythonWrapper::~PythonWrapper()
{

  
  //call xdecref in case any of these failed with a NULL return DECREF() requires a valid pointer
  Py_XDECREF(pFunc);
  Py_XDECREF(pModule);
  Py_XDECREF(pName);
  refCounter--;
  //WARNING, If I call py_finalize after numpy is imported and then py_initialize then the tool
  //segfaults, it should not work this way so we just never call py_finalize.
  //This is a well documented issue that just hasn't been fixed https://github.com/numpy/numpy/issues/8097
  if(refCounter == 0 && false){
    std::cout << "PythonWrapper:: last reference deleted, calling Py_Finalize()" << std::endl;
    Py_Finalize();
  }
}

//The python function is expected to accept a numpy array
//The argument should (*) splat expand the pyobject args
std::vector<double> PythonWrapper::CallDoubleVectorFunction(std::vector<double>& input)
{
  npy_intp inputsize = input.size();
  PyObject* numpyarray = (PyObject*)PyArray_SimpleNewFromData(1, &inputsize, PyArray_DOUBLE, (double*)input.data());
  PyObject* pArgs = PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, numpyarray);
   
  PyObject* pValue = PyObject_CallObject(pFunc, pArgs);

  PyArrayObject *np_ret = reinterpret_cast<PyArrayObject*>(pValue);
  int size = PyArray_DIMS(pValue)[0];
  std::cout << "Received array with dimensions " << size << std::endl;

  // Convert back to C++ array and print.
  double* c_out = reinterpret_cast<double*>(PyArray_DATA(pValue));
  std::vector<double> returnVec(size);
  //print and load up vector
  std::cout << "Printing output array" << std::endl;
  for (int i = 0; i < size; i++){
    std::cout << c_out[i] << ' ';
    returnVec[i] = c_out[i];
  }
  std::cout << std::endl;

  //free references
  Py_DECREF(pValue);
  Py_DECREF(pArgs);
  //Do Not delete numpyarray since it is simply a cast of the memory for the vector
  //Py_DECREF(numpyarray);
    
  return returnVec;
}

int PythonWrapper::CallTwoIntFunction(int a, int b)
{
  PyObject* pArgs;
  PyObject* pValue;
  int returnVal = -11111;

  pArgs = PyTuple_New(2);
  pValue = PyInt_FromLong(a);
  if(!pValue){
    Py_DECREF(pArgs);
    std::cerr << "Cannot convert argument in PythonWrapper" << std::endl;
    //TODO: call throwerror()
    return -999;
  }
  else{
    PyTuple_SetItem(pArgs, 0, pValue);
  }

  pValue = PyInt_FromLong(b);
  if(!pValue){
    Py_DECREF(pArgs);
    std::cerr << "Cannot convert argument in PythonWrapper" << std::endl;
    //TODO: call throwerror()
    return -998;
  }
  else{
    PyTuple_SetItem(pArgs, 1, pValue);
  }

  //Recall that if you need to unpack pArgs on the python side to use the (*) splat operator
  pValue = PyObject_CallObject(pFunc, pArgs);

  if(pValue == NULL){
    PyErr_Print();
    std::cerr << "PythonWrapper call did not return anything" << std::endl;
    return -997;
  }

  //we expect back a float
  returnVal = PyInt_AsLong(pValue);
  //returnVal = PyFloat_AsDouble(pValue);

  Py_DECREF(pValue);
  Py_DECREF(pArgs);
  
  return returnVal;
}

void PythonWrapper::CallBlank()
{
  std::cout << "Calling Blank Python module" << std::endl;
  PyObject* pValue;
  PyObject* pArgs = PyTuple_New(1);
  pValue = PyObject_CallObject(pFunc, pArgs);
}

void PythonWrapper::AddToPath(std::string path)
{
  //import the python path we need, append to current path
  std::string base("import sys\nsys.path.append('");
  std::string end("')\n");
  std::string command = base + path + end;

  PyRun_SimpleString(command.c_str());
  std::cout << "PythonWrapper:: added " << path << " to PYTHONPATH" << std::endl;
}

#endif
