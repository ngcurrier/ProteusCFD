#include "pythonInterface.h"

#ifdef _HAS_PYTHON

PythonWrapper::PythonWrapper(std::string fileRoot, std::string functionName)
{
  Py_Initialize();

  pName = PyString_FromString(fileRoot.c_str());
  pModule = PyImport_Import(pName);

  if(pModule != NULL){
    pFunc = PyObject_GetAttrString(pModule, functionName.c_str());

    if(pFunc && PyCallable_Check(pFunc)){
      //don't do anything
    }
    else{
      if(PyErr_Occurred()){
	PyErr_Print();
      }
      std::cerr << "Cannot load function in PythonWrapper " << functionName << " from file " << fileRoot << std::endl;
      //TODO: throwerror
    }
  }
  else{
    PyErr_Print();
    std::cerr << "PythonWrapper failed to load " << fileRoot << std::endl;
  }
  
  import_array(); //this is so we can use Numpy arrays
}

PythonWrapper::~PythonWrapper()
{
  //call xdecref in case any of these failed with a NULL return
  // DECREF() requires a valid pointer
  Py_XDECREF(pFunc);
  Py_XDECREF(pModule);
  Py_XDECREF(pName);
  Py_Finalize();
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
  Py_DECREF(numpyarray);
    
  return returnVec;
}

double PythonWrapper::CallTwoIntFunction(int a, int b)
{
  PyObject* pArgs;
  PyObject* pValue;
  double returnVal;

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
    return -999;
  }
  else{
    PyTuple_SetItem(pArgs, 1, pValue);
  }

  //Recall that if you need to unpack pArgs on the python side to use the (*) splat operator
  pValue = PyObject_CallObject(pFunc, pArgs);

  if(pValue == NULL){
    PyErr_Print();
    std::cerr << "PythonWrapper call did not return anything" << std::endl;
    return -999;
  }

  //we expect back a float
  //PyInt_AsLong(pValue);
  PyFloat_AsDouble(pValue);

  Py_DECREF(pValue);
  Py_DECREF(pArgs);
  
  return returnVal;
}

#endif
