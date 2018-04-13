#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL py_ARRAY_API

template <class Type>
void PythonWrapper::SetInitialConditions(Type* Qinf, int neqn, int nauxvars, Type* Qreturn, Type* coordsXYZ)
{
  //I'm not sure why but if I remove this import we segfault, it appears numpy is being
  //unloaded somehow from initialization
  import_array();

  //This routine assumes that Qinf and Qreturn are the same size... be warned
  //This routine is not yet generalized to deal with complex variables
  size_t neqns = neqn;
  size_t nauxvarss = nauxvars;
  size_t tot = neqn+nauxvars;
  npy_intp inputsize = 1;
  PyObject *npqinf,*npeqn,*npauxvars,*npq,*npxyz; 
  npqinf = (PyObject*)PyArray_SimpleNewFromData(1, &inputsize, NPY_DOUBLE, (double*)Qinf);
  npeqn = PyInt_FromSize_t(neqns);
  npauxvars = PyInt_FromSize_t(nauxvarss);
  npq = (PyObject*)PyArray_SimpleNewFromData(1, &inputsize, NPY_DOUBLE, (double*)Qreturn);
  npy_intp dim = 3;
  npxyz = (PyObject*)PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, coordsXYZ);
  PyObject* pArgs = PyTuple_New(5);
  PyTuple_SetItem(pArgs, 0, npqinf);
  PyTuple_SetItem(pArgs, 1, npeqn);
  PyTuple_SetItem(pArgs, 2, npauxvars);
  PyTuple_SetItem(pArgs, 3, npq);
  PyTuple_SetItem(pArgs, 4, npxyz);
   
  PyObject* pValue = PyObject_CallObject(pFunc, pArgs);
  //TODO: check that q was filled out
}

template <class Type>
void PythonWrapper::GetBoundaryVariables(Type* QL, Type* QR, int neqn, int nauxvars, Type* wallXYZ)
{
  //I'm not sure why but if I remove this import we segfault, it appears numpy is being
  //unloaded somehow from initialization
  import_array();

  //This routine assumes that Qinf and Qreturn are the same size... be warned
  //This routine is not yet generalized to deal with complex variables
  size_t neqns = neqn;
  size_t nauxvarss = nauxvars;
  size_t tot = neqn+nauxvars;
  npy_intp inputsize = 1;
  PyObject *npql,*npqr,*npeqn,*npauxvars,*npq,*npxyz; 
  npql = (PyObject*)PyArray_SimpleNewFromData(1, &inputsize, NPY_DOUBLE, (double*)QL);
  npqr = (PyObject*)PyArray_SimpleNewFromData(1, &inputsize, NPY_DOUBLE, (double*)QR);
  npeqn = PyInt_FromSize_t(neqns);
  npauxvars = PyInt_FromSize_t(nauxvarss);
  npy_intp dim = 3;
  npxyz = (PyObject*)PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, wallXYZ);
  PyObject* pArgs = PyTuple_New(5);
  PyTuple_SetItem(pArgs, 0, npql);
  PyTuple_SetItem(pArgs, 1, npqr);
  PyTuple_SetItem(pArgs, 2, npeqn);
  PyTuple_SetItem(pArgs, 3, npauxvars);
  PyTuple_SetItem(pArgs, 4, npxyz);
   
  PyObject* pValue = PyObject_CallObject(pFunc, pArgs);
}
