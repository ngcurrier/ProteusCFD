template <class Type>
void PythonWrapper::SetInitialConditions(const Type* Qinf, int neqn, int nauxvars, Type* Qreturn, const Type* coordsXYZ)
{
  //This routine assumes that Qinf and Qreturn are the same size... be warned
  //This routine is not yet generalized to deal with complex variables
  npy_intp inputsize = neqn+nauxvars;
  PyObject* npqinf = (PyObject*)PyArray_SimpleNewFromData(1, &inputsize, PyArray_DOUBLE, (double*)Qinf);
  PyObject* npeqn = PyInt_FromLong(neqn);
  PyObject* npauxvars = PyInt_FromLong(nauxvars);
  PyObject* npq = (PyObject*)PyArray_SimpleNewFromData(1, &inputsize, PyArray_DOUBLE, (double*)Qreturn);
  npy_intp dim = 3;
  PyObject* npxyz = (PyObject*)PyArray_SimpleNewFromData(1, &dim, PyArray_DOUBLE, (double*)coordsXYZ);
  PyObject* pArgs = PyTuple_New(5);
  PyTuple_SetItem(pArgs, 0, npqinf);
  PyTuple_SetItem(pArgs, 1, npeqn);
  PyTuple_SetItem(pArgs, 2, npauxvars);
  PyTuple_SetItem(pArgs, 3, npq);
  PyTuple_SetItem(pArgs, 4, npxyz);
   
  PyObject* pValue = PyObject_CallObject(pFunc, pArgs);

  //TODO: check that q was filled out
}
