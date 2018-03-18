template <class Type>
void PythonWrapper::SetInitialConditions(Type* Qinf, int neqn, int nauxvars, Type* Qreturn, Type* coordsXYZ)
{
  //This routine assumes that Qinf and Qreturn are the same size... be warned
  //This routine is not yet generalized to deal with complex variables
  size_t tot = neqn+nauxvars;
  std::cout << "Size input is: " << tot << std::endl;
  npy_intp inputsize = tot;
  PyObject *npqinf,*npeqn,*npauxvars,*npq,*npxyz; 
  std::cout << Qinf << " " << neqn << " " << nauxvars << " "  << Qreturn<< " "  << coordsXYZ << " " << &inputsize << std::endl;
  npqinf = (PyObject*)PyArray_SimpleNewFromData(1, &inputsize, PyArray_DOUBLE, (double*)Qinf);
  npeqn = PyInt_FromLong(neqn);
  npauxvars = PyInt_FromLong(nauxvars);
  npq = (PyObject*)PyArray_SimpleNewFromData(1, &inputsize, PyArray_DOUBLE, (double*)Qreturn);
  npy_intp dim = 3;
  npxyz = (PyObject*)PyArray_SimpleNewFromData(1, &dim, PyArray_DOUBLE, coordsXYZ);
  PyObject* pArgs = PyTuple_New(5);
  PyTuple_SetItem(pArgs, 0, npqinf);
  PyTuple_SetItem(pArgs, 1, npeqn);
  PyTuple_SetItem(pArgs, 2, npauxvars);
  PyTuple_SetItem(pArgs, 3, npq);
  PyTuple_SetItem(pArgs, 4, npxyz);
   
  PyObject* pValue = PyObject_CallObject(pFunc, pArgs);
  //TODO: check that q was filled out
}
