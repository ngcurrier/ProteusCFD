template <class Type>
BCObj<Type>::BCObj()
{
  neqn = 0;
  nvars = 0;
  Qref = NULL;
  QrefFromFile = 0;

  periodic = -1;
  Qmin = NULL;
  Qmax = NULL;
  period = -1.0; 
  flowDirection[0] = flowDirection[1] = flowDirection[2] = -1.0;
  rotationAxis[0] = rotationAxis[1] = rotationAxis[2] = -1.0;
  omega = 0.0;
  slipDirection[0] = slipDirection[2] = slipDirection[2] = -1.0;
  slipSpeed = 0.0;
  velocity = -1.0;
  backPressure = -1.0;
  density = -1.0;
  movement = false;
  massFractions = NULL;
  movingBC = -1;
  bleedSteps = 0;
  //set to constant temperature ratio of 1.0; twall < 0.0 is adiabatic
  twall = 1.0;
  flux = 0.0;

  //default to have the parallel factag of zero
  //this is done since we make a bcobj that never gets initialized
  //for this one
  factag = 0;
}

template <class Type>
BCObj<Type>::~BCObj()
{
  if(QrefFromFile){
    delete [] Qref;
  }
  delete [] Qmin;
  delete [] Qmax;
}

template <class Type>
void BCObj<Type>::GetQref(Real* QrefRet)
{
  Int i;
  for(i = 0; i < this->nvars; i++){
    QrefRet[i] = Qref[i];
  }
  return;

}

template <class Type>
void BCObj<Type>::GetQref(RCmplx* QrefRet)
{
  Int i;
  for(i = 0; i < this->nvars; i++){
    QrefRet[i] = Qref[i];
  }
  return;

}

template <class Type>
Int BCObj<Type>::IsMovingBC()
{
  return ((movingBC >= 0));
}
