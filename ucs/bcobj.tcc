template <class Type>
BCObj<Type>::BCObj() :
  name(BCs[Proteus_NULL]), bcType(Proteus_NULL)
{
  //NOTE: all values stored in bcobj are expected to be dimensional in SI
  // we non-dimensionalize these values elsewhere whenever they are used
  // b/c the non-dimensionalization scheme changes depending on what eqnset
  // we might be using
  // All Q values will be set and stored non-dimensionally by the eqnset class if used
  neqn = 0;
  nvars = 0;
  Qref = NULL;
  QrefFromFile = 0;

  Qmin = NULL;
  Qmax = NULL;

  flowDirection[0] = flowDirection[1] = flowDirection[2] = -1.0;
  rotationAxis[0] = rotationAxis[1] = rotationAxis[2] = -1.0;
  omega = 0.0;         // rad/s
  slipDirection[0] = slipDirection[2] = slipDirection[2] = -1.0;
  slipSpeed = 0.0;     // m/s
  velocity = -1.0;     // m/s
  backPressure = -1.0; // pa
  density = -1.0;      // kg/m^3
  movement = false;
  massFractions = NULL;
  movingBC = -1;
  bleedSteps = 0;
  //set to constant temperature ratio of 1.0; twall < 0.0 is adiabatic
  twall = 1.0;         // K
  flux = 0.0;          // W/m^2

  //default to have the parallel factag of negative number --> NULL BC
  factag = -1;
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
}

template <class Type>
void BCObj<Type>::GetQref(RCmplx* QrefRet)
{
  Int i;
  for(i = 0; i < this->nvars; i++){
    QrefRet[i] = Qref[i];
  }
}

template <class Type>
bool BCObj<Type>::IsMovingBC()
{
  return (movingBC >= 0);
}

template <class Type>
void BCObj<Type>::SetBCType(Int setType)
{
  if(setType < 0 || setType > NUM_BC_TYPES){
    std::stringstream ss;
    ss << "BCObj::SetBCType() requested to set invalid bc type: ";
    ss << "type requested is - " << setType <<  " and max available is " << NUM_BC_TYPES << std::endl;
    Abort << ss.str();
  }
  else{
    bcType = setType;
    std::cout << "Setting BC [" << factag << "] to type: " << BCs[setType] << std::endl;
  }
}
