//Equations of state - for differing types of flow
template <class Type>
EOS<Type>::EOS()
{
  return;
}

template <class Type>
EOS<Type>::~EOS()
{
  return;
}

template <class Type>
IdealGasEOS<Type>::IdealGasEOS()
{
  this->eosType = 0;
  return;
}

template <class Type>
Type IdealGasEOS<Type>::GetT(Type Rmix, Type rho, Type P)
{
  Type T = P/(rho*Rmix);
  return T;
}

template <class Type>
Type IdealGasEOS<Type>::GetRho(Type Rmix, Type P, Type T)
{
  Type rho = P/(Rmix*T);
  return rho;
}

template <class Type>
Type IdealGasEOS<Type>::GetP(Type Rmix, Type rho, Type T)
{
  return (rho*Rmix*T);
}

template <class Type>
Type IdealGasEOS<Type>::GetdT_dP(Type Rmix, Type rho, Type P, Type T)
{
  return (1.0/(rho*Rmix));
}

template <class Type>
Type IdealGasEOS<Type>::GetdT_dRho(Type Rmix, Type rho, Type P, Type T)
{
  return (-P/(Rmix*rho*rho));
}

template <class Type>
Type IdealGasEOS<Type>::GetdP_dR(Type Rmix, Type rho, Type P, Type T)
{
  return (rho*T);
}

template <class Type>
Type IdealGasEOS<Type>::GetdP_dRho(Type Rmix, Type rho, Type P, Type T)
{
  return (Rmix*T);
}

template <class Type>
Type IdealGasEOS<Type>::GetdP_dT(Type Rmix, Type rho, Type P, Type T)
{
  return (rho*Rmix);
}

template <class Type>
Type IdealGasEOS<Type>::GetCp(Type Cv, Type Rmix, Type /*rho*/, Type /*P*/, Type /*T*/)
{
  //Mayer's relation
  return (Cv + Rmix);
}

template <class Type>
Type IdealGasEOS<Type>::GetCv(Type Cp, Type Rmix, Type /*rho*/, Type /*P*/, Type /*T*/)
{
  //Mayer's relation
  Type Cv = Cp - Rmix;
  if(real(CAbs(Cv)) < 0.0) {
    std::stringstream ss;
    ss << "IdealGasEOS::GetCv() computed a negative cv\n";
    ss << "Rmixture: " << Rmix << "\n";
    ss << "Cp: " << Cp << "\n";
    Abort << ss.str();
  }
  return (Cp - Rmix);
}

template <class Type>
Type IdealGasEOS<Type>::GetdT_dR(Type Rmix, Type rho, Type P, Type T)
{
  return (-P/(rho*Rmix*Rmix));
}

template <class Type>
Type IdealGasEOS<Type>::GetRhoR(Type P, Type T)
{
  return (P/T);
}
