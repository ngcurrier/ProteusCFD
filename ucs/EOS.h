#ifndef EOS_H__
#define EOS_H__

//Equations of state - for differing types of flow
#include "general.h"

template <class Type>
class EOS
{
public:
  EOS();
  virtual ~EOS();


  //0 - perfect gas
  //1 - ...
  //2 - ... 
  Int eosType;

  
  virtual Type GetT(Type Rmix, Type rho, Type P) = 0;
  virtual Type GetRho(Type Rmix, Type P, Type T) = 0;
  virtual Type GetP(Type Rmix, Type rho, Type T) = 0;
  //FLUENT uses: drho_dT, drho_dp, dh_dT, dh_dP
  virtual Type GetdT_dP(Type Rmix, Type rho, Type P, Type T) = 0;
  virtual Type GetdT_dRho(Type Rmix, Type rho, Type P, Type T) = 0;
  virtual Type GetdT_dR(Type Rmix, Type rho, Type P, Type T) = 0;
  virtual Type GetdP_dRho(Type Rmix, Type rho, Type P, Type T) = 0;
  virtual Type GetdP_dR(Type Rmix, Type rho, Type P, Type T) = 0;
  virtual Type GetdP_dT(Type Rmix, Type rho, Type P, Type T) = 0;
  //Cp - Cv = T * (dP/dT)|_T * (dV/dT)|_T -- must be derived for every new EOS
  virtual Type GetCp(Type Cv, Type Rmix, Type rho, Type P, Type T) = 0;
  virtual Type GetCv(Type Cp, Type Rmix, Type rho, Type P, Type T) = 0;
  virtual Type GetRhoR(Type P, Type T) = 0;

private:
};

template <class Type>
class IdealGasEOS : public EOS<Type>
{

public:

  IdealGasEOS();

  Type GetT(Type Rmix, Type rho, Type P);
  Type GetRho(Type Rmix, Type P, Type T);
  Type GetP(Type Rmix, Type rho, Type T);
  Type GetdT_dP(Type Rmix, Type rho, Type P, Type T);
  Type GetdT_dRho(Type Rmix, Type rho, Type P, Type T);
  Type GetdT_dR(Type Rmix, Type rho, Type P, Type T);
  Type GetdP_dRho(Type Rmix, Type rho, Type P, Type T);
  Type GetdP_dR(Type Rmix, Type rho, Type P, Type T);
  Type GetdP_dT(Type Rmix, Type rho, Type P, Type T);
  Type GetCp(Type Cv, Type Rmix, Type rho, Type P, Type T);
  Type GetCv(Type Cp, Type Rmix, Type rho, Type P, Type T);
  Type GetRhoR(Type P, Type T);


private:

};

template <class Type>
void CreateEOS(EOS<Type>** eos, Int eosType)
{
  if(eosType == 0){
    *eos = new IdealGasEOS<Type>();
  }
  else{
    std::cerr << "WARNING: EOS type " << eosType << " not defined!" << std::endl;
  }
  return;
}

//include implementations
#include "EOS.tcc"


#endif
