template <class Type>
ChemModel<Type>::ChemModel(Param<Type>& param)
{
  Int i, j, k;
  Int eosType = 0; //ideal gas
  
  caseString = param.casestring;
  CreateEOS(&eos, eosType);

  nreactions = GetNumberOfReactionsFromFile();
  reactions = new Reaction<Type>[nreactions];

  std::string rxnfile = caseString + ".rxn";
  std::cout << "CHEM_RXN: Reading " << nreactions << " reactions from file " << rxnfile 
	    << std::endl;
  Int SpeciesCount = 0;
  Int CatCount = 0;
  for(i = 1; i <= nreactions; i++){
    reactions[i-1].ReadReactionFromFile(rxnfile,i);
    SpeciesCount += reactions[i-1].nspecies;
    CatCount += reactions[i-1].ncatalysts;
  }

  //from species information now contained in each reaction object
  //build a global list of all unique species
  std::string* speciesList = new std::string[SpeciesCount+CatCount];
  
  Int listSize = 0;
  Bool inList;
  //create a unique list using reactants
  for(i = 0; i < nreactions; i++){
    for(j = 0; j < reactions[i].nspecies; j++){
      inList = false;
      for(k = 0; k < listSize; k++){
	if(reactions[i].speciesSymbols[j] == speciesList[k]){
	  inList = true;
	  break;
	}
      }
      if(!inList){
	speciesList[listSize] = reactions[i].speciesSymbols[j];
	listSize++;
      }
    }
  }
  //add catalysts to the list
  for(i = 0; i < nreactions; i++){
    for(j = 0; j < reactions[i].ncatalysts+reactions[i].ntbodies; j++){
      inList = false;
      for(k = 0; k < listSize; k++){
	if(reactions[i].catSymbols[j] == speciesList[k]){
	  inList = true;
	  break;
	}
      }
      if(!inList){
	speciesList[listSize] = reactions[i].catSymbols[j];
	listSize++;
      }
    }
  }

  //set nspecies
  nspecies = listSize;
  species = new Species<Type>[nspecies];
  //initialize species we need
  for(i = 0; i < nspecies; i++){
    species[i].Init(speciesList[i], param.viscous, param.chemDB);
  }
  for(i = 0; i < nreactions; i++){
    reactions[i].SetSpeciesPointer(species, nspecies);
  }

  std::cout << "SPECIES IN MODEL: " << std::endl;
  std::cout << "================= " << std::endl;
  for(i = 0; i < nspecies; i++){
    std::cout << i << ": " << species[i].symbol << std::endl;
  }
  std::cout << std::endl;

  for(i = 0; i < nspecies; i++){
    species[i].Print();
  }

  for(i = 0; i < nreactions; i++){
    reactions[i].Print();
    std::cout << std::endl;
  }

  delete [] speciesList;

  return;
}

template <class Type>
ChemModel<Type>::ChemModel(Int modelId_, Param<Type>& param)
{
  Int eosType = 0; //ideal gas
  this->modelId = modelId_;
  //setup equation of state object
  CreateEOS(&eos, eosType);
  species = new Species<Type>[2];
  reactions = new Reaction<Type>[2];

  std::cerr << "WARNING: This is NOT finished... stubout only! ChemModel non-default constructor!" << std::endl;

  return;
}

template <class Type>
ChemModel<Type>::~ChemModel()
{
  delete [] reactions;
  delete [] species;
  delete eos;

  return;
}

template <class Type>
Int ChemModel<Type>::GetNumberOfReactionsFromFile()
{
  Int nreactions = 0;
  std::string fileName = this->caseString + ".rxn";
  std::ifstream fin;

  std::stringstream ss;
  std::string line, trash, subline;
  size_t loc;
  char c;


  std::string rxnKey = "rxnInModel";
  
  fin.open(fileName.c_str());
  if(fin.is_open()){
    while(!fin.eof()){
      c = fin.peek();
      //get rid of comments and blank lines
      if(c == '#' || c == ' ' || c == '\n'){
	getline(fin, trash);
	trash.clear();
	continue;
      }
      getline(fin, line);
      loc = line.find(rxnKey);
      if(loc != std::string::npos){
	loc = line.find('=');
	loc += 1;
	subline = line.substr(loc);
	ss << subline;
	ss >> nreactions;
	ss.clear();
	break;
      }
    }
  }
  else{
    std::cerr << "CHEM_RXN READFILE: Cannot open reaction file --> " << fileName 
	      << std::endl;
    return (0);

  }

  return nreactions;
}

template <class Type>
void ChemModel<Type>::GetFluidProperties(Type rho, Type* rhoi, Type* cvi, Type* cv, Type* R)
{
  Int i;
  *R = 0.0;
  *cv = 0.0;

  for(i = 0; i < this->nspecies; i++){
    *R += rhoi[i]*species[i].R;
    *cv += rhoi[i]*cvi[i];
  }

  (*cv) /= rho;
  (*R) /= rho;

  return;
}

template <class Type>
void ChemModel<Type>::GetMassProductionRates(Type* rhoi, Type T, Type* wdot)
{
  Int i;
  Int j;
  for(i = 0; i < this->nspecies; i++){
    wdot[i] = 0.0;
    for(j = 0; j < this->nreactions; j++){
      wdot[i] += this->reactions[j].GetMassProductionRate(rhoi, T, i);
    }
  }

  return;
}

template <class Type>
Type ChemModel<Type>::GetSpecificEnthalpy(Type* massfrac, Type T, Type * hi)
{
  Type h = 0.0;
  for(Int i = 0; i < nspecies; i++){
    hi[i] = species[i].GetH(T);
    h += hi[i]*massfrac[i];
  }
  return h;
}

template <class Type>
Type ChemModel<Type>::dEtdT_dEtdRhoi_FD(Type* rhoi, Type T, Type P, Type v2, Type* dEtdRhoi)
{
  Int i, j;
  Type dEtdT;
  Type Ht, Htpu, Htpd, Et, Etpu, Etpd;
  Type momV2 = 0.5*v2;
  Type h = 1.0e-8;
  Type rho;
  Type rhop;
  Type Rmix = 0.0;
  Type Rmixp;
  Type* rhoip = (Type*)alloca(sizeof(Type)*nspecies);
  Type* X = (Type*)alloca(sizeof(Type)*nspecies);
  Type* Xp = (Type*)alloca(sizeof(Type)*nspecies);
  Type* hi = (Type*)alloca(sizeof(Type)*nspecies);
  Type* hip = (Type*)alloca(sizeof(Type)*nspecies);
  Type Tpu, Tpd;
  Type Ppu, Ppd;
  Type dEidP;

  rho = 0.0;
  for(i = 0; i < nspecies; i++){
    Rmix += rhoi[i]*species[i].R;
    rho += rhoi[i];
  }
  Rmix /= rho;
  for(i = 0; i < nspecies; i++){
    X[i] = rhoi[i]/rho;
  }
  Ht = rho*GetSpecificEnthalpy(X, T, hi) + rho*momV2;
  Et = Ht - P;

  //calculate the derivative of TotalEnergy w.r.t partial density
  //at constant temperature
  for(i = 0; i < nspecies; i++){
    memcpy(rhoip, rhoi, sizeof(Type)*nspecies);
    rhoip[i] += h;
    rhop = 0.0;
    Rmixp = 0.0;
    //compute new mixture R and total density
    for(j = 0; j < nspecies; j++){
      Rmixp += rhoip[j]*species[j].R;
      rhop += rhoip[j];
    }
    Rmixp /= rhop;
    //compute new massfractions
    for(j = 0; j < nspecies; j++){
      Xp[j] = rhoip[j]/rhop;
    }
    //compute new P based on equation of state used
    Ppu = eos->GetP(Rmixp, rhop, T);
    //compute new total energy based on perturbed state
    Htpu = rhop*GetSpecificEnthalpy(Xp, T, hip) + rhop*momV2;
    Etpu = Htpu - Ppu;
    memcpy(rhoip, rhoi, sizeof(Type)*nspecies);
    rhoip[i] -= h;
    rhop = 0.0;
    Rmixp = 0.0;
    //compute new mixture R and total density
    for(j = 0; j < nspecies; j++){
      Rmixp += rhoip[j]*species[j].R;
      rhop += rhoip[j];
    }
    Rmixp /= rhop;
    //compute new massfractions
    for(j = 0; j < nspecies; j++){
      Xp[j] = rhoip[j]/rhop;
    }
    //compute new T based on equation of state used
    Ppd = eos->GetP(Rmixp, rhop, T);
    //compute new total energy based on perturbed state
    Htpd = rhop*GetSpecificEnthalpy(Xp, T, hip) + rhop*momV2;
    Etpd = Htpd - Ppd;
    //compute finite difference
    dEtdRhoi[i] = (Etpu - Etpd)/(2.0*h);
  }

  Tpu = T+h;
  Tpd = T-h;
  Htpu = rho*GetSpecificEnthalpy(X, Tpu, hip) + rho*momV2;
  Htpd = rho*GetSpecificEnthalpy(X, Tpd, hip) + rho*momV2;
  Ppu = eos->GetP(Rmix, rho, Tpu);
  Ppd = eos->GetP(Rmix, rho, Tpd);
  Etpu = Htpu - Ppu;
  Etpd = Htpd - Ppd;
  dEtdT = (Etpu - Etpd)/(2.0*h);

  return dEtdT;
}

template <class Type>
void ChemModel<Type>::dTdRhoi(Type* rhoi, Type P, Type* dTdRhoi)
{
  Type rho = 0.0;
  Type Rmix = 0.0;
  for(Int i = 0; i < nspecies; i++){
    rho += rhoi[i];
    Rmix += rhoi[i]*species[i].R;
  }
  Rmix /= rho;
  Type T = eos->GetT(Rmix, rho, P);
  
  Type dTdRho = eos->GetdT_dRho(Rmix, rho, P, T);
  Type dTdR = eos->GetdT_dR(Rmix, rho, P, T);
  
  for(Int i = 0; i < nspecies; i++){
    dTdRhoi[i] = dTdR*dRmixdRhoi(rhoi, rho, i) + dTdRho;
  }
}
   

template <class Type>
Type ChemModel<Type>::dEtdP_dEtdRhoi_FD(Type* rhoi, Type T, Type P, Type v2, Type* dEtdRhoi)
{
  Int i, j;
  Type dEtdP;
  Type Ht, Htpu, Htpd, Et, Etpu, Etpd;
  Type momV2 = 0.5*v2;
  Type h = 1.0e-8;
  Type rho;
  Type rhop;
  Type Rmix = 0.0;
  Type Rmixp;
  Type* rhoip = (Type*)alloca(sizeof(Type)*nspecies);
  Type* X = (Type*)alloca(sizeof(Type)*nspecies);
  Type* Xp = (Type*)alloca(sizeof(Type)*nspecies);
  Type* hi = (Type*)alloca(sizeof(Type)*nspecies);
  Type* hip = (Type*)alloca(sizeof(Type)*nspecies);
  Type Ppu, Ppd;
  Type Tpu, Tpd;
  Type dEidP;

  rho = 0.0;
  for(i = 0; i < nspecies; i++){
    Rmix += rhoi[i]*species[i].R;
    rho += rhoi[i];
  }
  Rmix /= rho;
  for(i = 0; i < nspecies; i++){
    X[i] = rhoi[i]/rho;
  }
  Ht = rho*GetSpecificEnthalpy(X, T, hi) + rho*momV2;
  Et = Ht - P;

  //calculate the derivative of TotalEnergy w.r.t partial density
  //at constant temperature
  for(i = 0; i < nspecies; i++){
    memcpy(rhoip, rhoi, sizeof(Type)*nspecies);
    rhoip[i] += h;
    rhop = 0.0;
    Rmixp = 0.0;
    //compute new mixture R and total density
    for(j = 0; j < nspecies; j++){
      Rmixp += rhoip[j]*species[j].R;
      rhop += rhoip[j];
    }
    Rmixp /= rhop;
    //compute new massfractions
    for(j = 0; j < nspecies; j++){
      Xp[j] = rhoip[j]/rhop;
    }
    //compute new P based on equation of state used
    Ppu = eos->GetP(Rmixp, rhop, T);
    //compute new total energy based on perturbed state
    Htpu = rhop*GetSpecificEnthalpy(Xp, T, hip) + rhop*momV2;
    Etpu = Htpu - Ppu;
    memcpy(rhoip, rhoi, sizeof(Type)*nspecies);
    rhoip[i] -= h;
    rhop = 0.0;
    Rmixp = 0.0;
    //compute new mixture R and total density
    for(j = 0; j < nspecies; j++){
      Rmixp += rhoip[j]*species[j].R;
      rhop += rhoip[j];
    }
    Rmixp /= rhop;
    //compute new massfractions
    for(j = 0; j < nspecies; j++){
      Xp[j] = rhoip[j]/rhop;
    }
    //compute new T based on equation of state used
    Ppd = eos->GetP(Rmixp, rhop, T);
    //compute new total energy based on perturbed state
    Htpd = rhop*GetSpecificEnthalpy(Xp, T, hip) + rhop*momV2;
    Etpd = Htpd - Ppd;
    //compute finite difference
    dEtdRhoi[i] = (Etpu - Etpd)/(2.0*h);
  }

  Ppu = P+h;
  Ppd = P-h;
  Tpu = eos->GetT(Rmix, rho, Ppu);
  Tpd = eos->GetT(Rmix, rho, Ppd);
  Htpu = rho*GetSpecificEnthalpy(X, Tpu, hip) + rho*momV2;
  Htpd = rho*GetSpecificEnthalpy(X, Tpd, hip) + rho*momV2;
  Etpu = Htpu - Ppu;
  Etpd = Htpd - Ppd;
  dEtdP = (Etpu - Etpd)/(2.0*h);

  return dEtdP;
}

template <class Type>
Type ChemModel<Type>::dEtdT_dEtdRhoi(Type* rhoi, Type T, Type P, Type v2, Type* dEtdRhoi)
{
  Type hv2 = 0.5*v2;
  Type rhomix = 0.0;
  Type Rmix = 0.0;
  for(Int i = 0; i < nspecies; i++){
    rhomix += rhoi[i];
    Rmix += rhoi[i]*species[i].R;
  }
  Rmix /= rhomix;
  Type dPdrhomix = eos->GetdP_dRho(Rmix, rhomix, P, T);
  Type dPdRmix = eos->GetdP_dR(Rmix, rhomix, P, T);

  //chain rule F.T.W. at constant T -- pick (rho, T)
  for(Int i = 0; i < nspecies; i++){
    dEtdRhoi[i] = hv2 + species[i].GetH(T) - (dPdrhomix + dPdRmix*dRmixdRhoi(rhoi, rhomix, i));
  }

  //now we are perturbing T -- pick (P, rho)
  Type dEtdT = 0.0;
  Type dPdT = eos->GetdP_dT(Rmix, rhomix, P, T);
  for(Int i = 0; i < nspecies; i++){
    dEtdT += rhoi[i]*(species[i].GetdHdT(T));
  }
  dEtdT -= dPdT;

  return dEtdT;
}

template <class Type>
Type ChemModel<Type>::dEtdP_dEtdRhoi(Type* rhoi, Type T, Type P, Type v2, Type* dEtdRhoi)
{
  Type hv2 = 0.5*v2;
  Type rhomix = 0.0;
  Type Rmix = 0.0;
  for(Int i = 0; i < nspecies; i++){
    rhomix += rhoi[i];
    Rmix += rhoi[i]*species[i].R;
  }
  Rmix /= rhomix;
  Type dPdrhomix = eos->GetdP_dRho(Rmix, rhomix, P, T);
  Type dPdRmix = eos->GetdP_dR(Rmix, rhomix, P, T);

  //chain rule F.T.W. at constant T -- pick (rho, T)
  for(Int i = 0; i < nspecies; i++){
    dEtdRhoi[i] = hv2 + species[i].GetH(T) - (dPdrhomix + dPdRmix*dRmixdRhoi(rhoi, rhomix, i));
  }
  
  //DEFINITION: Et = sum_i (rhoi * hi(T)) - P
  
  //now we are perturbing P -- pick (rho, T)
  Type dEtdP = 0.0;
  Type dTdP = eos->GetdT_dP(Rmix, rhomix, P, T);
  for(Int i = 0; i < nspecies; i++){
    dEtdP += rhoi[i]*(species[i].GetdHdT(T))*dTdP;
  }
  dEtdP -= 1.0;

  return dEtdP;
}

template <class Type>
Type ChemModel<Type>::dRmixdRhoi(Type* rhoi, Type rho, Int i)
{
  Type dRmixdRhoi = species[i].R*(rho - rhoi[i])/(rho*rho);
  for(Int j = 0; j < nspecies; j++){
    if(j == i) continue;
    dRmixdRhoi -= species[j].R*rhoi[j]/(rho*rho);
  }
  return dRmixdRhoi;
}

template <class Type>
Type ChemModel<Type>::WilkesMixtureRule(Type* rhoi, Type* property, Type T)
{
  Type mixtureProp = 0.0;
  Type sqrt8 = sqrt(8.0);

  Type rho = 0.0;
  for(Int i = 0; i < nspecies; i++){
    rho += rhoi[i];
  }

  //compute mol fraction following Cox's dissertation 
  Type* molfrac = (Type*)alloca(sizeof(Type)*nspecies);
	Type summ = 0.0;
  for(Int i = 0; i < nspecies; i++){
		molfrac[i] = (rhoi[i]/rho)/species[i].MW;
		summ += molfrac[i];
  }
  for(Int i = 0; i < nspecies; i++){
		molfrac[i] /= summ;
	}

  Type* visc = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; i++){
    visc[i] = species[i].GetViscosity(T);
  }

  for(Int i = 0; i < nspecies; i++){
    Type wi = 0.0;
    Type fraci = molfrac[i];
    Type MWi = species[i].MW;
    for(Int j = 0; j < nspecies; j++){
      Type fracj = molfrac[j];
      Type MWj = species[j].MW;
      Type temp = (1.0 + sqrt(visc[i]/visc[j])*pow(MWj/MWi, 0.25));
      Type phi = pow((1.0 + MWi/MWj), -0.5)*temp*temp/sqrt8;
      wi += fracj*phi;
    }
    mixtureProp += (fraci/wi)*property[i];
  }
  return mixtureProp;
}

template <class Type>
Type ChemModel<Type>::GetViscosity(Type* rhoi, Type T)
{
  Type* visc = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; i++){
    visc[i] = species[i].GetViscosity(T);
  }
  Type mu = WilkesMixtureRule(rhoi, visc, T);
  return mu;
}

template <class Type>
Type ChemModel<Type>::GetThermalConductivity(Type* rhoi, Type T)
{
  Type* k = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; i++){
    k[i] = species[i].GetThermalConductivity(T);
  }
  Type kt = WilkesMixtureRule(rhoi, k, T);
  return kt;
}
