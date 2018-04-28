#include "strings_util.h"

template <class Type>
ChemModel<Type>::ChemModel()
{
}

template <class Type>
ChemModel<Type>::ChemModel(std::string casestring, std::string databaseFile):
  caseString(casestring), nespecies(0), nspecies(0), reactions(0), databaseFile(databaseFile)
{
  Int i, j, k;
  Int eosType = 0; //ideal gas

  //check if we can open .rxn file first
  std::string rxnfile = caseString + ".rxn";
  std::ifstream rtest;
  rtest.open(rxnfile.c_str());
  if(rtest.is_open()){
    rtest.close();
    if(ReadReactionsFile(rxnfile)){
      std::stringstream ss;
      ss << "Error reading chemical reaction file " << rxnfile;
      Abort << ss.str();
    }
  }
  else{
    //check if we can open .ck file
    std::ifstream ctest;
    std::string filename = caseString + ".ck";
    ctest.open(filename.c_str());
    if(ctest.is_open()){
      ctest.close();
      if(ReadChemkinReactionsFile(filename)){
	std::stringstream ss;
	ss << "Error reading Chemkin chemical reaction file " << filename;
	Abort << ss.str();
      }
    }
  }
  
  //count total elements in model
  std::vector<std::string> nespeciesList = GetElementsInModel();
  nespecies = nespeciesList.size();
  
  std::cout << "SPECIES IN MODEL: " << std::endl;
  std::cout << "================= " << std::endl;
  for(i = 0; i < nspecies; i++){
    std::cout << i << ": " << species[i].symbol << std::endl;
  }
  std::cout << std::endl;

  for(i = 0; i < nspecies; i++){
    species[i].Print();
  }

  // setup one EOS model for each specie
  eos.resize(nspecies);
  std::cout << "Setting up EOSs:\n";
  std::cout << "-------------------------------------\n";
  for(i = 0; i < nspecies; ++i){
    EOS<Type>* eosTmp;
    CreateEOS(&eosTmp, eosType);
    eos[i] = eosTmp;
    std::cout << i << ": " << eos[i]->myType() << std::endl;
  }
  std::cout << "\n";
   
  for(i = 0; i < nreactions; i++){
    reactions[i].Print();
    std::cout << std::endl;
  }

}

template <class Type>
std::vector<std::string> ChemModel<Type>::GetElementsInModel()
{
  std::vector<std::string> vnespecies;
  //TODO: Fill in elemental species counter
  return vnespecies;
}

template <class Type>
Int ChemModel<Type>::ReadReactionsFile(std::string rxnfile)
{
  Int i, j, k;
  nreactions = ReadNumberOfReactionsFromFile(rxnfile);
  reactions = new Reaction<Type>[nreactions];

  std::cout << "CHEM_RXN: Reading " << nreactions << " reactions from file " << rxnfile 
	    << std::endl;

  for(Int i = 1; i <= nreactions; i++){
    reactions[i-1].ReadReactionFromFile(rxnfile,i);
  }
  
  BuildGlobalSpeciesList();
  
  return 0;

}

template <class Type>
void ChemModel<Type>::BuildGlobalSpeciesList()
{
  std::cout << "CHEM_RXN: Building global species list" << std::endl;
  
  Int SpeciesCount = 0;
  for(Int i = 1; i <= nreactions; i++){
    SpeciesCount += reactions[i-1].GetNspecies();
  }
  
  //from species information now contained in each reaction object
  //build a global list of all unique species
  std::vector<std::string> speciesList(SpeciesCount);
  
  Int listSize = 0;
  Bool inList;
  //create a unique list using reactants
  for(Int i = 0; i < nreactions; i++){
    for(Int j = 0; j < reactions[i].GetNspecies(); j++){
      inList = false;
      for(Int k = 0; k < listSize; k++){
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

  //set nspecies
  nspecies = listSize;
  species = new Species<Type>[nspecies];
  //initialize species we need
  for(Int i = 0; i < nspecies; i++){
    species[i].Init(speciesList[i],databaseFile);
  }
  for(Int i = 0; i < nreactions; i++){
    reactions[i].SetSpeciesPointer(species, nspecies);
  }

}

template <class Type>
Int ChemModel<Type>::ReadChemkinReactionsFile(std::string rxnfile)
{
  std::ifstream fin;
  std::string line, trash;
  size_t loc;

  nreactions = ReadNumberOfReactionsFromFileChemkin(rxnfile);
  std::cout << "CHEM_RXN CHEMKIN: Reading " << nreactions << " reactions from file " << rxnfile 
	    << std::endl;

  reactions = new Reaction<Type>[nreactions];
  
  //states for our state machine reader
  enum{
    stateFindingRXNS,
    stateReadingRXNSList,
    stateReadingRXNItem1,
    stateReadingRXNItem2,
    stateReadingRXNItem3,
    stateFoundEND,
    NSTATES
  };

  Int state = stateFindingRXNS;
  Int ireaction = 0;
  
  fin.open(rxnfile.c_str());
  if(fin.is_open()){
    //We simply look through all lines knowing that =, =>, AND <=> represent a reaction, count 'em
    //Ignore anything past a ! in each line
    while(!fin.eof()){
      getline(fin, trash);
      // tokenize based on !, take everything before it
      std::vector<std::string> tokens = Tokenize(trash, '!');
      std::vector<std::string> itokens, ctokens;
      std::vector<double> values;
      std::string LHS, RHS, coeffs, LOW, TROE;
      std::string proteusEquation;
      std::string proteusTBs;
      bool hasTB = false;
      bool haspdep = false;
      bool haslow, hastroe;
      trash.clear();
      std::string line = tokens[0];
      double A, b, Ea;
      std::stringstream ss;
      
      if(line.find("REACTIONS") != std::string::npos) state = stateReadingRXNSList;
      else if((line.find("END") != std::string::npos) && state == stateReadingRXNSList) state = stateFoundEND;
      else if(line.find("=") != std::string::npos) state = stateReadingRXNItem1;
      else if(line.find("=>") != std::string::npos) state = stateReadingRXNItem2;
      else if(line.find("<=>") != std::string::npos) state = stateReadingRXNItem3;
      else{
	//do nothing just skip
      }
      switch(state){
      case stateFindingRXNS:
	// just ignore and keep reading
	break;
      case stateReadingRXNSList:
	break;
      case stateReadingRXNItem1:
      case stateReadingRXNItem2: //we ignore irreversible reactions and treat as reversible
      case stateReadingRXNItem3:
	//All chemkin reactions are modified arrhenius form -- ModArrhenius k = A T^n exp(-Ea/RT)
	reactions[ireaction].rxnType = ModArrheniusCoeff;
	
	// jump down a state to read next reaction
	// Get left and right of equation
	itokens = Split(line, "=");
	StripFromString(itokens[0], "<", itokens[0]);
	StripFromString(itokens[1], ">", itokens[1]);
	LHS = itokens[0];
	RHS = itokens[1];
	//This is questionable but normally works since most users space the Arrhenius coefficients far right
	// and don't normaly space out the equation more than one space
	itokens = Split(RHS, "  ");
	RHS = itokens[0];
	coeffs = itokens[1];
	values = TokenizeToDoubles(coeffs, ' ');
	A = values[0];
	b = values[1];
	Ea = values[2];

	reactions[ireaction].A = A;
	reactions[ireaction].EA = Ea;
	reactions[ireaction].n = b; // difference in notation, temperature exponent
	std::cout << ireaction << ": " << LHS << " <==> " << RHS << std::endl;
	std::cout << "\tA: " << A << " b: " << b << " Ea: " << Ea << std::endl;
	
	RemoveWhitespace(LHS);
	RemoveWhitespace(RHS);
	if(LHS.find("(+M)") != std::string::npos){
	  haspdep = true;
	  hasTB = true;
	  //replace any (M+) values with M
	  Replace(LHS, "(+M)", "+M");
	  Replace(RHS, "(+M)", "+M");
	}
	else if(LHS.find("M") != std::string::npos){
	  hasTB = true;
	}

	//parser currently expects spaces separating operators, make it so!
	itokens = Tokenize(LHS, '+');
	LHS = itokens[0];
	for(Int i = 1; i < itokens.size(); ++i){
	  LHS += " + " + itokens[i];
	}
	itokens = Tokenize(RHS, '+');
	RHS = itokens[0];
	for(Int i = 1; i < itokens.size(); ++i){
	  RHS += " + " + itokens[i];
	}
	
	//create the equation proteus expects to see for parsing
	proteusEquation = LHS + "<==>" + RHS;
	reactions[ireaction].ParseReaction(proteusEquation);
	
	// If the letter M appears by itself in a reaction the next line must be thirdbodies
	// i.e. the rate constant is assumed to be in the low pressure limiting region
	if(hasTB){
	  std::cout << "WARNING: " << LHS << " has third bodies\n";
	  std::cerr << "WARNING: " << LHS << " has third bodies\n";
	  getline(fin, trash);
	  haslow = hastroe = false;
	  if(trash.find("LOW") != std::string::npos){
	    coeffs = GetStringBetween("/", "/", trash, LOW);
	    getline(fin,trash);
	    haslow = true;
	  }
	  if(trash.find("TROE") != std::string::npos){
	    coeffs = GetStringBetween("/", "/", trash, TROE);
	    getline(fin,trash);
	    hastroe = true;
	  }
	  if(trash.find("PLOG") != std::string::npos){
	    //should continue to getline until it fails.. TODO: fix PLOG reading
	  }
	  std::vector<std::string> TBs = Tokenize(trash, '/');
	  reactions[ireaction].thirdBodiesPresent = true;
	  //we need to put the thirdbodies in the following form for parsing
	  //M = NO(+)[1.0], O2(+)[1.0], N2(+)[1.0], O(+)[1.0], N(+)[1.0]
	  proteusTBs = "M =";
	  for(int i = 0; i < TBs.size(); i = i + 2){
	    proteusTBs += " " + TBs[i] + "[" + TBs[i+1] + "],";
	  }
	  proteusTBs = proteusTBs.substr(0, proteusTBs.size()-1);
	  reactions[ireaction].ParseThirdBodies(proteusTBs);
	}
	// If (+M) appears then the rate constant is in the fall off region
	// -- If only LOW/ follows the nthis is the Lindemann formula
	// -- If LOW/ then a line with TROE/ follows then this is Troe's formula
	if(haspdep){
	  std::cout << "\tWARNING: pressure dependent reaction\n";       
	  std::cerr << "\tWARNING: pressure dependent reaction\n";
	  std::cout << "\tWARNING: Proteus does not support pressure dependence at this time -- ignoring" << std::endl;
	  if(haslow)  std::cout << "\tLOW: " << LOW << std::endl;
	  if(hastroe) std::cout << "\tTROE: " << TROE << std::endl;
	}

	// If the next line contains the PLOG keyword then this may be used to express
	// pressure dependence at these tabuluar pressures
	ireaction++;
	state = stateReadingRXNSList;
	break;
      case stateFoundEND:
	// we're done, just close out the file
	fin.close();
	BuildGlobalSpeciesList();
	return 0;
	break;
      default:
	std::cerr << "State not defined reading Chemkin" << std::endl;
	break;
      };
      // wait until we get to line REACTIONS to start looking for reactions
      
      // read reactions until we see END as well
    }
  }
  else{
    std::stringstream ss;
    ss << "CHEM_RXN CHEMKIN READFILE: Cannot open reaction file --> " << rxnfile << std::endl;
    Abort << ss.str();
    return 1;
  }
  fin.close();

  BuildGlobalSpeciesList();

  return 0;
}

template <class Type>
ChemModel<Type>::~ChemModel()
{
  delete [] reactions;
  delete [] species;
  for(typename std::vector<EOS<Type>* >::iterator it = eos.begin(); it != eos.end(); ++it){
    delete *it;
  }

  return;
}

template <class Type>
Int ChemModel<Type>::ReadNumberOfReactionsFromFile(std::string rxnfile)
{
  Int nreactions = 0;
  std::string fileName = rxnfile;
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
    std::stringstream ss;
    ss << "CHEM_RXN READFILE: Cannot open reaction file --> " << fileName << std::endl;
    Abort << ss.str();
    return 1;
  }

  return nreactions;
}

//This function is just an empty state machine which serves as a counter for Chemkin files
template <class Type>
Int ChemModel<Type>::ReadNumberOfReactionsFromFileChemkin(std::string rxnfile)
{

  std::ifstream fin;
  std::string line, trash;

  //states for our state machine reader
  enum{
    stateFindingRXNS,
    stateReadingRXNSList,
    stateReadingRXNItem1,
    stateReadingRXNItem2,
    stateReadingRXNItem3,
    stateFoundEND,
    NSTATES
  };

  Int state = stateFindingRXNS;
  Int nreactions = 0;
  
  fin.open(rxnfile.c_str());
  if(fin.is_open()){
    //We simply look through all lines knowing that =, =>, AND <=> represent a reaction, count 'em
    //Ignore anything past a ! in each line
    while(!fin.eof()){
      getline(fin, trash);
      // tokenize based on !, take everything before it
      std::vector<std::string> tokens = Tokenize(trash, '!');
      trash.clear();
      std::string line = tokens[0];

      if(line.find("REACTIONS") != std::string::npos) state = stateReadingRXNSList;
      else if((line.find("END") != std::string::npos) && state == stateReadingRXNSList) state = stateFoundEND;
      else if(line.find("=") != std::string::npos) state = stateReadingRXNItem1;
      else if(line.find("=>") != std::string::npos) state = stateReadingRXNItem2;
      else if(line.find("<=>") != std::string::npos) state = stateReadingRXNItem3;
      else{
	//do nothing just skip
      }
      switch(state){
      case stateFindingRXNS:
	// just ignore and keep reading
	break;
      case stateReadingRXNSList:
	break;
      case stateReadingRXNItem1:
      case stateReadingRXNItem2: //NOTE: we ignore irreversible reactions and treat them as reversible
      case stateReadingRXNItem3:
	nreactions++;
	state = stateReadingRXNSList;
	break;
      case stateFoundEND:
	return nreactions;
	break;
      default:
	std::cerr << "State not defined reading Chemkin" << std::endl;
	break;
      };
      // wait until we get to line REACTIONS to start looking for reactions

      // read reactions until we see END as well
    }
  }
  else{
    std::stringstream ss;
    ss << "CHEM_RXN CHEMKIN READFILE: Cannot open reaction file --> " << rxnfile << std::endl;
    Abort << ss.str();
    return -1;
  }
  fin.close();
  return nreactions;
  
}

// returns fluid properties given rhoi's (kg/m^3), Temperature (K), and cvi's (J/kg.K)
// all properties passed in dimensional (SI)
// all properties passed back dimensional (SI)
template <class Type>
void ChemModel<Type>::GetFluidProperties(const Type* rhoi, const Type T, const Type* cvi,
					 Type& cv, Type& cp, Type& R, Type& P, Type& rho, Type& c2) const
{
  Int i;
  R = 0.0;
  P = 0.0;
  rho = 0.0;
  c2 = 0.0;

  for(i = 0; i < nspecies; ++i){
    rho += rhoi[i];
  }
  
  for(i = 0; i < this->nspecies; i++){
    Type Rs = species[i].R;
    R += rhoi[i]*Rs;
    // assumes that Dalton's law holds - P = sum_i (P_i)
    P += eos[i]->GetP(Rs, rhoi[i], T);
  }
  R /= rho;

  cp = GetCp(rhoi, T);
  cv = GetCv(rhoi, T);
  Type gamma = cp/cv;
  
  //TODO: generalize for non-ideal gas speed of sound
  //this simplification is done in Busby p.17, only valid for ideal gases
  c2 = gamma * R * T;
}

template <class Type>
void ChemModel<Type>::GetMassProductionRates(Type* rhoi, Type T, Type* wdot)
{
  Int i;
  Int j;
  for(i = 0; i < this->nspecies; i++){
    wdot[i] = 0.0;
    for(j = 0; j < this->nreactions; j++){
      wdot[i] += this->reactions[j].GetMassProductionRate(rhoi, this->nspecies, T, i);
    }
  }
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
  Type pert = 1.0e-8;
  Type rho;
  Type rhop;
  Type Rmix = 0.0;
  Type Rmixp;
  Type* rhoip = (Type*)alloca(sizeof(Type)*nspecies);
  Type* X = (Type*)alloca(sizeof(Type)*nspecies);
  Type* Xp = (Type*)alloca(sizeof(Type)*nspecies);
  Type* hi = (Type*)alloca(sizeof(Type)*nspecies);
  Type* hip = (Type*)alloca(sizeof(Type)*nspecies);
  Type Tpu = T + pert;
  Type Tpd = T - pert;
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
  //at constant temperature - dEt_drho(Tconst)
  for(i = 0; i < nspecies; i++){
    memcpy(rhoip, rhoi, sizeof(Type)*nspecies);
    rhoip[i] += pert;
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
    Ppu = eos[i]->GetP(Rmixp, rhop, T);
    //compute new total energy based on perturbed state
    Htpu = rhop*GetSpecificEnthalpy(Xp, T, hip) + rhop*momV2;
    Etpu = Htpu - Ppu;
    memcpy(rhoip, rhoi, sizeof(Type)*nspecies);
    rhoip[i] -= pert;
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
    Ppd = eos[i]->GetP(Rmixp, rhop, T);
    //compute new total energy based on perturbed state
    Htpd = rhop*GetSpecificEnthalpy(Xp, T, hip) + rhop*momV2;
    Etpd = Htpd - Ppd;
    //compute finite difference
    dEtdRhoi[i] = (Etpu - Etpd)/(2.0*pert);
  }

  Htpu = rho*GetSpecificEnthalpy(X, Tpu, hip) + rho*momV2;
  Htpd = rho*GetSpecificEnthalpy(X, Tpd, hip) + rho*momV2;
  Ppu = eos[0]->GetP(Rmix, rho, Tpu);
  Ppd = eos[0]->GetP(Rmix, rho, Tpd);
  Etpu = Htpu - Ppu;
  Etpd = Htpd - Ppd;
  dEtdT = (Etpu - Etpd)/(2.0*pert);

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
  Type T = eos[0]->GetT(Rmix, rho, P);
  
  Type dTdRho = eos[0]->GetdT_dRho(Rmix, rho, P, T);
  Type dTdR = eos[0]->GetdT_dR(Rmix, rho, P, T);
  
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
    Ppu = eos[i]->GetP(Rmixp, rhop, T);
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
    Ppd = eos[i]->GetP(Rmixp, rhop, T);
    //compute new total energy based on perturbed state
    Htpd = rhop*GetSpecificEnthalpy(Xp, T, hip) + rhop*momV2;
    Etpd = Htpd - Ppd;
    //compute finite difference
    dEtdRhoi[i] = (Etpu - Etpd)/(2.0*h);
  }

  Ppu = P+h;
  Ppd = P-h;
  Tpu = eos[0]->GetT(Rmix, rho, Ppu);
  Tpd = eos[0]->GetT(Rmix, rho, Ppd);
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
  Type dPdrhomix = eos[0]->GetdP_dRho(Rmix, rhomix, P, T);
  Type dPdRmix = eos[0]->GetdP_dR(Rmix, rhomix, P, T);

  //chain rule F.T.W. at constant T -- pick (rho, T)
  for(Int i = 0; i < nspecies; i++){
    dEtdRhoi[i] = hv2 + species[i].GetH(T) - (dPdrhomix + dPdRmix*dRmixdRhoi(rhoi, rhomix, i));
  }

  //now we are perturbing T -- pick (P, rho)
  Type dEtdT = 0.0;
  Type dPdT = eos[0]->GetdP_dT(Rmix, rhomix, P, T);
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
  Type dPdrhomix = eos[0]->GetdP_dRho(Rmix, rhomix, P, T);
  Type dPdRmix = eos[0]->GetdP_dR(Rmix, rhomix, P, T);

  //chain rule F.T.W. at constant T -- pick (rho, T)
  for(Int i = 0; i < nspecies; i++){
    dEtdRhoi[i] = hv2 + species[i].GetH(T) - (dPdrhomix + dPdRmix*dRmixdRhoi(rhoi, rhomix, i));
  }
  
  //DEFINITION: Et = sum_i (rhoi * hi(T)) - P
  
  //now we are perturbing P -- pick (rho, T)
  Type dEtdP = 0.0;
  Type dTdP = eos[0]->GetdT_dP(Rmix, rhomix, P, T);
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

  //compute mol fraction
  Type* massfrac = (Type*)alloca(sizeof(Type)*nspecies);
  Type* molefrac = (Type*)alloca(sizeof(Type)*nspecies);

  for(Int i = 0; i < nspecies; ++i){
    massfrac[i] = rhoi[i]/rho;
  }

  MassFractionToMoleFraction(massfrac, molefrac);
  
  Type* visc = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; i++){
    visc[i] = species[i].GetViscosity(T);
  }

  for(Int i = 0; i < nspecies; i++){
    Type wi = 0.0;
    Type fraci = molefrac[i];
    Type MWi = species[i].MW;
    for(Int j = 0; j < nspecies; j++){
      Type fracj = molefrac[j];
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

template <class Type>
void ChemModel<Type>::MassFractionToMoleFraction(Type* massfrac, Type* molefrac)
{
  //compute mol fraction following Cox's dissertation 
  Type summ = 0.0;
  for(Int i = 0; i < nspecies; ++i){
    molefrac[i] = (massfrac[i])/species[i].MW;
    summ += molefrac[i];
  }
  for(Int i = 0; i < nspecies; ++i){
    molefrac[i] /= summ;
  }
  //make sure we add up to one to rounding error
  summ = 0.0;
  for(Int i = 0; i < nspecies-1; ++i){
    summ += molefrac[i];
  }
  molefrac[nspecies-1] = 1.0 - summ;
}


template <class Type>
void ChemModel<Type>::MoleFractionToMassFraction(Type* molefrac, Type* massfrac)
{
  Type avgmolemass = 0.0;
  for(Int i = 0; i < nspecies; ++i){
    avgmolemass += molefrac[i]*species[i].MW;
  }
  for(Int i = 0; i < nspecies; ++i){
    massfrac[i] = molefrac[i]*species[i].MW/avgmolemass;
  }
  //make sure we add up to one to rounding error
  Type summ = 0.0;
  for(Int i = 0; i < nspecies-1; ++i){
    summ += massfrac[i];
  }
  massfrac[nspecies-1] = 1.0 - summ;
}

template <class Type>
void ChemModel<Type>::MassToMole(Type* rho, Type* moleConc)
{
  //converts from kg/m^3 to mol/m^3
  for(Int i = 0; i < nspecies; ++i){
    moleConc[i] = rho[i]/species[i].MW/1000.0; //gives units of mol/m^3
  }
}

template <class Type>
Type ChemModel<Type>::GetP(const Type* rhoi, const Type T) const
{
  Type P = 0.0;
  for(Int i = 0; i < this->nspecies; i++){
    Type Rs = species[i].R;
    // assumes that Dalton's law holds - P = sum_i (P_i)
    P += eos[i]->GetP(Rs, rhoi[i], T);
  }
  return P;
}

template <class Type>
Type ChemModel<Type>::GetCp(const Type* rhoi, const Type T) const
{
  Type rho = 0.0;
  Type* X = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; ++i){
    rho += rhoi[i];
  }
  for(Int i = 0; i < nspecies; ++i){
    X[i] = rhoi[i]/rho;
  }
  
  Type cp = 0.0;
  for(Int i = 0; i < nspecies; i++){
    cp += species[i].GetCp(T)*X[i];
  }
  return cp;
}

template <class Type>
Type ChemModel<Type>::GetCv(const Type* rhoi, const Type T) const
{
  Type rho = 0.0;
  Type* X = (Type*)alloca(sizeof(Type)*nspecies);
  Type* cvi = (Type*)alloca(sizeof(Type)*nspecies);
  for(Int i = 0; i < nspecies; ++i){
    rho += rhoi[i];
  }
  for(Int i = 0; i < nspecies; ++i){
    X[i] = rhoi[i]/rho;
  }
  
  Type cv = 0.0;
  Type P = 0.0; //dummy variable
  for(Int i = 0; i < nspecies; i++){
    Type cpi = species[i].GetCp(T);
    Type Ri = species[i].R;
    cv += X[i]*(eos[i]->GetCv(cpi, Ri, rhoi[i], P, T));
  }
  return cv;
}
  
