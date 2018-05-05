//used to create an individual reaction
template <class Type>
Reaction<Type>::Reaction() :
  backwardRateGiven(false), thirdBodiesPresent(false),
  rxnType(999), rxnTypeBackward(999), species(NULL)
{}

template <class Type>
Reaction<Type>::~Reaction()
{}

template <class Type>
Int Reaction<Type>::GetNspecies()
{
  return speciesSymbols.size();
}

template <class Type>
void Reaction<Type>::Print()
{
  Int i;
  Int indx;

  //print out data about the reaction which was just created
  std::cout << "REACTION INITIALIZED: " << std::endl;
  std::cout << "( " << GetFormattedReaction() << " )" << std::endl;
  if(rxnType >= NUM_RXN_TYPES){
    std::stringstream ss;
    ss << "Reaction type " << rxnType << " not valid" << std::endl;
    Abort << ss.str();
  }
  if(backwardRateGiven && rxnTypeBackward >= NUM_RXN_TYPES){
    std::stringstream ss;
    ss << "Backward reaction type " << rxnTypeBackward << " not valid" << std::endl;
    Abort << ss.str();
  }

  std::cout << "========================" << std::endl;
  std::cout << "nspecies: " << this->GetNspecies() << std::endl;
  for(i = 0; i < this->GetNspecies(); i++){
    indx = globalIndx[i];
    std::cout << "Species[" << i << "] = " << species[indx].symbol << "   \t\tNup: " << Nup[i] << "\tNupp: " << Nupp[i] << std::endl;
  }
  if(thirdBodiesPresent){
    for(i = 0; i < this->GetNspecies(); i++){
      indx = globalIndx[i];
      std::cout << "Thirdbody[" << i << "] = " << species[indx].symbol 
		<< "  \tTBEff: " << TBEff[i] << std::endl;
    }
  }
  
  std::cout << "Temp(K) \tKf \t\tKc \t\tKb" << std::endl;
  std::cout << "===============================================================" << std::endl;
  std::cout.setf(std::ios::scientific);
  std::cout.precision(8);
  for(Int T = 200; T < 6000; T+=1100){
    std::cout << T 
	      << "\t\t" << GetForwardReactionRate((Real)T) 
	      << "\t" << GetEquilibriumReactionRate((Real)T)
	      << "\t" << GetBackwardReactionRate((Real)T) << std::endl;
  }  
  std::cout << "***************************************************************" << std::endl;
  std::cout << "Forward rate given in " << rxnTypeNames[rxnType] << std::endl;
  if(backwardRateGiven){
    std::cout << "Backward rate given in " << rxnTypeNames[rxnTypeBackward] << std::endl;
  }
  std::cout << "***************************************************************" << std::endl;
  std::cout.unsetf(std::ios::floatfield);

  return;
}

// This function iterates through the files and counts for the looks for reactionId given
// and passes the reaction input block off to ParseReactionBlockLine() for parsing
// filename - reaction file name to open
// reactionid - number of the reaction to ingest
// units - SI, CGS, etc. see chem_constants header for values
template <class Type>
Int Reaction<Type>::ReadReactionFromFile(std::string fileName, Int reactionId, Int units)
{
  Int err = 0;
  Int err2 = 0;
  
  std::ifstream fin;
  std::string line, trash;
  size_t loc;
  char c;

  std::stringstream ss;
  std::string rxnKey = "REACTION ";
  std::string rxnStart = rxnKey;
  
  //build descriptor for the reaction we are looking for
  ss.clear();
  ss << reactionId;
  rxnStart += ss.str();

  fin.open(fileName.c_str());
  if(fin.is_open()){
    std::cout << "CHEM_RXN: Reading reaction file --> " << fileName << " for rxn " 
	      << reactionId << std::endl;

    while(!fin.eof() && err != 1){
      c = fin.peek();
      //get rid of comments and blank lines
      if(c == '#' || c == ' ' || c == '\n'){
	getline(fin, trash);
	trash.clear();
	continue;
      }
      //get data from real lines
      else{
	//continue to scan through file until we find the reaction
	//we are looking for or hit the end of the file
	do{
	  getline(fin, line);
	  loc = line.find(rxnStart);
	}
	while(loc == std::string::npos && !fin.eof());
	
	if(fin.eof()){
	  std::stringstream ss;
	  ss << "CHEM_RXN: Reached end of file before finding rxn " << reactionId << std::endl;
	  Abort << ss.str();
	  err = 2;
	  return err;
	}
	
	//continue to scan through file until we get to the next reaction
	//or hit the end of the file, this is wasting one line read
	//but the line parser won't do anything with it, so it's safe
	do{
	  c = fin.peek();
	  //get rid of comments and blank lines
	  if(c == '#' || c == ' ' || c == '\n'){
	    getline(fin, trash);
	    trash.clear();
	    continue;
	  }
	  getline(fin, line);
	  loc = line.find(rxnKey);
	  err = ParseReactionBlockLine(line);
	  if(err != 0 && err != 1){
	    Abort << "CHEM_RXN: could not parse line -- " + line + " in file " 
		  + fileName + "\nCHEM_RXN: check that correct keyword is used";
	    err2 = err;
	    line.clear();
	  }
	  if(err == 1){
	    //reached the next reaction...
	    break;
	  }
	}
	while(loc == std::string::npos && !fin.eof());

      }
    }
    
    fin.close();
  }
  else{
    Abort << "CHEM_RXN READFILE: Cannot open reaction file --> " + fileName;
    return (1);
  }

  if(units == UNITS::CGS){
    this->ConvertArrheniusCGSToSI();
  }
  
  return err2;
}

// This parses a general line found in a reaction block based on keyword
// An example of such a block is
// REACTION 1
// rxnType = 2
// thirdBodiesPresent = 1
// O2 + M <==> 2O + M
// A = 3.618E12
// EA = 5.94E+4
// n = -1.000
// M = N[1.0], N2[2.0], O2[9.0], NO[1.0], O[25.0]
// This block builds up a reaction
template <class Type>
Int Reaction<Type>::ParseReactionBlockLine(std::string& line)
{
  size_t loc;

  std::stringstream ss;
  std::string rxn = "REACTION ";
  std::string subline;

  if(line.length() == 0){
    return(0);
  }

  //if we find the reaction key (i.e. the next reaction) we've read far enough to get all the current reaction's data
  loc = line.find("REACTION");
  if(loc != std::string::npos && loc == 0){
    return(1);
  }

  loc = line.find("<==>");
  if(loc != std::string::npos){
    ParseReaction(line);
    return(0);
  }
  loc = line.find("rxnTypeBackward");
  if(loc != std::string::npos){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> rxnTypeBackward;
    ss.clear();
    backwardRateGiven = true;
    return(0);
  }
  loc = line.find("rxnType");
  if(loc != std::string::npos){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> rxnType;
    ss.clear();
    return(0);
  }
  loc = line.find("thirdBodiesPresent");
  if(loc != std::string::npos){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> thirdBodiesPresent;
    ss.clear();
    return(0);
  }
  loc = line.find("Ab");
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> Ab;
    ss.clear();
    return(0);
  }
  loc = line.find("nb");
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> nb;
    ss.clear();
    return(0);
  }
  loc = line.find("EAb");
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> EAb;
    ss.clear();
    return(0);
  }
  loc = line.find("A");
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> A;
    ss.clear();
    return(0);
  }
  loc = line.find("EA");
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> EA;
    ss.clear();
    return(0);
  }
  loc = line.find("n");
  if(loc != std::string::npos){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> n;
    ss.clear();
    return(0);
  }
  loc = line.find("M");
  if(loc != std::string::npos){
    ParseThirdBodies(line);
    ss.clear();
    return(0);
  }

  return (2);
}

// Function parses the reaction format which looks like (spaces allowed and also not)
// O2 + M <==> 2O + M
// OR
// O + M <==> O(+) + E(-) + M
// where (+) (-) indicate ionized species and M is a third-body catalyst
template <class Type>
void Reaction<Type>::ParseReaction(std::string& rxn)
{
  Int i, j, k;
  Int nlhs, nrhs;
  std::string lhs, rhs;

  Int ions = 0;
  std::string ionp = "(+)";
  std::string ionn = "(-)";
  std::stringstream ss;
  std::string iff = "<==>";
  std::string elem;

  std::vector<std::string> tokens = Split(rxn, iff);
  lhs = tokens[0];
  rhs = tokens[1];
  
  //count species on lhs
  nlhs = std::count(lhs.begin(), lhs.end(), '+') + 1;
  nlhs -= std::count(lhs.begin(), lhs.end(), 'M');
  //subtract off the (+) signs which are referencing ions
  nlhs -= CountSubstring(lhs, "(+)");

  //count species on rhs
  nrhs = std::count(rhs.begin(), rhs.end(), '+') + 1;
  nrhs -= std::count(rhs.begin(), rhs.end(), 'M');
  //subtract off the (+) signs which are referencing ions
  nrhs -= CountSubstring(rhs, "(+)");

  //this is an estimate we will get the real number by parsing
  //out duplicates later
  Int nspecies = nrhs+nlhs;

  //allocate nupp and nup temporarily
  std::vector<Type> NupTemp(nspecies);
  std::vector<Type> NuppTemp(nspecies);

  //allocate room for species symbols
  std::vector<std::string> speciesSymbolsTemp(nspecies);

  std::vector<std::string> tlhstokens = Tokenize(lhs, ' ');
  std::vector<std::string> trhstokens = Tokenize(rhs, ' ');

  //remove any '+' tokens
  tlhstokens.erase(std::remove(tlhstokens.begin(), tlhstokens.end(), "+"), tlhstokens.end());
  trhstokens.erase(std::remove(trhstokens.begin(), trhstokens.end(), "+"), trhstokens.end());
  //remove any third body tokens "M"
  tlhstokens.erase(std::remove(tlhstokens.begin(), tlhstokens.end(), "M"), tlhstokens.end());
  trhstokens.erase(std::remove(trhstokens.begin(), trhstokens.end(), "M"), trhstokens.end());
  
  k = 0;
  for(std::vector<std::string>::const_iterator it = tlhstokens.begin(); it != tlhstokens.end(); ++it){
    std::string token = *it;
    token.erase(remove_if(token.begin(), token.end(), isspace), token.end());
    std::stringstream stoich;
    std::stringstream chem;
    bool endstoich = false;
    for(Int i = 0; i < (Int)token.size(); ++i){
      if(!isalpha(token[i]) && !endstoich){
	stoich << token[i];
      }
      else{
	endstoich = true;
	chem << token[i];
      }
    }
    //if there is no stoich coefficient read set it to one
    if(stoich.str().size() == 0) stoich << 1;
    std::string elem = chem.str();
    //strip out our custom ion notation (+)
    ions = StripFromString(elem, ionp, elem);
    if(ions){
      elem += "+";
    }
    //strip out our custom ion notation (-)
    ions = StripFromString(elem, ionn, elem);
    if(ions){
      elem += "-";
    }
    stoich >> NupTemp[k];
    NuppTemp[k] = 0;
    speciesSymbolsTemp[k] = elem;
    k++;
  }
  for(std::vector<std::string>::const_iterator it = trhstokens.begin(); it != trhstokens.end(); ++it){ 
    std::string token = *it;
    token.erase(remove_if(token.begin(), token.end(), isspace), token.end());
    std::stringstream stoich;
    std::stringstream chem;
    bool endstoich = false;
    for(Int i = 0; i < (Int)token.size(); ++i){
      if(!isalpha(token[i]) && !endstoich){
	stoich << token[i];
      }
      else{
	endstoich = true;
	chem << token[i];
      }
    }
    if(stoich.str().size() == 0) stoich << 1;
    std::string elem = chem.str();
    //strip out our custom ion notation (+)
    ions = StripFromString(elem, ionp, elem);
    if(ions){
      elem += "+";
    }
    //strip out our custom ion notation (-)
    ions = StripFromString(elem, ionn, elem);
    if(ions){
      elem += "-";
    }
    stoich >> NuppTemp[k];
    NupTemp[k] = 0;
    speciesSymbolsTemp[k] = elem;
    k++;
  }

  //eliminate possible duplicates from species list
  Int flag[nspecies];
  Int unique;

  //initialize flag
  for(i = 0; i < nspecies; i++){
    flag[i] = -1;
  }
  //first value is unique
  unique = 0;
  for(i = 0; i < nspecies; i++){
    if(flag[i] < 0){
      flag[i] = unique;
      for(j = i+1; j < nspecies; j++){
	if(speciesSymbolsTemp[i] == speciesSymbolsTemp[j]){
	  flag[j] = unique;
	}
      }
      unique++;
    }
  }

  //allocate final value for species symbols
  speciesSymbols.resize(unique);
  Nup.resize(unique);
  Nupp.resize(unique);
  

  //move unique values to beginning of list
  //first value always checks as unique
  //also sum up nup, and nupp
  for(i = 0; i < unique; i++){
    Nup[i] = Nupp[i] = 0;
    for(j = 0; j < nspecies; j++){
      if(flag[j] == i){
	speciesSymbols[i] = speciesSymbolsTemp[j];
	Nup[i] += NupTemp[j];
	Nupp[i] += NuppTemp[j];
	flag[j] = -1;
      }
    }
  }

  //rename electrons to lower case e-
  for(i = 0; i < unique; i++){
    if(speciesSymbols[i] == "E-"){
      speciesSymbols[i] = "e-";
    }
  }

  //now set total number of species which are unique
  nspecies = unique;
}

//we need to put the third bodies in the following form for parsing
//M = NO(+)[1.0], O2(+)[1.0], N2(+)[1.0], O(+)[1.0], N(+)[1.0]
template <class Type>
void Reaction<Type>::ParseThirdBodies(std::string& line)
{
  std::string TBSym;
  std::stringstream ss;
  Int ions;
  std::string ionp = "(+)";
  std::string ionn = "(-)";

  line = Split(line, "=")[1];
  
  //count number of csv 
  Int ntb = std::count(line.begin(), line.end(), ',') + 1;

  //if third bodies are used, all unlisted species have an efficiency of unity, allocate memory
  TBEff.resize(this->GetNspecies(), 1.0);

  Int i = 0;
  while(i < (Int)line.size()){
    //if we find an upper case letter (chemical), parse the catalyst one letter at a time
    TBSym.clear();
    TBSym = "";
    if(isupper(line[i])){
      while((line[i] != '[') && (line[i] != ' ')){
	TBSym += line[i];
	i++;
      }
      //strip out our custom ion notation (+)
      ions = StripFromString(TBSym, ionp, TBSym);
      if(ions){
	TBSym += "+";
      }
      //strip out our custom ion notation (-)
      ions = StripFromString(TBSym, ionn, TBSym);
      if(ions){
	TBSym += "-";
      }
      //rename electrons to lower case e-
      if(TBSym == "E-"){
	TBSym = "e-";
      }
      
      //skip to next brackets and read out the values to TBefficiency array
      size_t loc = line.find('[', i-1);
      loc += 1;
      size_t loc2 = line.find(']', i-1);
      if(loc == std::string::npos || loc2 == std::string::npos){
	std::cerr << "Error reading thirdbodies, matching brackets: " + line << std::endl;
      }
      std::string val = line.substr(loc, loc2-loc);
      ss.clear();
      ss << val;
      Type vald;
      ss >> vald;
      
      //find the species in the species list for the reaction
      Int specielocalid = -99;
      for(int k = 0; k < speciesSymbols.size(); ++k){
	if(TBSym == speciesSymbols[k]){
	  specielocalid = k;
	  TBEff.at(specielocalid) = vald;
	  break;
	}
      }
      //if a new species shows up as a third body, extend both lists
      if(specielocalid == -99){
	speciesSymbols.push_back(TBSym);
	TBEff.push_back(vald);
	Nup.push_back(0.0);
	Nupp.push_back(0.0);
      }
      i = loc2;
      i++;
    }
    i++;
  }
  
  if(TBEff.size() != speciesSymbols.size()){
    Abort << "Third bodies and species symbols list sizes do not match";
  }
}

//looks through global species symbol list and sets the global indx locally
//globalSpecies - list of global species pointers (globally complete)
//globalCount - number of species in the list
template <class Type>
void Reaction<Type>::SetSpeciesPointer(Species<Type>* globalSpecies, Int globalCount)
{
  //copy a pointer to point at the global list
  species = globalSpecies;
  globalIndx.resize(this->GetNspecies());
  
  for(Int i = 0; i < this->GetNspecies(); i++){
    for(Int j = 0; j < globalCount; j++){
      if(speciesSymbols[i] == species[j].symbol){
	globalIndx[i] = j;
	break;
      }
    }
  }
}

//the function maps the species id from the global list (all reactions) to the
//list that is specific to this reaction only
//speciedId - local species indx from reaction
//returns - global species list indx
template <class Type>
Int Reaction<Type>::GetLocalSpeciesFromGlobal(Int speciesId)
{
  for(Int i = 0; i < this->GetNspecies(); i++){
    //the first occurence in the list is what we want
    if(globalIndx[i] == speciesId){
      return(i);
    }
  }
  return(-1);
}

//returns forward reaction rate Kf with units (kmol/m^3)^nu
template <class Type>
Type Reaction<Type>::GetForwardReactionRate(Type T)
{
  if(rxnType == ArrheniusCoeff){
    return Arrhenius(A,EA,T);
  }
  else if(rxnType == ModArrheniusCoeff){
    return ModArrhenius(A,EA,n,T);
  }
  else if(rxnType == GuptaModArrheniusCoeff){
    return GuptaModArrhenius(A, EA, n, T);
  }
  else if(rxnType == Power){
    return PowerRxn(A, n, T);
  }
  
  return(-999);
}

//returns equilibrium reaction rate Kc with no units (dimensionless iff number of moles of products == number of moles of reactants)
//will have units if the mole balance is not zero from products to reactants
template <class Type>
Type Reaction<Type>::GetEquilibriumReactionRate(Type T)
{
  Type dnu_i;
  Type a[7];
  Type nu = 0.0;
  Type dnua1 = 0.0;
  Type dnua2 = 0.0;
  Type dnua3 = 0.0;
  Type dnua4 = 0.0;
  Type dnua5 = 0.0;
  Type dnua6 = 0.0;
  Type dnua7 = 0.0;
  
  //calculate the delta in reference enthalpy and reference entropy
  for(Int i = 0; i < this->GetNspecies(); i++){
    Int ii = globalIndx[i];
    species[ii].GetThermoCoeff(T,a);
    dnu_i = Nupp[i] - Nup[i];
    nu += dnu_i;
    dnua1 += dnu_i*a[0];
    dnua2 += dnu_i*a[1];
    dnua3 += dnu_i*a[2];
    dnua4 += dnu_i*a[3];
    dnua5 += dnu_i*a[4];
    dnua6 += dnu_i*a[5];
    dnua7 += dnu_i*a[6];
  }
    
  //compute the equilibrium constant based on pressure
  //p. 502, Gardiner
  //This is DeltaS_rr/UnivR - DeltaH_rr/(UnivR*T)
  //p. 27 in my thesis(Nicholas Currier), exp(-Gibbs free energy)
  Type Kp = exp(dnua1*(log(T) - 1.0) + T*(dnua2/2.0 + T*(dnua3/6.0 + T*(dnua4/12.0 + dnua5/20.0*T))) 
		- dnua6/T + dnua7);
  
  // compute the equilibrium constant based on concentration.  Equation 9.91, Kee
  Type Kc = Kp;
  if(!(CAbs(real(nu)) < 1.0e-15)){
    Type Pref = 101325.0;
    //attempt to use the pow(Type, int) version if possible - much faster
    if(isWholeNumber(nu)){
      Int nuint = (Int)real(nu);
      Kc *= pow(Pref/(UNIV_R*T),nuint); // has units of (kmol/m^3)^nu, 
    }
    else{
      Kc *= pow(Pref/(UNIV_R*T), nu); // has units of (kmol/m^3)^nu
    }
  }
  
  return Kc;  //Unitless if number of moles of products == number of moles of reactants
}

//returns backwards reaction rate Kb with units (1/s)(kmol/m^3)^(1-zpp)
template <class Type>
Type Reaction<Type>::GetBackwardReactionRate(Type T)
{
  if(backwardRateGiven){
    if(rxnTypeBackward == ArrheniusCoeff){
      return Arrhenius(Ab,EAb,T);
    }
    else if(rxnTypeBackward == ModArrheniusCoeff){
      return ModArrhenius(Ab,EAb,nb,T);
    }
    else if(rxnTypeBackward == GuptaModArrheniusCoeff){
      return GuptaModArrhenius(Ab, EAb, nb, T);
    }
    else if(rxnTypeBackward == Power){
      return PowerRxn(Ab, nb, T);
    }
    return(-999);
  }
  
  //if a backward reaction form is not given assume we have to 
  //calculate it based on equilibrium rate constants
  Type Kf = GetForwardReactionRate(T);
  Type Kc = GetEquilibriumReactionRate(T);  
  return (Kf/Kc);
}

//returns reaction rate in units (1/s)(kmol/m^3)^(1-z')
//T - temperature in Kelvin
//Ea - activation energy in J/mol
//A - Prefactor in  units of (1/sec)(kmol/m^3)^(1-z') and kbr (1/sec)(kmol/m^3)^(1-z'') internally
//    -- z' is the sum of the left-hand side stoich coefficients(nu')
//    -- z'' is the sum of the right-hand size stoich coefficients (nu'')
template <class Type>
Type Reaction<Type>::Arrhenius(Type A, Type EA, Type T)
{
  //UNIV_R - J/(mol.K)
  return(A * exp(-EA/(UNIV_R * T)));
}


//returns reaction rate in units (1/s)(kmol/m^3)^(1-z')
//T - temperature in Kelvin
//Ea - activation energy in J/mol
//A - Prefactor in  units of (1/sec)(kmol/m^3)^(1-z') and kbr (1/sec)(kmol/m^3)^(1-z'') internally
//    -- z' is the sum of the left-hand side stoich coefficients(nu')
//    -- z'' is the sum of the right-hand size stoich coefficients (nu'')
template <class Type>
Type Reaction<Type>::ModArrhenius(Type A, Type EA, Type n, Type T)
{
  return(A * pow(T, n) * exp(-EA/(UNIV_R * T)));
}

//returns reaction rate in units (1/s)(kmol/m^3)^(1-z')
//T - temperature in Kelvin
//Ea - activation energy in J/mol
//A - Prefactor in  units of (1/sec)(kmol/m^3)^(1-z') and kbr (1/sec)(kmol/m^3)^(1-z'') internally
//    -- z' is the sum of the left-hand side stoich coefficients(nu')
//    -- z'' is the sum of the right-hand size stoich coefficients (nu'')
template <class Type>
Type Reaction<Type>::GuptaModArrhenius(Type A, Type EA, Type n, Type T)
{
  return(A * pow(T, n) * exp(-EA/T));
}

//returns reaction rate in units (1/s)(kmol/m^3)^(1-z')
//T - temperature in Kelvin
//Ea - activation energy in J/mol
//A - Prefactor in  units of (1/sec)(kmol/m^3)^(1-z') and kbr (1/sec)(kmol/m^3)^(1-z'') internally
//    -- z' is the sum of the left-hand side stoich coefficients(nu')
//    -- z'' is the sum of the right-hand size stoich coefficients (nu'')
template <class Type>
Type Reaction<Type>::PowerRxn(Type A, Type n, Type T){
  return(A * pow(T, n));
}

//returns production rate of given species in kg/s.m^3
//rhoi - dimensional vector of densities for all species (global) kg/m^3
//globalNspecies - number of entries in the rhoi vector
//T - temperature in Kelvin
//speciesId - global species indx (will be mapped to local species indx)
template <class Type>
Type Reaction<Type>::GetMassProductionRate(Type* rhoi, Int globalNspecies, Type T, Int speciesId)
{
  Type wdot = 0.0;
  Type dest,form,prod_dest,prod_form;
  Type nupp_kr,nup_kr;
  Type nu_ir;
  Type Kf,Kb;
  Type Mconc,net;
  Int nspecies = this->GetNspecies();
  std::vector<Type> X(globalNspecies); //precomputed concentrations for each species
  Int k,kk;	  
  
  kk = GetLocalSpeciesFromGlobal(speciesId);
  if(kk != -1){
    nu_ir = Nupp[kk] - Nup[kk];
  }
  else{
    //we don't have that species in this reaction
    nu_ir = 0;
  }

  Real zero = 0.0;
  Real diff = 1.0e-15;
  Real dnu_ir = real(nu_ir);
  //this kicks out so we don't do the expensive math below if coefficients are identical
  if (AlmostEqualRelative(dnu_ir, zero, diff)) return 0.0;
  
  //precompute the concentrations for all species in the gas mixture globally
  for(k = 0; k < globalNspecies; k++){
    X[k] = rhoi[k]/species[k].MW; //(kg/m^3)/(kg/kmol) => (kmol/m^3)
  }

  //Kfr has units of (1/sec)(kmol/m^3)^(1-z') and kbr (1/sec)(kmol/m^3)^(1-z'')
  //z' is the sum of the left-hand side stoich coefficients and z'' is the sum of the right-hand size stoich coefficients
  Kf = GetForwardReactionRate(T);
  Kb = GetBackwardReactionRate(T);

  prod_form = prod_dest = 1.0;
  
  //loop over the species that participate locally as reactants
  for(k = 0; k < nspecies; k++){
    //need this to reference the correction global concentration X[kk]
    kk = globalIndx[k];

    nupp_kr = Nupp[k];
    nup_kr  = Nup [k];

    //avoid using pow(Type, Type) if possible to use pow(Type, Int) - much cheaper
    if(isWholeNumber(nup_kr)){
      Int nup_kri = (Int)real(nup_kr);
      form = pow(X[kk],nup_kri); //(kmol/m^3)^(nu')
    }
    else{
      form = pow(X[kk], nup_kr); //(kmol/m^3)^(nu')
    }
    if(isWholeNumber(nupp_kr)){
      Int nupp_kri = (Int)real(nupp_kr);
      dest = pow(X[kk], nupp_kri); //(kmol/m^3)^(nu'')
    }
    else{
      dest = pow(X[kk],nupp_kr);   //(kmol/m^3)^(nu'')
    }

    prod_form *= form;
    prod_dest *= dest;
  }

  //if third body action is turned on we need to calculate the effective
  //concentration of the third bodies (equation 9.75, p. 384, Kee text)
  //Note that for radical reactions, all species are assumed to have an efficiency
  //of unity unless they show up in our input list as something else
  //this is b/c the third body collision provides the energy for the reaction
  //to proceed
  Mconc = 1.0;
  if(thirdBodiesPresent){
    Mconc = 0.0;
    //sum up concentration for all gases in the mixture globally - assum TBeff = 1.0
    for(k = 0; k < globalNspecies; ++k){
      Mconc += X[k];
    }
    //for gases that we have local information about, exchange their TEff*X concentrations instead
    for(k = 0; k < nspecies; ++k){
      kk = globalIndx[k];
      Mconc += (TBEff[k]-1.0)*X[kk]; //kmol/m^3
    }
  }
  //Mconc is the same as \Gamma in my disseration - Nick
  
  net  = Kf*prod_form - Kb*prod_dest;  // (kmol/m^3.s)
  wdot = (Type)nu_ir*Mconc*net;        // (kmol/m^3.s)(N.D.) in kmol/s.m^3
  wdot *= species[speciesId].MW;       // convert units (kmol/s.m^3)(kg/kmol) to (kg/s.m^3)

  return wdot; //(kg/s.m^3)
}

//returns a string with the formatted reaction
template <class Type>
std::string Reaction<Type>::GetFormattedReaction()
{
  Int i, indx; 
  Int first;
  std::string rxn = "";
  std::stringstream ss;
  Int nspecies = this->GetNspecies();
  
  first = 1;
  for(i = 0; i < nspecies; i++){
    indx = globalIndx[i];
    if(Nup[i] != (Type)0){
      if(!first){
	ss << "+ ";
      }
      if(Nup[i] == (Type)1){
	ss << species[indx].symbol << " ";
      }
      else{
	ss << Nup[i] << species[indx].symbol << " ";
      }
      first = 0;
    }
  }
  if(thirdBodiesPresent){
    ss << "+ M ";
  }
  ss << "<==> ";
  first = 1;
  for(i = 0; i < nspecies; i++){
    indx = globalIndx[i];
    if(Nupp[i] != (Type)0){
      if(!first){
	ss << "+ ";
      }
      if(Nupp[i] == (Type)1){
	ss << species[indx].symbol << " ";
      }
      else{
	ss << Nupp[i] << species[indx].symbol << " ";
      }
      first = 0;
    }
  }
  if(thirdBodiesPresent){
    ss << "+ M";
  }
  rxn = ss.str();

  return rxn;
}

//Function converts Arrhenius coefficients from cm, gram, moles to SI units for internal use
//Chemkin gives us rates in the form k = A T^n exp(-Ea/(Runiv.T))
//A - prefactor - units in cgs are (mol/(cm^3.s.K)) -- > we need something like (kmol/(m^3.s.K))
//    note that unimolecular reactions are 1/s
//              bimolecular reactions are mol/(cm^3.s.K)
//              trimolecular reactions are mol^2/(cm^6.s.K)
//n - power exponent
//Ea - activation energy (chemkin default cal/mol unless keyword Joule/mole is in file)
//zp_zpp - if forward rate constant is given pass z', if backward is given pass z''
//         z' or Sum(nu'), sum of left-hand side stoich coefficients
//         z''or Sum(nu''), sum of right-hand side stoich coefficients
template <class Type>
void Reaction<Type>::ConvertArrheniusCGSToSI()
{
  Type zp_zpp = 0.0;
  for(Int i = 0; i < Nup.size(); ++i){
    zp_zpp += Nup[i];
    //TODO: support backwards reactions specification with Nupp
  }
  //Ea is passed to us in calories/mole (not KCal) convert it to Joules/mole
  EA = EA * (Type)CAL_TO_JOULE;
  Real inverseOfcm3Tom3 = 1.0e-6;
  Real molTokmol = 1.0e-3;
  Real factor = molTokmol*inverseOfcm3Tom3;
  Type one = 1.0;
  A = A*pow(factor, (one - zp_zpp));
}

