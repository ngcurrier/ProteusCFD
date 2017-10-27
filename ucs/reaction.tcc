//used to create an individual reaction
template <class Type>
Reaction<Type>::Reaction() :
  ncatalysts(0), ntbodies(0), nspecies(0), backwardRateGiven(false),
  rxnType(999), rxnTypeBackward(999), TBEff(NULL), species(NULL), globalIndx(NULL), speciesSymbols(NULL),
  catSymbols(NULL), Nup(NULL), Nupp(NULL)
{
}

template <class Type>
Reaction<Type>::~Reaction()
{
  delete [] Nup;
  delete [] Nupp;
  delete [] speciesSymbols;
  delete [] catSymbols;
  if(TBEff != NULL){
    delete [] TBEff;
  }
  delete [] globalIndx;
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
  std::cout << "nspecies: " << nspecies << std::endl;
  for(i = 0; i < nspecies; i++){
    indx = globalIndx[i];
    std::cout << "Species[" << i << "] = " << species[indx].symbol << "   \t\tNup: " << Nup[i] << "\tNupp: " << Nupp[i] << std::endl;
  }
  std::cout << "ncatalysts: " << ncatalysts << std::endl;
  for(i = nspecies; i < ntbodies+ncatalysts+nspecies; i++){
    indx = globalIndx[i];
    std::cout << "Catalyst[" << i-nspecies << "] = " << species[indx].symbol 
	      << "  \tTBEff: " << TBEff[i-nspecies] << std::endl;
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

template <class Type>
Int Reaction<Type>::ReadReactionFromFile(std::string fileName, Int reactionId)
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

  return err2;
}

// This parse a general line found in a reaction block based on keyword
// An example of such a block is
// REACTION 1
// rxnType = 2
// catPresent = 1
// reactantsAsCat = 1
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
  std::string rxnStart = rxn;
  std::string iff = "<==>";
  std::string ArrCoeffA = "A";
  std::string ArrCoeffEA = "EA";
  std::string ArrCoeffn = "n";
  std::string CatalystSwitch = "catPresent";
  std::string TBodies = "reactantsAsCat";
  std::string rxnTypeK = "rxnType";
  std::string rxnTypeBackK = "rxnTypeBackward";
  std::string Catalyst = "M";
  std::string ArrCoeffAb = "Ab";
  std::string ArrCoeffnb = "nb";
  std::string ArrCoeffEAb = "EAb";

  std::string subline;

  if(line.length() == 0){
    return(0);
  }


  //if we find the reaction key we've read far enough to get all this
  //reaction's data
  loc = line.find(rxn);
  if(loc != std::string::npos && loc == 0){
    return(1);
  }

  loc = line.find(iff);
  if(loc != std::string::npos){
    ParseReaction(line);
    return(0);
  }
  loc = line.find(rxnTypeBackK);
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
  loc = line.find(rxnTypeK);
  if(loc != std::string::npos){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> rxnType;
    ss.clear();
    return(0);
  }
  loc = line.find(CatalystSwitch);
  if(loc != std::string::npos){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> catalysts;
    ss.clear();
    return(0);
  }
  loc = line.find(TBodies);
  if(loc != std::string::npos){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> thirdBodies;
    ss.clear();
    return(0);
  }
  loc = line.find(ArrCoeffAb);
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> Ab;
    ss.clear();
    return(0);
  }
  loc = line.find(ArrCoeffnb);
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> nb;
    ss.clear();
    return(0);
  }
  loc = line.find(ArrCoeffEAb);
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> EAb;
    ss.clear();
    return(0);
  }
  loc = line.find(ArrCoeffA);
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> A;
    ss.clear();
    return(0);
  }
  loc = line.find(ArrCoeffEA);
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> EA;
    ss.clear();
    return(0);
  }
  loc = line.find(ArrCoeffn);
  if(loc != std::string::npos){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> n;
    ss.clear();
    return(0);
  }
  loc = line.find(Catalyst);
  if(loc != std::string::npos){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ParseCatalysts(subline);
    ss.clear();
    return(0);
  }

  return (2);
}

// Function parses the reaction format which looks like
// O2 + M <==> 2O + M
// OR
// O + M <==> O(+) + E(-) + M
// where (+) (-) indicate ionized species and M is a third-body catalyst
template <class Type>
Int Reaction<Type>::ParseReaction(std::string& rxn)
{
  Int err = 0;
  Int i, j, k;
  size_t loc;
  Int nlhs, nrhs;
  std::string lhs, rhs;

  Int ions = 0;
  std::string ionp = "(+)";
  std::string ionn = "(-)";
  std::stringstream ss;
  std::string iff = "<==>";
  std::string elem;
  loc = rxn.find(iff);

  //count species on lhs
  lhs = rxn.substr(0, loc);
  nlhs = std::count(lhs.begin(), lhs.end(), '+') + 1;
  nlhs -= std::count(lhs.begin(), lhs.end(), 'M');
  //subtract off the (+) signs which are referencing ions
  nlhs -= CountSubstring(lhs, "(+)");

  //count species on rhs
  rhs = rxn.substr(loc+iff.size());
  nrhs = std::count(rhs.begin(), rhs.end(), '+') + 1;
  nrhs -= std::count(rhs.begin(), rhs.end(), 'M');
  //subtract off the (+) signs which are referencing ions
  nrhs -= CountSubstring(rhs, "(+)");

  //this is an estimate we will get the real number by parsing
  //out duplicates later
  nspecies = nrhs+nlhs;

  //allocate nupp and nup temporarily
  Type* NupTemp = new Type[nspecies];
  Type* NuppTemp = new Type[nspecies];

  //allocate room for species symbols
  std::string* speciesSymbolsTemp = new std::string[nspecies];

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
  speciesSymbols = new std::string[unique];
  Nup = new Type[unique];
  Nupp = new Type[unique];
  

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

  delete [] speciesSymbolsTemp;
  delete [] NupTemp;
  delete [] NuppTemp;

  return err;
}

template <class Type>
Int Reaction<Type>::ParseCatalysts(std::string& line)
{
  Int err = 0;
  Int i;
  size_t loc, loc2;
  std::string catSym;
  std::string str;
  std::stringstream ss;
  Int ions;
  std::string ionp = "(+)";
  std::string ionn = "(-)";
  
  //count number of csv 
  Int ncat = std::count(line.begin(), line.end(), ',') + 1;
  
  //strip off third bodies, i.e. rxn species acting as catalysts
  //we want a list of strict catatlysts ONLY
  ntbodies = 0;
  for(i = 0; i < nspecies; i++){
    //attempt to idiot proof the format in case spaces
    //are not put after the commas... this is easy to break
    //however, TODO: more robust parser here...
    loc = line.find(' '+speciesSymbols[i]+'[');
    if(loc != std::string::npos){
      ntbodies++;
    }
    loc = line.find(','+speciesSymbols[i]+'[');
    if(loc != std::string::npos){
      ntbodies++;
    }
  }
  if(!thirdBodies && ntbodies){
    Abort << "WARNING: third bodies switch is off but third bodies found in list!" +
      thirdBodies;
    err = 1;
  }

  //set ncatalysts equal to efficiencies found minus third bodies
  ncatalysts = ncat - ntbodies;
  catSymbols = new std::string[ncatalysts+ntbodies];
  TBEff = new Type[ncatalysts+ntbodies];

  i = 0;
  Int count = 0;
  while(i < (Int)line.size()){
    //if we find an upper case letter parse the catalyst
    if(isupper(line[i])){
      while((line[i] != '[') && (line[i] != ' ')){
	catSym += line[i];
	i++;
      }
      //strip out our custom ion notation (+)
      ions = StripFromString(catSym, ionp, catSym);
      if(ions){
	catSym += "+";
      }
      //strip out our custom ion notation (-)
      ions = StripFromString(catSym, ionn, catSym);
      if(ions){
	catSym += "-";
      }
      //rename electrons to lower case e-
      if(catSym == "E-"){
	catSym = "e-";
      }
      catSymbols[count] = catSym;
      catSym = "";
      
      //skip to next brackets and read out the values
      loc = line.find('[', i-1);
      loc += 1;
      loc2 = line.find(']', i-1);
      str = line.substr(loc, loc2-loc);
      ss.clear();
      ss << str;
      ss >> TBEff[count];
      i = loc2;

      count++;
      i++;
    }
    i++;
  }


  
  //for(i = 0; i < ncatalysts+ntbodies; i++){
  //  std::cout << catSymbols[i] << " " << TBEff[i] << std::endl;
  //}

  return err;
}

template <class Type>
void Reaction<Type>::SetSpeciesPointer(Species<Type>* globalSpecies, Int globalCount)
{
  Int i, j;
  
  species = globalSpecies;
  //this list will have some duplicates if ntbodies != 0
  globalIndx = new Int[nspecies+ncatalysts+ntbodies];
  
  for(i = 0; i < nspecies; i++){
    for(j = 0; j < globalCount; j++){
      if(speciesSymbols[i] == species[j].symbol){
	globalIndx[i] = j;
	break;
      }
    }
  }
  for(i = 0; i < ncatalysts+ntbodies; i++){
    for(j = 0; j < globalCount; j++){
      if(catSymbols[i] == species[j].symbol){
	globalIndx[i+nspecies] = j;
	break;
      }
    }
  }

  return;
}

//the function maps the species id from the global list (all reactions) to the
//list that is specific to this reaction only
template <class Type>
Int Reaction<Type>::GetLocalSpeciesFromGlobal(Int speciesId)
{
  Int i;
  for(i = 0; i < nspecies; i++){
    //the first occurence in the list is what we want
    //since reactants that act as catalysts are stored again
    //towards the end
    if(globalIndx[i] == speciesId){
      return(i);
    }
  }
  return(-1);
}

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

template <class Type>
Type Reaction<Type>::GetEquilibriumReactionRate(Type T)
{
  Int i,ii;
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
  for(i = 0; i < nspecies; i++){
    ii = globalIndx[i];
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
      Kc *= pow(Pref/(UNIV_R*T),nuint); // has units of (mol/m^3)^nu, 
    }
    else{
      Kc *= pow(Pref/(UNIV_R*T), nu); // has units of (mol/m^3)^nu
    }
  }
  
  return Kc;
}

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

template <class Type>
Type Reaction<Type>::Arrhenius(Type A, Type EA, Type T)
{
  return(A * exp(-EA/(UNIV_R * T)));
}

template <class Type>
Type Reaction<Type>::ModArrhenius(Type A, Type EA, Type n, Type T)
{
  return(A * pow(T, n) * exp(-EA/(UNIV_R * T)));
}

template <class Type>
Type Reaction<Type>::GuptaModArrhenius(Type A, Type EA, Type n, Type T)
{
  return(A * pow(T, n) * exp(-EA/T));
}

template <class Type>
Type Reaction<Type>::PowerRxn(Type A, Type n, Type T){
  return(A * pow(T, n));
}

template <class Type>
Type Reaction<Type>::GetMassProductionRate(Type* rhoi, Type T, Int speciesId)
{
  Type wdot = 0.0;
  Type dest,form,prod_dest,prod_form;
  Type nupp_kr,nup_kr;
  Type nu_ir;
  Type Kf,Kb;
  Type Mconc,net;
  Type* X = (Type*)alloca(sizeof(Type)*(nspecies+ncatalysts+ntbodies)); //precomputed concentrations for each species
  Type XX;
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
  if (AlmostEqualRelative(dnu_ir, zero, diff)) return 0.0;
  
  //precompute the concentrations
  for(k = 0; k < nspecies+ncatalysts+ntbodies; k++){
    kk = globalIndx[k];
    X[k] = rhoi[kk]/species[kk].MW; //(kg/m^3)/(kg/mol) => (mol/m^3)
  }
  
  Kf = GetForwardReactionRate(T);
  Kb = GetBackwardReactionRate(T);

  prod_form = prod_dest = 1.0;
  

  //loop over the species that participate as reactants
  for(k = 0; k < nspecies; k++){
    nupp_kr = Nupp[k];
    nup_kr  = Nup [k];

    //avoid using pow(Type, Type) if possible to use pow(Type, Int) - much cheaper
    if(isWholeNumber(nup_kr)){
      Int nup_kri = (Int)real(nup_kr);
      form = pow(X[k],nup_kri);
    }
    else{
      form = pow(X[k], nup_kr);
    }
    if(isWholeNumber(nupp_kr)){
      Int nupp_kri = (Int)real(nupp_kr);
      dest = pow(X[k], nupp_kri);
    }
    else{
      dest = pow(X[k],nupp_kr);
    }

    prod_form *= form;
    prod_dest *= dest;
  }

  //if catalyst action is turned on we need to calculate the effective
  //concentration of the catalysts (equation 9.75, p. 384, Kee text)
  Mconc = 1.0;
  if(catalysts){
    Mconc = 0.0;
    for(k = 0; k < ntbodies+ncatalysts; k++){
      kk = globalIndx[k+nspecies];
      XX = rhoi[kk]/species[kk].MW; //in mol/m^3
      Mconc += TBEff[k]*XX;
    }
  }
  
  net  = Kf*prod_form - Kb*prod_dest; 
  wdot = (Type)nu_ir*Mconc*net; //in mol/s.m^3
  wdot *= species[speciesId].MW;  //convert units to kg/s.m^3

  return wdot;
}

template <class Type>
std::string Reaction<Type>::GetFormattedReaction()
{
  Int i, indx; 
  Int first;
  std::string rxn = "";
  std::stringstream ss;

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
  if(catalysts){
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
  if(catalysts){
    ss << "+ M";
  }
  rxn = ss.str();

  return rxn;
}
