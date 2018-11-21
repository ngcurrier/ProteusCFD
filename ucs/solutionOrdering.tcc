template <class Type>
SolutionOrdering<Type>::~SolutionOrdering()
{
  for(typename std::vector<SolutionOperation<Type>*>::iterator it = opList.begin(); 
      it != opList.end(); ++it){
    delete *it;
  }
}

// returns the number of operations that we execute per time step
template <class Type>
UInt SolutionOrdering<Type>::size()
{
  return commands.size();
}

template <class Type>
Int SolutionOrdering<Type>::Insert(const std::string lineCommand)
{
  std::string tlineCommand = lineCommand;
  //eliminate any leading spaces
  size_t loc2 = 0;
  while(loc2 != std::string::npos){
    loc2 = tlineCommand.find(' ');
    if(loc2 != std::string::npos && loc2 == 0){
      tlineCommand = tlineCommand.substr(loc2+1);
      continue;
    }
    break;
  }
  //eliminate trailing spaces
  loc2 = 0;
  while(loc2 != std::string::npos){
    loc2 = tlineCommand.rfind(' ');
    if(loc2 != std::string::npos && loc2 == tlineCommand.size()-1){
      tlineCommand = tlineCommand.substr(0, loc2-1);
      continue;
    }
    break;
  }
  commands.push_back(tlineCommand);
  return 0;
}

template <class Type>
Int SolutionOrdering<Type>::Finalize(std::vector<SolutionSpaceBase<Type>*>& solSpaces)
{
  std::string iterateKey = "Iterate";
  std::string transferKey = "Transfer";

  for(std::vector<std::string>::iterator it = commands.begin(); it != commands.end(); ++it){
    //make a copy, we will modify this string
    std::string command = *it;
    std::vector<std::string> tokens = Tokenize(command, ' ');
    if(tokens.at(0) == iterateKey){
      //Command: Iterate solutionSpaceName
      for(typename std::vector<SolutionSpaceBase<Type>*>::iterator it = solSpaces.begin(); 
	  it != solSpaces.end(); ++it){
	if((*it)->name == tokens.at(1)){
	  SolutionOperation<Type>* op = new OperationUpdate<Type>(**it);
	  opList.push_back(op);
	}
      }
    }
    else if(tokens.at(0) == transferKey){
      //Command(field): Transfer solutionVariableName1 solutionSpaceName1 BC#n To 
      //                  solutionVariableName2 solutionSpaceName2 BC#m
      //OR Command(scalar): Transfer solutionVariableName1 solutionSpaceName1 To 
      //                      solutionVariableName2 solutionSpaceName2 
      size_t loc = 0;
      std::string bc1;
      std::string bc2; 
      std::string var1; 
      std::string var2;
      std::string space1;
      std::string space2;
      SolutionSpaceBase<Type>* s1 = NULL;
      SolutionSpaceBase<Type>* s2 = NULL;
      std::string bcKey = "BC#";
      Int transferType = 0;
      loc = command.find(bcKey);
      //if we don't find the bc key, we are transferring a scalar field
      if(loc != std::string::npos){
	bc1 = tokens.at(3);
	bc2 = tokens.at(7);
	var1 = tokens.at(1);
	var2 = tokens.at(5);
	space1 = tokens.at(2);
	space2 = tokens.at(6);
	transferType = 0;
      }
      else{
	bc1 = "";
	bc2 = "";
	var1 = tokens.at(1);
	var2 = tokens.at(4);
	space1 = tokens.at(2);
	space2 = tokens.at(5);
	transferType = 1;
      }

      //look for the first space (source)
      for(typename std::vector<SolutionSpaceBase<Type>*>::iterator it = solSpaces.begin(); 
	  it != solSpaces.end(); ++it){
	if((*it)->name == space1){
	  s1 = *it;
	  if(transferType == 0){
	    //check that solution field exists, must be able to perform this cast
	    SolutionSpace<Type>* space = dynamic_cast<SolutionSpace<Type>*>(s1);
	    SolutionField<Type>& field = space->GetField(var1);
	  }
	  else{
	    //check that scalar field exists
	    ScalarField<Type>& scalar = s1->GetScalarField(var1);
	  }
	  break;
	}
      }
      //look for the second space (destination)
      for(typename std::vector<SolutionSpaceBase<Type>*>::iterator it = solSpaces.begin(); 
	  it != solSpaces.end(); ++it){
	if((*it)->name == space2){
	  s2 = *it;
	  if(transferType == 0){
	    //create solution field, must be able to perform this cast
	    SolutionSpace<Type>* space = dynamic_cast<SolutionSpace<Type>*>(s2);
	    SolutionField<Type>& field = space->GetField(var2);
	  }
	  else{
	    //create scalar field
	    s2->AddScalarField(var2);
	  }
	  break;
	}
      }
      if(s1 == NULL || s2 == NULL){
	Abort << "WARNING: Transfer solution spaces not found.  Found space " + space1 
	  + " and space " + space2 + " instead. These are not registered. Operation - " + command;
	return (-1);
      }

      SolutionOperation<Type>* op = NULL;
      if(transferType == 0){
	op = new OperationTransfer<Type>(*s1, *s2, 
					 dynamic_cast<SolutionSpace<Type>*>(s1)->GetField(var1),
					 dynamic_cast<SolutionSpace<Type>*>(s2)->GetField(var2));
      }
      else{
	op = new OperationTransfer<Type>(*s1, *s2, 
					 s1->GetScalarField(var1), 
					 s2->GetScalarField(var2));
      }
      opList.push_back(op);
      //negative values indicate whole field transfer
      //opList.back()->bcSrc = -1;
      //opList.back()->bcDest = -1;
    }
    else{
      Abort << "WARNING: SolutionOrdering::Finalize() did not recognize operation line " + command;
    }
  }
  if(opList.size() != commands.size()){
    Abort << "WARNING: SolutionOrdering::Finalize() did not find the same number of operations as commands given in param file";
    return (-1);
  }
  return 0;
}

template <class Type>
void SolutionOrdering<Type>::Print()
{
  std::cout << "SOLUTION ORDERING\n";
  std::cout << "===================\n";
  for(std::vector<std::string>::iterator it = commands.begin(); it != commands.end(); ++it){
    std::cout << "\t+ " << *it << "\n";
  }
  std::cout << std::endl;
}

template <class Type>
void SolutionOrdering<Type>::Iterate()
{
  Int i = 0;
  for(typename std::vector<SolutionOperation<Type>*>::iterator it = opList.begin(); 
      it != opList.end(); ++it){
    SolutionOperation<Type>& op = **it;
    std::cout << " *** Applying operation: " << commands[i] << std::endl;
    op.Apply();
    i++;
  }
}


template <class Type> template <class Type2>
SolutionOrdering<Type>& SolutionOrdering<Type>::operator= (const SolutionOrdering<Type2>& orderingToCopy)
{
  for(typename std::vector<std::string>::const_iterator it = orderingToCopy.commands.begin(); 
      it != orderingToCopy.commands.end(); ++it){
    const std::string & command = *it;
    Insert(command);
  }
  return(*this);
}


//function which will parse operation ordering information from the param file
//casename - root of the case for instance cube.0.h5 files have a root of cube
//pathname - path to the case files
template <class Type>
Int SolutionOrdering<Type>::Read(std::string casename, std::string pathname)
{
  std::ifstream fin;
  std::string trash, temp;
  Int err = 0;
  Int err2 = 0;
  char c;

  std::string beginSolutionOrdering = "<<<BEGIN SOLUTION ORDERING>>>";
  std::string endSolutionOrdering = "<<<END SOLUTION ORDERING>>>";
  std::string filename = pathname+casename + ".param";
  
  fin.open(filename.c_str());
  if(fin.is_open()){
    std::cout << "PARAM: Reading solution ordering --> " << filename << std::endl;
    //parse all of the solution spaces from the file
    while(!fin.eof()){
      c = fin.peek();
      if(c == '#' || c == ' ' || c == '\n'){
	getline(fin, trash);
	trash.clear();
	continue;
      }
      else{
	size_t loc;
	getline(fin, temp);
	loc = temp.find(beginSolutionOrdering);
	if(loc != std::string::npos){
	  std::streampos beginOrderpos = fin.tellg();
	  std::streampos endOrderpos = fin.tellg();
	  //we have the beginning of the solution ordering definition, now find the end
	  while(!fin.eof()){
	    c = fin.peek();
	    if(c == '#' || c == ' ' || c == '\n'){
	      getline(fin, trash);
	      trash.clear();
	      continue;
	    }
	    else{
	      std::streampos preEndOrderpos = fin.tellg();
	      getline(fin, temp);
	      loc = temp.find(endSolutionOrdering);
	      if(loc != std::string::npos){
		endOrderpos = fin.tellg();
		ReadOrderingSegment(fin, beginOrderpos, preEndOrderpos);
		//seek back to where we were to continue, routine readspace may modify stream
		fin.seekg(endOrderpos);
		break;
	      }
	    }
	  }
	  if(beginOrderpos == endOrderpos){
	    Abort << "WARNING: Solution ordering end delimiter not found";
	    return(1);
	  }
	}
	temp.clear();
      }
    }
    fin.close();
  }
  else{
    Abort << "PARAM READFILE: Cannot open param file --> " + filename;
    return (1);
  }
  
  if(this->size() == 0){
    Abort << "PARAM READFILE: Did not find any defined solution spaces\n\tThis is defined by <<<BEGIN SOLUTION ORDERING>>> and <<<END SOLUTION ORDERING>>> delimiters";
    return (1);
  }

  if(err2 == 0){
    this->Print();
  }

  return err2;
}

template <class Type>
Int SolutionOrdering<Type>::ReadOrderingSegment(std::ifstream& fin, std::streampos locBegin, std::streampos locEnd)
{
  std::string trash, temp;
  Int err = 0;
  Int err2 = 0;
  char c;
  
  fin.clear();
  fin.seekg(locBegin);
  std::cout << "SOLUTION ORDERING: Reading solution ordering from parameter file" << std::endl;
  while(fin.tellg() < locEnd){
    c = fin.peek();
    if(c == '#' || c == ' ' || c == '\n'){
      getline(fin, trash);
      trash.clear();
      continue;
    }
    else{
      getline(fin, temp);
      err = Insert(temp);
      if(err != 0){
	Abort << "SOLUTION ORDERING: could not parse line -- " + temp +
	  "SOLUTION ORDERING: check that correct keyword is used";
      }
      err2 += err;
      temp.clear();
    }
  }
  return (err2);
};
