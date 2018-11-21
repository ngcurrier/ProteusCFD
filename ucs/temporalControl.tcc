template <class Type>
TemporalControl<Type>::TemporalControl()
{ 
  //SET THE REAL VALUED LIST PARAMETERS
  paramListReal.push_back(Parameter<Type>("timeStep", &dt, 0.0, 0.0, 10000.0));
  //SET THE INTEGER VALUED LIST PARAMETERS
  paramListInt.push_back(Parameter<Int>("numTimeSteps", &nSteps, 100, 0, 999999));
  paramListInt.push_back(Parameter<Int>("newtonIterations", &newtonIter, 1, 0, 999));
  //SET THE BOOLEAN VALUED LIST PARAMETERS
  //paramListBool.push_back(Parameter<Bool>("example", &reference, default_value(t/f));
    
}

template <class Type>
TemporalControl<Type>::~TemporalControl()
{ }
 
template <class Type>
void TemporalControl<Type>::Print()
{
  std::cout << "/***************************************************/\n";
  std::cout << "  TEMPORAL CONTROL: Runtime parameters from file" << "\n";
  std::cout << "/***************************************************/\n\n";
  std::cout << "Number of timesteps: " << nSteps << "\n";
  std::cout << "Number of Newton iterations: " << newtonIter << "\n\n";
}

//function which will parse temporal control things from the param file
//casename - root of the case for instance cube.0.h5 files have a root of cube
//pathname - path to the case files
template <class Type>
Int TemporalControl<Type>::Read(std::string casename, std::string pathname)
{
  std::ifstream fin;
  std::string trash, temp;
  Int err = 0;
  Int err2 = 0;
  char c;

  std::string beginTemporalControl = "<<<BEGIN TEMPORAL CONTROL>>>";
  std::string endTemporalControl = "<<<END TEMPORAL CONTROL>>>";
  std::string filename = pathname+casename + ".param";

  fin.open(filename.c_str());

  if(fin.is_open()){
    std::cout << "PARAM: Reading temporal control file --> " << filename << std::endl;
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
	loc = temp.find(beginTemporalControl);
	if(loc != std::string::npos){
	  std::streampos beginTemppos = fin.tellg();
	  std::streampos endTemppos = fin.tellg();
	  //we have the beginning of the temporal control section, now find the end
	  while(!fin.eof()){
	    c = fin.peek();
	    if(c == '#' || c == ' ' || c == '\n'){
	      getline(fin, trash);
	      trash.clear();
	      continue;
	    }
	    else{
	      std::streampos preEndTemppos = fin.tellg();
	      getline(fin, temp);
	      loc = temp.find(endTemporalControl);
	      if(loc != std::string::npos){
		endTemppos = fin.tellg();
		ReadTemporalControlSegment(fin, beginTemppos, preEndTemppos);
		//seek back to where we were to continue, routine readspace may modify stream
		fin.seekg(endTemppos);
		break;
	      }
	    }
	  }
	  if(beginTemppos == endTemppos){
	    std::cerr << "WARNING: Temporal control end delimiter not found";
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
  return err2;

}


template <class Type>
void TemporalControl<Type>::ReadTemporalControlSegment(std::ifstream& fin, std::streampos locBegin, std::streampos locEnd)
{
  std::string trash, temp;
  Int err = 0;
  Int err2 = 0;
  char c;
  
  std::cout << "TEMPORAL CONTROL: Reading temporal parameters" << std::endl;
  
  fin.clear();
  fin.seekg(locBegin);
  while(fin.tellg() < locEnd){
    c = fin.peek();
    if(c == '#' || c == ' ' || c == '\n'){
      getline(fin, trash);
      trash.clear();
      continue;
    }
    else{
      getline(fin, temp);
      err = ParseLine(temp);
      if(err != 0){
	std::cerr << "TEMPORAL CONTROL: could not parse line -- " + temp +  
	  "\nTEMPORAL CONTROL: check that correct keyword is used";
	abort();
      }
      err2 += err;
      temp.clear();
    } 
  }
  Print();
}

template <class Type>
Int TemporalControl<Type>::ParseLine(std::string& line)
{
  if(line.length() == 0){
    return (0);
  }
  for(std::vector<Parameter<Bool> >::iterator it = paramListBool.begin();
      it != paramListBool.end(); ++it){
    Parameter<Bool>& p = *it;
    if(p.ParseLine(line)) return 0;
  }
  for(std::vector<Parameter<Int> >::iterator it = paramListInt.begin();
      it != paramListInt.end(); ++it){
    Parameter<Int>& p = *it;
    if(p.ParseLine(line)) return 0;
  }
  for(typename std::vector<Parameter<Type> >::iterator it = paramListReal.begin();
      it != paramListReal.end(); ++it){
    Parameter<Type>& p = *it;
    if(p.ParseLine(line)) return 0;
  }

  return 1;
}

