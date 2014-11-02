template <class Type>
CompositeBody<Type>::CompositeBody()
{
  name = "default";
  cl = cd = cm = 0.0;
  nsurfs = 0;
  list = NULL;
  momentPt = new Type[3];
  momentAxis = new Type[3];
  surfArea = new Type[3];
  //pick origin as default location to compute moments
  momentPt[0] = momentPt[1] = momentPt[2] = 0.0;
  //pick z-axis as the default moment axis to project to
  momentAxis[2] = 1.0;
  momentAxis[0] = momentAxis[1] = 0.0;
  surfArea[0] = surfArea[1] = surfArea[2] = 0.0;

  forces = new Type[3];
  vforces = new Type[3];
  moments = new Type[3];
  vmoments = new Type[3];

  return;
}

template <class Type>
CompositeBody<Type>::~CompositeBody()
{
  delete [] list;
  delete [] momentPt;
  delete [] momentAxis;
  delete [] surfArea;

  delete [] forces;
  delete [] vforces;
  delete [] moments;
  delete [] vmoments;

  return;
}

template <class Type>
void CompositeBody<Type>::Init(Int* list, Int nsurfs)
{
  Int i;
  this->nsurfs = nsurfs;
  this->list = new Int[nsurfs];
  for(i = 0; i < nsurfs; i++){
    this->list[i] = list[i];
  }
  return;
}

template <class Type>
void CompositeBody<Type>::Print(Int id)
{
  Int i;
  if(id >= 0){
    std::cout << "COMPOSITE: Composite body " << id << " made of surfaces - [";
  }
  else{
    std::cout << "COMPOSITE: Composite body made of surfaces - [";
  }
  for(i = 0; i < nsurfs; i++){
    std::cout << " " << list[i];
  }
  std::cout << " ]";
  std::cout << " Moment CG: " << momentPt[0] << " " << momentPt[1] 
	    << " " << momentPt[2];
  std::cout << " Moment Axis: " << momentAxis[0] << " " << momentAxis[1] 
	    << " " << momentAxis[2] << std::endl;
  return;
}

template <class Type>
Int BodiesReadFile(CompositeBody<Type>* bodies, std::string casename)
{ 
  std::ifstream fin;
  std::string line;
  char c;
  size_t loc;

  std::string filename = casename + ".bc";
  fin.open(filename.c_str());
  std::string body_str = "body";

  if(fin.is_open()){
    std::cout << "COMPOSITE: Reading boundary conditions file --> " << filename << std::endl;

    while(!fin.eof()){
      c = fin.peek();
      getline(fin, line);
      if(c == '#' || c == ' ' || c == '\n'){
	line.clear();
	continue;
      }
      else{
	loc = line.find(body_str);
	if(loc != std::string::npos){
	  Int bodyId = ParseBody(bodies, line);
	  if(bodyId >= 0){
	    bodies[bodyId].Print(bodyId);
	  }
	  else{
	    std::cerr << "COMPOSITE: could not parse line -- " << line << " in file " 
		 << filename << std::endl;
	    std::cerr << "COMPOSITE: check that correct keyword is used" << std::endl;
	    return (1);
	  }
	}
	line.clear();
      }
    }
    
    fin.close();
  }
  else{
    std::cerr << "COMPOSITE READFILE: Cannot open bc file --> " << filename << std::endl;
    return (1);
  }
  
  return(0);
}

template <class Type>
Int ParseBody(CompositeBody<Type>* bodies, std::string& line)
{
  size_t loc;
  std::string subline;
  std::stringstream ss (std::stringstream::in | std::stringstream::out);

  Int* listSurfs = NULL;
  Int nsurfs;
  Int bodyId;

  loc = line.find('#');
  loc += 1;
  line = line.substr(loc);
  ss << line;
  ss >> bodyId;
  if(!GetStringBetween("[","]",line,subline)){
    nsurfs = StripCSV(subline, &listSurfs); 
    //initialize the composite body
    bodies[bodyId].Init(listSurfs, nsurfs);
    delete [] listSurfs;
  }
  else{
    std::cerr << "WARNING: Body " << bodyId << " found with no list of surfaces" << std::endl;
    return (-1);
  }
  
  return bodyId;
}

template <class Type>
Bool CompositeBody<Type>::SurfIsPart(Int surfId)
{
  Int i;
  for(i = 0; i < nsurfs; i++){
    if(surfId == list[i]){
      return true;
    }
  }
  return false;
}


