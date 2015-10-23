#include "dataInfo.h"
#include "h5layer.h"
#include "mesh.h"
#include "exceptions.h"

//TODO: if a field is declared everywhere and with a temporal state,
//we only need the most recent NP1 storage to have all the boundary data included, fix this -- expensive

template <class Type>
SolutionField<Type>::SolutionField(Mesh<Type> & mesh, DataInfo dataInfo, Int stateType, Int varLocation):
  mesh(mesh), dataInfo(dataInfo), mystate(stateType), varLoc(varLocation)
{
  Int mnnode = mesh.GetNumNodes();
  Int mgnode = mesh.GetNumParallelNodes();
  Int mnbnode = mesh.GetNumBoundaryNodes();
  Int mnallnode = mnnode + mgnode + mnbnode;
  Int narrays = 0;
  ndof = dataInfo.GetNdof();
  if(varLoc == FIELDS::VAR_EVERYWHERE){
    nentities = mnallnode;
  }
  else if(varLoc == FIELDS::VAR_INTERIOR){
    nentities = mnnode + mgnode;
  }
  if(mystate == FIELDS::STATE_NONE){
    narrays = 1;
  }
  else if(mystate == FIELDS::STATE_TIME){
    narrays = 3;
  }
  else{
    Abort << "State passed to solution field named " + dataInfo.GetName() + " not valid\n";
  }

  data = new Type[nentities*ndof*narrays];

#ifdef _OPENMP
  mutexes = new pthread_mutex_t[nentities];
  for(Int i = 0; i < nentities; i++){
    pthread_mutex_init(&mutexes[i], NULL);
  }
#endif
}

template <class Type>
SolutionField<Type>::SolutionField(Mesh<Type>& mesh, hid_t fileId, std::string path):
  mesh(mesh)
{
  Int mnnode = mesh.GetNumNodes();
  Int mgnode = mesh.GetNumParallelNodes();
  Int mnbnode = mesh.GetNumBoundaryNodes();
  Int mnallnode = mnnode + mgnode + mnbnode;
  Int nentities = 0;
  size_t loc = path.rfind("/");
  if(loc == std::string::npos){
    Abort << "SolutionField::SolutionField() " + path + " not a valid HDF file path";
  }
  std::string directory = path.substr(0, loc);
  std::string dataName = path.substr(loc+1, std::string::npos);
  dataInfo.SetFromHDF(fileId, directory, dataName);
  ndof = dataInfo.GetNdof();
  mystate = FIELDS::STATE_NONE;
  Int narrays = 1;

  //Read attribute to determine where data was written from (interior, boundaries, etc.)
  std::vector<std::string> location;
  HDF_ReadStringAttribute(fileId, directory, dataName, "variable_location", location);
  if(location.front().find("Interior and Boundaries") != std::string::npos){
    varLoc = FIELDS::VAR_EVERYWHERE;
    nentities = mnallnode;
  }
  else if(location.front().find("Interior Only") != std::string::npos){
    varLoc = FIELDS::VAR_INTERIOR;
    nentities = mnnode + mgnode;
  }
  else{
    Abort << "WARNING: In SolutionField::SolutionField() variable location undefined from " 
	  + path + " in HDF!" + "\n\tfound " + location.front() + "...";
    nentities = 0;
    return;
  }
  //we simply allocate enough space for the data, reading has to be done elsewhere
  Int size = nentities*ndof*narrays;
  data = new Type[size];

#ifdef _OPENMP
  mutexes = new pthread_mutex_t[nentities];
  for(Int i = 0; i < nentities; i++){
    pthread_mutex_init(&mutexes[i], NULL);
  }
#endif

  return;
}

template <class Type>
SolutionField<Type>::~SolutionField()
{
  delete [] data;

#ifdef _OPENMP
  delete [] mutexes;
#endif
}

template <class Type>
Type* SolutionField<Type>::GetData(Int state)
{
  if(mystate == FIELDS::STATE_NONE){
    return(data);
  }
  else if(mystate == FIELDS::STATE_TIME && state == FIELDS::STATE_NP1){ 
    return(data);
  }
  else if(mystate == FIELDS::STATE_TIME && state == FIELDS::STATE_N){
    return(&data[1*ndof*nentities]);
  }
  else if(mystate == FIELDS::STATE_TIME && state == FIELDS::STATE_NM1){
    return(&data[2*ndof*nentities]);
  }
  Abort << "WARNING: SolutionField::GetData() falling through if nest for field named " 
	+ dataInfo.GetName() + " ... returning NULL";
  return NULL;
}

template <class Type>
const Type* SolutionField<Type>::GetData(Int state) const
{
  if(mystate == FIELDS::STATE_NONE){
    return(data);
  }
  else if(mystate == FIELDS::STATE_TIME && state == FIELDS::STATE_NP1){ 
    return(data);
  }
  else if(mystate == FIELDS::STATE_TIME && state == FIELDS::STATE_N){
    return(&data[1*ndof*nentities]);
  }
  else if(mystate == FIELDS::STATE_TIME && state == FIELDS::STATE_NM1){
    return(&data[2*ndof*nentities]);
  }
  Abort << "WARNING: SolutionField::GetData() falling through if nest for field named " 
	+ dataInfo.GetName() + " ... returning NULL";
  return NULL;
}

template <class Type>
Int SolutionField<Type>::GetNdof() const
{
  return dataInfo.GetNdof();
}

template <class Type>
std::string SolutionField<Type>::GetDofName(Int dof) const
{
  return dataInfo.GetDofName(dof);
}

template <class Type>
Bool SolutionField<Type>::DofIsVector(Int dof) const
{
  return dataInfo.DofIsVector(dof);
}

template <class Type>
Bool SolutionField<Type>::DofIsScalar(Int dof) const
{
  return dataInfo.DofIsScalar(dof);
}

template <class Type>
void SolutionField<Type>::EvolveInTime()
{
  if(mystate == FIELDS::STATE_NONE){
    return;
  }
  else{
    memcpy(GetData(FIELDS::STATE_NM1), GetData(FIELDS::STATE_N), sizeof(Type)*nentities*ndof);
    memcpy(GetData(FIELDS::STATE_N), GetData(FIELDS::STATE_NP1), sizeof(Type)*nentities*ndof);
  }
}
template <class Type>
std::string SolutionField<Type>::GetName() const
{
  return dataInfo.GetName();
}

template <class Type>
Bool SolutionField<Type>::IsNamed(std::string testname) const
{
  if(testname == dataInfo.GetName()){
    return true;
  }
  return false;
}

template <class Type>
Bool SolutionField<Type>::IsTemporal() const
{
  if(mystate == FIELDS::STATE_TIME) return true;
  return false;
}

template <class Type>
void SolutionField<Type>::Reorder(Int* newOrder)
{
  Int i, j;
  Type* temp = new Type[nentities*ndof];
  if(mystate == FIELDS::STATE_NONE){
    Type* ptr = data[0];
    for(i = 0; i < nentities; i++){
      Int newloc = newOrder[i];
      memcpy(&temp[newloc*ndof + 0], &ptr[i*ndof + 0], sizeof(Type)*ndof);
    }
    memcpy(ptr, temp, sizeof(Type)*ndof*nentities);
    //now assume we've also reordered on all other processes, call parallel sync
    mesh.p->UpdateGeneralVectors(ptr, ndof);
  }
  else if(mystate == FIELDS::STATE_TIME){
    Int numStates = 3;
    for(j = 0; j < numStates; j++){
      Type* ptr = data[j*nentities*ndof];
      for(i = 0; i < nentities; i++){
	Int newloc = newOrder[i];
	memcpy(&temp[newloc*ndof + 0], &ptr[i*ndof + 0], sizeof(Type)*ndof);
      }
      memcpy(ptr, temp, sizeof(Type)*ndof*nentities);
      //now assume we've also reordered on all other processes, call parallel sync
      mesh.p->UpdateGeneralVectors(ptr, ndof);
    }
  }

  return;
}

template <class Type>
void SolutionField<Type>::Fill(Type value)
{
  for(Int i = 0; i < ndof*nentities; i++){
    data[i] = value;
  }
}


template <class Type>
void SolutionField<Type>::WriteH5(std::string filename, std::string directory, Int states)
{
  hid_t file_id = -1;
  Int narrays = 1;
  if(mystate == FIELDS::STATE_TIME && states == FIELDS::STATE_TIME){
    narrays = 3;
  }
  Int nvalues = ndof*nentities*narrays;
  file_id = HDF_OpenFile(filename, 1);
  if(file_id < 0){
    Abort << "SolutionField::WriteH5() could not open file -- " + filename;
    return;
  }
  HDF_WriteArray(file_id, directory, dataInfo.GetName(), data, nvalues);
  dataInfo.WriteHDFAttribute(file_id, directory);
  //write the variable location info to the attributes
  std::vector<std::string> location;
  if(varLoc == FIELDS::VAR_EVERYWHERE){
    location.push_back("Interior and Boundaries");
  }
  else if(varLoc == FIELDS::VAR_INTERIOR){
    location.push_back("Interior Only");
  }
  else{
    Abort << "WARNING: In SolutionField::WriteH5() variable location undefined!";
  }
  HDF_WriteStringAttribute(file_id, directory, dataInfo.GetName(), "variable_location", location);
  HDF_CloseFile(file_id);
}

template <class Type>
void SolutionField<Type>::ReadH5(std::string filename, std::string directory, Int states)
{
  hid_t file_id = -1;
  Int narrays = 1;
  if(mystate == FIELDS::STATE_TIME && states == FIELDS::STATE_TIME){
    narrays = 3;
  }
  Int nvalues = ndof*nentities*narrays;
  file_id = HDF_OpenFile(filename, 0);
  if(file_id < 0){
    Abort << "SolutionField::ReadH5() could not open file -- " + filename;
    return;
  }
  HDF_ReadArray(file_id, directory, dataInfo.GetName(), &data, &nvalues);
  HDF_CloseFile(file_id);
}

template <class Type>
PObj<Type>* SolutionField<Type>::GetPObj()
{
  return mesh.p;
}

template <class Type> template <class Type2>
SolutionField<Type> & SolutionField<Type>::operator= (const SolutionField<Type2>& fieldToCopy)
{
  //copy over the solution field, this is assuming the data descriptors, mesh object, etc.
  //have already been created in the destination field, if you need a deep copy of all
  //these things then this routine has to be re-done, this just copies values in the field
  
  //check ndofs and temporal status match
  if((GetNdof() != fieldToCopy.GetNdof()) || IsTemporal() != fieldToCopy.IsTemporal()){
    Abort << "In SolutionField::operator= NDOFs don't match, the probably means the scope of, this function has a need to expand, please re-write";
  }
 
  const Type2* dataToCopy = fieldToCopy.GetData(FIELDS::STATE_NP1);
  Type* dataDest = GetData(FIELDS::STATE_NP1);
  Int ndof = GetNdof();
  Int nvalues = 0;

  if(IsTemporal()){
    nvalues = ndof*3*nentities;
  }
  else{
    nvalues = ndof*nentities;
  }

  for(int i = 0; i < nvalues; i++){
    dataDest[i] = dataToCopy[i];
  }
  return(*this);
}

template <class Type>
Type SolutionField<Type>::GetMax() const
{
  MPI_Datatype mpit;
  Int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Int np = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  Type* max = (Type*)alloca(sizeof(Type)*np);

  max[rank] = CAbs(data[0]);
  //get mpi datatype to send
  mpit = MPI_GetType(max[0]);

  //reduce across all processes, MPI_MAX is not valid for complex types
  //we have to hack around it
  for(Int i = 0; i < nentities*ndof; i++){
    max[rank] = MAX(CAbs(data[i]), max[rank]);
  }
  MPI_Allgather(MPI_IN_PLACE, 1, mpit, max, 1, mpit, MPI_COMM_WORLD);

  for(Int i = 0; i < np; i++){
    max[0] = MAX(max[i], max[0]);
  }

  return max[0];
}

template <class Type>
void SolutionField<Type>::Normalize() 
{
  Type max = GetMax();
  if(max == (Type)(0.0)){
    Fill(Type(0.0));
  }
  else{
    for(Int i = 0; i < nentities*ndof; i++){
      data[i] /= max;
    }
  }
}
