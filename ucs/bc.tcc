#include "driver.h"
#include "strings_util.h"
#include "macros.h"
#include "composite.h"
#include "bc_defines.h"
#include "mesh.h"
#include "param.h"
#include "eqnset.h"
#include "solutionField.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <set>
#include <sstream>

template <class Type>
Int BoundaryConditions<Type>::GetBCType(Int factag)
{
  return GetBCObj(factag)->GetBCType();;
}

template <class Type>
BCObj<Type>* BoundaryConditions<Type>::GetBCObj(Int factag)
{
  BCObj<Type>* bco = bc_objs[factag];
  if(bco != NULL){
    return bco;
  }
  std::stringstream ss;
  ss << "GetBCObj() got a NULL pointer -- this shouldn't happen -- dying\n";
  ss << "Request made for factag " << factag;
  Abort << ss.str();
}

template <class Type>
BoundaryConditions<Type>::BoundaryConditions()
{
  num_bcs = 0;
}

template <class Type>
BoundaryConditions<Type>::~BoundaryConditions()
{
}

template <class Type>
void BoundaryConditions<Type>::AllocateBCObjects()
{
  std::cout << "BC: Allocating boundary condition objects" << std::endl;
  for(Int i = 0; i <= largest_bc_id; i++){
    bc_objs[i] = new BCObj<Type>;
    bc_objs[i]->SetBCFactag(i);
  }
  //setup the parallel object for internal use
  bc_objs[0]->SetName("PARALLEL_COMM_BOUNDARY");
  bc_objs[0]->SetBCType(Proteus_ParallelBoundary);
  std::cout << "BC: Finished allocating boundary condition objects" << std::endl;
}

template <class Type>
Int BoundaryConditions<Type>::InitInternals(EqnSet<Type>* eqnset)
{
  for(Int i = 0; i <= num_bcs; i++){
    BCObj<Type>* bcobj = GetBCObj(i);
    bcobj->neqn = eqnset->neqn;
    bcobj->nvars = bcobj->neqn + eqnset->nauxvars;
    //if we haven't already set the Qref state from file
    //go ahead and set it to point at the eqnset version
    if(bcobj->Qref == NULL){
      bcobj->Qref = eqnset->Qinf;
    }
    else{
      if(bcobj->QrefFromFile == true){
	eqnset->ComputeAuxiliaryVariables(bcobj->Qref);	
      }
      else{
	Abort << "BC: WARNING -- BCs not set from file and pointer not NULL";
	return (-1);
      }
    }
  }

  return 0;
}

template <class Type>
Int BoundaryConditions<Type>::ReadFile(std::string casename, Int maximumFactagExpected)
{
  largest_bc_id = maximumFactagExpected;
  std::ifstream fin;
  std::string temp;
  Int err = 0;
  Int err2 = 0;
  char c;
  Int temp_num_bcs = 0;

  std::string filename = casename + ".bc";
  fin.open(filename.c_str());

  if(fin.is_open()){
    std::cout << "BC: Reading boundary conditions file --> " << filename << std::endl;

    //count bcs in file
    CountBCsInFile(fin);
    //count and allocate the composite bodies in the file
    CountBodiesInFile(fin);
    //allocate all bc objects
    AllocateBCObjects();

    //set a temp value to check that we got all the bcs expected
    temp_num_bcs = num_bcs;
    //zero num_bcs to check later
    num_bcs = 0;

    while(!fin.eof()){
      c = fin.peek();
      getline(fin, temp);
      if(c == '#' || c == ' ' || c == '\n'){
	temp.clear();
	continue;
      }
      else{
	err = BCParseLine(temp);
	if(err != 0){
	  std::cerr << "BC: could not parse line -- " << temp << " in file " 
	       << filename << std::endl;
	  std::cerr << "BC: check that correct keyword is used" << std::endl;
	  err2 = err;
	}
	
	temp.clear();
      }
    }
    
    fin.close();
  }
  else{
    std::cerr << "BC READFILE: Cannot open bc file --> " << filename << std::endl;
    return (1);
  }
  
  //check that we actually read what we expected
  if(num_bcs != temp_num_bcs){
    std::cerr << "BC: number of BCs anticipated does not match number read!!" << std::endl;
    std::cerr << "BC: number of BCs anticipated " << temp_num_bcs << std::endl;
    std::cerr << "BC: number of BCs read " << num_bcs << std::endl;
    err = 1;
  }

  
  if(err2 == 0 && err == 0){
    PrintBCs();
  }

  return(err2||err);
}

template <class Type>
void BoundaryConditions<Type>::CountBCsInFile(std::ifstream& fin)
{
  std::string surf_str = "surface";
  size_t loc;
  char c;
  std::string line;
  std::string subline;
  Int surf_id;
  std::stringstream ss (std::stringstream::in | std::stringstream::out);
  num_bcs = 0;

  while(!fin.eof()){
    c = fin.peek();
    getline(fin, line);
    if(c == '#' || c == ' ' || c == '\n'){
      line.clear();
      continue;
    }
    else{
      if(line.length() == 0){ 
	line.clear();
	continue;
      }
      loc = line.find(surf_str);
      if(loc != std::string::npos){
	loc = line.find('#');
	loc += 1;
	subline = line.substr(loc);
	ss << subline;
	ss >> surf_id;
	//if we found a surface flag and a number assume the bc is complete
	num_bcs++;
	ss.clear();
	//empty std::stringstream
	ss.str("");
	line.clear();
	subline.clear();
	continue;
      }
      else{
	line.clear();
	continue;
      }
    }
  }

  std::cout << "BC: Surfaces found in file: " << num_bcs << std::endl;
  std::cout << "BC: Largest surface ID expected: " << largest_bc_id << std::endl;

  if(num_bcs == 0){
    Abort << "Did not read any BCs in boundary condition file ";
  }
  
  //clear failbits
  fin.clear();
  //return to file beginning
  fin.seekg(0, std::ios::beg);
}

template <class Type>
void BoundaryConditions<Type>::CountBodiesInFile(std::ifstream& fin)
{
  std::string body_str = "body";
  size_t loc;
  char c;
  std::string line;
  std::string subline;
  Int surf_id;
  std::stringstream ss (std::stringstream::in | std::stringstream::out);
  largest_body_id = -999;
  num_bodies = 0;

  while(!fin.eof()){
    c = fin.peek();
    getline(fin, line);
    if(c == '#' || c == ' ' || c == '\n'){
      line.clear();
      continue;
    }
    else{
      if(line.length() == 0){ 
	line.clear();
	continue;
      }
      loc = line.find(body_str);
      if(loc != std::string::npos){
	loc = line.find('#');
	loc += 1;
	subline = line.substr(loc);
	ss << subline;
	ss >> surf_id;
	largest_body_id = MAX(largest_body_id, surf_id);
	//if we found a body flag and a number assume the composite body is complete
	num_bodies++;
	ss.clear();
	//empty std::stringstream
	ss.str("");
	line.clear();
	subline.clear();
	continue;
      }
      else{
	line.clear();
	continue;
      }
    }
  }
  if(num_bodies){
    std::cout << "COMPOSITE: Composite bodies found in file: " << num_bodies << std::endl;
    std::cout << "COMPOSITE: Largest composite body ID found in file: " << largest_body_id << std::endl;
  }
  else{
    std::cout << "COMPOSITE: Composite bodies NOT found in file" << std::endl;
  }

  //clear failbits
  fin.clear();
  //return to file beginning
  fin.seekg(0, std::ios::beg);

  return;
}


template <class Type>
Int BoundaryConditions<Type>::BCParseLine(std::string& line)
{
  size_t loc, loc2;
  Int surf_id;
  std::string subline;
  std::string tempstring;
  std::stringstream ss (std::stringstream::in | std::stringstream::out);

  std::string surf_str = "surface";
  std::string body_str = "body";

  if(line.length() == 0){
    return (0);
  }

  loc = line.find(body_str);
  if(loc != std::string::npos){
    //just return, composite bodies are handled elsewhere
    return(0);
  }
  loc = line.find(surf_str);
  if(loc != std::string::npos){
    loc = line.find('#');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> surf_id;
    //set factag in bcobj for this surface
    BCObj<Type>* bcobj = GetBCObj(surf_id);
    bcobj->SetBCFactag(surf_id);
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    for(Int itype = 0; itype < NUM_BC_TYPES; itype++){
      loc = subline.find(BCs[itype]);
      if(loc != std::string::npos){
	bcobj->SetBCType(itype);
	subline = line.substr(loc);
	SetVarsFromLine(*bcobj, subline);
	break;
      }
      if(itype == NUM_BC_TYPES - 1){
	std::cerr << "BC: type not found" << std::endl;
	std::cerr << "Line read: " << line << std::endl;
	std::cerr << "VALID types: " << std::endl;
	for(Int j = 0; j < NUM_BC_TYPES; j++){
	  std::cerr << "\t" << BCs[j] << std::endl;
	}
	return (1);
      }
    }
    loc = line.find('"');
    loc2 = line.rfind('"');
    if(loc == std::string::npos){
      bcobj->SetName("Default_Face_Name");
    }
    else{
      loc2 --;
      loc2 -= loc;
      loc += 1;
      subline = line.substr(loc, loc2);
      bcobj->SetName(subline);
    }
    ss.clear();
    ss.str("");
    num_bcs++;
    return(0);
  }
  
  return (1);
}


template <class Type>
Int BoundaryConditions<Type>::SetVarsFromLine(BCObj<Type>& bcobj, std::string& line)
{
  std::string subline;
  std::stringstream ss (std::stringstream::in | std::stringstream::out);
  Int flagsfound[NUM_BC_VAR_IDS];
  
  //look for the flags.. we'll parse them later
  for(Int i = 0; i < NUM_BC_VAR_IDS; i++){
    size_t loc = line.find(BCVars[i]);
    if(loc != std::string::npos){
      flagsfound[i] = 1;
    }
    else{
      flagsfound[i] = 0;
    }
  }

  Int ObjId = bcobj.GetBCFactag();
  for(Int i = 0; i < NUM_BC_VAR_IDS; i++){
    size_t loc;
    if(flagsfound[i]){
      switch (i){
      case QREF:
	loc = line.find(BCVars[QREF]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    Int nitems = StripCSV(subline, &bcobj.Qref);
	    std::cout << "BC: BC #" << ObjId << " Qref found ( ";
	    for(Int j = 0; j < nitems; j++){
	      std::cout << bcobj.Qref[j] << " ";
	    }
	    std::cout << ")" << std::endl;
	    bcobj.QrefFromFile = true;
	  }
	  else{
	    std::cerr << "Found key " << BCVars[QREF]
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;
      case FlowDirection:
	loc = line.find(BCVars[FlowDirection]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.flowDirection[0];
	    ss >> bcobj.flowDirection[1];
	    ss >> bcobj.flowDirection[2];
	    //normalize the flow direction so we don't get momentum sources/sinks
	    Normalize(bcobj.flowDirection, bcobj.flowDirection);
	    std::cout << "BC: BC #" << ObjId << " flow direction found ( " 
		      <<  bcobj.flowDirection[0] << " " 
		      <<  bcobj.flowDirection[1] << " "  
		      <<  bcobj.flowDirection[2] << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	  }
	  else{
	    std::cerr << "Found key " << BCVars[FlowDirection] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	
	break;
      case Velocity:
	loc = line.find(BCVars[Velocity]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.velocity;
	    std::cout << "BC: BC #" << ObjId << " flow velocity magnitude found ( "
		      << bcobj.velocity << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	  }
	  else{
	    std::cerr << "Found key " << BCVars[Velocity] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;
      case BackPressure:
	loc = line.find(BCVars[BackPressure]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.backPressure;
	    std::cout << "BC: BC #" << ObjId << " back pressure found ( "
		      << bcobj.backPressure << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	  }
	  else{
	    std::cerr << "Found key " << BCVars[Velocity] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;
      case Density:
	loc = line.find(BCVars[Density]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.density;
	    std::cout << "BC: BC #" << ObjId << " density found ( "
		      << bcobj.density << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	  }
	  else{
	    std::cerr << "Found key " << BCVars[Density] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;
      case Axis:
	loc = line.find(BCVars[Axis]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.rotationAxis[0];
	    ss >> bcobj.rotationAxis[1];
	    ss >> bcobj.rotationAxis[2];
	    //normalize the rotation axis so we don't get momentum sources/sinks
	    Normalize(bcobj.rotationAxis, bcobj.rotationAxis);
	    std::cout << "BC: BC #" << ObjId << " rotation axis found ( " 
		      <<  bcobj.rotationAxis[0] << " " 
		      <<  bcobj.rotationAxis[1] << " "  
		      <<  bcobj.rotationAxis[2] << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	    bcobj.movement = 2;
	  }
	  else{
	    std::cerr << "Found key " << BCVars[Axis] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;
      case Point:
	loc = line.find(BCVars[Point]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.rotationPoint[0];
	    ss >> bcobj.rotationPoint[1];
	    ss >> bcobj.rotationPoint[2];
	    std::cout << "BC: BC #" << ObjId << " rotation point found ( " 
		      <<  bcobj.rotationPoint[0] << " " 
		      <<  bcobj.rotationPoint[1] << " "  
		      <<  bcobj.rotationPoint[2] << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	    bcobj.movement = 2;
	  }
	  else{
	    std::cerr << "Found key " << BCVars[Point] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;     
      case Omega:
	loc = line.find(BCVars[Omega]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.omega;
	    std::cout << "BC: BC #" << ObjId << " wall rotation speed found ( "
		      << bcobj.omega << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	    bcobj.movement = 2;
	  }
	  else{
	    std::cerr << "Found key " << BCVars[Omega] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;
      case SlipSpeed:
	loc = line.find(BCVars[SlipSpeed]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.slipSpeed;
	    std::cout << "BC: BC #" << ObjId << " wall velocity found ( "
		      << bcobj.slipSpeed << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	    bcobj.movement = 1;
	  }
	  else{
	    std::cerr << "Found key " << BCVars[SlipSpeed] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;
      case SlipDirection:
	loc = line.find(BCVars[SlipDirection]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.slipDirection[0];
	    ss >> bcobj.slipDirection[1];
	    ss >> bcobj.slipDirection[2];
	    //normalize the slip direction so we don't get momentum sources/sinks
	    Normalize(bcobj.slipDirection, bcobj.slipDirection);
	    std::cout << "BC: BC #" << ObjId << " wall slip direction found ( " 
		      <<  bcobj.slipDirection[0] << " " 
		      <<  bcobj.slipDirection[1] << " "  
		      <<  bcobj.slipDirection[2] << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	    bcobj.movement = 1;
	  }
	  else{
	    std::cerr << "Found key " << BCVars[SlipDirection] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;     
      case MassFractions:
	loc = line.find(BCVars[MassFractions]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    bcobj.nspecies = StripCSV(subline, &bcobj.massFractions);
	    //make sure that the mass fractions add up specifically to unity
	    //this will modify the each value slightly but is necessary for 
	    //consistency, normalize them
	    Type sum = 0.0;
	    for(Int j = 0; j < bcobj.nspecies; j++){
	      sum += bcobj.massFractions[j];
	    }
	    for(Int j = 0; j < bcobj.nspecies; j++){
	      bcobj.massFractions[j] = bcobj.massFractions[j]/sum;
	    }
	    std::cout << "BC: BC #" << ObjId << " " << bcobj.nspecies 
		      << " mass fractions found:" << std::endl;
	    std::cout << " ---------------------------" << std::endl;
	    for(Int j = 0; j < bcobj.nspecies; j++){
	      std::cout << "\tSpecies # " << j << ": " << bcobj.massFractions[j] << std::endl;
	    }
	  }
	  else{
	    std::cerr << "Found key " << BCVars[MassFractions] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;     	
      case Moving:
	loc = line.find(BCVars[Moving]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.movingBC;
	    std::cout << "BC: BC #" << ObjId << " movement flag found ( "
		      << bcobj.movingBC << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	  }
	  else{
	    std::cerr << "Found key " << BCVars[Moving] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;
      case BleedSteps:
	loc = line.find(BCVars[BleedSteps]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.bleedSteps;
	    std::cout << "BC: BC #" << ObjId << " bleed steps flag found ( "
		      << bcobj.bleedSteps << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	  }
	  else{
	    std::cerr << "Found key " << BCVars[BleedSteps] 
		      << " but no [...] for definition" << std::endl;
	  }
	}
	break;
      case Twall:
	loc = line.find(BCVars[Twall]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]",line,subline,loc)){
	    ss << subline;
	    ss >> bcobj.twall;
	    std::cout << "BC: BC #" << ObjId << " twall flag found ( "
		      << bcobj.twall << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	  }
	  else{
	    std::cerr << "Found key " << BCVars[Twall] 
		      << " but no [...] for definition" << std::endl;
	  }

	}
      case Flux:
	loc = line.find(BCVars[Flux]);
	if(loc != std::string::npos){
	  if(!GetStringBetween("[","]", line, subline, loc)){
	    ss << subline;
	    ss >> bcobj.flux;
	    std::cout << "BC: BC #" << ObjId << " heatFlux flag found ( "
		      << bcobj.flux << " )" << std::endl;
	    ss.clear();
	    ss.str("");
	  }
	}
	else{
	  std::cerr << "Found key " << BCVars[Flux]
		    << " but no [...] for defintion" << std::endl;
	}
      default:
	break;
      }
      
    }
  }
  return(0);
}

template <class Type>
void BoundaryConditions<Type>::PrintBCs()
{
  std::cout << "BC: Runtime boundary conditions from file" << std::endl;
  std::cout << "=========================================" << std::endl;
  
  //NOTE: zeroth face not used.... defaulted to volume factag
  for(Int i = 1; i <= this->largest_bc_id; i++){
    if(this->GetBCObj(i)->GetBCType() != Proteus_NULL){
      std::cout << "\tSurface " <<  i << "\t" << BCs[this->GetBCObj(i)->GetBCType()] << "\t : " << this->GetBCObj(i)->GetName() << std::endl;
    }
  }
  std::cout << std::endl;
}

template <class Type> template <class Type2>
Bool BoundaryConditions<Type>::IsNodeOnBC(Int nodeid, Mesh<Type2>* m, Int bcType)
{
  //Check surface elements surround point
  for(Int* indx = m->SelspBegin(nodeid); indx != m->SelspEnd(nodeid); ++indx){
    Int selemid = *indx;;
    Int factag = m->elementList[selemid]->GetFactag();
    Int bcTypeLocal = GetBCType(factag); 
    if(bcType == bcTypeLocal){
      return true;
    }
  }
  return false;
}

template <class Type>
void BC_Kernel(B_KERNEL_ARGS)
{
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;

  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* QL = &space->q[left_cv*nvars];    //boundary point
  Type* QR = &space->q[right_cv*nvars];   //ghost point
  BoundaryConditions<Real>* bc = space->bc;
  Int bcType = bc->GetBCType(factag); 
  BCObj<Real>* bcobj = bc->GetBCObj(factag);

  //TODO: how to handle this with the driver module
  //until then hack it out
  *ptrL = NULL;
  *size = 0;

  Type* Qref = (Type*)alloca(sizeof(Type)*nvars);
  bcobj->GetQref(Qref);
  CalculateBoundaryVariables(eqnset, m, space, QL, QR, Qref, avec, bcType, eid, bcobj, vdotn, velw);
}

template <class Type>
void Bkernel_BC_Jac_Modify(B_KERNEL_ARGS)
{
  Int i;

  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  BoundaryConditions<Real>* bc = space->bc;
  Param<Type>* param = space->param;

  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* QL = &space->q[left_cv*nvars];    //boundary point
  Type* QR = &space->q[right_cv*nvars];   //ghost point

  Int bcType = bc->GetBCType(factag); 
  BCObj<Real>* bcobj = bc->GetBCObj(factag);

  //TODO: how to handle this with the driver module
  //until then hack it out
  *ptrL = NULL;
  *size = 0;

  Type* Qref = (Type*)alloca(sizeof(Type)*nvars);
  bcobj->GetQref(Qref);

  if(bcType == Proteus_NoSlip){
    //find the most normal node to the wall
    Type* dx = (Type*)alloca(sizeof(Type)*3);
    Type* wallx = &m->xyz[left_cv*3];
    Int indx, indx1, indx2;
    Int normalNode = -1;
    indx1 = m->ipsp[left_cv];
    indx2 = m->ipsp[left_cv+1];
    Type dotmax = 0.0;
    for(indx = indx1; indx < indx2; indx++){
      Int pt = m->psp[indx];
      Type* ptx = &m->xyz[pt*3];
      Subtract(ptx, wallx, dx);
      Normalize(dx, dx);
      //this is to solve a templated type issue, standard DotProduct()
      //call doesn't like incompatible types
      Type dot = 0.0;
      for(i = 0; i < 3; i++){
	dot -= dx[i]*real(avec[i]);
      }
      //find the most normal point off the wall, use that distance for yplus
      if(real(dot) >= real(dotmax)){
	normalNode = pt;
	dotmax = dot;
      }
    }
    if(normalNode < 0){
      std::cerr << "WARNING: Normal node node found for point at "
		<< wallx[0] << " " << wallx[1] << " " << wallx[2] << std::endl;
    }
    //grid velocity at the wall
    Type vel[3];
    if(bcobj->movement == 1){
      //sliding movement
      vel[0] = bcobj->slipDirection[0]*bcobj->slipSpeed;
      vel[1] = bcobj->slipDirection[1]*bcobj->slipSpeed;
      vel[2] = bcobj->slipDirection[2]*bcobj->slipSpeed;
    }
    else if(bcobj->movement == 2){
      Type* pp = (Type*)alloca(sizeof(Type)*3);
      Type* pm = (Type*)alloca(sizeof(Type)*3);
      //get the point on the surface
      //find the point's location at the current timelevel
      Type* pt = &m->xyz[left_cv*3];
      //use a central difference to compute the velocity
      for(i = 0; i < 3; i++){
	pp[i] = pm[i] = pt[i];
      }
      Type* rotationPoint = (Type*)alloca(sizeof(Type)*3);
      Type* rotationAxis = (Type*)alloca(sizeof(Type)*3);
      Type omega = bcobj->omega;
      //we need to copy these over explicitly to take care of the complex templating issue
      for(i = 0; i < 3; i++){
	rotationPoint[i] = bcobj->rotationPoint[i];
	rotationAxis[i] = bcobj->rotationAxis[i];
      }

      RotatePoints(pp, rotationPoint, rotationAxis, -omega*param->dt, 1);
      RotatePoints(pm, rotationPoint, rotationAxis, +omega*param->dt, 1);
      
      //rotating movement
      for(i = 0; i < 3; i++){
	vel[i] = 0.5*(pp[i] - pm[i])/param->dt;
      }
    }
    else{
      //no movement
      vel[0] = vel[1] = vel[2] = 0.0;
    }
    //modify the velocities we expect via linear interpolation from 
    //freestream values, which is what we initialize to 
    if(space->iter < bcobj->bleedSteps){
      for(i = 0; i < 3; i++){
	Type V = param->GetVelocity(space->iter);
	Type Mach = V;
	Type fs = Mach*param->flowdir[i];
	Type dv = fs - vel[i]; 
	vel[i] = vel[i] + (fs - dv/(Type)bcobj->bleedSteps*(Type)space->iter);
      }
    }
    //add in the grid speed velocities
    vel[0] += velw[0];
    vel[1] += velw[1];
    vel[2] += velw[2];
    //add a gaussian source term velocity contribution to z-direction
    Bool holeCut = false;
    Type tmp_src = 0.0;
    Type* densities = new Type[param->massfractions.size()];
    if(param->gaussianSource && bc->GetBCType(factag) == param->gaussianBCid){
      Type* src = space->GetFieldData("gaussian", FIELDS::STATE_NONE);
      tmp_src = src[left_cv];
      vel[2] += src[left_cv]*param->gaussianVelAmpl;
      if(real(Abs(src[left_cv])) >= 1.0e-15){
	holeCut = true;
      }
    }
    if(!holeCut){
      eqnset->ModifyViscousWallJacobian(QL, QR, vel, cvid, space->crs, normalNode, bcobj->twall);
    }
    else{
      Type injected = tmp_src*param->gaussianAmpl;
      for(i = 0; i < param->massfractions.size(); i++){
	Type reduction = injected*(Type)param->massfractions.at(i);
	//if we get to the species we want to inject
	if(param->gaussianEqn == i){
	  densities[i] = injected;
	}
	//this reduces each species mass at the inlet to avoid adding
	//a mass dense plume to the flow.
	else{
	  densities[i] = param->massfractions.at(i) - reduction;
	}
      }
      Type flowDirection[3];
      //Type mag = Normalize(vel, flowDirection);
      flowDirection[0] = flowDirection[1] = 0.0;
      flowDirection[2] = 1.0;
      eqnset->ModifyInternalInflowJacobian(cvid, space->crs, densities, flowDirection, vel[2]);
    }
    delete [] densities;
  }

  return;
}

template <class Type>
void Bkernel_BC_Res_Modify(B_KERNEL_ARGS)
{
  Int i;

  Mesh<Type>* m = space->m;
  EqnSet<Type>* eqnset = space->eqnset;
  BoundaryConditions<Real>* bc = space->bc;
  Param<Type>* param = space->param;

  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  BCObj<Real>* bcobj = bc->GetBCObj(factag);
  Int bcType = bc->GetBCType(factag); 

  //TODO: how to handle this with the driver module
  //until then hack it out
  *ptrL = NULL;
  *size = 0;

  Type* Qref = (Type*)alloca(sizeof(Type)*nvars);
  bcobj->GetQref(Qref);

  if(bcType == Proteus_NoSlip){
    //find the most normal node to the wall
    Type* dx = (Type*)alloca(sizeof(Type)*3);
    Type* wallx = &m->xyz[left_cv*3];
    Int indx, indx1, indx2;
    Int normalNode = -1;
    Type dotmax = 0.0;
    indx1 = m->ipsp[left_cv];
    indx2 = m->ipsp[left_cv+1];
    for(indx = indx1; indx < indx2; indx++){
      Int pt = m->psp[indx];
      Type* ptx = &m->xyz[pt*3];
      Subtract(ptx, wallx, dx);
      Normalize(dx, dx);
      //this is to solve a templated type issue, standard DotProduct()
      //call doesn't like incompatible types
      Type dot = 0.0;
      for(i = 0; i < 3; i++){
	dot -= dx[i]*real(avec[i]);
      }
      //find the most normal point off the wall, use that distance for yplus
      if(real(dot) >= real(dotmax)){
	normalNode = pt;
	dotmax = dot;
      }
    }
    if(normalNode < 0){
      std::cerr << "WARNING: Normal node node found for point at "
		<< wallx[0] << " " << wallx[1] << " " << wallx[2] << std::endl;
    }
    //grid velocity at the wall
    Type vel[3];
    if(bcobj->movement == 1){
      //sliding movement
      vel[0] = bcobj->slipDirection[0]*bcobj->slipSpeed;
      vel[1] = bcobj->slipDirection[1]*bcobj->slipSpeed;
      vel[2] = bcobj->slipDirection[2]*bcobj->slipSpeed;
    }
    else if(bcobj->movement == 2){
      Type* pp = (Type*)alloca(sizeof(Type)*3);
      Type* pm = (Type*)alloca(sizeof(Type)*3);
      //get the point on the surface
      //find the point's location at the current timelevel
      Type* pt = &m->xyz[left_cv*3];
      //use a central difference to compute the velocity
      for(i = 0; i < 3; i++){
	pp[i] = pm[i] = pt[i];
      }
      Type* rotationPoint = (Type*)alloca(sizeof(Type)*3);
      Type* rotationAxis = (Type*)alloca(sizeof(Type)*3);
      Type omega = bcobj->omega;
      //we need to copy these over explicitly to take care of the complex templating issue
      for(i = 0; i < 3; i++){
	rotationPoint[i] = bcobj->rotationPoint[i];
	rotationAxis[i] = bcobj->rotationAxis[i];
      }

      RotatePoints(pp, rotationPoint, rotationAxis, -omega*param->dt, 1);
      RotatePoints(pm, rotationPoint, rotationAxis, +omega*param->dt, 1);
      
      //rotating movement
      for(i = 0; i < 3; i++){
	vel[i] = 0.5*(pp[i] - pm[i])/param->dt;
      }
    }
    else{
      //no movement
      vel[0] = vel[1] = vel[2] = 0.0;
    }
    //modify the velocities we expect via linear interpolation from 
    //freestream values, which is what we initialize to 
    if(space->iter < bcobj->bleedSteps){
      for(i = 0; i < 3; i++){
	Type V = param->GetVelocity(space->iter);
	Type Mach = V;
	Type fs = Mach*param->flowdir[i];
	Type dv = fs - vel[i]; 
	vel[i] = vel[i] + (fs - dv/(Type)bcobj->bleedSteps*(Type)space->iter);
      }
    }
    //add in the grid speed velocities
    vel[0] += velw[0];
    vel[1] += velw[1];
    vel[2] += velw[2];
    //add a gaussian source term velocity contribution to z-direction
    Bool holeCut = false;
    Type* densities = new Type[param->massfractions.size()];
    Type tmp_src = 0.0;
    if(param->gaussianSource && bc->GetBCType(factag) == param->gaussianBCid){
      Type* src = space->GetFieldData("gaussian", FIELDS::STATE_NONE);
      tmp_src = src[left_cv];
      vel[2] += src[left_cv]*param->gaussianVelAmpl;
      if(real(Abs(src[left_cv])) >= 1.0e-15){
	holeCut = true;
      }
    }
    if(!holeCut){
      eqnset->ModifyViscousWallResidual(&space->crs->b[cvid*neqn], vel, normalNode, bcobj->twall);
    }
    else{
      Type injected = tmp_src*param->gaussianAmpl;
      for(i = 0; i < param->massfractions.size(); i++){
	Type reduction = injected*(Type)param->massfractions.at(i);
	//if we get to the species we want to inject
	if(param->gaussianEqn == i){
	  densities[i] = injected;
	}
	//this reduces each species mass at the inlet to avoid adding
	//a mass dense plume to the flow.
	else{
	  densities[i] = param->massfractions.at(i) - reduction;
	}
      }
      Type flowDirection[3];
      //Type mag = Normalize(vel, flowDirection);
      Type pressure = 0.0;
      flowDirection[0] = flowDirection[1] = 0.0;
      flowDirection[2] = 1.0;
      eqnset->ModifyInternalInflowResidual(&space->crs->b[cvid*neqn], densities, flowDirection, vel[2]);
    }
    delete [] densities;
  }

}

template <class Type, class Type2, class Type3>
void CalculateBoundaryVariables(EqnSet<Type>* eqnset, Mesh<Type3>* m, SolutionSpace<Type3>* space, 
				Type* QL, Type* QR, Type* Qinf, Type* avec, Int bcType, 
				Int eid, BCObj<Type2>* bcobj, Type3 vdotn, Type3* velw)
{
  Int i;
  Int neqn = eqnset->neqn;
  Int nvars = eqnset->nauxvars + neqn;
  Param<Type>* param = eqnset->param;
  BoundaryConditions<Real>* bc = space->bc;
  Int left_cv = m->bedges[eid].n[0];
  Type3* beta = space->GetFieldData("beta", FIELDS::STATE_NONE);
  Type betaL = beta[left_cv];

  if(bcType == Proteus_ParallelBoundary){
    //do nothing this is a parallel updated edge/node
    return;
  }
  else if(bcType == Proteus_SonicInflow || bcType == Proteus_Dirichlet){
    for(i = 0; i < nvars; i++){
      //TODO: implement this correctly with param file
      Type* Qref = (Type*)alloca(sizeof(Type2)*nvars);
      bcobj->GetQref(Qref);
      QR[i] = QL[i] = Qref[i];
    }
  }
  else if(bcType == Proteus_SonicOutflow || bcType == Proteus_Neumann){
    for(i = 0; i < neqn; i++){
      QR[i] = QL[i];
    }
  }
  else if(bcType == Proteus_FarField){
    eqnset->GetFarfieldBoundaryVariables(QL, QR, Qinf, avec, vdotn, betaL);
  }
  else if(bcType == Proteus_FarFieldViscous){
    //modify qinf to follow the powerlaw assumption b/c it intersects
    //with a viscous surface, this aids stability of the solver
    Type3* dist = space->GetFieldData("wallDistance", FIELDS::STATE_NONE);
    Type d = dist[left_cv];
    Type one = 1.0;
    Type ubar = PowerLawU(one, d, param->Re);
    //now scale qinf velocity terms by ubar
    Int imom = eqnset->GetMomentumLocation();
    if(real(ubar) < 1.0){
      for(i = 0; i < 3; i++){
	Qinf[imom+i] = ubar*Qinf[imom+i];
      }
    }
    eqnset->GetFarfieldBoundaryVariables(QL, QR, Qinf, avec, vdotn, betaL);
  }
  else if(bcType == Proteus_ImpermeableWall || bcType == Proteus_Symmetry){
    eqnset->GetInviscidWallBoundaryVariables(QL, QR, Qinf, avec, vdotn, betaL);
  }
  else if(bcType == Proteus_InternalInflow || bcType == Proteus_InternalInflowDuctBL){
    Type pressure = bcobj->backPressure;
    Type* densities = NULL;
    if(bcobj->massFractions != NULL){
      densities = new Type[bcobj->nspecies];
      for(i = 0; i < bcobj->nspecies; i++){
	densities[i] = bcobj->density*bcobj->massFractions[i];
      }
    }
    Type flowDirection[3];
    //if the default settings haven't been changed use the normal to inlet condition
    if(bcobj->flowDirection[0] == -1.0 && bcobj->flowDirection[1] == -1.0){
      //our boundary edges all point outward... we need an inward pointing normal
      flowDirection[0] = -avec[0];
      flowDirection[1] = -avec[1];
      flowDirection[2] = -avec[2];
    }
    else{
      flowDirection[0] = bcobj->flowDirection[0];
      flowDirection[1] = bcobj->flowDirection[1];
      flowDirection[2] = bcobj->flowDirection[2];
    }
    Type vel = bcobj->velocity;
    Type ubar = 1.0;
    if(bcType == Proteus_InternalInflowDuctBL){
      //modify qinf to follow the powerlaw assumption b/c it intersects
      //with a viscous surface, this aids stability of the solver
      Type3* dist = space->GetFieldData("wallDistance", FIELDS::STATE_NONE);
      Type d = dist[left_cv];
      Type one = 1.0;
      ubar = PowerLawU(one, d, param->Re);
    }
    eqnset->GetInternalInflowBoundaryVariables(QL, QR, Qinf, pressure, densities, 
					       flowDirection, ubar*vel);
    delete [] densities;
  }
  else if(bcType == Proteus_InternalOutflow){
    Type pressure = bcobj->backPressure;
    //this only works for eqnsets where gamma is constant
    Type gamma = eqnset->param->gamma;
    eqnset->GetInternalOutflowBoundaryVariables(QL, QR, pressure, gamma);
  }
  else if(bcType == Proteus_TotalTempAndPressure){
    Type pressure = bcobj->backPressure;
    Type T = bcobj->twall;
    Type flowDirection[3];
    //if the default settings haven't been changed use the normal to inlet condition
    if(bcobj->flowDirection[0] == -1.0 && bcobj->flowDirection[1] == -1.0){
      //our boundary edges all point outward... we need an inward pointing normal
      flowDirection[0] = -avec[0];
      flowDirection[1] = -avec[1];
      flowDirection[2] = -avec[2];
    }
    else{
      flowDirection[0] = bcobj->flowDirection[0];
      flowDirection[1] = bcobj->flowDirection[1];
      flowDirection[2] = bcobj->flowDirection[2];
    }
    eqnset->GetTotalTempPressureBoundaryVariables(QL, QR, pressure, T, flowDirection);
  }
  else if(bcType == Proteus_NoSlip){
    //find the most normal node to the wall
    Type3* dx = (Type3*)alloca(sizeof(Type3)*3);
    Type3* wallx = &m->xyz[left_cv*3];
    Int indx, indx1, indx2;
    Int normalNode = -1;
    Type dotmax = 0.0;
    indx1 = m->ipsp[left_cv];
    indx2 = m->ipsp[left_cv+1];
    for(indx = indx1; indx < indx2; indx++){
      Int pt = m->psp[indx];
      Type3* ptx = &m->xyz[pt*3];
      Subtract(ptx, wallx, dx);
      Normalize(dx, dx);
      //this is to solve a templated type issue, standard DotProduct()
      //call doesn't like incompatible types
      Type3 dot = 0.0;
      for(i = 0; i < 3; i++){
	dot -= dx[i]*real(avec[i]);
      }
      //find the most normal point off the wall, use that distance for yplus
      if(real(dot) >= real(dotmax)){
	normalNode = pt;
	dotmax = dot;
      }
    }
    if(normalNode < 0){
      std::cerr << "WARNING: Normal node node found for point at "
		<< wallx[0] << " " << wallx[1] << " " << wallx[2] << std::endl;
    }
    //grid velocity at the wall
    Type vel[3];
    if(bcobj->movement == 1){
      //sliding movement
      vel[0] = bcobj->slipDirection[0]*bcobj->slipSpeed;
      vel[1] = bcobj->slipDirection[1]*bcobj->slipSpeed;
      vel[2] = bcobj->slipDirection[2]*bcobj->slipSpeed;
    }
    else if(bcobj->movement == 2){
      Type* pp = (Type*)alloca(sizeof(Type)*3);
      Type* pm = (Type*)alloca(sizeof(Type)*3);
      //get the point on the surface
      //find the point's location at the current timelevel
      Type3* pt = &m->xyz[left_cv*3];
      //use a central difference to compute the velocity
      for(i = 0; i < 3; i++){
	pp[i] = pm[i] = pt[i];
      }
      Type* rotationPoint = (Type*)alloca(sizeof(Type)*3);
      Type* rotationAxis = (Type*)alloca(sizeof(Type)*3);
      Type omega = bcobj->omega;
      //we need to copy these over explicitly to take care of the complex templating issue
      for(i = 0; i < 3; i++){
	rotationPoint[i] = bcobj->rotationPoint[i];
	rotationAxis[i] = bcobj->rotationAxis[i];
      }

      RotatePoints(pp, rotationPoint, rotationAxis, -omega*param->dt, 1);
      RotatePoints(pm, rotationPoint, rotationAxis, +omega*param->dt, 1);
      
      //rotating movement
      for(i = 0; i < 3; i++){
	vel[i] = 0.5*(pp[i] - pm[i])/param->dt;
      }
    }
    else{
      //no movement
      vel[0] = vel[1] = vel[2] = 0.0;
    }
    //modify the velocities we expect via linear interpolation from 
    //freestream values, which is what we initialize to 
    if(space->iter < bcobj->bleedSteps){
      for(i = 0; i < 3; i++){
	Type V = param->GetVelocity(space->iter);
	Type Mach = V;
	Type fs = Mach*param->flowdir[i];
	Type dv = fs - vel[i]; 
	vel[i] = vel[i] + (fs - dv/(Type)bcobj->bleedSteps*(Type)space->iter);
      }
    }
    //add in the grid speed velocities
    vel[0] += velw[0];
    vel[1] += velw[1];
    vel[2] += velw[2];
    Int nvars = eqnset->neqn + eqnset->nauxvars;
    Type nQ[nvars];
    for(Int kk = 0; kk < nvars; kk++){
      nQ[kk] = space->q[normalNode*nvars + kk];
    }
    Int factag = m->bedges[eid].factag;
    Bool holeCut = false;
    Type* densities = new Type[param->massfractions.size()];
    Type tmp_src = 0.0;
    if(param->gaussianSource && bc->GetBCType(factag) == param->gaussianBCid){
      Type3* src = space->GetFieldData("gaussian", FIELDS::STATE_NONE);
      tmp_src = src[left_cv];
      vel[2] += src[left_cv]*param->gaussianVelAmpl;
      if(real(Abs(src[left_cv])) >= 1.0e-15){
	holeCut = true;
      }
    }
    if(!holeCut){
      eqnset->GetViscousWallBoundaryVariables(QL, QR, vel, normalNode, nQ, bcobj->twall);
    }
    else{
      Type injected = tmp_src*param->gaussianAmpl;
      for(i = 0; i < param->massfractions.size(); i++){
	Type reduction = injected*(Type)param->massfractions.at(i);
	//if we get to the species we want to inject
	if(param->gaussianEqn == i){
	  densities[i] = injected;
	}
	//this reduces each species mass at the inlet to avoid adding
	//a mass dense plume to the flow.
	else{
	  densities[i] = param->massfractions.at(i) - reduction;
	}
      }
      Type flowDirection[3];
      //Type mag = Normalize(vel, flowDirection);
      Type pressure = 0.0;
      flowDirection[0] = flowDirection[1] = 0.0;
      flowDirection[2] = 1.0;
      eqnset->GetInternalInflowBoundaryVariables(QL, QR, Qinf, pressure, densities, 
						 flowDirection, vel[2]);
    }
    delete [] densities;
  }
  else if(bcType == Proteus_PitchingFarField){
    //we compute the farfield Qinf inflow angles externally to 
    //avoid repititve computation and simply apply the farfield BC
    //here, now compute the farfield BCs with the new AOA
    eqnset->GetFarfieldBoundaryVariables(QL, QR, Qinf, avec, vdotn, betaL);
  }
  else if(bcType == Proteus_HeatFlux){
    //find the most normal node to the wall
    Type3* dx = (Type3*)alloca(sizeof(Type3)*3);
    Type* ndx = (Type*)alloca(sizeof(Type)*3);
    Type3* wallx = &m->xyz[left_cv*3];
    Int indx, indx1, indx2;
    Int normalNode = -1;
    Type dotmax = 0.0;
    indx1 = m->ipsp[left_cv];
    indx2 = m->ipsp[left_cv+1];
    for(indx = indx1; indx < indx2; indx++){
      Int pt = m->psp[indx];
      Type3* ptx = &m->xyz[pt*3];
      Subtract(ptx, wallx, dx);
      Normalize(dx, dx);
      //this is to solve a templated type issue, standard DotProduct()
      //call doesn't like incompatible types
      Type3 dot = 0.0;
      for(i = 0; i < 3; i++){
	dot -= dx[i]*real(avec[i]);
      }
      //find the most normal point off the wall, use that distance for yplus
      if(real(dot) >= real(dotmax)){
	normalNode = pt;
	dotmax = dot;
	for(Int j = 0; j < 3; ++j){
	  Subtract(ptx, wallx, dx);
	  ndx[j] = dx[j];
	}
      }
    }
    if(normalNode < 0){
      std::cerr << "WARNING: Normal node node found for point at "
		<< wallx[0] << " " << wallx[1] << " " << wallx[2] << std::endl;
    }
    Int nvars = eqnset->neqn + eqnset->nauxvars;
    Type nQ[nvars];
    for(Int kk = 0; kk < nvars; kk++){
      nQ[kk] = space->q[normalNode*nvars + kk];
    }
    Type flux = bcobj->flux;
    eqnset->GetHeatFluxBoundaryVariables(QL, QR, nQ, ndx, flux);
  }
  else if(bcType == Proteus_Isothermal){
    Type Twall = bcobj->twall;
    eqnset->GetIsothermalBoundaryVariables(QL, QR, Twall);
  }
  else if(bcType == Proteus_PythonBC){
    Type wallx[3];
    wallx[0] = m->xyz[left_cv*3 + 0];
    wallx[1] = m->xyz[left_cv*3 + 1];
    wallx[2] = m->xyz[left_cv*3 + 2];
    eqnset->GetPythonBoundaryVariables(QL, QR, wallx, avec);
  }
  else{
    std::cout << "WARNING: EID for " << eid << " not found! " << std::endl;
    return;
  }

  //compute aux vars based on changes
  eqnset->ComputeAuxiliaryVariables(QR);
  //compute aux vars if we hard set the wall
  eqnset->ComputeAuxiliaryVariables(QL);
  
  return;
}

template <class Type>
void UpdateBCs(SolutionSpace<Type>* space)
{
  EqnSet<Type>* eqnset = space->eqnset;
  Param<Type>* param = space->param;
  BoundaryConditions<Real>* bc = space->bc;

  //check to see if our pitching BC is enabled, if it is then we have to 
  //modify the flow direction and update Qinf (this is done below)
  Int i;
  Bool pitching = false;
  for(i = 0; i <= bc->num_bcs; i++){
    Int factag = i;
    Int bcType = bc->GetBCType(factag);
    if(bcType == Proteus_PitchingFarField){
      pitching = true;
      std::cout << "RUNNING PITCHING MOTION" << std::endl;
      break;
    }
  }

  if(pitching){
    //the number of steps before we start pitching
    Int pitchingStart = 0;
    //the AOA in degrees about which we are pitching the airfoil
    Type alphaMean = 4.86;
    //the AOA in degrees (+/-) by which we oscillate
    Type alphaVary = 2.44;
    //this is the dimensional frequency omega (Hz - 1/s)
    Type dimFreq = 50.32;
    if(space->iter > pitchingStart){
      Type pi = 3.14159265359;
      //convert alphas to radians for sin/cos functions
      Type alphaMeanRad = alphaMean*2.0*pi/360.0;
      Type alphaVaryRad = alphaVary*2.0*pi/360.0;
      Type dimOmega = dimFreq*2.0*pi;      //(rad/s)
      //non-dimensional circular frequency
      Type nondimOmega = dimOmega*param->ref_time;
      //compute omega*t
      Type radPos = nondimOmega*(Type)(space->iter - pitchingStart)*(Type)param->dt;
      //since sine varies from 0 -> 1 -> 0 -> -1 over 2*pi
      //we can simply multiply by our alpha variation to get
      //a smooth pitching action
      Type alpha = alphaMeanRad + sin(radPos)*alphaVaryRad;
      param->flowdir[0] = real(cos(alpha));
      param->flowdir[1] = real(sin(alpha));
      std::cout << "Pitching angle: " << alpha*360.0/(2.0*pi) << std::endl;
    }
  }

  //update Qinf in case our farfield values are ramping or dynamic
  eqnset->UpdateQinf();

  Kernel<Type> BC(BC_Kernel);
  
  //call driver to loop over and set boundary conditions
  BdriverNoScatter(space, BC, eqnset->neqn+eqnset->nauxvars, NULL);
  return;
}


template <class Type>
Int GetNodesOnSymmetryPlanes(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int** nodes)
{
  Int ncount = 0;
  Int memChunk = 50;
  Int memSize = 50;

  *nodes = new Int[memChunk];

  Int nbedge = m->GetNumBoundaryEdges();
  Int ngedge = m->GetNumParallelEdges();
  for(Int eid = 0; eid < nbedge+ngedge; eid++){
    Int factag = m->bedges[eid].factag;
    Int bcType = bc->GetBCType(factag);
    Int node = m->bedges[eid].n[0];
    if(bcType == Proteus_Symmetry){
      //if node is connected to symmetry/parallel plane 
      //ensure it is not on another type of surface as well
      //this keeps us from moving solid boundaries around
      Int indx1 = m->ibesp[node];
      Int indx2 = m->ibesp[node+1];
      bool clean = true;
      for(Int i = indx1; i < indx2; i++){
	Int leid = m->besp[i];
	Int lfactag = m->bedges[leid].factag;
	Int lbcType = bc->GetBCType(lfactag);
	if(lbcType != Proteus_Symmetry){
	  //mark dirty and move to next node
	  clean = false;
	  break;
	}
      }
      if(clean){
	//if the node was only connected to parallel/symmetry planes
	//add it to the list
	if(memSize <= ncount){
	  MemResize(nodes, memSize, memSize+memChunk);
	  memSize += memChunk;
	}
	(*nodes)[ncount] = node;
	ncount++;
      }
    }
  }

  //do final memory resizing operation
  MemResize(nodes, memSize, ncount);

  return ncount;
}

template <class Type>
Int GetNodesOnMovingBC(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int** nodes)
{
  Int ncount = 0;
  Int memChunk = 50;
  Int memSize = 50;

  *nodes = new Int[memChunk];

  //check for duplicates as we go along
  Int nnode = m->GetNumNodes();
  Bool* dupeCheck = new Bool[nnode];
  for(Int i = 0; i < nnode; i++){
    dupeCheck[i] = false;
  }
  
  Int nbedge = m->GetNumBoundaryEdges();
  Int ngedge = m->GetNumParallelEdges();
  for(Int eid = 0; eid < nbedge+ngedge; eid++){
    Int factag = m->bedges[eid].factag;
    Int node = m->bedges[eid].n[0];
    BCObj<Real>* bcobj = bc->GetBCObj(factag);
    Int bcType = bc->GetBCType(factag);
    if( bcobj->IsMovingBC() && !dupeCheck[node]){
      //if the node was only connected to parallel/symmetry planes
      //add it to the list
      if(memSize <= ncount){
	MemResize(nodes, memSize, memSize+memChunk);
	memSize += memChunk;
      }
      (*nodes)[ncount] = node;
      dupeCheck[node] = true;
      ncount++;
    }
  }

  //do final memory resizing operation
  MemResize(nodes, memSize, ncount);

  delete [] dupeCheck;

  return ncount;
}
 


template <class Type>
Int GetNodesOnMovingBC(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int beta,
		       Int** nodes)
{
  Int i;
  Int ncount = 0;
  Int memChunk = 50;
  Int memSize = 50;

  *nodes = new Int[memChunk];

  //check for duplicates as we go along
  Int nnode = m->GetNumNodes();
  Bool* dupeCheck = new Bool[nnode];
  for(i = 0; i < nnode; i++){
    dupeCheck[i] = false;
  }
 
  Int nbedge = m->GetNumBoundaryEdges();
  Int ngedge = m->GetNumParallelEdges();
  for(Int eid = 0; eid < nbedge+ngedge; eid++){
    Int factag = m->bedges[eid].factag;
    int node = m->bedges[eid].n[0];
    BCObj<Real>* bcobj= bc->GetBCObj(factag);
    if( bcobj->movingBC >= 0 && !dupeCheck[node]){
      if(memSize <= ncount){
	MemResize(nodes, memSize, memSize+memChunk);
	memSize += memChunk;
      }
      (*nodes)[ncount] = node;
      dupeCheck[node] = true;
      ncount++;
    } 
  }

  //do final memory resizing operation
  MemResize(nodes, memSize, ncount);

  delete [] dupeCheck;

  return ncount;
}
 
template <class Type>
Int GetNodesOnBCType(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int bcType,
		   Int** nodes)
{
  Int ncount = 0;
  Int memChunk = 50;
  Int memSize = 50;

  *nodes = new Int[memChunk];

  //check for duplicates as we go along
  Int nnode = m->GetNumNodes();
  Bool* dupeCheck = new Bool[nnode];
  for(Int i = 0; i < nnode; i++){
    dupeCheck[i] = false;
  }
 
  Int nbedge = m->GetNumBoundaryEdges();
  Int ngedge = m->GetNumParallelEdges();
  for(Int eid = 0; eid < nbedge+ngedge; eid++){
    Int factag = m->bedges[eid].factag;
    Int node = m->bedges[eid].n[0];
    Int bcTypel = bc->GetBCType(factag);
    if(bcTypel== bcType && !dupeCheck[node]){
      if(memSize <= ncount){
	MemResize(nodes, memSize, memSize+memChunk);
	memSize += memChunk;
      }
      (*nodes)[ncount] = node;
      dupeCheck[node] = true;
      ncount++;
    } 
  }

  //do final memory resizing operation
  MemResize(nodes, memSize, ncount);

  delete [] dupeCheck;

  return ncount;
}

template <class Type>
void GetSElemsOnBCType(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int bcType, 
		     std::vector<Element<Type>*>& elementList)
{
  Int* bcnodes;
  Int nbcnodes = GetNodesOnBCType(m, bc, bcType, &bcnodes);

  std::set<Element<Type>*, ElementPtrComp<Type> > uniqueElems;

  for(Int i = 0; i < nbcnodes; ++i){
    Int node = bcnodes[i];
    for(const Int* indx = m->SelspBegin(node); indx != m->SelspEnd(node); ++indx){
      Int selemid = *indx;;
      Int factag = m->elementList[selemid]->GetFactag();
      Int bcTypeLocal = bc->GetBCType(factag);
      if(bcTypeLocal == bcType){
	//attempt to insert into the set, it will be rejected if non-unique
	uniqueElems.insert(m->elementList[selemid]);
      }
    }
  }

  delete [] bcnodes;

  for(typename std::set<Element<Type>*>::iterator it = uniqueElems.begin(); it != uniqueElems.end(); ++it){
    elementList.push_back(*it);
  }

}

template <class Type>
void SetBCStatFlags(Mesh<Type>* m, BoundaryConditions<Real>* bc)
{
  Int nbedge = m->GetNumBoundaryEdges();
  Int ngedge = m->GetNumParallelEdges();
  for(Int eid = 0; eid < nbedge+ngedge; eid++){
    Int factag = m->bedges[eid].factag;
    Int bcType = bc->GetBCType(factag);
    Int node = m->bedges[eid].n[0];
    if(bcType == Proteus_NoSlip){
      m->cvstat[node].SetViscous();
    }
    else if(bcType == Proteus_Dirichlet){
      m->cvstat[node].SetDirichlet();
    }
  }
}
