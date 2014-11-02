#include "octree.h"

template <class Type>
Sensors<Type>::Sensors(std::string filename, Int neqn, Int nstride, Type* Q, 
		       PObj<Type>* p, Int* ipsp, Int* psp, Type* mxyz, 
		       Int nnodes) : neqn(neqn), nstride(nstride), psp(psp), ipsp(ipsp), 
				     mxyz(mxyz), mnnode(nnodes), q(Q), xyz(NULL), 
				     values(NULL), domain(NULL), nodes(NULL), p(p)
  
{
  init = false;

  //this is a hardset of the bounding box tolerance before
  //an exhaustive search.
  tol = 5.0e-1;

  ReadFile(filename);

  //clear files from old data
  std::stringstream ss;
  std::ofstream fout;
  std::string file;
  
  //clear the end of the filename
  size_t pos = filename.rfind(".sensors");
  filename = filename.substr(0, pos);

  if(p->GetRank() == 0){
    for(Int i = 0; i < count; i++){
      ss.clear();
      ss.str("");
      ss << i;
      file = filename + ".Sensor-" + ss.str() + ".sens";
      fout.open(file.c_str(), std::ios::trunc);
      if(!fout.is_open()){
	std::cerr << "WARNING: SENSORS_IO: " << file << " is not found for writing!" << std::endl;
	return;
      }
      fout.close();
    }
  }  

  return;
}

template <class Type>
Sensors<Type>::~Sensors()
{
  delete [] xyz;
  delete [] values;
  delete [] domain;
  delete [] nodes;

  return;
}

template <class Type>
Int Sensors<Type>::ReadFile(std::string filename)
{
  Int i, j;
  Int err = 0;
  std::ifstream fin;

  std::cout << "SENSORS_IO: reading sensors file " << filename << std::endl;

  fin.open(filename.c_str());
  if(!fin.is_open()){
    std::cerr << "WARNING: SENSORS_IO: " << filename << " is not found for reading!" << std::endl;
    count = 0;
    return(-1);
  }

  fin >> count;

  std::cout << "SENSORS_IO: reading " << count << " sensors from file " << std::endl;
  
  //allocate xyz memory
  xyz = new Type[3*count];

  //allocate sampling memory
  values = new Type[count*neqn];

  //allocate domain list memory
  domain = new Int[count];

  //allocate node list memory
  nodes = new Int[count];

  //read in xyz location of each sensor
  for(i = 0; i < count; i++){
    for(j = 0; j < 3; j++){
      fin >> xyz[i*3 + j];
    }
  }

  std::cout << "SENSORS_IO: Read sensors" << std::endl;
  for(i = 0; i < count; i++){
    std::cout << i << ": ";
    for(j = 0; j < 3; j++){
      std::cout << xyz[i*3 + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  fin.close();

  init = true;

  return err;
}

template <class Type>
void Sensors<Type>::Search()
{
  Int j;
  RealInt* searchArray = (RealInt*)alloca(sizeof(RealInt)*count);

  //octree parameters
  Int minNodes = 10;
  Int maxDepth = 10;
  Octree<Type> octree(mxyz, mnnode, maxDepth, minNodes);

  //search for closest node at local level
  for(j = 0; j < count; j++){
    searchArray[j].rank = p->GetRank();
    searchArray[j].val = real(octree.FindDistance(&xyz[j*3 + 0], nodes[j]));
  }

  //sync and check for best global match
  MPI_Allreduce(MPI_IN_PLACE, searchArray, count, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD); 

  //now we need to reduce the list of local nodes such that if the local
  //process does not contain the smallest distance, i.e. the node 
  //we set the correct domain lookup and the node value to < 0
  for(j = 0; j < count; j++){
    if(searchArray[j].rank != p->GetRank()){
      nodes[j] = -999;
    }
    domain[j] = searchArray[j].rank;
  }

  return;
}

template <class Type>
void Sensors<Type>::Sample()
{
  Int j;
  Int indx1, indx2;
  std::vector<Int> stencil;
  Int nnode;
  MPI_Datatype mpit;

  //get mpi datatype to send
  mpit = MPI_GetType(q[0]);

  for(j = 0; j < count; j++){
    stencil.clear();
    if(domain[j] == p->GetRank()){
      indx1 = ipsp[nodes[j]];
      indx2 = ipsp[nodes[j]+1];
      for(Int indx = indx1; indx < indx2; indx++){
	stencil.push_back(psp[indx]);
      }
      //now eliminate ghost nodes since these live in psp array
      for(UInt indx = 0; indx < stencil.size(); indx++){
	if(stencil[indx] >= mnnode){
	  stencil.erase(stencil.begin()+indx);
	}
      }
      nnode = stencil.size();
      //allocate some temporary memory for interpolation
      Int* nodelist = (Int*)alloca(sizeof(Int)*nnode);
      for(Int i = 0; i < nnode; i++){
	nodelist[i] = stencil[i];
      }
      //this is the power parameter for Shepard's method
      Type p = 2.0;
      InverseWeighting(nnode, nodelist, mxyz, q, &xyz[j*3], &values[j*neqn], neqn, nstride, p);
    }
  }

  //sync the values across all processes, this is useful if we are doing
  //some kind of goal seek for sample data
  for(j = 0; j < count; j++){
    MPI_Bcast(&values[j*neqn], neqn, mpit, domain[j], MPI_COMM_WORLD);
  }

  return;
}

template <class Type>
Int Sensors<Type>::Write(std::string filename, Int stepTag)
{
  Int err = 0;
  std::ofstream fout;
  std::stringstream ss;
  std::string file;
  Int prec = 16;

  if(p->GetRank() == 0){
    for(Int i = 0; i < count; i++){
      ss.clear();
      ss.str("");
      ss << i;
      file = filename + ".Sensor-" + ss.str() + ".sens";
      fout.open(file.c_str(), std::ios::app);
      fout.setf(std::ios::scientific);
      fout.precision(prec);
      if(!fout.is_open()){
	std::cerr << "WARNING: SENSORS_IO: " << file << " is not found for writing!" << std::endl;
	return(-1);
      }
      fout << stepTag << " ";
      for(Int j = 0; j < neqn; j++){
	fout << values[i*neqn + j] << "\t";
      }
      fout << std::endl;
      fout.close();
    }
  }  

  return err;
}
