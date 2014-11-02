#include "etypes.h"
#include "mesh.h"
#include "crsmatrix.h"
#include "exceptions.h"

template <class Type>
PObj<Type>::PObj()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  
  //bad values in case something fails
  rank = -1;
  np = -1;
  //set object copies of np and rank
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  
  //allocate some memory we'll need for most things
  commCountsRecv = new Int[np];
  commCountsSend = new Int[np];
  commOffsetsRecv = new Int[np];
  commOffsetsSend = new Int[np];
  sRequest= new MPI_Request[np];
  rRequest= new MPI_Request[np];
  nodePackingList = new Int*[np];
  
  mapsBuilt = false;
  
  timers.InitList(2);
  timers.CreateTimer("CommTimer");
  timers.CreateTimer("CommTotalTimer");
}

template <class Type>
PObj<Type>::~PObj()
{
  if(mapsBuilt != 0){
    //delete allocated memory if commMaps were built
    delete [] nodePackingListData;
  }
  //delete standard allocated memory
  delete [] commCountsRecv;
  delete [] commCountsSend;
  delete [] commOffsetsRecv;
  delete [] commOffsetsSend;
  delete [] nodePackingList;
  delete [] sRequest;
  delete [] rRequest;
  
}

template <class Type>
Int PObj<Type>::TransposeCommCRS(CRSMatrix<Type>* crs)
{
  Int i, j, sumrecv, sumsend;
  MPI_Status status;
  MPI_Datatype mpit;
  MPI_Datatype mpitint;
  i = 0; //silence compiler warning

  //get mpi datatype to send
  Type check = 0;
  mpit = MPI_GetType(check);
  mpitint = MPI_GetType(i);


  //DESCRIPTION OF EXPECTED STORAGE PATTERN:
  //=========================================
  //
  //if process 0 needs a matrix swapped on row 200, col 300
  //then the ghost point is 300
  //this can be mapped to a process (say 3) with a local id (say 150)
  //
  //then the matrix process zero needs is on row 150 in process 3's CRS matrix
  // it will also be a ghost not local, since the request was a ghost
  // and all local blocks will have been transposed already
  //
  //use global node number to match the matrix we need
  //load it into the array to be sent
  //

  Int* globalToLocal;
  Int* localToGlobal;
  Int globalCount;
  globalCount = BuildGlobalMaps(&localToGlobal, &globalToLocal);

  //pick out the information we need from the CRS structure
  Int maxrow;
  Int maxconn;
  maxconn = crs->GetMaxConnectivity(&maxrow);
  Int* ghosts = new Int[maxconn];
  Int gcount;
  Int counter[np];
  Int countsSend[np];
  Int countsRecv[np];
  Int row, col;
  Int owner;
  Int OffsetsSend[np];
  Int OffsetsRecv[np];
  
  //set up indexes for receiving
  //we need to count the expected size of the buffer
  //b/c of each edge containing a jacobian this is by definition
  // >= the size of our normal q vector calls
  sumrecv = 0;
  for(i = 0; i < np; i++){
    countsSend[i] = 0;
    countsRecv[i] = 0;
    counter[i] = 0;
    OffsetsSend[i] = 0;
    OffsetsRecv[i] = 0;
  }
  for(row = 0; row < nnode; row++){
    gcount = crs->GetGhostsInRow(row, ghosts);
    sumrecv += gcount;
    for(i = 0; i < gcount; i++){
      owner = m->gNodeOwner[ghosts[i]-nnode];
      countsRecv[owner]++;
    }
  }
  //we now know the size of the buffer we need to recv data
  //because we just counted the blocks
  //setup the recv offset map
  OffsetsRecv[0] = 0;
  for(i = 1; i < np; i++){
    OffsetsRecv[i] = OffsetsRecv[i-1] + countsRecv[i-1];
  }

  //now we need to do a parallel comm so that we can count the number of 
  //blocks the will be requested of each process from other processes
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&countsRecv[i], 1, mpitint, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  //recv the counts from other processes
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&countsSend[i], 1, mpitint, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  //count the send size
  sumsend = 0;
  for(i = 0; i < np; i++){
    sumsend += countsSend[i];
  }
  //set the sending offsets
  OffsetsSend[0] = 0;
  for(i = 1; i < np; i++){
    OffsetsSend[i] = OffsetsSend[i-1] + countsSend[i-1];
  }

  //setup utility buffers for sends/recvs, we want to send/recv a contiguous block
  Type* sendBuff[np];
  Type* recvBuff[np];
  Type* sendBuffData;
  Type* recvBuffData;  

  
  //allocate memory for buffers
  sendBuffData = new Type[crs->neqn2*sumsend];
  recvBuffData = new Type[crs->neqn2*sumrecv];

  for(i = 0; i < np; i++){
    sendBuff[i] = &sendBuffData[OffsetsSend[i]*crs->neqn2];
    recvBuff[i] = &recvBuffData[OffsetsRecv[i]*crs->neqn2];
  }

  //now we need to sync a list of cols (by global id) that each process needs to load
  //since we are syncing a tranpose this is the current row of the matrix it will be 
  //swapped with

  Int* sendCol[np];
  Int* sendRow[np];
  Int* sendColData = new Int[sumsend];
  Int* sendRowData = new Int[sumsend];
  Int* recvCol[np];
  Int* recvColG[np];
  Int* recvColData = new Int[sumrecv];
  Int* recvColDataG = new Int[sumrecv];
  Int* recvRow[np];
  Int* recvRowG[np];
  Int* recvRowData = new Int[sumrecv];
  Int* recvRowDataG = new Int[sumrecv];

  for(i = 0; i < np; i++){
    sendCol[i] = &sendColData[OffsetsSend[i]];
    sendRow[i] = &sendRowData[OffsetsSend[i]];
    recvCol[i] = &recvColData[OffsetsRecv[i]];
    recvRow[i] = &recvRowData[OffsetsRecv[i]];
    recvColG[i] = &recvColDataG[OffsetsRecv[i]];
    recvRowG[i] = &recvRowDataG[OffsetsRecv[i]];
  }

  for(i = 0; i < np; i++){
    counter[i] = 0;
  }
  for(row = 0; row < nnode; row++){
    gcount = crs->GetGhostsInRow(row, ghosts);
    for(i = 0; i < gcount; i++){
      owner = m->gNodeOwner[ghosts[i]-nnode];
      recvCol[owner][counter[owner]] = ghosts[i];
      recvRow[owner][counter[owner]] = row;
      counter[owner]++;
    }
  }

  //sanity check
  for(i = 0; i < np; i++){
    if(counter[i] != countsRecv[i]){
      std::cerr << "PObj: TransposeCommCRS() counters don't match!" << std::endl;
      std::cerr << "PObj: process " << i << " counter[" << i << "] - " << counter[i] 
		<< " CountsRecv[" << i << "] - " << countsRecv[i] << std::endl;
      Abort << "Dead";
    }
  }

  //this list is still locally numbered, convert it to a global numbering
  for(i = 0; i < np; i++){
    for(j = 0; j < countsRecv[i]; j++){
      //this is NOT a bug the local row is the global column in transpose
      //and vice versa
      recvColG[i][j] = localToGlobal[recvRow[i][j]];
      recvRowG[i][j] = localToGlobal[recvCol[i][j]];
    }
  }

  //PARALLEL SYNC UP FOR THE COLS
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&recvColG[i][0], countsRecv[i], mpitint, 
		i, 0, MPI_COMM_WORLD, &sRequest[i]); 
    }
  }
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&sendCol[i][0], countsSend[i], mpitint, 
		i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  } 

  //PARALLEL SYNC UP FOR THE ROWS
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&recvRowG[i][0], countsRecv[i], mpitint, 
		i, 0, MPI_COMM_WORLD, &sRequest[i]); 
    }
  }
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&sendRow[i][0], countsSend[i], mpitint, 
		i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  } 

  //PACK THE BLOCKS
  //pack matrics from column and row as directed
  for(i = 0; i < np; i++){
    for(j = 0; j < countsSend[i]; j++){
      row = sendRow[i][j];
      col = sendCol[i][j];
      //these were sent with global indices, map 'em 
      row = globalToLocal[row];
      col = globalToLocal[col];
      Type* ptr = crs->GetPointer(row, col);
      memcpy(&sendBuff[i][j*crs->neqn2], ptr, sizeof(Type)*crs->neqn2);
    }
  }

  //PARALLEL SEND/RECV BLOCKS
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&sendBuff[i][0], countsSend[i]*crs->neqn2, mpit, 
		i, 0, MPI_COMM_WORLD, &sRequest[i]); 
    }
  }
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&recvBuff[i][0], countsRecv[i]*crs->neqn2, mpit, 
		i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  } 

  //UNPACK THE BLOCKS
  //unpack matrices from column and row as directed
  for(i = 0; i < np; i++){
    for(j = 0; j < countsRecv[i]; j++){
      row = recvRow[i][j];
      col = recvCol[i][j];
      Type* ptr = crs->GetPointer(row, col);
      memcpy(ptr, &recvBuff[i][j*crs->neqn2], sizeof(Type)*crs->neqn2);
    }
  }

  //And we're done, wrap it up!

  delete [] sendBuffData;
  delete [] recvBuffData;
  
  delete [] globalToLocal;
  delete [] localToGlobal;

  delete [] sendColData;
  delete [] sendRowData;
  delete [] recvColData;
  delete [] recvColDataG;
  delete [] recvRowData;
  delete [] recvRowDataG;

  delete [] ghosts;

  return 0;
}


template <class Type>
Int PObj<Type>::BuildGlobalMaps(Int** localToGlobal, Int** globalToLocal)
{
  Int i, j;
  MPI_Status status;
  Int localNodes[np];
  Int offsets[np];
  Int globalNodes;

  for(i = 0; i < np; i++){
    if(i == rank){
      localNodes[i] = nnode;
    }
    else{
      localNodes[i] = 0;
    }
  }

  //get number of real nodes on each process
  MPI_Allreduce(MPI_IN_PLACE, localNodes, np, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  //get number of real nodes in global mesh
  globalNodes = 0;
  for(i = 0; i < np; i++){
    globalNodes += localNodes[i];
  }
  
  //assign global numbers in order of process ids
  offsets[0] = 0;
  for(i = 1; i < np; i++){
    offsets[i] = offsets[i-1] + localNodes[i-1];
  }


  //allocate the lists
  *localToGlobal = new Int[nnode+gnode];
  //this allocation is a bit wasteful but saves arithmetic in lookups... we 
  //could also pass back the offsets and allocate this to length nnode
  *globalToLocal = new Int[globalNodes];
  
  //build what we can before the comm step
  for(i = 0; i < globalNodes; i++){
    (*globalToLocal)[i] = -1;
  }
  for(i = 0; i < nnode; i++){
    (*localToGlobal)[i] = offsets[rank] + i;
    (*globalToLocal)[offsets[rank] + i] = i;
  }

  //
  //now we need to pick up the gnodes in the local to global lists, rock on...
  //
  
  Int sum = 0;
  for(i = 0; i < np; i++){
    if(i != rank){
      sum += commCountsSend[i];
    }
  }
  Int* globalIdSendData = new Int[sum];
  Int* globalIdSend[np];
  for(i = 0; i < np; i++){
    globalIdSend[i] = &globalIdSendData[commOffsetsSend[i]];
  }

  //package nodes requested for each process
  for(i = 0; i < np; i++){
    if(i != rank){
      for(j = 0; j < commCountsSend[i]; j++){
	Int node = nodePackingList[i][j];
	globalIdSend[i][j] = (*localToGlobal)[node];
      }
    }
  }
  //send node data to processes
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(globalIdSend[i], commCountsSend[i], MPI_INT, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  //recv node data from processes directly to q vector 
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&((*localToGlobal)[nnode + commOffsetsRecv[i]]), 
		commCountsRecv[i], MPI_INT, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  //now add the ghost nodes to the global to local list... we may have to look
  //up matrices with one of their indices when using CRS
  for(i = nnode; i < nnode+gnode; i++){
    (*globalToLocal)[(*localToGlobal)[i]] = i;
  }

  //check list sanity
  for(i = 0; i < nnode; i++){
    Int gn = (*localToGlobal)[i];
    if((*globalToLocal)[gn] != i){
      std::cout << "ERROR: MAP SANITY FAILED! local: " << i << "global " << gn << std::endl;
    }
  }

  delete [] globalIdSendData;

  return globalNodes;
}

template <class Type>
Int PObj<Type>::BuildCommMaps(Mesh<Type>* m)
{
  //This routine allocates and leaves available this memory in PObj
  //--------------------------------------------------------------
  //nodePackingListData
  
  Int i, sum;
  
  this->m = m;
  
  gnode = m->gnode;
  nnode = m->nnode;
  
  MPI_Status status;
  
  for(i = 0; i < np; i++){
    commCountsRecv[i] = 0;
    commCountsSend[i] = 0;
    commOffsetsRecv[i] = 0;
    commOffsetsSend[i] = 0;
  }
  //count number of nodes to recv from each process
  for(i = 0; i < m->gnode; i++){
    commCountsRecv[m->gNodeOwner[i]]++;
  }
  
  //send number of entries expected to recv to other processes
  //not necessarily symmetric
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&commCountsRecv[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  //recv number of entries each process expects to recv from here
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&commCountsSend[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }

  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  //build utility offset mapping
  //for placing messages in correct segments
  for(i = 1; i < np; i++){
    commOffsetsRecv[i] = commCountsRecv[i-1] + commOffsetsRecv[i-1];
    commOffsetsSend[i] = commCountsSend[i-1] + commOffsetsSend[i-1];
  }

  //set up indexes for receiving
  sum = 0;
  for(i = 0; i < np; i++){
    if(i != rank){
      sum += commCountsSend[i]; 
    }
  }

  nodePackingListData = new Int[sum];

  for(i = 0; i < np; i++){
    nodePackingList[i] = &nodePackingListData[commOffsetsSend[i]];
  }

  //send list of nodes requested using ids local to owning process
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&m->gNodeLocalId[commOffsetsRecv[i]], commCountsRecv[i], MPI_INT, 
		i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  //recv list of nodes requested from each process
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&nodePackingList[i][0], commCountsSend[i], MPI_INT, 
		i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  mapsBuilt = 1;

  return 0;
}


template <class Type>
Int PObj<Type>::ReorderCommNodes(Mesh<Type>* m)
{
  //NOTE: BuildCommMaps() will not be called by this point for obvious reasons
  //this means that we can depend on none of the mesh specific internal data
  //to be correct... we have to build it here
  Int err = 0;
  Int i, j, sum;

  MPI_Status status;
  
  for(i = 0; i < np; i++){
    commCountsRecv[i] = 0;
    commCountsSend[i] = 0;
    commOffsetsRecv[i] = 0;
    commOffsetsSend[i] = 0;
  }
  //count number of nodes to recv from each process
  for(i = 0; i < m->gnode; i++){
    commCountsRecv[m->gNodeOwner[i]]++;
  }

  //send number of entries expected to recv to other processes
  //not necessarily symmetric
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&commCountsRecv[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  //recv number of entries each process expects to recv from here
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&commCountsSend[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }

  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  //build utility offset mapping
  //for placing messages in correct segments
  for(i = 1; i < np; i++){
    commOffsetsRecv[i] = commCountsRecv[i-1] + commOffsetsRecv[i-1];
    commOffsetsSend[i] = commCountsSend[i-1] + commOffsetsSend[i-1];
  }

  //set up indexes for receiving
  sum = 0;
  for(i = 0; i < np; i++){
    if(i != rank){
      sum += commCountsSend[i]; 
    }
  }
  
  nodePackingListData = new Int[sum];
  Int** orderingList = new Int*[np];
  Int* orderingListData = new Int[sum];
  //we will be recv'ing gnode new id's
  Int* recvOrdering = new Int[m->gnode];
  

  for(i = 0; i < np; i++){
    nodePackingList[i] = &nodePackingListData[commOffsetsSend[i]];
    orderingList[i] = &orderingListData[commOffsetsSend[i]];   
  }

  //send list of nodes requested for this process before remapping using ids local to owning process
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&m->gNodeLocalId[commOffsetsRecv[i]], commCountsRecv[i], MPI_INT, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  //recv list of nodes requested from each process
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&nodePackingList[i][0], commCountsSend[i], MPI_INT, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  //map requested ids to new local ordering
  for(i = 0; i < np; i++){
    if(i != rank){
      for(j = 0; j < commCountsSend[i]; j++){
  	orderingList[i][j] = m->ordering[nodePackingList[i][j]];
      }
    }
  }

  //send list of remapped nodes back to each process
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&orderingList[i][0], commCountsSend[i], MPI_INT, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }

  //receive list of remapped nodes
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&recvOrdering[commOffsetsRecv[i]], commCountsRecv[i], MPI_INT, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }

  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }
  
  //remap the gNodeLocalId list to reflect reordering on other processes
  for(i = 0; i < np; i++){
    if(i != rank){
      for(j = 0; j < commCountsRecv[i]; j++){
	m->gNodeLocalId[commOffsetsRecv[i]+j] = recvOrdering[commOffsetsRecv[i]+j];
      }
    }
  }

  delete [] nodePackingListData;
  delete [] orderingList;
  delete [] orderingListData;
  delete [] recvOrdering;

  return err;
}

template <class Type>
Int PObj<Type>::CheckSanityCoords(Type* xyz)
{
  Int i,j;
  Int err = 0;
  MPI_Status status;
  MPI_Datatype mpit;

  //set up indexes for receiving
  Int sum = 0;
  for(i = 0; i < np; i++){
    if(i != rank){
      sum += commCountsSend[i]; 
    }
  }

  //check to make sure that nodes that are being sent/recv'd are the
  //ones that we want -- use coordinates of nodes as test
  Type* recvnodes = new Type[gnode*3];
  Type* sendnodesdata = new Type[sum*3];
  Type* sendnodes[np];

  //get mpi datatype which is appropriate
  Type check = 0;
  mpit = MPI_GetType(check);

  for(i = 0; i < np; i++){
    sendnodes[i] = &sendnodesdata[commOffsetsSend[i]*3];
  }

  for(i = 0; i < np; i++){
    if(i != rank){
      for(j = 0; j < commCountsSend[i]; j++){
	Int node = nodePackingList[i][j];
	memcpy(&sendnodes[i][j*3], &xyz[node*3], sizeof(Type)*3);
      }
    }
  }

  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(sendnodes[i], commCountsSend[i]*3, mpit, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&recvnodes[commOffsetsRecv[i]*3], commCountsRecv[i]*3, mpit, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }

  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  //check that node coords are identical
  for(i = 0; i < gnode*3; i++){
    if(xyz[nnode*3 + i] != recvnodes[i]){
      Abort << "PARALLEL OBJ: Coordinates NOT matched";
      err = 1;
    }
  }

  if(err == 0){
    //std::cout << "PARALLEL OBJ: Sanity check of parallel object PASSED!!!" << std::endl;
  }
  else{
    Abort << "PARALLEL OBJ: Sanity check of parallel object FAILED!!!";
  }


  delete [] recvnodes;
  delete [] sendnodesdata;

  return err;
}


template <class Type>
Int PObj<Type>::UpdateGeneralVectors(Type* v, Int n)
{ 
  Int i,j;
  MPI_Status status;
  MPI_Datatype mpit;

  //start accumlating timer
  timers.StartAccumulate("CommTimer");
  timers.StartAccumulate("CommTotalTimer");

  //get mpi datatype to send
  Type check = 0;
  mpit = MPI_GetType(check);

  //set up indexes for receiving
  Int sum = 0;
  for(i = 0; i < np; i++){
    if(i != rank){
      sum += commCountsSend[i]; 
    }
  }

  Type* nodeVListSend[np];
  Type* nodeVListSendData = new Type[sum*n];

  for(i = 0; i < np; i++){
    nodeVListSend[i] = &nodeVListSendData[commOffsetsSend[i]*n];
  }

  //send list of nodes requested using ids local to owning process
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&m->gNodeLocalId[commOffsetsRecv[i]], commCountsRecv[i], MPI_INT, 
		i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  //recv list of nodes requested from each process
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&nodePackingList[i][0], commCountsSend[i], MPI_INT, 
		i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits for receive
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&rRequest[i], &status);
    }
  }

  //package nodes requested for each process
  for(i = 0; i < np; i++){
    if(i != rank){
      for(j = 0; j < commCountsSend[i]; j++){
	Int node = nodePackingList[i][j];
	memcpy(&nodeVListSend[i][j*n], &v[node*n], sizeof(Type)*n);
      }
    }
  }

  //clear waits for sending
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
    }
  }

  //send node data to processes
  for(i = 0; i < np; i++){
    if(i != rank && commCountsSend[i] != 0){
      MPI_Isend(nodeVListSend[i], commCountsSend[i]*n, mpit, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  //recv node data from processes directly to vector 
  for(i = 0; i < np; i++){
    if(i != rank && commCountsRecv[i] != 0){
      MPI_Irecv(&v[nnode*n + commOffsetsRecv[i]*n], commCountsRecv[i]*n, mpit, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  delete [] nodeVListSendData;

  //pause accumlating timer
  timers.PauseAccumulate("CommTimer");
  timers.PauseAccumulate("CommTotalTimer");

  return 0;
}

template <class Type>
Int PObj<Type>::GetRank()
{
  return rank;
}

template <class Type>
Int PObj<Type>::GetNp()
{
  return np;
}

template <class Type>
Int PObj<Type>::UpdateXYZ(Type* xyz)
{
  Int i,j;
  Int err = 0;
  MPI_Status status;
  MPI_Datatype mpit;

  //set up indexes for receiving
  Int sum = 0;
  for(i = 0; i < np; i++){
    if(i != rank){
      sum += commCountsSend[i]; 
    }
  }

  //check to make sure that nodes that are being sent/recv'd are the
  //ones that we want -- use coordinates of nodes as test
  Type* recvnodes = new Type[gnode*3];
  Type* sendnodesdata = new Type[sum*3];
  Type* sendnodes[np];

  //get mpi datatype which is appropriate
  Type check = 0;
  mpit = MPI_GetType(check);

  for(i = 0; i < np; i++){
    sendnodes[i] = &sendnodesdata[commOffsetsSend[i]*3];
  }

  for(i = 0; i < np; i++){
    if(i != rank){
      for(j = 0; j < commCountsSend[i]; j++){
	Int node = nodePackingList[i][j];
	memcpy(&sendnodes[i][j*3], &xyz[node*3], sizeof(Type)*3);
      }
    }
  }

  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(sendnodes[i], commCountsSend[i]*3, mpit, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&recvnodes[commOffsetsRecv[i]*3], commCountsRecv[i]*3, mpit, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }

  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  //copy sent coords over
  for(i = 0; i < gnode*3; i++){
    xyz[nnode*3 + i] = recvnodes[i];
  }

  if(err != 0){
    Abort << "UPDATE XYZ: FAILED!!!";
  }


  delete [] recvnodes;
  delete [] sendnodesdata;

  return err;
}
