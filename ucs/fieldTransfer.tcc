#include "../ucs/parallel.h"
#include "interpolation.h"
#include "octree.h"

template <class Type>
FieldTransfer<Type>::FieldTransfer(SolutionSpaceBase<Type>* source, SolutionSpaceBase<Type>* dest, 
				   SolutionField<Type>& sourceField, SolutionField<Type>& destField) :
  Transfer<Type>(source, dest), sourceArrayField(sourceField), sourceXyz(NULL), 
  nptsSource(0),  destArrayField(destField), destXyz(NULL), nptsDest(0), 
  ipspSource(NULL), pspSource(NULL), geometrySet(false)
{}

template <class Type>
void FieldTransfer<Type>::SetXyz(Type* sourceXyz, Int sourceNpts, Type* destXyz, Int destNpts)
{
  this->sourceXyz = sourceXyz;
  this->destXyz = destXyz;
  geometrySet = true;
}

template <class Type>
void FieldTransfer<Type>::SetPointToPointMap(Int* ipspSource, Int* pspSource)
{
  this->pspSource = pspSource;
  this->ipspSource = ipspSource;
}

template <class Type>
void FieldTransfer<Type>::DoTransfer()
{
  if(!geometrySet){
    std::cerr << "WARNING: FieldTransfer::DoTransfer() - transfer not fully defined, refusing to continue" 
	      << std::endl;
    return;
  }

  Int t2 = 0;
  MPI_Datatype mpitint = MPI_GetType(t2);

  //power parameter for inverse distance weighted interpolation
  Type pow = 1.0;

  //this is an assumption that we use the most recent data available
  //for transfer
  Type* sourceq= sourceArrayField.GetData(FIELDS::STATE_NONE);
  Type* destq = destArrayField.GetData(FIELDS::STATE_NONE);
  Int ndof = sourceArrayField.GetNdof();

  PObj<Type>* p = sourceArrayField.GetPObj();
  Int myrank = p->GetRank();
  Int np = p->GetNp();

  //step one: Build an octree around the source
  Int maxdepth = 10;
  Int minNodesPerLevel = 10;
  Octree<Type> octree(sourceXyz, nptsSource, maxdepth, minNodesPerLevel);
  
  if(destArrayField.GetNdof() != ndof){
    std::cerr << "WARNING: In Transfer::DofFieldTransfer() source and destination " 
	      << "degrees of freedom do not match -- failing" << std::endl;
    return;
  }

  //loop over all the points in the destination list and find the
  //nearest source node to each of these
  std::vector<TypeInt<Type> > nearestDist;
  nearestDist.reserve(nptsDest);
  for(Int i = 0; i < nptsDest; i++){
    Int nodeId;
    Type dist = octree.FindDistance(&destXyz[i*3 + 0], nodeId);
    TypeInt<Type> temp;
    temp.dist = dist;
    temp.rank = myrank;
    temp.nodeId = nodeId;
    nearestDist.push_back(temp);
  }

  //  THIS IS THE PARALLEL INTERPOLATION PROCESS
  //  --------------------------------------------------------------
  // -Communicate to all other nodes which might contain the point
  //  in question due to their bounding boxes
  // -If another proc contains the node in question, check distance to
  //  the nearest node, if closer than the current distance, correct
  //  that list
  // -Look at all the nodes and their closest points(procs), sync any
  //  necessary data across the processors
  // -Finally do the interpolating
  
  //sync everywhere the bounding box edges for each process
  Type box[6];
  Type* boxExtents = new Type[np*6];
  octree.GetExtents(box);
  MPI_Datatype mpit = MPI_GetType(box[0]);
  MPI_Alltoall(box, 6, mpit, boxExtents, 6, mpit, MPI_COMM_WORLD);

  //find candidate processors which might contain each point using simple extents, this
  //might fail if the box is aligned on a boundary which falls in a simple coordinate
  //direction, however, ghost node overlap should take care of correcting the interpolation
  Int* candidateNodes[np];
  Int count[np];
  for(Int ip = 0; ip < np; ip++){
    //count points for allocation
    count[ip] = 0;
    for(Int i = 0; i < nptsDest; i++){
      if((real(destXyz[i*3 + 0]) < real(boxExtents[ip*6 + 1])) && (real(destXyz[i*3 + 0]) > real(boxExtents[ip*6 + 0]))){
	if((real(destXyz[i*3 + 1]) < real(boxExtents[ip*6 + 3])) && (real(destXyz[i*3 + 1]) > real(boxExtents[ip*6 + 2]))){
	  if((real(destXyz[i*3 + 2]) < real(boxExtents[ip*6 + 5])) && (real(destXyz[i*3 + 2]) > real(boxExtents[ip*6 + 4]))){
	    count[ip]++;
	  }
	}
      }
    }
  }
  Int commCountsSend[np];
  Int commCountsRecv[np];
  for(Int ip = 0; ip < np; ip++){
    commCountsSend[ip] = count[ip];
  }
  //allocate and store the points
  for(Int ip = 0; ip < np; ip++){
    candidateNodes[ip] = new Int[commCountsSend[ip]];
    count[ip] = 0;
    for(Int i = 0; i < nptsDest; i++){
      if((real(destXyz[i*3 + 0]) < real(boxExtents[ip*6 + 1])) && (real(destXyz[i*3 + 0]) > real(boxExtents[ip*6 + 0]))){
	if((real(destXyz[i*3 + 1]) < real(boxExtents[ip*6 + 3])) && (real(destXyz[i*3 + 1]) > real(boxExtents[ip*6 + 2]))){
	  if((real(destXyz[i*3 + 2]) < real(boxExtents[ip*6 + 5])) && (real(destXyz[i*3 + 2]) > real(boxExtents[ip*6 + 4]))){
	    candidateNodes[ip][count[ip]] = i;
	    count[ip]++;
	  }
	}
      }
    }
  }
  MPI_Status status;
  MPI_Request* sRequest = new MPI_Request[np];
  MPI_Request* rRequest = new MPI_Request[np];

  //communicate the number of entries expected to recv to other processes
  //               ----   not necessarily symmetric -----
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Isend(&commCountsSend[i], 1, mpitint, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Irecv(&commCountsRecv[i], 1, mpitint, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  Type* sendXyzBuff[np];
  Type* servicedNodesXyz[np];
  for(Int ip = 0; ip < myrank; ++ip){
    if(ip != myrank){
      sendXyzBuff[ip] = new Type[commCountsSend[ip]*3];
      servicedNodesXyz[ip] = new Type[commCountsRecv[ip]*3];
    }
  }  
  //communicate nodes xyz coords to check distances on to other processes
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      for(Int j = 0; j < commCountsSend[i]; j++){
	for(Int k = 0; k < 3; k++){
	  sendXyzBuff[i][j*3 + k] = destXyz[candidateNodes[i][j]*3 + k];
	}
      }
      MPI_Isend(sendXyzBuff[i], commCountsSend[i]*3, mpit, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Irecv(&servicedNodesXyz[i][0], commCountsRecv[i]*3, mpit, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
      delete [] sendXyzBuff[i];
    }
  }

  Type* servicedNodesDist[np];
  Int* servicedNodesNodeId[np];
  Type* minDist[np];
  Int* minDistNodeId[np];
  for(Int ip = 0; ip < myrank; ++ip){
    if(ip != myrank){
      servicedNodesDist[ip] = new Type[commCountsRecv[ip]];
      servicedNodesNodeId[ip] = new Int[commCountsRecv[ip]];
      minDist[ip] = new Type[commCountsSend[ip]];
      minDistNodeId[ip] = new Int[commCountsSend[ip]];
    }
  }
  //check the nearest distance for each point requested from other processes
  //store the minimum distance and the local node id for communication
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      for(Int j = 0; j < commCountsRecv[i]; j++){
	Int nodeId;
	servicedNodesDist[i][j] = octree.FindDistance(&servicedNodesXyz[i][j*3 + 0], nodeId);
	servicedNodesNodeId[i][j] = nodeId;
      }
      delete [] servicedNodesXyz[i];
    }
  }

  //communicate the minimum distance to home process for each requested node
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Isend(servicedNodesDist[i], commCountsRecv[i], mpit, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Irecv(minDist[i], commCountsSend[i], mpit, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
      delete [] servicedNodesDist[i];
    }
  }

  //communicate the local id of the minimum distance node to the home process for each requested node
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Isend(servicedNodesNodeId[i], commCountsRecv[i], mpitint, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Irecv(minDistNodeId[i], commCountsSend[i], mpitint, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
      delete [] servicedNodesNodeId[i];
    }
  }

  //parse the list and find either the local node id of minimum distance of the local node id
  //on a non-local process with the minimum distance
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      for(Int j = 0; j < commCountsSend[i]; j++){
	Int localNode = candidateNodes[i][j];
	if(real(minDist[i][j]) < real(nearestDist[localNode].dist)){
	  nearestDist[localNode].dist = minDist[i][j];
	  nearestDist[localNode].rank = i;
	  nearestDist[localNode].nodeId = minDistNodeId[i][j];
	}
      }
      delete [] candidateNodes[i];
      delete [] minDist[i];
      delete [] minDistNodeId[i];
    }
  }

  std::vector<Int> unresolvedPoints;
  std::vector<Int> stencil;
  for(UInt i = 0; i < nptsDest; ++i){
    if(nearestDist[i].rank != myrank){
      //if the closest match is non-local, put it on the list to resolve
      //via parallel methods
      unresolvedPoints.push_back(i);
    }
    else{
      //point is local, do the interpolation
      stencil.clear();
      Int nearestNode = nearestDist[i].nodeId;
      stencil.push_back(nearestNode);
      for(Int j = ipspSource[nearestNode]; j < ipspSource[nearestNode + 1]; j++){
	stencil.push_back(pspSource[j]);
      }
      InverseWeighting((Int)stencil.size(), stencil.data(), sourceXyz, sourceq, 
		       &destXyz[i*3 + 0], &destq[i*ndof], ndof, ndof, pow);
    }
  }

  //parse through the list of uresolved points and separate them out according to rank
  for(Int i = 0; i < np; i++){
    commCountsSend[i] = 0;
    commCountsRecv[i] = 0;
  }
  for(UInt i = 0; i < unresolvedPoints.size(); i++){
    Int pt = unresolvedPoints[i];
    Int rank = nearestDist[pt].rank;
    commCountsSend[i]++;
  }
  //communicate the number of entries expected to recv to other processes
  //               ----   not necessarily symmetric -----
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Isend(&commCountsSend[i], 1, mpitint, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Irecv(&commCountsRecv[i], 1, mpitint, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }
  
  Int* servicedNodesSendBuff[np];
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      sendXyzBuff[i] = new Type[commCountsSend[i]*3];
      servicedNodesXyz[i] = new Type[commCountsRecv[i]*3];
      servicedNodesSendBuff[i] = new Int[commCountsSend[i]];
      servicedNodesNodeId[i] = new Int[commCountsRecv[i]];
    }
  }

  //load up all arrays with request lists and xyz coords. for interpolation
  for(Int i = 0; i < np; i++){
    count[i] = 0;
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      for(UInt j = 0; j < unresolvedPoints.size(); j++){
	Int pt = unresolvedPoints[j];
	Int rank = nearestDist[pt].rank;
	servicedNodesSendBuff[rank][count[rank]] = nearestDist[pt].nodeId;
	for(Int k = 0; k < 3; k++){
	  sendXyzBuff[rank][count[rank]*3 + k] = destXyz[j*3 + k];
	}
	count[rank]++;
      }
    }
  }

  //request interpolation for all points which have a non-local minimum distance on another rank
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Isend(servicedNodesSendBuff[i], commCountsSend[i], mpitint, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Irecv(servicedNodesNodeId[i], commCountsRecv[i], mpitint, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
      delete [] servicedNodesSendBuff[i];
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Isend(sendXyzBuff[i], commCountsSend[i]*3, mpit, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Irecv(servicedNodesXyz[i], commCountsRecv[i]*3, mpit, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
      delete [] sendXyzBuff[i];
    }
  }
  
  //do interpolation as a service for all non-local points whose minimum distance host node is here
  Type* qSendBuff[np];
  Type* qRecvBuff[np];
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      qSendBuff[i] = new Type[ndof*commCountsRecv[i]];
      qRecvBuff[i] = new Type[ndof*commCountsSend[i]];
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      for(Int j = 0; j < commCountsRecv[i]; j++){
	stencil.clear();
	Int nearestNode = nearestDist[i].nodeId;
	stencil.push_back(nearestNode);
	for(Int k = ipspSource[nearestNode]; k < ipspSource[nearestNode + 1]; k++){
	  stencil.push_back(pspSource[k]);
	}
	InverseWeighting((Int)stencil.size(), stencil.data(), sourceXyz, sourceq, 
			 &servicedNodesXyz[i][j*3 + 0], &qSendBuff[i][j*ndof], ndof, ndof, pow);
      }
    }
  }

  //communicate the interpolation back to the owning process
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Isend(qSendBuff[i], commCountsRecv[i]*ndof, mpit, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Irecv(qRecvBuff[i], commCountsSend[i]*ndof, mpit, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
      delete [] qSendBuff[i];
      delete [] servicedNodesXyz[i];
    }
  }

  //Put the received interpolated values back in the destination solution field array 
  for(Int i = 0; i < np; i++){
    count[i] = 0;
  }
  for(Int i = 0; i < np; i++){
    if(i != myrank){
      for(UInt j = 0; j < unresolvedPoints.size(); j++){
	Int pt = unresolvedPoints[j];
	Int rank = nearestDist[pt].rank;
	for(Int k = 0; k < ndof; k++){
	 destq[pt*ndof + k] = qRecvBuff[rank][count[rank]*ndof + k];
	}
	count[rank]++;
      }
      delete [] qRecvBuff[i];
    }
  }


  delete [] sRequest;
  delete [] sRequest;
  delete [] boxExtents;
}
