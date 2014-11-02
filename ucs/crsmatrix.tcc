#include "matrix.h"
#include "parallel.h"
#include "mem_util.h"
#include "matrix.h"
#include "threaded.h"
#include "exceptions.h"
#include <list>

template <class Type>
CRSMatrix<Type>::CRSMatrix()
{
  //DO NOT allocate any memory in here... That
  //should always be handled in Init() since we want to 
  //be able to control memory allocations directly, while
  //still being able to declare things for future use

  M = NULL;
  pv = NULL;

  ia = NULL;
  ja = NULL;
  iau = NULL;

  ludiag = false;
}

template <class Type>
CRSMatrix<Type>::~CRSMatrix()
{
  delete [] M;
  delete [] pv;
  delete [] ia;
  delete [] ja;
  delete [] iau;
}

template <class Type>
Int CRSMatrix<Type>::GetNnodes()
{
  return nnodes;
}

template <class Type>
Int CRSMatrix<Type>::GetGnodes()
{
  return gnodes;
}

template <class Type>
PObj<Type>* CRSMatrix<Type>::GetPObj()
{
  return p;
}

template <class Type>
void CRSMatrix<Type>::Init(Int nnodes_, Int gnodes_, Int neqn_, Int* ipsp, Int* psp, PObj<Type>* pobj)
{
  Int i;
  Int indx, indx1, indx2, count;

  nnodes = nnodes_;
  gnodes = gnodes_;
  p = pobj;
  neqn = neqn_;
  neqn2 = neqn_*neqn_;

  iau = new Int[nnodes+1];
  ia = new Int[nnodes+1];
  
  //setup ia
  ia[0] = 0;
  for(i = 0; i < nnodes; i++){
    //the entries on the row are the points connected plus the node itself
    ia[i+1] = ia[i] + (ipsp[i+1] - ipsp[i]) + 1;
  }
  nblocks = ia[nnodes];

  ja = new Int[nblocks];

  //setup ja and iau
  for(i = 0; i < nnodes; i++){
    indx1 = ipsp[i];
    indx2 = ipsp[i+1];
    //the diagonal element is the first in the CRS row
    count = ia[i];
    iau[i] = count;
    ja[count] = i;
    count++;
    //add in entries for the points connected to node i
    for(indx = indx1; indx < indx2; indx++){
      ja[count] = psp[indx];
      count++;
    }
  }

  //allocate system memory
  M = new Type[nblocks*neqn2];
  if(M == NULL){
    Abort << "WARNING in CRSMatrix<>::Init() matrix allocation failed";
  }
  pv = new Int[neqn*nnodes];
  if(pv == NULL){
    Abort << "WARNING in CRSMatrix<>::Init() permutation vector allocation failed";
  }

  return;
}

//this requires that ipsp contains all ghost nodes out to a second order stencil level
//i.e. we must have a double layer of ghost nodes in the grid otherwise the 2nd order
//stencil cannot be resolved
template <class Type>
void CRSMatrix<Type>::Init2(Int nnodes_, Int gnodes_, Int neqn_, Int* ipsp, Int* psp, PObj<Type>* pobj)
{

  nnodes = nnodes_;
  gnodes = gnodes_;
  p = pobj;
  neqn = neqn_;
  neqn2 = neqn_*neqn_;

  iau = new Int[nnodes+1];
  ia = new Int[nnodes+1];

  //The first step here is to build the psp second order stencil lists
  //This is accomplished by traversing all first neighbors in psp

  //setup ia
  ia[0] = 0;
  std::list<Int> list;
  for(Int i = 0; i < nnodes; i++){
    list.clear();
    Int indx1 = ipsp[i];
    Int indx2 = ipsp[i+1];
    for(Int indx = indx1; indx < indx2; indx++){
      Int pt = psp[indx];
      list.push_back(pt);
      //we now look at all of the points attached to node i's neighbor
      //and add them to the list, if we encounter a ghost node here
      //that is fine because we can resolve all of this locally
      Int indxd1 = ipsp[pt];
      Int indxd2 = ipsp[pt+1];
      for(Int indxd = indxd1; indxd < indxd2; indxd++){
	Int ptd = psp[indxd];
	if(ptd != i){
	  list.push_back(ptd);
	}
      }
    }
    //now we have to remove duplicates, list must be sorted first
    list.sort();
    list.unique();
    Int size = list.size();
    //the entries on the row are the points connected (2nd order) plus the node itself (diagonal)
    ia[i+1] = ia[i] + size + 1;
    //the first block is the diagonal
    iau[i] = ia[i];
  }

  nblocks = ia[nnodes];
  ja = new Int[nblocks];

  //now setup ja, same as above only load up ja this pass
  for(Int i = 0; i < nnodes; i++){
    list.clear();
    Int indx1 = ipsp[i];
    Int indx2 = ipsp[i+1];
    for(Int indx = indx1; indx < indx2; indx++){
      Int pt = psp[indx];
      list.push_back(pt);
      Int indxd1 = ipsp[pt];
      Int indxd2 = ipsp[pt+1];
      for(Int indxd = indxd1; indxd < indxd2; indxd++){
	Int ptd = psp[indxd];
	if(ptd != i){
	  list.push_back(ptd);
	}
      }
    }
    //again, remove dupes
    list.sort();
    list.unique();
    Int size = list.size();
    for(Int indx = ia[i]+1; indx < ia[i+1]; indx++){
      ja[indx] = list.front();
      list.pop_front();
    }
  }
  M = new Type[nblocks*neqn2];
  if(M == NULL){
    Abort << "WARNING in CRSMatrix<>::Init2() matrix allocation failed";
  }
  pv = new Int[neqn*nnodes];
  if(pv == NULL){
    Abort << "WARNING in CRSMatrix<>::Init() permutation vector allocation failed";
  }

  return;
}

template <class Type> 
void CRSMatrix<Type>::BuildBlockDiagPrecond(CRSMatrix<Type>* A)
{
  Int i, count;
  Type* ptr;

  nnodes = A->GetNnodes();
  gnodes = A->GetGnodes();
  neqn = A->neqn;
  neqn2 = neqn*neqn;
  p = A->GetPObj();

  iau = new Int[nnodes+1];
  ia = new Int[nnodes+1];

  //This is a block diagonal structure... therefore, ia will be mostly zeroed
  ia[0] = 0;
  for(i = 0; i < nnodes; i++){
    //the entries on the row are the points connected plus the node itself
    ia[i+1] = ia[i] + 1;
  }
  nblocks = ia[nnodes];

  ja = new Int[nblocks];
  
  //setup ja and iau
  for(i = 0; i < nnodes; i++){
    //the diagonal element is the first in the CRS row
    count = ia[i];
    iau[i] = count;
    ja[count] = i;
    //no other entries b/c there are no points connected to this node
  }

  //allocate system memory
  M = new Type[nblocks*neqn2];
  pv = new Int[neqn*nnodes];

  //now that we've allocate the memory, copy over the diagonal blocks
  for(i = 0; i < nblocks; i++){
    ptr = A->GetPointer(i, i);
    memcpy(&M[neqn2*i], ptr, sizeof(Type)*neqn2);
  }

  return;
}

template <class Type>
void CRSMatrix<Type>::CopyMatrixStructure(CRSMatrix<Type>* A)
{
  Int i, indx;

  nnodes = A->GetNnodes();
  gnodes = A->GetGnodes();
  neqn = A->neqn;
  neqn2 = neqn*neqn;
  p = A->GetPObj();

  iau = new Int[nnodes+1];
  ia = new Int[nnodes+1];

  //This is the full matrix structure
  for(i = 0; i <= nnodes; i++){
    //the entries on the row are the points connected plus the node itself
    ia[i] = A->ia[i];
  }
  nblocks = ia[nnodes];

  ja = new Int[nblocks];
  
  //copy over ja
  for(i = 0; i < nblocks; i++){
    ja[i] = A->ja[i];
  }
  
  //copy over iau
  for(i = 0; i <= nnodes; i++){
    iau[i] = A->iau[i];
  }

  //allocate system memory
  M = new Type[nblocks*neqn2];
  pv = new Int[neqn*nnodes];

  //copy over all the values in original matrix
  memcpy(M, A->M, sizeof(Type)*nblocks*neqn2);

  return;
}

template <class Type>
void CRSMatrix<Type>::BuildILU0Local(CRSMatrix<Type>* A)
{
  Int i, j, k, l;
  Int localrow, localcol; 
  Int blockrow, blockcol;
  Int localn;
  Type* diag;
  Type* block;
  Type* block2;
  Type* block3;
  Type* block4;
  Type pivot;
  
  Type smallnum;
  //decide on a smallnumber... this is a bit of a hack..
  //there is likely something better
  if(typeid(Type) == typeid(double)){
    smallnum = 1.0e-15;
  }
  else if(typeid(Type) == typeid(float)){
    smallnum = 1.0e-6;
  }
  else{
    smallnum = 1.0e-15;
  }

  CopyMatrixStructure(A);

  Type* temp = new Type[neqn2];
  
  //Assume this current matrix currently has the value of the original matrix

  //total number of rows
  Int n = nnodes*neqn;

  //loop over rows 
  for(i = 0; i < n; i++){
    //find the current row of blocks we're in
    //we know that diagonal blocks are always non-zero
    blockrow = i / neqn;
    blockcol = blockrow;
    //find the remainder of division using mod operator
    localrow = i % neqn;
    localcol = localrow;
    diag = GetPointer(blockrow, blockrow);
    //pick the pivot as the current diagonal, i.e. don't pivot
    pivot = diag[localrow*neqn + localcol];

    //check for ill-conditioning
    if(real(pivot) < real(smallnum)){
      std::cerr << "Ill Conditioned matrix.. continuing but be warned" << std::endl;
    }
    pivot = 1.0/pivot;
    //loop down column and multiply by the pivot's inverse
    //first do the portion on the block we are looking at
    for(j = localrow+1; j < neqn; j++){
      diag[j*neqn + localcol] *= pivot;
    }
    //we know that any blocks we encounter after ia[blockrow]
    //will not be above the current row, make use of that
    for(j = ia[blockrow+1]; j < nblocks; j++){
      //check for blocks directly below our pivot
      blockcol = ja[j];
      if(blockcol == blockrow){
	block = &M[j*neqn2];
	//loop over all the rows in the block
	for(k = 0; k < neqn; k++){
	  block[k*neqn + localcol] *= pivot;
	}
      }
    }
    //loop across the rest of the submatrix and subtract off the outer product
    //of the pivot column and the pivot row
    //
    //  IMPORTANT!!!
    //
    //  we know that our matrices are going to be symmetric in structure
    //  but not in value take advantage of this fact
    //
    //
    for(j = ia[blockrow]; j < ia[blockrow+1]; j++){
      //this is a block on the row which we know is nonzero
      blockcol = ja[j];
      if(blockcol >= nnodes){
	//we do not want to include ghost information
	//since this is local ILU0
	continue;
      }
      //this is the block on the pivot row
      block2 = &M[j*neqn2];
      //this is the symmetric block on the pivot column
      block3 = GetPointer(blockcol, blockrow, true);
      if(block3 != NULL){
	localn = i % neqn;
	//get the outer product for these two blocks
	//that is, product of appropriate row in block 2 and column in block 3
	for(k = 0; k < neqn; k++){
	  for(l = 0; l < neqn; l++){
	    temp[k*neqn + l] = block2[localn*neqn + l]*block3[k*neqn + localn];
	  }
	}
	//this outer product will only be subtracted from the diagonal due
	//to matrix structure
	block4 = GetPointer(blockrow, blockrow);
	for(k = 0; k < neqn; k++){
	  for(l = 0; l < neqn; l++){
	    block4[k*neqn + l] -= temp[k*neqn + l];
	  }
	}
      }

      //we also need the outer product subtracted from every block in the
      //pivot row
      for(k = localn+1; k < neqn; k++){
	for(l = 0; l < neqn; l++){
	  temp[k*neqn + l] = block2[localn*neqn + l]*diag[k*neqn + localn];
	}
      }
      //this is a block on the pivot row
      block4 = block2;
      for(k = localn+1; k < neqn; k++){
	for(l = 0; l < neqn; l++){
	  block4[k*neqn + l] -= temp[k*neqn + l];
	}
      }

      //we also need the outer product subtracted from every block in the
      //pivot column
      for(k = 0; k < neqn; k++){
	for(l = localn+1; l < neqn; l++){
	  temp[k*neqn + l] = block2[localn*neqn + l]*block3[k*neqn + localn];
	}
      }
      //this is a block on the pivot column
      block4 = block3;
      for(k = 0; k < neqn; k++){
	for(l = localn+1; l < neqn; l++){
	  block4[k*neqn + l] -= temp[k*neqn + l];
	}
      }
    }
  }

  //zero blocks which are related to parallel nodes
  for(i = 0; i < nblocks; i++){
    blockcol = ja[i];
    if(blockcol >= nnodes){
      MemBlank(&M[i*neqn2] , neqn2);
    }
  }

  delete [] temp;
  
  return;
}

template <class Type>
void CRSMatrix<Type>::ILU0BackSub(Type* x, Type* b)
{
  Int i, j, k, l;
  Type* temp = new Type[neqn];
  Type* temp2 = new Type[neqn];
  Int col;

  //NOTE: we cannot destroy b in this process

  //zero x since we will accumulate here
  BlankVector(x);

  //pass one Ly = b
  //sweep down
  for(i = 0; i < nnodes; i++){
    MemBlank(temp, neqn);
    for(j = ia[i]; j < ia[i+1]; j++){
      col = ja[j];
      //check to make sure we are below the diagonal
      if(col < i){
	MatVecMult(&M[j*neqn2], &x[col*neqn], temp2, neqn);
	for(k = 0; k < neqn; k++){
	  temp[k] += temp2[k];
	}
      }
      //check to see if we are on the diagonal
      //if we are just use the lower half of the block, i.e. L
      else if(col == i){
	for(k = 0; k < neqn; k++){
	  for(l = 0; l < k; l++){
	    temp[k] += M[j*neqn2 + k*neqn + l] * x[i*neqn + l];
	  }
	}
      }
    }
    for(k = 0; k < neqn; k++){
      x[i*neqn + k] = b[i*neqn + k] - temp[k];
    }
  }

  //second pass Ux = y
  //seep up
  for(i = nnodes-1; i >= 0; i--){
    MemBlank(temp, neqn);
    for(j = ia[i]; j < ia[i+1]; j++){
      col = ja[j];
      //check to make sure we are above the diagonal
      if(col > i){
	MatVecMult(&M[j*neqn2], &x[col*neqn], temp2, neqn);
	for(k = 0; k < neqn; k++){
	  temp[k] += temp2[k];
	}
      }
      //check to see if we are on the diagonal
      //if we are just use the upper half of the block, i.e. U
      else if(col == i){
	for(k = neqn-1; k >= 0; k--){
	  for(l = neqn-1; l > k; l--){
	    temp[k] += M[j*neqn2 + k*neqn + l] * b[i*neqn + l];
	  }
	}
      }
    }
    for(k = 0; k < neqn; k++){
      Int diag = iau[i];
      x[i*neqn + k] = (x[i*neqn + k] - temp[k])/M[diag*neqn2 + k*neqn + k];
    }
  }

  delete [] temp;
  delete [] temp2;

  return;
}


template <class Type>
void CRSMatrix<Type>::Blank()
{
  MemBlank(M, nblocks*neqn2);
  ludiag = false;
  return;
}

template <class Type>
void CRSMatrix<Type>::BlankRow(Int row)
{
  Int indx, indx1, indx2;
  indx1 = ia[row];
  indx2 = ia[row+1];
  for(indx = indx1; indx < indx2; indx++){
    MemBlank(&M[indx*neqn2], neqn2);
  }
  //what to do with ludiag?... assume the user isn't dumb 
  return;
}

template <class Type>
void CRSMatrix<Type>::BlankSubRow(Int blockrow, Int subrow)
{
  Int indx, indx1, indx2;
  Int i;
  indx1 = ia[blockrow];
  indx2 = ia[blockrow+1];
  for(indx = indx1; indx < indx2; indx++){
    for(i = 0; i < neqn; i++){
      M[indx*neqn2 + subrow*neqn + i] = 0.0;
    }
  }
  //get the diagonal block to put a one on the diagonal
  Int diag = iau[blockrow];
  M[diag*neqn2 + subrow*neqn + subrow] = 1.0;

  //what to do with ludiag?... assume the user isn't dumb 
  return;
}


template <class Type>
void CRSMatrix<Type>::PrintRow(Int row)
{
  Int indx, indx1, indx2;
  Int col;
  indx1 = ia[row];
  indx2 = ia[row+1];
  //diagonal is stored first
  for(indx = indx1; indx < indx2; indx++){
    col = ja[indx];
    std::cout << "(" << row << "," << col << ")" << std::endl;
    MatPrint(&M[indx*neqn2], neqn);
  }
  return;
}

template <class Type>
void CRSMatrix<Type>::PrintMatrix()
{
  Int i;
  for(i = 0; i < nnodes; i++){
    PrintRow(i);
  }
  return;
}

template <class Type>
void CRSMatrix<Type>::CRSTranspose()
{
  Int i, j;
  Int indx, indx1, indx2;
  Type temp[neqn2];

  //Transpose all blocks in place
  for(i = 0; i < nblocks; i++){
    Transpose(&M[i*neqn2], neqn);
  }
  //Physically tranpose blocks
  for(i = 0; i < nnodes; i++){
    indx1 = ia[i];
    indx2 = ia[i+1];
    //skip the diagonals in the first row entry.. they stay in place
    for(indx = indx1+1; indx < indx2; indx++){
      j = ja[indx];
      //since the matrix is symmetric in structure we only want to do this once
      //also, we need to leave the non-symmetric parallel pieces out of here
      //for the parallel op which is coming
      if((i < j) && j < nnodes){
	memcpy(temp, GetPointer(i,j), sizeof(Type)*neqn2);
	memcpy(GetPointer(i,j), GetPointer(j,i), sizeof(Type)*neqn2);
	memcpy(GetPointer(j,i), temp, sizeof(Type)*neqn2);
      }
    }
  }
  //Do parallel sync for the tranposed matrix blocks
  p->TransposeCommCRS(this);
  
  return;
}


template <class Type>
void CRSMatrix<Type>::TestTranspose()
{
  Int row, col, indx, indx1, indx2;
  Int rowg, colg;
  Int roffset, coffset;
  Int i, j;

  //get global indices
  Int* globalToLocal;
  Int* localToGlobal;
  Int globalCount;
  globalCount = p->BuildGlobalMaps(&localToGlobal, &globalToLocal);

  //fill the lower half of our matrix with the negative global col indx
  //fill the upper half of our matrix with the positive global row indx
  //put 999 on the absolute diagonal
  for(row = 0; row < nnodes; row++){
    rowg = localToGlobal[row];
    roffset = rowg*neqn;
    indx1 = ia[row];
    indx2 = ia[row+1];
    //do the diagonal separately, here the row == col
    indx = iau[row];
    for(i = 0; i < neqn; i++){
      for(j = 0; j < neqn; j++){
	if(i == j){
	  M[indx*neqn2 + i*neqn + j] = 999;
	}
	else if(j > i){
	  M[indx*neqn2 + i*neqn + j] = (roffset+j);
	}
	else{
	  M[indx*neqn2 + i*neqn + j] = -(roffset+i);
	}
      }
    }
    std::cout << "DIAG: " << row << std::endl;
    MatPrint(&M[indx*neqn2], neqn);
    
    //skip the diagonals in the first row entry
    for(indx = indx1+1; indx < indx2; indx++){
      col = ja[indx];
      colg = localToGlobal[col];
      coffset = colg*neqn;
      //if the col is greater than the row, this is an upper block
      if(colg > rowg){
	for(i = 0; i < neqn; i++){
	  for(j = 0; j < neqn; j++){
	    M[indx*neqn2 + i*neqn + j] = (roffset+j);
	  }
	}
      }
      else{
	for(i = 0; i < neqn; i++){
	  for(j = 0; j < neqn; j++){
	    M[indx*neqn2 + i*neqn + j] = -(coffset+i);
	  }
	}
      }
    }
  }

  //tranpose matrix
  CRSTranspose();

  //check that the expected output is there, if not throw error
  for(row = 0; row < nnodes; row++){
    rowg = localToGlobal[row];
    roffset = rowg*neqn;
    indx1 = ia[row];
    indx2 = ia[row+1];
    //do the diagonal separately, here the row == col
    indx = iau[row];
    for(i = 0; i < neqn; i++){
      for(j = 0; j < neqn; j++){
	if(i == j){
	  if(M[indx*neqn2 + i*neqn + j] != 999){
	    std::cerr << "NO MATCH IN TRANSPOSE CHECK DIAG" << std::endl;
	    std::cerr << "A: " << M[indx*neqn2 + i*neqn + j] << " vs 999" << std::endl;
	  }
	}
	else if(j > i){
	  if(M[indx*neqn2 + i*neqn + j] != -(roffset+j)){
	    std::cerr << "NO MATCH IN TRANSPOSE CHECK DIAG" << std::endl;
	    std::cerr << "A: " << M[indx*neqn2 + i*neqn + j] 
		      << " vs " << -(roffset+j) << std::endl;
	  }
	}
	else{
	  if(M[indx*neqn2 + i*neqn + j] != (roffset+i)){
	    std::cerr << "NO MATCH IN TRANSPOSE CHECK DIAG" << std::endl;
	    std::cerr << "A: " << M[indx*neqn2 + i*neqn + j] 
		      << " vs " << (roffset+i) << std::endl;

	  }
	}
      }
    }
    
    //skip the diagonals in the first row entry
    for(indx = indx1+1; indx < indx2; indx++){
      col = ja[indx];
      colg = localToGlobal[col];
      coffset = colg*neqn;
      //if the col is greater than the row, this is an upper block
      if(colg > rowg){
	for(i = 0; i < neqn; i++){
	  for(j = 0; j < neqn; j++){
	    if(M[indx*neqn2 + i*neqn + j] != -(roffset+j)){
	      std::cerr << "NO MATCH IN TRANSPOSE CHECK UPPER" << std::endl;
	      std::cerr << "A: " << M[indx*neqn2 + i*neqn + j] 
			<< " vs " << -(roffset+j) << std::endl;
	    }
	  }
	}
      }
      else{
	for(i = 0; i < neqn; i++){
	  for(j = 0; j < neqn; j++){
	    if(M[indx*neqn2 + i*neqn + j] != (coffset+i)){
	      std::cerr << "NO MATCH IN TRANSPOSE CHECK LOWER" << std::endl;
	      std::cerr << "A: " << M[indx*neqn2 + i*neqn + j] 
			<< " vs " << (coffset+i) << std::endl;

	    }
	  }
	}
      }
    }
  }

  //clean up
  delete [] localToGlobal;
  delete [] globalToLocal;
}

template <class Type>
Int CRSMatrix<Type>::GetIndex(Int row, Int col)
{
  //TODO: implement some fancy binary searching based on edge ordering...
  Int indx, indx1, indx2;

  //if we are getting the diagonal element return immediately
  //its index is stored in an adjacent list
  if(row == col){
    return iau[row];
  }

  indx1 = ia[row];
  indx2 = ia[row+1];
  //diagonal is stored first... skip checking that element due to above code
  for(indx = indx1+1; indx < indx2; indx++){
    //ja[indx] contains the column id of the block on the row at location indx
    if(ja[indx] == col){
      return indx;
    }
  }

  //return an error
  return -1;
}

template <class Type>
Type* CRSMatrix<Type>::GetPointer(Int row, Int col, Bool silent)
{
  Int indx = GetIndex(row, col);
  Type* pointer;
  if(indx != -1){
    pointer = &M[indx*neqn2];
  }
  else{
    pointer = NULL;
    if(!silent){
      std::cerr << "WARNING: ROW: " << row << " COL: " << col << " not found!" << std::endl;
      std::cerr << "WARNING: CRS pointer returned is NULL. Invalid row/column pair!" << std::endl;
    }
  }
  return (pointer);
}

template <class Type>
Int CRSMatrix<Type>::GetMaxConnectivity(Int* row)
{
  Int maxconn = 0;
  Int conn;
  Int i;
  Int indx1, indx2;

  for(i = 0; i < nnodes; i++){
    indx1 = ia[i];
    indx2 = ia[i+1];
    //skip diagonal
    conn = indx2 - (indx1+1);
    if(conn > maxconn){
      *row = i;
      maxconn = conn;
    }
  }
  return maxconn;
}

template <class Type>
Int CRSMatrix<Type>::GetGhostsInRow(Int row, Int* ghosts)
{
  Int gcount = 0;
  Int indx, indx1, indx2;
  Int col;

  if(row < nnodes){
    indx1 = ia[row];
    indx2 = ia[row+1];
    //don't look at the diagonal, it's not a ghost, trust me, BOO!
    for(indx = indx1+1; indx < indx2; indx++){
      col = ja[indx];
      if(col >= nnodes){
	ghosts[gcount] = col;
	gcount++;
      }
    }
    return gcount;
  }
  else{
    return -1;
  }
}

template <class Type>
Int CRSMatrix<Type>::GetMemUsage()
{
  Int tsize = sizeof(Type);
  Int isize = sizeof(Int);
  //This is for M, and pv only
  Int size = nblocks*neqn2*tsize + neqn*nnodes*isize;
  return size;
}

template <class Type>
void CRSMatrix<Type>::PrepareSGS()
{
  Int i;

  if(!ludiag){
#ifdef _OPENMP
    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel default(none) \
  private(i)
#endif
    {
      //compute LU decomposition for the diagonal jacobians
      if(neqn == 1){
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif    
	for(i = 0; i < nnodes; i++){
	  Type* d = GetPointer(i,i);
	  *d = 1.0/(*d);
	}
      }
      else{
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif    
	for(i = 0; i < nnodes; i++){
	  LU(GetPointer(i,i), &pv[i*neqn], neqn);
	}
      }
    }
    ludiag = true;
  }
  return;
}


//NOTE: PA = LU for diagonal blocks, we have to undo this
template <class Type>
void CRSMatrix<Type>::UndoPrepareSGS()
{
  Int i;
  Type* temp = new Type[neqn2];

  if(ludiag){
#ifdef _OPENMP
    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel default(none) shared(neqn, neqn2)	\
  private(i)
#endif
    {
      //undo the LU decomposition of the diagonal blocks
      if(neqn == 1){
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif    
	for(i = 0; i < nnodes; i++){
	  Type* d = GetPointer(i,i);
	  *d = 1.0/(*d);
	}
      }
      else{
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif    
	for(i = 0; i < nnodes; i++){
	  UndoLU(GetPointer(i,i), temp, &pv[i*neqn], neqn);
	}
      }
    }
    ludiag = false;
  }
  delete [] temp;

  return;
}

template <class Type>
void CRSMatrix<Type>::BlankVector(Type* v)
{
  MemBlank(v, neqn*(nnodes+gnodes));
  return;
}

template <class Type>
Int CRSMatrix<Type>::WriteMatrix(std::string casename)
{
  Int err = 0;
  std::ofstream fout;
  Int rank = p->GetRank();
  Int np = p->GetNp();
  std::stringstream ss;
  ss << rank;
  casename += ss.str()+".crsmat";

  err = fout.open(casename.c_str());

  //write number of processes
  fout.write((char*)&np, sizeof(Int));

  //write number of nodes
  fout.write((char*)&nnodes, sizeof(Int));

  //write number of ghost nodes
  fout.write((char*)&gnodes, sizeof(Int));

  //write number of equations
  fout.write((char*)&neqn, sizeof(Int));

  //write number of blocks
  fout.write((char*)&nblocks, sizeof(Int));

  //write ia
  fout.write((char*)ia, sizeof(Int)*(nnodes+1));

  //write iau
  fout.write((char*)iau, sizeof(Int)*(nnodes+1));

  //write ja
  fout.write((char*)ja, sizeof(Int)*(nblocks));

  //write matrix data
  fout.write((char*)M, sizeof(Type)*nblocks*neqn*neqn);

  fout.close();

  return err;
}

template <class Type>
Int CRSMatrix<Type>::ReadMatrix(std::string casename, Bool allocate)
{
  Int err = 0;
  std::ifstream fin;
  Int rank = p->GetRank();
  Int np = p->GetNp();
  Int npcheck;
  std::stringstream ss;
  ss << rank;
  casename += ss.str()+".crsmat";

  err = fin.open(casename.c_str());

  //read number of processes
  fin.read((char*)&npcheck, sizeof(Int));

  if(npcheck != np){
    std::cerr << "WARNING: Matrix being read has different number of processes than current job!!" 
	      << std::endl;
  }

  //read number of nodes
  fin.read((char*)&nnodes, sizeof(Int));

  //read number of ghost nodes
  fin.read((char*)&gnodes, sizeof(Int));

  //read number of equations
  fin.read((char*)&neqn, sizeof(Int));

  //read number of blocks
  fin.read((char*)&nblocks, sizeof(Int));

  if(allocate){
    iau = new Int[nnodes+1];
    ia = new Int[nnodes+1];
    ja = new Int[nblocks];
    M = new Type[nblocks*neqn*neqn];
    pv = new Int[neqn*nnodes];
  }

  //read ia
  fin.read((char*)ia, sizeof(Int)*(nnodes+1));

  //read iau
  fin.read((char*)iau, sizeof(Int)*(nnodes+1));

  //read ja
  fin.read((char*)ja, sizeof(Int)*(nblocks));

  //read matrix data
  fin.read((char*)M, sizeof(Type)*nblocks*neqn*neqn);

  fin.close();

  return err;
}
