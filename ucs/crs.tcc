#include "parallel.h"
#include "matrix.h"
#include "mem_util.h"
#include "crsmatrix.h"

template <class Type>
CRS<Type>::CRS()
{
  A = NULL;
  b = NULL;
  x = NULL;

  p = NULL;
}

template <class Type>
CRS<Type>::~CRS()
{
  delete A;
  delete [] b;
  delete [] x;
}

template <class Type>
//NOTE: psp map must contain nodes for parallel boundary edges
//      that is, looping over ipsp gets parallel updates included
void CRS<Type>::Init(Int nnodes_, Int gnodes_, Int neqn_, Int* ipsp, Int* psp, PObj<Type>* pobj)
{
  nnodes = nnodes_;
  gnodes = gnodes_;
  neqn = neqn_;
  neqn2 = neqn*neqn;
  p = pobj;

  A = new CRSMatrix<Type>();
  A->Init(nnodes, gnodes, neqn, ipsp, psp, pobj);

  //allocate system memory
  x = new Type[neqn*(nnodes+gnodes)];
  b = new Type[neqn*nnodes];

  return;
}

template <class Type>
Int CRS<Type>::GetSystemMemUsage(Bool print)
{
  Int mbyte = 1048576;  //since we are referring to memory
  Int tsize = sizeof(Type);
  Int isize = sizeof(Int);
  //This is for A, x, b, and pv only
  Int size = neqn*(nnodes+gnodes)*tsize + neqn*nnodes*tsize;
  size += A->GetMemUsage();
  if(print){
    std::cout << "CRS: Using " << size/mbyte << " MB for main CRS system" << std::endl;
  }
  return size;
}


template <class Type>
Type CRS<Type>::SGS(Int nSgs, CRSMatrix<Type>* A, Type* x, Type* b, Bool trueResidual)
{
  Int i,k;
  Int node, node2, diag;
  Int isgs;
  Int indx, indx1, indx2;
  Type* rhs = (Type*)alloca(sizeof(Type)*neqn);
  Type* vout = (Type*)alloca(sizeof(Type)*neqn);
  Type* temp = (Type*)alloca(sizeof(Type)*neqn);

  Type xOld = 0.0;
  Type xNorm = 0.0;

  //use internal system if another is not given
  if(A == NULL){
    A = this->A;
  }
  if(x == NULL){
    x = this->x;
  }
  if(b == NULL){
    b = this->b;
  }

  //assume initial guess x0 is passed in via the x vector we are using
  //parallel update of vectors, this allows us to take a non-zero initial guess
  p->UpdateGeneralVectors(x, neqn);

  for(isgs = 0; isgs < nSgs; isgs++){
    //forward sweep
    for(i = 0; i < nnodes; i++){
      node = i;
      memcpy(rhs, &b[node*neqn], neqn*sizeof(Type));
      indx1 = A->ia[node];
      indx2 = A->ia[node+1];
      //don't use first entry it is the diagonal of the row
      for(indx = indx1+1; indx < indx2; indx++){
	node2 = A->ja[indx];
	MatVecMult(&A->M[indx*neqn2], &x[node2*neqn], vout, neqn);
	for(k = 0; k < neqn; k++){
	  rhs[k] -= vout[k];
	}
      }
      //solve diagonal block
      diag = A->iau[node];
      if(neqn == 1){
	//if prepareSGS() was called on matrix A then the diagonal simply
	//contains the inverse, LuSolve() doesn't do the correct thing
	x[node*neqn] = A->M[diag*neqn2] * rhs[0];
      }
      else{
	LuSolve(&A->M[diag*neqn2], rhs, &A->pv[node*neqn], temp, neqn); 
	//copy solution to correct location
	memcpy(&x[node*neqn], rhs, neqn*sizeof(Type));
      }
    }
    
    //backward sweep
    for(i = nnodes-1; i >= 0; i--){
      node = i;
      memcpy(rhs, &b[node*neqn], neqn*sizeof(Type));
      indx1 = A->ia[node];
      indx2 = A->ia[node+1];
      //don't use first entry it is the diagonal of the row
      for(indx = indx1+1; indx < indx2; indx++){
	node2 = A->ja[indx];
	MatVecMult(&A->M[indx*neqn2], &x[node2*neqn], vout, neqn);
	for(k = 0; k < neqn; k++){
	  rhs[k] -= vout[k];
	}
      }
      //solve diagonal block
      diag = A->iau[node];
      if(neqn == 1){
	//if prepareSGS() was called on matrix A then the diagonal simply
	//contains the inverse, LuSolve doesn't do the correct thing
	x[node*neqn] = A->M[diag*neqn2] * rhs[0];
      }
      else{
	LuSolve(&A->M[diag*neqn2], rhs, &A->pv[node*neqn], temp, neqn); 
	//copy solution to correct location
	memcpy(&x[node*neqn], rhs, neqn*sizeof(Type));
      }
    }
    p->UpdateGeneralVectors(x, neqn);
    xOld = xNorm;
    //xNorm = ParallelL2Norm(p, x, (nnodes+gnodes)*neqn);
    xNorm = ParallelL2Norm(p, x, nnodes*neqn);
  }

  if(trueResidual){
    Type* prod = new Type[neqn*(nnodes+gnodes)];
    //performs A * x0 computation
    MatVecMultiply(A, x, prod);
    
    //find residual r0 = b - A * x0
    //set this as first search vector
    for(i = 0; i < neqn*nnodes; i++){
      prod[i] = b[i] - prod[i];
    }
    
    Type dot1_g = ParallelL2NormLittle(p, prod, neqn*nnodes);
    
    delete [] prod;

    //return the true residual
    return dot1_g;
  }

  //return the change in the x vector norm
  return (CAbs(xOld-xNorm));
}

template <class Type>
Type CRS<Type>::GMRES(Int restarts, Int nSearchDir, Int precondType, 
		      CRSMatrix<Type>* A, Type* x, Type* b)
{
  Int i, ii, j, jj;
  Int irestart, idir;
  Type dqNorm;

  //set small number for type of floating point arithmetic we are doing
  Type smallnum;
  if(typeid(Type) == typeid(double)){
    smallnum = 1.0e-15;
  }
  else if(typeid(Type) == typeid(float)){
    smallnum = 1.0e-6;
  }
  else{
    smallnum = 1.0e-15;
  }

  //use internal system if another is not given
  if(A == NULL){
    A = this->A;
  }
  if(x == NULL){
    x = this->x;
  }
  if(b == NULL){
    b = this->b;
  }

  //arnoldi search vector space
  //create doubly indexed linear array
  Type** v = new Type*[nSearchDir+1];
  Int vstride = neqn*(nnodes+gnodes);
  //allocate vectors with room for parallel syncs
  Type* vdat = new Type[(nSearchDir+1)*vstride];
  for(i = 0; i <= nSearchDir; i++){
    v[i] = &vdat[i*vstride];
  }
  
  //current search vector temp space
  Type* vtemp;
  AllocateBlankVector(&vtemp);
  //working vector
  Type* uk;
  AllocateBlankVector(&uk);
  //This is the rhs in our problem (GMRES)
  Type* g = new Type[nSearchDir+2];
 
  //Variables for solution iteration
  std::vector<Type> Q; //stores sin and cos of givens rotations
  Type cos, sin;  //temp variables
  Type a1, a2, alpha; //temp variables
  Type temp1, temp2;  //temp variables
  Type residNorm; //norm of the residual
  Type dot1_g; //global norm temp location

  //Hessenberg storage
  std::vector<Type> H;
  //stores the beginning of each col in H's index
  std::vector<Int> Hoffset;
  
  //compute preconditioner, we will be using a right preconditioner
  //That is, A*(Ninv * N)*x = b
  //         Abar = A*Ninv
  //         xbar = N*x
  //         Abar*xbar = b
  //         r0bar = b - Abar*xbar = b - A*Ninv*N*x = b - Ax = r0
  CRSMatrix<Type>* N = new CRSMatrix<Type>; 
  Preconditioner(precondType, &N, A);

  //assume initial guess x0 is passed in via the x vector we are using
  //parallel update of vectors, this allows us to take a non-zero initial guess
  p->UpdateGeneralVectors(x, neqn);

  //outer restart loop
  for(irestart = 0; irestart < restarts; irestart++){

    //performs A * x0 computation
    MatVecMultiply(A, &x[0], &v[0][0]);
    
    //find residual r0 = b - A * x0
    //set this as first search vector
    for(i = 0; i < neqn*nnodes; i++){
      v[0][i] = b[i] - v[0][i];
    }
    
    //find ||r0||
    dot1_g = ParallelL2NormLittle(p, &v[0][0], neqn*nnodes);

    //check if the residual is zero, if so, bump out of loop
    //this implies x is the solution
    if(real(dot1_g) < real(smallnum)){
      //set idir to zero to avoid segfaults :0
      idir = 0;
      break;
    }

    //normalize v0 = r0/||r0||
    for(i = 0; i < neqn*nnodes; i++){
      v[0][i] /= dot1_g;
    }
 
    //init g = ||r0||e1, i.e. rhs in our problem
    for(i = 0; i < (nSearchDir+2); i++){
      g[i] = 0.0;
    }
    g[0] = dot1_g;

    //inner search direction loop
    Int hpos = 0;
    for(idir = 0; idir < nSearchDir; idir++){
      //allocate more memory for Hess. matrix
      H.resize(hpos+(idir+2), 0.0);
      Hoffset.push_back(hpos);
      Q.resize((idir+1)*2, 0.0);
	
      //do solve for N(yk) = vk for preconditioning
      //That is, solve for yk
      PrecondBackSolve(precondType, N, vtemp, &v[idir][0]);
      
      //sync vector vk to every process
      p->UpdateGeneralVectors(vtemp, neqn);

      //init uk vector for next pass
      BlankVector(uk);

      //matrix vector multiply
      MatVecMultiply(A, vtemp, uk);
      
      //update search space
      for(j = 0; j <= idir; j++){
	dot1_g = 0.0;
	//parallel dot product
	dot1_g = ParallelDotProduct(p, uk, &v[j][0], nnodes*neqn); 
	H[hpos] = dot1_g;
	hpos++;
	//orthogonlize uk with current column
	for(i = 0; i < neqn*nnodes; i++){
	  uk[i] -= (dot1_g * v[j][i]);
	}
      }
      
      //take ||uk||_2
      dot1_g = ParallelL2NormLittle(p, uk, neqn*nnodes);
      H[hpos] = dot1_g;
      hpos++;

      //check if the residual is zero, if so, bump out of loop
      //this we can build up the solution
      if(real(dot1_g) < real(smallnum)){
	break;
      }
      
      //create next search vector
      for(i = 0; i < neqn*nnodes; i++){
	v[idir+1][i] = uk[i] / dot1_g;
      }
      
      //apply previous Givens rotations to new column in H
      for(jj = 0; jj < idir; jj++){
	cos = Q[jj*2 + 0];
	sin = Q[jj*2 + 1];
	temp1 = H[Hoffset[idir] + jj];
	temp2 = H[Hoffset[idir] + jj + 1];
	H[Hoffset[idir] + jj] = cos * temp1 + sin * temp2;
	H[Hoffset[idir] + jj + 1] = -sin * temp1 + cos * temp2;
      }	  
      //term to zero
      a2 = H[Hoffset[idir] + idir + 1];
      //term on diagonal above target
      a1 = H[Hoffset[idir] + idir];
      
      alpha = sqrt(a1*a1 + a2*a2);
      cos = a1/alpha;
      sin = a2/alpha;
      //store s & c 
      Q[idir*2 + 0] = cos;
      Q[idir*2 + 1] = sin;
      
      //apply new Givens rotation to new column in H
      H[Hoffset[idir] + idir] = alpha;
      H[Hoffset[idir] + idir + 1] = 0.0;
      
      //apply new Givens rotation to the rhs vector
      temp1 = g[idir];
      temp2 = g[idir+1];
      g[idir] = cos * temp1 + sin * temp2;
      g[idir+1] = -sin * temp1 + cos * temp2;
      
    }
    //At this point GMRES has done its thing and we are ready to 
    //compute the final solution vector
   
    //solve Ha = g
    //store a in place of g
    for(jj = idir-1; jj >=0; jj--){
      temp1 = 0.0;
      for(ii = jj+1; ii <= idir-1; ii++){
	temp1 += H[Hoffset[ii] + jj] * g[ii];
      }
      g[jj] -= temp1;
      g[jj] /= H[Hoffset[jj] + jj];
    }

    //zk = sum(j = 0, k-1) {aj * vj}
    //"build up" the solution correction
    //use vector allocation uk to store it
    for(ii = 0; ii < neqn*nnodes; ii++){
      uk[ii] = 0.0;
    }
    for(jj = 0; jj <= idir-1; jj++){
      for(ii = 0; ii < neqn*nnodes; ii++){
	//sum up zk
	uk[ii] += (v[jj][ii]*g[jj]);
      }
    }
    //solve for x_(k+1) - x_0
    //N *  (x_(k+1) - x_0) = zk
    //x_(k+1) = (x_(k+1) - x_0) + x0
    PrecondBackSolve(precondType, N, vtemp, uk);
    for(i = 0; i < neqn*nnodes; i++){
      x[i] += vtemp[i];
    }
    p->UpdateGeneralVectors(x, neqn);
  }

  //get the residual
  dqNorm = g[idir-1+1];
  
  //clean up
  delete [] v;
  delete [] vdat;
  delete [] uk;
  delete [] vtemp;
  delete [] g;
  delete N;

  return Abs(dqNorm);
}

template <class Type>
void CRS<Type>::BlankSystem()
{   
  A->Blank();
  MemBlank(x, neqn*(nnodes+gnodes));
  MemBlank(b, neqn*nnodes);

  return;
}

template <class Type>
void CRS<Type>::BlankMatrix()
{
  A->Blank();
  return;
}

template <class Type>
void CRS<Type>::BlankMatrixRow(Int row)
{
  A->BlankRow(row);
  return;
}

template <class Type>
void CRS<Type>::BlankVector(Type* v)
{
  MemBlank(v, neqn*(nnodes+gnodes));
  return;
}

template <class Type>
void CRS<Type>::BlankX()
{
  MemBlank(x, neqn*(nnodes+gnodes));
  return;
}

template <class Type>
void CRS<Type>::BlankB()
{
  MemBlank(b, neqn*nnodes);
  return;
}

template <class Type>
Int CRS<Type>::AllocateVector(Type** v)
{
  Int err = 0;
  *v = NULL;
  //this is enough space to perform a parallel update on the ghost nodes as well
  *v = new Type[neqn*(nnodes+gnodes)];
  if(*v == NULL){
    err = 1;
  }
  return err;
}

template <class Type>
Int CRS<Type>::AllocateBlankVector(Type** v)
{
  Int err = AllocateVector(v);
  MemBlank(*v, neqn*(nnodes+gnodes));
  return err;
}

template <class Type>
void CRS<Type>::MatVecMultiply(CRSMatrix<Type>* M, Type* vin, Type* vout)
{
  Int i, j;
  Int node2;
  Int indx, indx1, indx2;
  Type temp[neqn];

  //zero vout
  MemBlank(vout, neqn*(nnodes+gnodes));

  //if diagonal block is not LU decomposed, do the standard thing
  if(!M->ludiag){
    //loop over the rows
    for(i = 0; i < nnodes; i++){
      indx1 = M->ia[i];
      indx2 = M->ia[i+1];
      for(indx = indx1; indx < indx2; indx++){
	node2 = M->ja[indx];
	MatVecMult(&M->M[indx*neqn2], &vin[node2*neqn], temp, neqn);
	for(j = 0; j < neqn; j++){
	  vout[i*neqn + j] += temp[j];
	}
      }
    }
  }
  //if the diagonal block is LU decomposed, do a little more work
  else{
    //loop over the rows
    Type* diag = new Type[neqn2];
    Int dblock;
    for(i = 0; i < nnodes; i++){
      indx1 = M->ia[i];
      indx2 = M->ia[i+1];
      //skip the diagonal block
      for(indx = indx1+1; indx < indx2; indx++){
	node2 = M->ja[indx];
	MatVecMult(&M->M[indx*neqn2], &vin[node2*neqn], temp, neqn);
	for(j = 0; j < neqn; j++){
	  vout[i*neqn + j] += temp[j];
	}
      }
      //now undo the LU decomposition and compute the diagonal contribution
      dblock = M->iau[i];
      node2 = M->ja[dblock];  //this should return same value as the row id
      memcpy(diag, &M->M[dblock*neqn2], sizeof(Type)*neqn2);
      UndoLU(&M->M[dblock*neqn2], diag, &M->pv[i*neqn], neqn);
      MatVecMult(diag, &vin[node2*neqn], temp, neqn);
      for(j = 0; j < neqn; j++){
	vout[i*neqn + j] += temp[j];
      }
    }
    delete [] diag;
  }
  
  //perform a parallel sync in case we need to use vout for another parallel op
  p->UpdateGeneralVectors(vout, neqn);

  return;
}

template <class Type> template <class Type2>
void CRS<Type>::CopyMatrixSlow(CRSMatrix<Type>* Mout, CRSMatrix<Type2>* Min)
{
  Int i;
  for(i = 0; i < Min->nblocks*neqn2; i++){
    Mout[i] = Min[i];
  }
  return;
}

template <class Type>
void CRS<Type>::Preconditioner(Int type, CRSMatrix<Type>** N, CRSMatrix<Type>* A)
{
  switch(type){
  case 0:
    //This is no preconditioner
    //Don't build N just let the pointer hang out... we shouldn't use it
    break;
  case 1:
    //This is a diagonal preconditioning 
    //we just build the full block diagonal here.. though we only use the strict diagonal
    (*N)->BuildBlockDiagPrecond(A);
   break;
  case 2:
    //This is a block diagonal preconditioning
    (*N)->BuildBlockDiagPrecond(A);
    //we need to LU decomp the diagonal for backsubstitution
    (*N)->PrepareSGS();
    break;
  case 3:
    //This is local ILU(0)
    (*N)->BuildILU0Local(A);
    break;
  case 4:
    //This is SGS preconditioning
    (*N)->CopyMatrixStructure(A);
    (*N)->PrepareSGS();
    break;
  default:
    //error
    std::cerr << "WARNING: preconditioner " << type << " not available" << std::endl;
    break;
  }
  return;
}

template <class Type>
void CRS<Type>::PrecondBackSolve(Int type, CRSMatrix<Type>* N, Type* x, Type* b)
{
  //This solves the following problem
  //N*x = b
  Int i, j;
  Type* ptr;
  Type* temp = (Type*)alloca(sizeof(Type)*neqn);
  //subiterations if SGS preconditioning is used
  Int subit = 6;

  switch(type){
  case 0:
    //This is no preconditioner, simply copy over solution
    memcpy(x, b, sizeof(Type)*nnodes*neqn);
    break;
  case 1:
    //This is a diagonal preconditioning
    for(i = 0; i < nnodes; i++){
      ptr = N->GetPointer(i, i);
      for(j = 0; j < neqn; j++){
	x[i*neqn + j] = b[i*neqn + j]/ptr[j*neqn + j];
      }
    }
    break;
  case 2:
    //This is a block diagonal preconditioning
    for(i = 0; i < nnodes; i++){
      ptr = N->GetPointer(i, i);
      //copy rhs vector into the temporary memory provided by the solution
      //memory space, we do this to avoid destroying b during LuSolve() call
      //since we need it later
      memcpy(&x[i*neqn + 0], &b[i*neqn + 0], sizeof(Type)*neqn);
      LuSolve(ptr, &x[i*neqn + 0], &N->pv[i*neqn], temp, neqn);
      //solution vector is in correct location after exiting
    }
    break;
  case 3:
    //This is localized ILU0
    N->ILU0BackSub(x, b);
    break;
  case 4:
    //This is SGS preconditioning
    SGS(subit, N, x, b);
    break;
  default:
    //error
    std::cerr << "WARNING: preconditioner " << type << " not available" << std::endl;
    break;
  }
  return;
}


