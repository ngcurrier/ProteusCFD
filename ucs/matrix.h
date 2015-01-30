#ifndef MATRIX_H__
#define MATRIX_H__

#include <iostream>
#include <typeinfo>
#include <fstream>
#include <cmath>
#include "macros.h"
#include "mem_util.h"

//Perform matrix matrix multiply for
//square matrices (n x n)
template<class theType, class theType2>
void MatMatMult(const theType* a1, const theType* a2, theType* result, 
		const theType2 n){
  int i, j, k;
  theType sum;
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      sum = 0.0;
      for(k = 0; k < n; k++){
	sum += a1[i*n + k] * a2[k*n + j];
      }
      result[i*n + j] = sum;
    }
  }
  return;
}

//check for matrices being inverse of each other
//return 1 if true 0 if false
template<class theType, class theType2>
int MatInvCheck(const theType* a1, const theType* a2, const theType2 n){
  theType I[n*n];
  theType2 i, j;
  theType smallnum = 1.0e-12;
  MatMatMult(a1, a2, I, n);
  //MatPrint(I, n);
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      if(j == i){
	//check diagonal
	theType diff = 1.0 - I[i*n + i];
	if(Abs(real(diff)) > Abs(real(smallnum))){
	  std::cerr << "Diff[" << i << "]: " << diff << std::endl;
	  return (0);
	}
      }
      else{
	//check off-diagonal
	if(Abs(real(I[i*n + j])) > Abs(real(smallnum))){
	  std::cerr << "Non-zero[" << i << "," << j << "]: " << I[i*n+j] << std::endl;
	  return (0);
	}
      }
    }
  }
  return (1);
}

//Perform matrix vector multiply for
//square matrices (n x n)
template<class theType, class theType2>
void MatVecMult(const theType* __restrict__ a1, const theType* __restrict__ v1, theType* __restrict__ vout, 
		const theType2 n){
  int i,j;
  for(i = 0; i < n; i++){
    vout[i] = a1[i*n + 0] * v1[0]; 
    for(j = 1; j < n; j++){
      vout[i] += a1[i*n + j] * v1[j];
    }
  }
  return;
}

//Matrix transpose routine for a square matrix
//a - pointer to square matrix (n x n)
//n - size of matrix
//performed in place
template<class theType, class theType2>
void Transpose(theType* a, theType2 n){
  theType temp;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      temp = a[i*n + j];
      a[i*n + j] = a[j*n + i];
      a[j*n + i] = temp;
    }
  }
}

//Performs and L2 Norm calculation on vector
template<class theType, class theType2>
theType VecL2Norm(const theType *vin, const theType2 n)
{
  int i;
  theType l2norm = 0.0;
  for(i = 0; i < n; i++){
    l2norm += vin[i]*vin[i];
  }
  return (sqrt(l2norm)/(theType)n);
}

//LU decomposition routine... in place, with partial pivoting
//works on a ribbon vector style matrix i.e. a[i][j] = a[i*nj + j] 
//
//a - pointer to the square matrix (ni x ni)
//p - integer vector to hold permutations (ni)
//n - size of the matrix
template<class theType, class theType2>
int LU(theType* a, theType2* p, theType2 n){

  int i,j,k;
  int condBad = 0;
  theType large;
  int row = 0;
  int temp;
  theType smallnum;
  
  //decide on a smallnumber... this is a bit of a hack..
  //there is likely something better
  if(typeid(theType) == typeid(double)){
    smallnum = 1.0e-15;
  }
  else if(typeid(theType) == typeid(float)){
    smallnum = 1.0e-6;
  }
  else{
    smallnum = 1.0e-15;
  }
  //std::cout << smallnum << std::endl;

  //load up the default ordering of the permutations
  for(i = 0; i < n; i++){
    p[i] = i;
  }

  //loop through the matrix
  for(i = 0; i < n; i++){
    large = 0.0;
    for(j = i; j < n; j++){
      //find largest pivot in current column
      if(real(Abs(a[p[j]*n + i])) > real(Abs(large))){
	large = a[p[j]*n + i];
	row = j;
      }
    }
    //error check for a singular matrix
    if(large == 0.0){
      //std::cout << "ERROR: Singular matrix in LuRowwise" << std::endl;
      //MatPrint(a, n);
      condBad++;
    }
    else if(real(Abs(large)) < real(smallnum)){
      condBad++;
    }
    if(condBad == 1){
      //std::cout << "WARNING: Ill conditioned matrix... continuing but be warned" 
      //<< std::endl;
      //std::cout << "A Badly Conditioned by " << i << "th column: " << std::endl;
      //MatPrint(a, n);
    }
    //swap rows in permutation matrix
    temp = p[i];
    p[i] = p[row];
    p[row] = temp;
  
    //loop down the column and multiply by the new pivot's inverse
    large = 1.0/large;
    for(j = i+1; j < n; j++){
      a[p[j]*n + i] *= large; 
    }

    //loop across each subrow and subcolumn subtracting off 
    //multiplication of pivot column and pivot row
    for(j = i+1; j < n; j++){
      for(k = i+1; k < n; k++){
	a[p[j]*n + k] -= a[p[j]*n + i] * a[p[i]*n + k];
      }
    } 
  }  

  if(condBad > 0){
    std::cerr << "Ill Conditioned matrix.. continuing but be warned" << std::endl;
    //MatPrint(a,n);
    return 1;
  }

  return 0;
}

//This reverses an LU decomposition with partial pivoting so that
//we arrive back at the original matrix
//
//a - pointer to the square matrix (ni x ni)
//temp - pointer to the matrix (ni x ni) for the result
//p - integer vector to hold permutations (ni)
//n - size of the matrix
//NOTE: PA = LU
//WARNING: If the conditioning of A was bad to begin with we aren't likely
//         to get back anything meaningful
template <class theType, class theType2>
void UndoLU(const theType* a, theType* temp, const theType2* p, const theType2 n)
{
  int i, j, l;
  theType sum;

  MemBlank(temp, n*n);
  
  //loop over rows in matrix L
  for(i = 0; i < n; i++){
    //loop over columns in matrix L and rows in matrix U
    for(j = 0; j < n; j++){
      sum = 0.0;
      for(l = 0; l < i; l++){
	sum += a[p[i]*n + l]*a[p[l]*n + j];
	//no data in U past this point, break
	if(l >= j) break;
      }
      temp[p[i]*n + j] += sum;
      //add the unity contributions from the diagonal of L*U
      if(i <= j){
	temp[p[i]*n + j] += a[p[i]*n + j];
      }
    }
  }

  return;
}

//algorithm which solves a linear system given as an LU decomposition
//a - prefactored LU matrix
//b - rhs vector (will be destroyed), resultant x is stored here
//p - integer vector containing matrix permutations
//x - scratch memory needed of at least length n... can be bigger
//n - size of matrix a (n x n)
template<class theType, class theType2>
void LuSolve(theType* __restrict__ a, theType* __restrict__ b, theType2* __restrict__ p, 
	     theType* __restrict__ x, theType2 n){
  int i, j;
  theType sum;
  
  //pass one - Ly = Pb
  //sweep down
  for(i = 0; i < n; i++){
    sum = 0.0;
    for(j = 0; j < i; j++){
      sum += a[p[i]*n + j] * x[j];
    }
    x[i] = b[p[i]] - sum;
  }
  
  //second pass - Ux = y
  //sweep up
  for(i = n-1; i >= 0; i--){
    sum = 0.0;
    for(j = n-1; j > i; j--){
      sum += a[p[i]*n + j] * b[j];
    }
    b[i] = (x[i] - sum) / a[p[i]*n + i];
  }
  //at this point b holds the solution vector for the system
  return;
}

//prints out an n x n  matrix in rowwise format
//a - matrix
//n - size
template<class theType, class theType2>
void MatPrint(theType* a, theType2 n, theType2 prec = 6){
  int i,j;

  //print out matrix
  std::cout.setf(std::ios::scientific);
  //std::cout.setf(std::ios::fixed);
  std::cout.precision(prec);
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      std::cout << a[i*n + j] << "\t " ;
    }
    std::cout << std::endl;
  }
}

//prints out an n x n  matrix in rowwise format
//a - matrix
//n - size
template<class theType, class theType2>
void MatPrintStream(theType* a, theType2 n, std::ostream& str, theType2 prec = 6){
  int i,j;

  //print out matrix
  str.setf(std::ios::scientific);
  //std::cout.setf(std::ios::fixed);
  str.precision(prec);
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      str << a[i*n + j] << "\t " ;
    }
    str << std::endl;
  }
}

//prints out an m x n  matrix in rowwise format
//NS is for non-square to be non-ambiguous
//a - matrix
//m - rows
//n - columns
template<class theType, class theType2>
void MatPrintNS(theType* a, theType2 m, theType2 n, theType2 prec = 6){
  int i,j;

  //print out matrix
  std::cout.setf(std::ios::fixed);
  std::cout.precision(prec);
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      std::cout << a[i*n + j] << "\t " ;
    }
    std::cout << std::endl;
  }
}

//prints out an n x n  matrix in rowwise format
//take a permutation vector as well
//a - matrix
//n - size
template<class theType, class theType2>
void MatPrint(theType* a, theType2* p, theType2 n, theType2 prec = 6){
  
  int i,j;

  //print out matrix
  std::cout.setf(std::ios::fixed);
  std::cout.precision(prec);
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      std::cout << a[p[i]*n + j] << "\t " ;
    }
    std::cout << std::endl;
  }
}

//print out a vector of length n
//a - vector
//n - size
//prec - number of digits to print (default 6)
template<class theType, class theType2>
void VecPrint(theType* a, theType2 n, theType2 prec = 6){

  int i;
  
  std::cout.setf(std::ios::fixed);
  std::cout.precision(prec);
  for(i = 0; i < n; i++){
    std::cout << a[i] << "\t";
  }
  std::cout << std::endl;
}

//print out a permuted vector of length n
//a - vector
//p - permutation vector
//n - size
//prec - number of digits to print (default 6)
template<class theType, class theType2>
void VecPrint(theType* a, theType2* p, theType2 n, theType2 prec = 6){
  int i;
  
  std::cout.setf(std::ios::fixed);
  std::cout.precision(prec);
  for(i = 0; i < n; i++){
    std::cout << a[p[i]] << "\t";
  }
  std::cout << std::endl;
}

//write out an n x n matrix in MATLAB readable format
//filename - full filename should be appended .dat
//a - vector
//n - size
//prec - number of digits to print out (default 6)
template<class theType, class theType2, class nameType>
void MatMatlabWrite(nameType filename, theType* a, theType2 n, theType2 prec = 6){
  int i, j;
  
  std::ofstream fout(filename, std::ios::out);

  fout.setf(std::ios::fixed);
  fout.precision(prec);
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      fout << a[i*n + j] << "\t";
    }
    fout << std::endl;
  }
  fout.close();
}

//write out an m x n matrix in MATLAB readable format
//NS stands for non-square
//filename - full filename should be appended .dat
//a - vector
//m -rows
//n - columns
//prec - number of digits to print out (default 6)
template<class theType, class theType2, class nameType>
void MatMatlabWriteNS(nameType filename, theType* a, theType2 m, theType2 n, theType2 prec = 6){

  int i, j;
  
  std::ofstream fout(filename, std::ios::out);

  fout.setf(std::ios::fixed);
  fout.precision(prec);
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      fout << a[i*n + j] << "\t";
    }
    fout << std::endl;
  }
  fout.close();
}

//write out an n length vector in MATLAB readable format
//filename - full filename should be appended .dat
//a - vector
//n - size
//prec - number of digits to print out (default 6)
template<class theType, class theType2, class nameType>
void VecMatlabWrite(nameType filename, theType* a, theType2 n, theType2 prec = 6){

  int i;
  
  std::ofstream fout(filename, std::ios::out);

  fout.setf(std::ios::fixed);
  fout.precision(prec);
  for(i = 0; i < n; i++){
    fout << a[i] << "\n";
  }
  fout.close();
}


//Solve a least squares problem with Gram-Schmidt orthogonalization
//i.e. get the QR factorization
//a - matrix (will contain values of Q on return (ni x nj) )
//r - matrix (will contain values of R on return (nj x nj) )
//b - rhs vector (will contain values c1 on return)
//ni - number of rows in a
//nj - number of columns in a
//Note: QRx = b  &  Rx = Q'b  &  Rx = c1
template<class theType, class theType2>
void LsqModGS(theType* a, theType* r, theType* b, theType2 ni, theType2 nj){
  int i, j, k;
  theType norm2, norm2inv;
  theType dotprod, dotb;
  theType c1[nj];
  
  //cycle over the columns of a
  for(j = 0; j < nj; j++){
    
    //find norm 2 of the column
    norm2 = 0.0;
    for(i = 0; i < ni; i++){
      norm2 += a[i*nj + j] * a[i*nj + j];
    }
    norm2 = sqrt(norm2);
    
    //set correct value in R on diagonal
    r[j*nj + j] = norm2;
    //take the inverse of norm2
    norm2 = 1.0/norm2;
    
    //normalize the column
    for(i = 0; i < ni; i++){
      a[i*nj + j] *= norm2;
    }
    
    //orthogonalize b with new column
    dotb = 0.0;
    for(i = 0; i < ni; i++){
      dotb += a[i*nj + j] * b[i];
    }
    for(i = 0; i < ni; i++){
      b[i] -= dotb * a[i*nj + j];
    }
    //save the coefficient, but we still need all of b at this point
    c1[j] = dotb;

    for(k = j+1; k < nj; k++){
      
      //get the dot product of qT and the next column
      dotprod = 0.0;
      for(i = 0; i < ni; i++){
	dotprod += a[i*nj + j] * a[i*nj + k];
      }
      
      //save the dotproduct in our R matrix
      r[j*nj + k] = dotprod;
      
      //update the trailing columns
      for(i = 0; i < ni; i++){
	a[i*nj + k] -= dotprod * a[i*nj + j];
      }
    }
  }
  //push c1 onto b
  for(j = 0; j < nj; j++){
    b[j] = c1[j];
  }
  return;
}


//Back substitution routine
//a - matrix which is upper triangular (n x n)
//    (not strictly necessary, lower is just ignored)
//b - rhs of the system (will contain solution on exit)
//n - size of a
template<class theType, class theType2>
void BackSub(theType* a, theType* b, theType2 n){
  int i, j;
  theType sum;

  //Ax = b
  //sweep up
  for(i = n-1; i >= 0; i--){
    sum = 0.0;
    for(j = n-1; j > i; j--){
      sum += a[i*n + j] * b[j];
    }
    b[i] = (b[i] - sum) / a[i*n + i];
  }
  //at this point b holds the solution vector for the system
  return;
}


//Solve a least squares problem with Householder transformation
//i.e. get the QR factorization
//a - matrix (will contain values of R on return (ni x nj) )
//    also strictly lower triangular section contains 
//    the householder vectors for reconstruction
//b - rhs vector (will contain values c1 on return)
//aux - returns with diagonal components of Householder vectors
//      must be of minimum length nj
//ni - number of rows in a
//nj - number of columns in a
//
//TODO: implement version which does not require vector aux
//      i.e. we never want to reconstruct Q
template<class theType, class theType2>
void LsqHH(theType* a, theType* b, theType* aux, theType2 ni, theType2 nj){
  int i, j, k;
  theType alpha;
  theType sign;
  theType frac; // 2.0 * (vTu / vTv)
  theType dot1, dot2;
  //dot1 is for vTv - denominator
  //dot2 is for vTu - numerator

  //loop over the columns
  for(j = 0; j < nj; j++){
    if(a[j*nj + j] > 0.0){
      sign = -1.0;
    }
    else{
      sign = 1.0;
    }
    alpha = 0.0;
    //loop down working column to get norm2 
    //for below the diagonal
    for(i = j; i < ni; i++){
      alpha += a[i*nj + j] * a[i*nj + j];
    }
    alpha = sign * sqrt(alpha);
    
    //subtract alpha off of the diagonal
    //vector v is now in column j but
    //with terms above the diagonal
    //being non-zero.. hybrid
    a[j*nj + j] -= alpha;

    //calculate dot product vTv
    //ignore terms above j as they
    //should be zero anyway
    dot1 = 0.0;
    for(i = j; i < ni; i++){
      dot1 += a[i*nj + j] * a[i*nj + j];
    }

    //update b vector 
    dot2 = 0.0;
    for(i = j; i < ni; i++){
      dot2 += b[i] * a[i*nj + j];
    }
    frac = 2.0 * (dot2/dot1);
    for(i = j; i < ni; i++){
      b[i] -= (frac * a[i*nj + j]);
    }
   
    //loop across submatrix
    for(k = j+1; k < nj; k++){
      
      //calculate dot product vTu
      dot2 = 0.0;
      for(i = j; i < ni; i++){
	dot2 += a[i*nj + j] * a[i*nj + k];
      }
      //update column
      frac = 2.0 * (dot2/dot1);
      for(i = j; i < ni; i++){
	a[i*nj + k] -= (frac * a[i*nj + j]);
      }
    }
    //save off diagonal term of Householder vector
    aux[j] = a[j*nj + j];
    //set diagonal term explicitly
    a[j*nj + j] = alpha;
    //leave sub diagonals containing the rest of the HH vector
    //for later reconstruction of Q
  }
  return;
}

//Reconstruct the Q factor using the HH vector submatrix
//and the auxiliary vector which holds the diagonal vector entries
//
//a - matrix which has R in the upper triangular section and
//    the partial householder vectors in the strictly lower
//    triangular region (contains Q on return)
//aux - contains diagonal components of Householder vectors
//      must be of minimum length nj
//q - matrix which will contain Q on return (ni x ni)
//ni - number of rows in a (and size of q)
//nj - number of columns in a
template<class theType, class theType2>
void HHReconstructQ(theType* a, theType* aux, theType* q, theType2 ni, theType2 nj){
  int i, j, k;
  theType dot1; //vTu
  theType dot2; //vTv
  theType frac; // 2 x (vTu/vTv)

  //initialize q to the identity matrix
  for(i = 0; i < ni; i++){
    for(j = 0; j < ni; j++){
      if(i == j){
	q[i*ni + j] = 1.0;
      }
      else{
	q[i*ni + j] = 0.0;
      }
    }
  }

  //apply the HH vectors in sequence onto the identity matrix
  //Hu = u - 2((vTu)/(vTv))v
  //
  for(j = 0; j < nj; j++){
    //vTv
    dot1 = 0.0;
    dot1 += aux[j] * aux[j];
    for(i = j+1; i < ni; i++){
      dot1 += a[i*nj + j] * a[i*nj + j];
    }
    
    //apply to the matrix
    for(k = 0; k < ni; k++){
      //vTu
      dot2 = 0.0;
      dot2 += aux[j] * q[j*ni + k];
      for(i = j+1; i < ni; i++){
	dot2 += a[i*nj + j] * q[i*ni + k];
      }
      //apply to the column 
      frac = 2.0 * (dot2 / dot1);
      q[j*ni + k] -= frac * aux[j];
      for(i = j+1; i < ni; i++){
	q[i*ni + k] -= frac * a[i*nj + j];
      }
    }
  }
  //transpose matrix since Hn .. H1 = QT
  Transpose(q, ni);
  return;
}

//Cholesky decomposition A = L* L^T
//Given an SPD matrix a of size n
//On input only the upper triangle need be given, it is not modified
//The Cholesky factor L is returned in the lower triangle of a,
//except for its diagonal elements which are returned in p
template <class Type>
void Cholesky(Type* a, Type* p, int n)
{
  int i, j, k;
  Type sum;

  for(i = 0; i < n; i++){
    for(j = i; j < n; j++){
      sum = a[i*n + j];
      for(k = i-1; k >= 0; k--){
	sum -= a[i*n + k] * a[j*n + k];
      }
      if(i == j){
	if(sum <= 0.0){
	  std::cerr << "WARNING: Cholesky failed, not SPD" << std::endl;
	}
	p[i] = sqrt(sum);
      }
      else{
	a[j*n + i] = sum/p[i];
      }
    }
  }
  return;
}

//solves Ax = b given the items returned from the cholesky factorization
//routine
//a - matrix
//p - L diagonal vector
//x - solution vector
//b - right hand side
//n - size of matrix/vector
template <class Type>
void CholeskySolve(Type* a, Type* p, Type* x, Type* b, int n)
{
  int i, k;
  Type sum;

  //solve Ly = b
  //store y in x
  for(i = 0; i < n; i++){
    sum = b[i];
    for(k = i-1; k >= 0; k--){
      sum -= a[i*n + k]*x[k];
    }
    x[i] = sum/p[i];
  }
  
  //solve L^T x = y
  for(i = n-1; i >= 0; i--){
    sum = x[i];
    for(k = i+1; k < n; k++){
      sum -= a[k*n + i]*x[k];
    }
    x[i] = sum/p[i];
  }

  return;
}

//performs a matrix vector multiply with a matrix that has been decomposed
//with Cholesky routines above; Ax = b
//a - matrix
//x - vector
//n - size of matrix/vector
template <class Type>
void CholeskyMatVec(Type* a, Type* x, Type* b, int n)
{
  int i, j;

  for(i = 0; i < n; i++){
    b[i] = 0.0;
  }

  for(i = 0; i < n; i++){
    //the data is in the correct location
    for(j = i; j < n; j++){
      b[i] += a[i*n + j] * x[j];
    }
    //data is in the opposite half of matrix
    //take advantage of symmetricness
    for(j = 0; j < i; j++){
      b[i] += a[j*n + i] * x[j];
    }
  }
  return;
}

#endif
