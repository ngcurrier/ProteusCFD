#include "bc.h"
#include "solutionField.h"
#include "parallel.h"

template <class Type>
GaussianSource<Type>::GaussianSource() :
  sigma2(0.01)
{}

template <class Type>
GaussianSource<Type>::GaussianSource(Int eqnid, Int bcid, Type xloc, Type yloc, Type ampl, 
				     SolutionSpace<Type>& space):
  sigma2(0.01), xloc(xloc), yloc(yloc), ampl(ampl), eqnid(eqnid), bcid(bcid), space(space)
{}

template <class Type>
void GaussianSource<Type>::SetXY(Type xloc, Type yloc)
{
  this->xloc = xloc;
  this->yloc = yloc;
}

template <class Type>
void GaussianSource<Type>::SetAmplitude(Type ampl)
{
  this->ampl = ampl;
}

template <class Type>
void GaussianSource<Type>::ApplyToResidual()
{
  Mesh<Type>* m = space.m;
  EqnSet<Type>* eqnset = space.eqnset;
  Param<Type>* param = space.param;
  Int neqn = eqnset->neqn;

  Int velz = eqnset->GetMomentumLocation()+2;

  SolutionField<Type>& fsrc = space.GetField("gaussian");
  fsrc.Fill(0.0);
  Type* src = fsrc.GetData(FIELDS::STATE_NONE);

  std::vector<Element<Type>*> selemList;
  GetSElemsOnBCType(m, space.bc, bcid, selemList);

  SpatialFunctor<Type>* GFunc = new GaussianFunctor<Type>(*this);

  Int* nodes;
  Type elemXYZ[4*3];
  Type I = 0.0;
  for(Int nelem = 0; nelem < selemList.size(); nelem++){
    Int nnodes = selemList[nelem]->GetNodes(&nodes);
    for(Int i = 0; i < nnodes; ++i){
      Int node = nodes[i];
      elemXYZ[i*3 + 0] = m->xyz[node*3 + 0];
      elemXYZ[i*3 + 1] = m->xyz[node*3 + 1];
      elemXYZ[i*3 + 2] = m->xyz[node*3 + 2];
    }
    //our normals face outward in the XY plane, flip the sign to match with standard quadrature rules
    //here we use order 4 b/c that is as high as our quadrature rules go for triangles at this time
    Type val = -selemList[nelem]->IntegrateFunction(GFunc, elemXYZ, 4);
    I += val;
    Type contrib = val/(Type)nnodes;
    Type w[4];
    Type wsum = 0.0;
    //prevent divide by zero errors
    if(real(val) >= 1.0e-16){
      //this has problems when evaluate returns 0.0 on all nodes, the source is no longer constant
      //happens if the mesh is poorly resolved and the source is centered in the boundary element0
#if 0
      //build a weighting function assuming linear function, helps smooth source on course grids
      for(Int i = 0; i < nnodes; i++){
	w[i] = Evaluate(elemXYZ[i*3 + 0], elemXYZ[i*3 + 1], elemXYZ[i*3 + 2]);
	wsum += w[i]*w[i];
      }
      wsum = sqrt(wsum);
      for(Int i = 0; i < nnodes; i++){
	w[i] /= wsum;
      }
#endif
      for(Int i = 0; i < nnodes; i++){
	w[i] = 1.0/(Type)nnodes;
      }
      for(Int i = 0; i < nnodes; i++){
	Int node = nodes[i];
	//check that we don't try to contribute to a ghost node
	if(node < m->GetNumNodes()){
	  //we actually want to control the precise mass being added, not the density integration
	  //this is currently commented out b/c the mass source addition on viscous
	  //boundaries is unstable, all of the supporting code is useful, however, for 
	  //computing a strict inflow condition which is mass conservative when moved
	  //space.crs->b[node*neqn + eqnid] += w[i]*val;
	}
	src[node] += w[i]*val;
      }
    }
  }
  MPI_Datatype mpit;
  //get mpi datatype to send
  mpit = MPI_GetType(I);
  //we have to work around the absence of MPI_SUM for complex numbers
  Int rank, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  Type* values = new Type[np];
  values[rank] = I;
  MPI_Allgather(MPI_IN_PLACE, 1, mpit, values, 1, mpit, MPI_COMM_WORLD);
  I = 0.0;
  for(Int i = 0; i < np; i++){
    I += values[i];
  }
  delete [] values;

  std::cout << "Integrated Source: " << I << std::endl;
  Type ref_mass = param->ref_density*param->ref_length*param->ref_length*param->ref_length;
  std::cout << "IS (kg/s): " << I*ref_mass/param->ref_time << std::endl;

  //normalize this here b/c we reuse the gaussian distribution to 
  //specify slight inflow velocity at the plume location
  fsrc.Normalize();

  delete GFunc;
}

template <class Type>
Type GaussianSource<Type>::Evaluate(Type x, Type y , Type z)
{
  Type sigmax2, sigmay2;
  sigmax2 = sigmay2 = sigma2;
  
  Type dx = x - xloc;
  Type dy = y - yloc;
  Type dx2 = dx*dx;
  Type dy2 = dy*dy;
  
  Type s = ampl*exp(-1.0*(dx2/(2.0*sigmax2) + dy2/(2.0*sigmay2)));
  return s;
}

template <class Type>
GaussianFunctor<Type>::GaussianFunctor(GaussianSource<Type>& gaussianSource) :
  gaussianSource(gaussianSource)
{
}

template <class Type>
Type GaussianFunctor<Type>::Evaluate(Type x, Type y, Type z)
{
  return gaussianSource.Evaluate(x, y, z);
}
