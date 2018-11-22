#include <gtest/gtest.h>
#include "mesh.h"
#include "parallel.h"
#include "gradient.h"
#include "solutionSpace.h"
#include "create_functions.h"
#include <math.h>

class MeshTestGradient : public testing::Test
{
 protected:
  MeshTestGradient(){};
  ~MeshTestGradient(){};

  void SetUp(){};
  void TearDown(){};
};


TEST_F(MeshTestGradient, testLowFiMesh)
{
  std::string pathname = "../unitTest/meshResources/cubeStructuredSeries/";
  std::string casestring = "cube";
  
  PObj<Real> pobj;
  std::vector<Param<Real>* > paramList; 
  SolutionOrdering<Real> operations;
  TemporalControl<Real> temporalControl;

  ReadParamFile(paramList, casestring, pathname);
  operations.Read(casestring, pathname);
  temporalControl.Read(casestring, pathname);

  std::vector<SolutionSpaceBase<Real>*> solSpaces = CreateSolutionSpaces(paramList, pobj, temporalControl);

  //setup the solution operations (advanements, transfers, etc.)
  //TODO: this should work but it's not finding the correct # of entries
  //operations.Finalize(solSpaces);
 
  SolutionSpace<Real>* lowFiSpace = static_cast<SolutionSpace<Real>*>(solSpaces[0]);
  SolutionSpace<Real>* medFiSpace = static_cast<SolutionSpace<Real>*>(solSpaces[1]);
  SolutionSpace<Real>* hiFiSpace = static_cast<SolutionSpace<Real>*>(solSpaces[2]);
    
  
  //sanity check before proceeding
  EXPECT_NEAR(1.0, lowFiSpace->m->GetVolumeTotal(), 1.0e-14);
  EXPECT_NEAR(1.0, medFiSpace->m->GetVolumeTotal(), 1.0e-14);
  EXPECT_NEAR(1.0, hiFiSpace->m->GetVolumeTotal(), 1.0e-14);

  
  // now we set each direction x,y,z individually to sin(4*pi*x)
  // The result is that we get two full periods within the 1 x 1 x 1 grid
  // The derivative of this is 4*pi*cos(4*pi*x) and we can check that analytically
  // on a series of refined grids to determin order of grid convergence

  Mesh<Real>* m;
  EqnSet<Real>* eqnset = lowFiSpace->eqnset;
  Int neqn = eqnset->neqn;
  Int nauxvars = eqnset->nauxvars;
  Int ntotal = neqn+nauxvars;
  Real* q;

  std::cout << "Manipulating the field solution" << std::endl;
  
  m = lowFiSpace->m;
  q = lowFiSpace->q;
  for(int i = 0; i < m->GetNumNodes(); ++i){
    Real x = m->xyz[i*3 + 0];
    Real y = m->xyz[i*3 + 1];
    Real z = m->xyz[i*3 + 2];
    q[i*ntotal + 0] = sin(4.0*PI*x);
  }

  m = medFiSpace->m;
  q = medFiSpace->q;
  for(int i = 0; i < m->GetNumNodes(); ++i){
    Real x = m->xyz[i*3 + 0];
    Real y = m->xyz[i*3 + 1];
    Real z = m->xyz[i*3 + 2];
    q[i*ntotal + 0] = sin(4.0*PI*x);
  }
  m = hiFiSpace->m;
  q = hiFiSpace->q;
  for(int i = 0; i < m->GetNumNodes(); ++i){
    Real x = m->xyz[i*3 + 0];
    Real y = m->xyz[i*3 + 1];
    Real z = m->xyz[i*3 + 2]; 
    q[i*ntotal + 0] = sin(4.0*PI*x);
  }

  std::cout << "Computing gradients" << std::endl;
  
  // we'll do it manually
  lowFiSpace->grad->Compute();
  medFiSpace->grad->Compute();
  hiFiSpace->grad->Compute();

  // write out solutions modified for debugging purposes
  lowFiSpace->WriteSolution();
  medFiSpace->WriteSolution();
  hiFiSpace->WriteSolution();

  Real* qgrad = lowFiSpace->qgrad;
  std::vector<Int> gradientList;
  Int nterms = eqnset->GetGradientsLocation(gradientList);

# if 0
  for(int i = 0; i < m->GetNumNodes(); ++i){
    std::cout << qgrad[i*nterms + 0] << std::endl;
    std::cout << qgrad[i*nterms + 3] << std::endl;
    std::cout << qgrad[i*nterms + 6] << std::endl;
    std::cout << "---------" << std::endl;
  }
#endif
  
  
  //this is to clean up our solution spaces, etc.
  for(std::vector<SolutionSpaceBase<Real>*>::iterator it = solSpaces.begin(); it != solSpaces.end(); ++it){
    SolutionSpaceBase<Real>* space = *it;
    delete space;
  }
  solSpaces.clear();

}
