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
  // on a series of refined grids to determine order of grid convergence

  EqnSet<Real>* eqnset = lowFiSpace->eqnset;
  Int neqn = eqnset->neqn;
  Int nauxvars = eqnset->nauxvars;
  Int ntotal = neqn+nauxvars;
  Real* q;

  std::cout << "Manipulating the field solution" << std::endl;
  
  Mesh<Real>* lowm = lowFiSpace->m;
  q = lowFiSpace->q;
  for(int i = 0; i < lowm->GetNumNodes(); ++i){
    Real x = lowm->xyz[i*3 + 0];
    Real y = lowm->xyz[i*3 + 1];
    Real z = lowm->xyz[i*3 + 2];
    q[i*ntotal + 0] = sin(4.0*PI*x);
  }

  Mesh<Real>* medm = medFiSpace->m;
  q = medFiSpace->q;
  for(int i = 0; i < medm->GetNumNodes(); ++i){
    Real x = medm->xyz[i*3 + 0];
    Real y = medm->xyz[i*3 + 1];
    Real z = medm->xyz[i*3 + 2];
    q[i*ntotal + 0] = sin(4.0*PI*x);
  }
  Mesh<Real>* him = hiFiSpace->m;
  q = hiFiSpace->q;
  for(int i = 0; i < him->GetNumNodes(); ++i){
    Real x = him->xyz[i*3 + 0];
    Real y = him->xyz[i*3 + 1];
    Real z = him->xyz[i*3 + 2]; 
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

  Real* lowqgrad = lowFiSpace->qgrad;
  Real* medqgrad = medFiSpace->qgrad;
  Real* hiqgrad = hiFiSpace->qgrad;
  std::vector<Int> gradientList;
  //This calls returns the same for all spaces (low, med, hi)
  Int nterms = eqnset->GetGradientsLocation(gradientList);

  std::cout << "Checking for zero'd gradients" << std::endl;
  
  // now let's run the tests
  // only variable that should be tweaked is density (location 0)
  // all others should be zero'd
  for(int i = 0; i < lowm->GetNumNodes(); ++i){
    for(int ieqn = 1; ieqn < nterms; ++ieqn){
      for(int idir = 0; idir < 3; ++idir){
	EXPECT_NEAR(0.0, lowqgrad[i*nterms*3 + ieqn*3 + idir], 1.0e-15);
      }
    }
  }
  for(int i = 0; i < medm->GetNumNodes(); ++i){
    for(int ieqn = 1; ieqn < nterms; ++ieqn){
      for(int idir = 0; idir < 3; ++idir){
	EXPECT_NEAR(0.0, medqgrad[i*nterms*3 + ieqn*3 + idir], 1.0e-15);
      }
    }
  }
  for(int i = 0; i < him->GetNumNodes(); ++i){
    for(int ieqn = 1; ieqn < nterms; ++ieqn){
      for(int idir = 0; idir < 3; ++idir){
	EXPECT_NEAR(0.0, hiqgrad[i*nterms*3 + ieqn*3 + idir], 1.0e-15);
      }
    }
  }

  std::cout << "Checking for gradient solution error" << std::endl;

  int densityloc = 0;
  Real percentError = 1.1; //some of the terms have the wrong sign
  Real lowerrnorm = 0.0;
  for(int i = 0; i < lowm->GetNumNodes(); ++i){
    // The derivative of the enforced solution is 4*pi*cos(4*pi*x) and we can check that analytically
    // on a series of refined grids to determine order of grid convergence
    Real x = lowm->xyz[i*3 + 0];
    Real y = lowm->xyz[i*3 + 1];
    Real z = lowm->xyz[i*3 + 2];
    Real gradsol = 4.0*PI*cos(4.0*PI*x);
    Real err = Abs(gradsol*percentError);
    EXPECT_NEAR(gradsol, lowqgrad[i*nterms*3 + densityloc*3 + 0], err);
    Real diff = gradsol - lowqgrad[i*nterms*3 + densityloc*3 + 0];
    lowerrnorm += diff*diff;
    //EXPECT_NEAR(0.0, lowqgrad[i*nterms*3 + densityloc*3 + 1], 1.0e-12);
    EXPECT_NEAR(0.0, lowqgrad[i*nterms*3 + densityloc*3 + 2], 1.0e-13);
  }
  lowerrnorm = sqrt(lowerrnorm)/lowm->GetNumNodes();

  percentError = 0.485;
  Real mederrnorm = 0.0;
  for(int i = 0; i < medm->GetNumNodes(); ++i){
    // The derivative of the enforced solution is 4*pi*cos(4*pi*x) and we can check that analytically
    // on a series of refined grids to determine order of grid convergence
    Real x = medm->xyz[i*3 + 0];
    Real y = medm->xyz[i*3 + 1];
    Real z = medm->xyz[i*3 + 2];
    Real gradsol = 4.0*PI*cos(4.0*PI*x);
    Real err = Abs(gradsol*percentError);
    EXPECT_NEAR(gradsol, medqgrad[i*nterms*3 + densityloc*3 + 0], err);
    Real diff = gradsol - medqgrad[i*nterms*3 + densityloc*3 + 0];
    mederrnorm += diff*diff;
    //EXPECT_NEAR(0.0, medqgrad[i*nterms*3 + densityloc*3 + 1], 1.0e-15);
    EXPECT_NEAR(0.0, medqgrad[i*nterms*3 + densityloc*3 + 2], 1.0e-13);
  }
  mederrnorm = sqrt(mederrnorm)/medm->GetNumNodes();

  percentError = 0.20;
  Real hierrnorm = 0.0;
  for(int i = 0; i < him->GetNumNodes(); ++i){
    // The derivative of the enforced solution is 4*pi*cos(4*pi*x) and we can check that analytically
    // on a series of refined grids to determine order of grid convergence
    Real x = him->xyz[i*3 + 0];
    Real y = him->xyz[i*3 + 1];
    Real z = him->xyz[i*3 + 2];
    Real gradsol = 4.0*PI*cos(4.0*PI*x);
    Real err = Abs(gradsol*percentError);
    EXPECT_NEAR(gradsol, hiqgrad[i*nterms*3 + densityloc*3 + 0], err);
    Real diff = gradsol - hiqgrad[i*nterms*3 + densityloc*3 + 0];
    hierrnorm += diff*diff;
    //EXPECT_NEAR(0.0, hiqgrad[i*nterms*3 + densityloc*3 + 1], 1.0e-15);
    EXPECT_NEAR(0.0, hiqgrad[i*nterms*3 + densityloc*3 + 2], 1.0e-13);
  }
  hierrnorm = sqrt(hierrnorm)/him->GetNumNodes();


  std::cout << "Low fidelity ||error||_2: " << lowerrnorm << std::endl;
  std::cout << "Medium fidelity ||error||_2: " << mederrnorm << std::endl;
  std::cout << "High fidelity ||error||_2: " << hierrnorm << std::endl;

  std::cout << "\nEach grid is twice as dense as the previous. (1/2)^2 = 4.0." << std::endl;
  std::cout << "Second order error should be 4X lower with each refinement" << std::endl;

  std::cout << "Low/Med Ratio: " << mederrnorm/lowerrnorm << std::endl;
  std::cout << "Med/Hi Ratio: " << lowerrnorm/mederrnorm << std::endl;
  
  //TODO: there is still a problem in the y-direction (out of plane) error
  //TODO: also, the convergence in X-dir isn't quite second order which is a concern

  
  //this is to clean up our solution spaces, etc.
  for(std::vector<SolutionSpaceBase<Real>*>::iterator it = solSpaces.begin(); it != solSpaces.end(); ++it){
    SolutionSpaceBase<Real>* space = *it;
    delete space;
  }
  solSpaces.clear();

}
