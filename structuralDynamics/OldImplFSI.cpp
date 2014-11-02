  //include fluid_structure interation code
#include "../structuralDynamics/fluid_structure.h"

SectionMovement sectMove;   //stores angle and plunge for typical section analysis
  double* dxyz = NULL;        //stores explicit dxyz for full FSI integration
  int bodyId = 1;  //composite surface which is moving, for multiphysics

  std::string multPhysFileName = param.casestring + ".multphys";
  std::ofstream mpFout;

      
if((param.movement && newtonit == 1) || (param.movement && param.multiPhysics)){ 
  void* custom = NULL;
  if(param.movement == 3 && param.multiPhysics == 1){
    custom = (void*)&sectMove;
  }
  else if(param.movement == 6 && param.multiPhysics == 1){
    custom = (void*)&sectMove;
  }
  else if(param.movement == 7 && param.multiPhysics == 2){
    custom = (void*)&dxyz;
  }
  //move grid, update grid velocities, etc.
  MoveMesh(solSpaces[0], custom);
  //we need to update the points used for sampling if the mesh
  //is moving around
  solSpaces[0]->sensors->Search();
  //if we've moved the mesh we need to recompute the surface areas
  //for accurate lift/drag, etc. calculations
  ComputeSurfaceAreas(solSpaces[0], 1);
 }
  


  else if(param.multiPhysics == 2){
    Int i, j, node;
    SolidMech2->Init();
    Int* nodelist = NULL;
    Int nodecount = GetNodesOnMovingBC(m, solSpaces[0]->bc, &nodelist);
    double* xyz = new double[nodecount*3];
    for(i = 0; i < nodecount; i++){
      node = nodelist[i];
      for(j = 0; j < 3; j++){
	xyz[node*3 + j] = real(m->xyz[node*3 + j]);
      }
    }
    //build the internal connectivity rays for CFD/FEA interaction
    SolidMech2->BuildFieldConnectivity(xyz, nodecount);
    delete [] nodelist;
    delete [] xyz;
  }
  
//do the coupled physics problem solve here
// -- this is typical section airfoil analysis
if(param.multiPhysics == 1){
  //we have to compute the forces inside the Newton loop
  //if doing fluid-structure interactions, we need cl, and cm
  solSpaces[0]->forces->Compute();
  std::cout << "TYPICAL_SECTION:\tcl = " << m->bodies[bodyId].cl << "\tcm = " 
	    << m->bodies[bodyId].cm << std::endl;
  solSpaces[0]->SolidMech->BuildRHS((double)real(m->bodies[bodyId].cl), (double)real(m->bodies[bodyId].cm));
  solSpaces[0]->SolidMech->ComputeDx(&sectMove.alpha, &sectMove.h);
 }
// -- this is fully coupled FSI
 else if(param.multiPhysics == 2){
   Int* nodelist = NULL;
   Int nodecount = GetNodesOnMovingBC(m, solSpaces[0]->bc, &nodelist);
   double* forceslist = new double[nodecount*3];
   //compute the surface forces and transfer them to the FEA model
   forces->DiscreteSurfaceForces(nodelist, nodecount, forceslist);
   solSpaces[0]->SolidMech2->ApplyCFDForce(forceslist);
   //compute the FEA solution
   solSpaces[0]->SolidMech2->ComputeDx();
   //compute the displacement of the wetted surface points
   solSpaces[0]->SolidMech2->GetPointsDisplacements(dxyz);
   delete [] nodelist;
   delete [] forceslist;
 }
