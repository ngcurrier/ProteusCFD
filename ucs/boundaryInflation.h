#include "general.h"
#include "exceptions.h"
#include "move.h"
#include "bc.h"
#include "mesh.h"
#include <vector>

// NOTE: this set of functions are designed for taking a volume mesh that does not have
//       a boundary layer and systematically inflating it one boundary at a time
//       This follows the general methods outlined in "Unstructured Viscous Layer
//       Insertion Using Linear-Elastic Smoothing", by Karman et.al.

void GenerateBoundaryLayers(std::vector<int> boundaryFactagList, std::vector<Real> boundaryThicknesses,
			    std::vector<int> numberOfLayers, Mesh<Real>* m);
void InflateBoundary(int boundaryFactag, Real inflationDistance, Mesh<Real>* m);


// Generates new boundary layers on select surfaces given a volume mesh
// boundaryFactagList - vector of factags to inflate a boundary layer from
// boundaryThickness - vector of thickness (one per factag) to use for inflation
// numberOfLayers - vector of layer # (one per factag) to use for inflation
// m - mesh pointer to volume mesh to inflate
void GenerateBoundaryLayers(std::vector<int> boundaryFactagList, std::vector<Real> boundaryThicknesses,
			    std::vector<int> numberOfLayers, Mesh<Real>* m)
{
  //sanity checking
  if((boundaryFactagList.size() != boundaryThicknesses.size()) || (numberOfLayers.size() != boundaryThicknesses.size())){
    Abort << "GenerateBoundaryLayers: arguments not matching in length";
  }
  
  // Procedure:
  for(int i = 0; i < boundaryFactagList.size(); ++i){
    int factag = boundaryFactagList[i];
    Real dist = 0.01;
    // 1) Displace a single boundary from the list in the list using linear elastic smoothing
    // compute the next layer's distance
    InflateBoundary(factag, dist, m);
    // continue to next boundary..
  }
  // end when all layers have been inserted
  
}

// Insert a layer of prism/hexes using the surface mesh as a connecting region
// boundaryFactag - factag of boundary we'd like to inflate
// inflationDistance - distance in mesh coordinates (dimensional/non-dimensional) to inflate
// m - volume mesh to insert a layer into on boundary
void InflateBoundary(int boundaryFactag, Real inflationDistance, Mesh<Real>* m)
{
  Int iter = 10;
  Real* dx = new Real[m->GetNumNodes()*3];

  Int* pts;
  Int npts = m->FindPointsWithFactag(&pts, boundaryFactag);

  // save off the current location of these points
  Real* oldxyz = new Real[npts*3];
  m->GetCoordsForPoints(oldxyz, pts, npts);

  // TODO: check to see if any of these points are connected on multiple factags
  // If a point is connected to multiple factags that are being extruded it's distance
  // should be split between those two surfaces, if it is connected to another that is
  // not being extruded then the adjacent boundary should be treated as a symmetry surface
  // (that is) the movement normal is in the plane with the symmetry wall
  
  // compute the normals at each point using surrounding geometry
  // and then compute the total displacement at that node
  for(Int i = 0; i < npts; ++i){
    int ptid = pts[i];
    std::vector<Real> normal;
    m->GetNodeNeighborhoodNormal(ptid, normal);
    dx[ptid*3 + 0] = normal[0]*inflationDistance;
    dx[ptid*3 + 1] = normal[1]*inflationDistance;
    dx[ptid*3 + 2] = normal[2]*inflationDistance;
  }

  // find the surface elements which are on the appropriate boundary
  std::vector<Int> elementIds;
  m->FindSurfaceElementsWithFactag(elementIds, boundaryFactag);
  
  //TODO: read in BC file so we can do the appropriate thing on symmetry planes (slipping nodes, etc.)
  BoundaryConditions<Real>* bc = NULL;
  MoveMeshLinearElastic(m, bc, dx, iter);

  // append the old nodes back and use them to reset the boundary elements and create an
  // interstitial layer of volume elements
  m->AppendNodes(npts, oldxyz);
  
  // loop over all the elements on that factag
  for(Int i = 0; i < elementIds.size(); ++i){
    Int ielem = elementIds[i];
    Element<Real>* elem = m->elementList[ielem];
    Int* oldnodes;
    Int nnodes = elem->GetNnodes();
    elem->GetNodes(&oldnodes);

    // --- reset nodes on surface elements to the new nodes at old location
    Int newnodes[8];
//*********** TODO: compute the nodes in the new list
    elem->Init(newnodes);

    Int etype = elem->GetType();
    Element<Real>* tempe = NULL;
      
    // --- create new volume element to stitch the layer
    // TRI -> PRISM : QUAD -> HEX keywords
    switch(etype) {
    case TRI :
      tempe = new Prism<Real>;
      break;
    case QUAD :
      tempe = new Hexahedron<Real>;
      break;
    default:
      std::cerr << "Type not inflatable in InflateBoundary()" << std::endl;
    }
//*********** TODO: compute the nodes in the new element list
    tempe->Init(newnodes);
    tempe->SetFactag(boundaryFactag);
    m->elementList.push_back(tempe);
  }
  
  delete [] pts;
  delete [] oldxyz;

  // generate the maps required to take the next layer insertion
  m->BuildMapsDecomp();

  std::cout << "LAYER INFLATION SUCCESSFUL ON FACTAG: " << boundaryFactag << std::endl;
}


