#include "general.h"
#include "exceptions.h"
#include "move.h"
#include "bc.h"
#include "mesh.h"
#include <vector>
#include <map>

// NOTE: this set of functions are designed for taking a volume mesh that does not have
//       a boundary layer and systematically inflating it one boundary at a time
//       This follows the general methods outlined in "Unstructured Viscous Layer
//       Insertion Using Linear-Elastic Smoothing", by Karman et.al.

void GenerateBoundaryLayers(std::vector<int> boundaryFactagList, std::vector<Real> boundaryThicknesses,
			    std::vector<int> numberOfLayers, Mesh<Real>* m, Real growthRate);
void InflateBoundary(int boundaryFactag, Real inflationDistance, Mesh<Real>* m);


// Generates new boundary layers on select surfaces given a volume mesh
// boundaryFactagList - vector of factags to inflate a boundary layer from
// boundaryThicknesses - vector of thickness (one per factag) to use for inflation
// numberOfLayers - vector of layer # (one per factag) to use for inflation
// m - mesh pointer to volume mesh to inflate
// growthRate - geometric growth rate of layers
void GenerateBoundaryLayers(std::vector<int> boundaryFactagList, std::vector<Real> boundaryThicknesses,
			    std::vector<int> numberOfLayers, Mesh<Real>* m, Real growthRate)
{
  //sanity checking
  if((boundaryFactagList.size() != boundaryThicknesses.size()) || (numberOfLayers.size() != boundaryThicknesses.size())){
    Abort << "GenerateBoundaryLayers: arguments not matching in length";
  }

  // Procedure:
  for(int i = 0; i < boundaryFactagList.size(); ++i){
    int factag = boundaryFactagList[i];
    int layers = numberOfLayers[i];
    int power = layers - 1;
    Real factor = 0.0;
    for(int j = 0; j < layers; ++j){
      factor += pow(growthRate,j);
    }
    Real dist = boundaryThicknesses[i]/factor;
    for(int j = 1; j <= numberOfLayers[i]; ++j){
      // 1) Displace a single boundary from the list in the list using linear elastic smoothing
      // compute the next layer's distance
      InflateBoundary(factag, dist*pow(growthRate,layers-j), m);
      // continue to next boundary..
    }
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

  std::cout << "MESH UTILITY: Inflating " << npts << " boundary nodes by " << inflationDistance << std::endl;
   
  if(npts == 0){
    std::cout << "WARNING: factag " << boundaryFactag << " does not seem to have points associated" << std::endl;
    return;
  }
  

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

  std::cout << "MESH UTILITY: " << elementIds.size() << " surface elements ready for extrusion " << std::endl;
  
  //TODO: read in BC file so we can do the appropriate thing on symmetry planes (slipping nodes, etc.)
  BoundaryConditions<Real>* bc = NULL;
  MoveMeshLinearElastic(m, bc, dx, iter);

  // append the old nodes back and use them to reset the boundary elements and create an
  // interstitial layer of volume elements
  Int oldnnode = m->GetNumNodes();
  m->AppendNodes(npts, oldxyz);
  // map goes from oldNodeId --> adjacent inserted node
  std::map<int,int> nodemap;
  for(Int jpt = 0; jpt < npts; ++jpt){
    Int nodeid = pts[jpt];
    //new nodes are appended at the end of the mesh
    nodemap[nodeid] = jpt+oldnnode;
  }

  std::cout << "MESH UTILITY: pre-extrusion mesh has " << m->elementList.size() << " elements" << std::endl;
  std::cout << "MESH UTILITY: pre-extrusion mesh has " << oldnnode << " nodes" << std::endl;

  // loop over all the elements on that factag
  for(Int i = 0; i < elementIds.size(); ++i){
    Int ielem = elementIds[i];
    Element<Real>* elemsurf = m->elementList[ielem];
    Int etypesurf = elemsurf->GetType();
    Int* oldsurfnodes;
    Int enodes = elemsurf->GetNnodes();
    elemsurf->GetNodes(&oldsurfnodes);

    // --- get nodes on surface elements from the new nodes list
    Int newsurfnodes[4];
    for(Int j = 0; j < enodes; ++j){
      // access the map like
      newsurfnodes[j] = nodemap.at(oldsurfnodes[j]);
    }

    // we now have, for each surface element, the old node numbers which are
    // now pushed into the volume field and the new node numbers which define the extrusion
    // at the surface stitch a new volume element from that information
    Element<Real>* tempe = NULL;
    // allocate maximum needed space for new node list for volume element
    Int newvolnodes[8];
    // --- create new volume element to stitch the layer
    // TRI -> PRISM : QUAD -> HEX keywords
    switch(etypesurf) {
    case TRI :
      tempe = new Prism<Real>;
      for(Int k = 0; k < 3; ++k){
	newvolnodes[k] = oldsurfnodes[k];
	newvolnodes[k+3] = newsurfnodes[k];
      }
      break;
    case QUAD :
      tempe = new Hexahedron<Real>;
      for(Int k = 0; k < 4; ++k){
	newvolnodes[k] = oldsurfnodes[k];
	newvolnodes[k+4] = newsurfnodes[k];
      }
      break;
    default:
      std::cerr << "Type not inflatable in InflateBoundary()" << std::endl;
    }
    tempe->Init(newvolnodes);
    tempe->SetFactag(-1);
    m->elementList.push_back(tempe);
    
    // set the old surface element to map to its new nodes on the boundary
    elemsurf->Init(newsurfnodes);
  }

  std::cout << "MESH UTILITY: extruded mesh has " << m->elementList.size() << " elements" << std::endl;
  std::cout << "MESH UTILITY: extruded mesh has " << m->GetNumNodes() << " nodes" << std::endl;

  m->UpdateElementCounts();
  delete [] pts;
  delete [] oldxyz;

  // generate the maps required to take the next layer insertion
  m->BuildMaps();
  m->CalcMetrics();
  
  std::cout << "LAYER INFLATION SUCCESSFUL ON FACTAG: " << boundaryFactag << std::endl;
}


