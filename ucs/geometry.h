#ifndef GEOMETRY_H__
#define GEOMETRY_H__

#include <cmath>
#include <iostream>
#include "general.h"
#include "macros.h"
#include "etypes.h"

template <class theType>
// subtracts pt2 from pt1 -> r
void Subtract(const theType pt1[3], const theType pt2[3], theType r[3])
{
  r[0] = pt1[0] - pt2[0];
  r[1] = pt1[1] - pt2[1];
  r[2] = pt1[2] - pt2[2];
}

template <class theType>
void Add(const theType pt1[3], const theType pt2[3], theType r[3])
{
  r[0] = pt1[0] + pt2[0];
  r[1] = pt1[1] + pt2[1];
  r[2] = pt1[2] + pt2[2];
}

//when we are calling distance, we are almost always talking about the
//distance magnitude between the real parts of the geometry
template <class theType>
Real Distance(const theType pt1[3], const theType pt2[3])
{
  Real dx = real(pt1[0] - pt2[0]);
  Real dy = real(pt1[1] - pt2[1]);
  Real dz = real(pt1[2] - pt2[2]);
  return( sqrt(dx*dx + dy*dy + dz*dz) );
}

//computes the magnitude of a vector
//returns - magnitude
template <class theType>
theType Magnitude(const theType v[3])
{
  return( sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) );
}

//normalizes a vector
//v - 3 space vector
//rv - normalize vector returned
//returns - magnitude of normalization
template <class theType>
theType Normalize(const theType v[3], theType rv[3])
{
  theType mag = Magnitude(v);
  rv[0] = v[0]/mag;
  rv[1] = v[1]/mag;
  rv[2] = v[2]/mag;
  return mag;
}

template <class theType>
theType DotProduct(const theType v1[3], const theType v2[3])
{
  return( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

template <class theType>
//computes v1 x v2 = rv
void CrossProduct(const theType v1[3], const theType v2[3], theType rv[3])
{
  rv[0] = v1[1]*v2[2] - v2[1]*v1[2]; 
  rv[1] = v1[2]*v2[0] - v2[2]*v1[0];
  rv[2] = v1[0]*v2[1] - v2[0]*v1[1];
}

//vector scalar multiplication
//ni - entries in vector
template <class theType>
void Scale(const theType s, theType r[], int ni)
{
  Int i;
  for(i = 0; i < ni; i++){
    r[i] = s * r[i];
  }
}

//v should be normalized
//removes components of v1 in direction of norv
template <class theType>
void RemoveNormalComponent(const theType v1[3], const theType norv[3], theType rv[3])
{
  theType dot = DotProduct(v1,norv);
  
  rv[0] = v1[0] - dot*norv[0];
  rv[1] = v1[1] - dot*norv[1];
  rv[2] = v1[2] - dot*norv[2];
}


//require vectors to be mutually orthonormal i.e. n = v1 x v2
template <class theType>
void PerpVectors(const theType n[3], theType v1[3], theType v2[3])
{
  theType dot;

  v1[0] = v1[1] = v1[2] = 0.0;

  if (real(CAbs(dot = n[0])) < 0.95)
    {
      v1[0] = 1.0;
    }
  else if (real(CAbs(dot = n[1])) < 0.95)
    {
      v1[1] = 1.0;
    }
  else 
    {
      dot = n[2];
      v1[2] = 1.0;
    }

  //subtract off portion of v1 lying in direction n
  v1[0] -= dot*n[0];
  v1[1] -= dot*n[1];
  v1[2] -= dot*n[2];
  Normalize(v1, v1);

  CrossProduct(n,v1,v2);
  Normalize(v2, v2);
}

//centroid for beam
template <class theType>
void Centroid(const theType pt1[3], const theType pt2[3], theType rpt[3])
{
  rpt[0] = (pt1[0] + pt2[0]);
  rpt[1] = (pt1[1] + pt2[1]);
  rpt[2] = (pt1[2] + pt2[2]);
  rpt[0] *= 0.5;
  rpt[1] *= 0.5;
  rpt[2] *= 0.5;
}

//centroid for triangle
template <class theType>
void Centroid(const theType pt1[3], const theType pt2[3], const theType pt3[3], 
		     theType rpt[3])
{
  theType frac = 1.0/3.0;
  rpt[0] = (pt1[0] + pt2[0] + pt3[0]);
  rpt[1] = (pt1[1] + pt2[1] + pt3[1]);
  rpt[2] = (pt1[2] + pt2[2] + pt3[2]);
  rpt[0] *= frac;
  rpt[1] *= frac;
  rpt[2] *= frac;
}


//centroid for quads/tets
template <class theType>
void Centroid(const theType pt1[3], const theType pt2[3], const theType pt3[3],
		     const theType pt4[3], theType rpt[3])
{
  rpt[0] = (pt1[0] + pt2[0] + pt3[0] + pt4[0]);
  rpt[1] = (pt1[1] + pt2[1] + pt3[1] + pt4[1]);
  rpt[2] = (pt1[2] + pt2[2] + pt3[2] + pt4[2]);
  rpt[0] *= 0.25;
  rpt[1] *= 0.25;
  rpt[2] *= 0.25;
}

//centroid for pyramid
template <class theType>
void Centroid(const theType pt1[3], const theType pt2[3], const theType pt3[3],
		     const theType pt4[3], const theType pt5[3], theType rpt[3])
{
  rpt[0] = (pt1[0] + pt2[0] + pt3[0] + pt4[0] + pt5[0]);
  rpt[1] = (pt1[1] + pt2[1] + pt3[1] + pt4[1] + pt5[1]);
  rpt[2] = (pt1[2] + pt2[2] + pt3[2] + pt4[2] + pt5[2]);
  rpt[0] *= 0.2;
  rpt[1] *= 0.2;
  rpt[2] *= 0.2;
}

//centroid for prism
template <class theType>
void Centroid(const theType pt1[3], const theType pt2[3], const theType pt3[3],
		     const theType pt4[3], const theType pt5[3], const theType pt6[3], 
		     theType rpt[3])
{
  theType frac = 1.0/6.0;
  rpt[0] = (pt1[0] + pt2[0] + pt3[0] + pt4[0] + pt5[0] + pt6[0]);
  rpt[1] = (pt1[1] + pt2[1] + pt3[1] + pt4[1] + pt5[1] + pt6[1]);
  rpt[2] = (pt1[2] + pt2[2] + pt3[2] + pt4[2] + pt5[2] + pt6[2]);
  rpt[0] *= frac;
  rpt[1] *= frac;
  rpt[2] *= frac;
}

//centroid for hexahedral elements
template <class theType>
void Centroid(const theType pt1[3], const theType pt2[3], const theType pt3[3],
		     const theType pt4[3], const theType pt5[3], const theType pt6[3],
		     const theType pt7[3], const theType pt8[3], theType rpt[3])
{
  theType frac = 1.0/8.0;
  rpt[0] = (pt1[0] + pt2[0] + pt3[0] + pt4[0] + pt5[0] + pt6[0] + pt7[0] + pt8[0]);
  rpt[1] = (pt1[1] + pt2[1] + pt3[1] + pt4[1] + pt5[1] + pt6[1] + pt7[1] + pt8[1]);
  rpt[2] = (pt1[2] + pt2[2] + pt3[2] + pt4[2] + pt5[2] + pt6[2] + pt7[2] + pt8[2]);
  rpt[0] *= frac;
  rpt[1] *= frac;
  rpt[2] *= frac;
}

//compute centroid based on element type id
//nodes is a list containing id's for nodes in element
template <class theType, class theType2>
void ElementCentroid(const theType* nodes, const theType2 *xyz, theType etype, 
		     theType2* centroid)
{
  switch(etype){
  case TRI:
    Centroid(&xyz[nodes[0]*3], &xyz[nodes[1]*3], &xyz[nodes[2]*3], centroid);
    break;
  case QUAD:
    Centroid(&xyz[nodes[0]*3], &xyz[nodes[1]*3], &xyz[nodes[2]*3], &xyz[nodes[3]*3], 
	     centroid);
    break;
  case TET:
    Centroid(&xyz[nodes[0]*3], &xyz[nodes[1]*3], &xyz[nodes[2]*3], &xyz[nodes[3]*3], 
	     centroid);
    break;
  case PYRAMID:
    Centroid(&xyz[nodes[0]*3], &xyz[nodes[1]*3], &xyz[nodes[2]*3], &xyz[nodes[3]*3], 
	     &xyz[nodes[4]*3], centroid);
    break;
  case PRISM:
    Centroid(&xyz[nodes[0]*3], &xyz[nodes[1]*3], &xyz[nodes[2]*3], &xyz[nodes[3]*3], 
	     &xyz[nodes[4]*3], &xyz[nodes[5]*3], centroid);
    break;
  case HEX:
    Centroid(&xyz[nodes[0]*3], &xyz[nodes[1]*3], &xyz[nodes[2]*3], &xyz[nodes[3]*3], 
	     &xyz[nodes[4]*3], &xyz[nodes[5]*3], &xyz[nodes[6]*3], &xyz[nodes[7]*3], 
	     centroid);
    break;
  default:
    std::cerr << "WARNING!!! Geometry.h -- element" 
	      << " type id unknown, cannot calculate centroid" << std::endl;
  }
} 

//Computes the area of a triangle
//pt1 - xyz values of 1st point
//pt2 - xyz values of 2nd point
//pt3 - xyz values of 3rd point
//rpt - return value, direction area (x,y,z) of the resultant triangle
template <class theType>
void CalcTriArea(const theType pt1[3], const theType pt2[3], const theType pt3[3], theType rpt[3])
{
  theType v1[3];
  theType v2[3];
  
  Subtract (pt2, pt1, v1);
  Subtract (pt3, pt1, v2);
  
  CrossProduct(v1, v2, rpt);

  Scale((theType)0.5, rpt, 3);
}

//computes the volume of a tetrahedron
//pt1 - xyz values of 1st point
//pt2 - xyz values of 2nd point
//pt3 - xyz values of 3rd point
//pt4 - xyz values of 3rd point
//returns volume of the tet
template <class theType>
theType CalcTetVol(const theType pt1[3], const theType pt2[3], const theType pt3[3], 
			const theType pt4[3])
{
  theType v1[3];
  theType v2[3];
  theType v3[3];
  theType tmp[3];
  theType vol;
 
  Subtract (pt2, pt1, v1);
  Subtract (pt3, pt1, v2);
  Subtract (pt4, pt1, v3);

  CrossProduct(v1, v2, tmp);
  vol = DotProduct(v3, tmp)/6.0;

  return vol;
}

//takes three points, returns one point and a normal vector
//return value is the magnitude of the normal before normalization
template <class theType>
theType GetPlaneDefinition(const theType pt1[3], const theType pt2[3], const theType pt3[3], theType rpt[3], theType rn[3])
{
  theType diff1[3];
  theType diff2[3];
  theType mag;

  // use centroid
  rpt[0] = (pt1[0] + pt2[0] + pt3[0])/3.0;
  rpt[1] = (pt1[1] + pt2[1] + pt3[1])/3.0;
  rpt[2] = (pt1[2] + pt2[2] + pt3[2])/3.0;
  
  // pt2 - pt1
  diff1[0] = pt2[0] - pt1[0];
  diff1[1] = pt2[1] - pt1[1];
  diff1[2] = pt2[2] - pt1[2];

  // pt3 - pt1
  diff2[0] = pt3[0] - pt1[0];
  diff2[1] = pt3[1] - pt1[1];
  diff2[2] = pt3[2] - pt1[2];

  CrossProduct(diff1, diff2, rn);
  mag = Normalize(rn, rn);

  return mag;
}

//computes the intersection of two planes, expects normalized vectors for normals
//pt1 - (xyz) of a point on the first plane
//n1 - normal (xyz) of first plane
//pt2 - (xyz) of a point on the second plane
//n2 - normal (xyz) of second plane
//rpt - (xyz) of a point on the intersection line
//rn - normal vector which follows the direction of the intersection
//return - angle between planes in radians
template <class theType>
theType FindPlanesIntersection(const theType pt1[3], const theType n1[3], const theType pt2[3], const theType n2[3],
			       theType rpt[3], theType rn[3])
{
  theType alpha;
  theType abs, max;
  Int indx;
  theType normMin;
  theType h1,h2;

  if(sizeof(theType) == sizeof(double)){
    normMin = 1.0e-15;
  }
  else if(sizeof(theType) == sizeof(float)){
    normMin = 1.0e-7;
  }

  //get direction of intersection
  CrossProduct(n1, n2, rn);

  if(Normalize(rn,rn) < normMin){
    //return zero radians between planes
    std::cerr << "WARNING: Planes are parallel!!" << std::endl;
    return (0);
  }

  //find parametric value of h which is valid for each plane
  //use the point we are given to define each plane
  h1 = DotProduct(n1, pt1);
  h2 = DotProduct(n2, pt2);

  //Find the intersection with the most stable coordinate planes
  max = CAbs(rn[0]);
  indx = 0;
  if(CAbs(rn[1]) > max){
    indx = 1;
    max = Abs(rn[1]);
  }
  if(CAbs(rn[2]) > max){
    indx = 2;
    max = CAbs(rn[2]);
  }
  //pick correct analytic solution
  switch(indx){
  case 0:
    //x-plane
    rpt[0] = 0.0;
    rpt[2] = (h2*n1[1] - n2[1]*h1) / (n2[2]*n1[1] - n2[1]*n1[2]);
    rpt[1] = (h1 - n1[2]*rpt[2]) / n1[1];
    break;
  case 1:
    //y-plane
    rpt[1] = 0.0;
    rpt[2] = (h2*n1[0] - n2[0]*h1) / (n2[2]*n1[0] - n2[0]*n1[2]);
    rpt[0] = (h1 - n1[2]*rpt[2]) / n1[0];
    break;
  case 2:
    //z-plane
    rpt[2] = 0.0;
    rpt[1] = (h2*n1[0] - n2[0]*h1) / (n2[1]*n1[0] - n2[0]*n1[1]);;
    rpt[0] = (h1 - n1[1]*rpt[1]) / n1[0];
    break;
  default:
    break;
  }
  
  //find cos (alpha) of the intersection
  alpha = DotProduct(n1, n2);
  //return alpha in radians
  return (acos(alpha));
}

//translate list of points by given distance
template <class theType>
void TranslatePoints(theType* pts, theType dx[3], Int n)
{
  Int i;
  theType* pt;
  for(i = 0; i < n; i++){
    pt = pts + i*3;

    pt[0] += dx[0];
    pt[1] += dx[1];
    pt[2] += dx[2];
  }
}

//rotate list of points by theta given
//pts - list of points (x1 y1 z1, x2 y2 z2, ...)
//apt - axis point. center of rotation
//axis - (xyz) normal of axis to rotate around
//theta (rads) - angle to rotate by
//n - number of points in list pts
template <class theType>
void RotatePoints(theType* pts, const theType apt[3], const theType axis[3], const theType theta, const Int n)
{
  Int i, j;
  theType* pt;
  theType r[3];
  theType v[3];
  theType mag;
  theType dv, dr;
  theType sn = sin(theta);
  theType cs1 = 1.0 - cos(theta);

  for(i = 0; i < n; i++){
    pt = pts + i*3;
    
    //get the radius vector
    Subtract(apt, pt, r);
    //remove components of r in direction of the rotation axis
    RemoveNormalComponent(r, axis, r);

    mag = Normalize(r,r);
    
    //no need to rotate if point is already on the rotation axis
    if(real(mag) > 1.0e-14){
      //get the other basis
      CrossProduct(axis,r,v);

      //get components along both basis vectors
      dv = mag*sn;
      dr = -mag*cs1;

      //move the point by dx
      for (j = 0; j < 3; j++){
	pt[j] += dv*v[j] + dr*r[j];
      }
  
    }
  }
} 

//Mirror a vector about a plane defined by normalized vector avec
template <class theType>
inline void MirrorVector(const theType* vin, const theType* avec, theType* vout)
{
  theType dot = 2.0*DotProduct(vin,avec);

  vout[0] = vin[0] - dot*avec[0];
  vout[1] = vin[1] - dot*avec[1];
  vout[2] = vin[2] - dot*avec[2];
}

//reverse a vector
template <class theType>
inline void ReverseVector(const theType* vin, theType* vout)
{
  vout[0] = -vin[0];
  vout[1] = -vin[1];
  vout[2] = -vin[2];
}

//distance from a line to a point
//returns linear distance and vector from the nearest point to the line
template <class theType>
theType DistanceLinePoint(const theType* linePt1, const theType* linePt2, const theType* oPt, 
			  theType* posVec)
{
  theType dist;
  theType r1[3];
  theType r2[3];
  theType rCross[3];
  theType rCMag;
  theType s[3];
  theType sMag;

  //find the minimum distance
  Subtract(oPt, linePt1, r1);
  Subtract(oPt, linePt2, r2);
  CrossProduct(r1, r2, rCross);
  rCMag = Magnitude(rCross);

  Subtract(linePt2, linePt1, s);
  sMag = Magnitude(s);

  dist = rCMag / sMag;

  //now find the parametric value of t along the line to get back the point
  Subtract(linePt1, oPt, r1);
  Subtract(linePt2, linePt1, r2);
  theType dot = DotProduct(r1, r2);
  theType s2 = sMag*sMag;
  //this parameter minimized the distance
  theType t = -dot/s2;

  //compute the actual point on the line
  theType newPt[3];
  newPt[0] = linePt1[0] + (linePt2[0]-linePt1[0])*t;
  newPt[1] = linePt1[1] + (linePt2[1]-linePt1[1])*t;
  newPt[2] = linePt1[2] + (linePt2[2]-linePt1[2])*t;

  //compute vector from the minimizing point to the query point
  Subtract(oPt, newPt, posVec);

  return dist;
}

template <class theType>
// computes the velocity of a point given a rotating field
// rotationPoint - center of rotation
// movingPoint - point to compute velocity of
// omegaRate - rotation rate of the field (i.e. rad/s)
// v - return value 
void VelocityFromRotation(const theType rotationPoint[3], const theType movingPoint[3], const theType omegaRate, theType v[3])
{
  // compute the position vector, r
  theType r[3];
  Subtract(rotationPoint, movingPoint, r);
  // velocity is w x r = v
  CrossProduct(omegaRate, r, v);
}

#endif
