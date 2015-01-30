#ifndef INTERPOLATION_H__
#define INTERPOLATION_H__

#include <iostream>
#include <cmath>

#include "macros.h"
#include "etypes.h"
#include "general.h"
#include "geometry.h"

//interpolate a value for a point inside a cell
//=============================================
//etypes defined in enumeration in uns_base.h
//x_new are coordinates for interpolation point
//q is vector of nodes solutions indexed by nodes[i]*nvars
//q_new is array to return value in 
//neqn is the number of variables we actually care about sampling ()
//nvars is number of variables present for each node in q vector, the stride
//nodes contains node ids to look up in lists xyz and q
//Methods:
//0 - inverse distance weighting
template <class theType, class theType2>
void InterpolateElement(const theType etype, const theType* nodes, const theType2* xyz,
			const theType2* q, const theType2* x_new, theType2* q_new, const theType neqn,
			const theType nvars, const Int method = 0)
{
  if(method == 0){
    switch(etype){
    case TRI:
      InverseWeighting(3, nodes, xyz, q, x_new, q_new, neqn, nvars);
      break;
    case QUAD:
      InverseWeighting(4, nodes, xyz, q, x_new, q_new, neqn, nvars);
      break;
    case TET:
      InverseWeighting(5, nodes, xyz, q, x_new, q_new, neqn, nvars);
      break;
    case PYRAMID:
      InverseWeighting(6, nodes, xyz, q, x_new, q_new, neqn, nvars);
      break;
    case PRISM:
      InverseWeighting(7, nodes, xyz, q, x_new, q_new, neqn, nvars);
      break;
    case HEX:
      InverseWeighting(8, nodes, xyz, q, x_new, q_new, neqn, nvars);
      break;
    default:
      std::cerr << "WARNING!!! interpolation.h -- element" 
		<< " type id unknown, cannot calculate interpolation" << std::endl;
    }
  }
  return;
}

//uses simple Shepard's inverse method
//nnodes - number of nodes we are weighting from
//nodes - linear array with node ids in it (i.e. 2, 5, 7, 8)
//xyz - linear array with coords in it for nodes in order (i.e. xyz-n1, xyz-n2, xyz-n3, ...)
//q - linear array with interpolated value in it (i.e. q-n1, q-n2, q-n3, ...)
//x_new - xyz coords of point to interpolate to
//q_new - value of the interpolation
//neqn - number of values in q to interpolate
//nvars - size of stride to next q sub vector ( q[3*(neqn + ...)] = q[3*nvars] )
//p - power to weight nearest points, for R->inf and p <= ndim weighting is 
//    dominated by farther points, this may not be desirable if nnodes is large
template <class theType, class theType2>
void InverseWeighting(const theType nnodes, const theType* nodes, const theType2* xyz,
		      const theType2* q, const theType2* x_new, theType2* q_new, 
		      const theType neqn, const theType nvars, const theType2 p = 2.0)
{
  Int i,j;
  theType2* w = (theType2*)alloca(sizeof(theType2)*nnodes);
  theType2 wsum = 0.0;
  //calculate weighting function for each node
  for(i = 0; i < nnodes; i++){
    theType2 d = Distance(&xyz[nodes[i]*3], x_new);
    //avoid divide by zero errors
    if(real(d) > 1.0e-18){
      if(p == 1.0){
	w[i] = 1.0 / d;
      }
      else{
	w[i] = 1.0 / pow(d,p);
      }
      wsum += w[i];
    }
    else{
      //we are on the point all weight goes here
      w[i] = 1.0;
      for(j = 0; j < nnodes; j++){
	if(i != j){
	  w[j] = 0.0;
	}
      }
      wsum = 1.0;
      break;
    }
  }
  //calculate the weighted variables
  for(j = 0; j < neqn; j++){
    q_new[j] = 0.0;
    for(i = 0; i < nnodes; i++){
      q_new[j] += (w[i] / wsum) * q[nodes[i]*nvars + j];
    }
  }
  return;
}
#endif
