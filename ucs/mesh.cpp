#include "mesh.h"
#include "etypes.h"

void TranslateWinding(Int* nodes, Int translation[6][8], Int num_nodes, Int etype, Int to_other_format){
  Int tempnode[num_nodes];
  Int k;

  memcpy(tempnode,nodes,num_nodes*sizeof(Int));
  
  if(to_other_format){
    for (k = 0; k < num_nodes; k++){
      nodes[k] = tempnode[translation[etype][k]];
    }
  }
  else{
    for (k = 0; k < num_nodes; k++){
      // translation[k] gives our local node number
      // k = other format local node number
      nodes[translation[etype][k]] = tempnode[k];
    }
  }
}

