
template <class Type>
void InitIC(Int icid, SolutionSpace<Type>* space){
  Int i;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Param<Type>* p = space->param;
  Type* q = space->q;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int nnode = m->GetNumNodes();
  Type grad;

  std::cout << "CUSTOM IC: using custom ic: " << icid << std::endl;
  std::cout << "CUSTOM IC: please see implementation in source code" << std::endl;

  switch(icid)
    {
    case 0:
      //shouldn't ever be called b/c it should default to farfield in EqnSet
      break;
    case 1:
      //shock tube for compressible regime (1.0 x ?? x ??) grid
      for(i = 0; i < nnode; i++){
	if(real(m->xyz[i*3 + 0]) <= 0.5){
	  q[i*nvars + 0] = 1.0;
	  q[i*nvars + 1] = 0.0;
	  q[i*nvars + 2] = 0.0;
	  q[i*nvars + 3] = 0.0;	   
	  q[i*nvars + 4] = 1.0/p->gamma;
	}
	else{
	  q[i*nvars + 0] = 1.0/8.0;
	  q[i*nvars + 1] = 0.0;
	  q[i*nvars + 2] = 0.0;
	  q[i*nvars + 3] = 0.0;	   
	  q[i*nvars + 4] = (1.0/p->gamma)/10.0;
	}
      }
      break;
    case 2:
      //shock box -- 2D version (1.0 x 1.0 x ??) grid
      for(i = 0; i < nnode; i++){
	if((real(m->xyz[i*3 + 0]) <= 0.5) && (real(m->xyz[i*3 + 1]) <= 0.5)){
	  q[i*nvars + 0] = 1.0;
	  q[i*nvars + 1] = 0.0;
	  q[i*nvars + 2] = 0.0;
	  q[i*nvars + 3] = 0.0;	   
	  q[i*nvars + 4] = 1.0/p->gamma;
	}
	else{
	  q[i*nvars + 0] = 1.0/8.0;
	  q[i*nvars + 1] = 0.0;
	  q[i*nvars + 2] = 0.0;
	  q[i*nvars + 3] = 0.0;	   
	  q[i*nvars + 4] = (1.0/p->gamma)/10.0;
	}
	eqnset->ComputeAuxiliaryVariables(&q[i*nvars+0]);
      }	
      break;
    case 3:
      //shock box -- 3D version (1.0 x 1.0 x 1.0) grid
      for(i = 0; i < nnode; i++){
	if((real(m->xyz[i*3 + 0]) <= 0.5) && (real(m->xyz[i*3 + 1]) <= 0.5) && (real(m->xyz[i*3 + 2]) <= 0.5)){
	  q[i*nvars + 0] = 1.0;
	  q[i*nvars + 1] = 0.0;
	  q[i*nvars + 2] = 0.0;
	  q[i*nvars + 3] = 0.0;	   
	  q[i*nvars + 4] = 1.0/p->gamma;
	}
	else{
	  q[i*nvars + 0] = 1.0/8.0;
	  q[i*nvars + 1] = 0.0;
	  q[i*nvars + 2] = 0.0;
	  q[i*nvars + 3] = 0.0;	   
	  q[i*nvars + 4] = (1.0/p->gamma)/10.0;
	}
	eqnset->ComputeAuxiliaryVariables(&q[i*nvars+0]);
      }
      break;
    case 4:
      //linear progression in the x direction - gradient testing case
      grad = 0.5;
      for(i = 0; i < nnode; i++){
	q[i*nvars + 0] = real(m->xyz[i*3 + 0])*grad;
	q[i*nvars + 1] = real(m->xyz[i*3 + 0])*grad;
	q[i*nvars + 2] = real(m->xyz[i*3 + 0])*grad;
	q[i*nvars + 3] = real(m->xyz[i*3 + 0])*grad;
	q[i*nvars + 4] = real(m->xyz[i*3 + 0])*grad;
      }
      break;
    case 5:
      //2D shock with box in the middle of a 1x1 domain
      for(i = 0; i < nnode; i++){
	if((real(m->xyz[i*3 + 0]) <= 0.75) && (real(m->xyz[i*3 + 1]) <= 0.75)
	   && (real(m->xyz[i*3 + 0]) >= 0.25) && (real(m->xyz[i*3 + 1]) >= 0.25) ){
	  q[i*nvars + 0] = 1.0;
	  q[i*nvars + 1] = 0.0;
	  q[i*nvars + 2] = 0.0;
	  q[i*nvars + 3] = 0.0;	   
	  q[i*nvars + 4] = 1.0/p->gamma;
	}
	else{
	  q[i*nvars + 0] = 1.0/8.0;
	  q[i*nvars + 1] = 0.0;
	  q[i*nvars + 2] = 0.0;
	  q[i*nvars + 3] = 0.0;	   
	  q[i*nvars + 4] = (1.0/p->gamma)/10.0;
	}
	eqnset->ComputeAuxiliaryVariables(&q[i*nvars+0]);
      }	
      break;
    case 6:
      //set for quiescent flow (still)
      Int vloc = eqnset->GetMomentumLocation();
      for(i = 0; i < nnode; i++){
	q[i*nvars + vloc+0] = 0.0;
	q[i*nvars + vloc+1] = 0.0;
	q[i*nvars + vloc+2] = 0.0;
      }
      eqnset->ComputeAuxiliaryVariables(&q[i*nvars+0]);
      break;
    }
  
  return;
}
