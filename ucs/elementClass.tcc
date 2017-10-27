#include "etypes.h"
#include "mem_util.h"
#include "macros.h"
#include "strings_util.h"
#include "functor.h"
#include <string>
#include <sstream>
#include <iomanip>

template <class Type>
Element<Type>::Element() :
  nodes(NULL), nnodes(0), active(false), typeset(false), factag(-1)
{}

template <class Type>
Element<Type>::~Element()
{
  delete [] nodes;
}

template <class Type>
Int Element<Type>::GetNodes(Int** nodes)
{
  *nodes = this->nodes;
  return nnodes;
}

template <class Type>
Int Element<Type>::GetNodes(Int* nodes) const
{
  for(Int i = 0; i < nnodes; ++i){
    nodes[i] = this->nodes[i];
  }
  return nnodes;
}

template <class Type>
Int Element<Type>::GetNnodes() const
{
  return nnodes;
}

template <class Type>
void Element<Type>::Init(Int* nodes)
{
  if(!typeset){
    std::cerr << "Element type not yet declared, cannot Init" << std::endl;
    return;
  }  
  memcpy(this->nodes, nodes, sizeof(Int)*nnodes);
}

template <class Type>
Int Element<Type>::GetType() const
{
  return (-1);
}

template <class Type>
void Element<Type>::SetFactag(Int factag)
{
  this->factag = factag;
}

template <class Type>
Int Element<Type>::GetFactag() const
{
  return factag;
}

template <class Type>
bool Element<Type>::operator==(const Element<Type>& e2) const
{
  if(GetType() != e2.GetType()){
    return false;
  }
  if(nnodes != e2.GetNnodes()){
    return false;
  }
  if(factag != e2.GetFactag()){
    return false;
  }
  Int* nodes2 = new Int[e2.GetNnodes()];
  e2.GetNodes(nodes2);
  Bool check;
  //check to see if each element contains the same nodes
  for(Int i = 0; i < nnodes; i++){
    check = false;
    Int n1 = nodes[i];
    for(Int j = 0; j < nnodes; j++){
      Int n2 = nodes2[j];
      if(n1 == n2){
	check = true;
	break;
      }
    }
    if(check == false){
      delete [] nodes2;
      return false;
    }
  }
  delete [] nodes2;
  return true;
}

template <class Type>
bool Element<Type>::operator!=(const Element& e2) const
{
  return !(*this == e2);
}

template <class Type>
bool Element<Type>::operator<(const Element& e2) const
{
  std::stringstream e1ss;
  std::stringstream e2ss;

  Int c = 99999999;
  Int w = 16;

  e1ss << std::setfill('0');
  e2ss << std::setfill('0');
  e1ss << std::setw(w) << GetType();
  e2ss << std::setw(w) << e2.GetType();
  e1ss << std::setw(w) << GetNnodes();
  e2ss << std::setw(w) << e2.GetNnodes();
  e1ss << std::setw(w) << GetFactag();
  e2ss << std::setw(w) << e2.GetFactag();
  e1ss << std::setw(w) << GetNodeSum();
  e2ss << std::setw(w) << e2.GetNodeSum();

  Int maxNodes = 9;
  for(Int i = 0; i < maxNodes; i++){
    if(i < GetNnodes()){
      e1ss << std::setw(w) << nodes[i];
    }
    else{
      e1ss << std::setw(w) << 0;
    }
  }
  Int* nodes2 = new Int[e2.GetNnodes()];
  e2.GetNodes(nodes2);
  for(Int i = 0; i < maxNodes; i++){
    if(i < e2.GetNnodes()){
      e2ss << std::setw(w) << nodes2[i];
    }
    else{
      e2ss << std::setw(w) << 0;
    }
  }
  delete [] nodes2;

  return (e1ss.str() < e2ss.str());
}

template <class Type>
bool Element<Type>::IsActive() const
{
  return active;
}
 
template <class Type>
Int Element<Type>::GetNodeSum() const
{
  Int sum = 0;
  for(Int i = 0; i < nnodes; i++){
    sum += nodes[i];
  }
  return sum;
}

template <class Type>
Int Element<Type>::GetMaxNode() const
{
  Int max = nodes[0];
  for(Int i = 1; i < nnodes; i++){
    if(nodes[i] > max) max = nodes[i];
  }
  return max;
}

template <class Type>
Int Element<Type>::GetMinNode() const
{
  Int min = nodes[0];
  for(Int i = 1; i < nnodes; i++){
    if(nodes[i] < min) min = nodes[i];
  }
  return min;
}

template <class Type> template <class Type2>
Element<Type>::Element(const Element<Type2>& elemToCopy) :
  nodes(NULL), nnodes(elemToCopy.GetNnodes()), active(elemToCopy.IsActive()), 
  typeset(false), factag(elemToCopy.GetFactag())
{ 
  nodes = new Int[nnodes];
}

template <class Type>
std::ostream& operator<<(std::ostream& os, const Element<Type>& e)
{
  os << "Type: " << e.GetType() << " Nnodes: " << e.GetNnodes() << 
    " Factag: " << e.GetFactag() << std::endl;
  Int* nodes = new Int[e.GetNnodes()];
  e.GetNodes(nodes);
  for(Int i = 0; i < e.GetNnodes(); i++){
    os << "Node[" << i << "]: " << nodes[i] << std::endl;
  }
  os << "NodeSum: " << e.GetNodeSum() << std::endl;
  delete [] nodes;
  return os;
}

template <class Type>
Triangle<Type>::Triangle()
{
  this->nnodes = 3;
  this->active = true;
  this->typeset = true;
  this->nodes = new Int[3];
}

template <class Type>
Int Triangle<Type>::GetType() const
{
  return TRI;
}

template <class Type> template <class Type2>
Triangle<Type>::Triangle(const Triangle<Type2>& elemToCopy) :
  Element<Type>(elemToCopy)
{
  elemToCopy.GetNodes(this->nodes);
  this->typeset = true;
}

template <class Type>
Type Triangle<Type>::IntegrateFunction(SpatialFunctor<Type>* f, Type* elemXYZ, Int degree) const
{ 
  std::vector<Type> gpoints;
  std::vector<Type> gweights;
  Type sfunc[this->nnodes];

  //get standard tabulated gauss weights and locations
  GetGaussWeights(gpoints, gweights, degree);

  Type I = 0.0;

  //loop over gauss points
  for(Int i = 0; i < gweights.size(); i++){
    Type gpt_xi = gpoints[i*2 + 0];
    Type gpt_eta = gpoints[i*2 + 1];
    Type det;
    EvaluateShapeFunctions(elemXYZ, gpt_xi, gpt_eta, det, sfunc);
    Type K = det*gweights[i];
    
    //get contributions from each element nodal point
    Type x = 0.0;
    Type y = 0.0;
    Type z = 0.0;
    for(Int nn = 0; nn < this->nnodes; nn++){
      x += elemXYZ[nn*3 + 0]*sfunc[nn];
      y += elemXYZ[nn*3 + 1]*sfunc[nn];
    }

    //evaluate function at the gauss point
    Type contrib = f->Evaluate(x, y, z);
    for(Int nn = 0; nn < this->nnodes; nn++){
      I += 0.5*K*contrib*sfunc[nn];
    }
  }
  return I;
}  

template <class Type>
void Triangle<Type>::GetGaussWeights(std::vector<Type>& gpoints, std::vector<Type>& gweights, Int degree) const
{
  std::vector<Type>& x = gpoints;
  std::vector<Type>& w = gweights;

  if(degree == 1){
    Int n = 1;
    x.resize(n*2);
    w.resize(n);

    x[0*2 + 0] = 0.3333333333333333;
    x[0*2 + 1] = 0.3333333333333333;
    
    w[0] = 1.0000000000000000000;
  }
  else if(degree == 2){
    Int n = 3;
    x.resize(n*2);
    w.resize(n);

    x[0*2 + 0] = 0.1666666666666667;
    x[0*2 + 1] = 0.1666666666666667;
    
    x[1*2 + 0] = 0.6666666666666667;
    x[1*2 + 1] = 0.1666666666666667;

    x[2*2 + 0] = 0.1666666666666667;
    x[2*2 + 1] = 0.6666666666666667;

    w[0] = 0.3333333333333333;
    w[1] = 0.3333333333333333;
    w[2] = 0.3333333333333333;
  }
  else if(degree == 3){
    Int n = 4;
    x.resize(n*2);
    w.resize(n);

    x[0*2 + 0] = 0.3333333333333333;
    x[0*2 + 1] = 0.3333333333333333;

    x[1*2 + 0] = 0.20;
    x[1*2 + 1] = 0.20;

    x[2*2 + 0] = 0.60;
    x[2*2 + 1] = 0.20;

    x[3*2 + 0] = 0.20;
    x[3*2 + 1] = 0.60;
    
    w[0] = -0.5625;
    w[1] =  0.5208333333333333;
    w[2] =  0.5208333333333333;
    w[3] =  0.5208333333333333;
  }
  else if(degree == 4){
    Int n = 6;
    x.resize(n*2);
    w.resize(n);

    x[0*2 + 0] = 0.091576213509771;
    x[0*2 + 1] = 0.091576213509771;

    x[1*2 + 0] = 0.445948490915965;
    x[1*2 + 1] = 0.108103018168070;

    x[2*2 + 0] = 0.816847572980459;
    x[2*2 + 1] = 0.091576213509771;

    x[3*2 + 0] = 0.108103018168070;
    x[3*2 + 1] = 0.445948490915965;

    x[4*2 + 0] = 0.445948490915965;
    x[4*2 + 1] = 0.445948490915965;

    x[5*2 + 0] = 0.091576213509771;
    x[5*2 + 1] = 0.816847572980459;

    w[0] = 0.109951743655322;
    w[1] = 0.223381589678011;
    w[2] = 0.109951743655322;
    w[3] = 0.223381589678011;
    w[4] = 0.223381589648011;
    w[5] = 0.109951743655322;
  }
  else{
    std::cerr << "STOP.... greater than 4th degree (6 point) triangle rule not allowed..." << std::endl;
  }
}

template <class Type>
void Triangle<Type>::EvaluateShapeFunctions(Type* elemXYZ, Type xi, Type eta, Type& det, Type* sfunc) const
{
  Type J[2][2];
  Type dsf[this->nnodes][2];

#if 0
  //node locations in "master" quad
  Type xnode[this->nnodes][2];

  xnode[0][0] = 0.0;
  xnode[0][1] = 0.0;

  xnode[1][0] =  1.0;
  xnode[1][1] =  0.0;

  xnode[2][0] =  0.0;
  xnode[2][1] =  1.0;
#endif
  
  //using linear lagrange interpolation on three nodes
  sfunc[0] = 1.0 - xi - eta;
  sfunc[1] = xi;
  sfunc[2] = eta;

  //derivative of the shape function inode, for jacobian
  dsf[0][0] = -1.0;
  dsf[0][1] = -1.0;

  dsf[1][0] = 1.0;
  dsf[1][1] = 0.0;

  dsf[2][0] = 0.0;
  dsf[2][1] = 1.0;

  //compute the jacobian matrix (J) and its inverse (Jinv) and determinant (det)
  for(Int i = 0; i < 2; i++){
    for(Int j = 0; j < 2; j++){
      J[i][j] = 0.0;
      for(Int inode = 0; inode < this->nnodes; inode++){
	J[i][j] = J[i][j] + dsf[inode][i]*elemXYZ[inode*3 + j];
      }
    }
  }
  
  det = J[0][0]*J[1][1] - J[0][1]*J[1][0];

#if 0
  //turned off to reduce computation
  Type Jinv[2][2];
  Jinv[0][0] =  J[1][1]/det;
  Jinv[1][1] =  J[0][0]/det;
  Jinv[0][1] = -J[0][1]/det;
  Jinv[1][0] = -J[1][0]/det;
#endif
}

template <class Type>
Quadrilateral<Type>::Quadrilateral()
{
  this->nnodes = 4;
  this->active = true;
  this->typeset = true;
  this->nodes = new Int[4];
}
 
template <class Type>
Int Quadrilateral<Type>::GetType() const
{
  return QUAD;
}

template <class Type> template <class Type2>
Quadrilateral<Type>::Quadrilateral(const Quadrilateral<Type2>& elemToCopy) :
  Element<Type>(elemToCopy)
{
  elemToCopy.GetNodes(this->nodes);
  this->typeset = true;
}

template <class Type>
Type Quadrilateral<Type>::IntegrateFunction(SpatialFunctor<Type>* f, Type* elemXYZ, Int degree) const
{
  std::vector<Type> gpoints;
  std::vector<Type> gweights;
  Type sfunc[this->nnodes];

  //get standard tabulated gauss weights and locations
  GetGaussWeights(gpoints, gweights, degree);

  Type I = 0.0;

  //loop over grid of gauss points
  for(Int ni = 0; ni < degree; ni++){
    Type gpt_xi = gpoints[ni];
    for(Int nj = 0; nj < degree; nj++){
      Type gpt_eta = gpoints[nj];
      Type det;
      EvaluateShapeFunctions(elemXYZ, gpt_xi, gpt_eta, det, sfunc);
      Type K = det*gweights[ni]*gweights[nj];

      //get contributions from each element nodal point
      Type x = 0.0;
      Type y = 0.0;
      Type z = 0.0;
      for(Int nn = 0; nn < this->nnodes; nn++){
	x += elemXYZ[nn*3 + 0]*sfunc[nn];
	y += elemXYZ[nn*3 + 1]*sfunc[nn];
      }

      //evaluate function at the gauss point
      Type contrib = f->Evaluate(x, y, z);
      for(Int nn = 0; nn < this->nnodes; nn++){
	I += K*contrib*sfunc[nn];
      }
    }
  }
  return I;
}

template <class Type>
void Quadrilateral<Type>::EvaluateShapeFunctions(Type* elemXYZ, Type xi, Type eta, Type& det, Type* sfunc) const
{
  Type J[2][2];
  Type Jinv[2][2];
  Type dsf[this->nnodes][2];

  //node locations in "master" quad
  Type xnode[this->nnodes][2];

  xnode[0][0] = -1.0;
  xnode[0][1] = -1.0;

  xnode[1][0] =  1.0;
  xnode[1][1] = -1.0;

  xnode[2][0] =  1.0;
  xnode[2][1] =  1.0;

  xnode[3][0] = -1.0;
  xnode[3][1] =  1.0;
  
  //using linear lagrange interpolation on four nodes
  for(Int inode = 0; inode < this->nnodes; inode++){
    //these are just used to get the sign right for the basis functions
    Type xnatural = xnode[inode][0];
    Type ynatural = xnode[inode][1];

    //really, 1 - xi, 1 + xi, 1 - eta, 1 + eta as appropriate
    Type xi_basis = (1.0 + xi*xnatural);
    Type eta_basis = (1.0 + eta*ynatural);

    sfunc[inode] = 0.25 * xi_basis * eta_basis;
    //derivative of the shape function inode, for jacobian
    dsf[inode][0] = 0.25 * xnatural * eta_basis;
    dsf[inode][1] = 0.25 * ynatural * xi_basis;
  }

  //compute the jacobian matrix (J) and its inverse (Jinv) and determinant (det)
  for(Int i = 0; i < 2; i++){
    for(Int j = 0; j < 2; j++){
      J[i][j] = 0.0;
      for(Int inode = 0; inode < this->nnodes; inode++){
	J[i][j] = J[i][j] + dsf[inode][i]*elemXYZ[inode*3 + j];
      }
    }
  }
  
  det = J[0][0]*J[1][1] - J[0][1]*J[1][0];

  Jinv[0][0] =  J[1][1]/det;
  Jinv[1][1] =  J[0][0]/det;
  Jinv[0][1] = -J[0][1]/det;
  Jinv[1][0] = -J[1][0]/det;
}

template <class Type>
void Quadrilateral<Type>::GetGaussWeights(std::vector<Type>& gpoints, std::vector<Type>& gweights, Int degree) const
{
  std::vector<Type>& x = gpoints;
  std::vector<Type>& w = gweights;
  Int n = degree;
  x.resize(n);
  w.resize(n);

  //...quadrature rules.... from Jim Newman's 2D integration code
  if (n  ==  1 ){
    w[0] = 2.0; 
    x[0] = 0.0; 
  }
  else if (n  ==  2 ){
    w[0] =  1.0; 
    w[1] =  w[0]; 
    
    x[0] = -0.577350269189625764509148780502;
    x[1] =  0.577350269189625764509148780502;
  }
  else if (n  ==  3 ){
    w[0] =  5.0 / 9.0;  
    w[1] =  8.0 / 9.0; 
    w[2] =  5.0 / 9.0; 
     
    x[0] =  -0.774596669241483377035853079956; 
    x[1] =   0.000000000000000; 
    x[2] =   0.774596669241483377035853079956; 
  }
  else if (n  ==  4 ){
    w[0] = 0.347854845137453857373063949222; 
    w[1] = 0.652145154862546142626936050778; 
    w[2] = 0.652145154862546142626936050778; 
    w[3] = 0.347854845137453857373063949222; 
     
    x[0] = - 0.861136311594052575223946488893; 
    x[1] = - 0.339981043584856264802665759103; 
    x[2] =   0.339981043584856264802665759103; 
    x[3] =   0.861136311594052575223946488893; 
  }
  else if (n  ==  5 ){
    w[0] = 0.236926885056189087514264040720; 
    w[1] = 0.478628670499366468041291514836; 
    w[2] = 0.568888888888888888888888888889; 
    w[3] = 0.478628670499366468041291514836; 
    w[4] = 0.236926885056189087514264040720; 
     
    x[0] =  -0.906179845938663992797626878299; 
    x[1] =  -0.538469310105683091036314420700; 
    x[2] =   0.0;
    x[3] =   0.538469310105683091036314420700; 
    x[4] =   0.906179845938663992797626878299; 
  }
  else if (n  ==  6 ) {
    w[0] = 0.171324492379170345040296142173; 
    w[1] = 0.360761573048138607569833513838; 
    w[2] = 0.467913934572691047389870343990; 
    w[3] = 0.467913934572691047389870343990; 
    w[4] = 0.360761573048138607569833513838; 
    w[5] = 0.171324492379170345040296142173; 
     
    x[0] =  -0.932469514203152027812301554494; 
    x[1] =  -0.661209386466264513661399595020; 
    x[2] =  -0.238619186083196908630501721681; 
    x[3] =   0.238619186083196908630501721681; 
    x[4] =   0.661209386466264513661399595020; 
    x[5] =   0.932469514203152027812301554494; 
  }
  else if (n  ==  7 ) {
    w[0] = 0.129484966168869693270611432679; 
    w[1] = 0.279705391489276667901467771424;
    w[2] = 0.381830050505118944950369775489; 
    w[3] = 0.417959183673469387755102040816; 
    w[4] = 0.381830050505118944950369775489; 
    w[5] = 0.279705391489276667901467771424; 
    w[6] = 0.129484966168869693270611432679; 
    
    x[0] = - 0.949107912342758524526189684048; 
    x[1] = - 0.741531185599394439863864773281; 
    x[2] = - 0.405845151377397166906606412077; 
    x[3] =   0.0; 
    x[4] =   0.405845151377397166906606412077; 
    x[5] =   0.741531185599394439863864773281; 
    x[6] =   0.949107912342758524526189684048; 
  }
  else if (n  ==  8 ){
    w[0] = 0.101228536290376259152531354310; 
    w[1] = 0.222381034453374470544355994426; 
    w[2] = 0.313706645877887287337962201987; 
    w[3] = 0.362683783378361982965150449277; 
    w[4] = 0.362683783378361982965150449277; 
    w[5] = 0.313706645877887287337962201987; 
    w[6] = 0.222381034453374470544355994426; 
    w[7] = 0.101228536290376259152531354310; 
     
    x[0] = - 0.960289856497536231683560868569; 
    x[1] = - 0.796666477413626739591553936476; 
    x[2] = - 0.525532409916328985817739049189; 
    x[3] = - 0.183434642495649804939476142360; 
    x[4] =   0.183434642495649804939476142360; 
    x[5] =   0.525532409916328985817739049189; 
    x[6] =   0.796666477413626739591553936476; 
    x[7] =   0.960289856497536231683560868569; 
  }
  else if (n  ==  9 ){
    w[0] = 0.0812743883615744119718921581105; 
    w[1] = 0.180648160694857404058472031243; 
    w[2] = 0.260610696402935462318742869419; 
    w[3] = 0.312347077040002840068630406584; 
    w[4] = 0.330239355001259763164525069287; 
    w[5] = 0.312347077040002840068630406584; 
    w[6] = 0.260610696402935462318742869419; 
    w[7] = 0.180648160694857404058472031243; 
    w[8] = 0.0812743883615744119718921581105; 
     
    x[0] = - 0.968160239507626089835576202904; 
    x[1] = - 0.836031107326635794299429788070; 
    x[2] = - 0.613371432700590397308702039341; 
    x[3] = - 0.324253423403808929038538014643; 
    x[4] =   0.0;
    x[5] =   0.324253423403808929038538014643; 
    x[6] =   0.613371432700590397308702039341; 
    x[7] =   0.836031107326635794299429788070; 
    x[8] =   0.968160239507626089835576202904; 
  }
  else if (n  ==  10 ){
    w[0] =  0.0666713443086881375935688098933; 
    w[1] =  0.149451349150580593145776339658; 
    w[2] =  0.219086362515982043995534934228; 
    w[3] =  0.269266719309996355091226921569; 
    w[4] =  0.295524224714752870173892994651; 
    w[5] =  0.295524224714752870173892994651; 
    w[6] =  0.269266719309996355091226921569; 
    w[7] =  0.219086362515982043995534934228; 
    w[8] =  0.149451349150580593145776339658; 
    w[9] = 0.0666713443086881375935688098933; 
     
    x[0] =  - 0.973906528517171720077964012084; 
    x[1] =  - 0.865063366688984510732096688423; 
    x[2] =  - 0.679409568299024406234327365115; 
    x[3] =  - 0.433395394129247290799265943166; 
    x[4] =  - 0.148874338981631210884826001130; 
    x[5] =    0.148874338981631210884826001130; 
    x[6] =    0.433395394129247290799265943166; 
    x[7] =    0.679409568299024406234327365115; 
    x[8] =    0.865063366688984510732096688423; 
    x[9] =   0.973906528517171720077964012084; 
  }
  else if (n == 11){
    w[0] = 0.5566856711617350E-01; 
    w[1] = 0.1255803694649042E+00; 
    w[2] = 0.1862902109277336E+00; 
    w[3] = 0.2331937645919900E+00; 
    w[4] = 0.2628045445102460E+00; 
    w[5] = 0.2729250867779003E+00; 
    w[6] = 0.2628045445102467E+00; 
    w[7] = 0.2331937645919913E+00; 
    w[8] = 0.1862902109277353E+00; 
    w[9] = 0.1255803694649053E+00; 
    w[10] = 0.5566856711617379E-01; 

    x[0] = -.9782286581460569E+00; 
    x[1] = -.8870625997680954E+00; 
    x[2] = -.7301520055740491E+00; 
    x[3] = -.5190961292068117E+00; 
    x[4] = -.2695431559523447E+00; 
    x[5] = 0.0000000000000000E+00; 
    x[6] = 0.2695431559523451E+00; 
    x[7] = 0.5190961292068120E+00; 
    x[8] = 0.7301520055740494E+00; 
    x[9] = 0.8870625997680953E+00; 
    x[10] = 0.9782286581460570E+00; 
  }
  else if (n == 12){
    w[0] = 0.4717533638651189E-01; 
    w[1] = 0.1069393259953186E+00; 
    w[2] = 0.1600783285433465E+00; 
    w[3] = 0.2031674267230662E+00; 
    w[4] = 0.2334925365383550E+00; 
    w[5] = 0.2491470458134032E+00; 
    w[6] = 0.2491470458134031E+00; 
    w[7] = 0.2334925365383548E+00; 
    w[8] = 0.2031674267230660E+00; 
    w[9] = 0.1600783285433459E+00; 
    w[10] = 0.1069393259953177E+00; 
    w[11] = 0.4717533638651136E-01; 
  
    x[0] = -.9815606342467189E+00; 
    x[1] = -.9041172563704747E+00; 
    x[2] = -.7699026741943045E+00; 
    x[3] = -.5873179542866176E+00; 
    x[4] = -.3678314989981801E+00; 
    x[5] = -.1252334085114690E+00; 
    x[6] = 0.1252334085114690E+00; 
    x[7] = 0.3678314989981801E+00; 
    x[8] = 0.5873179542866173E+00; 
    x[9] = 0.7699026741943045E+00; 
    x[10] = 0.9041172563704748E+00; 
    x[11] = 0.9815606342467192E+00; 
  }
  else if (n == 13){
    w[0] = 0.4048400476531592E-01; 
    w[1] = 0.9212149983772858E-01; 
    w[2] = 0.1388735102197876E+00; 
    w[3] = 0.1781459807619462E+00; 
    w[4] = 0.2078160475368887E+00; 
    w[5] = 0.2262831802628974E+00; 
    w[6] = 0.2325515532308739E+00; 
    w[7] = 0.2262831802628975E+00; 
    w[8] = 0.2078160475368894E+00; 
    w[9] = 0.1781459807619455E+00; 
    w[10] = 0.1388735102197863E+00; 
    w[11] = 0.9212149983772767E-01; 
    w[12] = 0.4048400476531557E-01; 
    
    x[0] = -.9841830547185881E+00; 
    x[1] = -.9175983992229779E+00; 
    x[2] = -.8015780907333099E+00; 
    x[3] = -.6423493394403401E+00; 
    x[4] = -.4484927510364467E+00; 
    x[5] = -.2304583159551348E+00; 
    x[6] = 0.0000000000000000E+00; 
    x[7] = 0.2304583159551348E+00; 
    x[8] = 0.4484927510364468E+00; 
    x[9] = 0.6423493394403401E+00; 
    x[10] = 0.8015780907333099E+00; 
    x[11] = 0.9175983992229779E+00; 
    x[12] = 0.9841830547185881E+00; 
  }
  else if (n == 14){
    w[0] = 0.3511946033175159E-01; 
    w[1] = 0.8015808715975993E-01; 
    w[2] = 0.1215185706879026E+00; 
    w[3] = 0.1572031671581927E+00; 
    w[4] = 0.1855383974779368E+00; 
    w[5] = 0.2051984637212945E+00; 
    w[6] = 0.2152638534631576E+00; 
    w[7] = 0.2152638534631579E+00; 
    w[8] = 0.2051984637212967E+00; 
    w[9] = 0.1855383974779394E+00; 
    w[10] = 0.1572031671581947E+00; 
    w[11] = 0.1215185706879038E+00; 
    w[12] = 0.8015808715976029E-01; 
    w[13] = 0.3511946033175129E-01; 
     
    x[0] = -.9862838086968122E+00; 
    x[1] = -.9284348836635739E+00; 
    x[2] = -.8272013150697646E+00; 
    x[3] = -.6872929048116854E+00; 
    x[4] = -.5152486363581539E+00; 
    x[5] = -.3191123689278896E+00; 
    x[6] = -.1080549487073434E+00; 
    x[7] = 0.1080549487073439E+00; 
    x[8] = 0.3191123689278900E+00; 
    x[9] = 0.5152486363581541E+00; 
    x[10] = 0.6872929048116853E+00; 
    x[11] = 0.8272013150697648E+00; 
    x[12] = 0.9284348836635734E+00; 
    x[13] = 0.9862838086968123E+00; 
  }
  else if (n == 15){
    w[0] = 0.3075324199611686E-01; 
    w[1] = 0.7036604748810761E-01; 
    w[2] = 0.1071592204671713E+00; 
    w[3] = 0.1395706779261536E+00; 
    w[4] = 0.1662692058169937E+00; 
    w[5] = 0.1861610000155617E+00; 
    w[6] = 0.1984314853271112E+00; 
    w[7] = 0.2025782419255614E+00; 
    w[8] = 0.1984314853271118E+00; 
    w[9] = 0.1861610000155627E+00; 
    w[10] = 0.1662692058169951E+00; 
    w[11] = 0.1395706779261555E+00; 
    w[12] = 0.1071592204671726E+00; 
    w[13] = 0.7036604748810819E-01; 
    w[14] = 0.3075324199611669E-01; 
     
    x[0] = -.9879925180204854E+00; 
    x[1] = -.9372733924007057E+00; 
    x[2] = -.8482065834104273E+00; 
    x[3] = -.7244177313601701E+00; 
    x[4] = -.5709721726085388E+00; 
    x[5] = -.3941513470775635E+00; 
    x[6] = -.2011940939974344E+00; 
    x[7] = -.2220446049250313E-15; 
    x[8] = 0.2011940939974346E+00; 
    x[9] = 0.3941513470775635E+00; 
    x[10] = 0.5709721726085387E+00; 
    x[11] = 0.7244177313601700E+00; 
    x[12] = 0.8482065834104271E+00; 
    x[13] = 0.9372733924007058E+00; 
    x[14] = 0.9879925180204854E+00; 
  }
  else if (n == 16){
    w[0] = 0.2715245941175374E-01; 
    w[1] = 0.6225352393864749E-01; 
    w[2] = 0.9515851168249226E-01; 
    w[3] = 0.1246289712555333E+00; 
    w[4] = 0.1495959888165762E+00; 
    w[5] = 0.1691565193950020E+00; 
    w[6] = 0.1826034150449232E+00; 
    w[7] = 0.1894506104550686E+00; 
    w[8] = 0.1894506104550693E+00; 
    w[9] = 0.1826034150449242E+00; 
    w[10] = 0.1691565193950039E+00; 
    w[11] = 0.1495959888165779E+00; 
    w[12] = 0.1246289712555344E+00; 
    w[13] = 0.9515851168249294E-01; 
    w[14] = 0.6225352393864737E-01; 
    w[15] = 0.2715245941175337E-01; 
     
    x[0] = -.9894009349916499E+00; 
    x[1] = -.9445750230732328E+00; 
    x[2] = -.8656312023878314E+00; 
    x[3] = -.7554044083550031E+00; 
    x[4] = -.6178762444026438E+00; 
    x[5] = -.4580167776572275E+00; 
    x[6] = -.2816035507792587E+00; 
    x[7] = -.9501250983763732E-01; 
    x[8] = 0.9501250983763765E-01; 
    x[9] = 0.2816035507792590E+00; 
    x[10] = 0.4580167776572276E+00; 
    x[11] = 0.6178762444026438E+00; 
    x[12] = 0.7554044083550030E+00; 
    x[13] = 0.8656312023878316E+00; 
    x[14] = 0.9445750230732325E+00; 
    x[15] = 0.9894009349916499E+00; 
  }
  else if (n == 17){
    w[0] = 0.2414830286854770E-01; 
    w[1] = 0.5545952937398704E-01; 
    w[2] = 0.8503614831717891E-01; 
    w[3] = 0.1118838471934039E+00; 
    w[4] = 0.1351363684685252E+00; 
    w[5] = 0.1540457610768101E+00; 
    w[6] = 0.1680041021564503E+00; 
    w[7] = 0.1765627053669930E+00; 
    w[8] = 0.1794464703562071E+00; 
    w[9] = 0.1765627053669935E+00; 
    w[10] = 0.1680041021564503E+00; 
    w[11] = 0.1540457610768105E+00; 
    w[12] = 0.1351363684685254E+00; 
    w[13] = 0.1118838471934039E+00; 
    w[14] = 0.8503614831717914E-01; 
    w[15] = 0.5545952937398678E-01; 
    w[16] = 0.2414830286854736E-01; 
     
    x[0] = -.9905754753144176E+00; 
    x[1] = -.9506755217687675E+00; 
    x[2] = -.8802391537269865E+00; 
    x[3] = -.7815140038968011E+00; 
    x[4] = -.6576711592166906E+00; 
    x[5] = -.5126905370864769E+00; 
    x[6] = -.3512317634538762E+00; 
    x[7] = -.1784841814958478E+00; 
    x[8] = 0.0000000000000000E+00; 
    x[9] = 0.1784841814958479E+00; 
    x[10] = 0.3512317634538763E+00; 
    x[11] = 0.5126905370864770E+00; 
    x[12] = 0.6576711592166908E+00; 
    x[13] = 0.7815140038968014E+00; 
    x[14] = 0.8802391537269858E+00; 
    x[15] = 0.9506755217687677E+00; 
    x[16] = 0.9905754753144174E+00;
  } 
  else if (n == 18){
    w[0] = 0.2161601352648302E-01; 
    w[1] = 0.4971454889496964E-01; 
    w[2] = 0.7642573025488883E-01; 
    w[3] = 0.1009420441062869E+00; 
    w[4] = 0.1225552067114782E+00; 
    w[5] = 0.1406429146706500E+00; 
    w[6] = 0.1546846751262646E+00; 
    w[7] = 0.1642764837458321E+00; 
    w[8] = 0.1691423829631434E+00; 
    w[9] = 0.1691423829631430E+00; 
    w[10] = 0.1642764837458328E+00; 
    w[11] = 0.1546846751262664E+00; 
    w[12] = 0.1406429146706515E+00; 
    w[13] = 0.1225552067114795E+00; 
    w[14] = 0.1009420441062878E+00; 
    w[15] = 0.7642573025488943E-01; 
    w[16] = 0.4971454889496965E-01; 
    w[17] = 0.2161601352648310E-01; 

    x[0] = -.9915651684209310E+00; 
    x[1] = -.9558239495713980E+00; 
    x[2] = -.8926024664975556E+00; 
    x[3] = -.8037049589725231E+00; 
    x[4] = -.6916870430603528E+00; 
    x[5] = -.5597708310739478E+00; 
    x[6] = -.4117511614628426E+00; 
    x[7] = -.2518862256915055E+00; 
    x[8] = -.8477501304173507E-01; 
    x[9] = 0.8477501304173540E-01; 
    x[10] = 0.2518862256915056E+00; 
    x[11] = 0.4117511614628427E+00; 
    x[12] = 0.5597708310739475E+00; 
    x[13] = 0.6916870430603532E+00; 
    x[14] = 0.8037049589725229E+00; 
    x[15] = 0.8926024664975557E+00; 
    x[16] = 0.9558239495713977E+00; 
    x[17] = 0.9915651684209309E+00; 
  }
  else if (n == 19){
    w[0] = 0.1946178822972596E-01; 
    w[1] = 0.4481422676569922E-01; 
    w[2] = 0.6904454273764057E-01;
    w[3] = 0.9149002162244917E-01; 
    w[4] = 0.1115666455473329E+00; 
    w[5] = 0.1287539625393350E+00; 
    w[6] = 0.1426067021736055E+00; 
    w[7] = 0.1527660420658583E+00; 
    w[8] = 0.1589688433939536E+00; 
    w[9] = 0.1610544498487835E+00; 
    w[10] = 0.1589688433939550E+00; 
    w[11] = 0.1527660420658612E+00; 
    w[12] = 0.1426067021736085E+00; 
    w[13] = 0.1287539625393376E+00; 
    w[14] = 0.1115666455473350E+00; 
    w[15] = 0.9149002162245097E-01; 
    w[16] = 0.6904454273764168E-01; 
    w[17] = 0.4481422676569959E-01; 
    w[18] = 0.1946178822972679E-01; 
   
    x[0] = -.9924068438435840E+00; 
    x[1] = -.9602081521348302E+00; 
    x[2] = -.9031559036148180E+00; 
    x[3] = -.8227146565371426E+00; 
    x[4] = -.7209661773352292E+00; 
    x[5] = -.6005453046616811E+00; 
    x[6] = -.4645707413759608E+00; 
    x[7] = -.3165640999636299E+00; 
    x[8] = -.1603586456402253E+00; 
    x[9] = 0.1110223024625157E-15; 
    x[10] = 0.1603586456402255E+00; 
    x[11] = 0.3165640999636298E+00; 
    x[12] = 0.4645707413759610E+00; 
    x[13] = 0.6005453046616810E+00; 
    x[14] = 0.7209661773352294E+00; 
    x[15] = 0.8227146565371428E+00; 
    x[16] = 0.9031559036148178E+00; 
    x[17] = 0.9602081521348300E+00; 
    x[18] = 0.9924068438435845E+00; 
  }
  else if (n == 20){
    w[0] = 0.1761400713915191E-01; 
    w[1] = 0.4060142980038711E-01; 
    w[2] = 0.6267204833410912E-01; 
    w[3] = 0.8327674157670478E-01; 
    w[4] = 0.1019301198172406E+00; 
    w[5] = 0.1181945319615187E+00; 
    w[6] = 0.1316886384491771E+00; 
    w[7] = 0.1420961093183825E+00; 
    w[8] = 0.1491729864726041E+00; 
    w[9] = 0.1527533871307253E+00; 
    w[10] = 0.1527533871307241E+00; 
    w[11] = 0.1491729864726025E+00; 
    w[12] = 0.1420961093183820E+00; 
    w[13] = 0.1316886384491767E+00; 
    w[14] = 0.1181945319615187E+00; 
    w[15] = 0.1019301198172408E+00; 
    w[16] = 0.8327674157670528E-01; 
    w[17] = 0.6267204833410926E-01; 
    w[18] = 0.4060142980038687E-01; 
    w[19] = 0.1761400713915236E-01; 
   
    x[0] = -.9931285991850947E+00; 
    x[1] = -.9639719272779137E+00; 
    x[2] = -.9122344282513259E+00; 
    x[3] = -.8391169718222187E+00; 
    x[4] = -.7463319064601510E+00; 
    x[5] = -.6360536807265151E+00; 
    x[6] = -.5108670019508268E+00; 
    x[7] = -.3737060887154198E+00; 
    x[8] = -.2277858511416453E+00; 
    x[9] = -.7652652113349734E-01; 
    x[10] = 0.7652652113349734E-01; 
    x[11] = 0.2277858511416453E+00; 
    x[12] = 0.3737060887154198E+00; 
    x[13] = 0.5108670019508271E+00; 
    x[14] = 0.6360536807265151E+00; 
    x[15] = 0.7463319064601508E+00; 
    x[16] = 0.8391169718222189E+00; 
    x[17] = 0.9122344282513259E+00; 
    x[18] = 0.9639719272779138E+00; 
    x[19] = 0.9931285991850949E+00;
  }
  else{
    std::cerr << "STOP.... greater than 20-point rule not allowed..." << std::endl;
  }
}

template <class Type>
Tetrahedron<Type>::Tetrahedron()
{
  this->nnodes = 4;
  this->active = true;
  this->typeset = true;
  this->nodes = new Int[4];
}

template <class Type>
Int Tetrahedron<Type>::GetType() const
{
  return TET;
}

template <class Type> template <class Type2>
Tetrahedron<Type>::Tetrahedron(const Tetrahedron<Type2>& elemToCopy) :
  Element<Type>(elemToCopy)
{
  elemToCopy.GetNodes(this->nodes);
  this->typeset = true;
}

template <class Type>
Pyramid<Type>::Pyramid()
{
  this->nnodes = 5;
  this->active = true;
  this->typeset = true;
  this->nodes = new Int[5];
}

template <class Type>
Int Pyramid<Type>::GetType() const
{
  return PYRAMID;
}

template <class Type> template <class Type2>
Pyramid<Type>::Pyramid(const Pyramid<Type2>& elemToCopy) :
  Element<Type>(elemToCopy)
{
  elemToCopy.GetNodes(this->nodes);
  this->typeset = true;
}

template <class Type>
Prism<Type>::Prism()
{
  this->nnodes = 6;
  this->active = true;
  this->typeset = true;
  this->nodes = new Int[6];
}

template <class Type>
Int Prism<Type>::GetType() const
{
  return PRISM;
}

template <class Type> template <class Type2>
Prism<Type>::Prism(const Prism<Type2>& elemToCopy) :
  Element<Type>(elemToCopy)
{
  elemToCopy.GetNodes(this->nodes);
  this->typeset = true;
}

template <class Type>
Hexahedron<Type>::Hexahedron() 
{
  this->nnodes = 8;
  this->active = true;
  this->typeset = true;
  this->nodes = new Int[8];
}

template <class Type>
Int Hexahedron<Type>::GetType() const 
{
  return HEX;
}

template <class Type> template <class Type2>
Hexahedron<Type>::Hexahedron(const Hexahedron<Type2>& elemToCopy) :
  Element<Type>(elemToCopy)
{
  elemToCopy.GetNodes(this->nodes);
  this->typeset = true;
}

