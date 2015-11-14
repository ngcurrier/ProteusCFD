#ifndef ELEMENT_CLASS_H__
#define ELEMENT_CLASS_H__

#include "general.h"
#include "etypes.h"
#include "mem_util.h"
#include <iostream>
#include <vector>

template <class Type> class SpatialFunctor;

template <class Type>
class Element
{
 public:

  Element<Type>();
  template <class Type2>
  Element<Type>(const Element<Type2>& elemToCopy);
  virtual ~Element<Type>();
  
  void Init(Int* nodes);
  Int GetNodes(Int** nodes);
  Int GetNodes(Int* nodes) const;
  Int GetNnodes() const;
  virtual Int GetType() const;
  void SetFactag(Int factag);
  Int GetFactag() const;
  bool IsActive() const;
  bool operator==(const Element<Type>& e2) const;
  bool operator!=(const Element<Type>& e2) const;
  bool operator<(const Element<Type>& e2) const;
  template <class Type2>
  friend std::ostream& operator<<(std::ostream& os, const Element<Type2>& e);
  //simple indicator of likeness without directly comparing nodes
  Int GetNodeSum() const;
  Int GetMinNode() const;
  Int GetMaxNode() const;
  virtual void GetGaussWeights(std::vector<Type>& gpoints, std::vector<Type>& gweights, Int degree) const
  { std::cerr << "WARNING: weights not defined for element" << std::endl;};
  virtual void EvaluateShapeFunctions(Type* elemXYZ, Type xi, Type eta, Type& det, Type* sfunc) const
  { std::cerr << "WARNING: shape functions not defined for element" << std::endl;};
  //pass a function pointer which takes x,y,z and returns a value
  virtual Type IntegrateFunction(SpatialFunctor<Type>* f, Type* elemXYZ, Int degree) const
  { std::cerr << "WARNING: integrate function not defined for element" << std::endl; return 0.0;};

protected:
  Int* nodes;
  Int nnodes;
  Bool active;
  Bool typeset;
  Int factag;
private:
};

template <class Type>
class Triangle : public Element<Type>
{
public:
  Triangle();
  template <class Type2>
  Triangle(const Triangle<Type2>& elemToCopy);
  Int GetType() const;
  
  void GetGaussWeights(std::vector<Type>& gpoints, std::vector<Type>& gweights, Int degree) const;
  void EvaluateShapeFunctions(Type* elemXYZ, Type xi, Type eta, Type& det, Type* sfunc) const;
  //elem xyz is a list of points for the element
  Type IntegrateFunction(SpatialFunctor<Type>* f, Type* elemXYZ, Int degree) const;
protected:
private:
};

template <class Type>
class Quadrilateral : public Element<Type>
{
public:
  Quadrilateral();
  template <class Type2>
  Quadrilateral(const Quadrilateral<Type2>& elemToCopy);
  Int GetType() const;
  void GetGaussWeights(std::vector<Type>& gpoints, std::vector<Type>& gweights, Int degree) const;
  void EvaluateShapeFunctions(Type* elemXYZ, Type xi, Type eta, Type& det, Type* sfunc) const;
  //elem xyz is a list of points for the element
  Type IntegrateFunction(SpatialFunctor<Type>* f, Type* elemXYZ, Int degree) const;
 protected:
private:
};

template <class Type>
class Tetrahedron : public Element<Type>
{
public:
  Tetrahedron();
  template <class Type2>
  Tetrahedron(const Tetrahedron<Type2>& elemToCopy);
  Int GetType() const;
protected:
private:
};

template <class Type>
class Pyramid : public Element<Type>
{
public:
  Pyramid();
  template <class Type2>
  Pyramid(const Pyramid<Type2>& elemToCopy);
  Int GetType() const;
protected:
private:
};

template <class Type>
class Prism : public Element<Type>
{
public:
  Prism();
  template <class Type2>
  Prism(const Prism<Type2>& elemToCopy);
  Int GetType() const;
protected:
private:
};

template <class Type>
class Hexahedron : public Element<Type>
{
public:
  Hexahedron();
  template <class Type2>
  Hexahedron(const Hexahedron<Type2>& elemToCopy);
  Int GetType() const;
protected:
private:
};

//used for set based insertion so compares are deep not at ptr level
template <class Type>
struct ElementPtrComp
{
  bool operator()(const Element<Type>* e1, const Element<Type>* e2) const
  {
    return ((*e1) < (*e2));
  }
};

//include implementations
#include "elementClass.tcc"

#endif
