#ifndef POST_PLUGIN_H__
#define POST_PLUGIN_H__

#include <iostream>
#include "kernel.h"
#include "driver.h"
#include "solutionSpace.h"

//post plugin defines a virtual interface class for
//post-processing simulation results
template <class Type>
class PostPlugin
{
 public:
  PostPlugin(SolutionSpace<Type>& space);
  virtual ~PostPlugin();
  virtual void Compute() = 0;
  virtual void WriteTabularHeader() = 0;
  virtual void Report() const = 0;
  
 protected:
  SolutionSpace<Type>& space;
};

// post surface integral is a derived class of post plugin
// which provides a framework for computing surface integrals
// on a particular solution space
template <class Type>
class PostSurfaceIntegral : public PostPlugin<Type>
{
public:
  PostSurfaceIntegral(SolutionSpace<Type>& space);
  virtual ~PostSurfaceIntegral();
  virtual void Compute();
  virtual void WriteTabularHeader() = 0;
  virtual void Report() const = 0; 
private:
  //this is the kernel where the core of the work happens for the plugin
  void myBoundaryKernel(B_KERNEL_ARGS) = 0;
};

// post volume integral is a derived class of post plugin
// which provides a framework for computing volume integrals
// on a particular solution space
template <class Type>
class PostVolumeIntegral : public PostPlugin<Type>
{
public:
  PostVolumeIntegral(SolutionSpace<Type>& space);
  virtual void Compute();
  virtual void WriteTabularHeader() = 0;
  virtual void Report() const = 0; 
private:
  //this is the kernel where the core of the work happens for the plugin
  void myVolumeKernel(B_KERNEL_ARGS) = 0;
};

//include implementations
#include "postPlugin.tcc"

#endif
