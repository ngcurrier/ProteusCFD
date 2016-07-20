#ifndef POST_PLUGIN_H__
#define POST_PLUGIN_H__

#include <iostream>
#include "kernel.h"
#include "driver.h"
#include "solutionSpace.h"
#include <string>

//post plugin defines a virtual interface class for
//post-processing simulation results
template <class Type>
class PostPlugin
{
 public:
  PostPlugin(SolutionSpace<Type>& space, std::string name);
  virtual ~PostPlugin(){};

  //these are pure virtuals, must be defined as inherited
  virtual void Compute() = 0;
  virtual void WriteTabularHeader() = 0;
  virtual void Report() const = 0;

  //these are default implementations (utilities)
  std::string getName() const {return name;};
  Bool isName(std::string name__) const {return (name__ == name);};
  
 protected:
  SolutionSpace<Type>& space;
  std::string name;
};

// post surface integral is a derived class of post plugin
// which provides a framework for computing surface integrals
// on a particular solution space
template <class Type>
class PostSurfaceIntegral : public PostPlugin<Type>
{
public:
  PostSurfaceIntegral(SolutionSpace<Type>& space, std::string name, void(*bkernel)(B_KERNEL_ARGS));
  virtual ~PostSurfaceIntegral(){};
  virtual void Compute();
  virtual void WriteTabularHeader() = 0;
  virtual void Report() const = 0; 
protected:
  //this is the kernel where the core of the work happens for the plugin
  Kernel<Type> myKernel;
};

// post volume integral is a derived class of post plugin
// which provides a framework for computing volume integrals
// on a particular solution space
template <class Type>
class PostVolumeIntegral : public PostPlugin<Type>
{
public:
  PostVolumeIntegral(SolutionSpace<Type>& space, std::string name, void(*kernel)(KERNEL_ARGS));
  virtual ~PostVolumeIntegral(){};
  virtual void Compute();
  virtual void WriteTabularHeader() = 0;
  virtual void Report() const = 0; 
protected:
  //this is the kernel where the core of the work happens for the plugin
  Kernel<Type> myKernel;
};

//include implementations
#include "postPlugin.tcc"

#endif
