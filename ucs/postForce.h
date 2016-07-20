#ifndef POST_FORCE_H__
#define POST_FORCE_H__

#include "postPlugin.h"
#include "solutionSpace.h"

//forward declarations
template <class Type> class SolutionSpace;

//post plugin which computes forces applied to a composite body
template <class Type>
class PostForcePlugin : public PostSurfaceIntegral<Type>
{
 public:
  PostForcePlugin(SolutionSpace<Type>& space);
  virtual ~PostForcePlugin();
  void WriteTabularHeader();
  void Report() const;

  //TODO: make these things private once this gets cleaned up (bodies,cp)
  CompositeBody<Type>* bodies; //list of the composite bodies
  Type* cp;    //coefficient of pressure
  
 protected:
  
 private:

  //things we store only on the surface
  Type* cf;    //skin friction coefficient
};

template <class Type>
void Post_Force_Kernel(B_KERNEL_ARGS);

//include implementations
#include "postForce.tcc"

#endif
