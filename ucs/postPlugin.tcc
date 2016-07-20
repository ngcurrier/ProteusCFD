#include "kernel.h"

template <class Type>
PostPlugin<Type>::PostPlugin(SolutionSpace<Type>& space, std::string name):
  space(space), name(name)
{
}

template <class Type>
PostSurfaceIntegral<Type>::PostSurfaceIntegral(SolutionSpace<Type>& space, std::string name, void(*bkernel)(B_KERNEL_ARGS)):
  PostPlugin<Type>(space, name), myKernel(bkernel)
{
}

template <class Type>
void PostSurfaceIntegral<Type>::Compute()
{
  BdriverNoScatter(&(this->space), myKernel, 0, NULL);
}

template <class Type>
PostVolumeIntegral<Type>::PostVolumeIntegral(SolutionSpace<Type>& space, std::string name, void(*kernel)(KERNEL_ARGS)):
  PostPlugin<Type>(space, name), myKernel(kernel)
{
}

template <class Type>
void PostVolumeIntegral<Type>::Compute()
{
  DriverNoScatter(&(this->space), myKernel, 0, NULL);
}
