template <class Type>
PostPlugin<Type>::PostPlugin(SolutionSpace<Type>& space):
  space(space)
{
}

template <class Type>
PostSurfaceIntegral<Type>::PostSurfaceIntegral(SolutionSpace<Type>& space):
  PostPlugin<Type>(space)
{
}

template <class Type>
void PostSurfaceIntegral<Type>::Compute()
{
  Kernel<Type> surfaceKernel(myBoundaryKernel);
  BdriverNoScatter(&(this->space), surfaceKernel, 0, NULL);
}

template <class Type>
PostVolumeIntegral<Type>::PostVolumeIntegral(SolutionSpace<Type>& space):
  PostPlugin<Type>(space)
{
  Kernel<Type> volumeKernel(myVolumeKernel);
  DriverNoScatter(&(this->space), volumeKernel, 0, NULL);
}

template <class Type>
void PostVolumeIntegral<Type>::Compute()
{

}
