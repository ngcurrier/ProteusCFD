#ifndef POST_PLUGIN_H__
#define POST_PLUGIN_H__

//post plugin defines a virtual interface class for
//post-processing simulation results
template <class Type>
class PostPlugin
{
 public:
  PostPlugin();
  virtual ~PostPlugin();
  virtual void Compute() = 0;
  virtual void WriteTabularHeader() = 0;
  virtual void Report() const = 0;
  
 private:
  
};



#endif
