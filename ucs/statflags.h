#include "bitflags.h"

// This class is designed to hold status information in boolean form
// for a particular control volume in the mesh

class CVStat
{
public:
  CVStat() :
    numflags(4), bitfield(numflags)
  { };

  void Clear()
  {
    bitfield.clear();
  }

  void SetViscous()
  {
    bitfield.set(VISCOUS_LOC, true);
  };

  bool IsViscous()
  {
    return bitfield.get(VISCOUS_LOC);
  };

  void SetDirichlet()
  {
    bitfield.set(DIRICHLET_LOC, true);
  };

  bool IsDirichlet()
  {
    return bitfield.get(DIRICHLET_LOC);
  };

  ~CVStat()
  {  };


private:
  int numflags;
  BitField bitfield;
  
  enum{
    VISCOUS_LOC,
    DIRICHLET_LOC
  };

};
