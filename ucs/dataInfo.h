#ifndef DATA_INFO_HH__
#define DATA_INFO_HH__

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "general.h"
#include "h5layer.h"

class DataInfo
{
public:
  DataInfo();
  //dumb constructor w/o vector, scalar, or naming information
  DataInfo(Int ndof, std::string bulkName);
  ~DataInfo();
  void AddVector(Int dof, std::string name);
  void AddScalar(Int dof, std::string name);
  void Verify();
  void Print() const;
  void WriteBinary(std::ofstream & fout, Int mode = 0);
  void ReadBinary(std::ifstream & fin);
  Int GetNdof() const;
  std::string GetName() const;
  std::string GetDofName(Int dof) const;
  const std::vector<std::string>& GetNames() const;
  void WriteHDFAttribute(hid_t fileId, std::string directory);
  Bool DofIsScalar(Int dof) const;
  Bool DofIsVector(Int dof) const;
  //use HDF file and a path to a dataset to set interior
  void SetFromHDF(hid_t fileId, std::string directory, std::string dataName);
  friend std::ostream& operator<<(std::ostream& os, const DataInfo& obj);
  
private:
  Int ndof;
  std::string name;
  std::vector<std::string> names;
  //store scalar locations
  std::vector<Int> descriptorS;
  //store vector locations
  std::vector<Int> descriptorV;

  Int nvector;
  Int nscalar;
};



#endif
