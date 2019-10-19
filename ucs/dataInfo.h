#ifndef DATA_INFO_HH__
#define DATA_INFO_HH__

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "general.h"
#include "h5layer.h"

template <class Type>
class DataInfo
{
public:
  DataInfo();
  //dumb constructor w/o vector, scalar, or naming information
  DataInfo(Int ndof, std::string bulkName);
  ~DataInfo();
  void AddVector(Int dof, std::string name, Type refValue);
  void AddScalar(Int dof, std::string name, Type refValue);
  void Verify();
  void Print() const;
  void WriteBinary(std::ofstream & fout, Int mode = 0);
  void ReadBinary(std::ifstream & fin);
  Int GetNdof() const;
  std::string GetName() const;
  std::string GetDofName(Int dof) const;
  Type GetDofReferenceValue(Int dof) const;
  const std::vector<std::string>& GetNames() const;
  void WriteHDFAttribute(hid_t fileId, std::string directory);
  Bool DofIsScalar(Int dof) const;
  Bool DofIsVector(Int dof) const;
  //use HDF file and a path to a dataset to set interior
  void SetFromHDF(hid_t fileId, std::string directory, std::string dataName);
  template <class Type2>
  friend std::ostream& operator<<(std::ostream& os, const DataInfo<Type2>& obj);
  
private:
  Int ndof;
  std::string name;
  std::vector<std::string> names;
  //store scalar locations
  std::vector<Int> descriptorS;
  //store vector locations
  std::vector<Int> descriptorV;
  // reference value to dimensionalize  the variable
  std::vector<Type> refValues;

  Int nvector;
  Int nscalar;
};

// include implementation
#include "dataInfo.tcc"


#endif
