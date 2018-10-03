#include "dataInfo.h"
#include "h5layer.h"
#include "exceptions.h"
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

DataInfo::DataInfo():
  name("NULL"), ndof(0), nvector(0), nscalar(0)
{
}

DataInfo::DataInfo(Int ndof, std::string bulkName):
  ndof(ndof), name(bulkName), nvector(0), nscalar(0)
{
  descriptorS.resize(ndof);
  descriptorV.resize(ndof);
  names.resize(ndof);
  for(Int i = 0; i < ndof; i++){
    descriptorS[i] = -1;
    descriptorV[i] = -1;
  }
}

void DataInfo::SetFromHDF(hid_t fileId, std::string directory, std::string dataName)
{
  HDF_ReadArrayAttribute(fileId, directory, dataName, "scalars", descriptorS);
  HDF_ReadArrayAttribute(fileId, directory, dataName, "vectors", descriptorV);
  ndof = descriptorS.size();
  names.resize(ndof);
  HDF_ReadStringAttribute(fileId, directory, dataName, "variable_names", names);
}

DataInfo::~DataInfo()
{
}

void DataInfo::AddVector(Int dof, std::string name)
{
  Int i;
  if(dof > ndof){
    std::stringstream ss;
    ss << "WARNING: degree of freedom " << dof << " greater than allocated space of "
       << ndof << " -- FAILING! "
       << " in data descriptor of name " << name << " !" << std::endl;
    Abort << ss.str();
    return;
  }
  for(i = dof; i < dof+3; i++){
    if(descriptorS[i] != -1){
      std::stringstream ss;
      ss << "WARNING: degree of freedom " << dof
	 << " has already been declared a scalar -- FAILING! " 
	 << " in data descriptor of name " << name << " !" << std::endl;
      Abort << ss.str();
      return;
    }
    descriptorV[i] = nvector;
    names[i] = name;
  }
  nvector++;
}

void DataInfo::AddScalar(Int dof, std::string name)
{
  if(dof > ndof){
    std::stringstream ss;
    ss << "WARNING: degree of freedom " << dof << " greater than allocated space of "
       << ndof << " -- FAILING! " 
       << " in data descriptor of name " << name << " !" << std::endl;
    Abort << ss.str();
    return;
  }
  if(descriptorV[dof] != -1){
    std::stringstream ss;
    ss << "WARNING: degree of freedom " << dof
       << " has already been declared a vector -- FAILING! " 
       << " in data descriptor of name " << name << " !" << std::endl;
    Abort << ss.str();
    return;
  }
  descriptorS[dof] = nscalar;
  names[dof] = name;
  nscalar++;
}

void DataInfo::Verify()
{
  Int i;
  for(i = 0; i < ndof; i++){
    if((descriptorS[i] == -1) && (descriptorV[i] == -1)){
      std::stringstream ss;
      ss << "WARNING: degree of freedom " <<  i  << " in data descriptor of name " 
	 << name << " is not set!" << std::endl;
      Abort << ss.str();
    }
  }
}

void DataInfo::Print() const
{
  Int i;
  std::cout << "BULK DATA DESCRIPTOR " << name << std::endl;
  for(i = 0; i < ndof; i++){
    std::cout << "Dof[" << i << "] Type: ";
    if(descriptorS[i] != -1){
      std::cout << "Scalar ";
    }
    else{
      std::cout << "Vector ";
    }
    std::cout << names[i] << std::endl;
  }
}

void DataInfo::WriteBinary(std::ofstream & fout, Int mode)
{
  //okay doing this sucks but we are having some huge problems
  //with the binary data which follows this... oi!

  Int i;
  std::string::size_type sz;
  if(mode == 1){
    //TODO: roll this write function out soon
    fout.write((char*)&ndof, sizeof(Int));
    sz = names.size();
    fout.write(reinterpret_cast<char*>(&sz), sizeof(std::string::size_type));
    fout.write(name.data(), sz);
  }
  for(i = 0; i < ndof; i++){
    sz = names[i].size();
    fout.write(reinterpret_cast<char*>(&sz), sizeof(std::string::size_type));
    fout.write(names[i].data(), sz);
  }
  for(i = 0; i < ndof; i++){
    fout.write((char*)&descriptorS[i], sizeof(Int));
  }
  for(i = 0; i < ndof; i++){
    fout.write((char*)&descriptorV[i], sizeof(Int));
  }
}

void DataInfo::ReadBinary(std::ifstream & fin)
{
  //okay doing this sucks but we are having some huge problems
  //with the binary data which follows this... oi!

  Int i;
  std::string::size_type sz;
  char temp[128];

  fin.read((char*)&ndof, sizeof(Int));
  names.resize(ndof);
  fin.read(reinterpret_cast<char*>(&sz), sizeof(std::string::size_type));
  fin.read(temp, sz);
  name = std::string(temp, temp + sz);

  for(i = 0; i < ndof; i++){
    fin.read(reinterpret_cast<char*>(&sz), sizeof(std::string::size_type));
    fin.read(temp, sz);
    names[i] = std::string(temp, temp + sz);
  }
  for(i = 0; i < ndof; i++){
    fin.read((char*)&descriptorS[i], sizeof(Int));
  }
  for(i = 0; i < ndof; i++){
    fin.read((char*)&descriptorV[i], sizeof(Int));
  }
  Verify();
}

Int DataInfo::GetNdof() const
{
  return ndof;
}

std::string DataInfo::GetName() const
{
  return name;
}

std::string DataInfo::GetDofName(Int dof) const
{
  return names[dof];
}

Bool DataInfo::DofIsVector(Int dof) const
{
  if(descriptorV[dof] != -1){
    return true;
  }
  return false;
}

Bool DataInfo::DofIsScalar(Int dof) const
{
  if(descriptorS[dof] != -1){
    return true;
  }
  return false;
}

const std::vector<std::string>& DataInfo::GetNames() const
{
  return names;
}

void DataInfo::WriteHDFAttribute(hid_t file_id, std::string directory)
{
  HDF_WriteStringAttribute(file_id, directory, name, "variable_names", names);
  HDF_WriteArrayAttribute(file_id, directory, name, "scalars", descriptorS);
  HDF_WriteArrayAttribute(file_id, directory, name, "vectors", descriptorV);
}

std::ostream& operator<<(std::ostream& os, const DataInfo& obj){
  os << "\t" << obj.name << "\n";
  os << "\t-----------------------------------------------\n";
  for(int i = 0; i < obj.names.size(); ++i){
    os << "\t" << obj.names[i] << "\t"
       << "\t" << obj.descriptorS[i] << "\t"
       << "\t" << obj.descriptorV[i] << "\n";
  }
  return os;
}

