#include "parameterParser.h"
#include "strings_util.h"
#include <sstream>

ParameterEnum::ParameterEnum(std::string key, std::vector<std::string> enumKeys, Int* valRef, Int defaultValue) :
  key(key), enumKeys(enumKeys), valRef(valRef)
{
  *valRef = defaultValue;
}

void ParameterEnum::Copy(const ParameterEnum& toCopy)
{
  *valRef = toCopy.GetValue();
}

void ParameterEnum::Print() const
{
  std::cout << key << " <Enum>" << std::endl;
}

Bool ParameterEnum::ParseLine(std::string line)
{
  size_t loc;
  std::string subline;
  loc = line.find(key);
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    Int i = 0;
    //eliminate any leading spaces
    size_t loc2 = 0;
    while(loc2 != std::string::npos){
      loc2 = subline.find(' ');
      if(loc2 != std::string::npos && loc2 == 0){
	subline = subline.substr(loc2+1);
	continue;
      }
      break;
    }
    //eliminate any trailing spaces
    while(loc2 != std::string::npos){
      loc2 = subline.rfind(' ');
      if(loc2 != std::string::npos && loc2 == subline.size()-1){
	subline = subline.substr(0, subline.size()-1);
	continue;
      }
      break;
    }
    for(std::vector<std::string>::iterator it = enumKeys.begin(); it != enumKeys.end(); ++it){
      std::string& key2 = *it;
      //look for an exact match or we might match similar keywords line compressible for in[compressible]
      if(subline == key2){
	*valRef = i;
    	return true;
      }
      i++;
    }
    std::cerr << "WARNING: key found " << key << " for enumerated values, parsed invalid option \'" 
	      << subline << "\'" << std::endl;
    std::cerr << "\tEnumerated match was not found, possible values are: " << std::endl;
    for(std::vector<std::string>::iterator it = enumKeys.begin(); it != enumKeys.end(); ++it){
      std::string& key2 = *it;
      std::cerr << "\t* " << key2 << "\n";
    }
    return true;
  }
  return false;
}

Int ParameterEnum::GetValue() const
{
  return *valRef;
}

std::string ParameterEnum::GetKey() const
{
  return key;
}

ParameterBool::ParameterBool(std::string key, Bool* valRef, Bool defaultValue) :
  key(key), valRef(valRef)
{
  *valRef = defaultValue;
}

void ParameterBool::Copy(const ParameterBool& toCopy)
{
  *valRef = toCopy.GetValue();
}

void ParameterBool::Print() const
{
  std::cout << key << " <boolean>" << std::endl;
}

Bool ParameterBool::ParseLine(std::string line)
{
  size_t loc;
  std::string subline;
  std::stringstream ss;
  loc = line.find(key);
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> *valRef;
    ss.clear();
    return true;
  }
  return false;
}

Bool ParameterBool::GetValue() const
{
  return *valRef;
}

std::string ParameterBool::GetKey() const
{
  return key;
}

ParameterString::ParameterString(std::string key, std::string* valRef, std::string defaultValue) :
  key(key), valRef(valRef)
{
  *valRef = defaultValue;
}

void ParameterString::Copy(const ParameterString& toCopy)
{
  *valRef = toCopy.GetValue();
}

void ParameterString::Print() const
{
  std::cout << key << " kind(string)" << std::endl;
}

Bool ParameterString::ParseLine(std::string line)
{
  size_t loc;
  std::string subline;
  std::stringstream ss;
  loc = line.find(key);
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    ss << subline;
    ss >> *valRef;
    ss.clear();
    return true;
  }
  return false;
}

std::string ParameterString::GetValue() const 
{
  return *valRef;
}

std::string ParameterString::GetKey() const
{
  return key;
}

ParameterStringList::ParameterStringList(std::string key, std::vector<std::string>* valRef, 
					 std::string defaultValue) :
  key(key), valRef(valRef)
{
  (*valRef).push_back(defaultValue);
}

void ParameterStringList::Copy(const ParameterStringList& toCopy)
{
  std::vector<std::string> vals = toCopy.GetValue();
  (*valRef).clear();
  for(UInt i = 0; i < vals.size(); i++){
    (*valRef).push_back(vals[i]);
  }
}

void ParameterStringList::Print() const
{
  std::cout << key << " <string - list>" << std::endl;
}

Bool ParameterStringList::ParseLine(std::string line)
{
  size_t loc;
  std::string subline;
  std::vector<std::string> temp;
  loc = line.find(key);
  if(loc != std::string::npos && loc == 0){
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    GetStringBetween("[", "]", subline, subline);
    temp = Tokenize(subline, ',');
    for(std::vector<std::string>::iterator it = temp.begin(); it != temp.end(); ++it){
      (*valRef).push_back(*it);
    }
    for(std::vector<std::string>::iterator it = (*valRef).begin();  it != (*valRef).end(); ++it){
      std::string& value = *it;
    }
    return true;
  }
  return false;
}

std::vector<std::string> ParameterStringList::GetValue() const
{
  return *valRef;
}

std::string ParameterStringList::GetKey() const
{
  return key;
}
