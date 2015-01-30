#include "strings_util.h"
#include <sstream>
#include <iostream>
#include <vector>

template <class Type>
Parameter<Type>::Parameter(std::string key, Type* valRef, Type defaultValue, Type lowRange, Type highRange) :
  key(key), valRef(valRef), lowRange(lowRange), highRange(highRange)
{
  *valRef = defaultValue;
}

template <class Type> template <class Type2>
void Parameter<Type>::Copy(const Parameter<Type2>& toCopy)
{
  *valRef = toCopy.GetValue();
}

template <class Type>
Bool Parameter<Type>::ParseLine(std::string line)
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

    if(*valRef > highRange){
      std::cerr << "WARNING: Parameter " << key << " above high range -- " << highRange << std::endl;
    }
    if(*valRef < lowRange){
      std::cerr << "WARNING: Parameter " << key << " below low range -- " << lowRange << std::endl;
    }
    ss.clear();
    return true;
  }
  return false;
}

template <class Type>
Type Parameter<Type>::GetValue() const
{
  return *valRef;
}

template <class Type>
std::string Parameter<Type>::GetKey() const
{
  return key;
}

template <class Type>
void Parameter<Type>::Print() const
{
  std::cout << key << " <int, double>" << std::endl;
}

template <class Type>
ParameterList<Type>::ParameterList(std::string key, std::vector<Type>* valRef, std::vector<Type> defaultValue, 
				   Type lowRange, Type highRange, Int maxValues, Int minvalues) :
  key(key), valRef(valRef), lowRange(lowRange), highRange(highRange), maxValues(maxValues), minValues(minValues)
{
  (*valRef) = defaultValue;
}

template <class Type> template <class Type2>
void ParameterList<Type>::Copy(const ParameterList<Type2>& toCopy)
{
  const std::vector<Type2> srcVect = toCopy.GetValue();
  (*valRef).clear();
  for(Int i = 0; i < srcVect.size(); i++){
    (*valRef).push_back(srcVect[i]);
  }
}

template <class Type>
void ParameterList<Type>::Print() const
{
  std::cout << key << " <int, double - list>" << std::endl;
}

template <class Type>
Bool ParameterList<Type>::ParseLine(std::string line)
{
  size_t loc;
  std::string subline;
  std::stringstream ss;
  std::vector<std::string> temp;
  Type value;
  loc = line.find(key);
  if(loc != std::string::npos && loc == 0){
    (*valRef).clear();
    loc = line.find('=');
    loc += 1;
    subline = line.substr(loc);
    GetStringBetween("[", "]", subline, subline);
    temp = Tokenize(subline, ',');
    for(std::vector<std::string>::iterator it = temp.begin(); it != temp.end(); ++it){
      ss.clear();
      std::string& str = *it;
      ss << str;
      ss >> value;
      (*valRef).push_back(value);
    }
    for(typename std::vector<Type>::iterator it = (*valRef).begin(); it != (*valRef).end(); ++it){
      Type value = *it;
      if(value > highRange){
	std::cerr << "WARNING: Parameter " << key << " above high range -- " << highRange << std::endl;
      }
      if(value < lowRange){
	std::cerr << "WARNING: Parameter " << key << " below low range -- " << lowRange << std::endl;
      }
    }
    if((*valRef).size() > maxValues){
      std::cerr << "WARNING: Parameter " << key << " more than " << maxValues << " found" << std::endl;
    }
    if((*valRef).size() > maxValues){
      std::cerr << "WARNING: Parameter " << key << " less than " << minValues << " found" << std::endl;
    }
    ss.clear();
    return true;
  }
  return false;
}

template <class Type>
std::vector<Type> ParameterList<Type>::GetValue() const
{
  std::vector<Type> values;
  for(UInt i = 0; i < (*valRef).size(); i++){
    values.push_back( (*valRef)[i] );
  }
  return values;
}

template <class Type>
std::string ParameterList<Type>::GetKey() const
{
  return key;
}
