#ifndef PARAM_PARSER_H__
#define PARAM_PARSER_H__

#include "general.h"
#include <string>
#include <vector>


//this class defines a keyword based parameter to search
//for in a parameter file
template <class Type>
class Parameter
{
public:
  Parameter(std::string key, Type* valRef, Type defaultValue, Type lowRange, Type highRange);
  ~Parameter(){};
  template <class Type2>
  void Copy(const Parameter<Type2>& toCopy);
  void Print() const;
  Bool ParseLine(std::string line);
  Type GetValue() const;
  std::string GetKey() const;
private:
  std::string key;
  Type* valRef;
  Type lowRange;
  Type highRange;
};

//parses a list of numerical values like key = [ 1, 2, 3, 4, ...]
template <class Type>
class ParameterList
{
public:
  ParameterList(std::string key, std::vector<Type>* valRef, std::vector<Type> defaultValue, 
		Type lowRange, Type highRange, Int maxValues = 999, Int minValues = 1);
  ~ParameterList(){};
  template <class Type2>
  void Copy(const ParameterList<Type2>& toCopy);
  void Print() const;
  Bool ParseLine(std::string line);
  std::vector<Type> GetValue() const;
  std::string GetKey() const;
private:
  std::string key;
  std::vector<Type>* valRef;
  Type lowRange;
  Type highRange;
  Int maxValues;
  Int minValues;
};

class ParameterEnum
{
public:
  ParameterEnum(std::string key, std::vector<std::string> enumKeys, Int* valRef, Int defaultValue);
  ~ParameterEnum(){};
  void Copy(const ParameterEnum& toCopy);
  void Print() const;
  Bool ParseLine(std::string line);
  Int GetValue() const;
  std::string GetKey() const;
private:
  std::string key;
  std::vector<std::string> enumKeys;
  Int* valRef;
};

class ParameterBool
{
public: 
  ParameterBool(std::string key, Bool* valRef, Bool defaultValue);
  ~ParameterBool(){};
  void Copy(const ParameterBool& toCopy);
  void Print() const;
  Bool ParseLine(std::string line);
  Bool GetValue() const;
  std::string GetKey() const;
private:
  std::string key;
  Bool* valRef;
};

class ParameterString
{
 public:
  ParameterString(std::string key, std::string* valRef, std::string defaultValue);
  ~ParameterString(){};
  void Copy(const ParameterString& toCopy);
  void Print() const;
  Bool ParseLine(std::string line);
  std::string GetValue() const;
  std::string GetKey() const;
 private:
  std::string key;
  std::string* valRef;
};

//parses a list of comma separated strings in a list like key = [ string1, string2, string3, ...]
class ParameterStringList
{
public:
  ParameterStringList(std::string key, std::vector<std::string>* valRef, std::string defaultValue);
  ~ParameterStringList(){};
  void Copy(const ParameterStringList& toCopy);
  void Print() const;
  Bool ParseLine(std::string line);
  std::vector<std::string> GetValue() const;
  std::string GetKey() const;
private:
  std::string key;
  std::vector<std::string>* valRef;
};

//include implementation
#include "parameterParser.tcc"

#endif
