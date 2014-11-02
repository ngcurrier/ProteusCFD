#ifndef STRINGS_UTIL_H__
#define STRINGS_UTIL_H__

#include "general.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>

std::string UpperToLower(std::string str);
std::string LowerToUpper(std::string str);

//returns string between first set of left, right strings... 
//if pos is given search starts from there
Int GetStringBetween(const std::string left, const std::string right, const std::string strIn, std::string& strOut, size_t pos = 0);

//returns string with section between first set of left, right strings removed... 
//if pos is given search starts from there
Int StripStringBetween(const std::string left, const std::string right, const std::string strIn, std::string& strOut, size_t pos = 0);

//strips the string 'stripped' from base and returns the result in strOut
//return value is the number of times the string 'stripped' was found
Int StripFromString(std::string base, std::string stripped, std::string& strOut);

//count the occurences of a substring in a string
Int CountSubstring(const std::string& str, const std::string& sub);

Int CountCSV(const std::string str);

//takes a string and a delimiter and returns a vector of string tokens
std::vector<std::string> Tokenize(const std::string str, const char delimiter);

//strips comma separated data, allocates memory and stores it
template <class theType>
Int StripCSV(const std::string str, theType** items)
{
  Int i;
  size_t loc;
  std::string value, rest;
  std::stringstream ss (std::stringstream::in | std::stringstream::out);

  Int nitems = CountCSV(str);		       
  theType* temp = new theType[nitems];

  rest = str;
  for(i = 0; i < nitems; i++){
    loc = rest.find(',');
    value = rest.substr(0, loc);

#if 0
    Int j;
    //NOTE: this normally shouldn't be necessary, left intact for 
    //strange compiler issues if they arise
    //strip off leading and trailing spaces
    for(j = 0; j < value.size(); j++){
      if(value[j] != ' '){
	value = value.substr(j, value.size()-j);
	break;
      }
    }
    for(j = value.size()-1; j >=0; j--){
      if(value[j] != ' '){
	value = value.substr(0, j+1);
	break;
      }
    }
#endif

    ss.clear();
    ss.str("");
    ss << value;
    ss >> temp[i];
    loc += 1;
    rest = rest.substr(loc, rest.size()- loc);
  }
  //just in case someone was goofy enough to pass allocated memory
  if(*items != NULL){
    std::cerr << "WARNING: StripCSV() ... passing a non-null pointer may result in a memory leak" << std::endl;;
    std::cerr << "WARNING: DON'T DO IT!!!" << std::endl;
    delete [] *items;
  }
  //point at new array
  *items = temp;
  return nitems;
}

//generates hash of a string
unsigned int stringHash(std::string str); 

#endif
