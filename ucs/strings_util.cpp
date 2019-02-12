#include "strings_util.h"

std::string UpperToLower(std::string str) {
  unsigned int i;
  for (i = 0; i < str.length(); i++){
    //do some bit twiddling
    //another option is using the c std function
    //tolower
    //if(str[i] >= 0x41 && str[i] <= 0x5A){
    //  str[i] = str[i] + 0x20;
    //}
    str[i] = tolower(str[i]);
  }
  return str;
}

std::string LowerToUpper(std::string str) {
  unsigned int i;
  for (i = 0; i < str.length(); i++){
    str[i] = toupper(str[i]);
  }
  return str;
}

//returns string between first set of left, right strings... 
//if pos is given search starts from there
Int GetStringBetween(const std::string left, const std::string right, const std::string strIn, std::string& strOut, size_t pos)
{
  Int err = 0;
  size_t loc1,loc2;
  loc1 = strIn.find(left, pos);
  loc2 = strIn.find(right, loc1+1);
  if(loc1 == std::string::npos || loc2 == std::string::npos){
    strOut = "";
    err = 1;
  }
  else{
    //eliminate the delimiters on each end
    loc1 += left.size();
    strOut = strIn.substr(loc1, loc2-loc1);
  }
  return err;
}

//returns string with section between first set of left, right strings removed... 
//if pos is given search starts from there
Int StripStringBetween(const std::string left, const std::string right, const std::string strIn, std::string& strOut, size_t pos)
{
  Int err = 0;
  size_t loc1,loc2;
  loc1 = strIn.find(left, pos);
  loc2 = strIn.find(right, loc1+1);
  if(loc1 == std::string::npos || loc2 == std::string::npos){
    strOut = "";
    err = 1;
  }
  else{
    //eliminate the delimiters on each end
    strOut = strIn.substr(0, loc1) + strIn.substr(loc2+right.size());
  }
  return err;
}

//strips the string 'stripped' from base and returns the result in strOut
//return value is the number of times the string 'stripped' was found
Int StripFromString(std::string base, std::string stripped, std::string& strOut)
{
  Int nitems = 0;
  size_t loc1;
  size_t s = stripped.size();

  //initialize return string to the input
  strOut = base;

  //look for the string to strip out
  loc1 = strOut.find(stripped);
  while(loc1 != std::string::npos){
    //eliminate the string 'stripped'
    strOut = strOut.substr(0, loc1) + strOut.substr(loc1+s);
    nitems++;
    loc1 = strOut.find(stripped);
  }

  return nitems;
}

//count the occurences of a substring in a string
Int CountSubstring(const std::string& str, const std::string& sub)
{
  if (sub.length() == 0) return 0;
  Int count = 0;
  for (size_t offset = str.find(sub); offset != std::string::npos;
       offset = str.find(sub, offset + sub.length())){
    ++count;
  }
  return count;
}


//Assumes form of 1,2,3,etc.
//will return number of commas found plus one for last value
Int CountCSV(const std::string str)
{
  Int count;
  size_t loc;
  std::string temp;

  count = 0;
  temp = str;
  loc = temp.find(',');
  while(loc != std::string::npos){
    count++;
    temp = temp.substr(loc+1, temp.size()-1);
    loc = temp.find(',');
  }

  return count+1;
}

//Removes whitespace from a string
void RemoveWhitespace(std::string& s)
{
  std::string out;
  for(size_t i = 0; i < s.size(); ++i){
    if(s[i] != ' '){
      out.push_back(s[i]);
    }
  }
  s = out;
}

//You should not pass in ' ' as the split feature, use tokenize if you'd like to do that here
std::vector<std::string> Split(const std::string str, const std::string splitOn)
{
  std::vector<std::string> v;
  size_t loc = str.find(splitOn);
  size_t l = splitOn.size();
  if(loc != std::string::npos){
    v.push_back(str.substr(0,loc));
    v.push_back(str.substr(loc+l,str.size()-loc)); 
  }
  else{
    std::stringstream ss;
    ss << "Split() could not split on values passed : " << str << " -- " << splitOn << std::endl;
    throw std::invalid_argument(ss.str());
  }
		     
  return v;
}

void Replace(std::string& in, const std::string replaceThis, const std::string withThis)
{
  std::string tmp("");
  size_t loc = in.find(replaceThis);
  size_t l = replaceThis.size();
  if(loc != std::string::npos){
    tmp += in.substr(0, loc);
    tmp += withThis;
    tmp += in.substr(loc+l, in.size() - loc);
    in = tmp;
  }
  else{
    std::stringstream ss;
    ss << "Replace() could not split on values passed : " << in << " -- " << replaceThis << std::endl;
    throw std::invalid_argument(ss.str());
  }
}

std::vector<std::string> Tokenize(const std::string str, const char delimiter)
{
  std::vector<std::string> tokens;
  std::string token;
  std::string workingStr = str;
  size_t loc1;

  //get rid of the first character if it is a delimiter
  if(workingStr[0] == delimiter){
    workingStr = workingStr.substr(1);
  }
  
  //place a delimiter at the end of the string if it is not there
  //the back() member function is not portable yet, sigh...
  size_t last = workingStr.size()-1;
  if(workingStr[last] != delimiter){
    workingStr += delimiter;
  } 

  while(workingStr.size() != 0){
    loc1 = workingStr.find(delimiter, 0);
    token = workingStr.substr(0, loc1);
    workingStr = workingStr.substr(loc1+1);
    //do not return empty tokens
    bool allwhite = true;
    for(size_t i = 0; i < token.size(); ++i){
      if(token[i] != ' ') allwhite = false;
    }
    if(allwhite) continue;
    //attempt to eliminate any leading spaces
    size_t loc2 = 0;
    while(loc2 != std::string::npos){
      loc2 = token.find(' ');
      if(loc2 != std::string::npos && loc2 == 0){
	token = token.substr(loc2+1);
	continue;
      }
      break;
    }
    tokens.push_back(token);
  }
  return tokens;
}

std::vector<double> TokenizeToDoubles(const std::string str, const char delimiter){
  std::stringstream ss;
  std::vector<std::string> sv = Tokenize(str, delimiter);
  std::vector<double> v(sv.size(), 0.0);
  for(size_t i = 0; i < sv.size(); ++i){
    ss << sv[i];
    ss >> v[i];
    ss.clear();
  }
  return v;
}

unsigned int stringHash(std::string str){
  unsigned int hash = 0;
  for(size_t i = 0; i < str.size(); ++i) 
    hash = 65599 * hash + str[i];
  return hash ^ (hash >> 16);
}



//Remove surrounding whitespace from a std::string.
// @param s The string to be modified.
// @param t The set of characters to delete from each end
// of the string.
// @return The same string passed in as a parameter reference.
//
std::string& trim(std::string& s, const char* t)
{
	s.erase(0, s.find_first_not_of(t));
	s.erase(s.find_last_not_of(t) + 1);
	return s;
}
