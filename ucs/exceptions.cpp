#include "exceptions.h"

std::string GetStack(std::string errorstr)
{				    
  void* message[500];					    
  char** strings;					    
  int nptrs;						    
  nptrs = backtrace(message, 500);			    
  strings = backtrace_symbols(message, nptrs);		    
  std::stringstream ERROR_MESSAGE;			    
  for (int II = 0; II < nptrs; II++) {			    
    ERROR_MESSAGE << std::endl << strings[II];		    
  }									
  ERROR_MESSAGE << std::endl << __FILE__<<": " << __LINE__ << ": " << errorstr; 
  return ERROR_MESSAGE.str();
};
