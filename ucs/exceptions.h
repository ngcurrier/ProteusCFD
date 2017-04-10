#ifndef ABORT_H__
#define ABORT_H__

#include <mpi.h>
#include "general.h"
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <stdexcept>
#include <execinfo.h>

#define GetStack(errorstr){                                 \
    void* message[500];					    \
    char** strings;					    \
    int nptrs;                                              \
    nptrs = backtrace(message, 500);                        \
    strings = backtrace_symbols(message, nptrs);	    \
    std::stringstream ERROR_MESSAGE;			    \
    for (int II = 0; II < nptrs; II++) {		    \
      ERROR_MESSAGE << std::endl << strings[II];	    \
    }							    \
    ERROR_MESSAGE << std::endl << __FILE__<<": " << __LINE__ << ": " << errorstr; \
    std::cout << ERROR_MESSAGE.str() << std::endl;		                  \
  }

class Abort
{
public:
  std::string rootDirectory;
  std::string softAbortMessage;
  bool softAbort;
  
  Abort(){
    rootDirectory = "";
    softAbortMessage = "";
    softAbort = false;
  };

  void operator << (const std::string message)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::ofstream fout;
    std::string filename = rootDirectory + "Abort.Log";

    fout.open(filename.c_str());
    fout << "Aborting (Process Id - " << rank << "): " << message << std::endl;
    fout.close();

    MPI_Abort(MPI_COMM_WORLD, -999);
  };

  void SetSoftAbort(const std::string message)
  {
    softAbortMessage += message +"\n";
    softAbort = true;
  };

  void CheckForSoftAbort(const std::string message)
  {
    if(softAbort){
      std::string tmp = message + softAbortMessage;
      *this << tmp;
    }
  };

};

static Abort Abort;

#endif
