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
#include <unistd.h>

//returns a stack trace
std::string GetStack(std::string errorstr);

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

    std::string stack = GetStack(message);
    
    std::cerr << "Aborting (Process Id - " << rank << "): " << message << stack << std::endl << std::flush;
    std::cout << "Aborting (Process Id - " << rank << "): " << message << stack << std::endl << std::flush;
    
    fout.open(filename.c_str());
    fout << "Aborting (Process Id - " << rank << "): " << message << stack << std::endl;
    fout.close();
    sleep(1); //wait 1 sec, this sometimes prevents a race condition with the file I/O

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
