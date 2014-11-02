#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "portFileio.h"
#include "portFunctions.h"

void function_(FUNC_ARGS)
{
  int err = 0;
  double* bounds = new double[*ndv * 2];
  double* grad = new double[*ndv];
  int* dvType = new int[*ndv];

  double* trashx = new double[*ndv];

  //trash the x values since we don't won't to kill the ones passed in
  ReadDesignFile(casename, ndv, trashx, bounds, f, grad, dvType);
  //write out design file for movement using new x values
  WriteDesignFile(casename, *ndv, x, bounds, *f, grad, dvType);

  //perturb points according to x
  MoveMesh(nprocs);

#if 0
  //this section used to make movies only
  std::string recomp = "./urecomp moved >> recomp.out";
  std::string mvbase = "mv moved-0.vtk moved-";
  std::stringstream ss;
  ss << *nf;
  mvbase += (ss.str() + ".vtk");
    
  //recomp mesh
  system(recomp.c_str());
  //mv mesh to place in series
  system(mvbase.c_str());
#endif

  //run flow solver
  ComputeObjectiveFunction(nprocs);

  //Read design file again to retrieve the "new" objective function value
  ReadDesignFile(casename, ndv, trashx, bounds, f, grad, dvType);

  if(err){
    *nf = 0;
  }
  
  urparm[0] = *f;

  WriteDesignFile(casename, *ndv, x, bounds, *f, grad, dvType);

  WriteLog(casename, *ndv, x, grad, *f);

  delete [] bounds;
  delete [] grad;
  delete [] dvType;
  delete [] trashx;

  return;
}

void gradient_(GRAD_ARGS)
{
  int err = 0;  
  
  double* bounds = new double[(*ndv) * 2];
  int* dvType = new int[*ndv];
  double f = 0.0;

  double* trashx = new double[*ndv];

  ReadDesignFile(casename, ndv, trashx, bounds, &f, grad, dvType);
  //write out design file for movement using new x values
  WriteDesignFile(casename, *ndv, x, bounds, f, grad, dvType);

  //move points in case betas have changed
  MoveMesh(nprocs);
  //compute grid sensitivities
  ComputeMeshSensitivities(nprocs);
  //compute gradients
  ComputeGradient(nprocs, GRAD_METHOD);
  //read gradients from file
  ReadDesignFile(casename, ndv, trashx, bounds, &f, grad, dvType);
  //write design file out with most recent data
  WriteDesignFile(casename, *ndv, x, bounds, f, grad, dvType);
  
  //if error occurs let port know
  if(err){
    *nf = 0;
  }

  WriteLog(casename, *ndv, x, grad, f);

  delete [] trashx;
  delete [] bounds;
  delete [] dvType;
}


int GetNptsPointFile(std::string filename)
{
  int npts;

  std::ifstream fin;

  fin.open(filename.c_str());
  if(fin.is_open()){
    fin >> npts;
  }
  else{
    std::cerr << "POINT READFILE: Cannot open point file --> " << filename << std::endl;
    return(-1);
  }
  fin.close();

  return npts;
}

void ReadPointMoveFile(std::string filename, int npts, int ndv, 
		       int* pts, double* xyz, double* dxyz, int* betaid)
{
  int i, j;
  int tnpts;
  std::ifstream fin;
  fin.open(filename.c_str());
  if(!fin.is_open()){
    std::cerr << "WARNING: opening file failed --> " << filename << std::endl;
    return;
  }

  fin >> tnpts;
  if(tnpts != npts){
    std::cerr << "WARNING: npts does not match points expected ReadPointMoveFile()" << std::endl;
    return;
  }
  for(i = 0; i < npts; i++){
    fin >> pts[i];
    for(j = 0; j < 3; j++){
      fin >> xyz[i*3 + j];
    }
    for(j = 0; j < 3; j++){
      fin >> dxyz[i*3 + j];
    }
  }
  if(betaid != NULL){
    for(i = 0; i < npts; i++){
      fin >> betaid[i];
    }
  }


  fin.close();
  return;
}

void WritePointMoveFile(std::string filename, int npts, int ndv, 
			int* pts, double* xyz, double* dxyz, int* betaid)
{
  int i, j;
  std::ofstream fout;
  fout.open(filename.c_str());
  if(!fout.is_open()){
    std::cerr << "WARNING: opening file failed --> " << filename << std::endl;
    return;
  }

  fout.setf(std::ios::scientific);
  fout.precision(15);

  fout << npts << std::endl;
  for(i = 0; i < npts; i++){
    fout << pts[i] << " ";
    for(j = 0; j < 3; j++){
      fout << xyz[i*3 + j] << " ";
    }
    for(j = 0; j < 3; j++){
      fout << dxyz[i*3 + j] << " ";
    }
    fout << std::endl;
  }
  for(i = 0; i < npts; i++){
    fout << betaid[i] << std::endl;
  }

  fout << std::endl << std::endl;
  fout << "#This file format is as follows" << std::endl;
  fout << "#npts" << std::endl;
  fout << "#ptid      x y z	deltax	deltay	deltaz" << std::endl;
  fout << "#ptid      x y z	deltax	deltay	deltaz" << std::endl;
  fout << "#beta id associated with point 0" << std::endl;
  fout << "#beta id associated with point 1" << std::endl;

  fout.close();
  return;
}

void ComputeObjectiveFunction(int np)
{
  std::stringstream ss;
  std::string command;
  ss << np;
  command = "mpirun -np " + ss.str() + " ./ucs " + casename + " 1";
  if(locationSwitch == 0){
    system(command.c_str());
  }
  else if(locationSwitch == 1){
    ExecuteRemoteAndWait("jobsub.py -n " + ss.str() + " -j PortObj \"./ucs " + casename + " 1\"");
  }
  return;
}

void MoveMesh(int np)
{
  std::stringstream ss;
  std::string command;
  ss << np;

  //copy original (non-moved) partitions down
  CopyBasePartitionFilesDown(np);

  //now perturb the copies
  command = "mpirun -np " + ss.str() + " ./ucs " + casename + " 5";
  if(locationSwitch == 0){
    system(command.c_str());
  }
  else if(locationSwitch == 1){
    ExecuteRemoteAndWait("jobsub.py -n " + ss.str() + " -j PortMove \"./ucs " + casename + " 5\"");
  }

  return;
}

void SaveOriginalMesh(int np)
{
  std::stringstream ss;
  std::string command;
  ss << np;
  std::cout << "SAVING ORIGINAL MESH TO " + pathname + "originalMesh.* \n";

  for(int i = 0; i < np; i++){
    ss.str("");
    ss.clear();
    ss << i;
    command = "cp " + casename + "." + ss.str() + ".h5 " + pathname + "originalMesh." + ss.str() + ".h5";
    system(command.c_str());
  }
}

void SaveOriginalDesignFile()
{
  std::cout << "SAVING ORIGINAL DESIGN FILE TO " + pathname + "original.design \n";
  std::string command = "cp " + casename + ".design " + pathname + "original.design";
  system(command.c_str());
}

void CopyBasePartitionFilesDown(int np)
{
  int i;
  std::stringstream ss;
  std::string command;
  for(i = 0; i < np; i++){
    ss.str("");
    ss.clear();
    ss << i;
    command = "cp " + pathname + "originalMesh." + ss.str() + ".h5 " + casename + "." + ss.str() + ".h5";
    system(command.c_str());
  }
  return;
}

void ComputeMeshSensitivities(int np)
{
  std::stringstream ss;
  std::string command;
  ss << np;

  command = "mpirun -np " + ss.str() + " ./ucs " + casename + " 6";
  if(locationSwitch == 0){
    system(command.c_str());
  }
  else if (locationSwitch == 1){
    ExecuteRemoteAndWait("jobsub.py -n " + ss.str() + " -j PortSens \"./ucs " + casename + " 6\"");
  }
  return;
}

void ComputeGradient(int np, int method)
{
  std::stringstream ss1, ss2;
  std::string command;
  ss1 << np;
  
  command = "mpirun -np " + ss1.str() + " ./ucs " + casename + " "; 

  //now add the correct method of gradient computation
  //2 - direct
  //3 - adjoint
  //4 - CTSE
  //7 - FD
  ss2.clear();
  ss2.str("");
  ss2 << method;
  command += ss2.str();
  if(locationSwitch == 0){
    system(command.c_str());
  }
  else if(locationSwitch == 1){
    ExecuteRemoteAndWait("jobsub.py -n " + ss1.str() + " -j PortGrad \"./ucs " + casename + " " + ss2.str() + "\"");
  }

  return;
}

std::string GetOutputFromCommand(std::string cmd)
{
  char buffer[MAX_BUFFER];
  std::string data;
  FILE* stream;
  
  cmd += " 2>&1";

  stream = popen(cmd.c_str(), "r");
  while ( fgets(buffer, MAX_BUFFER, stream) != NULL){
    data.append(buffer);
  }

  pclose(stream);

  return data;
}


void ExecuteRemoteAndWait(std::string cmd)
{
  std::cout << "COMMAND: " << cmd << std::endl;

  //stream the output from our queue submission to a string
  std::string data = GetOutputFromCommand(cmd);
  std::stringstream ss;
  //push the output into a string stream
  ss << data;

  //this string is immediately before the job id which we need
  std::string jobFlag = "Your job";
  size_t pos;
  //look for the job id indicator
  pos = data.find(jobFlag.c_str());

  //this string is present when there is either a failure or a completion
  std::string finishedFlag = "Following jobs do not exist:";

  //seek to the location of the jobid and stream it out
  ss.seekg(pos+jobFlag.size());
  int jobid;
  ss >> jobid;

  std::cout << "JOB ID: " << jobid << std::endl;

  //clear the stream and create the command to check for the job status
  ss.clear();
  ss.str("");
  ss << jobid;
  cmd = "qstat -j " + ss.str(); 

  bool tmp;

  //this loop continues to poll the queueing system to look for a
  //finished job, this is frought with peril but the best 
  //I have come up with at this point
  do{
    //query the queueing system
    data = GetOutputFromCommand(cmd);

    //look for our string we know indicates either a failure 
    //or a finished job
    pos = data.find(finishedFlag);

    //wait for 30 seconds before polling again
    sleep(30);

    tmp = (pos == std::string::npos);

  } while(tmp);

  std::cout << "THE JOB " << jobid << "HAS FINISHED RUNNING" << std::endl;

  return;
}

