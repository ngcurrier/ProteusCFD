#ifndef PORT_FILEIO_H__
#define PORT_FILEIO_H__

//NOTE: both of these functions assume memory is previously allocated to appropriate size
int ReadDesignFile(std::string casename, int* ndv, double* x, double* bounds, double* f, 
		   double* grad, int* dvType);
int WriteDesignFile(std::string casename, int ndv, double* x, double* bounds, double f, 
		    double* grad, int* dvType);
int WriteGradToDesignFile(std::string casename, double* grad);
int WriteFunctionToDesignFile(std::string casename, double f);
//returns number of design variables in the file
int GetNdvDesignFile(std::string casename);

//clear log file
void ClearLog(std::string casename);
//write state to log file
void WriteLog(std::string casename, int ndv, double* x, double* grad, double f);

#endif
