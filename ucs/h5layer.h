#ifndef H5_LAYER_H__
#define H5_LAYER_H__

#include <hdf5.h>
#include <typeinfo>
#include "general.h"
#include "strings_util.h"
#include <iostream>
#include <string>
#include <vector>

//helper struct for HDF_GetAllDatasetsInFile()
class RecurseDirectory
{
public:
  RecurseDirectory()
  {
    level = 0;
  };
  Int level;
  std::vector<std::pair<std::string, Int> > list;
private:
};

//returns the hdf type for the passed in value
template <class Type>
hid_t HDF_GetType(Type a)
{
  if (typeid(a) == typeid(int)) return H5T_NATIVE_INT;
  else if (typeid(a) == typeid(unsigned int)) return H5T_NATIVE_UINT;
  else if (typeid(a) == typeid(short int)) return H5T_NATIVE_SHORT;
  else if (typeid(a) == typeid(unsigned short int)) return H5T_NATIVE_USHORT;
  else if (typeid(a) == typeid(long int)) return H5T_NATIVE_LONG;
  else if (typeid(a) == typeid(unsigned long int)) return H5T_NATIVE_ULONG;
  else if (typeid(a) == typeid(long long int)) return H5T_NATIVE_LLONG;
  else if (typeid(a) == typeid(unsigned long long int)) return H5T_NATIVE_ULLONG;
  else if (typeid(a) == typeid(float)) return H5T_NATIVE_FLOAT;
  else if (typeid(a) == typeid(double)) return H5T_NATIVE_DOUBLE;
  else if (typeid(a) == typeid(long double)) return H5T_NATIVE_LDOUBLE;
  else if (typeid(a) == typeid(char)) return H5T_NATIVE_CHAR;
  else if (typeid(a) == typeid(unsigned char)) return H5T_NATIVE_UCHAR;
  else
    {
      std::cerr << "HDF: unable to map type!!\n" << std::endl;
    }

  return(-1);
};

//Turns off HDF specific error handling, internal error handling still on
void HDF_TurnOffErrorHandling();
//opens an hdf file, if writable is true file will also be created/writable if it doesn't exist
hid_t HDF_OpenFile(std::string filename, Int writable);
//close hdf file by its handle
hid_t HDF_CloseFile(hid_t fileId);
//create an hdf file
hid_t HDF_CreateFile(std::string filename, Int overwrite);
//will parse a directory structure i.e. /root/rootplusone/data/ and open that group
hid_t HDF_GetDirectory(hid_t fileId, std::string directory);
//will test if a certain group exists and return true if it does, false if not
hid_t HDF_TestDirectory(hid_t fileId, std::string directory);

//reads a single scalar value from directory and dataset given
template <class Type>
hid_t HDF_ReadScalar(hid_t fileId, std::string directory, std::string dataName, Type* scalar);

//reads an array from directory and dataset given
//if n == -1 and data == NULL this routine will simply return the size of the array in n
template <class Type>
hid_t HDF_ReadArray(hid_t fileId, std::string directory, std::string dataName, Type** data, Int* n, 
		    Int* ncols = NULL);

//write a single scalar to directory and dataset given
template <class Type>
hid_t HDF_WriteScalar(hid_t fileId, std::string directory, std::string dataName, Type* scalar);

//write an array to directory and dataset given
template <class Type>
hid_t HDF_WriteArray(hid_t fileId, std::string directory, std::string dataName, Type* data, Int nrows, 
		    Int ncols = 1);

//write an array to an attribute
template <class Type>
hid_t HDF_WriteArrayAttribute(hid_t fileId, std::string directory, std::string dataName, 
			      std::string attrName, std::vector<Type>& attrArray);
//read an array from an attribute
template <class Type>
hid_t HDF_ReadArrayAttribute(hid_t fileId, std::string directory, std::string dataName, 
			     std::string attrName, std::vector<Type>& attrArray);

//write a list of strings to a comma delimited attribute
hid_t HDF_WriteStringAttribute(hid_t fileId, std::string directory, std::string dataName, 
			       std::string attrName, const std::vector<std::string>& strings);

//read a list of strings from a comma delimited attribute
hid_t HDF_ReadStringAttribute(hid_t fileId, std::string directory, std::string dataName, 
			      std::string attrName, std::vector<std::string>& strings);

//function to parse the directory structure for every dataset in an hdf file
void HDF_GetAllDatasetsInFile(hid_t fileId, std::vector<std::string>& pathList);
//helper function for HDF_GetAllDatasetsInFile()
void HDF_GetAllDatasetsBelowGroupInternal(hid_t group, RecurseDirectory& pathList);
//helper function for HDF_GetAllDatasetsInFile()
herr_t HDF_AddToList(hid_t group, const char* name, const H5L_info_t* info, void* op_data);

//include implementations
#include "h5layer.tcc"

#endif
