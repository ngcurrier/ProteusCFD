#include "h5layer.h"

void HDF_TurnOffErrorHandling()
{
  //turn off error handling
  H5Eset_auto1(NULL, NULL);

  return;
}

hid_t HDF_OpenFile(std::string filename, Int writable)
{
  hid_t fid = -1;
  unsigned int flags;

  //get readonly access
  if(writable == 0){
    flags = H5F_ACC_RDONLY;
  }
  else{
    flags = H5F_ACC_RDWR;
  }

  fid = H5Fopen(filename.c_str(), flags, H5P_DEFAULT);
  
  //if the file didn't open, create it
  if(fid < 0 && writable){
    fid = HDF_CreateFile(filename, writable);
  }

  return fid;
}

hid_t HDF_CreateFile(std::string filename, Int overwrite)
{
  hid_t fid;
  hid_t flags;
  
  if (overwrite){
    flags = H5F_ACC_TRUNC;
  }
  else{
    flags = H5F_ACC_EXCL;
  }
  fid = H5Fcreate(filename.c_str(),flags,H5P_DEFAULT,H5P_DEFAULT);
  if(fid < 0){
    std::cerr << "HDF_IO: Cannot create file " << filename << std::endl;
  }

  return fid;
}

hid_t HDF_CloseFile(hid_t fileId)
{
  hid_t err = H5Fclose(fileId);
  return (err);
}

hid_t HDF_GetDirectory(hid_t fileId, std::string directory)
{
  hid_t group = -1;

  //transfer plist to set up temp buffer
  hid_t xfer_plist = H5P_DEFAULT;

  std::string path = "/";

  group = H5Gopen(fileId, directory.c_str(), xfer_plist);
  //if group returns a failure iterate through the directory structure and build up the
  //correct group paths
  if(group < 0){
    std::vector<std::string> tokens = Tokenize(directory, '/');
    std::vector<std::string>::iterator it;
    for(it = tokens.begin(); it != tokens.end(); it++){
      path += *it + "/";
      group = H5Gopen(fileId, path.c_str(), xfer_plist);
      if(group >= 0) continue;
      group = H5Gcreate(fileId, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if(group < 0){
	std::cerr << "HDF_IO: Cannot create group " << path << " in HDF_GetDirectory()" << std::endl;
	return (-1);
      }
      H5Gclose(group);
    }
    group = H5Gopen2(fileId, path.c_str(), xfer_plist);
  }

  if(group < 0){
    std::cerr << "HDF_IO: Group could not be built " << path << " in HDF_GetDirectory()" << std::endl;
  }

  return group;
}

hid_t HDF_TestDirectory(hid_t fileId, std::string directory)
{
  herr_t status = H5Gget_objinfo(fileId, directory.c_str(), 0, NULL);
  if(status == 0){
    return true;
  }
  else{
    return false;
  }
}

hid_t HDF_WriteStringAttribute(hid_t fileId, std::string directory, std::string dataName, 
			       std::string attrName, const std::vector<std::string>& strings)
{
  hid_t err = 0;
  hid_t group, aspace, attribute, dset;
  hid_t a_type;
  Int stringSize = 0;
  std::string concat = "";
  hid_t xfer_plist = H5P_DEFAULT;

  //concatenate the string list to a comma separated list
  for(std::vector<std::string>::const_iterator it = strings.begin(); it != strings.end(); ++it){
    const std::string temp = *it;
    concat += temp;
    //add a comma unless this is the last value
    if(it != strings.end() - 1){
      concat += ",";
    }
  }

  //copy the type we need
  a_type = H5Tcopy(H5T_C_S1);

  //set the size to write
  H5Tset_size(a_type, concat.size());
  
  group = HDF_GetDirectory(fileId, directory);
  if(group < 0){
    std::cerr << "HDF_IO: Could not open group " << directory 
	      << " in HDF_WriteStringAttribute()" << std::endl;
  }

  //create the attribute space
  aspace = H5Screate(H5S_SCALAR);
  if(aspace < 0){
    std::cerr << "HDF_IO: Could not create attribute space in HDF_WriteStringAttribute()" << std::endl;
  }

  dset = H5Dopen(group, dataName.c_str(), H5P_DEFAULT);
  if(dset < 0){
    dset = H5Dcreate(group, dataName.c_str(), a_type, aspace, xfer_plist, xfer_plist, xfer_plist);
    if(dset < 0){
      std::cerr << "HDF_IO: Cannot create dataset " << directory+dataName 
		<< " in HDF_WriteStringAttribute()!" << std::endl;
      H5Gclose(group);
      return (-1);
    }
  }
  
  attribute = H5Acreate(dset, attrName.c_str(), a_type, aspace, H5P_DEFAULT, H5P_DEFAULT);
  if(attribute < 0){
    //if fail, delete then recreate it, not elegant but it works
    err = H5Adelete(dset, attrName.c_str());
    if(err < 0){
      std::cerr << "HDF_IO: Could not delete attribute " << directory+dataName 
		<< " - " << attrName << " failed in HDF_WriteStringAttribute()" << std::endl;
    }
    attribute = H5Acreate(dset, attrName.c_str(), a_type, aspace, H5P_DEFAULT, H5P_DEFAULT);
    if(attribute < 0){
      std::cerr << "HDF_IO: Could not create string attribute" << directory+dataName 
		<< " - " << attrName << " failed in HDF_WriteStringAttribute()" << std::endl;
    }
  }

  err = H5Awrite(attribute, a_type, concat.c_str());
  if(err < 0){
    std::cerr << "HDF_IO: Writing string attribute " << directory+dataName << " - " << attrName 
	      << " failed in HDF_WriteStringAttribute()" << std::endl;
  }
  
  //free the resources
  H5Gclose(group);
  H5Dclose(dset);
  H5Aclose(attribute);
  H5Sclose(aspace);
  H5Tclose(a_type);

  return err;
}

hid_t HDF_ReadStringAttribute(hid_t fileId, std::string directory, std::string dataName, 
			      std::string attrName, std::vector<std::string>& strings)
{
  hid_t err = 0;
  hid_t group, aspace, attribute, dset;
  hid_t xfer_plist = H5P_DEFAULT;

  //get the group
  group = HDF_GetDirectory(fileId, directory);

  //open the dataset
  dset = H5Dopen(group, dataName.c_str(), xfer_plist);
  if(dset < 0){
    H5Gclose(group);
    std::cerr << "HDF_IO: Cannot open dataset " << directory+dataName 
	      << " in HDF_ReadStringAttribute()" << std::endl;
    return (-1);
  }

  //open attribute
  attribute = H5Aopen(dset, attrName.c_str(), xfer_plist);
  if(attribute < 0){
    H5Gclose(group);
    H5Dclose(dset);
    std::cerr << "HDF_IO: Cannot open attribute " << directory+dataName << " - "
	      << attrName << " in HDF_ReadStringAttribute() " << std::endl;
    return (-1);
  }

  //get the data type
  hid_t atype2 = H5Aget_type(attribute);
  hsize_t size = H5Tget_size(atype2);
  char* concat = new char[size+1];
  
  //get the dataspace
  aspace = H5Aget_space(attribute);

  //read the attribute
  err = H5Aread(attribute, atype2, concat);
  if(err < 0){
    std::cerr << "HDF_IO: Read failed for attribute " << directory+dataName << " - "
	      << attrName << " in HDF_ReadStringAttribute() " << std::endl;
  }

  //add null terminator so string initialization doesn't overrun buffer
  concat[size] = '\0';
  std::string sconcat = concat;
  strings = Tokenize(sconcat, ',');

  //free the resources
  H5Gclose(group);
  H5Dclose(dset);
  H5Aclose(attribute);
  H5Sclose(aspace);
  H5Tclose(atype2);

  delete [] concat;

  return err;
}

//This function takes an hdf file handle and will return a vector of paths to
//all the datasets in that file
void HDF_GetAllDatasetsInFile(hid_t fileId, std::vector<std::string>& pathList)
{
  std::string groupname = "/";
  hid_t group = HDF_GetDirectory(fileId, groupname);
  RecurseDirectory tempList;
  HDF_GetAllDatasetsBelowGroupInternal(group, tempList);
  H5Gclose(group);

  for(std::vector<std::pair<std::string, Int> >::iterator it = tempList.list.begin(); it != tempList.list.end(); ++it){
    std::pair<std::string, Int>& name = *it;
    std::string namestring = name.first;
    std::string stripstring;
    Int namelevel = name.second;
    std::string groupcheck = "GROUP-";
    std::string dsetcheck = "DSET-";
    std::string concat;
    //check for a group
    if(namestring.substr(0,6) == groupcheck){
      stripstring = namestring.substr(6, std::string::npos);
    }
    //check for a dataset - if found, build up the path to here
    else if(namestring.substr(0,5) == dsetcheck){
      stripstring = namestring.substr(5, std::string::npos);
      concat = "/" + stripstring;
      //if not at the root level, build up the path
      Int levelit = namelevel;
      if(namelevel != 1){
	for(std::vector<std::pair<std::string, Int> >::iterator itn = it; itn != tempList.list.begin(); --itn){
	  std::pair<std::string, Int>& namein = *itn;
	  std::string namestringin = namein.first;
	  Int namelevelin = namein.second;
	  if(namelevelin == levelit -1){
	    //check for a group
	    if(namestringin.substr(0,6) == groupcheck){
	      namestringin = namestringin.substr(6, std::string::npos);
	      concat = "/" + namestringin + concat;
	      levelit--;
	    }
	    //check for a dataset
	    else if(namestringin.substr(0,5) == dsetcheck){
	      namestringin = namestringin.substr(5, std::string::npos);
	    }
	    else{
	      std::cerr << "Type found in list not identified in HDF_GetAllDatasetsInFile()" << std::endl;
	    }
	  }
	  if(levelit == 1){
	    continue;
	  }
	}
      }
      pathList.push_back(concat);
    }
    else{
      std::cerr << "Type found in list not identified in HDF_GetAllDatasetsInFile()" << std::endl;
    }
  }
  return;
}

//helper function for HDF_GetAllDatasetsInFile()
void HDF_GetAllDatasetsBelowGroupInternal(hid_t group, RecurseDirectory& pathList)
{
  herr_t err = 0;
  hsize_t* idx = NULL;
  //make sure that we increment the level b/c this function is recursive
  pathList.level++;
  err = H5Literate(group, H5_INDEX_NAME, H5_ITER_INC, idx, HDF_AddToList, (void*)&pathList);
  //when function returns move the level counter back up again
  pathList.level--;
}

//helper function for HDF_GetAllDatasetsInFile() and HDF_GetAllDatasetsBelowGroupInternal()
herr_t HDF_AddToList(hid_t group, const char* name, const H5L_info_t* info, void* op_data)
{
  H5G_stat_t statbuf;
  RecurseDirectory& pathList = *static_cast<RecurseDirectory*>(op_data);
  hid_t xfer_plist = H5P_DEFAULT;
  hid_t group2 = 0;
  std::pair<std::string, Int> paired;

  H5Gget_objinfo(group, name, false, &statbuf);
  
  switch(statbuf.type) {
  case H5G_GROUP:
    paired.first = name;
    paired.first = "GROUP-" + paired.first;
    paired.second = pathList.level;
    pathList.list.push_back(paired);
    group2 = H5Gopen(group, name, xfer_plist);
    HDF_GetAllDatasetsBelowGroupInternal(group2, pathList);
    H5Gclose(group2);
    break;
  case H5G_DATASET:
    paired.first = name;
    paired.first = "DSET-" + paired.first;
    paired.second = pathList.level;
    pathList.list.push_back(paired);
    break;
  case H5G_TYPE:
    break;
  default:
    std::cerr << "HDF_AddToList() Unable to identify an object" << std::endl;
    return (-1);
  }


  //from this function, returning 0 causes iteration to continue
  //                    returning + causes iterator to short circuit positive success
  //                    returning - causes iterator to short circuit negative failure
  return 0;
}

