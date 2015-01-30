template <class Type>
hid_t HDF_ReadScalar(hid_t fileId, std::string directory, std::string dataName, Type* scalar)
{
  hid_t xfer_plist = H5P_DEFAULT;
  hid_t group, dset, space;
  hid_t t_id, ms_id;
  herr_t err = 0;

  //open the group
  group = HDF_GetDirectory(fileId, directory);

  //open the dataset
  dset = H5Dopen(group, dataName.c_str(), xfer_plist);
  if(dset < 0){
    H5Gclose(group);
    std::cerr << "HDF_IO: Cannot open dataset " << dataName << std::endl;
    return (-1);
  }
  //get type id from dataset
  t_id = H5Dget_type(dset);
  //det dataspace id from dataset
  space = H5Dget_space(dset);
  //create hdf memspace
  ms_id = H5Screate(H5S_SCALAR);

  //get type that was passed in
  Type check = 0;
  t_id = HDF_GetType(check);

  //read in the scalar
  err = H5Dread(dset, t_id, ms_id, space, H5P_DEFAULT, scalar);
  if(err < 0){
    std::cerr << "HDF_IO: Read failed from " << directory+dataName << " in HDF_ReadScalar()!" << std::endl;
  }

  //clean up
  H5Dclose(dset);
  H5Sclose(space);
  H5Gclose(group);
  H5Tclose(t_id);

  return err;
}

template <class Type>
hid_t HDF_ReadArray(hid_t fileId, std::string directory, std::string dataName, Type** data, Int* n, Int* ncols)
{
  hid_t group = -1;
  hid_t d_id = -1;
  hid_t t_id = -1;
  hid_t ht_id = -1;
  hid_t ds_id = -1;
  hid_t ms_id = -1;
  
  hsize_t rank;
  hsize_t dims[2];
  hsize_t sstride[2], count[2], block[2], start[2];
  
  //transfer plist to set up temp buffer
  hid_t xfer_plist = H5P_DEFAULT;
  herr_t err = 0;
  
  //open the group
  group = HDF_GetDirectory(fileId, directory);
  
  //open the dataset
  d_id = H5Dopen(group, dataName.c_str(), xfer_plist);
  if(d_id < 0){
    H5Gclose(group);
    std::cerr << "HDF_IO: Cannot open dataset " << dataName << std::endl;
    return (-1);
  }

  //get type from dataset
  ht_id = H5Dget_type(d_id);
  if(ht_id < 0){
    std::cerr << "HDF_IO: Cannot get type from dataset in HDF_ReadArray()" << std::endl;
    H5Gclose(group);
    H5Dclose(d_id);
    return(-1);
  }

  //det dataspace from dataset
  ds_id = H5Dget_space(d_id);
  if(ds_id < 0){
    std::cerr << "HDF_IO: Cannot get dataset from dataset in HDF_ReadArray()" << std::endl;
  }
  
  //query the dataspace for the dimensions
  rank = H5Sget_simple_extent_ndims(ds_id);
  if(rank != 1 && rank != 2){
    H5Gclose(group);
    H5Sclose(ds_id);
    std::cerr << "WARNING: Rank == 1 or Rank==2 expected.  Dataset rank = " << (int)rank << std::endl;
    return(-1);
  }
  err = H5Sget_simple_extent_dims(ds_id,dims,NULL);

  //if the input size is requested i.e. n == -1 set it
  //otherwise make sure they match
  if(*n == -1){
    *n = dims[0];
  }
  else{
    if((hsize_t)(*n) != dims[0]){
      std::cerr << "WARNING: size of array expected does NOT match in HDF_ReadArray()";
      std::cerr << "\texpected " << dims[0] << " passed in " << *n << std::endl;
    }
  }

  if(ncols){
    if(*ncols == -1){
      *ncols = dims[1];
    }
    else{
      if((hsize_t)(*ncols) != dims[1]){
	std::cerr << "WARNING: column count of array expected does NOT match!" << std::endl;
      }
    }
  }

  if(*data == NULL){
    //return after getting number of cols and rows for allocations
    H5Tclose(ht_id);
    H5Dclose(d_id);
    H5Sclose(ds_id);
    H5Gclose(group);
    return (0);
  }

  //get type that was passed in
  Type check = 0;
  t_id = HDF_GetType(check);
  if(t_id < 0){
    std::cerr << "HDF_IO: Cannot get data type" << std::endl;
  }
  
  
  //define the hyperslab
  start[0] = 0;
  sstride[0] = 1;
  count[0] = dims[0];
  block[0] = 1;
  
  if (rank == 2){
    start[1] = 0;
    sstride[1] = 1;
    count[1] = dims[1];
    block[1] = 1;
  }
  
  err = H5Sselect_hyperslab(ds_id,H5S_SELECT_SET,start,sstride,count,block);
  //create hdf memspace
  ms_id = H5Screate_simple(rank, dims, NULL);
  //define memory hyperslab
  H5Sselect_hyperslab(ms_id,H5S_SELECT_SET,start,sstride,count,block);
  
  //in case the type in the file is different that what we request
  //increase the buffer size over 1MB so conversion can proceed at a decent speed
  //if needed
  if(t_id != ht_id){
    xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    //set to read n doubles... decent guess
    size_t size = (*n)*(8);
    H5Pset_buffer(xfer_plist, size, NULL, NULL);
  } 
  
  err = H5Dread(d_id, t_id, ms_id, ds_id, xfer_plist, *data);
  if(err < 0){
    std::cerr << "HDF_IO: Could not read data " << directory+dataName << std::endl;
  }
  
  //clean up
  H5Pclose(xfer_plist);
  H5Sclose(ms_id);
  H5Tclose(ht_id);
  H5Dclose(d_id);
  H5Sclose(ds_id);
  H5Gclose(group);
 
  return err;
}

template <class Type>
hid_t HDF_WriteScalar(hid_t fileId, std::string directory, std::string dataName, Type* scalar)
{
  hid_t xfer_plist = H5P_DEFAULT;
  hid_t group = -1;
  hid_t dset = -1;
  hid_t space = -1;
  hid_t t_id = -1;
  
  Int err = 0;

  hsize_t dims[] = {1};

  //get type that was passed in
  Type check = 0;
  t_id = HDF_GetType(check);

  //open the group
  group = HDF_GetDirectory(fileId, directory);

  //create the dataspace
  space = H5Screate_simple(1, dims, NULL);
  if(space < 0){
    std::cerr << "HDF_IO: Could not create dataspace <" << directory + dataName 
	      << "> in HDF_WriteScalar()!" << std::endl;
  }

  //open the dataset
  dset = H5Dopen(group, dataName.c_str(), xfer_plist);
  if(dset < 0){
    dset = H5Dcreate(group, dataName.c_str(), t_id, space, xfer_plist, xfer_plist, xfer_plist);
    if(dset < 0){
      std::cerr << "HDF_IO: Cannot create dataset " << directory+dataName 
		<< " in HDF_WriteScalar()!" << std::endl;
      H5Gclose(group);
      return (-1);
    }
  }

  //write the scalar
  err = H5Dwrite(dset, t_id, H5S_ALL, space, xfer_plist, scalar);
  if(err < 0){
    std::cerr << "HDF_IO: Write failed to " << directory+dataName << " in HDF_WriteScalar()!" << std::endl;
  }

  //clean up
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(group);

  return err;
}


template <class Type>
hid_t HDF_WriteArray(hid_t fileId, std::string directory, std::string dataName, Type* data, Int nrows, Int ncols)
{
  hid_t space, group, dset, t_id;

  hsize_t rank = 2;
  hsize_t dims[rank];
  dims[0] = nrows;
  dims[1] = ncols;
  hsize_t sstride[2], start[2];

  //transfer plist to set up temp buffer
  hid_t xfer_plist = H5P_DEFAULT;

  herr_t err = 0;

  Type check = 0;
  //get type that was passed in
  t_id = HDF_GetType(check);
  if(t_id < 0){
    std::cerr << "HDF_IO: Cannot determine H5 datatype -- this is going to fail!" << std::endl;
  }
  
  //setting maximum size to NULL sets the maximum size to the current size
  space = H5Screate_simple(2, dims, NULL);
  if(space < 0){
    std::cerr << "HDF_IO: Could not create dataspace <" << directory + dataName 
	      << "> in HDF_WriteArray()!" << std::endl;
  }

  //select the entire region of the dataspace
  start[0] = 0;
  start[1] = 0;
  sstride[0] = 1;
  sstride[1] = 1;

  //we don't use blocks here, pass NULL
  err = H5Sselect_hyperslab(space, H5S_SELECT_SET, start, sstride, dims, NULL);
  if(err < 0){
    std::cerr << "HDF_IO: Select hyperslab failed for writing in HDF_WriteArray()!" << std::endl;
    return (-1);
  }

  //open the group
  group = HDF_GetDirectory(fileId, directory);

  //attempt to open dataset, if fail, create it
  dset = H5Dopen(group, dataName.c_str(), xfer_plist);
  if(dset < 0){
    dset = H5Dcreate(group, dataName.c_str(), t_id, space, xfer_plist, xfer_plist, xfer_plist);
    if(dset < 0){
      std::cerr << "HDF_IO: Cannot create dataset " << directory+dataName << " in HDF_WriteArray()!" << std::endl;
      H5Gclose(group);
      return (-1);
    }
  }

  err = H5Dwrite(dset, t_id, H5S_ALL, space, xfer_plist, data);
  if(err < 0){
    std::cerr << "HDF_IO: Write failed to " << directory+dataName << " in HDF_WriteArray()!" << std::endl;
  }

  //clean up
  H5Dclose(dset);
  H5Gclose(group);
  H5Sclose(space);

  return err;
}

template <class Type>
hid_t HDF_WriteArrayAttribute(hid_t fileId, std::string directory, std::string dataName, 
			      std::string attrName, std::vector<Type>& attrArray)
{
  hid_t err = 0;
  hid_t group, aspace, attribute, dset;
  hid_t a_type;
  hsize_t nrows = attrArray.size();
  hsize_t ncols = 1;
  hsize_t rank = 2;
  hsize_t dims[rank];
  dims[0] = nrows;
  dims[1] = ncols;
  hid_t xfer_plist = H5P_DEFAULT;

  //get type that was passed in
  Type check = 0;
  a_type = HDF_GetType(check);
  if(a_type < 0){
    std::cerr << "HDF_IO: Cannot determine H5 datatype -- this is going to fail!" << std::endl;
  }

  aspace = H5Screate_simple(rank, dims, NULL);
  if(aspace < 0){
    std::cerr << "HDF_IO: Could not create attribute space in HDF_WriteArrayAttribute()" << std::endl;
  }
  group = HDF_GetDirectory(fileId, directory);
  
  dset = H5Dopen(group, dataName.c_str(), H5P_DEFAULT);
  if(dset < 0){
    dset = H5Dcreate(group, dataName.c_str(), a_type, aspace, xfer_plist, xfer_plist, xfer_plist);
    if(dset < 0){
      std::cerr << "HDF_IO: Cannot create dataset " << directory+dataName 
		<< " in HDF_WriteArrayAttribute()!" << std::endl;
      H5Gclose(group);
      return (-1);
    }
  }

  attribute = H5Acreate(dset, attrName.c_str(), a_type, aspace, H5P_DEFAULT, H5P_DEFAULT);
  if(attribute < 0){
    //if fail, delete then recreate it, not elegant but it works
    err = H5Adelete(dset, attrName.c_str());
    if(err < 0){
      std::cerr << "HDF_IO: Could not delete attribute " 
		<< directory+dataName << " - " << attrName << std::endl;
    }
    attribute = H5Acreate(dset, attrName.c_str(), a_type, aspace, H5P_DEFAULT, H5P_DEFAULT);
    if(attribute < 0){
      std::cerr << "HDF_IO: Could not create vector attribute " 
		<< directory+dataName << " - " << attrName << std::endl;
    }
  }

  err = H5Awrite(attribute, a_type, attrArray.data());
  if(err < 0){
    std::cerr << "HDF_IO: Writing vector attribute " << directory+dataName << " - " 
	      << attrName << " failed " << std::endl;
  }
  
  //free the resources
  H5Gclose(group);
  H5Dclose(dset);
  H5Aclose(attribute);
  H5Sclose(aspace);

  return err;
}

template <class Type>
hid_t HDF_ReadArrayAttribute(hid_t fileId, std::string directory, std::string dataName, 
			     std::string attrName, std::vector<Type>& attrArray)
{
  hid_t err = 0;
  hid_t group, aspace, attribute, dset;
  hsize_t rank = 2;
  hsize_t dims[rank];
  hid_t xfer_plist = H5P_DEFAULT;

  //get type of data expected
  Type check = 0;
  hid_t atype = HDF_GetType(check);
  if(atype < 0){
    std::cerr << "HDF_IO: Cannot determine H5 datatype -- this is going to fail!" << std::endl;
  }

  //get the group
  group = HDF_GetDirectory(fileId, directory);

  //open the dataset
  dset = H5Dopen(group, dataName.c_str(), xfer_plist);
  if(dset < 0){
    H5Gclose(group);
    std::cerr << "HDF_IO: Cannot open dataset " << directory+dataName 
	      << " in HDF_ReadArrayAttribute()" << std::endl;
    return (-1);
  }

  //open attribute
  attribute = H5Aopen(dset, attrName.c_str(), xfer_plist);
  if(attribute < 0){
    H5Gclose(group);
    H5Dclose(dset);
    std::cerr << "HDF_IO: Cannot open attribute " << directory+dataName << " - "
	      << attrName << " in HDF_ReadArrayAttribute() " << std::endl;
    return (-1);
  }

  //get the data type
  hid_t atype2 = H5Aget_type(attribute);

  //get the dataspace
  aspace = H5Aget_space(attribute);

  hsize_t rank2 = H5Sget_simple_extent_ndims(aspace);
  if(rank2 != rank){
    H5Gclose(group);
    H5Dclose(dset);
    H5Aclose(attribute);
    H5Sclose(aspace);
    std::cerr << "HDF_IO: Ranks non-matching for attribute " << directory+dataName << " - "
	      << attrName << " in HDF_ReadArrayAttribute() " << std::endl;
    return (-1);
  }
  err = H5Sget_simple_extent_dims(aspace, dims, NULL);
  Type* data = new Type[dims[0]*dims[1]];

  err = H5Aread(attribute, atype2, data);
  if(err < 0){
    std::cerr << "HDF_IO: Read failed for attribute " << directory+dataName << " - "
	      << attrName << " in HDF_ReadArrayAttribute() " << std::endl;
  }

  attrArray.reserve(dims[0]*dims[1]);
  for(hsize_t i = 0; i < dims[0]*dims[1]; i++){
    attrArray.push_back(data[i]);
  }

  //clean up
  H5Gclose(group);
  H5Dclose(dset);
  H5Sclose(aspace);
  H5Aclose(attribute);
  H5Tclose(atype2);
  delete [] data;

  return err;
}
