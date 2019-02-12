#ifndef MEM_UTIL_H__
#define MEM_UTIL_H__

#include <cstring>
#include <iostream>

/*******************************************************/
//Function to either shrink or extend a heap memory block
//while maintaining contents if extending and truncating
//end entries if shrinking
/*******************************************************/
template <class theType>
void MemResize(theType** data, size_t old_size, size_t new_size){
  theType* temp = NULL;
  temp = new theType[new_size];
  if(temp == NULL){
    std::cerr << "WARNING: MemResize() just failed" << std::endl;
  }
  if(old_size < new_size){
    memcpy(temp, *data, old_size*sizeof(theType)); 
  }
  else{
    memcpy(temp, *data, new_size*sizeof(theType));
  }
  //free old mem block
  delete [] *data;
  //reassign pointer to new mem block
  *data = temp;
}

//this has similar functionality to the above but operates on a list
//of pointers to an object of type theType*
template <class theType>
void MemResize(theType*** data, size_t old_size, size_t new_size){
  theType** temp = NULL;
  temp = new theType*[new_size];
  if(temp == NULL){
    std::cerr << "WARNING: MemResize() just failed" << std::endl;
  }
  if(old_size < new_size){
    memcpy(temp, *data, old_size*sizeof(theType*)); 
  }
  else{
    memcpy(temp, *data, new_size*sizeof(theType*));
  }
  //free old mem block
  delete [] *data;
  //reassign pointer to new mem block
  *data = temp;
}

//Blanks memory to zero
template <class Type>
void MemBlank(Type* x, int n){
  int i;
  for(i = 0; i < n; i++) x[i] = 0.0;
};

//sets memory to value
template <class Type>
void MemSet(Type* x, Type val, int n)
{
  int i;
  for(i = 0; i < n; i++) x[i] = val;
};

template <class Type, class Type2>
void DuplicateArray(Type** dataNew, const Type2* dataToCopy, size_t n)
{
  *dataNew = new Type[n];
  if(*dataNew == NULL){
    std::cerr << "WARNING: MemResize() just failed" << std::endl;
  }
  for(int i = 0; i < n; i++){
    (*dataNew)[i] = dataToCopy[i];
  }
}

#endif
