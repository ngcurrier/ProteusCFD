#include <iostream>

class BitField
{
  
public:
  BitField() :
  bits(NULL), size(0), bitsize(0), allocated(false)
  { };
  
  BitField(int length) :
    bitsize(length), size(length / 8)
  {
    if(bitsize % 8){
       size++;
    }
    bits = new unsigned char[size];
    allocated = true;
    clear();
  };

  ~BitField()
  {
    if (allocated){
      delete [] bits;
    }
    size = 0;
    bitsize = 0;
  };

  void print()
  {
    std::cout << "BITS: ";
    for(int i = 0; i < bitsize; i++){
      std::cout << get(i) << " ";
    }
    std::cout << std::endl;
  };
    

  void clear()
  {
    for(int i = 0; i < size; i++){
      bits[i] = 0;
    }
  };

  void set(int loc, bool val){
    unsigned char* byte = bits + loc/8;
    unsigned char tgt = (!!val) << (loc % 8);
    unsigned char mask = ~(0x01 << (loc % 8));
    *byte &= mask;
    *byte |= tgt;
  };

  bool get(int loc){
    unsigned char byte = bits[loc/8];
    unsigned char mask = 0x01 << (loc % 8);
    return (bool)(!!(byte & mask));
  };

private:
  unsigned char * bits;
  int size;
  int bitsize;
  bool allocated;
};
