
#include "endian_util.h"

typedef union {
  double value;
  char   bytes[8];
} double_endian;


typedef union {
  float value;
  char   bytes[4];
} float_endian;


int IsLittleEndian(void)
{
  short w = 0x0001;
  char * byte = (char*)&w;
  return (*byte ? 1 : 0 );
}


double ReverseEndian_Double8(double in)
{
  int i;
  double_endian out;
  double_endian temp;

  temp.value = in;

  //
  // the fixing of this was truly an AAAARGH! moment.  On the Intel 9.0 compiler,
  // the original version of the code (hardcoded, no loop) just did not work at 
  // all; endian reversals failed to happen.
  //
  // This turned out to be a compiler optimization issue.  Normally we compile
  // with -O3 -xW -Ob2 flags, and these interacted (remember 9.0 ONLY) in some
  // very strange ways.  If I downgraded any one of the three (-O3 replace with -O2,
  // removing -xW altogether, or -Ob1 instead of -Ob2), then it did not matter what
  // I did with the other two, and it would work with the Intel 9.0 compiler.
  //
  // The suspicion is that if the compiler inlines this function, it did not 
  // properly take into account the dependencies (ref. ReverseEndian_Double8Array).
  //
  // The other thing that works is replacing the hardcoding with a loop, as is done
  // below.  If we do this, we can use any of the optimization flags we want.
  // Bizarre compiler bug.  As stated before, this strange behavior only happened 
  // with Intel 9.0.  So, we will call it fixed, and leave this as a momento of
  // the amazingness of compiler bugs.
  // 
  // Q.E.D.
  //
#if 1
  for (i = 0; i < 8; i++) out.bytes[i] = temp.bytes[7-i];
#else
  out.bytes[0] = temp.bytes[7];
  out.bytes[1] = temp.bytes[6];
  out.bytes[2] = temp.bytes[5];
  out.bytes[3] = temp.bytes[4];
  out.bytes[4] = temp.bytes[3];
  out.bytes[5] = temp.bytes[2];
  out.bytes[6] = temp.bytes[1];
  out.bytes[7] = temp.bytes[0];
#endif

  return(out.value);
}


void ReverseEndian_Double8Array(double *in,int n)
{
    int i;

  for (i = 0; i < n; i++)
    {
      in[i] = ReverseEndian_Double8(in[i]);
    }

  return;
}


void ReverseEndian_Int4Array(int * in,int n)
{
  int i;

  for (i = 0; i < n; i++)
    {
      in[i] = ReverseEndian_Int4(in[i]);
    }

  return;
}

int ReverseEndian_Int4(int in)
{
  int retval;

  retval  = ((in & 0x000000FF) << 24);
  retval |= ((in & 0x0000FF00) <<  8);
  retval |= ((in & 0x00FF0000) >>  8);
  retval |= ((in & 0xFF000000) >> 24);

  return(retval);
}


short int ReverseEndian_Int2(short int in)
{
  short int retval;

  retval  = ((in & 0x00FF) << 8);
  retval |= ((in & 0xFF00) >> 8);

  return(retval);
}

void ReverseEndian_Int2Array(short int * in,int n)
{
  int i;

  for (i = 0; i < n; i++)
    {
      in[i] = ReverseEndian_Int2(in[i]);
    }

  return;
}

void ReverseEndian_Float4Array(float * in,int n)
{
  int i;

  for (i = 0; i < n; i++)
    {
      in[i] = ReverseEndian_Float4(in[i]);
    }

  return;
}

float ReverseEndian_Float4(float in)
{ 
  float_endian out;
  float_endian temp;

  temp.value = in;
  
  out.bytes[0] = temp.bytes[3];
  out.bytes[1] = temp.bytes[2];
  out.bytes[2] = temp.bytes[1];
  out.bytes[3] = temp.bytes[0];

  return(out.value);
}

