#include "parallel.h"

Int MPI_TypeMap(int a) {return TYPE_INT;};
Int MPI_TypeMap(float a) {return TYPE_FLOAT;};
Int MPI_TypeMap(double a) {return TYPE_DOUBLE;};
Int MPI_TypeMap(std::complex<double> a) {return TYPE_COMPLEX_DOUBLE;};
Int MPI_TypeMap(std::complex<float> a) {return TYPE_COMPLEX_FLOAT;};
