#include "global.h"
#include <complex>
#define  MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "ndarray.h"

ndarray<double> transpose(ndarray<double> &in_array);
ndarray<MKL_Complex16> transpose(ndarray<MKL_Complex16> &in_array);
ndarray<MKL_Complex16> conj_transpose(ndarray<MKL_Complex16> &in_array);
ndarray<double> multi_2darray(ndarray<double> &M1, ndarray<double> &M2);