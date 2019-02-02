#include "util.h"

ndarray<double> transpose(ndarray<double> &in_array)
{
  // Get dimensions of arrays
  if (in_array.get_n_dim() == 2)
  {
    size_t dim_0 = in_array.get_dim(0);
    size_t dim_1 = in_array.get_dim(1);

    ndarray<double> transpose{dim_1, dim_0};
    mkl_domatcopy('C', 'T', dim_1, dim_0,
                  1., in_array.get_ptr(), dim_0, transpose.get_ptr(), dim_1);
    return transpose;
  }
  else
  {
    Fatal_Error("ERROR: Array transpose function only accepts a 2-dimensional square ndarray");
  }
}

ndarray<MKL_Complex16> transpose(ndarray<MKL_Complex16> &in_array)
{
  // Get dimensions of arrays
  if (in_array.get_n_dim() == 2)
  {
    size_t dim_0 = in_array.get_dim(0);
    size_t dim_1 = in_array.get_dim(1);

    ndarray<MKL_Complex16> transpose{dim_1, dim_0};
    mkl_zomatcopy('C', 'T', dim_1, dim_0,
                  1., in_array.get_ptr(), dim_0, transpose.get_ptr(), dim_1);
    return transpose;
  }
  else
  {
    Fatal_Error("ERROR: Array transpose function only accepts a 2-dimensional square ndarray");
  }
}

ndarray<MKL_Complex16> conj_transpose(ndarray<MKL_Complex16> &in_array)
{
  // Get dimensions of arrays
  if (in_array.get_n_dim() == 2)
  {
    size_t dim_0 = in_array.get_dim(0);
    size_t dim_1 = in_array.get_dim(1);

    ndarray<MKL_Complex16> transpose{dim_1, dim_0};
    mkl_zomatcopy('C', 'C', dim_1, dim_0,
                  MKL_Complex16(1.,0.), in_array.get_ptr(), dim_0, transpose.get_ptr(), dim_1);
    return transpose;
  }
  else
  {
    Fatal_Error("ERROR: Array transpose function only accepts a 2-dimensional square ndarray");
  }
}

// Multiply M1(L*M) by M2(M*N)
ndarray<double> multi_2darray(ndarray<double> &M1, ndarray<double> &M2)
{
  // Get dimensions of arrays
  size_t dim_1_0 = M1.get_dim(0);
  size_t dim_1_1 = M1.get_dim(1);

  size_t dim_2_0 = M2.get_dim(0);
  size_t dim_2_1 = M2.get_dim(1);

  // Only 2D arrays
  if (M1.get_n_dim() == 2 && M2.get_n_dim() == 2)
  {
    // Ensure consistent inner dimensions
    if (dim_1_1 == dim_2_0)
    {
      ndarray<double> product{dim_1_0, dim_2_1};
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, dim_1_0, dim_2_1, dim_1_1,
                  1.0, M1.get_ptr(), dim_1_0, M2.get_ptr(), dim_2_0, 0.0, product.get_ptr(), dim_1_0);
      return product;
    }
    else
    {
      Fatal_Error("ERROR: ndarray dimensions are not compatible in multiplication function");
    }
  }
  else
  {
    Fatal_Error("ERROR: Array multiplication function can only multiply 2-dimensional arrays together");
  }
}
