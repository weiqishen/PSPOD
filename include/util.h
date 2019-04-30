/**
 * @file util.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2019-02-02
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once
#include <complex>
#include <cmath>
#define  MKL_Complex16 std::complex<double>
#define  MKL_Complex8 std::complex<float>
#include "mkl.h"
#include "ndarray.h"

extern double pi;
ndarray<double> transpose(ndarray<double> &in_array);
ndarray<MKL_Complex16> transpose(ndarray<MKL_Complex16> &in_array);
ndarray<MKL_Complex16> conj_transpose(ndarray<MKL_Complex16> &in_array);
ndarray<double> multi_2darray(ndarray<double> &M1, ndarray<double> &M2);

void hann_window(ndarray<double> &out_hann);
