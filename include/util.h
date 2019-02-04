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

inline double hann_window(size_t in_n, size_t in_len)
{
#ifdef _DEBUG
    if (in_n < in_len) //in_n from 0 to len-1
#endif
        return 0.5 * (1. - cos(2 * pi * in_n / (double)(in_len - 1)));
#ifdef _DEBUG
    else
    {
        Fatal_Error("Hann window out of bound")
    }
#endif
}