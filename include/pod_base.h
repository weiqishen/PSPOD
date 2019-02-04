
/**
 * @file pod_base.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2019-02-02
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once
#include "global.h"
#include "util.h"//contain array, mkl complex definition and tool functions

class pod_base
{
  public:
    //constructors/destructors
    pod_base();
    pod_base(size_t in_n_probe, size_t in_n_snap);
    ~pod_base();

    //computational methods
    virtual void calc_mode(); //!<calculate mode from real_data using eig(real_data'*real_data)
    void calc_mean();         //!<calculate mean of real_data along axis=0
    void subtract_mean();     //!<subtract mean from real_data along axis=0

    //output methods
    virtual void write_results(); //!<write modal energy, modes, etc. to hdf5 file
    //others
    double *get_data_ptr(); //!<get pointer to data_array

  protected:
    //meta data to describe local data array
    size_t n_probe_local;
    size_t n_realization_local;
    size_t n_fields;

    //heavy data for pod calculation
    ndarray<double> real_data; //time*space
    ndarray<double> mean_data; //space
    ndarray<double> U, D;      //!<POD modes(space*n), modal energy(n)
    ndarray<double> a;         //!<temporal coefficient(time*n)
};
