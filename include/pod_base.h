#pragma once
#include "global.h"
#include "ndarray.h"

class pod_base
{
    public:
    pod_base();
    pod_base(size_t n_probe,size_t n_snap,size_t n_fields);
    ~pod_base();

    virtual void calc_mode();
    void calc_mean();
    void subtract_mean();

    double* get_data_ptr();
    protected:
    ndarray<double> real_data;
    ndarray<double> mean_data;
    ndarray<double> U,D;//right eigenvectors/left singularvectors and eigenvalues
};