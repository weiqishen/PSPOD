#pragma once
#include "global.h"
#include "ndarray.h"
#include "pod_base.h"
class pod_spectral:public pod_base
{
    public:
    pod_spectral();
    pod_spectral(size_t n_probe,size_t n_snap,size_t n_fields);
    ~pod_spectral();
    void calc_fft();
    void calc_mode();

};