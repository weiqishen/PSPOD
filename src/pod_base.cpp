
#include "pod_base.h"
#include "mkl.h"
#include <numeric>

using namespace std;
pod_base::pod_base()
{
}
pod_base::pod_base(size_t n_probe, size_t n_snap, size_t n_fields)
{
    real_data.setup({n_snap, n_probe * n_fields});
}

pod_base::~pod_base()
{
}
void pod_base::calc_mean()
{
    mean_data.setup(real_data.get_dim(1));
    for (size_t i = 0; i < real_data.get_dim(1); i++)
    {
        mean_data(i) = accumulate(real_data.get_ptr({0, i}), real_data.get_ptr({0, i + 1}), 0.) / real_data.get_dim(0);
    }
}

void pod_base::subtract_mean()
{
    for (size_t i = 0; i < real_data.get_dim(1); i++)
        for (size_t j = 0; j < real_data.get_dim(0); j++)
            real_data({j, i}) -= mean_data(i);
}

void pod_base::calc_mode()
{
    //setup eigenvector and eigenvalue arrays
    U.setup({real_data.get_dim(1), real_data.get_dim(1)});//space*space
    D.setup(real_data.get_dim(1));//space
    //do matrix multiplication
    U = 0.;
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, real_data.get_dim(1), real_data.get_dim(1), real_data.get_dim(0), 1.0,
                real_data.get_ptr(), real_data.get_dim(0), real_data.get_ptr(), real_data.get_dim(0), 0.,
                U.get_ptr(), real_data.get_dim(1));
                
    //delete upper triangular elements
    for (size_t i = 0; i < U.get_dim(1); i++) //row for dim0
        for (size_t j = 0; j < i; j++)        //column for dim0-1
            U({j, i}) = 0.;

    //compute eigenvalue decomposition
    LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'L', real_data.get_dim(1), U.get_ptr(), real_data.get_dim(1), D.get_ptr());
}

double *pod_base::get_data_ptr()
{
    return real_data.get_ptr();
}
