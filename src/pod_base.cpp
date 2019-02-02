
#include "pod_base.h"
#include <numeric>

using namespace std;
pod_base::pod_base()
{
}
pod_base::pod_base(size_t in_n_probe, size_t in_n_snap, size_t in_n_fields)
{
    //copy data for local
    n_realization_local = in_n_snap;
    n_probe_local = in_n_probe;
    n_fields = in_n_fields;
    //declare local array
    real_data.setup({n_realization_local, n_probe_local * n_fields});
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
    U.setup({real_data.get_dim(1), real_data.get_dim(1)});
    U = 0.;                        //initialize to 0.
    D.setup(real_data.get_dim(1)); //space
    //do matrix multiplication
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, real_data.get_dim(1), real_data.get_dim(1), real_data.get_dim(0), 1.0,
                real_data.get_ptr(), real_data.get_dim(0), real_data.get_ptr(), real_data.get_dim(0), 0.,
                U.get_ptr(), real_data.get_dim(1));

    //delete upper triangular elements
    for (size_t i = 0; i < U.get_dim(1); i++) //loop over row
        for (size_t j = 0; j < i; j++)        //loop over 0->i-1
            U({j, i}) = 0.;

    //compute eigenvalue decomposition
    LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'L', U.get_dim(0), U.get_ptr(), U.get_dim(0), D.get_ptr());
}

void pod_base::write_results()
{
    //create hdf5 file

    //write metadata

    //write modal energy

    //write modes
}

double *pod_base::get_data_ptr()
{
    return real_data.get_ptr();
}
