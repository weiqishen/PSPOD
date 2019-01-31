
#include "pod_snap.h"
#include "mkl.h"
#include <numeric>
#include "util.h"
using namespace std;

pod_snap::pod_snap(size_t n_probe, size_t n_snap, size_t n_fields) : pod_base(n_probe, n_snap, n_fields)
{
}

pod_snap::~pod_snap()
{
}

void pod_snap::calc_mode()
{
    real_data.trans();
    U.setup({real_data.get_dim(0), real_data.get_dim(1)});
    D.setup(real_data.get_dim(1));
    //declare dummy array
    ndarray<double>vt_dumm(1);
    ndarray<double>superb{min(real_data.get_dim(0), real_data.get_dim(1))-1};

    //compute svd
    LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'N', real_data.get_dim(0),  real_data.get_dim(1),
                   real_data.get_ptr(), real_data.get_dim(0), D.get_ptr(), U.get_ptr(),
                   U.get_dim(0), vt_dumm.get_ptr(), real_data.get_dim(1), superb.get_ptr());
}
