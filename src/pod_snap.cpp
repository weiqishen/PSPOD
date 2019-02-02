
#include "pod_snap.h"
#include <numeric>
using namespace std;

pod_snap::pod_snap(size_t in_n_probe, size_t in_n_snap, size_t in_n_fields) : pod_base(in_n_probe, in_n_snap, in_n_fields)
{
}

pod_snap::~pod_snap()
{
}

void pod_snap::calc_mode()
{
    //transpose real_data array
    real_data.trans();//(time*space)->(sapce*time)
    U.setup({real_data.get_dim(0), real_data.get_dim(1)});
    D.setup(real_data.get_dim(1));
    //declare dummy array
    ndarray<double>vt_dumm(1);
    ndarray<double>superb{min(real_data.get_dim(0), real_data.get_dim(1))-1};//min(m,n)-1
    //compute svd
    LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'N', real_data.get_dim(0),  real_data.get_dim(1),
                   real_data.get_ptr(), real_data.get_dim(0), D.get_ptr(), U.get_ptr(),
                   U.get_dim(0), vt_dumm.get_ptr(), real_data.get_dim(1), superb.get_ptr());
    //rescale singlular value to be correct eigenvalue
    vdSqr(D.get_len(),D.get_ptr(),D.get_ptr());
}
