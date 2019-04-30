/**
 * @file pod_snap.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2019-02-02
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include <numeric>
#include "pod_snap.h"

using namespace std;

pod_snap::pod_snap(size_t in_n_probe, size_t in_n_snap) : pod_base(in_n_probe, in_n_snap)
{
}

pod_snap::~pod_snap()
{
}

void pod_snap::calc_mode()
{
    //square root of weight
    vdSqrt(w.get_len(), w.get_ptr(), w.get_ptr());
    //multiply weight(probe*field)
    //for catesian coord multiply const
    cblas_dscal(real_data.get_len(), w(0), real_data.get_ptr(), 1);
    U.setup({real_data.get_dim(0), min(real_data.get_dim(0),real_data.get_dim(1))});
    D.setup(min(real_data.get_dim(0),real_data.get_dim(1)));
    //declare dummy array
    ndarray<double>vt_dumm(1);
    ndarray<double>superb{min(real_data.get_dim(0), real_data.get_dim(1))-1};//min(m,n)-1
    //compute svd
    LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'N', real_data.get_dim(0),  real_data.get_dim(1),
                   real_data.get_ptr(), real_data.get_dim(0), D.get_ptr(), U.get_ptr(),
                   U.get_dim(0), vt_dumm.get_ptr(), U.get_dim(1), superb.get_ptr());
    //rescale singlular value to be correct eigenvalue
    vdSqr(D.get_len(),D.get_ptr(),D.get_ptr());
}
