/**
 * @file pod_base.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2019-02-02
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include <numeric>
#include "pod_base.h"
#include "hdf5.h"

using namespace std;

pod_base::pod_base(size_t in_n_probe, size_t in_n_snap, double in_dt)
{
    //copy data for local
    n_realization = in_n_snap;
    n_probe = in_n_probe;
    dt = in_dt;
    //declare data array space*time
    real_data.setup({n_probe * run_input.fields_pod.get_len(), n_realization});
}

void pod_base::calc_mean()
{
    mean_data.setup(real_data.get_dim(0));
    mean_data = 0;
    for (size_t j = 0; j < real_data.get_dim(1); j++)
        for (size_t i = 0; i < real_data.get_dim(0); i++)
            mean_data(i) += real_data({i, j});
    cblas_dscal(mean_data.get_len(), 1. / real_data.get_dim(1), mean_data.get_ptr(), 1);
}

void pod_base::subtract_mean()
{
    for (size_t j = 0; j < real_data.get_dim(1); j++)
        for (size_t i = 0; i < real_data.get_dim(0); i++)
            real_data({i, j}) -= mean_data(i);
}

void pod_base::calculateWeight(double *in_coord)
{
    double dw = 1.;

    //calculate product of increment
    for (size_t i = 0; i < run_input.d_xyz.get_len(); i++)
        dw *= run_input.d_xyz(i);

    if (run_input.coord_sys == CARTESIAN) //only one number is needed
    {
        w.setup(1);
        w = dw;
    }
    else if (run_input.coord_sys == CYLINDRICAL) //one weight for each probe
    {
        w.setup(n_probe);
        for (size_t i = 0; i < n_probe; i++)
        {
            w(i) = dw * in_coord[run_input.d_xyz.get_len() * i]; //rdrdthetadz
        }
    }
    else
    {
        Fatal_Error("Unsupported coordinate system!")
    }
}

void pod_base::calc_mode()
{
    U.setup({real_data.get_dim(0), min(real_data.get_dim(0), real_data.get_dim(1))});
    D.setup(min(real_data.get_dim(0), real_data.get_dim(1)));
    a.setup({real_data.get_dim(1), U.get_dim(1)});

    //multiply sqrt(w) matrix
    if (run_input.coord_sys == CARTESIAN) //scale uniformly
        cblas_dscal(real_data.get_len(), sqrt(w(0)), real_data.get_ptr(), 1);
    else if (run_input.coord_sys == CYLINDRICAL)
    {
        for (size_t j = 0; j < run_input.fields_pod.get_len(); j++)                                                           //loop over each field
            for (size_t i = 0; i < n_probe; i++)                                                                              //loop over each probe
                cblas_dscal(real_data.get_dim(1), sqrt(w(i)), real_data.get_ptr({i + j * n_probe, 0}), real_data.get_dim(0)); //rescale each row
    }
    //scale by sqrt(1/n_realization)
    cblas_dscal(real_data.get_len(), sqrt(1. / n_realization), real_data.get_ptr(), 1);

    //declare dummy array
    ndarray<double> vt_dumm(1);
    ndarray<double> superb{min(real_data.get_dim(0), real_data.get_dim(1)) - 1}; //min(m,n)-1
    //compute svd
    ndarray<double> temp_real = real_data; //copy real data
    LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'N', real_data.get_dim(0), real_data.get_dim(1),
                   real_data.get_ptr(), real_data.get_dim(0), D.get_ptr(), U.get_ptr(),
                   U.get_dim(0), vt_dumm.get_ptr(), U.get_dim(1), superb.get_ptr());
    //rescale singlular value to be correct eigenvalue
    vdSqr(D.get_len(), D.get_ptr(), D.get_ptr());
    //calc coefficient X^T*U
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, temp_real.get_dim(1), U.get_dim(1), temp_real.get_dim(0), sqrt(n_realization), temp_real.get_ptr(), temp_real.get_dim(0), U.get_ptr(), U.get_dim(0), 0, a.get_ptr(), a.get_dim(0));
    //multiply 1/sqrt(w) matrix
    if (run_input.coord_sys == CARTESIAN) //scale uniformly
        cblas_dscal(U.get_len(), 1. / sqrt(w(0)), U.get_ptr(), 1);
    else if (run_input.coord_sys == CYLINDRICAL)
    {
        for (size_t j = 0; j < run_input.fields_pod.get_len(); j++)                                        //loop over each field
            for (size_t i = 0; i < n_probe; i++)                                                           //loop over each probe
                cblas_dscal(U.get_dim(1), 1. / sqrt(w(i)), U.get_ptr({i + j * n_probe, 0}), U.get_dim(0)); //rescale each row
    }
}

void pod_base::write_results()
{
    //create hdf5 file and dataset
    hid_t of_id, attr_id, dataspace_id, dataset_id, datatype;
    hsize_t dim[2];
    int n_modes = U.get_dim(1);

    of_id = H5Fcreate(run_input.output_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //write metadata to root
    dataspace_id = H5Screate(H5S_SCALAR);
    attr_id = H5Acreate(of_id, "dt", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &dt);
    H5Aclose(attr_id);
    H5Sclose(dataspace_id);

    //fields
    dim[0] = run_input.fields_pod.get_len();
    const char **temp_field = new const char *[dim[0]];
    for (size_t i = 0; i < dim[0]; i++)
        temp_field[i] = run_input.fields_pod(i).c_str();
    dataspace_id = H5Screate_simple(1, dim, NULL);
    datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, H5T_VARIABLE);
    attr_id = H5Acreate(of_id, "fields", datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, datatype, temp_field);
    delete[] temp_field;
    H5Aclose(attr_id);
    H5Tclose(datatype);
    H5Sclose(dataspace_id);

    //write modal energy
    dim[0] = n_modes;
    dataspace_id = H5Screate_simple(1, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "modal_energy", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D.get_ptr());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

    //write modes
    dim[1] = U.get_dim(0);
    dim[0] = n_modes;
    dataspace_id = H5Screate_simple(2, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "modes", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, U.get_ptr());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

    //write coefficients
    dim[1] = a.get_dim(0);
    dim[0] = a.get_dim(1);
    dataspace_id = H5Screate_simple(2, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "coeff", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, a.get_ptr());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

    //write mean if needed
    if (run_input.write_mean)
    {
        dim[0] = mean_data.get_dim(0);
        dataspace_id = H5Screate_simple(1, dim, NULL);
        dataset_id = H5Dcreate2(of_id, "mean", H5T_NATIVE_DOUBLE, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mean_data.get_ptr());
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
    }
    //close file
    H5Fclose(of_id);
}

void pod_base::write_coord(ndarray<double> &in_coord)
{
    //create hdf5 file and dataset
    hid_t of_id, dataspace_id, dataset_id;
    hsize_t dim[2];
    of_id = H5Fopen(run_input.output_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (of_id < 0)
        Fatal_Error("Failed to open output file!");
    dim[1] = in_coord.get_dim(0);
    dim[0] = in_coord.get_dim(1);
    dataspace_id = H5Screate_simple(2, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "coord", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, in_coord.get_ptr());
    //close file
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(of_id);
}