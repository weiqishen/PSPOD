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
pod_base::pod_base()
{
}
pod_base::pod_base(size_t in_n_probe, size_t in_n_snap)
{
    //copy data for local
    n_realization_local = in_n_snap;
    n_probe_local = in_n_probe;
    n_fields = run_input.fields.get_len();
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
    //rescale modal energy with weights(assume uniform weight)
    mkl_dimatcopy ('C', 'N', D.get_len(), 1, run_input.dw, D.get_ptr(), D.get_len(), D.get_len());

}

void pod_base::write_results()
{
    //create hdf5 file and dataset
    hid_t of_id, attr_id, dataspace_id, dataset_id, datatype;
    hsize_t *dim = new hsize_t[2];
    int n_modes=U.get_dim(1);

    of_id = H5Fcreate(run_input.output_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //write metadata to root
    dataspace_id = H5Screate(H5S_SCALAR);
    //number of modes
    attr_id = H5Acreate(of_id, "n_modes", H5T_NATIVE_INT32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT32, &n_modes);
    H5Aclose(attr_id);
    //sampling interval
    attr_id = H5Acreate(of_id, "dt", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &run_input.dt);
    H5Aclose(attr_id);
    //total number of probes
    attr_id = H5Acreate(of_id, "n_probe", H5T_NATIVE_INT32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT32, &run_input.total_n_probe);
    H5Aclose(attr_id);
    H5Sclose(dataspace_id);

    dim[0]=run_input.xyz_0.get_len();
    dataspace_id = H5Screate_simple(1, dim, NULL);
    //coordinate of the first probe
    attr_id = H5Acreate(of_id, "xyz_0", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, run_input.xyz_0.get_ptr());
    H5Aclose(attr_id);
    //number of probes in each direction
    attr_id = H5Acreate(of_id, "np_xyz", H5T_NATIVE_INT32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT32, run_input.np_xyz.get_ptr());
    H5Aclose(attr_id);
    //seperation of probes in each direction
    attr_id = H5Acreate(of_id, "d_xyz", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, run_input.d_xyz.get_ptr());
    H5Aclose(attr_id);
    H5Sclose(dataspace_id);

    //fields
    dim[0]=run_input.fields.get_len();
    const char **temp_field=new const char*[run_input.fields.get_len()];
    for (size_t i=0;i<run_input.fields.get_len();i++)
        temp_field[i]=run_input.fields(i).c_str();
    dataspace_id=H5Screate_simple(1,dim,NULL);
    datatype = H5Tcopy (H5T_C_S1);
    H5Tset_size (datatype, H5T_VARIABLE);
    attr_id = H5Acreate(of_id, "fields", datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, datatype, temp_field);
    H5Aclose(attr_id);
    H5Tclose(datatype);
    H5Sclose(dataspace_id);
    //write modal energy
     dim[0]= n_modes;
    dataspace_id = H5Screate_simple(1, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "/modal_energy", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D.get_ptr());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    //write modes
    dim[1] = U.get_dim(0);
    dim[0] = n_modes;
    dataspace_id = H5Screate_simple(2, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "/modes", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, U.get_ptr());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    //write mean if needed
    if (run_input.write_mean)
    {
        dim[0] = mean_data.get_dim(0);
        dataspace_id = H5Screate_simple(1, dim, NULL);
        dataset_id = H5Dcreate2(of_id, "/mean", H5T_NATIVE_DOUBLE, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mean_data.get_ptr());
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
    }
    //close file
    H5Fclose(of_id);
    delete[] dim;
    delete[] temp_field;
}

double *pod_base::get_data_ptr()
{
    return real_data.get_ptr();
}
