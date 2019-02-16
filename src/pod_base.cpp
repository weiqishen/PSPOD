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
    //declare local array
    real_data.setup({n_realization_local, in_n_probe * run_input.fields.get_len()});
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
                U.get_ptr(), real_data.get_dim(1)); //S=Q*Q^T
    //multiply weight(probe*field)            
    for (size_t i = 0; i < U.get_dim(1); i++)
        mkl_dimatcopy('C', 'N', U.get_dim(0), 1, w(i), U.get_ptr({0, i}), U.get_dim(0), U.get_dim(0));

    //delete upper triangular elements
    for (size_t i = 0; i < U.get_dim(1); i++) //loop over row
        for (size_t j = 0; j < i; j++)        //loop over 0->i-1
            U({j, i}) = 0.;

    //compute eigenvalue decomposition
    LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'L', U.get_dim(0), U.get_ptr(), U.get_dim(0), D.get_ptr());

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
    H5Awrite(attr_id, H5T_NATIVE_INT32, &run_input.n_probe);
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

void pod_base::write_coord(ndarray<double> &in_coord)
{
    //create hdf5 file and dataset
    hid_t of_id, dataspace_id, dataset_id;
    hsize_t *dim = new hsize_t[2];
    of_id = H5Fopen(run_input.output_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (of_id < 0)
        Fatal_Error("Failed to open output file!");
    dim[1] = in_coord.get_dim(0);
    dim[0] = in_coord.get_dim(1);
    dataspace_id = H5Screate_simple(2, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "/coord", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, in_coord.get_ptr());
    //close file
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(of_id);
    delete[] dim;
}