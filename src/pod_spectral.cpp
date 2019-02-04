/**
 * @file pod_spectral.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2019-02-02
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include <numeric>
#include "pod_specteal.h"
#include "hdf5.h"

using namespace std;

pod_spectral::pod_spectral(size_t in_n_probe, size_t in_block_size, size_t in_n_blocks)
{
    n_realization_local = in_n_blocks;
    n_probe_local = in_n_probe;
    n_fields = run_input.fields.get_len();
    block_size_local = in_block_size;
    real_data.setup({block_size_local, n_probe_local * n_fields});
    fft_data.setup({block_size_local / 2 + 1, n_probe_local * n_fields, n_realization_local});
}

pod_spectral::~pod_spectral()
{
}

void pod_spectral::calc_fft(size_t block_id)
{
    //create and set descriptor
    DFTI_DESCRIPTOR_HANDLE desc_1;
    DftiCreateDescriptor(&desc_1, DFTI_DOUBLE, DFTI_REAL, 1, block_size_local);
    DftiSetValue(desc_1, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(desc_1, DFTI_NUMBER_OF_TRANSFORMS, real_data.get_dim(1));
    DftiSetValue(desc_1, DFTI_FORWARD_SCALE, 1. / block_size_local);
    DftiSetValue(desc_1, DFTI_INPUT_DISTANCE, block_size_local);
    DftiSetValue(desc_1, DFTI_OUTPUT_DISTANCE, block_size_local / 2 + 1);
    DftiSetValue(desc_1, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiCommitDescriptor(desc_1);

    //compute fft
    MKL_LONG status;
    status = DftiComputeForward(desc_1, real_data.get_ptr(), fft_data.get_ptr({0, 0, block_id}));
    //check state
    if (status && !DftiErrorClass(status, DFTI_NO_ERROR))
    {
        printf("Error: %s\n", DftiErrorMessage(status));
        exit(1);
    }
    DftiFreeDescriptor(&desc_1);
}

void pod_spectral::calc_mode()
{
    //scale the amplitude of half spectrum
    for (size_t i = 0; i < fft_data.get_len() / fft_data.get_dim(0); i++)//loop over space*block
    {
        transform(fft_data.get_ptr({1, i}), fft_data.get_ptr({block_size_local - block_size_local / 2, i}),
                  fft_data.get_ptr({1, i}), [](MKL_Complex16 x) { return x * 2.; });
    }
    //reshape fft data ->space*block*frequency
    fft_data.reshape({block_size_local / 2 + 1, n_probe_local * n_fields * n_realization_local});
    fft_data.trans();
    fft_data.reshape({n_probe_local * n_fields, n_realization_local, block_size_local / 2 + 1});
    //calculate svd
    U_spectral.setup({fft_data.get_dim(0), min(fft_data.get_dim(0), fft_data.get_dim(1)), fft_data.get_dim(2)});
    D.setup({min(fft_data.get_dim(0), fft_data.get_dim(1)), fft_data.get_dim(2)});
    //declare dummy array
    ndarray<MKL_Complex16> vt_dumm(1);
    ndarray<double> superb{min(fft_data.get_dim(0), fft_data.get_dim(1)) - 1};//min(m,n)-1
    //compute svd for each frequency
    for (size_t i = 0; i < fft_data.get_dim(2); i++)
    {
        LAPACKE_zgesvd(LAPACK_COL_MAJOR, 'S', 'N', fft_data.get_dim(0), fft_data.get_dim(1),
                       fft_data.get_ptr({0, 0, i}), fft_data.get_dim(0), D.get_ptr({0, i}), U_spectral.get_ptr({0, 0, i}),
                       U_spectral.get_dim(0), vt_dumm.get_ptr(), fft_data.get_dim(1), superb.get_ptr());
    }
    //scale singular value to modal energy
    vdSqr(D.get_len(), D.get_ptr(), D.get_ptr());
    mkl_dimatcopy('C', 'N', D.get_dim(0), D.get_dim(1), run_input.dw, D.get_ptr(), D.get_dim(0), D.get_dim(0));
}

void pod_spectral::write_results()
{
    //create hdf5 file and dataset
    hid_t of_id, attr_id, dataspace_id, dataset_id, datatype;
    hsize_t *dim = new hsize_t[3];

    of_id = H5Fcreate(run_input.output_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //write metadata to root
    dataspace_id = H5Screate(H5S_SCALAR);
    //number of modes
    int n_modes_pf = U_spectral.get_dim(1);
    attr_id = H5Acreate(of_id, "n_modes", H5T_NATIVE_INT32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT32, &n_modes_pf);
    H5Aclose(attr_id);
    //sampling interval
    double df = 1 / (run_input.dt * run_input.block_size);
    attr_id = H5Acreate(of_id, "df", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &df);
    H5Aclose(attr_id);
    //number of frequency components
    int n_freq = U_spectral.get_dim(2);
    attr_id = H5Acreate(of_id, "n_freq", H5T_NATIVE_INT32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT32, &n_freq);
    H5Aclose(attr_id);
    //total number of probes
    attr_id = H5Acreate(of_id, "n_probe", H5T_NATIVE_INT32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT32, &run_input.total_n_probe);
    H5Aclose(attr_id);
    H5Sclose(dataspace_id);

    dim[0] = run_input.xyz_0.get_len();
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
    dim[0] = run_input.fields.get_len();
    const char **temp_field = new const char *[run_input.fields.get_len()];
    for (size_t i = 0; i < run_input.fields.get_len(); i++)
        temp_field[i] = run_input.fields(i).c_str();
    dataspace_id = H5Screate_simple(1, dim, NULL);
    datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, H5T_VARIABLE);
    attr_id = H5Acreate(of_id, "fields", datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, datatype, temp_field);
    H5Aclose(attr_id);
    H5Tclose(datatype);
    H5Sclose(dataspace_id);
    //write modal energy
    dim[1] = n_modes_pf;
    dim[0] = n_freq;
    dataspace_id = H5Screate_simple(2, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "/modal_energy", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D.get_ptr());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    //write modes
    dim[2] = U_spectral.get_dim(0);
    dim[1] = n_modes_pf;
    dim[0] = n_freq;
    dataspace_id = H5Screate_simple(3, dim, NULL);
    U.setup(U_spectral.get_len());
    for (size_t i = 0; i < U_spectral.get_len(); i++)
        U(i) = U_spectral(i).real();
    dataset_id = H5Dcreate2(of_id, "/modes_real", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, U_spectral.get_ptr());
    H5Dclose(dataset_id);

    for (size_t i = 0; i < U_spectral.get_len(); i++)
        U(i) = U_spectral(i).imag();
    dataset_id = H5Dcreate2(of_id, "/modes_imag", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, U_spectral.get_ptr());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    //close file
    H5Fclose(of_id);
    delete[] dim;
    delete[] temp_field;
}
