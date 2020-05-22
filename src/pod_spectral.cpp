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
#include <sys/stat.h>
using namespace std;

pod_spectral::pod_spectral(size_t in_n_probe, size_t in_block_size, size_t in_n_blocks, double in_dt)
{
    n_realization = in_n_blocks;
    block_size = in_block_size;
    n_probe = in_n_probe;
    dt = in_dt;
    real_data.setup({n_probe * run_input.fields_pod.get_len(), block_size});        //space*time
    fft_data.setup({block_size / 2 + 1, in_n_probe * run_input.fields_pod.get_len()}); //freq*space
    fft_comp.setup({fft_data.get_dim(0), fft_data.get_dim(1)});                        //same size as fft_data to write fft_data to file
    //create window array
    if (run_input.window)
    {
        hann_array.setup(block_size);
        hann_window(hann_array);
        hann_sqr = pow(cblas_dnrm2(block_size, hann_array.get_ptr(), 1), 2);
    }
}

void pod_spectral::calc_fft(size_t block_id)
{
    //inplace transpose
    real_data.trans(); //time*space
    //subtract mean
    double temp_mean;
    for (size_t i = 0; i < real_data.get_dim(1); i++)//loop over space
    {
        temp_mean = accumulate(real_data.get_ptr({0, i}), real_data.get_ptr({0, i + 1}), 0.);
        temp_mean /= real_data.get_dim(0);
        for (size_t j = 0; j < real_data.get_dim(0);j++)//loop over time
            real_data({j, i}) -= temp_mean;
    }
    //apply window
    if (run_input.window)
    {
        for (size_t i = 0; i < real_data.get_dim(1); i++)
            vdMul(block_size, real_data.get_ptr({0, i}), hann_array.get_ptr(), real_data.get_ptr({0, i}));
    }
    //create and set descriptor
    DFTI_DESCRIPTOR_HANDLE desc_1;
    DftiCreateDescriptor(&desc_1, DFTI_DOUBLE, DFTI_REAL, 1, block_size);
    DftiSetValue(desc_1, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(desc_1, DFTI_NUMBER_OF_TRANSFORMS, real_data.get_dim(1)); //space
    if (run_input.window)
        DftiSetValue(desc_1, DFTI_FORWARD_SCALE, sqrt((double)block_size / (double)(hann_sqr * n_realization)) / (double)block_size); //sqrt(N_f/w^2)/N_fenergy correction 1.63
    else
        DftiSetValue(desc_1, DFTI_FORWARD_SCALE, sqrt(1. / (double)(n_realization)) / (double)block_size); //sqrt(1/N_b)/N_f
    DftiSetValue(desc_1, DFTI_INPUT_DISTANCE, real_data.get_dim(0)); //time
    DftiSetValue(desc_1, DFTI_OUTPUT_DISTANCE, fft_data.get_dim(0)); //freq
    DftiSetValue(desc_1, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiCommitDescriptor(desc_1);

    //compute fft
    MKL_LONG status;
    status = DftiComputeForward(desc_1, real_data.get_ptr(), fft_data.get_ptr());
    //check state
    if (status && !DftiErrorClass(status, DFTI_NO_ERROR))
    {
        printf("Error: %s\n", DftiErrorMessage(status));
        exit(1);
    }
    DftiFreeDescriptor(&desc_1);
    dump_fft(block_id);                                              //write to disk
    real_data.reshape({real_data.get_dim(1), real_data.get_dim(0)}); //change shape back to space*time
}

void pod_spectral::calc_mode()
{
    //release memory
    real_data.setup(1);

    fft_data.setup({n_probe* run_input.fields_pod.get_len(), n_realization}); //space*block
    fft_comp.setup({fft_data.get_dim(0), fft_data.get_dim(1)}); //same as fft_data to load fft_data from file
    ndarray<MKL_Complex16> fft_temp;                            //temp array for coeff calculation in case fft_data is destroyed
    U_spectral.setup({fft_data.get_dim(0), min(fft_data.get_dim(0), fft_data.get_dim(1))});//space*block
    U.setup({U_spectral.get_dim(0),U_spectral.get_dim(1)});//same as U_spectral to write U_Spectral to file
    D.setup({min(fft_data.get_dim(0), fft_data.get_dim(1))});//block
    a_spectral.setup({fft_data.get_dim(1), U_spectral.get_dim(1)}); //block*block
    a.setup({a_spectral.get_dim(0), a_spectral.get_dim(1)}); //block*block
    ndarray<MKL_Complex16> vt_dumm(1);
    ndarray<double> superb{min(fft_data.get_dim(0), fft_data.get_dim(1)) - 1}; //min(m,n)-1
    //create result file
    create_result();

    //compute svd for each frequency
    for (size_t i = 0; i < block_size / 2 + 1; i++)
    {
        freq_id = i;
        cout << "calulating modes ... for freq: " << i + 1 << " of " << block_size / 2 + 1 << '\r' << flush;
        //load fft from file
        load_fft();
            //multiply sqrt(w) matrix
        if (run_input.coord_sys == CARTESIAN) //scale uniformly
            cblas_zdscal(fft_data.get_len(), sqrt(w(0)), fft_data.get_ptr(), 1);
        else if (run_input.coord_sys == CYLINDRICAL)
        {
            for (size_t j = 0; j < run_input.fields_pod.get_len(); j++) //loop over each field
                for (size_t i = 0; i < n_probe; i++) //loop over each probe
                    cblas_zdscal(fft_data.get_dim(1), sqrt(w(i)), fft_data.get_ptr({i + j * n_probe, 0}), fft_data.get_dim(0)); //rescale each row
        }
        //calc svd
        fft_temp = fft_data; //copy to temporary storage
        LAPACKE_zgesvd(LAPACK_COL_MAJOR, 'S', 'N', fft_data.get_dim(0), fft_data.get_dim(1),
                       fft_data.get_ptr(), fft_data.get_dim(0), D.get_ptr(), U_spectral.get_ptr(), U_spectral.get_dim(0), vt_dumm.get_ptr(), U_spectral.get_dim(1), superb.get_ptr());
        // calc coefficient X^T*U
        MKL_Complex16 alpha(sqrt(n_realization), 0);
        MKL_Complex16 beta(0, 0);
        cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, fft_temp.get_dim(1), U_spectral.get_dim(1), fft_temp.get_dim(0), &alpha, fft_temp.get_ptr(), fft_temp.get_dim(0), U_spectral.get_ptr(), U_spectral.get_dim(0), &beta, a_spectral.get_ptr(), a_spectral.get_dim(0));
        //multiply 1/sqrt(w) matrix
        if (run_input.coord_sys == CARTESIAN) //scale uniformly
            cblas_zdscal(U_spectral.get_len(), 1. / sqrt(w(0)), U_spectral.get_ptr(), 1);
        else if (run_input.coord_sys == CYLINDRICAL)
        {
            for (size_t j = 0; j < run_input.fields_pod.get_len(); j++) //loop over each field
                for (size_t i = 0; i < n_probe; i++) //loop over each probe
                    cblas_zdscal(U_spectral.get_dim(1), 1. / sqrt(w(i)), U_spectral.get_ptr({i + j * n_probe, 0}), U_spectral.get_dim(0)); //rescale each row
        }
        vdSqr(D.get_len(), D.get_ptr(), D.get_ptr());
        write_results();
    }
}

void pod_spectral::dump_fft(size_t block_id)
{
    hid_t f_id, dataspace_id, dataset_id, memspace_id;
    hsize_t dim[3], offset[3], count[3];
    string dump_fname("dump_fft.h5");
    //if file not exist
    if (block_id == 0) //for the first time
    {
        //create file, if exist,then overwrite
        f_id = H5Fcreate(dump_fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        //create datasets
        dim[2] = fft_data.get_dim(0); //freq
        dim[1] = fft_data.get_dim(1); //space
        dim[0] = n_realization;       //block
        dataspace_id = H5Screate_simple(3, dim, NULL);

        dataset_id = H5Dcreate2(f_id, "fft_real", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dataset_id);

        dataset_id = H5Dcreate2(f_id, "fft_imag", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dataset_id);

        H5Sclose(dataspace_id);
        H5Fclose(f_id);
    }

    //write to dataset
    f_id = H5Fopen(dump_fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    dim[2] = fft_data.get_dim(0); //freq
    dim[1] = fft_data.get_dim(1); //space
    dim[0] = 1;                   //1

    offset[2] = 0;        //0
    offset[1] = 0;        //0
    offset[0] = block_id; //block

    count[2] = dim[2]; //freq
    count[1] = dim[1]; //space
    count[0] = dim[0]; //1
    memspace_id = H5Screate_simple(3, dim, NULL);

    dataset_id = H5Dopen2(f_id, "fft_real", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    for (size_t i = 0; i < fft_data.get_len(); i++)
        fft_comp(i) = fft_data(i).real();
    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                            NULL, count, NULL) < 0)
        Fatal_Error("Failed to get hyperslab");
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, fft_comp.get_ptr());
    if (H5Dflush(dataset_id) < 0)
        Fatal_Error("Error writing fft_real");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    dataset_id = H5Dopen2(f_id, "fft_imag", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    for (size_t i = 0; i < fft_data.get_len(); i++)
        fft_comp(i) = fft_data(i).imag();
    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                            NULL, count, NULL) < 0)
        Fatal_Error("Failed to get hyperslab");
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, fft_comp.get_ptr());
    if (H5Dflush(dataset_id) < 0)
        Fatal_Error("Error writing fft_imag");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    H5Sclose(memspace_id);
    H5Fclose(f_id);
}

void pod_spectral::load_fft()
{
    hid_t f_id, dataset_id, memspace_id, dataspace_id;
    hsize_t dim[3]; //dimension for date of each snapshot
    hsize_t count[3];
    hsize_t offset[3];
    string dump_fname("dump_fft.h5");

    f_id = H5Fopen(dump_fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    //set subset dataset to read
    dim[2] = 1;                   //1
    dim[1] = fft_data.get_dim(0); //space
    dim[0] = fft_data.get_dim(1); //block

    offset[2] = freq_id;
    offset[1] = 0;
    offset[0] = 0;

    count[2] = dim[2];
    count[1] = dim[1];
    count[0] = dim[0];
    memspace_id = H5Screate_simple(3, dim, NULL);

    //read dataset
    dataset_id = H5Dopen2(f_id, "fft_real", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                            NULL, count, NULL) < 0)
        Fatal_Error("Failed to get hyperslab");
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, fft_comp.get_ptr());
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    for (size_t i = 0; i < fft_data.get_len(); i++)
        fft_data(i).real(fft_comp(i));

    dataset_id = H5Dopen2(f_id, "fft_imag", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                            NULL, count, NULL) < 0)
        Fatal_Error("Failed to get hyperslab");
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, fft_comp.get_ptr());
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    for (size_t i = 0; i < fft_data.get_len(); i++)
        fft_data(i).imag(fft_comp(i));

    H5Sclose(memspace_id);
    H5Fclose(f_id);
}

void pod_spectral::create_result()
{
    //create hdf5 file and dataset
    hid_t of_id, attr_id, dataspace_id, dataset_id, datatype;
    hsize_t dim[3];

    of_id = H5Fcreate(run_input.output_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //write metadata to root
    dataspace_id = H5Screate(H5S_SCALAR);
    //sampling interval
    double df = 1. / (dt * block_size);
    attr_id = H5Acreate(of_id, "df", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &df);
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
    H5Aclose(attr_id);
    H5Tclose(datatype);
    H5Sclose(dataspace_id);
    delete[] temp_field;

    //create modal energy
    dim[1] = D.get_dim(0);
    dim[0] = block_size / 2 + 1;
    dataspace_id = H5Screate_simple(2, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "modal_energy", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    //create modes
    dim[2] = U_spectral.get_dim(0); //space
    dim[1] = U_spectral.get_dim(1); //block
    dim[0] = block_size / 2 + 1;    //freq
    dataspace_id = H5Screate_simple(3, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "modes_real", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dataset_id);

    dataset_id = H5Dcreate2(of_id, "modes_imag", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    //create modal coeff
    dim[2] = a_spectral.get_dim(0); //space
    dim[1] = a_spectral.get_dim(1); //block
    dim[0] = block_size / 2 + 1;    //freq
    dataspace_id = H5Screate_simple(3, dim, NULL);
    dataset_id = H5Dcreate2(of_id, "coeff_real", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dataset_id);

    dataset_id = H5Dcreate2(of_id, "coeff_imag", H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    //close file
    H5Fclose(of_id);
}

void pod_spectral::write_results()
{
    hid_t f_id, dataset_id, dataspace_id, memspace_id;
    hsize_t dim[3], offset[3], count[3];

    f_id = H5Fopen(run_input.output_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    //write modal energy
    dim[1] = D.get_dim(0); //n_mode
    dim[0] = 1;

    offset[1] = 0;
    offset[0] = freq_id;

    count[1] = dim[1];
    count[0] = dim[0];
    memspace_id = H5Screate_simple(2, dim, NULL);

    dataset_id = H5Dopen2(f_id, "modal_energy", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                            NULL, count, NULL) < 0)
        Fatal_Error("Failed to get hyperslab");
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, D.get_ptr());
    if (H5Dflush(dataset_id) < 0)
        Fatal_Error("Error writing modes_real");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);

    //write modes
    dim[2] = U_spectral.get_dim(0); //space
    dim[1] = U_spectral.get_dim(1); //block
    dim[0] = 1;                     //1
    offset[2] = 0;
    offset[1] = 0;
    offset[0] = freq_id; //freq
    count[2] = dim[2];   //space
    count[1] = dim[1];   //block
    count[0] = dim[0];

    memspace_id = H5Screate_simple(3, dim, NULL);

    dataset_id = H5Dopen2(f_id, "modes_real", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    for (size_t i = 0; i < U_spectral.get_len(); i++)
        U(i) = U_spectral(i).real();
    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                            NULL, count, NULL) < 0)
        Fatal_Error("Failed to get hyperslab");
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, U.get_ptr());
    if (H5Dflush(dataset_id) < 0)
        Fatal_Error("Error writing modes_real");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    dataset_id = H5Dopen2(f_id, "modes_imag", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    for (size_t i = 0; i < U_spectral.get_len(); i++)
        U(i) = U_spectral(i).imag();
    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                            NULL, count, NULL) < 0)
        Fatal_Error("Failed to get hyperslab");
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, U.get_ptr());
    if (H5Dflush(dataset_id) < 0)
        Fatal_Error("Error writing modes_imag");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);

    //write coeff
    dim[2] = a_spectral.get_dim(0); //space
    dim[1] = a_spectral.get_dim(1); //block
    dim[0] = 1;                     //1
    offset[2] = 0;
    offset[1] = 0;
    offset[0] = freq_id; //freq
    count[2] = dim[2];   //space
    count[1] = dim[1];   //block
    count[0] = dim[0];
    memspace_id = H5Screate_simple(3, dim, NULL);

    //open dataset
    dataset_id = H5Dopen2(f_id, "coeff_real", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    for (size_t i = 0; i < a_spectral.get_len(); i++)
        a(i) = a_spectral(i).real();
    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                            NULL, count, NULL) < 0)
        Fatal_Error("Failed to get hyperslab");
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, a.get_ptr());
    if (H5Dflush(dataset_id) < 0)
        Fatal_Error("Error writing coeff_real");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    dataset_id = H5Dopen2(f_id, "coeff_imag", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    for (size_t i = 0; i < a_spectral.get_len(); i++)
        a(i) = a_spectral(i).imag();
    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                            NULL, count, NULL) < 0)
        Fatal_Error("Failed to get hyperslab");
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, a.get_ptr());
    if (H5Dflush(dataset_id) < 0)
        Fatal_Error("Error writing coeff_imag");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    //close file
    H5Fclose(f_id);
}
