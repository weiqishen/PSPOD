/**
 * @file hdf_rw.cpp
 * @author Weiqi Shen (weiqishen1994@ufl.edu)
 * @brief 
 * @version 0.1
 * @date 2019-01-29
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "snapshot_reader.h"

snapshot_reader::snapshot_reader(string in_filename)
{
    if (H5Fis_hdf5(in_filename.c_str()) > 0)
        Fatal_Error("Given snapshot file format is not HDF5");
    file_nameS = in_filename;
}

snapshot_reader::~snapshot_reader()
{
    close_file();
}

void snapshot_reader::open_file()
{
    file_id = H5Fopen(file_nameS.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0)
        Fatal_Error("Failed to open file")
}

void snapshot_reader::close_file()
{
    if (H5Fclose(file_id) < 0)
        Fatal_Error("Failed to close file");
}

void snapshot_reader::read_metadata()
{

}

void snapshot_reader::partial_load_data(int p_start,int n_p,int s_start,int n_s,double *out_data)
{
    hid_t dataset_id,memspace_id,dataspace_id;
    herr_t status;
    ndarray<hsize_t>dimsm{3};
    hsize_t count[3];  // number of blocks
    hsize_t offset[3]; // offset of the sub data
    hsize_t stride[3]; // stride
    hsize_t block[3];   //block size

    //set subset dataset to read
    dimsm(0) = n_s;
    dimsm(1) = fields.get_len();
    dimsm(2) = n_p;

    offset[0] = s_start; //snapshot offset
    offset[1] = 0;       //field offset according to metadata
    offset[2] = p_start; //probe offset

    count[0] = dimsm(0);
    count[1] = dimsm(1);
    count[2] = dimsm(2);

    stride[0] = 1;
    stride[1] = 1;
    stride[2] = 1;

    block[0] = 1;
    block[1] = 1;
    block[2] = 1;

    //open dataset
    hid_t dataset_id = H5Dopen2(file_id, "/snapshots", H5P_DEFAULT);
    //set subset to read
    memspace_id = H5Screate_simple (3, dimsm.get_ptr(), NULL);
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset,
                                  stride, count, block);

    //read from dataset
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     out_data);

}