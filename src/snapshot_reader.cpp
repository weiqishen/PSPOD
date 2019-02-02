/**
 * @file snapshot_reader.cpp
 * @author Weiqi Shen (weiqishen1994@ufl.edu)
 * @brief 
 * @version 0.1
 * @date 2019-01-29
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "snapshot_reader.h"

snapshot_reader::snapshot_reader()
{
}

snapshot_reader::snapshot_reader(string in_filename)
{
    if (in_filename.compare(in_filename.size() - 2, 2, "h5"))
        Fatal_Error("Given snapshot file not a HDF5 file");
    file_nameS = in_filename;
    open_flag = false; //set open flag to false
    open_file();
    read_attr();
    close_file();
}

snapshot_reader::~snapshot_reader()
{
    close_file();
}

void snapshot_reader::open_file()
{
    if (open_flag == true) //opened
        return;

    file_id = H5Fopen(file_nameS.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    if (file_id >= 0) //success
        open_flag = true;
    else //failed
        Fatal_Error("Failed to open file");
}

void snapshot_reader::close_file()
{
    if (open_flag) //opened
    {
        if (H5Fclose(file_id) >= 0) //success
            open_flag = false;
        else //failed
            Fatal_Error("Failed to close file");
    }
}

void snapshot_reader::read_attr()
{
    hid_t data_id, dt_id, n_probe_id, n_snaps_id, fields_id, n_fields_id;
    hid_t ftype, type;
    int n_fields;
    char **temp;
    //open dataset
    data_id = H5Dopen2(file_id, "snapshots", H5P_DEFAULT);
    if (data_id < 0)
        Fatal_Error("Failed to open dataset");

    //open attributions
    dt_id = H5Aopen(data_id, "dt", H5P_DEFAULT);
    n_probe_id = H5Aopen(data_id, "num_probe", H5P_DEFAULT);
    n_snaps_id = H5Aopen(data_id, "n_snaps", H5P_DEFAULT);
    n_fields_id = H5Aopen(data_id, "n_fields", H5P_DEFAULT);
    fields_id = H5Aopen(data_id, "fields", H5P_DEFAULT);

    //create memory space for fields

    ftype = H5Aget_type(fields_id);
    type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
    H5Aread(n_fields_id, H5T_NATIVE_INT32, &n_fields);
    fields.setup(n_fields); //setup fields array
    temp = new char *[n_fields];

    //read attributes
    H5Aread(dt_id, H5T_NATIVE_DOUBLE, &dt);
    H5Aread(n_probe_id, H5T_NATIVE_INT32, &total_n_probe);
    H5Aread(n_snaps_id, H5T_NATIVE_INT32, &total_n_snaps);
    H5Aread(fields_id, type, temp);

    //copy back to fields
    for (size_t i = 0; i < size_t(n_fields); i++)
    {
        fields(i).assign(temp[i]);
        delete[] temp[i];
    }
    //close objects
    H5Tclose(ftype);
    H5Tclose(type);
    H5Aclose(dt_id);
    H5Aclose(n_probe_id);
    H5Aclose(n_snaps_id);
    H5Aclose(fields_id);
    delete[] temp;
}

void snapshot_reader::partial_load_data(size_t p_start, size_t n_p, size_t s_start, size_t n_s, double *out_data)
{
    hid_t dataset_id, memspace_id, dataspace_id;
    ndarray<hsize_t> dimsm{3};
    hsize_t count[3];  // number of blocks
    hsize_t offset[3]; // offset of the sub data
    hsize_t stride[3]; // stride
    hsize_t block[3];  //block size

    //set subset dataset to read
    dimsm(0) = fields.get_len();
    dimsm(1) = n_p;
    dimsm(2) = n_s;

    offset[0] = 0;       //field offset
    offset[1] = p_start; //probe offset
    offset[2] = s_start; //snapshot offset

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
    dataset_id = H5Dopen2(file_id, "/snapshots", H5P_DEFAULT);
    //set subset to read
    memspace_id = H5Screate_simple(3, dimsm.get_ptr(), NULL);
    dataspace_id = H5Dget_space(dataset_id);
    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                            stride, count, block) < 0)
        Fatal_Error("Failed to get hyperslab");

    //read from dataset
    if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT,
                out_data) < 0)
        Fatal_Error("Failed to read subset of data");

    //close objects
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}
