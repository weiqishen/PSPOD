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

snapshot_reader::snapshot_reader(string in_filename)
{
    if (in_filename.compare(in_filename.size() - 2, 2, "h5"))
        Fatal_Error("Given snapshot file not a HDF5 file");
    file_nameS = in_filename;
    open_flag = false; //set open flag to false
    cout<<"Loading snapshot file info..."<<flush;
    open_file();
    read_attr();
    close_file();
    cout << "done." << endl;
    cout << endl;
    cout << "Sampling interval: " << run_input.dt << "sec" << endl;
    cout << "Total number of probes: " << run_input.n_probe << endl;
    cout << "Total number of snapshots: " << run_input.n_snap_global << endl;
    cout<<"Fields: "<<run_input.fields<<endl;
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
    hid_t data_id, dt_id, n_probe_id, n_snap_id, fields_id, space_id;
    hid_t str_ftype, str_type;
    hsize_t n_fields;
    char **temp_field;//char* array to hold field
    //open dataset
    data_id = H5Dopen2(file_id, "/snapshots", H5P_DEFAULT);
    if (data_id < 0)
        Fatal_Error("Failed to open dataset");

    //open attributions
    dt_id = H5Aopen(data_id, "dt", H5P_DEFAULT);
    n_probe_id = H5Aopen(data_id, "n_probe", H5P_DEFAULT);
    n_snap_id = H5Aopen(data_id, "n_snaps", H5P_DEFAULT);
    fields_id = H5Aopen(data_id, "fields", H5P_DEFAULT);

    //create memory space for fields
    space_id = H5Aget_space(fields_id);
    H5Sget_simple_extent_dims(space_id, &n_fields, NULL);
    run_input.fields.setup(n_fields); //setup fields array
    temp_field = new char *[n_fields];
    str_ftype = H5Aget_type(fields_id);
    str_type = H5Tget_native_type(str_ftype, H5T_DIR_ASCEND);

    //read attributes
    H5Aread(dt_id, H5T_NATIVE_DOUBLE, &run_input.dt);
    H5Aread(n_probe_id, H5T_NATIVE_INT32, &run_input.n_probe);
    H5Aread(n_snap_id, H5T_NATIVE_INT32, &run_input.n_snap_global);
    H5Aread(fields_id, str_type, temp_field);

    //copy back to fields
    for (size_t i = 0; i < size_t(n_fields); i++)
    {
        run_input.fields(i).assign(temp_field[i]);
        delete[] temp_field[i];
    }
    //close objects
    H5Tclose(str_ftype);
    H5Tclose(str_type);
    H5Sclose(space_id);
    H5Aclose(dt_id);
    H5Aclose(n_probe_id);
    H5Aclose(n_snap_id);
    H5Aclose(fields_id);
    H5Dclose(data_id);
    delete[] temp_field;
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
    dimsm(0) = run_input.fields.get_len();
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

void snapshot_reader::read_coord()
{
    int n_dim;            //to varify coordinate array dimension
    ndarray<hsize_t> dim; //coordinate array dimension
    hid_t dataset_id, dataspace_id;

    dataset_id = H5Dopen2(file_id, "/coord", H5P_DEFAULT);
    if (dataset_id < 0)
        Fatal_Error("Failed to open /coord");

    dataspace_id = H5Dget_space(dataset_id);
    n_dim = H5Sget_simple_extent_ndims(dataspace_id);
    if (n_dim == 2)
        dim.setup(2);
    else
        Fatal_Error("Invalid coordinate array!");
    H5Sget_simple_extent_dims(dataspace_id, dim.get_ptr(), NULL);
    coord.setup({(size_t)dim(1), (size_t)dim(0)}); //dimension inverse

    if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord.get_ptr()) < 0)
        Fatal_Error("Failed to read subset of data");
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

void snapshot_reader::calc_weight(ndarray<double> &weight)
{
    weight.setup(run_input.n_probe * run_input.fields.get_len()); //probe*field
    if (run_input.coord_type == CARTESIAN_COORD)
    {
        double dw = 1.;
        for (size_t i = 0; i < run_input.d_xyz.get_len(); i++)
            dw *= run_input.d_xyz(i);
        weight = dw;
    }
    else if (run_input.coord_type == CYLINDRICAL_COORD)
    {
        Fatal_Error("Cylindrical coordinate weight hasn't been implemented");
    }
}