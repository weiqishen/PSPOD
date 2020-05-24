/**
 * @file data_loader.cpp
 * @author Weiqi Shen (weiqishen1994@ufl.edu)
 * @brief 
 * @version 0.1
 * @date 2019-01-29
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "data_loader.h"
#include <sys/stat.h>

data_loader::data_loader(string in_filename)
{
    struct stat st = {0};
    if (stat(in_filename.c_str(), &st) == -1)
        Fatal_Error("data file not exist");

    file_nameS = in_filename;
    open_flag = false; //set open flag to false

    cout << "Loading snapshot file info..." << flush;
    open_file();
    read_attr();
    read_coord();
    close_file();
    cout << "done." << endl;
    cout << endl;
    cout << "Sampling interval: " << dt << " s" << endl;
    cout << "Total number of probes: " << n_probe_data << endl;
    cout << "Total number of snapshots: " << n_snap_data << endl;
    cout << "Fields from data: " << fields_data << endl;
}

data_loader::~data_loader()
{
    close_file();
}

void data_loader::open_file()
{
    if (open_flag == true) //opened
        return;

    file_id = H5Fopen(file_nameS.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    if (file_id >= 0) //success
        open_flag = true;
    else //failed
        Fatal_Error("Failed to open file");
}

void data_loader::close_file()
{
    if (open_flag) //opened
    {
        if (H5Fclose(file_id) >= 0) //success
            open_flag = false;
        else //failed
            Fatal_Error("Failed to close file");
    }
}

void data_loader::read_attr()
{
    hid_t data_id, dt_id, fields_id, space_id;
    hid_t str_ftype, str_type;
    hsize_t dim[3];
    char **temp_field;//char* array to hold field

    //open attributions
    dt_id = H5Aopen(file_id, "dt", H5P_DEFAULT);
    fields_id = H5Aopen(file_id, "fields", H5P_DEFAULT);

    //create memory space for fields
    space_id = H5Aget_space(fields_id);
    H5Sget_simple_extent_dims(space_id, dim, NULL);
    H5Sclose(space_id);
    fields_data.setup(dim[0]); //setup fields array
    temp_field = new char *[dim[0]];
    str_ftype = H5Aget_type(fields_id);
    str_type = H5Tget_native_type(str_ftype, H5T_DIR_ASCEND);

    //read attributes
    H5Aread(dt_id, H5T_NATIVE_DOUBLE, &dt);
    H5Aread(fields_id, str_type, temp_field);
    H5Tclose(str_ftype);
    H5Tclose(str_type);
    H5Aclose(dt_id);
    H5Aclose(fields_id);
    //copy back to fields
    for (size_t i = 0; i < size_t(dim[0]); i++)
    {
        fields_data(i).assign(temp_field[i]);
        delete[] temp_field[i];
    }
    delete[] temp_field;

    //set id of fields to read
    if (run_input.norm_pod == SPECIFIC_TOTAL_ENTHALPY)
    {
        if (fields_data(0) == "rho")
            field_data_id.push_back(0);
        else
            Fatal_Error("The first field in data file must be density!");
        for (size_t i = 0; i < fields_data.get_len(); i++)
        {
            if (fields_data(i) == "pressure")
            {
                pressure_id = i;
                break;
            }
            if (i == fields_data.get_len() - 1)
                Fatal_Error("Cant find pressure field in the data!");
        }
    }

    for (size_t i = field_data_id.size(); i < run_input.fields_pod.get_len(); i++) //loop over the pod fields
    {
        for (size_t j = 0; j < fields_data.get_len(); j++) //loop over the data fields
        {
            if (run_input.fields_pod(i) == fields_data(j)) //found in the data file
            {
                field_data_id.push_back(j);
                break;
            }
            else
            {
                if (j == fields_data.get_len() - 1)
                    Fatal_Error("Can't find the POD field in the data file!");
            }
        }
    }
    if (!is_sorted(field_data_id.begin(), field_data_id.end()))
        Fatal_Error("Field for POD should in the same order of field in data file!");

    //open dataset and read dimensions
    data_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
    if (data_id < 0)
        Fatal_Error("Failed to open dataset");
    space_id = H5Dget_space(data_id);
    H5Sget_simple_extent_dims(space_id, dim, NULL);
    n_probe_data = dim[1];
    n_snap_data = dim[2];
    H5Sclose(space_id);
    H5Dclose(data_id);
}

void data_loader::partial_load_data(size_t p_start, size_t n_p, size_t s_start, size_t n_s, double *out_data)
{
    hid_t dataset_id, memspace_id, dataspace_id;
    hsize_t dim[3]; //dimension for date of each snapshot
    hsize_t count[3]; 
    hsize_t offset[3];
    ndarray<double> temp_pressure;

    if (run_input.norm_pod == SPECIFIC_TOTAL_ENTHALPY)
        temp_pressure.setup(n_p);

    //set subset dataset to read
    dim[0] = run_input.fields_pod.get_len();//dimension for each read field*point*1
    dim[1] = n_p;
    dim[2] = 1;

    offset[1] = p_start; //offset of each slice field_i*p_0*t_i

    count[0] = 1;//count of each slice 1*point*1
    count[1] = n_p;
    count[2] = 1;
    //open dataset
    dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
    //set subset to read
    memspace_id = H5Screate_simple(3, dim, NULL);
    dataspace_id = H5Dget_space(dataset_id);

    for (size_t i = 0; i < n_s; i++)//for each snapshot
    {
        H5Sselect_none(dataspace_id);
        offset[2] = s_start + i;
        for (auto id : field_data_id)
        {
            offset[0] = id;
            if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_OR, offset,
                                    NULL, count, NULL) < 0)
                Fatal_Error("Failed to get hyperslab");
        }
        //read from dataset
        if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT,
                    out_data + dim[0]*dim[1] * i) < 0)
            Fatal_Error("Failed to read subset of data");
        
        //calculate speed of sound if pod field has "a" or norm=2
        if (run_input.norm_pod == SPECIFIC_TOTAL_ENTHALPY)
        {
            hid_t memspace_id2;
            offset[0] = pressure_id;
            if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                                    NULL, count, NULL) < 0)
                Fatal_Error("Failed to get hyperslab");
            //read from dataset
            memspace_id2 = H5Screate_simple(3, count, NULL);
            if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id2, dataspace_id, H5P_DEFAULT,
                        temp_pressure.get_ptr()) < 0)
                Fatal_Error("Failed to read subset of data");
            for (size_t j = 0; j < n_p; j++)
                out_data[j + dim[0] * dim[1] * i] = sqrt(run_input.gamma * temp_pressure(j) / out_data[j + dim[0] * dim[1] * i]);
            H5Sclose(memspace_id2);
        }
    }
    //close objects
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

void data_loader::read_coord()
{
    int n_dim;            //to varify coordinate array dimension
    ndarray<hsize_t> dim; //coordinate array dimension
    hid_t dataset_id, dataspace_id;

    dataset_id = H5Dopen2(file_id, "coord", H5P_DEFAULT);
    if (dataset_id < 0)
        Fatal_Error("Failed to open \'coord\'");

    dataspace_id = H5Dget_space(dataset_id);
    n_dim = H5Sget_simple_extent_ndims(dataspace_id);
    if (n_dim == 2)
        dim.setup(2);
    else
        Fatal_Error("Invalid coordinate array!");
    H5Sget_simple_extent_dims(dataspace_id, dim.get_ptr(), NULL);
    coord.setup({(size_t)dim(1), (size_t)dim(0)}); //dim*n_points

    if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord.get_ptr()) < 0)
        Fatal_Error("Failed to read subset of data");
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}
