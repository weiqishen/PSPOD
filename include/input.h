/**
 * @file input.h
 * @author Weiqi Shen (weiqishen1994@ufl.edu)
 * @brief 
 * @version 0.1
 * @date 2019-01-28
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once
#include <string>
#include "ndarray.h"

class input
{
    public:
    //methods
    
    /**
     * @brief Construct a new input object
     * 
     */
    input();

    /**
     * @brief Destroy the input object
     * 
     */
    ~input();

    /**
     * @brief setup input
     * 
     * @param input_fname name of input file 
     */
    void setup(char* input_fname);

    /**
     * @brief read parameters from input file
     * 
     */
    void read_param(void);


    //---------------data--------------------

    //Common params
    int task;//!< task to do 0[default]- classic POD; 1-Snapshot POD; 2-spectral POD; 3-DMD
    string snap_filename;//!< snapshot file name
    string output_filename;//!< output file name
    ndarray<double> d_xyz;//!< space between each probe in each direction
    double dw;//!< integration weight, calculated from d_xyz
    ndarray<int> np_xyz;//!< number of probes in each direction
    ndarray<double> xyz_0;//!< position of the first probe
    // Classic POD or Snapshot POD
    int write_mean;
    //Spectral POD
    int window;
    size_t overlap;//!<number of snapshots which overlap between blocks
    size_t block_size;//!< size of each block
    //meta data from snapshot files
    double dt;
    int total_n_probe;
    int total_n_snap;
    ndarray<string> fields; //fields read from metadata
  protected:
    std::string file_nameS;
};
