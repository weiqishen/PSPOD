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
    size_t n_snap_read;//!< number of snapshot to read from 

    // Classic POD or Snapshot POD


    //Spectral POD
    size_t overlap;//!<number of snapshots which overlap between blocks
    size_t block_size;//!< size of each block
    protected:
    std::string file_nameS;
};
