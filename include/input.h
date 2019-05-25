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

enum TASK
{
  CLASSIC_POD = 0,
  SNAPSHOT_POD = 1,
  SPECTRAL_POD = 2,
};

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
  void setup(char *input_fname);

  /**
     * @brief read parameters from input file
     * 
     */
  void read_param(void);

  //---------------data--------------------

  //***Common params***
  int task;               //!< task to do 0[default]- classic POD; 1-Snapshot POD; 2-spectral POD; 3-DMD
  string data_filename;   //!< data file name
  string output_filename; //!< output file name
  ndarray<double> d_xyz; //!< space between each probe in each direction
  ndarray<string> fields_pod;//!< fields need to perform pod

  // Classic POD or Snapshot POD
  int write_mean;     //!< if write mean data to output
  //Spectral POD
  int window;         //!<whether apply windowing before fft
  size_t overlap;    //!<number of snapshots which overlap between blocks
  size_t block_size; //!< size of each block
  int from_dump;

  //meta data from snapshot files
  double dt; //!<time step size
  int n_probe_global; //!< total number of probes
  int n_snap_global;  //!< total number of snapshots
  ndarray<string> fields_data; //!<fields in data file

protected:
  std::string file_nameS;
};
