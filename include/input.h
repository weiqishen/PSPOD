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

enum COORD
{
  CARTESIAN_COORD = 0,
  CYLINDRICAL_COORD = 1
};

enum TASK
{
  CLASSIC_POD = 0,
  SNAPSHOT_POD = 1,
  SPECTRAL_POD = 2,
  DMD = 3
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
  string snap_filename;   //!< snapshot file name
  string output_filename; //!< output file name

  int coord_type; //!< type of coordinate system 0: cartesian; 1: cylindrical
  //cartesian
  ndarray<double> d_xyz; //!< space between each probe in each direction
  //cylindrical
  ndarray<double> z_axis;//!< axis of z direction in cylindrical coord

  // Classic POD or Snapshot POD
  int write_mean;
  //Spectral POD
  int window;
  size_t overlap;    //!<number of snapshots which overlap between blocks
  size_t block_size; //!< size of each block

  //meta data from snapshot files
  double dt;
  int n_probe;
  int n_snap_global;
  ndarray<string> fields;

protected:
  std::string file_nameS;
};
