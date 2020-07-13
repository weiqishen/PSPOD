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
  SNAPSHOT_POD = 0,
  SPECTRAL_POD = 1,
};

enum COORD_SYS
{
  CARTESIAN = 0,
  CYLINDRICAL = 1,
};

enum NORM_TYPE
{
  DEFAULT = 0,
  SPECIFIC_KINETIC_ENERGY = 1,
  SPECIFIC_TOTAL_ENTHALPY = 2,
};

class input
{
public:
  //methods

  /**
     * @brief Construct a new input object
     * 
     */
  input() = default;

  /**
     * @brief Destroy the input object
     * 
     */
  ~input() = default;

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
  int task;                   //!< task to do 0[default]-Snapshot POD; 1-spectral POD; 2-azimuthal decomposed spectral POD
  string data_filename;       //!< data file name
  string output_filename;     //!< output file name
  ndarray<double> d_xyz;      //!< space between each probe in each direction
  ndarray<string> fields_pod; //!< fields need to perform pod
  int norm_pod;               //!< type of norm to perform pod. 0[default]-manually specify fields; 1-specific kinetic energy; 2-specific total enthalpy
  ndarray<double> w_field;    //!<weight for each pod field variable
  int coord_sys;               //!< coordinated system 0: cartesian; 1: cylindrical
  double gamma;               //ratio of specific heat
  //Snapshot POD
  int write_mean; //!< if write mean data to output
  //Spectral POD
  int window;        //!< whether apply windowing before fft
  size_t overlap;    //!< number of snapshots which overlap between blocks
  size_t block_size; //!< size of each block
  int from_dump;     //!< if load fft data from dumped temp file

protected:
  std::string file_nameS;
};
