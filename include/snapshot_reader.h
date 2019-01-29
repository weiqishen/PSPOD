/**
 * @file snapshot_reader.h
 * @author Weiqi Shen (weiqishen1994@ufl.edu)
 * @brief 
 * @version 0.1
 * @date 2019-01-29
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once
#include "global.h"
#include "hdf5.h"
#include "ndarray.h"
/**
 * @brief snapshot file reader
 * @class snapshot_reader
 */
class snapshot_reader
{
public:
  /**
     * @brief Construct a new snapshot_reader object
     * 
     * @param in_fname name of snapshot file
     */
  snapshot_reader(string in_fname);

  /**
    * @brief Destroy the snapshot object
    * 
    */
  ~snapshot_reader();

  /**
 * @brief Open snapshot file
 * 
 */
  void open_file(void);

  /**
    * @brief Close snapshot file
    * 
    */
  void close_file(void);

  void read_metadata(void);
  /**
   * @brief partial load data from snapshot file
   * 
   * @param p_start start index of probe
   * @param n_p number of probe to read
   * @param s_start start index of snapshot
   * @param n_s number of snapshot to read
   * @param out_data data pointer to write
   */
  void partial_load_data(int p_start, int n_p, int s_start, int n_s, double *out_data);

private:
  string file_nameS;
  hid_t file_id;


  //metadata
  int dt;
  int total_n_probe;
  int total_n_snaps;
  ndarray<string> fields;

};
