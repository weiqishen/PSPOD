/**
 * @file data_loader.h
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
#include <vector>
/**
 * @brief snapshot file reader
 * @class data_loader
 */
class data_loader
{
public:
  /**
     * @brief Construct a new data_loader object
     * 
     * @param in_fname name of snapshot file
     */
  data_loader(string in_fname);

  /**
    * @brief Destroy the snapshot object
    * 
    */
  ~data_loader();

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

  /**
 * @brief read attribution from snapshot dataset
 * 
 */
  void read_attr(void);

  /**
   * @brief read coordinates from snapshot file. 
   * If using cylindrical coordinate, transform to cylindrical coordinate system.
   * 
   */
  void read_coord(void);
  
  /**
   * @brief calculate quadrature weight
   * 
   * @param weight weight array to hold quadrature weight
   */
  void calc_weight(ndarray<double> &weight);

  /**
   * @brief partial load data from snapshot file
   * 
   * @param p_start start index of probe
   * @param n_p number of probe to read
   * @param s_start start index of snapshot
   * @param n_s number of snapshot to read
   * @param out_data data pointer to write
   */
  void partial_load_data(size_t p_start, size_t n_p, size_t s_start, size_t n_s, double *out_data);

  ndarray<double> coord;
  
private:
  //data for file reading
  string file_nameS;
  hid_t file_id;
  bool open_flag;
  vector<size_t> field_data_id;
};
