
/**
 * @file pod_base.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2019-02-02
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once
#include "global.h"
#include "util.h"//contain array, mkl complex definition and tool functions

class pod_base
{
public:
  //constructor/destructor
  pod_base() = default;
  pod_base(size_t in_n_probe, size_t in_n_snap, double in_dt);
  virtual ~pod_base() = default;

  //computation methods
  virtual void calc_mode(void); //!<calculate mode from real_data using svd(real_data)
  void calc_mean(void);  //!<calculate mean of real_data along axis=0
  void subtract_mean(void);     //!<subtract mean from real_data along axis=0
  void calculateWeight(double *in_coord);  //!<calculate the quadrature weight array 

  //output methods
  virtual void write_results(); //!<write modal energy, modes, etc. to hdf5 file
  void write_coord(ndarray<double> &in_coord); //!< write coordinates to hdf5 file

  //public data member
  ndarray<double> real_data; //space*time
  size_t n_realization, n_probe;
  double dt;

protected:
  //heavy data for pod calculation
  ndarray<double> w;         //!< quadrature weight
  ndarray<double> mean_data; //space
  ndarray<double> U, D;      //!<POD modes(space*n), modal energy(n)
  ndarray<double> a;         //!<modal coefficient(time*time)
};
