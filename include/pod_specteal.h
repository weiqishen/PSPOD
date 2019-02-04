/**
 * @file pod_specteal.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2019-02-02
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once
#include "pod_base.h"
class pod_spectral : public pod_base
{
public:
  pod_spectral();
  pod_spectral(size_t in_n_probe, size_t in_block_size, size_t in_n_blocks);
  ~pod_spectral();
  
  void calc_fft(size_t block_id);//!<calculate fft of real_data(time*space)
  void calc_mode();//!<calculate mode from fft_data using svd(fft_data)

  void write_results();//!<write modal energy, modes, etc. to hdf5 file

protected:
  //meta data to describe local data array
  size_t block_size_local;
  //heavy data for pod calculation
  ndarray<MKL_Complex16> fft_data;//!<(frequency*space*block)->(space*block*frequency)
  ndarray<MKL_Complex16> U_spectral;//spectral POD modes
};
