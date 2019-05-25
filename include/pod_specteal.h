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
  void dump_fft(size_t block_id);
  void load_fft(size_t freq_id);
  void create_result();//write attributions to result file
  void write_mode(size_t freq_id);//!<write modes to file
  void write_energy();//!<write modal energy
protected:
  //meta data to describe local data array
  size_t block_size;
  //heavy data for pod calculation
  ndarray<double> hann_array;//!< hann window array
  double hann_sqr;
  ndarray<MKL_Complex16> fft_data;//!<(frequency*space)
  ndarray<MKL_Complex16> U_spectral;//spectral POD modes
};
