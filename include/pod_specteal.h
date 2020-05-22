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
  //constructor/destructor
  ~pod_spectral() = default;
  pod_spectral(size_t in_n_probe, size_t in_block_size, size_t in_n_blocks, double in_dt);

  void calc_fft(size_t block_id); //!<calculate fft of real_data(time*space)
  void dump_fft(size_t block_id);
  void load_fft();
  void calc_mode(); //!<calculate mode from fft_data using svd(fft_data)
  void create_result(); //write attributions to result file
  void write_results(); //!<write results to file

  //meta data to describe local data array
  size_t block_size;
protected:
  //heavy data for pod calculation
  ndarray<double> hann_array; //!< hann window array
  double hann_sqr;
  ndarray<MKL_Complex16> fft_data;   //!<(frequency*space)
  ndarray<MKL_Complex16> U_spectral; //!<spectral POD modes
  ndarray<MKL_Complex16> a_spectral; //!<(block*block)
  ndarray<double> fft_comp;          //!<temp array to load fft data

private:
  size_t freq_id; //a counter to determine which frequency to write to
};
