/**
 * @file pod_snap.h
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

class pod_snap : public pod_base
{
public:
  pod_snap();
  pod_snap(size_t in_n_probe, size_t in_n_snap);
  ~pod_snap();

  void calc_mode();//!<calculate mode from real_data using svd(real_data')

protected:
};
