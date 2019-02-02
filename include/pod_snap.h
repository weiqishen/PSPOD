#pragma once
#include "pod_base.h"
class pod_snap : public pod_base
{
public:
  pod_snap();
  pod_snap(size_t in_n_probe, size_t in_n_snap, size_t in_n_fields);
  ~pod_snap();

  void calc_mode();//!<calculate mode from real_data using svd(real_data')

protected:
};
