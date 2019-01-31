#pragma once
#include "global.h"
#include "ndarray.h"
#include "pod_base.h"
class pod_snap : public pod_base
{
  public:
    pod_snap();
    pod_snap(size_t n_probe, size_t n_snap, size_t n_fields);
    ~pod_snap();
    void calc_mode();

  protected:
};