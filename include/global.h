
/**
 * @file global.h
 * @author Weiqi Shen (weiqishen1994@ufl.edu)
 * @brief 
 * @version 0.1
 * @date 2019-01-28
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <cstring>
#include <string>
#include <cmath>
#include "err.h"
#include "input.h"

//-------------Defines---------------
enum TASK
{
    CLASSIC_POD=0,
    SNAPSHOT_POD=1,
    SPECTRAL_POD=2,
    DMD=3
};
//-------
extern double pi;
extern input run_input;