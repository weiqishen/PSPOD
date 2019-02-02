/**
 * @file funcs.h
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

/**
 * @brief calculate classic pod
 * @fn calc_classic_pod
 * @return int flag
 */
void calc_classic_pod(void);

/**
 * @brief calculate snapshot pod
 * @fn calc_snapshot_pod
 * @return int flag
 */
void calc_snapshot_pod(void);

/**
 * @brief calculate spectral pod
 * @fn calc_spectral_pod
 * @return int flag
 */
void calc_spectral_pod(void);
