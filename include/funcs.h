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
 * @brief calculate snapshot pod
 */
void CalculateSpatialPOD(void);

/**
 * @brief calculate spectral pod
 */
void CalculateSpectralPOD(void);

/**
 * @brief calculate azimuthal decomposed spectral pod
 */
void CalculateAzimuthalSpectralPOD(void);
