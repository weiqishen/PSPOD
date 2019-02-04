/**
 * @file err.h
 * @author Weiqi Shen (weiqishen1994@ufl.edu)
 * @brief error handler
 * @version 0.1
 * @date 2019-01-28
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once
#include <stdio.h>

#define Fatal_Error(ss)                                                \
    {                                                                  \
        printf("Fatal error '%s' at %s:%d\n", ss, __FILE__, __LINE__); \
        exit(1);                                                       \
    }
    