/**
 * @file input.cpp
 * @class input
 * @author Weiqi Shen weiqishen1994@ufl.edu
 * @brief input reader
 * @version 0.1
 * @date 2019-01-28
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include "global.h"
#include "param_reader.h"

using namespace std;

input::input()
{
}

input::~input()
{
}

void input::setup(char *input_fname)
{
    file_nameS.assign(input_fname);
    read_param();
}

void input::read_param(void)
{
    param_reader pr(file_nameS);
    pr.openFile();
    //--------------Read params------------------
    pr.getScalarValue("task", task, int(CLASSIC_POD));
    pr.getScalarValue("snap_file", snap_filename);
    pr.getScalarValue("n_snap_read", n_snap_read);

    if (task == SPECTRAL_POD)
    {
        pr.getScalarValue("overlap", overlap);
        pr.getScalarValue("block_size", block_size);
    }
    pr.closeFile();
}