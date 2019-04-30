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
    cout << "Reading input file..." << flush;
    read_param();
    cout << "done." << endl;
    cout << "Field to perform POD" << fields_pod << endl;
}

void input::read_param(void)
{
    param_reader pr(file_nameS);
    pr.openFile();
    //--------------Read params------------------
    pr.getScalarValue("task", task, int(CLASSIC_POD));
    pr.getScalarValue("data_file", data_filename);
    pr.getScalarValue("output_file", output_filename, string("output.h5"));
    if (run_input.output_filename.find(".") == string::npos)
        run_input.output_filename += ".h5";
    pr.getVectorValue("fields", fields_pod);
    pr.getVectorValue("d_xyz", d_xyz);

    if (task == SPECTRAL_POD)
    {
        pr.getScalarValue("window", window);
        pr.getScalarValue("overlap", overlap);
        pr.getScalarValue("block_size", block_size);
    }
    else
    {
        pr.getScalarValue("write_mean", write_mean);
    }
    pr.closeFile();
}
