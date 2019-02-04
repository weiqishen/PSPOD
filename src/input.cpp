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
    cout<<"Reading input file..."<<flush;
    read_param();
    if (d_xyz.get_len() != np_xyz.get_len() || np_xyz.get_len() != xyz_0.get_len())
        Fatal_Error("Inconsistent dimension for probe point parameters")
    cout<<"done."<<endl;
}

void input::read_param(void)
{
    param_reader pr(file_nameS);
    pr.openFile();
    //--------------Read params------------------
    pr.getScalarValue("task", task, int(CLASSIC_POD));
    pr.getScalarValue("snap_file", snap_filename);
    pr.getScalarValue("output_file", output_filename, string("output.h5"));
    if (run_input.output_filename.find(".") == string::npos)
        run_input.output_filename += ".h5";

    pr.getVectorValue("d_xyz", d_xyz);
    dw = 1.;
    for (size_t i = 0; i < d_xyz.get_len(); i++)
        dw *= d_xyz(i);

    pr.getVectorValue("np_xyz", np_xyz);
    pr.getVectorValue("xyz_0", xyz_0);

    if (task == SPECTRAL_POD)
    {
        pr.getScalarValue("window", window);
        pr.getScalarValue("overlap", overlap);
        pr.getScalarValue("block_size", block_size);
    }
    pr.closeFile();
}
