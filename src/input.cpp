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

    pr.getScalarValue("coord_type", coord_type, int(CARTESIAN_COORD));
    if (coord_type == CARTESIAN_COORD)
    {
        pr.getVectorValue("d_xyz", d_xyz);
    }
    else if (coord_type == CYLINDRICAL_COORD)
    {
        pr.getVectorValue("axis", z_axis);
    }
    else
    {
        Fatal_Error("Coordinate system not supported!");
    }

    if (task == SPECTRAL_POD)
    {
        pr.getScalarValue("window", window);
        pr.getScalarValue("overlap", overlap);
        pr.getScalarValue("block_size", block_size);
    }
    pr.closeFile();
}
