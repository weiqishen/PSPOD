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
    pr.getScalarValue("task", task, int(SNAPSHOT_POD));
    pr.getScalarValue("data_file", data_filename);
    pr.getScalarValue("output_file", output_filename, string("output.h5"));
    if (run_input.output_filename.find(".") == string::npos)
        run_input.output_filename += ".h5";
    pr.getVectorValue("fields", fields_pod);
    pr.getScalarValue("coord_sys",coord_sys,int(CARTESIAN));
    if (run_input.coord_sys == CARTESIAN)
        pr.getVectorValue("d_xyz", d_xyz);
    else if (run_input.coord_sys == CYLINDRICAL)
        pr.getVectorValue("d_rtz", d_xyz);
    else
        Fatal_Error("Unsupported coordinate system");

    if (task ==SNAPSHOT_POD)
    {
        pr.getScalarValue("write_mean", write_mean);
    }
    else if (task == SPECTRAL_POD)
    {
        pr.getScalarValue("window", window);
        pr.getScalarValue("overlap", overlap);
        pr.getScalarValue("block_size", block_size);
        pr.getScalarValue("from_dump", from_dump, 0);
    }
    else if (task == AZIMUTHAL_SPOD)
    {
        if (coord_sys == CARTESIAN)
            Fatal_Error("Azimuthal decomposed spectral POD only support cylindrical coordinate data.");
    }

    pr.closeFile();
}
