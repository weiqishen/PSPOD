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
    pr.getScalarValue("gamma", gamma, 1.4);
    pr.getScalarValue("output_file", output_filename, string("output.h5"));
    if (output_filename.find(".") == string::npos)
        output_filename += ".h5";

    pr.getScalarValue("norm", norm_pod, int(DEFAULT));
    if (norm_pod == DEFAULT)
    {
        pr.getVectorValue("fields", fields_pod);
        pr.getVectorValue("w_field", w_field);
        if (fields_pod.get_len() != w_field.get_len())
            Fatal_Error("w_fields must has the same length as fields!");
    }
    else if (norm_pod == SPECIFIC_KINETIC_ENERGY)
    {
        fields_pod.setup(3);
        w_field.setup(3);
        w_field = 1.;
        fields_pod(0) = "u";
        fields_pod(1) = "v";
        fields_pod(2) = "w";
    }
    else if (norm_pod == SPECIFIC_TOTAL_ENTHALPY)
    {
        fields_pod.setup(4);
        w_field.setup(4);
        fields_pod(0) = "a";
        w_field(0) = 2. / (gamma - 1.);
        fields_pod(1) = "u";
        w_field(1) = 1.;
        fields_pod(2) = "v";
        w_field(2) = 1.;
        fields_pod(3) = "w";
        w_field(3) = 1.;
    }
    else if (norm_pod == COMPRESSIBLE_ENERGY)
    {
    pr.getScalarValue("R_gas", R_gas, 286.9);
    pr.getScalarValue("Mach", Mach);
        fields_pod.setup(5);
        fields_pod(0) = "rho";
        fields_pod(1) = "u";
        fields_pod(2) = "v";
        fields_pod(3) = "w";
        fields_pod(4) = "T";
    }

    pr.getScalarValue("coord_sys", coord_sys, int(CARTESIAN));
    if (coord_sys == CARTESIAN)
        pr.getVectorValue("d_xyz", d_xyz);
    else if (coord_sys == CYLINDRICAL)
        pr.getVectorValue("d_rtz", d_xyz);
    else
        Fatal_Error("Unsupported coordinate system");

    if (task == SNAPSHOT_POD)
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
    else
        Fatal_Error("Unsupported task!");
    pr.closeFile();
}
