/**
 * @file funcs.cpp
 * @author Weiqi Shen (weiqishen1994@ufl.edu)
 * @brief 
 * @version 0.1
 * @date 2019-01-29
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "funcs.h"
#include "ndarray.h"
#include "data_loader.h"
#include "pod_base.h"
#include "pod_spectral.h"

using namespace std;

void CalculateSpatialPOD()
{
    cout << "Start computing Snapshot POD" << endl;
    //read snapshot file info
    data_loader sr(run_input.data_filename);
    //initialize snapshotpod
    pod_base pod(sr.n_probe_data, sr.n_snap_data, sr.dt);
    //load snapshots
    sr.open_file();
    sr.partial_load_data(0, pod.n_probe, 0, pod.n_realization, pod.real_data.get_ptr());
    sr.close_file();

    //calculate weight, mean and subtract from mean
    pod.calculateWeight(sr.coord.get_ptr());
    pod.calc_mean();
    pod.subtract_mean();
    //calculate svd
    cout << "calulating modes ... " << flush;
    pod.calc_mode();
    cout << "done." << endl;
    //write modes to hdf5
    cout << "Writing data to " << run_input.output_filename << " ... " << flush;
    pod.write_results();
    pod.write_coord(sr.coord);
    cout << "done." << endl;
}

void CalculateSpectralPOD()
{
    size_t n_blocks;
    cout << "Start computing Spectral POD" << endl;
    //read snapshot file info
    data_loader sr(run_input.data_filename);
    //setup blocks
    n_blocks = (sr.n_snap_data - run_input.overlap) / (run_input.block_size - run_input.overlap); //calculate maximum number of blocks
    //initialize spectral pod
    pod_spectral pod(sr.n_probe_data, run_input.block_size, n_blocks, sr.dt);
    if (!run_input.from_dump) //calculate fft
    {
        //load snapshot one block at a time
        sr.open_file();
        for (size_t i = 0; i < n_blocks; i++)
        {
            cout << "Loading data... block " << i + 1 << " of " << n_blocks << "  \r" << flush;
            sr.partial_load_data(0, pod.n_probe, i * (pod.block_size - run_input.overlap), pod.block_size, pod.real_data.get_ptr());
            //calculate fft and store it in fft_data
            pod.calc_fft(i);
        }
        sr.close_file();
    }
    cout << endl;
    //calculate modes
    pod.calculateWeight(sr.coord.get_ptr());
    pod.calc_mode();
    cout << "done." << endl;
    //write modes to hdf5
    cout << "Writing data to " << run_input.output_filename << " ... " << flush;
    pod.write_coord(sr.coord);
    cout << "done." << endl;
}