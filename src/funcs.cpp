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
#include "pod_snap.h"
#include "pod_specteal.h"

using namespace std;

void calc_classic_pod()
{
    cout << "Start computing Classic POD" << endl;
    //read snapshot file info
    data_loader sr(run_input.data_filename);
    //initialize classic pod
    pod_base pb(run_input.n_probe_global, run_input.n_snap_global);
    //load snapshots
    sr.open_file();
    sr.read_coord();
    sr.calc_weight(pb.w);
    sr.partial_load_data(0, run_input.n_probe_global, 0, run_input.n_snap_global, pb.real_data.get_ptr());
    sr.close_file();
    //calculate mean and subtract from mean
    pb.calc_mean();
    pb.subtract_mean();
    //calculate eigenvalues and eigenvectors
    cout << "calulating modes ... " << flush;
    pb.calc_mode();
    cout << "done." << endl;
    //write modes to hdf5
    cout << "Writing data to " << run_input.output_filename << " ... " << flush;
    pb.write_results();
    pb.write_coord(sr.coord);
    cout << "done." << endl;
}

void calc_snapshot_pod()
{
    cout << "Start computing Snapshot POD" << endl;
    //read snapshot file info
    data_loader sr(run_input.data_filename);
    //initialize snapshotpod
    pod_snap psn(run_input.n_probe_global, run_input.n_snap_global);
    //load snapshots
    sr.open_file();
    sr.read_coord();
    sr.calc_weight(psn.w);
    sr.partial_load_data(0, run_input.n_probe_global, 0, run_input.n_snap_global, psn.real_data.get_ptr());
    sr.close_file();
    //calculate mean and subtract from mean
    psn.calc_mean();
    psn.subtract_mean();
    //calculate svd
    cout << "calulating modes ... " << flush;
    psn.calc_mode();
    cout << "done." << endl;
    //write modes to hdf5
    cout << "Writing data to " << run_input.output_filename << " ... " << flush;
    psn.write_results();
    psn.write_coord(sr.coord);
    cout << "done." << endl;
}

void calc_spectral_pod()
{
    size_t n_blocks;
    cout << "Start computing Spectral POD" << endl;
    //read snapshot file info
    data_loader sr(run_input.data_filename);
    //setup blocks
    n_blocks = (run_input.n_snap_global - run_input.overlap) / (run_input.block_size - run_input.overlap); //calculate maximum number of blocks
    //initialize spectral pod
    pod_spectral psp(run_input.n_probe_global, run_input.block_size, n_blocks);
    //load snapshot one block at a time
    sr.open_file();
    sr.read_coord();
    sr.calc_weight(psp.w);
    if (!run_input.from_dump)
    {
        for (size_t i = 0; i < n_blocks; i++)
        {
            cout << "Loading data... block " << i + 1 << " of " << n_blocks << "  \r" << flush;
            sr.partial_load_data(0, run_input.n_probe_global, i * (run_input.block_size - run_input.overlap), run_input.block_size, psp.real_data.get_ptr());
            //calculate fft and store it in fft_data
            psp.calc_fft(i);
        }
    }
    sr.close_file();
    cout << endl;
    //calculate modes
    psp.calc_mode();
    cout << "done." << endl;
    //write modes to hdf5
    cout << "Writing data to " << run_input.output_filename << " ... " << flush;
    psp.write_coord(sr.coord);
    cout << "done." << endl;
}
