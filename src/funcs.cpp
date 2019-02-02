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
#include "snapshot_reader.h"
#include "pod_base.h"
#include "pod_snap.h"
#include "pod_specteal.h"

void calc_classic_pod()
{
    //read snapshot file info
    snapshot_reader sr(run_input.snap_filename);
    //initialize classic pod
    pod_base pb(sr.total_n_probe, sr.total_n_snaps, sr.fields.get_len());
    //load snapshots
    sr.open_file();
    sr.partial_load_data(0, sr.total_n_probe, 0, sr.total_n_snaps, pb.get_data_ptr()); //read snapce*time
    sr.close_file();
    //calculate mean and subtract from mean
    pb.calc_mean();
    pb.subtract_mean();
    //calculate eigenvalues and eigenvectors
    pb.calc_mode();
    //write modes to hdf5
    pb.write_results();
}

void calc_snapshot_pod()
{
    //read snapshot file info
    snapshot_reader sr(run_input.snap_filename);
    //initialize snapshotpod
    pod_snap psn(sr.total_n_probe, sr.total_n_snaps, sr.fields.get_len());
    //load snapshots
    sr.open_file();
    sr.partial_load_data(0, sr.total_n_probe, 0, sr.total_n_snaps, psn.get_data_ptr()); //read snapce*time
    sr.close_file();
    //calculate mean and subtract from mean
    psn.calc_mean();
    psn.subtract_mean();
    //calculate svd
    psn.calc_mode();
    //write modes to hdf5
    psn.write_results();

}

void calc_spectral_pod()
{
    size_t n_blocks, block_size_global, overlap;
    block_size_global = run_input.block_size;
    overlap = run_input.overlap;
    //setup blocks
    n_blocks = (run_input.n_snap_read - overlap) / (block_size_global - overlap); //calculate maximum number of blocks
    //read snapshot file info
    snapshot_reader sr(run_input.snap_filename);
    //initialize spectral pod
    pod_spectral psp(sr.total_n_probe, block_size_global, sr.fields.get_len(), n_blocks);
    //load snapshot one block at a time
    sr.open_file();
    for (size_t i = 0; i < n_blocks; i++)
    {
        sr.partial_load_data(0, sr.total_n_probe, i * (block_size_global - overlap), block_size_global, psp.get_data_ptr());
        //calculate fft and store it in fft_data
        psp.calc_fft(i);
    }
    sr.close_file();
    //calculate modes
    psp.calc_mode();
    
}
