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
    //load data
    snapshot_reader sr(run_input.snap_filename);
    pod_base pb(sr.total_n_probe,sr.total_n_snaps,sr.fields.get_len());
    sr.open_file();
    sr.partial_load_data(0, sr.total_n_probe, 0, sr.total_n_snaps, pb.get_data_ptr()); //read snap*pts
    sr.close_file();
    //calculate mean and subtract from mean
    pb.calc_mean();
    pb.subtract_mean();
    //calculate eigenvalues and eigenvectors
    pb.calc_mode();
    //write modes to hdf5
}

void calc_snapshot_pod()
{
    //load data
    snapshot_reader sr(run_input.snap_filename);
    pod_snap psn(sr.total_n_probe,sr.total_n_snaps,sr.fields.get_len());
    sr.open_file();
    sr.partial_load_data(0, sr.total_n_probe, 0, sr.total_n_snaps, psn.get_data_ptr()); //read snap*pts
    sr.close_file();
    //calculate mean and subtract from mean
    psn.calc_mean();
    psn.subtract_mean();
    //calculate svd
    psn.calc_mode();
    //write modes to hdf5

}

void calc_spectral_pod()
{
    size_t num_blocks, block_size, overlap;
    block_size = run_input.block_size;
    overlap = run_input.overlap;
    //setup blocks
    num_blocks = (run_input.n_snap_read - overlap) / (block_size - overlap); //maximum number of blocks
    //load all points to compose the blocks
    snapshot_reader sr(run_input.snap_filename);

    ndarray<double> temp_real{block_size,size_t(sr.total_n_probe)};
    sr.open_file();
    for (size_t i = 0; i < num_blocks; i++)
    {
        sr.partial_load_data(0, sr.total_n_probe, i * (block_size - overlap), block_size, temp_real.get_ptr());
        //do fft

        //save to hdf5
    }
        sr.close_file();

    //load 
    //do svd/array multiplication

}
