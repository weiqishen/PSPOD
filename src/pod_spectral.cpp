#include <numeric>
#include "pod_specteal.h"

using namespace std;

pod_spectral::pod_spectral(size_t in_n_probe, size_t in_block_size, size_t in_n_fields, size_t in_n_blocks)
{
    n_realization_local = in_n_blocks;
    n_probe_local = in_n_probe;
    n_fields = in_n_fields;
    block_size_local = in_block_size;
    real_data.setup({block_size_local, n_probe_local * n_fields});
    fft_data.setup({block_size_local / 2 + 1, n_probe_local * n_fields, n_realization_local}); 
}

pod_spectral::~pod_spectral()
{
}

void pod_spectral::calc_fft(size_t block_id)
{
    //create and set descriptor
    DFTI_DESCRIPTOR_HANDLE desc_1;
    DftiCreateDescriptor(&desc_1, DFTI_DOUBLE, DFTI_REAL, 1, block_size_local);
    DftiSetValue(desc_1, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(desc_1, DFTI_NUMBER_OF_TRANSFORMS, real_data.get_dim(1));
    DftiSetValue(desc_1, DFTI_FORWARD_SCALE, 1. / block_size_local);
    DftiSetValue(desc_1, DFTI_INPUT_DISTANCE, block_size_local);
    DftiSetValue(desc_1, DFTI_OUTPUT_DISTANCE, block_size_local / 2 + 1);
    DftiSetValue(desc_1, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiCommitDescriptor(desc_1);

    //compute fft
    MKL_LONG status;
    status = DftiComputeForward(desc_1, real_data.get_ptr(), fft_data.get_ptr({0, 0, block_id}));
    //check state
    if (status && !DftiErrorClass(status, DFTI_NO_ERROR))
    {
        printf("Error: %s\n", DftiErrorMessage(status));
        exit(1);
    }
    DftiFreeDescriptor(&desc_1);
}

void pod_spectral::calc_mode()
{
    //scale the amplitude of half spectrum
    for (size_t i = 0; i < fft_data.get_len() / fft_data.get_dim(0); i++)//loop over space*block
    {
        transform(fft_data.get_ptr({1, i}), fft_data.get_ptr({block_size_local - block_size_local / 2, i}),
                  fft_data.get_ptr({1, i}), [](MKL_Complex16 x) { return x * 2.; });
    }
    //reshape fft data ->space*block*frequency
    fft_data.reshape({block_size_local, n_probe_local * n_fields * n_realization_local});
    fft_data.trans();
    fft_data.reshape({n_probe_local * n_fields, n_realization_local, block_size_local});
    //calculate svd
    U_spectral({fft_data.get_dim(0), fft_data.get_dim(1),fft_data.get_dim(2)});
    D.setup({fft_data.get_dim(1),fft_data.get_dim(2)});
    //declare dummy array
    ndarray<MKL_Complex16> vt_dumm(1);
    ndarray<double> superb{min(fft_data.get_dim(0), fft_data.get_dim(1)) - 1};//min(m,n)-1
    //compute svd for each frequency
    for (size_t i = 0; i < block_size_local; i++)
    {
        LAPACKE_zgesvd(LAPACK_COL_MAJOR, 'S', 'N', fft_data.get_dim(0), fft_data.get_dim(1),
                       fft_data.get_ptr({0, 0, i}), fft_data.get_dim(0), D.get_ptr({0, i}), U_spectral.get_ptr({0, 0, i}),
                       U_spectral.get_dim(0), vt_dumm.get_ptr(), fft_data.get_dim(1), superb.get_ptr());
    }
    //scale singular value to modal energy
    vdSqr(D.get_len(), D.get_ptr(), D.get_ptr());
}

void pod_base::write_results()
{
    //create hdf5 file

    //write metadata

    //write modal energy

    //write modes
}
