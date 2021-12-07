import numpy as np
from scipy import fft
from scipy.signal import convolve
import os
import EWSPy as ews
import configuration_weekly as cfg
from pcraster import numpy2pcr, report, Scalar

### Null models adapted from (Dakos et al. 2008) ###

# TODO - method 2 does not return the right mean, other values are A-OK

## First method ##
def method1_(dataset, realizations=1, path='./1/', file_name='xxx', replace=False):
    generated_number_length = 4
    if len(str(realizations)) > 4:
        generated_number_length = len(str(realizations))

    data_shape = dataset[0].shape
    steps = np.arange(cfg.interval_map_snapshots, cfg.numberOfTimeSteps + cfg.interval_map_snapshots,
                      cfg.interval_map_snapshots)

    for k, realization in enumerate(range(realizations)):
        for data in dataset:
            data_new = data.copy()
            data_1d = data_new.ravel()
            generated_dataset_numpy = np.random.choice(data_1d, len(data_1d), replace=replace).reshape(data_shape)
            generated_dataset = numpy2pcr(Scalar, generated_dataset_numpy, np.NaN)

            generated_number_string = 'm1g' + str(realization).zfill(generated_number_length)
            dir_name = os.path.join(path + generated_number_string)

            if os.path.isdir(dir_name) == False:
                os.makedirs(dir_name)

            fname = ews.file_name_str(file_name, steps[k])
            fpath = os.path.join(dir_name, fname)
            #np.savetxt(fpath + '.numpy.txt', generated_dataset)
            report(generated_dataset, fpath)

## Second method ## TODO - Still returns real and imag number
def method2_(dataset, realizations=1, path='./1/', file_name='xxx', replace=False):
    generated_number_length = 4
    if len(str(realizations)) > 4:
        generated_number_length = len(str(realizations))

    steps = np.arange(cfg.interval_map_snapshots, cfg.numberOfTimeSteps + cfg.interval_map_snapshots,
                      cfg.interval_map_snapshots)

    for k, realization in enumerate(range(realizations)):
        for data in dataset:
            fft2_ = fft.fft2(data)
            fft2_mag = np.abs(fft2_)
            fft2_phases = np.angle(fft2_)
            fft2_phases_shape = fft2_phases.shape

            fft2_phases_new = fft2_phases.copy()
            fft2_phases_1d = fft2_phases_new.ravel()
            fft2_phases_new = np.random.choice(fft2_phases_1d, len(fft2_phases_1d), replace=replace).reshape(fft2_phases_shape)

            fft2_sym = fft2_mag * (np.cos(fft2_phases_new) + 1j * np.sin(fft2_phases_new))
            generated_dataset_numpy = fft.ifft2(fft2_sym)
            generated_dataset = numpy2pcr(Scalar, generated_dataset_numpy, np.NaN)

            generated_number_string = 'm2g' + str(realization).zfill(generated_number_length)
            dir_name = os.path.join(path + generated_number_string)

            if os.path.isdir(dir_name) == False:
                os.makedirs(dir_name)

            fname = ews.file_name_str(file_name, steps[k])
            fpath = os.path.join(dir_name, fname)
            #np.savetxt(fpath + '.numpy.txt', generated_dataset)
            report(generated_dataset, fpath)


## Third method ##
"""
Combination of Jon Yearsley (2021). Generate AR1 spatial data (https://www.mathworks.com/matlabcentral/fileexchange/5099-generate-ar1-spatial-data), MATLAB Central File Exchange. Retrieved November 30, 2021.
and Dakos et al. 10.1073/pnas.0802430105
"""

rook_neighborhood = np.array([
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 0]
])

def spatial_corr(numpy_matrix): # Moran's I, same method as used in EWSPy
    mean = np.nanmean(numpy_matrix)
    var = np.nanvar(numpy_matrix)

    numpy_matrix_mmean = np.copy(numpy_matrix)
    numpy_matrix_mmean -= mean

    is_nan = np.isnan(numpy_matrix_mmean)
    is_not_nan = ~ is_nan
    is_not_nan_as_nr = is_not_nan.astype(float)

    numpy_matrix_var = np.copy(is_not_nan_as_nr)
    numpy_matrix_var *= var

    numpy_matrix_copy = np.copy(numpy_matrix)
    numpy_matrix_copy[is_nan] = 0

    sum_neighbours = convolve(numpy_matrix_copy, rook_neighborhood, mode='same')

    n_neighbours_times_avg = convolve(is_not_nan_as_nr * mean, rook_neighborhood, mode='same')
    n_neighbours_times_avg[is_nan] = 0

    P1 = np.nansum(numpy_matrix_mmean * (sum_neighbours - n_neighbours_times_avg))
    P2 = np.nansum(sum_neighbours * numpy_matrix_var)
    return P1 / P2

def method3_(dataset, realizations=1, path='./1/', file_name='xxx', stdev_error=1.0):
    generated_number_length = 4
    if len(str(realizations)) > 4:
        generated_number_length = len(str(realizations))

    steps = np.arange(cfg.interval_map_snapshots, cfg.numberOfTimeSteps + cfg.interval_map_snapshots,
                      cfg.interval_map_snapshots)

    for k, realization in enumerate(range(realizations)):
        for data in dataset:
            Morans_I = spatial_corr(data)
            alpha0 = np.nanmean(data) * (1 - Morans_I)
            sig2 = np.nanvar(data) * (1 - Morans_I**2)

            dim = data.shape
            N = np.prod(dim)

            W = np.zeros((N, N))
            np.fill_diagonal(W[1:], Morans_I/4)
            np.fill_diagonal(W[:, 1:], Morans_I/4)

            M = np.identity(N) - W
            inv_M = np.linalg.inv(M)

            random_error = np.random.normal(loc=0.0, scale=stdev_error, size=N)
            generated_dataset_numpy = np.dot(inv_M, random_error * np.sqrt(sig2)).reshape(dim) + alpha0
            generated_dataset = numpy2pcr(Scalar, generated_dataset_numpy, np.NaN)

            generated_number_string = 'm3g' + str(realization).zfill(generated_number_length)
            dir_name = os.path.join(path + generated_number_string)

            if os.path.isdir(dir_name) == False:
                os.makedirs(dir_name)

            fname = ews.file_name_str(file_name, steps[k])
            fpath = os.path.join(dir_name, fname)
            #np.savetxt(fpath + '.numpy.txt', generated_dataset)
            report(generated_dataset, fpath)


## TESTING ##
# data = np.random.random((100, 100)) * 10
#
# #print(data)
# print(np.nanmean(data), np.nanvar(data), spatial_corr(data))
#
# bata = method1_(data)[0]
# #print(bata)
# print(np.nanmean(bata), np.nanvar(bata), spatial_corr(bata))
#
# bata = method2_(data)[0]
# #print(bata
# print(np.nanmean(bata), np.nanvar(bata), spatial_corr(bata))
#
# bata = method3_(data)[0]
# #print(bata)
# print(np.nanmean(bata), np.nanvar(bata), spatial_corr(bata))