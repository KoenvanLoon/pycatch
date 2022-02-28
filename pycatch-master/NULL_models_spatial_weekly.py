import numpy as np
from scipy import fft
from scipy.signal import convolve
from scipy import ndimage
import os
import EWSPy as ews
import EWS_configuration as cfg
from pcraster import numpy2pcr, report, Scalar


### Null models adapted from (Dakos et al. 2008) ###

# TODO - method 2 did not always return the right mean --> check solution -, check method 3

# Detrend dataset
"""
Detrends the given dataset making use of either Gaussian filtering or linear detrending, as specified in the 
configuration. Optionally, this data is saved.

Args:
-----

data : numpy array, the spatial datasets.

realizations : int, the number of datasets generated. Used for folder name in the case of detrending.

path : str, the filepath where the original dataset can be found.

file_name : str, name of the variable.

Return:
-----

detrended_data : The detrended spatial datasets.

"""


def detrend_(dataset, realizations=1, path='./1/', variable='xxx'):
    generated_number_length = ews.generated_number_length(realizations)

    steps = np.arange(cfg.interval_map_snapshots, cfg.number_of_timesteps_weekly + cfg.interval_map_snapshots,
                      cfg.interval_map_snapshots)

    detrended_dataset = [0.0] * len(dataset)

    for k, data in enumerate(dataset):
        detrended_data = data
        if cfg.detrended_method == 'Gaussian':
            gaussian_filter = ndimage.gaussian_filter(data, cfg.detrended_sigma)
            detrended_data -= gaussian_filter
        elif cfg.detrended_method == 'Linear':
            mean = np.nanmean(data)
            detrended_data -= mean
        elif cfg.detrended_method is not 'None':
            print("Invalid input for detrending_spat in generate_datasets (EWS_weekly.py). No detrending done.")

        detrended_dataset[k] = detrended_data

        if cfg.save_detrended_data:
            generated_number_string = 'dtr' + str(0).zfill(generated_number_length)
            dir_name = os.path.join(path + generated_number_string)

            if os.path.isdir(dir_name) == False:
                os.makedirs(dir_name)

            fname1 = ews.file_name_str(variable.name, steps[k])
            fpath1 = os.path.join(dir_name, fname1)
            # np.savetxt(fpath1 + '.numpy.txt', detrended_data)
            detrended_data_pcr = numpy2pcr(Scalar, detrended_data, np.NaN)
            report(detrended_data_pcr, fpath1)

            if cfg.detrended_method == 'Gaussian':
                fname2 = ews.file_name_str(variable.name + 'g', steps[k])
                fpath2 = os.path.join(dir_name, fname2)
                # np.savetxt(fpath2 + '.numpy.txt', gaussian_filter)
                gaussian_filter_pcr = numpy2pcr(Scalar, gaussian_filter, np.NaN)
                report(gaussian_filter_pcr, fpath2)

    return detrended_dataset


# Generate datasets method 1
"""
Generates dataset(s) with similar mean and variance by randomly picking values from the original dataset. In the case
of replace==False, this is similar to shuffling the dataset.

Args:
-----

data : numpy array, the spatial datasets.

realizations : int, the number of datasets generated.

path : str, the filepath where the original dataset can be found.

variable : str, name of the variable.

replace : bool, selects whether new values are picked from the original dataset or the original dataset minus previously
    picked values. Usually set to False to ensure similar mean and variance for smaller datasets.

"""


def method1_(dataset, realizations=1, path='./1/', variable='xxx', replace=False):
    generated_number_length = ews.generated_number_length(realizations)

    data_shape = dataset[0].shape
    steps = np.arange(cfg.interval_map_snapshots, cfg.number_of_timesteps_weekly + cfg.interval_map_snapshots,
                      cfg.interval_map_snapshots)

    for k, data in enumerate(dataset):
        data_new = data.copy()
        data_1d = data_new.ravel()
        for realization in range(realizations):
            generated_dataset_numpy = np.random.choice(data_1d, len(data_1d), replace=replace).reshape(data_shape)
            generated_dataset = numpy2pcr(Scalar, generated_dataset_numpy, np.NaN)

            generated_number_string = 'm1g' + str(realization).zfill(generated_number_length)
            dir_name = os.path.join(path + generated_number_string)

            if os.path.isdir(dir_name) == False:
                os.makedirs(dir_name)

            fname = ews.file_name_str(variable.name, steps[k])
            fpath = os.path.join(dir_name, fname)
            # np.savetxt(fpath + '.numpy.txt', generated_dataset)
            report(generated_dataset, fpath)


# Generate datasets method 2
"""
Generates dataset(s) with similar autocorrelation, mean and variance by generating datasets with the same Fourier 
spectrum and amplitudes.

Args:
-----

data : numpy array, the spatial datasets.

realizations : int, the number of datasets generated.

method : str, either 'None' or 'Detrending', if detrended data is used as input, no further detrending is needed. If
    not-detrended data is used, linear detrending is applied before the Fourier spectrum and amplitudes are calculated,
    with the linear detrend added after the generation of datasets. 

path : str, the filepath where the original dataset can be found.

variable : name of the variable.

replace : bool, selects whether new values are picked from the original dataset or the original dataset minus previously
    picked values.

"""


def method2_(dataset, realizations=1, method='Detrending', path='./1/', variable='xxx', replace=False):
    generated_number_length = ews.generated_number_length(realizations)

    steps = np.arange(cfg.interval_map_snapshots, cfg.number_of_timesteps_weekly + cfg.interval_map_snapshots,
                      cfg.interval_map_snapshots)

    for k, data in enumerate(dataset):
        if method == 'Detrending':
            mean = np.nanmean(data)
            data = data - mean

        fft2_ = fft.fft2(data)
        fft2_mag = np.abs(fft2_)
        fft2_phases = np.angle(fft2_)
        fft2_phases_shape = fft2_phases.shape

        fft2_phases_new = fft2_phases.copy()
        fft2_phases_1d = fft2_phases_new.ravel()
        for realization in range(realizations):
            fft2_phases_new = np.random.choice(fft2_phases_1d, len(fft2_phases_1d), replace=replace).reshape(
                fft2_phases_shape)

            fft2_sym = fft2_mag * (np.cos(fft2_phases_new) + 1j * np.sin(fft2_phases_new))
            generated_dataset_numpy = fft.ifft2(fft2_sym)

            generated_dataset_numpy = np.absolute(generated_dataset_numpy)  # TODO - Check if this is correct

            generated_dataset = numpy2pcr(Scalar, generated_dataset_numpy, np.NaN)

            generated_number_string = 'm2g' + str(realization).zfill(generated_number_length)
            dir_name = os.path.join(path + generated_number_string)

            if os.path.isdir(dir_name) == False:
                os.makedirs(dir_name)

            fname = ews.file_name_str(variable.name, steps[k])
            fpath = os.path.join(dir_name, fname)
            # np.savetxt(fpath + '.numpy.txt', generated_dataset)
            report(generated_dataset, fpath)


# Third method note
"""
Combination of Jon Yearsley (2021). Generate AR1 spatial data (https://www.mathworks.com/matlabcentral/fileexchange/5099-generate-ar1-spatial-data), MATLAB Central File Exchange. Retrieved November 30, 2021.
and Dakos et al. 10.1073/pnas.0802430105
"""

# Generate datasets method 3
"""
Generates dataset(s) with similar autocorrelation, mean and variance by generating datasets with an AR(1) model trained
on the original dataset.

Args:
-----

data : numpy array, the spatial datasets.

realizations : int, the number of datasets generated.

method : str, either 'Normal' or 'Adjusted'. For 'Normal', the standard AR(1) format is used. For 'Adjusted', the AR(1)
    format of Jon Yearsley (2021) is used.

path : str, the filepath where the original dataset can be found.

variable : name of the variable.

stdev_error : int/float, the standard deviation of the white noise process.

"""

def method3_(dataset, realizations=1, method='Normal', path='./1/', variable='xxx', stdev_error=1.0):
    generated_number_length = ews.generated_number_length(realizations)

    steps = np.arange(cfg.interval_map_snapshots, cfg.number_of_timesteps_weekly + cfg.interval_map_snapshots,
                      cfg.interval_map_snapshots)

    spatial_correlation = ews.spatial_corr(dataset)

    for k, data in enumerate(dataset):
        Morans_I = spatial_correlation[k]
        alpha0_1 = np.nanmean(data) * (1 - Morans_I)
        alpha0_2 = np.nanmean(data)
        sig2 = np.nanvar(data) * (1 - Morans_I ** 2)

        dim = data.shape
        N = np.prod(dim)

        W = np.zeros((N, N))
        np.fill_diagonal(W[1:], Morans_I / 4)
        np.fill_diagonal(W[:, 1:], Morans_I / 4)

        M = np.identity(N) - W
        inv_M = np.linalg.inv(M)

        for realization in range(realizations):
            random_error = np.random.normal(loc=0.0, scale=stdev_error, size=N)

            if method == 'Adjusted':
                if np.isnan(np.sqrt(sig2)):
                    generated_dataset_numpy = np.dot(inv_M, random_error * 0.0).reshape(dim) + alpha0_1
                else:
                    generated_dataset_numpy = np.dot(inv_M, random_error * np.sqrt(sig2)).reshape(dim) + alpha0_1
            elif method == 'Normal':
                generated_dataset_numpy = np.dot(inv_M, random_error).reshape(dim) + alpha0_2

            generated_dataset = numpy2pcr(Scalar, generated_dataset_numpy, np.NaN)
            generated_number_string = 'm3g' + str(realization).zfill(generated_number_length)
            dir_name = os.path.join(path + generated_number_string)

            if os.path.isdir(dir_name) == False:
                os.makedirs(dir_name)

            fname = ews.file_name_str(variable.name, steps[k])
            fpath = os.path.join(dir_name, fname)
            # np.savetxt(fpath + '.numpy.txt', generated_dataset)
            report(generated_dataset, fpath)
