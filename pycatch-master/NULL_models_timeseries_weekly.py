import numpy as np
from scipy import fft
import statsmodels.api
import os
from scipy import ndimage
import EWSPy as ews
import EWS_configuration as cfg


### Null models timeseries (Dakos et al. 2008) ###

# TODO - method 2 did not return the right mean - check solution -, check method 3

# Detrend dataset
"""
Detrends the given dataset making use of either Gaussian filtering or linear detrending, as specified in the 
configuration. Optionally, this data is saved.

Args:
-----

data : numpy array, the timeseries data.

realizations : int, the number of datasets generated. Used for folder name in the case of detrending.

path : str, the filepath where the original dataset can be found.

file_name : str, name of the variable.

Return:
-----

detrended_data : The detrended timeseries data.

"""


def detrend_(data, realizations=1, path='./1/', file_name='xxx'):
    detrended_data = data

    if cfg.detrended_method == 'Gaussian':
        gaussian_filter = ndimage.gaussian_filter1d(data, cfg.detrended_sigma)
        detrended_data -= gaussian_filter
    elif cfg.detrended_method == 'Linear':
        mean = np.nanmean(data)
        detrended_data -= mean
    elif cfg.detrended_method is not 'None':
        print("Invalid input for detrending_temp in generate_datasets (EWS_weekly.py). No detrending done.")

    if cfg.save_detrended_data:
        generated_number_length = ews.generated_number_length(realizations)
        generated_number_string = 'dtr' + str(0).zfill(generated_number_length)
        dir_name = os.path.join(path + generated_number_string)

        if os.path.isdir(dir_name) == False:
            os.makedirs(dir_name)

        fname1 = ews.file_name_str(file_name, cfg.number_of_timesteps_weekly)
        fpath1 = os.path.join(dir_name, fname1)
        np.savetxt(fpath1 + '.numpy.txt', detrended_data)

        if cfg.detrended_method == 'Gaussian':
            fname2 = ews.file_name_str(file_name + 'g', cfg.number_of_timesteps_weekly)
            fpath2 = os.path.join(dir_name, fname2)
            np.savetxt(fpath2 + '.numpy.txt', gaussian_filter)

    return detrended_data


# Generate datasets method 1
"""
Generates dataset(s) with similar mean and variance by randomly picking values from the original dataset. In the case
of replace==False, this is similar to shuffling the dataset.

Args:
-----

data : numpy array, the spatial datasets.

realizations : int, the number of datasets generated.

path : str, the filepath where the original dataset can be found.

file_name : str, name of the variable.

replace : bool, selects whether new values are picked from the original dataset or the original dataset minus previously
    picked values. Usually set to False to ensure similar mean and variance for smaller datasets.

"""


def method1_(data, realizations=1, path='./1/', file_name='xxx', replace=False):
    generated_number_length = ews.generated_number_length(realizations)

    for realization in range(realizations):
        generated_dataset = np.random.choice(data, len(data), replace=replace)

        generated_number_string = 'm1g' + str(realization).zfill(generated_number_length)
        dir_name = os.path.join(path + generated_number_string)

        if os.path.isdir(dir_name) == False:
            os.makedirs(dir_name)

        fname = ews.file_name_str(file_name, cfg.number_of_timesteps_weekly)
        fpath = os.path.join(dir_name, fname)
        np.savetxt(fpath + '.numpy.txt', generated_dataset)


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

file_name : str, name of the variable.

replace : bool, selects whether new values are picked from the original dataset or the original dataset minus previously
    picked values.

"""


def method2_(data, realizations=1, method='None', path='./1/', file_name='xxx', replace=False):
    generated_number_length = ews.generated_number_length(realizations)
    if method == 'Detrending':
        y = data
        x = np.arange(len(data))
        coef = np.polyfit(x, y, 1)
        poly1d_fn = np.poly1d(coef)  # Function which takes in x and returns an estimate for y
        lin_detr = poly1d_fn(x)
        data = data - lin_detr  # Remove trend from data

    fft_ = fft.fft(data)
    fft_mag = np.abs(fft_)
    fft_phases = np.angle(fft_)

    for realization in range(realizations):
        fft_phases_new = fft_phases.copy()
        if data.size % 2 == 0:
            i = int(len(fft_phases_new) / 2)
            fft_phases_left_half = fft_phases[1:i]
            fft_shuffled_phases_lh = np.random.choice(fft_phases_left_half, len(fft_phases_left_half), replace=replace)
            fft_shuffled_phases_rh = - fft_shuffled_phases_lh[::-1]
            fft_phases_new = np.concatenate((np.array((fft_[0],)), fft_shuffled_phases_lh, np.array((fft_phases[i],)),
                                             fft_shuffled_phases_rh))
        else:
            i = int(len(fft_phases_new) / 2 + 1)
            fft_phases_left_half = fft_phases[1:i]
            fft_shuffled_phases_lh = np.random.choice(fft_phases_left_half, len(fft_phases_left_half), replace=replace)
            fft_shuffled_phases_rh = - fft_shuffled_phases_lh[::-1]
            fft_phases_new = np.concatenate((np.array((fft_[0],)), fft_shuffled_phases_lh,
                                             fft_shuffled_phases_rh))

        fft_sym = fft_mag * (np.cos(fft_phases_new) + 1j * np.sin(fft_phases_new))
        generated_dataset = fft.ifft(fft_sym)

        generated_dataset = np.absolute(generated_dataset)  # TODO - Check if this is correct
        if method == 'Detrending':
            generated_dataset = generated_dataset + lin_detr

        generated_number_string = 'm2g' + str(realization).zfill(generated_number_length)
        dir_name = os.path.join(path + generated_number_string)

        if os.path.isdir(dir_name) == False:
            os.makedirs(dir_name)

        fname = ews.file_name_str(file_name, cfg.number_of_timesteps_weekly)
        fpath = os.path.join(dir_name, fname)
        np.savetxt(fpath + '.numpy.txt', generated_dataset)


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

file_name : str, name of the variable.

stdev_error : int/float, the standard deviation of the white noise process.

"""


def method3_(data, realizations=1, method='Normal', path='./1/', file_name='xxx', stdev_error=1):
    generated_number_length = ews.generated_number_length(realizations)

    alpha1 = statsmodels.api.tsa.acf(data, nlags=1)
    sig2 = np.nanvar(data) * (1 - alpha1[1] ** 2)
    alpha0_1 = np.nanmean(data) * (1 - alpha1[1])
    alpha0_2 = np.nanmean(data)

    for realization in range(realizations):
        e = np.random.normal(loc=0.0, scale=stdev_error, size=len(data))
        AR1m = np.ones(len(data)) * alpha0_2
        if method == 'Adjusted':
            for i in range(len(data)):
                AR1m[i] = alpha1[1] * AR1m[i - 1] + alpha0_1 + np.sqrt(sig2) * e[i]
        elif method == 'Normal':
            for i in range(len(data)):
                AR1m[i] = alpha1[1] * AR1m[i - 1] + alpha0_2 + e[i]
        generated_dataset = AR1m

        generated_number_string = 'm3g' + str(realization).zfill(generated_number_length)
        dir_name = os.path.join(path + generated_number_string)

        if os.path.isdir(dir_name) == False:
            os.makedirs(dir_name)

        fname = ews.file_name_str(file_name, cfg.number_of_timesteps_weekly)
        fpath = os.path.join(dir_name, fname)
        np.savetxt(fpath + '.numpy.txt', generated_dataset)
