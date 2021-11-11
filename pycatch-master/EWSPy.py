import pcraster as pcr
from pcraster import pcr2numpy
from pcraster._pcraster import readmap
import pcraster.framework as pcrfw
import scipy.stats
from scipy import sparse
from scipy.spatial.distance import pdist

from statsmodels.tsa.ar_model import AutoReg
from statsmodels.stats.diagnostic import het_arch

import numpy as np
import skgstat
import gstools
import os
import operator
import glob
import subprocess
import math
import rpy2
from libpysal.weights import Queen, lat2W, KNN, lat2SW
from esda.moran import Moran, Moran_Local
from scipy.signal import convolve, convolve2d
from scipy import signal
from scipy.spatial.distance import pdist, squareform
from scipy import fftpack
from scipy import fft
import pylab as py
import matplotlib.pyplot as plt

import sys
sys.path.append("./pcrasterModules/")

from collections import deque
# from PCRaster.NumPy import *
import random

################################################
### Class StateVariable for Variable objects ###
################################################

class StateVariable:
    def __init__(self, name, meanmax='mean', window_size=100, snapshot_interval=100, spatial=True, temporal=True, datatype='map'):
        self.name = name
        self.meanmax = meanmax
        self.window_size = window_size
        self.snapshot_interval = snapshot_interval
        self.spatial = spatial
        self.temporal = temporal
        self.datatype = datatype

#####################################
### Spatial early-warning signals ###
#####################################
"""
--------------------------------
Methods/indicators per phenomena
--------------------------------

Rising memory:
    Spatial correlation: Done & tested
    Return time*: NOT INCLUDED - Only applied in a single paper (Wissel, 1984), only mentioned by name in other works.
                                Furthermore, this return time (to equilibrium) is a function of time --> temporal
                                Or, it should be seen as return distance (spatial distance to mean state, bit vague)
    Discrete Fourier Transform (DFT): Done - Needs testing & conformation of method! -> how2plot?
Rising variability:
    Spatial variance: Done & tested
    Spatial skewness: Done & tested
Patchiness*:
    Spatial variance & skewness: Done & tested (as above)
    Patch-size distributions: NOT INCLUDED - Only relevant for biomass
    Regular spotted patterns: NOT INCLUDED - Only relevant for biomass
    Power spectrum: Done - Needs testing & conformation of method! -> only works for square matrices
    
*= not sure if necessary for this study
--------------------------------
"""
#####################################

def time_series2snapshots(numpy_matrix, interval):
    # return numpy_matrix[::interval]
    return np.array([numpy_matrix[i] for i in range(len(numpy_matrix)) if i % interval == 0])

def spatial_mean(numpy_matrix):
    # return np.array([np.nanmean(array) for array in numpy_matrix])
    return np.nanmean(numpy_matrix, axis=(1,2))

def spatial_std(numpy_matrix):
    return np.nanstd(numpy_matrix, axis=(1,2))

def spatial_var(numpy_matrix):
    return np.nanvar(numpy_matrix, axis=(1,2))

def spatial_skw(numpy_matrix):
    return [scipy.stats.skew(np.nditer(array), nan_policy='omit') for array in numpy_matrix]

def spatial_krt(numpy_matrix):
    return [scipy.stats.kurtosis(np.nditer(array), nan_policy='omit') for array in numpy_matrix]

# Rook-neighborhood for spatial correlation analogue to lag-1 autocorrelation in time
rook_neighborhood = np.array([
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 0]
])

queen_neighborhood = np.array([
    [1, 1, 1],
    [1, 0, 1],
    [1, 1, 1]
])

rook_neighborhood_2 = np.array([
    [0, 0, 1, 0, 0],
    [0, 1, 1, 1, 0],
    [1, 1, 0, 1, 1],
    [0, 1, 1, 1, 0],
    [0, 0, 1, 0, 0]
])

queen_neighborhood_2 = np.array([
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 0, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1]
])

def spatial_corr(numpy_matrix): # Moran's I
    mean = spatial_mean(numpy_matrix)

    var = spatial_var(numpy_matrix)

    numpy_matrix_mmean = np.copy(numpy_matrix)
    numpy_matrix_mmean -= mean[:, None, None]

    is_nan = np.isnan(numpy_matrix_mmean) # missing values in map are assumed to be np.NaN
    is_not_nan = ~ is_nan
    is_not_nan_as_nr = is_not_nan.astype(float)

    numpy_matrix_copy = np.copy(numpy_matrix)
    numpy_matrix_copy[is_nan] = 0

    sum_neighbours = convolve(numpy_matrix_copy, rook_neighborhood[None, :, :], mode='same')
    # sum_neighbours = np.array([convolve2d(array, rook_neighborhood, mode='same') for array in numpy_matrix_copy])
    # sum_neighbours = convolve3D???

    n_neighbours_times_avg = convolve(is_not_nan_as_nr * mean[:, None, None], rook_neighborhood[None, :, :], mode='same')
    # n_neighbours_times_avg = np.array([convolve2d(is_not_nan_as_nr[i], rook_neighborhood * mean[i], mode='same') for i in range(len(is_not_nan_as_nr))])
    n_neighbours_times_avg[is_nan] = 0

    P1 = np.nansum(numpy_matrix_mmean * (sum_neighbours - n_neighbours_times_avg), axis=(1,2))
    P2 = (4 * var * is_not_nan_as_nr.sum(axis=(1, 2)))
    return P1 / P2

# def spatial_corr_2D_inside(numpy_matrix): # not fit for numpy matrices that contain np.NaN, use spatial_corr instead!
#     mean = np.nanmean(numpy_matrix)
#     var = np.nanvar(numpy_matrix)
#     storage = 0.0
#     for m in range(1, numpy_matrix.shape[0]-1):
#         for n in range(1, numpy_matrix.shape[1]-1):
#             storage += (numpy_matrix[m, n] - mean) * ((numpy_matrix[m, n-1] + numpy_matrix[m, n+1] + numpy_matrix[m-1, n]
#                                                       + numpy_matrix[m+1, n]) - (4 * mean))
#     spatial_corr = (storage) / (4 * var * ((numpy_matrix.shape[0]-2) * (numpy_matrix.shape[1]-2)))
#     return spatial_corr
#
# def spatial_corr_2D(numpy_matrix):
#     mean = np.nanmean(numpy_matrix)
#
#     var = np.nanvar(numpy_matrix)
#
#     numpy_matrix_mmean = np.copy(numpy_matrix)
#     numpy_matrix_mmean -= mean
#
#     is_nan = np.isnan(numpy_matrix_mmean) # missing values in map are assumed to be np.NaN
#     is_not_nan = ~ is_nan
#     is_not_nan_as_nr = is_not_nan.astype(float)
#
#     numpy_matrix_copy = numpy_matrix.copy()
#     numpy_matrix_copy[is_nan] = 0
#
#     sum_neighbours = convolve2d(numpy_matrix_copy, rook_neighborhood, mode='same')
#
#     n_neighbours_times_avg = convolve2d(is_not_nan_as_nr, rook_neighborhood * mean, mode='same')
#     n_neighbours_times_avg[is_nan] = 0
#
#     P1 = np.nansum(np.nansum(numpy_matrix_mmean * (sum_neighbours - n_neighbours_times_avg)))
#     P2 = (4 * var * np.count_nonzero(is_not_nan_as_nr))
#     return P1 / P2

def spatial_DFT(numpy_matrix):
    return fft.fft2(numpy_matrix, axes=(-2,))

def spatial_power_spec(numpy_matrix): # Only works for square matrices! Power spectrum as function of wave number (P(k))
    n = numpy_matrix.shape[0]

    fourier_image = fft.fft2(numpy_matrix) # fourier 'image'
    fourier_amplitudes = np.abs(fourier_image) ** 2 # fourier amplitudes

    kfreq = fft.fftfreq(n) * n # 1D array containing the wave vectors in k space
    kfreq2D = np.meshgrid(kfreq, kfreq) # convertion to 2D array matching the layout of the fourier 'image'
    knorm = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2) # norm of wave vectors
    knorm = knorm.flatten()
    fourier_amplitudes = fourier_amplitudes.flatten()

    kbins = np.arange(0.5, n//2+1, 1.) # start & end points of all bins
    kvals = 0.5 * (kbins[1:] + kbins[:-1]) # corresponding k values

    Abins, _, _ = scipy.stats.binned_statistic(knorm, fourier_amplitudes, statistic = "mean", bins = kbins) # average Fourier amplitude (**2) in each bin
    Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2) # total power --> multiply by area in each bin (eq5)

    return fourier_image, kvals, Abins

# def spatial_spec(numpy_matrix):
    # nr = numpy_matrix.shape[0]
    # nc = numpy_matrix.shape[1]
    #
    # n0x = np.floor(nc/2)+1
    # n0y = np.floor(nr/2)+1
    #
    # f1 = np.tile(range(1,nc+1) - n0x, (nc, 1)).T
    # f2 = np.tile(range(1,nr+1) - n0y, (nr, 1))
    #
    # DIST = np.sqrt(f1**2 + f2**2)
    #
    # ANGLE = np.arctan2(-f2, f1) * 180 / np.pi

#     fourier = np.fft.fft2(numpy_matrix)
#     fourier_shift = np.fft.fftshift(fourier, axes=(1,))
#     # F1 = fftpack.fftshift( fftpack.fft2(numpy_matrix) )
#     F2 = fft.fftshift(fft.fft2(numpy_matrix), axes=(1,))
#     psd2D = np.abs(fourier_shift) ** 2
#
#     return fourier_shift, F2, fourier_shift - F2

# data = (np.array(range(1,101)).reshape(10,10))
# print(data)
# print(spatial_power_spec(data))

#########################################
### Time series early-warning signals ###
#########################################
"""
--------------------------------
Methods/indicators per phenomena
--------------------------------
Rising memory:
    Autocorrelation lag 1: Done & tested
    Autoregressive coefficient of AR(1) model: Done - Needs testing & further improvements
    Return rate (inverse of AR(1) coefficient): Done - Needs testing & further improvements
    Detrended fluctuation analysis indicator: Done - Needs testing & further improvements
    Spectral density: TODO
    Spectral ratio (of low to high frequencies): TODO
    Spectral exponent: TODO
Rising variability & flickering:
    Standard deviation/variance: Done & tested
    Coefficient of variation: Done & tested
    Skewness: Done & tested
    Kurtosis: Done & tested
    Conditional heteroskedasticity: Done - Needs testing & further improvements
    BDS test**: TODO

*Models are not included, metrics are.
**= Can help to avoid false detections due to model misspecification.
***= alternative ways to measure autocorrelation lag-1
--------------------------------
"""
#########################################

def time_series2time_windows(time_series, window_size=10):
    # return time_series[::window_size]
    return np.array([time_series[i:i + window_size] for i in range(0, len(time_series), window_size)])

def mean_time_series(stack_of_maps_as_list): # needs testing, improvements
    mean_time_series = [0.0] * len(stack_of_maps_as_list)
    for k, map in enumerate(stack_of_maps_as_list):
        mean_time_series[k] += np.nanmean(map)
    return mean_time_series

def max_time_series(stack_of_maps_as_list): # needs testing, improvements
    max_time_series = [0.0] * len(stack_of_maps_as_list)
    for k, map in enumerate(stack_of_maps_as_list):
        max_time_series[k] += np.nanmax(map)
    return max_time_series

#########################################

def temporal_AR1(stack_of_windows):
    mean = temporal_mean(stack_of_windows)

    stack_of_windows_mmean = np.copy(stack_of_windows)
    stack_of_windows_mmean -= mean[:, None]

    AR1_params = []
    for numpy_array in stack_of_windows_mmean:
        mod = AutoReg(numpy_array, 1, trend='n').fit()
        AR1_params = np.append(AR1_params, mod.params)
    return AR1_params

def temporal_returnrate(stack_of_windows):
    # returns inf when division by zero; NaN is used for division by zero
    return np.reciprocal(temporal_AR1(stack_of_windows))

def temporal_cond_het(stack_of_windows, n_lags=4, ddof=1):
    mean = temporal_mean(stack_of_windows)

    stack_of_windows_mmean = np.copy(stack_of_windows)
    stack_of_windows_mmean -= mean[:, None]

    cond_het = [0.0] * len(stack_of_windows)
    for k, numpy_array in enumerate(stack_of_windows_mmean):
        residuals = AutoReg(numpy_array, 1, trend='n').fit().resid
        cond_het[k] = np.array(het_arch(residuals, nlags=n_lags, ddof=ddof))
    return cond_het

# def autocovariance(numpy_array, lag=1):
#     number = len(numpy_array)
#     mean = sum(numpy_array) / number
#     autocovariance = sum([(numpy_array[i] - mean) * (numpy_array[i+lag] - mean) for i in range(number - lag)]) / number
#     return autocovariance

#np.seterr(invalid='ignore') # would like to do without this
def temporal_autocorrelation(numpy_array, lag=1):
    return np.true_divide(temporal_autocovariance(numpy_array, lag=lag), temporal_var(numpy_array))

def temporal_autocovariance(numpy_array, lag=1):
    auto_cov = [0.0] * len(numpy_array)
    for k, window in enumerate(numpy_array):
        N = len(window)
        mean = np.nanmean(window)

        end_padded_series = np.zeros(N+lag)
        end_padded_series[:N] = window - mean
        start_padded_series = np.zeros(N+lag)
        start_padded_series[lag:] = window - mean

        auto_cov[k] = np.sum(start_padded_series * end_padded_series) / N
    return auto_cov

# def temporal_spectrum(numpy_array):
#     freqs, times, spectrogram = signal.spectrogram(numpy_array)
#     return spectrogram
#
# def temporal_PSD(numpy_array):
#     # freqs, psd = signal.welch(numpy_array)
#     # freqs, psd = signal.periodogram(numpy_array)
#     # return freqs, psd
#     return None

def temporal_mean(numpy_array):
    return np.nanmean(numpy_array, axis=1) # return np.array([np.nanmean(array) for array in numpy_array])

def temporal_std(numpy_array):
    return np.nanstd(numpy_array, axis=1)

def temporal_var(numpy_array):
    return np.nanvar(numpy_array, axis=1)

np.seterr(invalid='ignore') # would like to do without this
def temporal_cv(numpy_array):
    return np.true_divide(temporal_std(numpy_array), temporal_mean(numpy_array))

def temporal_skw(numpy_array):
    return scipy.stats.skew(numpy_array, axis=1, nan_policy='omit')

def temporal_krt(numpy_array):
    return scipy.stats.kurtosis(numpy_array, axis=1, nan_policy='omit')

def calc_rms(numpy_array, scale): # windowed Root Mean Square with linear detrending
    # Making of an array with data divided into segments
    shape = (numpy_array.shape[0] // scale, scale)
    segments = np.array([numpy_array[i:i + scale] for i in range(0, len(numpy_array), scale)]).reshape(shape)

    # x = np.copy(numpy_array)
    # segments = np.lib.stride_tricks.as_strided(x, shape=shape) # - 'Dangerous' function; excluded

    # Vector of x-axis
    scale_ax = np.arange(scale)
    rms = np.zeros(segments.shape[0])
    for i, segment in enumerate(segments):
        coeff = np.polyfit(scale_ax, segment, 1)
        xfit = np.polyval(coeff, scale_ax)
        # Detrending & computing RMS of each window
        rms[i] = np.sqrt(np.nanmean((segment - xfit)**2))
    return rms

def temporal_dfa(stack_of_windows, scales=np.array([10])):
    # TODO - Works for a single time_window --> needs to be working *nicely* for array of time_windows
    fluct = []
    coeff = []

    for numpy_array in stack_of_windows:
        # Cumulative sum of a single window with subtracted offset
        y = np.nancumsum(numpy_array - np.nanmean(numpy_array)) # noise like time series to random walk time series

        # RMS for each segment
        fluctuations = np.zeros(len(scales)) # a.ii.1
        for i, sc in enumerate(scales):
            fluctuations[i] = np.sqrt(np.mean(calc_rms(y, sc)**2))

        coefficients = np.polyfit(np.log2(scales), np.log2(fluctuations), 1)

        fluct = np.append(fluct, fluctuations)
        coeff = np.append(coeff, coefficients[0])

    fluct = np.array_split(fluct, len(scales))

    return scales, fluct, coeff

#########################################
