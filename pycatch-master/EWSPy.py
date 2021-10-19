import pcraster as pcr
from pcraster import pcr2numpy
from pcraster._pcraster import readmap
import pcraster.framework as pcrfw
import scipy.stats
from scipy import sparse
from scipy.spatial.distance import pdist

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
from scipy.signal import convolve2d
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

#####################################
### Spatial early-warning signals ###
#####################################
"""
--------------------------------
Methods/indicators per phenomena
--------------------------------

Rising memory:
    Spatial correlation: Done & tested
    Return time*: TODO - Only applied in a single paper (Wissel, 1984), only mentioned by name in other works.
                         Furthermore, this return time (to equilibrium) is a function of time --> temporal
    Discrete Fourier Transform (DFT): TODO
Rising variability:
    Spatial variance: Done & tested
    Spatial skewness: Done & tested
Patchiness*:
    Spatial variance & skewness: Done & tested (as above)
    Patch-size distributions: NOT INCLUDED - Only relevant for biomass
    Regular spotted patterns: NOT INCLUDED - Only relevant for biomass
    Power spectrum: TODO
    
*= not sure if necessary for this study
--------------------------------
"""
#####################################

def time_series2snapshots(numpy_matrix, interval):
    return numpy_matrix[::interval]
    # return np.array([numpy_matrix[i] for i in range(len(numpy_matrix)) if i % window_size == 0])

def spatial_mean(numpy_matrix):
    # return np.array([np.nanmean(array) for array in numpy_matrix])
    return np.nanmean(numpy_matrix, axis=(1,2))

def spatial_std(numpy_matrix):
    return np.nanstd(numpy_matrix, axis=(1,2))

def spatial_var(numpy_matrix):
    return np.nanvar(numpy_matrix, axis=(1,2))

def spatial_skw(numpy_matrix):
    return scipy.stats.skew(np.nditer(numpy_matrix), axis=(1,2), nan_policy='omit')

def spatial_krt(numpy_matrix):
    return scipy.stats.kurtosis(np.nditer(numpy_matrix), axis=(1,2), nan_policy='omit')

# Rook-neighborhood for spatial correlation analogue to lag-1 autocorrelation in time
rook_neighborhood = np.array([
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 0]
])

def spatial_corr(numpy_matrix):
    mean = np.nanmean(numpy_matrix)
    var = np.nanvar(numpy_matrix)
    numpy_matrix -= mean

    is_nan = np.isnan(numpy_matrix) # missing values in map are assumed to be np.NaN
    is_not_nan = ~ is_nan
    is_not_nan_as_nr = is_not_nan.astype(float)

    numpy_matrix_copy = numpy_matrix.copy()
    numpy_matrix_copy[is_nan] = 0
    sum_neighbours = convolve2d(numpy_matrix_copy, rook_neighborhood, mode='same')
    n_neighbours_times_avg = convolve2d(is_not_nan_as_nr, rook_neighborhood * mean, mode='same')
    n_neighbours_times_avg[is_nan] = 0
    P1 = sum(sum(numpy_matrix * (sum_neighbours - n_neighbours_times_avg)))

    return P1 / (4 * var * np.count_nonzero(is_not_nan_as_nr))

def spatial_corr_inside(numpy_matrix): # not fit for numpy matrices that contain np.NaN, use spatial_corr instead!
    mean = spatial_mean(numpy_matrix)
    var = spatial_var(numpy_matrix)
    storage = 0.0
    for m in range(1, numpy_matrix.shape[0]-1):
        for n in range(1, numpy_matrix.shape[1]-1):
            storage += (numpy_matrix[m, n] - mean) * ((numpy_matrix[m, n-1] + numpy_matrix[m, n+1] + numpy_matrix[m-1, n]
                                                      + numpy_matrix[m+1, n]) - (4 * mean))
    spatial_corr = (storage) / (4 * var * ((numpy_matrix.shape[0]-2) * (numpy_matrix.shape[1]-2)))
    return spatial_corr

# def spatial_corr_package(numpy_matrix):
#     weights = lat2W(nrows=numpy_matrix.shape[0], ncols=numpy_matrix.shape[1], rook=True)
#     mi = Moran_Local(numpy_matrix, weights)
#     return sum(mi.Is)

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
#
# data = (np.array(range(1,101)).reshape(10,10))
# print(data)
# print(spatial_spec(data))

# psd2D = spatial_spec(data)
# print(psd2D)

# py.figure(1)
# py.clf()
# py.imshow(np.log10(psd2D), cmap=py.cm.jet)
# py.show()

#########################################
### Time series early-warning signals ###
#########################################
"""
--------------------------------
Methods/indicators per phenomena
--------------------------------
Rising memory:
    Autocorrelation lag 1: Done & tested
    Autoregressive coefficient of AR(1) model* ***: NOT INCLUDED
    Return rate (inverse of AR(1) coefficient)* ***: NOT INCLUDED
    Detrended fluctuation analysis indicator: TODO
    Spectral density: TODO
    Spectral ratio (of low to high frequencies): TODO
    Spectral exponent: TODO
Rising variability & flickering:
    Standard deviation/variance: Done & tested
    Coefficient of variation: Done & tested
    Skewness: Done & tested
    Kurtosis: Done & tested
    Conditional heteroskedasticity: TODO
    BDS test**: TODO

*Models are not included, metrics are.
**= Can help to avoid false detections due to model misspecification.
***= alternative ways to measure autocorrelation lag-1
--------------------------------
"""
#########################################

def time_series2time_windows(time_series, window_size=10):
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

###

def autocovariance(numpy_array, lag=1):
    number = len(numpy_array)
    mean = sum(numpy_array) / number
    autocovariance = sum([(numpy_array[i] - mean) * (numpy_array[i+lag] - mean) for i in range(number - lag)]) / number
    return autocovariance

np.seterr(invalid='ignore') # would like to remove this ofc
def temporal_autocorrelation(numpy_array, lag=1):
    return np.true_divide(autocovariance(numpy_array, lag=lag), temporal_var(numpy_array))

def temporal_spectrum(numpy_array):
    freqs, times, spectrogram = signal.spectrogram(numpy_array)
    return spectrogram

def temporal_PSD(numpy_array):
    freqs, psd = signal.welch(numpy_array)
    return freqs, psd

def temporal_mean(numpy_array):
    return np.nanmean(numpy_array, axis=1) # return np.array([np.nanmean(array) for array in numpy_array])

def temporal_std(numpy_array):
    return np.nanstd(numpy_array, axis=1)

def temporal_var(numpy_array):
    return np.nanvar(numpy_array, axis=1)

np.seterr(invalid='ignore') # would like to remove this ofc
def temporal_cv(numpy_array):
    return np.true_divide(temporal_std(numpy_array), temporal_mean(numpy_array))

def temporal_skw(numpy_array):
    return scipy.stats.skew(numpy_array, axis=1, nan_policy='omit')

def temporal_krt(numpy_array):
    return scipy.stats.kurtosis(numpy_array, axis=1, nan_policy='omit')

#########################################

# np.random.seed(42)
# # data = np.random.rand(100, 100) * 100
# # time_series_stack = np.random.rand(10, 1000) * 100
# data = np.random.normal(10, 5, (100, 100))
# #data = np.array(range(1,101)).reshape(10,10)
# time_series_stack = np.random.normal(10, 5, (10, 1000))
#
# import time
#
# ### Spatial tests for single map ###
# print("\n---=<#>=--- start of spatial tests ---=<#>=---")
# start_time = time.time()
#
# print("spatial mean:", spatial_mean(data))
# print("spatial std:", spatial_std(data))
# print("spatial var:", spatial_var(data))
# print("spatial skw:", spatial_skw(data))
# print("spatial krt:", spatial_krt(data))
# print("spatial corr (inside):", spatial_corr_inside(data))
# print("spatial corr (full):", spatial_corr(data))
# #print("rspec:", spatial_spec(data))
#
# finished_time = time.time() - start_time
# print(f"-=- runtime is {finished_time} seconds -=-")
#
# ### Temporal tests for 10 time series ###
# print("\n---=<#>=--- start of temporal tests ---=<#>=---")
# start_time = time.time()
#
# tstorage = [0.0] * time_series_stack.shape[0]
#
# for k, time_series in enumerate(time_series_stack):
#     tstorage[k] += temporal_mean(time_series)
# print("temporal mean:", tstorage)
#
# for k, time_series in enumerate(time_series_stack):
#     tstorage[k] += temporal_autocorrelation(time_series)
# print("temporal autocorrelation lag-1:", tstorage)
#
# for k, time_series in enumerate(time_series_stack):
#     tstorage[k] += temporal_std(time_series)
# print("temporal std:", tstorage)
#
# for k, time_series in enumerate(time_series_stack):
#     tstorage[k] += temporal_skw(time_series)
# print("temporal skw:", tstorage)
#
# for k, time_series in enumerate(time_series_stack):
#     tstorage[k] += temporal_krt(time_series)
# print("temporal krt:", tstorage)
#
# for k, time_series in enumerate(time_series_stack):
#     tstorage[k] += temporal_cv(time_series)
# print("temporal CV:", tstorage)
#
# finished_time = time.time() - start_time
# print(f"-=- runtime is {finished_time} seconds -=-")

#####################################
