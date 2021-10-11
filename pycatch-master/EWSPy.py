import pcraster as pcr
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

import sys
sys.path.append("./pcrasterModules/")

##
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()  # nodig bij nieuwere rpy2 versie
import rpy2.robjects as robjects
##

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

def spatial_mean(numpy_matrix):
    return np.nanmean(numpy_matrix)

def spatial_std(numpy_matrix):
    return np.nanstd(numpy_matrix)

def spatial_var(numpy_matrix):
    return np.nanvar(numpy_matrix)

def spatial_skw(numpy_matrix):
    return scipy.stats.skew(np.nditer(numpy_matrix), nan_policy='omit')


def spatial_krt(numpy_matrix):
    return scipy.stats.kurtosis(np.nditer(numpy_matrix), nan_policy='omit')

# Rook-neighborhood for spatial correlation analogue to lag-1 autocorrelation in time
rook_neighborhood = np.array([
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 0]
])

def spatial_corr(numpy_matrix):
    # P1 = (numpy_matrix[m, n] - mean) * (numpy_matrix[m, n-1] + numpy_matrix[m, n+1] + numpy_matrix[m-1, n]
    #                                     + numpy_matrix[m+1, n]) - (4 * mean)
    mean = spatial_mean(numpy_matrix)
    var = spatial_var(numpy_matrix)
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

    spatial_corr = P1 / (4 * var * np.count_nonzero(is_not_nan_as_nr))
    return spatial_corr

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

#########################################
### Time series early-warning signals ###
#########################################
"""
--------------------------------
Methods/indicators per phenomena
--------------------------------
Rising memory:
    Autocorrelation lag 1: Done & tested
    Autoregressive coefficient of AR(1) model: TODO
    Return rate (inverse of AR(1) coefficient): TODO
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
--------------------------------
"""
#########################################

def stack_of_maps2time_series(stack_of_maps_as_list): # needs testing, improvements
    time_series = [0.0] * len(stack_of_maps_as_list)
    for k, map in enumerate(stack_of_maps_as_list):
        time_series[k] += np.nanmean(pcr2numpy(map, np.NaN))
    return time_series

def autocovariance(numpy_array, lag=1):
    number = len(numpy_array)
    mean = sum(numpy_array) / number
    autocovariance = sum([(numpy_array[i] - mean) * (numpy_array[i+lag] - mean) for i in range(number - lag)]) / number
    return autocovariance

def temporal_autocorrelation(numpy_array, lag=1):
    return autocovariance(numpy_array, lag=lag) / temporal_var(numpy_array)

def temporal_mean(numpy_array):
    return np.nanmean(numpy_array)

def temporal_std(numpy_array):
    return np.nanstd(numpy_array)

def temporal_var(numpy_array):
    return np.nanvar(numpy_array)

def temporal_cv(numpy_array):
    return temporal_std(numpy_array) / temporal_mean(numpy_array)

def temporal_skw(numpy_array):
    return scipy.stats.skew(np.nditer(numpy_array), nan_policy='omit')

def temporal_krt(numpy_array):
    return scipy.stats.kurtosis(np.nditer(numpy_array), nan_policy='omit')

def conditional_heteroskedacity(numpy_array):
    return None

#########################################

np.random.seed(42)
data = np.random.rand(100, 100) * 100
time_series_stack = np.random.rand(10, 1000) * 100

import time

### Spatial tests for single map ###
print("\n---=<#>=--- start of spatial tests ---=<#>=---")
start_time = time.time()

print("spatial mean:", spatial_mean(data))
print("spatial std:", spatial_std(data))
print("spatial var:", spatial_var(data))
print("spatial skw:", spatial_skw(data))
print("spatial krt:", spatial_krt(data))
print("spatial corr (inside):", spatial_corr_inside(data))
print("spatial corr (full):", spatial_corr(data))

finished_time = time.time() - start_time
print(f"-=- runtime is {finished_time} seconds -=-")

### Temporal tests for 10 time series ###
print("\n---=<#>=--- start of temporal tests ---=<#>=---")
start_time = time.time()

tstorage = [0.0] * time_series_stack.shape[0]

for k, time_series in enumerate(time_series_stack):
    tstorage[k] += temporal_mean(time_series)
print("temporal mean:", tstorage)

for k, time_series in enumerate(time_series_stack):
    tstorage[k] += temporal_autocorrelation(time_series)
print("temporal autocorrelation lag-1:", tstorage)

for k, time_series in enumerate(time_series_stack):
    tstorage[k] += temporal_std(time_series)
print("temporal std:", tstorage)

for k, time_series in enumerate(time_series_stack):
    tstorage[k] += temporal_skw(time_series)
print("temporal skw:", tstorage)

for k, time_series in enumerate(time_series_stack):
    tstorage[k] += temporal_krt(time_series)
print("temporal krt:", tstorage)

for k, time_series in enumerate(time_series_stack):
    tstorage[k] += temporal_cv(time_series)
print("temporal CV:", tstorage)

finished_time = time.time() - start_time
print(f"-=- runtime is {finished_time} seconds -=-")

#####################################
