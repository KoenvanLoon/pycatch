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
# Datasets should be composed of snapshots,
# each of which is a two-dimensional space
# discretized into M and N units (in x and y
# directions).
# Therefore, a total of M x N units of equal
# size is expected for each snapshot.
# z[m,n] is the value of the local state
# variable at location p = (m,n), with
# m or n = 1, 2, ..., M or N.
#####################################

def spatial_mean(numpy_matrix):
    return numpy_matrix.nanmean()

def spatial_std(numpy_matrix):
    return numpy_matrix.nanstd()

def spatial_var(numpy_matrix):
    return numpy_matrix.nanvar()

# def spatial_skw(numpy_matrix):
#     number = numpy_matrix.size
#     mean = spatial_mean(numpy_matrix)
#     std = spatial_std((numpy_matrix))
#     return (sum([(((z - mean) ** 3) / (std ** 3)) for z in np.nditer(numpy_matrix)]) / number)

def spatial_skw(numpy_matrix):
    return scipy.stats.skew(np.nditer(numpy_matrix), nan_policy='omit')

# def spatial_krt(numpy_matrix):
#     number = numpy_matrix.size
#     mean = spatial_mean(numpy_matrix)
#     std = spatial_std((numpy_matrix))
#     return (sum([(((z - mean) ** 4) / (std ** 4)) for z in np.nditer(numpy_matrix)]) / number - 3)

def spatial_krt(numpy_matrix):
    return scipy.stats.kurtosis(np.nditer(numpy_matrix), nan_policy='omit')

### Experimental functions ###

## spatial correlation ##

rook_neighborhood = np.array([
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 0]
])

def spatial_corr_complete(numpy_matrix):
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
    numpy_matrix *= sum_neighbours
    n_neighbours_times_avg = convolve2d(is_not_nan_as_nr, rook_neighborhood * mean, mode='same')
    n_neighbours_times_avg[is_nan] = 0
    P1 = sum(sum(numpy_matrix - n_neighbours_times_avg))

    # spatial_corr = storage / (4 * var * ((numpy_matrix.shape[0]-2) * (numpy_matrix.shape[1]-2)))
    spatial_corr = P1 / (4 * var * np.count_nonzero(is_not_nan_as_nr))
    return spatial_corr

## end of spatial corr. ##

def spatial_corr_inside(numpy_matrix): # Moran correlation from "EWS of Ecological Transitions: Methods for spatial patterns"
                                # - Does not seem to work properly (?), spatial_corr_other seems to peform better (but
                                # takes forever).
    mean = spatial_mean(numpy_matrix)
    var = spatial_var(numpy_matrix)
    storage = 0.0
    for m in range(1, numpy_matrix.shape[0]-1):
        for n in range(1, numpy_matrix.shape[1]-1):
            storage += (numpy_matrix[m, n] - mean) * (numpy_matrix[m, n-1] + numpy_matrix[m, n+1] + numpy_matrix[m-1, n]
                                                      + numpy_matrix[m+1, n]) - (4 * mean)
    spatial_corr = storage / (4 * var * ((numpy_matrix.shape[0]-2) * (numpy_matrix.shape[1]-2)))
    return spatial_corr

def spatial_corr_package(numpy_matrix):
    weights = lat2W(nrows=numpy_matrix.shape[0], ncols=numpy_matrix.shape[1], rook=True)
    mi = Moran_Local(numpy_matrix, weights)
    return sum(mi.Is)

#####################################

# from pcraster.framework import *
# setclone("inputs_weekly/clone.map")
# clone = boolean("inputs_weekly/clone.map")
# print(pcr2numpy(clone, np.nan))

np.random.seed(42)
data = np.random.rand(1000, 1000)

import time
start_time = time.time()

print(spatial_mean(data))
print(spatial_std(data))
print(spatial_var(data))
print(spatial_skw(data))
print(spatial_krt(data))

finished_time = time.time() - start_time
print(f"--- runtime is {finished_time} seconds ---")

start_time = time.time()
print(spatial_corr_inside(data))
finished_time = time.time() - start_time
print(f"--- runtime is {finished_time} seconds ---")

start_time = time.time()
print(spatial_corr_complete(data))
finished_time = time.time() - start_time
print(f"--- runtime is {finished_time} seconds ---")

# start_time = time.time()
# print(spatial_corr_other(data)) # Very slow (1h15m for matrix of 1000x1000 (one-million) values)
# finished_time = time.time() - start_time
# print(f"--- runtime is {finished_time} seconds ---")

#####################################

def mean(data):
    number = len(data)
    mean = sum(data) / number
    return mean

def variance(data):
    number = len(data)
    mean = sum(data) / number
    variance = sum([(x - mean) ** 2 for x in data]) / number
    return variance

def varianceSample(data):
    number = len(data)
    mean = sum(data) / number
    variance = sum([(x - mean) ** 2 for x in data]) / (number - 1)
    error = (2 * (variance ** 2) / (number - 1)) ** (1/2)
    return variance, error

def skewness(data):
    number = len(data)
    mean = sum(data) / number
    numerator = sum([(x - mean) ** 3 for x in data]) / number
    denominator = sum([(x - mean) ** 2 for x in data]) / number
    skewness = numerator / (denominator ** (3/2))
    return skewness

def skewnessSample(data):
    number = len(data)
    mean = sum(data) / number
    numerator = sum([(x - mean) ** 3 for x in data]) / number
    denominator = sum([(x - mean) ** 2 for x in data]) / (number - 1)
    skewness = numerator / (denominator ** (3/2))
    error = (6 * (number - 2) / ((number + 1) * (number + 3))) ** (1/2)
    return skewness, error

def kurtosis(data):
    number = len(data)
    mean = sum(data) / number
    numerator = sum([(x - mean) ** 4 for x in data]) / number
    denominator = sum([(x - mean) ** 2 for x in data]) / number
    kurtosis = numerator / (denominator ** 2) - 3
    return kurtosis

def kurtosisSample(data):
    number = len(data)
    mean = sum(data) / number
    numerator = sum([(x - mean) ** 4 for x in data]) / number
    denominator = sum([(x - mean) ** 2 for x in data]) / number
    kurtosis = numerator / (denominator ** 2) - 3
    error = (
            (24 * (number - 2) * (number - 3))
            / (((number + 1) ** 2) * (number + 3) * (number + 5))
    ) ** (1/2)
    return kurtosis, error

def autocovariance(data, lag=1):
    number = len(data)
    mean = sum(data) / number
    autocovariance = sum([(data[i] - mean) * (data[i+lag] - mean) for i in range(number - lag)]) / number
    return autocovariance

def autocorrelation(data, lag=1):
    return autocovariance(data, lag=lag) / variance(data)
