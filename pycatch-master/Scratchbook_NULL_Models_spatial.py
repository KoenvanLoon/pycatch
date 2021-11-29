import numpy as np
import random
import matplotlib.pyplot as plt
from scipy import interpolate, ndimage, fft
import statsmodels.api
import math

data_ = np.random.random((10,12))

## First method ##
# Save shape of dataset & flatten to 1D array
data_new = data_.copy()
data_shape = data_new.shape
data_1d = data_new.ravel()

# Detrending & residual time series
sigma = 1 # Estimate on data? TODO how2optimize, detrending of the original timeseries even necessary?
data_g1d = ndimage.gaussian_filter1d(data_1d, sigma)
data_resid = data_1d - data_g1d

data_shuffled = np.random.choice(data_resid, len(data_resid), replace=False).reshape(data_shape)
data_shuffled_replace = np.random.choice(data_resid, len(data_resid)).reshape(data_shape)

## Second method ##
#
fft2_ = fft.fft2(data_)
fft2_mag = np.abs(fft2_)
fft2_phases = np.angle(fft2_)

fft2_phases_new = fft2_phases.copy()
fft2_phases_shape = fft2_phases_new.shape
fft2_phases_1d = fft2_phases_new.ravel()
fft2_phases_new = np.random.choice(fft2_phases_1d, len(fft2_phases_1d), replace=False).reshape(fft2_phases_shape)

# Symmetrize the phases
fft2_sym = fft2_mag * (np.cos(fft2_phases_new) + 1j * np.sin(fft2_phases_new))

# Invert the DFT
ifft2_ = fft.ifft2(fft2_sym)

## Third method ##


####

start = 0
end = 5200
step = 100

start_shift, end_shift = 750, 901
num = 10

s_s = start_shift/step
e_s = end_shift/step
d_ss_es = e_s - s_s
print(d_ss_es)

if d_ss_es < num:
    num = math.ceil(d_ss_es)

print(num)

arange = np.arange(start, end+1, step)
linspace1 = np.linspace(start_shift, end_shift, num, dtype='int')
linspace2 = np.linspace(math.floor(start_shift/step), math.ceil(end_shift/step), num, dtype='int') * 100

print(arange)
print(linspace1)
print(linspace2)
