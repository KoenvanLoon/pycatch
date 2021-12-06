import numpy as np
from scipy import ndimage, fft
import EWSPy as ews

data_ = np.random.random((100, 100))

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
"""
Combination of Jon Yearsley (2021). Generate AR1 spatial data (https://www.mathworks.com/matlabcentral/fileexchange/5099-generate-ar1-spatial-data), MATLAB Central File Exchange. Retrieved November 30, 2021.
and Dakos et al. 10.1073/pnas.0802430105
"""

data_ = np.random.random((100, 100)) * 10
I = 0 #ews.spatial_corr(data_) # Variable controlling autocorr (Moran's I)
e_stdev = 1. # Variance of the normally distributed error (usually 1.0)

#
sig2 = np.nanvar(data_) * (1 - I**2)
alpha0 = np.nanmean(data_) * (1 - I)

#
dim = data_.shape
N = np.prod(dim)

W = np.zeros((N, N)) # matrix with neighbours (rook neighborhood)
np.fill_diagonal(W[1:], I/4)
np.fill_diagonal(W[:, 1:], I/4)

M = np.identity(N) - W # (inv(M) == generator matrix)
inv_M = np.linalg.inv(M)

e = np.random.normal(loc=0.0, scale=e_stdev, size=N)

x = np.linalg.solve(M, e).reshape(dim) # + np.mean(data_) # analog to Jon Yearsley (2021). Generate AR1 spatial data (https://www.mathworks.com/matlabcentral/fileexchange/5099-generate-ar1-spatial-data), MATLAB Central File Exchange. Retrieved December 3, 2021.

x_test1 = np.dot(inv_M, e * np.sqrt(sig2)).reshape(dim) + alpha0

print("Mean, var, and I from data:", np.mean(data_), np.var(data_)) #, I)
print("Mean, var, and I from Yearsley (2021):", np.mean(x), np.var(x)) #, ews.spatial_corr(x))
print("Mean, var, and I from this research:", np.mean(x_test1), np.var(x_test1)) #, ews.spatial_corr(x_test1))
