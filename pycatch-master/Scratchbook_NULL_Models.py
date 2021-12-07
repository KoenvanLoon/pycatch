import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage, fft
import scipy.stats
import statsmodels.api
from statsmodels.nonparametric.kernel_regression import KernelReg
import EWSPy as ews

####

# Method 1
def method1_(data, replace=False):
    return np.random.choice(data, len(data), replace=replace)

# Method 2
def method2_(data, replace=False):
    fft_ = fft.fft(data_)
    fft_mag = np.abs(fft_)
    fft_phases = np.angle(fft_)

    fft_phases_new = fft_phases.copy()
    if len(data_) % 2 == 0:
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

    ifftm_ = fft.ifft(fft_sym)

    return ifft_m

# Method 3
def method3_(data):
    alpha1 = statsmodels.api.tsa.acf(data, nlags=1)
    sig2 = np.nanvar(data) * (1 - alpha1[1] ** 2)
    alpha0 = np.nanmean(data) * (1 - alpha1[1])
    e = np.random.normal(loc=0.0, scale=1.0, size=len(data))

    AR1m = np.zeros(len(data))
    for i in range(len(data)):
        AR1m[i] = alpha1[1] * AR1m[i - 1] + alpha0 + np.sqrt(sig2) * e[i]

    return AR1m

####

time_step = 0.02
period = 5.
time_vec = np.arange(0, 20, time_step)
sig = (np.sin(2 * np.pi / period * time_vec)
       + 0.5 * np.random.randn(time_vec.size))
data_ = sig

data_ = np.loadtxt('./1/bioA0005.200.numpy.txt')
data_org = np.copy(data_)

# kde = KernelReg(endog=data_, exog=np.arange(data_.shape[0]), var_type='c')
# estimator = kde.fit(np.arange(data_.shape[0]))
# estimator = np.reshape(estimator[0], data_.shape[0])
#
# data_ = data_ - estimator

sigma = 30
data_g1d = ndimage.gaussian_filter1d(data_, sigma)
data_ = data_ - data_g1d

#data_ = np.random.normal(loc=0, size=1000)

mean_ = np.mean(data_) # actual mean
var_ = np.var(data_) # actual variance


### Null models (Dakos et al. 2008) ###

## First method ##
"Similar probability distribution (mean and var)"
# Detrending & residual time series
# sigma = 1 # Estimate on data? TODO how2optimize, detrending of the original timeseries even necessary?
# data_g1d = ndimage.gaussian_filter1d(data_, sigma)
data_resid = data_ #- data_g1d

# Shuffle the (detrended) original time series
data_g1d_shuffled = np.random.choice(data_resid, len(data_resid), replace=False)
data_g1d_shuffled_replace = np.random.choice(data_resid, len(data_resid))

# # Bootstrap sampling
# data_g1_shuffled_list = data_g1d_shuffled.tolist()
# data_bootstrap = []
# for _ in range(50): # TODO which range, sample size?
#     data_bootstrap.append(np.mean(random.sample(data_g1_shuffled_list, 5)))
# mean_sample = np.mean(data_bootstrap)

## Second method ##
"Same autocorrelations and the same probability distribution" # TODO check lit; not same mean?
"Same Fourier spectrum and amplitudes as the original sets"
# Compute DFT, split into magnitude and phase
fft_ = fft.fft(data_)
fft_mag = np.abs(fft_)
fft_phases = np.angle(fft_)

# Randomize phases
"A real valued time domain signal has a conjugate symmetric freq. domain signal;" \
"hence, the right half of the freq. domain is determined as a sign-inverted of" \
"the left half in reverse order"

fft_phases_new = fft_phases.copy()
if len(data_) % 2 == 0:
    i = int(len(fft_phases_new)/2)
    fft_phases_left_half = fft_phases[1:i]
    # np.random.shuffle(fft_phases_left_half)
    # fft_shuffled_phases_lh = fft_phases_left_half
    fft_shuffled_phases_lh = np.random.choice(fft_phases_left_half, len(fft_phases_left_half), replace=False)
    fft_shuffled_phases_rh = - fft_shuffled_phases_lh[::-1]
    fft_phases_new = np.concatenate((np.array((fft_[0],)), fft_shuffled_phases_lh, np.array((fft_phases[i],)),
                                     fft_shuffled_phases_rh))
else:
    i = int(len(fft_phases_new)/2 + 1)
    fft_phases_left_half = fft_phases[1:i]
    fft_shuffled_phases_lh = np.random.choice(fft_phases_left_half, len(fft_phases_left_half), replace=False)
    fft_shuffled_phases_rh = - fft_shuffled_phases_lh[::-1]
    fft_phases_new = np.concatenate((np.array((fft_[0],)), fft_shuffled_phases_lh,
                                     fft_shuffled_phases_rh))

# Symmetrize the phases
fft_sym = fft_mag * (np.cos(fft_phases_new) + 1j * np.sin(fft_phases_new))

# Invert the DFT
ifft_ = fft.ifft(fft_sym)

## ADDITION - might not be needed
# if not np.allclose(ifft_.imag, np.zeros(ifft_.shape)):
#         max_imag = (np.abs(ifft_.imag)).max()
#         imag_str = '\nNOTE: a non-negligible imaginary component was discarded.\n\tMax: {}'
#         print(imag_str.format(max_imag))
# ifft_ = ifft_.real

## Third method ##
alpha1 = statsmodels.api.tsa.acf(data_, nlags=1)
sig2 = np.nanvar(data_) * (1 - alpha1[1]**2)
alpha0 = np.nanmean(data_) * (1 - alpha1[1])
e = np.random.normal(loc=0.0, scale=1.0, size=len(data_))

z = np.zeros(len(data_))
for i in range(len(data_)):
    z[i] = alpha1[1] * z[i-1] + alpha0 + np.sqrt(sig2) * e[i]

## Plots ##
fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2,2, figsize=(10, 10))

#ax1.plot(data_ + data_g1d)
ax1.plot(data_)
#ax1.plot(data_g1d)
#ax1.plot(data_ + estimator)
#ax1.plot(data_g1d)

#ax2.plot(data_)
ax2.plot(data_g1d_shuffled)
#ax2.plot(data_g1d_shuffled + data_g1d)
#ax2.plot(data_g1d_shuffled + estimator)

#ax3.plot(data_)
ax3.plot(ifft_)
#ax3.plot(ifft_ + data_g1d)
#ax3.plot(ifft_ + estimator)

# ax4.plot(fft.fft(data_))
# ax4.plot(fft.ifft(fft.fft(data_))) # ! ifft(fft(a)) == a to within numerical accuracy
#ax4.plot(data_)
ax4.plot(z)
#ax4.plot(z + data_g1d)
#ax4.plot(z + estimator)

plt.show()

a = ews.time_series2time_windows(data_)
a = ews.temporal_var(a)
print(a)

b = ews.time_series2time_windows(data_g1d_shuffled)
b = ews.temporal_var(b)
print(b)

kendall_tau_shuffle = scipy.stats.kendalltau(a, b, nan_policy='omit')
print(kendall_tau_shuffle)
