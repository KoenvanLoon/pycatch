import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
from statsmodels.graphics.tsaplots import plot_acf

# Seed the random number generator for dummy data
np.random.seed(0)

time_step = .1
time_vec = np.arange(0, 100, time_step)
time_window = 100

# Generate a signal with a small frequency chirp; no (critical) transition
signal = np.sin(2 * np.pi * 6 * time_vec) + np.random.randn(len(time_vec))
#plt.plot(time_vec, signal)
#plt.show()

# Calculate mean, variance, skewness, kurtosis
mean = np.mean(signal)
variance = np.var(signal)
skewness = scipy.stats.skew(signal)
kurtosis = scipy.stats.kurtosis(signal)
print(mean, variance, skewness, kurtosis)

# Calculate lag-1 autocorrelation in defined datawindow
data_windows = np.array_split(signal, time_window)

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)

mean_windowed = []
variance_windowed = []
skewness_windowed = []
kurtosis_windowed = []
lag1_autocor_windowed = []

for data_window in data_windows:
    mean_windowed.append(np.mean(data_window))

    variance_windowed.append(np.var(data_window))

    skewness_windowed.append(scipy.stats.skew(data_window))

    kurtosis_windowed.append(scipy.stats.kurtosis(data_window))

    pd_data = pd.Series(data_window)
    lag1_autocor_window = pd_data.autocorr(lag=1)
    lag1_autocor_windowed.append(lag1_autocor_window)

ax1.plot(time_vec, signal)
ax1.set_title('Data')
ax2.plot(mean_windowed)
ax2.set_title('Mean, window = 100 steps')
ax3.plot(variance_windowed)
ax3.set_title('Variance, window = 100 steps')
ax4.plot(skewness_windowed)
ax4.set_title('Skewness, window = 100 steps')
ax5.plot(kurtosis_windowed)
ax5.set_title('Kurtosis, window = 100 steps')
ax6.plot(lag1_autocor_windowed)
ax6.set_title('Lag-1 autocorrelation, window = 100 steps')
plt.tight_layout()
plt.show()

# Low-frequency power spectrum (LFPS)

## Compute Fourier transform of x
fourier_transform_signal = np.fft.rfft(signal)
abs_fourier_transform = np.abs(fourier_transform_signal)
power_spectrum_signal = np.square(abs_fourier_transform)
freq = np.linspace(0, (1/time_step)/2, len(power_spectrum_signal))
plt.plot(freq, power_spectrum_signal)
plt.show()

#plot_acf(signal)
#plt.show()
