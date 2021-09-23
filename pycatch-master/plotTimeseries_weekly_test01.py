import matplotlib.pyplot as plt
import numpy
from statsmodels.graphics.tsaplots import plot_acf
import pandas as pd

# script to create plots for timeseries output of main_weekly.py

# first item: numpy array; first index time steps, second index samples, third row number, fourth column number
# second item: unit for the samples inside the numpy array
unit_dict = {
    'gA': 'unit of grazing rate',  # GRAZING RATE # Added
    'bioA': 'unit of biomass',  # BIOMASS # Added
    'bioS': 'unit of biomass var-S',
    'bioT': 'unit of biomass var-T',
    'regA': 'unit of regolith thickness',  # REGOLITH THICKNESS # Added
    'sfA': 'unit of moisture content',  # moisture content
    'qA': 'unit of discharge',  # DISCHARGE # Added
    'gpA': 'unit of growth part',  # growth part
    'grA': 'unit of grazing',  # grazing
    'grNA': 'unit of net growth',  # net growth
    'depA': 'unit of net deposition',  # net deposition
    'weaA': 'unit of net weathering',  # net weathering
    'creA': 'unit of net deposition due to creep',  # net deposition due to creep
    'demA': 'unit of DEM'  # DEM # Added
}

var_dict = {
    filename: (numpy.load(f'{filename}.npy'), desc)
    for filename, desc in unit_dict.items()
}


def plot_func(ax, var_str, location=None):
    ax.set_xlabel("Weeks")  # hardcoded in main_weekly
    ax.set_ylabel(f"{var_dict.get(var_str[1])}")  # changes per var; edit in unit_dict
    if location is None:
        ax.set_title(f'{var_str} semivariance for two max_lags')
        ax.plot(var_dict.get(var_str)[0][:, sample, row])
    else:
        ax.set_title(f'{var_str} at location {location}')
        ax.plot(var_dict.get(var_str)[0][:, sample, row, location])


# Making the plots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

# time steps (= all; per definition), samples (= one; tbd),
# row number (= one; tbd), column number (= up & down; tbd in func call)
sample = 0
row = 0

# plot_func(subplot, variable name, location in area int; None by default)
plot_func(ax1, 'bioA', 0)
plot_func(ax2, 'bioA', 7)
plot_func(ax3, 'bioS')
plot_func(ax4, 'bioT')

plt.tight_layout()
#plt.show()

# Correlogram & lag-1 autocorrelation
plot_acf(var_dict.get('bioA')[0][:, sample, row, 7])
plt.show()

data_windows = numpy.array_split(var_dict.get('bioA')[0][:, sample, row, 7], 5)
print(len(var_dict.get('bioA')[0][:, sample, row, 7]), len(data_windows))

lag1_autocor_windowed = []
for data_window in data_windows:
    pd_data = pd.Series(data_window)
    lag1_autocor_window = pd_data.autocorr(lag=1)
    lag1_autocor_windowed.append(lag1_autocor_window)
plt.plot(lag1_autocor_windowed)

#pd.plotting.autocorrelation_plot(var_dict.get('bioA')[0][:, sample, row, 7])
plt.show()

fig.savefig("plot_timeseries_weekly_test01.pdf", format="pdf")
