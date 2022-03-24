"""
EWS - Early Warning Signals
EWS Tests

@authors: KoenvanLoon & TijmenJanssen
"""

import numpy as np
import os
import scipy.stats
import matplotlib.pyplot as plt
from scipy import ndimage, signal
import time as timeit

import EWSPy as ews
import EWS_configuration as cfg
import EWS_StateVariables as ews_sv

# State variables for EWS
"""
State variables present in EWS_StateVariables.py can be added through EWS_configuration.py

Args:
----

variables : list of either the hourly or weekly variables from EWS_StateVariables as specified in EWS_configuration.

names : list of the names (full and shortened) of the variables specified in the hourly or weekly variables.

"""

variables = ews_sv.variables_weekly

names = []
for variable in variables:
    names.append([f'{variable.full_name} as {variable.name}'])

# Early warning signals names
"""
Early warning signals (both temporal and spatial) which are included in EWSPy.py

Args:
----

ews_temporal_signals : dict of shorthand notation and name of the temporal early warning signals.

ews_spatial_signals : dict of shorthand notation and name of the spatial early warning signals.

"""

ews_temporal_signals = {'mn': "mean", 'std': "standard deviation", 'var': "variance",
                        'cv': "coefficient of variation", 'skw': "skewness", 'krt': "kurtosis",
                        'dfa': "detrended fluctuation analysis", 'acr': "autocorrelation", 'AR1': "AR1",
                        'rr': "return rate", 'coh': "conditional heteroskedasticity", 'timeseries': "timeseries",
                        'gauss': "gauss"}
ews_spatial_signals = {'std': "standard deviation", 'var': "variance", 'skw': "skewness", 'krt': "kurtosis",
                       'mI': "Moran's I"}


# Kendall tau stats
"""
Returns the Kendall tau value and its significance (p-value).

Args:
----

state_variable : The state variable of interest.

sum_stat : str, the summary statistic for which the Kendall tau value is calculated.
 
comp2 : str, either 'Same' or 'Forcing', sets the comparison to the mean of the same state variable or the forcing 
    (grazing) rate

path : str, path where inputs from the hourly/weekly model are stored.

Returns:
----

tau : float, the Kendall tau value (rank correlation coefficient).

p : float, the p-value (significance) of the calculated Kendall tau value.

"""


def kendalltau_stats(state_variable, sum_stat, comp2='Same', path='./1/'):
    if state_variable.temporal:
        dim = '.t.'
    elif state_variable.spatial:
        dim = '.s.'

    tau, p = np.NaN, np.NaN
    fdict = os.path.join(path + state_variable.name + dim)
    if os.path.exists(fdict + sum_stat + '.numpy.txt'):
        X = np.loadtxt(fdict + sum_stat + '.numpy.txt')

        if comp2 == 'Same':  # Dakos et al, 2008
            Y = np.loadtxt(fdict + 'mn.numpy.txt')
        elif comp2 == 'Forcing':  # Dakos et al, 2011 - Does not work if window sizes are different for the forcing & SV
            Y = np.loadtxt('./1/gA.t.mn.numpy.txt')

        tau, p = scipy.stats.kendalltau(X, Y, nan_policy='propagate')

    return tau, p


# Kendall tau stats for dummy data
"""
Returns the Kendall tau value and its significance (p-value).

Args:
----

state_variable : The state variable of interest.

sum_stat : str, the summary statistic for which the Kendall tau value is calculated.

method : str, either 'm1g', 'm2g', or 'm3g', the dummy data generation method for which the Kendall tau value is 
    calculated.
 
comp2 : str, either 'Same' or 'Forcing', sets the comparison to the mean of the same state variable or the forcing 
    (grazing) rate

path : str, path where inputs from the hourly/weekly model are stored.

Returns:
----

tau : array, the Kendall tau values (rank correlation coefficient) for each realization of dummy data.

p : array, the p-value (significance) of the calculated Kendall tau values for each realization of dummy data.

"""


def kendalltau_stats_dummy(state_variable, sum_stat, method='m1g', comp2='Same', path='./1/'):
    if state_variable.temporal:
        dim = '.t.'
    elif state_variable.spatial:
        dim = '.s.'

    generated_number_length = ews.generated_number_length(cfg.nr_generated_datasets)

    taurray = [np.NaN] * cfg.nr_generated_datasets
    parray = [np.NaN] * cfg.nr_generated_datasets
    for realization in range(cfg.nr_generated_datasets):
        fdict = os.path.join(path + method + str(realization).zfill(generated_number_length) + '/'
                             + state_variable.name + dim)
        if os.path.exists(fdict + sum_stat + '.numpy.txt'):
            X = np.loadtxt(fdict + sum_stat + '.numpy.txt')

            if comp2 == 'Same':  # Dakos et al, 2008
                Y = np.loadtxt(fdict + 'mn.numpy.txt')
            elif comp2 == 'Forcing':  # Dakos et al, 2011 - Does not work if window sizes are different for the forcing & SV, detrending == Gaus
                Y = np.loadtxt('./1/gA.t.mn.numpy.txt')

            taurray[realization], parray[realization] = scipy.stats.kendalltau(X, Y, nan_policy='propagate')

    return taurray, parray


# Histogram plot maker
"""
Returns a histogram plot of the Kendall tau value(s) of the modelled dataset and dummy datasets.

Args:
----

variable : The state variable of interest.

statistic : str, the summary statistic for which the Kendall tau value is calculated.
 
method : str, either 'm1g', 'm2g', or 'm3g', the dummy data generation method for which the Kendall tau value is 
    calculated.
    
tau_values : array, the Kendall tau values (rank correlation coefficient) for each realization of dummy data.

tau_original : float, the Kendall tau value for the modelled dataset.

chance_cor : float, the p-value of the Kendall tau value for the modelled dataset.

path : str, path where inputs from the hourly/weekly model are stored.

Returns:
----

A histogram plot of the Kendall tau value(s) of the modelled dataset and dummy datasets. Optionally saved to disk.

"""


def histogram_plot(variable, statistic, method, tau_values, tau_original=np.NaN, chance_cor=np.NaN, path='./1/'):
    bins = np.linspace(-1, 1, num=41)
    histogram = np.histogram(tau_values, bins)
    print("Histogram values ([values],[bins]):", histogram)

    plt.tight_layout()
    plt.hist(tau_values, bins, ec='black')
    if tau_original != np.NaN:
        plt.axvline(x=tau_original, color='r')
        if variable.temporal:
            plt.title(f"Probability distribution of the trend statistic ({ews_temporal_signals[statistic]}) under {method} of {variable.full_name}")
        elif variable.spatial:
            plt.title(f"Probability distribution of the trend statistic ({ews_spatial_signals[statistic]}) under {method} for  of {variable.full_name}")
        plt.xlabel("K \u03C4")
        plt.ylabel("Frequency")
        plt.text(0.7, 0, f" \u03C4 = {round(tau_original, 3)} \n p = {round(chance_cor, 3)} \n")

    print("Save the plot as a .pdf? [Y/n]")
    save_plot = input()
    if save_plot == 'Y' or save_plot == 'y':
        plt.savefig(path + f"{variable.name}_{statistic}_{method}_histogram.pdf", format="pdf")

    plt.show()
    plt.close()


# Chance value
"""
The chance of the Kendall tau value exceeding the value found in the dummy datasets.

Args:
----

values : array, Kendall tau values of dummy datasets.

value : float, Kendall tau value of modelled dataset.

sign_ini : str, 'None' by default. Whether the value should be greater ('<') or smaller ('>') than values.

"""


def chance_value(values, value, sign_ini='None'):
    sign = sign_ini
    if sign == 'None':
        if value > 0.:
            sign = '>'
        elif value < 0.:
            sign = '<'
        else:
            print(f"Tau value of {value} encountered while sign is not specified. Defaulted to 'tau values' > 'tau value'.")
            sign = '>'

    if sign == '>':
        p = np.nansum(np.array(values) > value) / len(values)
    elif sign == '<':
        p = np.nansum(np.array(values) < value) / len(values)
    else:
        print(f"{sign} is not supported")

    return p


# Histogram plot maker user input
"""
Returns a histogram plot of the Kendall tau value(s) of the modelled dataset and dummy datasets by running the
    histogram_plot() function.

Args:
----

path : str, path where inputs from the hourly/weekly model are stored.

Returns:
----

optionally : A histogram plot of the Kendall tau value(s) of the modelled dataset and dummy datasets. Optionally saved 
    to disk.
    
optionally : Kendall tau and p-values otherwise shown in the histogram plot.

"""


def kendall_tau_valhist(path='./1/'):
    print("Variables present in the current run are:")
    for name in names:
        print(name)

    print("Enter the short name for the state variable:")
    state_variable_input = input()
    state_variable = [variable for variable in variables if variable.name == state_variable_input][0]

    if state_variable.temporal:
        print("EW signals present are:", ews_temporal_signals)
    elif state_variable.spatial:
        print("EW signals present are:", ews_spatial_signals)

    print("Enter the short name for the EW signal:")
    summary_statistic = input()
    tau_m, p_m = kendalltau_stats(state_variable=state_variable, sum_stat=summary_statistic)

    print("Do you want to make a histogram comparison graph? [Y/n]")
    histogram = input()
    if histogram == 'Y' or histogram == 'y':
        print("For which null model do you want to test? [m1g, m2g, m3g]")
        method = input()
        tau_values, p_values = kendalltau_stats_dummy(state_variable=state_variable, sum_stat=summary_statistic, method=method)
        chance_cor = chance_value(tau_values, tau_m)

        histogram_plot(variable=state_variable, statistic=summary_statistic, method=method, tau_values=tau_values, tau_original=tau_m, chance_cor=chance_cor, path=path)
    else:
        print("Print Kendall tau and p values? [Y/n]")
        print_on = input()
        if print_on == 'Y' or print_on == 'y':
            print(tau_m, p_m)

    print("Would you like to make another plot? [Y/n]")
    answer = input()
    if answer == 'Y' or answer == 'y':
        kendall_tau_valhist(path=path)
    else:
        print("Terminated histogram maker. Goodbye.")


# TODO - Variables sorted as in weekly_plots, move other things to cfg
# TODO - Change method names


# Time series to time windows
"""
Divides a time series (2D numpy array) into an array of evenly sized time windows (2D numpy arrays). If remaining data-
points do not fill the last time window, they are dropped from the stack of time windows.

Args:
-----

timeseries : A 2D numpy array containing data points of a early-warning signal.

window_size : The size (int) of the windows into which the time series is to be divided.

window_overlap : The number (int) of data points in the window equal to the last data points of the previous time 
    window.

Returns:
-----

view : A 3D numpy array containing evenly sized time windows (2D numpy arrays).

    ! - Note that the amount of data points in 'view' does not need to be equal to the amount of data points in 
    'timeseries' due to the possibility of dropping data points if they do not fill the last time window completely.

"""


def window(timeseries, window_size, window_overlap):
    actual_window_overlap = window_size - window_overlap
    sh = (timeseries.size - window_size + 1, window_size)
    st = timeseries.strides * 2
    if window_overlap != 0:
        view = np.lib.stride_tricks.as_strided(timeseries, strides=st, shape=sh)[::actual_window_overlap]
    elif window_overlap == 0:
        view = np.lib.stride_tricks.as_strided(timeseries, strides=st, shape=sh)[::window_size]
    return view


# Windowsize tests
"""
Makes a contourplot with Kendall tau and p-values for different window sizes and overlaps.

Args:
----

path : str, path where inputs from the hourly/weekly model are stored.

method : str, either 'None' or 'Linear'. If 'Linear', a linear detrending is performed for each window.

Returns:
----

Two contourplots with Kendall tau and p-values for different window sizes and overlaps. Optionally saved to disk.

"""


def test_windowsize(path='./1/', method='None'):
    # Loading files
    fname = ews.file_name_str('bioA', cfg.number_of_timesteps_weekly)
    fpath = os.path.join(path + fname)
    biomass_timeseries = np.loadtxt(fpath + '.numpy.txt')
    if cfg.cutoff:
        biomass_timeseries = biomass_timeseries[:cfg.cutoff_point]

    # Window sizes & overlap

    # # Normal
    # window_overlaps = np.arange(0, 1000, 10)
    # if cfg.cutoff:
    #     window_sizes = np.arange(1000, cfg.cutoff_point // 2 + 1, 10)
    # else:
    #     window_sizes = np.arange(1000, cfg.number_of_timesteps_weekly // 2 + 1, 10)

    # Zoom
    window_overlaps = np.arange(0, 50, 1)
    if cfg.cutoff:
        window_sizes = np.arange(52, 1560, 1)
    else:
        window_sizes = np.arange(52, 1560, 1)

    # X and Y coords
    x, y = np.meshgrid(window_overlaps, window_sizes)

    # Calculating and storing tau and p values
    tau_arr = np.zeros((len(window_sizes), len(window_overlaps)))
    p_arr = np.zeros((len(window_sizes), len(window_overlaps)))
    for i, window_size in enumerate(window_sizes):
        print(f"Moved to windowsize {window_size}")
        for j, window_overlap in enumerate(window_overlaps):
            stack_of_windows = window(biomass_timeseries, window_size, window_overlap)

            if method == 'Linear':
                stack_of_windows = signal.detrend(stack_of_windows)

            mean = np.nanmean(stack_of_windows, axis=1)
            stat = ews.temporal_autocorrelation(stack_of_windows)

            tau, p = scipy.stats.kendalltau(stat, mean, nan_policy='propagate')
            # tau_arr[i, j] = tau
            # p_arr[i, j] = p

            if np.abs(tau) >= 0.3:
                tau_arr[i, j] = tau
            else:
                tau_arr[i, j] = np.nan

            if p <= 0.05:
                p_arr[i, j] = p
            else:
                p_arr[i, j] = np.nan

    # Making the plots
    z_array = [tau_arr, p_arr]
    fig, axs = plt.subplots(ncols=2)
    for ax, z in zip(axs, z_array):
        # Making the contour plot
        if z.all() == tau_arr.all():
            ax.set_title('Biomass autocorrelation windowtest \u03C4-value')
            # # Max value(s)
            # max_value = np.max(z)
            # local_max_index = np.where(z == max_value)
            # max_x, max_y = x[local_max_index[0], local_max_index[1]], y[local_max_index[0], local_max_index[1]]
            # ax.plot(max_x, max_y, color='red', marker="v", zorder=10, markersize=10, clip_on=False)
        elif z.all() == p_arr.all():
            ax.set_title('Biomass autocorrelation windowtest p-value')
        cs = ax.contourf(x, y, z, 10)
        #ax.contour(cs, colors='k')

        # Creating z coord for plot
        x_flat, y_flat, z_flat = x.flatten(), y.flatten(), z.flatten()
        def fmt(x, y):
            dist = np.linalg.norm(np.vstack([x_flat - x, y_flat - y]), axis=0)
            index = np.argmin(dist)
            z = z_flat[index]
            return 'x={x:.5f}   y={y:.5f}   z={z:.5f}'.format(x=x, y=y, z=z)
        fig.gca().format_coord = fmt

        # Labels, colorbar and grid
        ax.set_ylabel('Window size')
        ax.set_xlabel('Window overlap')
        ax.grid(c='k', alpha=0.3)
        fig.colorbar(cs, ax=ax)

        # # Min values
        # min_value = np.min(z)
        # local_min_index = np.where(z == min_value)
        # min_x, min_y = x[local_min_index[0], local_min_index[1]], y[local_min_index[0], local_min_index[1]]
        # ax.plot(min_x, min_y, color='blue', marker="v", zorder=10, markersize=10, clip_on=False)

    print("Save the plot as a .pdf? [Y/n]")
    save_plot = input()
    if save_plot == 'Y' or save_plot == 'y':
        plt.savefig(path + f"{method}_windowsize_test_var.pdf", format="pdf")

    plt.show()


# Windowsize and gaussian filtering tests
"""
Makes a contourplot with Kendall tau and p-values for different window sizes and Gaussian filter sizes.

Args:
----

path : str, path where inputs from the hourly/weekly model are stored.

Returns:
----

Two contourplots with Kendall tau and p-values for different window sizes and Gaussian filter sizes. Optionally saved to 
    disk.

"""


def test_windowgauss(path='./1/'):
    # Loading files
    fname = ews.file_name_str('bioA', cfg.number_of_timesteps_weekly)
    fpath = os.path.join(path + fname)
    biomass_timeseries = np.loadtxt(fpath + '.numpy.txt')
    if cfg.cutoff:
        biomass_timeseries = biomass_timeseries[:cfg.cutoff_point]

    # Window sizes & overlap

    # # Normal - Don't even try to run this, it takes ages
    # gaussian_filter = np.arange(0, 1000, 10)
    # if cfg.cutoff:
    #     window_sizes = np.arange(1000, cfg.cutoff_point // 2 + 1, 10)
    # else:
    #     window_sizes = np.arange(1000, cfg.number_of_timesteps_weekly // 2 + 1, 10)

    # Zoom
    gaussian_filter = np.arange(0, 1000, 10)
    if cfg.cutoff:
        window_sizes = np.arange(100, 15000, 10)
    else:
        window_sizes = np.arange(100, 15000, 10)

    # X and Y coords
    x, y = np.meshgrid(gaussian_filter, window_sizes)

    # Calculating and storing tau and p values
    tau_arr = np.zeros((len(window_sizes), len(gaussian_filter)))
    p_arr = np.zeros((len(window_sizes), len(gaussian_filter)))

    for j, gauss_filter in enumerate(gaussian_filter):
        print(f"Moved to Gaussian filter-size {gauss_filter}")
        filter = ndimage.gaussian_filter1d(biomass_timeseries, gauss_filter)
        biomass_timeseries_detrended = biomass_timeseries - filter
        for i, window_size in enumerate(window_sizes):

            stack_of_windows = window(biomass_timeseries_detrended, window_size, 0)

            mean = np.nanmean(stack_of_windows, axis=1)
            stat = ews.temporal_var(stack_of_windows)

            tau, p = scipy.stats.kendalltau(stat, mean, nan_policy='propagate')
            # tau_arr[i, j] = tau
            # p_arr[i, j] = p

            if np.abs(tau) >= 0.2:
                tau_arr[i, j] = tau
            else:
                tau_arr[i, j] = np.nan

            if p <= 0.1:
                p_arr[i, j] = p
            else:
                p_arr[i, j] = np.nan

    # Making the plots
    z_array = [tau_arr, p_arr]
    fig, axs = plt.subplots(ncols=2)
    for ax, z in zip(axs, z_array):
        # Making the contour plot
        if z.all() == tau_arr.all():
            ax.set_title('Biomass variance window-gaussian-test \u03C4-value')
            # # Max value(s)
            # max_value = np.max(z)
            # local_max_index = np.where(z == max_value)
            # max_x, max_y = x[local_max_index[0], local_max_index[1]], y[local_max_index[0], local_max_index[1]]
            # ax.plot(max_x, max_y, color='red', marker="v", zorder=10, markersize=10, clip_on=False)
        elif z.all() == p_arr.all():
            ax.set_title('Biomass variance window-gaussian-test p-value')
        cs = ax.contourf(x, y, z, 10)
        #ax.contour(cs, colors='k')

        # Creating z coord for plot
        x_flat, y_flat, z_flat = x.flatten(), y.flatten(), z.flatten()
        def fmt(x, y):
            dist = np.linalg.norm(np.vstack([x_flat - x, y_flat - y]), axis=0)
            index = np.argmin(dist)
            z = z_flat[index]
            return 'x={x:.5f}   y={y:.5f}   z={z:.5f}'.format(x=x, y=y, z=z)
        fig.gca().format_coord = fmt

        # Labels, colorbar and grid
        ax.set_ylabel('Window size')
        ax.set_xlabel('Gaussian filter size')
        ax.grid(c='k', alpha=0.3)
        fig.colorbar(cs, ax=ax)

        # # Min values
        # min_value = np.min(z)
        # local_min_index = np.where(z == min_value)
        # min_x, min_y = x[local_min_index[0], local_min_index[1]], y[local_min_index[0], local_min_index[1]]
        # ax.plot(min_x, min_y, color='blue', marker="v", zorder=10, markersize=10, clip_on=False)

    print("Save the plot as a .pdf? [Y/n]")
    save_plot = input()
    if save_plot == 'Y' or save_plot == 'y':
        plt.savefig(path + f"window_gauss_test_var.pdf", format="pdf")

    plt.show()


#test_windowgauss(path='./1/')
test_windowsize(path='./1/')
#kendall_tau_valhist()
