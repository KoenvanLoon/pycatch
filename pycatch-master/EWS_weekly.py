import EWSPy as ews
from pcraster import *
import numpy as np
import os
import time
import scipy.stats
from scipy import ndimage

import configuration_weekly as cfg
import NULL_models_timeseries as temp_NULL
import NULL_models_spatial as spat_NULL

### User input ###

## State variables for EWS ##
"ews.StateVariable(file name, spatial ews, temporal ews, snapshot interval, window size, window overlap, datatype 'map'/'numpy'(.txt))"

sfM = ews.StateVariable("sfM", spatial=True)
bioM = ews.StateVariable("bioM", spatial=True)
bioA = ews.StateVariable("bioA", temporal=True, datatype='numpy')
regM = ews.StateVariable("regM", spatial=True)
demM = ews.StateVariable("demM", spatial=True)
qA = ews.StateVariable("qA", temporal=True, datatype='numpy')
gM = ews.StateVariable("gM", spatial=True)

#variables = [sfM, bioM, regM, demM, qA, gM]
variables = [bioA, bioM]

### End user input ###

## Realizations/MC runs ##
realizations = cfg.nrOfSamples
# realizations = 1 # for test cases

## Timesteps, intervals ##
spatial_ews_present = cfg.map_data
spatial_ews_interval = np.arange(cfg.interval_map_snapshots, cfg.numberOfTimeSteps + cfg.interval_map_snapshots,
                                 cfg.interval_map_snapshots)

temporal_ews_present = cfg.mean_timeseries_data
temporal_ews_interval = cfg.numberOfTimeSteps

## Functions ##

def time_series2time_windows(time_series, window_size=100, window_overlap=0):
    return np.array([time_series[i:i + window_size] for i in range(0, len(time_series), window_size - window_overlap)])


def generate_datasets(variable, path='./1/', nr_realizations=1, detrending_temp='None', sigma=50, method1=False,
                      method2=False, method3=False):

    if variable.temporal:
        ## Load data ##
        state_variable_timeseries = []
        if variable.datatype == 'numpy':
            file_name = ews.file_name_str(variable.name, temporal_ews_interval)
            state_variable_timeseries = np.loadtxt(path + file_name + ".numpy.txt")
        else:
            print(f"Datatype for {variable.name} currently not supported.")

        ## Detrending: 'None', 'Gaussian' ##
        if detrending_temp == 'None':
            state_variable_timeseries = state_variable_timeseries
        if detrending_temp == 'Gaussian': # TODO - Multiple sigmas?
            state_variable_timeseries = state_variable_timeseries - ndimage.gaussian_filter1d(state_variable_timeseries, sigma)

        ## Generate dummy datasets ##
        if method1:
            temp_NULL.method1_(state_variable_timeseries, realizations=nr_realizations, path=path,
                               file_name=variable.name)
        if method2:
            temp_NULL.method2_(state_variable_timeseries, realizations=nr_realizations, path=path,
                               file_name=variable.name)
        if method3:
            temp_NULL.method3_(state_variable_timeseries, realizations=nr_realizations, path=path,
                               file_name=variable.name)

    if variable.spatial:
        ## Load data ##
        state_variable_snapshots = [0.0] * len(spatial_ews_interval)
        if variable.datatype == 'numpy':
            for k, timestep in enumerate(spatial_ews_interval):
                file_name = ews.file_name_str(variable.name, timestep)
                state_variable_snapshots[k] = np.loadtxt(path + file_name + 'numpy.txt')
        if variable.datatype == 'map':
            for k, timestep in enumerate(spatial_ews_interval):
                file_name = ews.file_name_str(variable.name, timestep)
                state_variable_snapshots[k] = pcr2numpy(readmap(path + file_name), np.NaN)
        else:
            print(f"Datatype for {variable.name} currently not supported.")
        state_variable_snapshots = np.asarray(state_variable_snapshots)

        ## Generate dummy datasets ##
        if method1:
            spat_NULL.method1_(state_variable_snapshots, realizations=nr_realizations, path=path,
                               file_name=variable.name)
        if method2:
            spat_NULL.method2_(state_variable_snapshots, realizations=nr_realizations, path=path,
                               file_name=variable.name)
        if method3:
            spat_NULL.method3_(state_variable_snapshots, realizations=nr_realizations, path=path,
                               file_name=variable.name)


def ews_calculations(variable, path='./1/', timer_on=False):
    ## Temporal EWS calculations ##
    if variable.temporal:
        if temporal_ews_present:
            print("Started temporal EWS calculations")

            ## Start timer if set to True ##
            if timer_on:
                start_time = time.time()

            ## Timeseries file loading ##
            state_variable_timeseries = []
            if variable.datatype == 'numpy':
                file_name = ews.file_name_str(variable.name, temporal_ews_interval)
                state_variable_timeseries = np.loadtxt(path + file_name + ".numpy.txt")
            else:
                print(f"Datatype for {variable.name} currently not supported.")

            ## Splitting timeseries into (overlapping) windows ##
            stack_of_windows = time_series2time_windows(state_variable_timeseries, variable.window_size,
                                                        variable.window_overlap)

            ## EWS calculations ###
            print("temporal mean", ews.temporal_mean(stack_of_windows))

            # print(
            #     "temporal mean", ews.temporal_mean(stack_of_windows), '\n',
            #     "temporal std", ews.temporal_std(stack_of_windows), '\n',
            #     "temporal cv", ews.temporal_cv(stack_of_windows), '\n',
            #     "temporal skw", ews.temporal_skw(stack_of_windows), '\n',
            #     "temporal krt", ews.temporal_krt(stack_of_windows), '\n',
            #     "temporal dfa", ews.temporal_dfa(stack_of_windows, scales=np.array([10, 5])), '\n',
            #     "temporal autocorr", ews.temporal_autocorrelation(stack_of_windows), '\n',
            #     "temporal AR1", ews.temporal_AR1(stack_of_windows), '\n',
            #     "temporal return rate", ews.temporal_returnrate(stack_of_windows), '\n',
            #     "temporal cond het", ews.temporal_cond_het(stack_of_windows)
            # )

            ## End timer if set to True##
            if timer_on:
                end_time = time.time()
                print("Elapsed time for temporal EWS calculations equals:", end_time - start_time, '\n')

        elif temporal_ews_present == False:
            print(f"Mean timeseries data == False in configuration_weekly.py, could not calculate EWS for {variable.name}.")


# for realization in range(1, realizations + 1):
#     for variable in variables:
#         ews_calculations(variable, path=f'./{realization}/')

for realization in range(1, realizations + 1):
    for variable in variables:
        generate_datasets(variable, path=f'./{realization}/', detrending_temp='Gaussian', method1=True, method2=False, method3=True)
