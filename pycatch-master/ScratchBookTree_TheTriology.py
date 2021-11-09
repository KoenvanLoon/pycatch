import EWSPy as ews
from pcraster import *
import numpy as np
import os
import time
import itertools
import threading
import sys
import configuration_weekly as cfg

### User input ### TODO - Only single variable/realization; needs automation for multiple variables/realizations

variables = ['bioM', 'demM']
# variable = 'bioM'

realizations = cfg.nrOfSamples
# realizations = 1

spatial_ews = True
temporal_ews = True

window_size = 100
snapshot_interval = window_size

### End of user input ###

### Loading files animation ###

def animated_loading():
    chars = "/-\|"
    for char in chars:
        sys.stdout.write('\r'+'loading '+char)
        time.sleep(.1)
        sys.stdout.flush()

### Loading files ###

def file_loader(variable, path='./1/', timer_on=False):
    if timer_on == True:
        start_time = time.time()

    stack_of_maps_as_list = [pcr2numpy(readmap(path + i), np.NaN) for i in os.listdir(path)
                             if os.path.isfile(os.path.join(path, i)) and variable in i]

    if timer_on == True:
        print('\n',"Loading finished with loading time:", time.time() - start_time,'\n')

    return stack_of_maps_as_list

### EWS  calculations ###

def ews_calculations(stack_of_maps_as_list, window_size=100, snapshot_interval=100, time_series='mean', timer_on=True):
    ### Temporal EWS ###
    if temporal_ews == True:
        print("Started temporal EWS calculations")
        if timer_on == True:
            start_time = time.time()

        input_time_series = []
        if time_series == 'mean':
            input_time_series = ews.mean_time_series(stack_of_maps_as_list)
        if time_series == 'max':
            input_time_series = ews.max_time_series(stack_of_maps_as_list)

        stack_of_windows = ews.time_series2time_windows(input_time_series, window_size)

        print(
            " temporal mean", ews.temporal_mean(stack_of_windows),'\n',
            "temporal std", ews.temporal_std(stack_of_windows),'\n',
            "temporal cv", ews.temporal_cv(stack_of_windows),'\n',
            "temporal skw", ews.temporal_skw(stack_of_windows),'\n',
            "temporal krt", ews.temporal_krt(stack_of_windows),'\n',
            "temporal dfa", ews.temporal_dfa(stack_of_windows, scales=np.array([10, 5])),'\n',
            "temporal autocor", ews.temporal_autocorrelation(stack_of_windows),'\n',
            "temporal AR1", ews.temporal_AR1(stack_of_windows),'\n',
            "temporal return rate", ews.temporal_returnrate(stack_of_windows),'\n',
            "temporal cond het", ews.temporal_cond_het(stack_of_windows)
        )

        if timer_on == True:
            end_time = time.time()
            print("Elapsed time for temporal EWS calculations equals:", end_time - start_time,'\n')

    ### Spatial EWS ###
    if spatial_ews == True:
        print("Started spatial EWS calculations")
        if timer_on == True:
            start_time = time.time()

        stack_of_snapshots = ews.time_series2snapshots(stack_of_maps_as_list, snapshot_interval)

        print(
            " spatial mean", ews.spatial_mean(stack_of_snapshots),'\n',
            "spatial corr", ews.spatial_corr(stack_of_snapshots),'\n',
            # ews.spatial_DFT(stack_of_snapshots),'\n', # TODO - large output; how to plot?
            "spatial var", ews.spatial_var(stack_of_snapshots),'\n',
            "spatial skw", ews.spatial_skw(stack_of_snapshots),'\n',
            # ews.temporal_krt(stack_of_snapshots) #,'\n', # most of the time not included
            # ews.spatial_power_spec(stack_of_snapshots) # TODO - only works for square matrices; representative?
        )

        if timer_on == True:
            end_time = time.time()
            print("Elapsed time for spatial EWS calculations equals:", end_time - start_time,'\n')

    return 'EWS calculated succesfully'

###

for realization in range(1, realizations + 1):
    for variable in variables:
        process = threading.Thread(name='process', target=file_loader(variable))
        process.start()
        while process.isAlive():
            animated_loading()

        stack_of_maps_as_list = file_loader(variable, path=f'./{realization}/', timer_on=True)
        run = ews_calculations(stack_of_maps_as_list, window_size=window_size, snapshot_interval=snapshot_interval,
                               time_series='mean', timer_on=True)

# done_loading = False
# threading.Thread(target=animate_loading(variable)).start()
# stack_of_maps_as_list = file_loader(variable, './1/', timer_on=True)
# done_loading = True
# run = ews_calculations(stack_of_maps_as_list, window_size=window_size, snapshot_interval=snapshot_interval,
#                        time_series='mean', timer_on=True)

###
