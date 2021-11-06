import EWSPy as ews
from pcraster import *
import numpy as np
import os

rising_memory = True
rising_variability = True

window_size = 100
snapshot_interval = 100

###

path = './1/'
stack_of_maps_as_list = [pcr2numpy(readmap(path + i), np.NaN) for i in os.listdir(path) if os.path.isfile(os.path.join(path, i)) and 'bioT' in i]
# print(stack_of_maps_as_list)

###

max_time_series = ews.max_time_series(stack_of_maps_as_list)
mean_time_series = ews.mean_time_series(stack_of_maps_as_list)

stack_of_windows = ews.time_series2time_windows(max_time_series, window_size)

###

stack_of_snapshots = ews.time_series2snapshots(stack_of_maps_as_list, snapshot_interval)

###

print(
    "spatial mean", ews.spatial_mean(stack_of_snapshots),'\n',
    "spatial corr", ews.spatial_corr(stack_of_snapshots),'\n',
    # ews.spatial_DFT(stack_of_snapshots),'\n', # TODO - large output; how to plot?
    "spatial var", ews.spatial_var(stack_of_snapshots),'\n',
    "spatial skw", ews.spatial_skw(stack_of_snapshots),'\n',
    # ews.temporal_krt(stack_of_snapshots) #,'\n', # most of the time not included
    # ews.spatial_power_spec(stack_of_snapshots) # TODO - only works for square matrices; representative?
)

print(
    "temporal mean", ews.temporal_mean(stack_of_windows),'\n',
    "temporal std", ews.temporal_std(stack_of_windows),'\n',
    "temporal cv", ews.temporal_cv(stack_of_windows),'\n',
    # ews.temporal_skw(stack_of_windows),'\n',
    "temporal dfa", ews.temporal_dfa(stack_of_windows, scales=np.array([10, 5])),'\n',
    "temporal autocor", ews.temporal_autocorrelation(stack_of_windows),'\n',
    "temporal AR1", ews.temporal_AR1(stack_of_windows),'\n',
    "temporal return rate", ews.temporal_returnrate(stack_of_windows),'\n',
    "temporal cond het", ews.temporal_cond_het(stack_of_windows)
)
