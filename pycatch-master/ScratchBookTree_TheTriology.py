import EWSPy as ews
from pcraster import *
import numpy as np
import os

rising_memory = True
rising_variability = True

window_size = 10
snapshot_interval = window_size

###

path = './1/'
stack_of_maps_as_list = [pcr2numpy(readmap(path + i), np.NaN) for i in os.listdir(path) if os.path.isfile(os.path.join(path, i)) and 'Mfs' in i]
# print(stack_of_maps_as_list)

###

max_time_series = ews.max_time_series(stack_of_maps_as_list)
mean_time_series = ews.mean_time_series(stack_of_maps_as_list)

stack_of_windows = ews.time_series2time_windows(max_time_series, window_size)

###

stack_of_snapshots = ews.time_series2snapshots(stack_of_maps_as_list, snapshot_interval)

###

print(
    ews.temporal_mean(stack_of_windows),'\n',
    ews.temporal_std(stack_of_windows),'\n',
    ews.temporal_cv(stack_of_windows),'\n',
    ews.temporal_skw(stack_of_windows)
)

print(
    # stack_of_snapshots,'\n',
    ews.spatial_mean(stack_of_snapshots),'\n',
    #len(ews.spatial_mean(stack_of_snapshots)),'\n',
    #type(ews.spatial_mean(stack_of_snapshots)),'\n',
    #ews.spatial_corr(stack_of_snapshots),'\n', # needs work!
    ews.spatial_var(stack_of_snapshots)
)

for i in stack_of_snapshots:
    print(ews.spatial_corr(i))
