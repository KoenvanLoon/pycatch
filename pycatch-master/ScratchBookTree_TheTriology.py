import EWSPy as ews
from pcraster import *
import numpy as np
import os

rising_memory = True
rising_variability = True

window_size = 10
snapshot_interval = 10

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
    ews.spatial_mean(stack_of_snapshots),'\n',
    ews.spatial_corr(stack_of_snapshots),'\n',
    # ews.spatial_DFT(stack_of_snapshots),'\n', # TODO - large output; how2plot?
    ews.spatial_var(stack_of_snapshots),'\n',
    ews.spatial_skw(stack_of_snapshots),'\n',
    # ews.temporal_krt(stack_of_snapshots) #,'\n', # most of the time not included
    # ews.spatial_power_spec(stack_of_snapshots)
)

print(
    ews.temporal_mean(stack_of_windows),'\n',
    ews.temporal_std(stack_of_windows),'\n',
    ews.temporal_cv(stack_of_windows),'\n',
    ews.temporal_skw(stack_of_windows),'\n',
    ews.temporal_dfa(stack_of_windows, scales=np.array([5,10]))
)

# print(stack_of_windows)
# print(len(stack_of_windows))
# print(len(stack_of_windows[0]))

# for i in stack_of_snapshots:
#     print(ews.spatial_power_spec(i))
