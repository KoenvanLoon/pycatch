import EWSPy as ews
from pcraster import *
import numpy as np
import os
import time
import itertools
import threading
import sys

### User input ### TODO - Only single variable; needs automation for multiple variables

variable = 'bioM'
variables = ['bioM', 'demM']

spatial_ews = True
temporal_ews = True

window_size = 100
snapshot_interval = 100

### End of user input ###

### Loading files animation ###

def animate_loading():
    for i in itertools.cycle(['.', '..', '...']): # ['|', '/', '-', '\\']
        if done:
            break
        sys.stdout.write(f'\rLoading {variable} files ' + i)
        sys.stdout.flush()
        time.sleep(0.1)

### Loading files ###

start_time = time.time()
done = False

t = threading.Thread(target=animate_loading).start()

path = './1/'
stack_of_maps_as_list = [pcr2numpy(readmap(path + i), np.NaN) for i in os.listdir(path) if os.path.isfile(os.path.join(path, i)) and variable in i]

done = True
end_time = time.time() - start_time

print('\n',"Loading finished with loading time:", time.time() - start_time,'\n')

### Temporal EWS ###

if temporal_ews == True:
    print("Started temporal EWS calculations")
    start_time = time.time()

    max_time_series = ews.max_time_series(stack_of_maps_as_list)
    mean_time_series = ews.mean_time_series(stack_of_maps_as_list)

    stack_of_windows = ews.time_series2time_windows(mean_time_series, window_size)

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

    end_time = time.time()
    print("Elapsed time for temporal EWS calculations equals:", end_time - start_time,'\n')

### Spatial EWS ###

if spatial_ews == True:
    print("Started spatial EWS calculations")
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

    end_time = time.time()
    print("Elapsed time for spatial EWS calculations equals:", end_time - start_time,'\n')
