import numpy as np
from pcraster import *
import configuration_weekly as cfg
import os
import matplotlib.pyplot as plt
import threading
import queue

def file_loader(variable, path='./1/', datatype='map'):
    data_stack = []
    if datatype=='map':
        data_stack = [pcr2numpy(readmap(path + i), np.NaN) for i in os.listdir(path)
                      if os.path.isfile(os.path.join(path, i)) and variable in i]
    if datatype=='numpy':
        data_stack = [np.loadtxt(path + i) for i in os.listdir(path)
                      if os.path.isfile(os.path.join(path, i)) and variable in i]
    return data_stack

def snapshot_loader(variable, timestep_str, path='./1/', datatype='map'):
    snapshot_stack = []
    for timestep_string in timestep_str:
        if datatype=='map':
            snapshot_stack = [pcr2numpy(readmap(path + i), np.NaN) for i in os.listdir(path)
                          if os.path.isfile(os.path.join(path, i)) and (variable in i) and (timestep_string in i)]
        if datatype=='numpy':
            snapshot_stack = [np.loadtxt(path + i) for i in os.listdir(path)
                          if os.path.isfile(os.path.join(path, i)) and (variable in i) and (timestep_string in i)]
    return snapshot_stack

def mean_time_series(stack_of_maps_as_list):
    mean_time_series = [0.0] * len(stack_of_maps_as_list)
    for k, map in enumerate(stack_of_maps_as_list):
        mean_time_series[k] += np.nanmean(map)
    return mean_time_series

def timesteps(start_shift, end_shift, steps_in_shift=10, steps_ba_shift=10):
    before_shift = np.linspace(start=0, stop=start_shift, num=steps_ba_shift, endpoint=False, dtype=int)
    during_shift = np.linspace(start=start_shift, stop=end_shift, num=steps_in_shift, endpoint=False, dtype=int)
    after_shift = np.linspace(start=end_shift, stop=cfg.numberOfTimeSteps, num=steps_ba_shift, dtype=int)
    timesteps = np.concatenate((before_shift, during_shift, after_shift))
    return timesteps

def timestep_str(timesteps):
    timestep_strings = []
    for time in timesteps:
        if len(str(time)) == 1:
            timestep_strings += "".join(".00" + str(time))
        elif len(str(time)) == 2:
            timestep_strings += "".join(".0" + str(time))
        elif len(str(time)) == 3:
            timestep_strings += "".join("." + str(time))
        else:
            time = list(str(time))
            time.insert(-3, ".")
            timestep_strings += "".join(time)
    return timestep_strings

### INPUTS & CALCULATIONS ###

variables = ["bioM"]

# Load bio mass and soil thickness timeseries
timeseries_bioM = file_loader("bioM", path="./1/", datatype="map")
mean_time_series_bioM = mean_time_series(timeseries_bioM)

# Plot & user input functions
## Threading to enter start & end point while plot is open
my_queue = queue.Queue()
def store_in_queue(f):
    def wrapper(*args):
        my_queue.put(f(*args))
    return wrapper

def plot():
    plt.plot(mean_time_series_bioM)
    plt.show()

@store_in_queue
def user_inputs():
    print("The total number of timesteps equals: ", cfg.numberOfTimeSteps, "steps.")
    start_shift = input("Enter starting point of shift:")
    end_shift = input("Enter end point of shift:")
    plt.close()
    return start_shift, end_shift

mthread = threading.Thread(target=user_inputs)
mthread.start()
plot()
shift_data = my_queue.get()
mthread.join()

## end of threading

start_shift = shift_data[0]
end_shift = shift_data[1]

# Calculations
timesteps_ = timesteps(int(start_shift), int(end_shift))
timestep_str_ = timestep_str(timesteps_)

# Loop over variables
#for variable in variables: TODO - Assign snapshots of different state variables to different arrays, or calculate EWS and save results
print("Started loading snapshots, this might take a while...")
snapshots_bioM = snapshot_loader("bioM", timestep_str_, path='./1/', datatype='map')
print(snapshots_bioM)
