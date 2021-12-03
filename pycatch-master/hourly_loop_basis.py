import numpy as np
from pcraster import *
import configuration_weekly as cfg
import os
import matplotlib.pyplot as plt
import threading
import queue
import math

# def file_loader(variable, path='./1/', datatype='map'):
#     data_stack = []
#     if datatype=='map':
#         data_stack = [pcr2numpy(readmap(path + i), np.NaN) for i in os.listdir(path)
#                       if os.path.isfile(os.path.join(path, i)) and variable in i]
#     if datatype=='numpy':
#         data_stack = [np.loadtxt(path + i) for i in os.listdir(path)
#                       if os.path.isfile(os.path.join(path, i)) and variable in i]
#     return data_stack
#
# def snapshot_loader(variable, timestep_str, path='./1/', datatype='map'):
#     snapshot_stack = []
#     for timestep_string in timestep_str:
#         if datatype=='map':
#             snapshot_stack = [pcr2numpy(readmap(path + i), np.NaN) for i in os.listdir(path)
#                           if os.path.isfile(os.path.join(path, i)) and (variable in i) and (timestep_string in i)]
#         if datatype=='numpy':
#             snapshot_stack = [np.loadtxt(path + i) for i in os.listdir(path)
#                           if os.path.isfile(os.path.join(path, i)) and (variable in i) and (timestep_string in i)]
#     return snapshot_stack
#
# def mean_time_series(stack_of_maps_as_list):
#     mean_time_series = [0.0] * len(stack_of_maps_as_list)
#     for k, map in enumerate(stack_of_maps_as_list):
#         mean_time_series[k] += np.nanmean(map)
#     return mean_time_series

def timesteps(start_shift, end_shift, steps_in_shift=10, steps_ba_shift=10):
    first_state_duration = (start_shift/cfg.interval_map_snapshots) - 1
    shift_duration = (end_shift/cfg.interval_map_snapshots) - (start_shift/cfg.interval_map_snapshots)
    second_state_duration = (cfg.numberOfTimeSteps/cfg.interval_map_snapshots) - (end_shift/cfg.interval_map_snapshots)

    num1, num2, num3 = steps_ba_shift, steps_in_shift, steps_ba_shift
    if first_state_duration < num1:
        num1 = math.floor(first_state_duration)
        print("Timesteps before transition set to:", num1)
    if shift_duration < num2:
        num2 = math.ceil(shift_duration)
        print("Timesteps during transition set to:", num2)
    if second_state_duration < num3:
        num3 = math.floor(second_state_duration)
        print("Timesteps after transition set to:", num3)

    before_shift = np.linspace(1,
                               math.floor(start_shift/cfg.interval_map_snapshots) - 1,
                               num1, dtype='int') * cfg.interval_map_snapshots
    during_shift = np.linspace(math.floor(start_shift/cfg.interval_map_snapshots),
                               math.ceil(end_shift/cfg.interval_map_snapshots),
                               num2, dtype='int') * cfg.interval_map_snapshots
    after_shift = np.linspace(math.ceil(end_shift/cfg.interval_map_snapshots) + 1,
                               cfg.numberOfTimeSteps/cfg.interval_map_snapshots,
                               num3, dtype='int') * cfg.interval_map_snapshots
    timesteps = np.concatenate((before_shift, during_shift, after_shift))
    return timesteps

def file_name_str(name, timestep):
    file_name_str = ["0"]*11
    name_list = list(name)
    timestep_list = reversed(str(timestep))
    for k, letter in enumerate(name_list):
        file_name_str[k] = letter
    for k, number in enumerate(timestep_list):
        file_name_str[-(k+1)] = number
    file_name_str.insert(-3, '.')
    file_name_str = [''.join(file_name_str)]
    return file_name_str

timesteps1748_4249 = timesteps(1748,4249)
for step in timesteps1748_4249:
    print(file_name_str("demM", step))

### INPUTS & CALCULATIONS ###

# variables = ["bioM", "dem"]
#
# # Load bio mass and soil thickness timeseries
# timeseries_bioM = file_loader("bioM", path="./1/", datatype="map") # TODO timeseries only, not from maps
# mean_time_series_bioM = mean_time_series(timeseries_bioM)

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

# mthread = threading.Thread(target=user_inputs)
# mthread.start()
# plot()
# shift_data = my_queue.get()
# mthread.join()
#
# ## end of threading
#
# start_shift = shift_data[0]
# end_shift = shift_data[1]
#
# # Calculations
# timesteps_ = timesteps(int(start_shift), int(end_shift))
# timestep_str_ = timestep_str(timesteps_)
#
# # Loop over variables
# #for variable in variables: TODO - Assign snapshots of different state variables to different arrays, or calculate EWS and save results
# print("Started loading snapshots, this might take a while...")
# snapshots_bioM = snapshot_loader("bioM", timestep_str_, path='./1/', datatype='map')
# print(snapshots_bioM)
