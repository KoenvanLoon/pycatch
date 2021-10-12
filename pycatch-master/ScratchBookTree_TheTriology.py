import EWSPy as ews
import numpy as np
import os
path = 'C:/Users/koenv/Bureaublad/Thesis/pycatch-master/pycatch-master/1'

stack_of_maps_as_list = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path,i)) and 'rQ' in i]
print(stack_of_maps_as_list)

time_window = 100
stacks_of_maps = ews.hap_klare_hapjes(stack_of_maps_as_list, time_window)
print(stacks_of_maps)
a = np.loadtxt(stacks_of_maps)
print(a)

# time_series_stack = ews.stack_of_maps2max_time_series(stacks_of_maps, time_window)
#
# tstorage = [0.0] * time_series_stack.shape[0]
#
# for k, time_series in enumerate(time_series_stack):
#     tstorage[k] += ews.temporal_autocorrelation(time_series)
# print("temporal autocorrelation lag-1:", tstorage)
