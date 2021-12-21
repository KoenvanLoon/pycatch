import numpy as np

nr_of_weeks = 10
hours_per_week = 7 #* 24
probabilityOfRainstorm = 0.4
rainstormDuration = 2

rairray = np.random.rand(nr_of_weeks)

rairray[rairray >= (1-probabilityOfRainstorm)] = 1
rairray[rairray < (1-probabilityOfRainstorm)] = 0

rairray_bool = rairray.astype(bool)

print(rairray_bool)

rairray_hours = np.zeros((nr_of_weeks, hours_per_week))

rairray_hours[:, 0] = rairray
print(rairray_hours, '\n')

rairray_hours_shuffle = np.copy(rairray_hours)

idx = np.random.rand(*rairray_hours_shuffle.shape).argsort(axis=1)
rairray_hours_shuffle = np.take_along_axis(rairray_hours_shuffle, idx, axis=1)

rairray_hours_shuffle_flat = rairray_hours_shuffle.flatten()
rairray_hours_shuffle_flat_rolled = np.roll(rairray_hours_shuffle_flat, 1)

rairray_timeseries = rairray_hours_shuffle_flat + rairray_hours_shuffle_flat_rolled

print(rairray_hours_shuffle_flat)
print(rairray_hours_shuffle_flat_rolled)
print(rairray_timeseries)
