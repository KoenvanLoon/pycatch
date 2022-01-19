import numpy as np
import math
from pcraster import *
import matplotlib.pyplot as plt


# # Configuration file inputs
# #nr_of_weeks = 52
# #hours_per_week = 7 * 24
# nr_of_weeks = 10  # For testing
# hours_per_week = 10  # For testing
# probabilityOfRainstorm = 0.4
# rainstormDuration = 2
#
# # 1 if this week contains rain, else 0 (uniform distribution)
# rairray = np.random.rand(nr_of_weeks)
# rairray[rairray >= (1-probabilityOfRainstorm)] = 1
# rairray[rairray < (1-probabilityOfRainstorm)] = 0
# #print(rairray)
#
# # 2D array of size nr_of_weeks x hours_per_week (= number of hours in the hourly model run)
# rairray_2D = np.zeros((nr_of_weeks, hours_per_week))
#
# # First day of the week equals 1 if it rains in this week
# rairray_2D[:, 0] = rairray
# #print(rairray_2D)
#
# # Shuffles the hours in the week such that the rain does (most likely) not fall on the same day every week
# rairray_hours = np.copy(rairray_2D)
# idx = np.random.rand(*rairray_hours.shape).argsort(axis=1)
# rairray_hours = np.take_along_axis(rairray_hours, idx, axis=1)
#
# # Flatten the 2D array to a 1D array with the length of the hourly model run
# rairray_hours_flat = rairray_hours.flatten()
# #print(rairray_hours_flat)
#
# # For the set duration of the rainstorm, a 1 is added after an already given 1
# rairray_timeseries = rairray_hours_flat
# for i in range(1, rainstormDuration):
#     rairray_timeseries += np.roll(rairray_hours_flat, i)
#
# #print(rairray_timeseries)
#
# # Boolean timeseries for rainstorm
# rairray_timeseries_bool = rairray_timeseries.astype(bool)
# print(rairray_timeseries_bool)


# def time_series2time_windows(time_series, window_size=100, window_overlap=0):
#     return np.array([time_series[i:i + window_size] for i in range(0, len(time_series), window_size - window_overlap)])
#
#
# demini = pcr2numpy(readmap('./inputs_weekly/demini.map'), np.NaN)
# locations2report = pcr2numpy(readmap('./inputs_weekly/mlocs.map'), 0)
# print(locations2report)
# a = demini[locations2report.astype(bool)]
#
# # print(np.loadtxt('./1/bioL0005.000.numpy.txt'))
# # print(np.loadtxt('./1/bioA0005.000.numpy.txt'))
# print(type(np.loadtxt('./1/bioL0005.000.numpy.txt')), np.loadtxt('./1/bioL0005.000.numpy.txt').shape)
# print(type(np.loadtxt('./1/bioA0005.000.numpy.txt')), np.loadtxt('./1/bioA0005.000.numpy.txt').shape)
#
# state_variable_timeseries = np.loadtxt('./1/bioL0005.000.numpy.txt')
#
# if state_variable_timeseries.ndim == 1:
#     stack_of_windows = time_series2time_windows(state_variable_timeseries, 1000, 0)
# else:
#     stack_of_windows = [0] * state_variable_timeseries.shape[0]
#     for k, timeseries in enumerate(state_variable_timeseries):
#         stack_of_windows[k] = time_series2time_windows(timeseries, 1000, 0)
#     stack_x, stack_y, stack_z = np.asarray(stack_of_windows).shape
#     stack_of_windows = np.asarray(stack_of_windows).reshape(-1, stack_z)
#
# print(stack_of_windows)
# print(stack_of_windows.shape)
# print(stack_x, stack_y, int(stack_z / (100 - 0)))
# print(state_variable_timeseries.ndim)
#
# if state_variable_timeseries.ndim > 1:
#     stack_of_windows = stack_of_windows.reshape(stack_x, stack_y, stack_z)
#
# print(stack_of_windows)
# print(stack_of_windows.shape)

# x = np.linspace(0, 4, 5)
# y = np.asarray([[0.32,1.25,2.36,3.36,3.52],[0.2,1.5,2.6,2.3,1.5]])
# print(x, y)
#
# plt.plot(x, y.T)
# plt.show()
#
# import numpy as np
#
# fpath = 'C:/Users/koenv/Downloads/IndexFund/'
# index_fund_data = np.genfromtxt(fpath + 'VFINX.csv', delimiter=',')
# index_fund_closing_value = index_fund_data.T[4][1:10501]
#
# np.savetxt(fpath + 'INDF0010.500.numpy.txt', index_fund_closing_value)


# import numpy as np
#
# def dfa_propagator(alpha, c_guess=0.1):
#     # 0.91 * (c ** 3) - 0.37 * (c ** 2) + 0.49 * c + c - alpha = 0
#     # a*x**3 + b*x**2 + c*x + d = 0.
#
#     c1 = c_guess
#     count = 0
#     while count < 10:
#
#         if 0 < c1 <= 0.936:
#             a = 0.91
#             b = - 0.37
#             c = 0.49
#             d = 0.52 - alpha
#
#         elif 0.936 < c1 <= 0.967:
#             a = 0
#             b = -12.38
#             c = 25.14
#             d = 11.28 - alpha
#
#         # elif 0.967 < c1 < 1:
#         else:
#             a = 0
#             b = 0
#             c = 0.72
#             d = 0.75 - alpha
#
#         q = (3*a*c - b**2) / (9 * a**2)
#         r = (9*a*b*c - 27*a**2*d - 2*b**3) / (54*a**3)
#
#         print("q = ", q)
#         print("r = ", r)
#
#         delta = q**3 + r**2
#
#         print("delta = ", delta)
#
#         # here delta is less than zero so we use the second set of equations from the article:
#
#         rho = (-q**3)**0.5
#
#         # For c1 the imaginary part is unimportant since it cancels out
#         theta = np.arccos(r/rho)
#         s_real = rho**(1./3.) * np.cos(theta/3)
#         t_real = rho**(1./3.) * np.cos(-theta/3)
#
#         print("s [real] = ", s_real)
#         print("t [real] = ", t_real)
#
#         c1 = s_real + t_real - b / (3. * a)
#
#         print("c = ", c1)
#
#         print("should be zero: ", a*c1**3+b*c1**2+c*c1+d)
#
#         count += 1
#
#     return c1
#
#
# dfa_propagator(alpha=1)

import numpy as np

def dfa_propagator(alpha, c_guess=0.5):
    # 0.91 * (c ** 3) - 0.37 * (c ** 2) + 0.49 * c + c - alpha = 0
    # a*x**3 + b*x**2 + c*x + d = 0.

    x1 = c_guess
    count = 0
    while count < 5:

        if 0 < x1 <= 0.936:
            a = 0.91
            b = - 0.37
            c = 0.49
            d = 0.52 - alpha

        elif 0.936 < x1 <= 0.967:
            a = 0
            b = -12.38
            c = 25.14
            d = 11.28 - alpha

        # elif 0.967 < x1 < 1:
        else:
            a = 0
            b = 0
            c = 0.72
            d = 0.75 - alpha

        # Constructs the polynomial a*x**3 + b*x**2 + c*x + d
        poly = np.poly1d([a, b, c, d])
        roots = np.roots(poly)
        print(roots)
        x1 = roots[-1].real

        count += 1

    return x1


# dfap = dfa_propagator(0.5)
# print(dfap)

# a = np.linspace(0, 0.9, 9)
# b = np.linspace(0.1, 1, 9)
# c = np.linspace(0, 1, 10)
#
# print(a)
# print(b)
# print(c)

# import scipy.stats
#
# chi = scipy.stats.chi2.ppf((1 - 0.1), df=1)
# print(chi)
