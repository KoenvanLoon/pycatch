import numpy as np
import os
import scipy.stats

import EWS_main_configuration as cfg

# variable_window_size = 1000
#
# step_size = int(cfg.number_of_timesteps_weekly * cfg.rel_start_grazing)
# timeseries_index = np.arange(0, cfg.number_of_timesteps_weekly, step_size, dtype=int)
# range_index = timeseries_index // variable_window_size
#
# biomass_mean = np.loadtxt('./1/dtr0000/bioA.t.mn.numpy.txt')
# biomass_var = np.loadtxt('./1/dtr0000/bioA.t.var.numpy.txt')
#
# biomass_mean_sliced = np.array([biomass_mean[i:i + step_size] for i in range_index])
# biomass_var_sliced = np.array([biomass_var[i:i + step_size] for i in range_index])
#
# print(scipy.stats.kendalltau(biomass_mean, biomass_var))
#
# for i in range(len(biomass_mean_sliced)):
#     kt_coef, kt_p = scipy.stats.kendalltau(biomass_mean_sliced[i], biomass_var_sliced[i])
#     sm_coef, sm_p = scipy.stats.spearmanr(biomass_mean_sliced[i], biomass_var_sliced[i])
#     print("Timestep ", timeseries_index[i], "to", timeseries_index[i] + step_size, ":", '\n',
#           "Correlation:", kt_coef, sm_coef, '\n',
#           "p-value:", kt_p, sm_p)

variable_window_size = 200
start_index = int(cfg.rel_start_grazing * cfg.number_of_timesteps_weekly / variable_window_size)
end_index = int(4000 / variable_window_size)

# Detrended data
biomass_mean_dtr = np.loadtxt('./1/dtr0000/bioA.t.mn.numpy.txt')
biomass_var_dtr = np.loadtxt('./1/dtr0000/bioA.t.var.numpy.txt')

# Trended (?) data
biomass_mean_tr = np.loadtxt('./1/dtr0000/bioA.t.mn.numpy.txt')

# M1G: Data shuffled
biomass_mean_m1g = np.loadtxt('./1/m1g0000/bioA.t.mn.numpy.txt')
biomass_var_m1g = np.loadtxt('./1/m1g0000/bioA.t.var.numpy.txt')

# M2G: FFT
biomass_mean_m2g = np.loadtxt('./1/m2g0000/bioA.t.mn.numpy.txt')
biomass_var_m2g = np.loadtxt('./1/m2g0000/bioA.t.var.numpy.txt')

# M3G: AR
biomass_mean_m3g = np.loadtxt('./1/m3g0000/bioA.t.mn.numpy.txt')
biomass_var_m3g = np.loadtxt('./1/m3g0000/bioA.t.var.numpy.txt')

dtr_tau, dtr_p = scipy.stats.kendalltau(biomass_mean_dtr[start_index:end_index], biomass_var_dtr[start_index:end_index])
# tr_tau, tr_p = scipy.stats.kendalltau(biomass_mean_tr[start_index:end_index], biomass_var_tr[start_index:end_index])
m1g_tau, m1g_p = scipy.stats.kendalltau(biomass_mean_m1g[start_index:end_index], biomass_var_m1g[start_index:end_index])
m2g_tau, m2g_p = scipy.stats.kendalltau(biomass_mean_m2g[start_index:end_index], biomass_var_m2g[start_index:end_index])
m3g_tau, m3g_p = scipy.stats.kendalltau(biomass_mean_m3g[start_index:end_index], biomass_var_m3g[start_index:end_index])

print('detrended', dtr_tau, dtr_p)
# print('trend', tr_tau, tr_p)
print('m1', m1g_tau, m1g_p)
print('m2', m2g_tau, m2g_p)
print('m3', m3g_tau, m3g_p)
