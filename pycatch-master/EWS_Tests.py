import numpy as np
import os
import scipy.stats
import matplotlib.pyplot as plt

import EWSPy as ews
import EWS_main_configuration as cfg
import EWS_StateVariables as ews_sv

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

# variable_window_size = 100
# start_index = int(cfg.rel_start_grazing * cfg.number_of_timesteps_weekly / variable_window_size)
# end_index = int(4000 / variable_window_size)
#
# # Biomass data
# biomass_state_var = np.loadtxt('./1/bioA.t.mn.numpy.txt')
# biomass_sum_stat = np.loadtxt('./1/bioA.t.acr.numpy.txt')
#
# # Grazing rate
# grazing_mean = np.loadtxt('./1/gA.t.mn.numpy.txt')

# # M1G: Data shuffled
# biomass_mean_m1g = np.loadtxt('./1/m1g0000/bioM.s.mn.numpy.txt')
# biomass_var_m1g = np.loadtxt('./1/m1g0000/bioM.s.var.numpy.txt')
#
# # M2G: FFT
# biomass_mean_m2g = np.loadtxt('./1/m2g0000/bioM.s.mn.numpy.txt')
# biomass_var_m2g = np.loadtxt('./1/m2g0000/bioM.s.var.numpy.txt')
#
# # M3G: AR
# biomass_mean_m3g = np.loadtxt('./1/m3g0000/bioM.s.mn.numpy.txt')
# biomass_var_m3g = np.loadtxt('./1/m3g0000/bioM.s.var.numpy.txt')

# bsv_tau, bss_p = scipy.stats.kendalltau(biomass_state_var[start_index:end_index], biomass_sum_stat[start_index:end_index],
#                                         nan_policy='propagate')
#
# grb_tau, grb_p = scipy.stats.kendalltau(grazing_mean[start_index:end_index], biomass_sum_stat[start_index:end_index],
#                                         nan_policy='propagate')

# m1g_tau, m1g_p = scipy.stats.kendalltau(biomass_mean_m1g[start_index:end_index], biomass_var_m1g[start_index:end_index])
# m2g_tau, m2g_p = scipy.stats.kendalltau(biomass_mean_m2g[start_index:end_index], biomass_var_m2g[start_index:end_index])
# m3g_tau, m3g_p = scipy.stats.kendalltau(biomass_mean_m3g[start_index:end_index], biomass_var_m3g[start_index:end_index])

# print('Biomass state variable & summary statistic', bsv_tau, bss_p)
# print('Grazing rate & biomass summary statistic', grb_tau, grb_p)
# print('m1', m1g_tau, m1g_p)
# print('m2', m2g_tau, m2g_p)
# print('m3', m3g_tau, m3g_p)

variables = ews_sv.variables_weekly

generate_dummy_datasets = True
save_detrended_data = True  # Temporal only, and only relevant when detrending != None
method_1 = True
method_2 = True
method_3 = True
nr_generated_datasets = 1000


def kendalltau_stats(state_variable, sum_stat, comp2='Same', path='./1/'):
    if state_variable.temporal:
        dim = '.t.'
    elif state_variable.spatial:
        dim = '.s.'

    tau, p = np.NaN, np.NaN
    fdict = os.path.join(path + state_variable.name + dim)
    if os.path.exists(fdict + sum_stat + '.numpy.txt'):
        X = np.loadtxt(fdict + sum_stat + '.numpy.txt')

        if comp2 == 'Same':  # Dakos et al, 2008
            Y = np.loadtxt(fdict + 'mn.numpy.txt')
        elif comp2 == 'Forcing':  # Dakos et al, 2011 - Does not work if window sizes are different for the forcing & SV
            Y = np.loadtxt('./1/gA.t.mn.numpy.txt')

        tau, p = scipy.stats.kendalltau(X, Y, nan_policy='propagate')

    return tau, p


def kendalltau_stats_dummy(state_variable, sum_stat, comp2='Same', path='./1/'):
    if state_variable.temporal:
        dim = '.t.'
    elif state_variable.spatial:
        dim = '.s.'

    generated_number_length = 4
    if len(str(nr_generated_datasets)) > 4:
        generated_number_length = len(str(nr_generated_datasets))

    method = 'm2g'

    taurray = [np.NaN] * nr_generated_datasets
    parray = [np.NaN] * nr_generated_datasets
    for realization in range(nr_generated_datasets):
        fdict = os.path.join(path + method + str(realization).zfill(generated_number_length) + '/'
                             + state_variable.name + dim)
        if os.path.exists(fdict + sum_stat + '.numpy.txt'):
            X = np.loadtxt(fdict + sum_stat + '.numpy.txt')

            if comp2 == 'Same':  # Dakos et al, 2008
                Y = np.loadtxt(fdict + 'mn.numpy.txt')
            elif comp2 == 'Forcing':  # Dakos et al, 2011 - Does not work if window sizes are different for the forcing & SV, detrending == Gaus
                Y = np.loadtxt('./1/gA.t.mn.numpy.txt')

            taurray[realization], parray[realization] = scipy.stats.kendalltau(X, Y, nan_policy='propagate')

    return taurray, parray


def histogram_plot(variable, statistic, method, tau_values, tau_original=np.NaN, chance_cor=np.NaN):
    bins = np.linspace(-1, 1, num=41)
    histogram = np.histogram(tau_values, bins)
    print(histogram)

    _ = plt.hist(tau_values, bins, ec='black')
    if tau_original != np.NaN:
        plt.axvline(x=tau_original, color='r')
        plt.title(f"Probability distribution of the trend statistic under {method} for {statistic} of {variable.full_name}")
        plt.xlabel("K \u03C4")
        plt.ylabel("Frequency")
        plt.text(0.7, 0, f" \u03C4 = {round(tau_original, 3)} \n p = {round(chance_cor, 3)} \n")
    plt.show()
    return None

def chance_value(values, value, sign='None'):
    if sign == 'None':
        if value > 0.:
            sign = '>'
        elif value < 0.:
            sign = '<'
        else:
            print(f"Tau value of {value} encountered while sign is not specified. Defaulted to 'tau values' > 'tau value'.")
            sign = '>'

    if sign == '>':
        p = np.nansum(np.array(values) > value) / len(values)
    elif sign == '<':
        p = np.nansum(np.array(values) < value) / len(values)
    else:
        print(f"{sign} is not supported")

    return p


for variable in variables:
    tau, p = kendalltau_stats(variable, 'AR1')
    print(variable.name, tau, p)

for variable in variables:
    # tau, _ = kendalltau_stats(variable, 'AR1', path='./1/dtr0000/')
    tau, _ = kendalltau_stats(variable, 'AR1')
    tau_values, p = kendalltau_stats_dummy(variable, 'AR1')
    chance_cor = chance_value(tau_values, tau)
    print(chance_cor)
    histogram_plot(variable, statistic='var', method='m1g', tau_values=tau_values, tau_original=tau, chance_cor=chance_cor)
    print(variable.name, tau, p, np.nansum(np.array(p) > 0.05))
