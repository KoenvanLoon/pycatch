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

## Variables ##
variables = ews_sv.variables_weekly

names = []
for variable in variables:
    names.append([f'{variable.full_name} as {variable.name}'])

## Statistical EWS ##
ews_temporal_signals = {'std': "standard deviation", 'var': "variance", 'cv': "coefficient of variation",
                        'skw': "skewness", 'krt': "kurtosis", 'dfa': "detrended fluctuation analysis",
                        'acr': "autocorrelation", 'rr': "return rate", 'coh': "conditional heteroskedasticity"}
ews_spatial_signals = {'std': "standard deviation", 'var': "variance", 'skw': "skewness", 'krt': "kurtosis",
                       'mI': "Moran's I"}

## Null models ##
generate_dummy_datasets = True
save_detrended_data = False  # Temporal only, and only relevant when detrending != None
method_1 = True
method_2 = True
method_3 = True
nr_generated_datasets = 100


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


def histogram_plot(variable, statistic, method, tau_values, tau_original=np.NaN, chance_cor=np.NaN, path='./1/'):
    bins = np.linspace(-1, 1, num=41)
    histogram = np.histogram(tau_values, bins)
    print("Histogram values ([values],[bins]):", histogram)

    plt.hist(tau_values, bins, ec='black')
    if tau_original != np.NaN:
        plt.axvline(x=tau_original, color='r')
        if variable.temporal:
            plt.title(f"Probability distribution of the trend statistic ({ews_temporal_signals[statistic]}) under {method} of {variable.full_name}")
        elif variable.spatial:
            plt.title(f"Probability distribution of the trend statistic ({ews_spatial_signals[statistic]}) under {method} for  of {variable.full_name}")
        plt.xlabel("K \u03C4")
        plt.ylabel("Frequency")
        plt.text(0.7, 0, f" \u03C4 = {round(tau_original, 3)} \n p = {round(chance_cor, 3)} \n")
    plt.tight_layout()
    plt.show()

    print("Save the plot as a .pdf? [Y/n]")
    save_plot = input()
    if save_plot == 'Y' or save_plot == 'y':
        plt.savefig(path + f"{variable.name}_{statistic}_{method}_histogram.pdf", format="pdf")
    plt.close()


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


# for variable in variables:
#     tau, p = kendalltau_stats(variable, 'var')
#     print(variable.name, tau, p)
#
# for variable in variables:
#     # tau, _ = kendalltau_stats(variable, 'AR1', path='./1/dtr0000/')
#     tau, _ = kendalltau_stats(variable, 'AR1')
#     tau_values, p = kendalltau_stats_dummy(variable, 'AR1')
#     chance_cor = chance_value(tau_values, tau)
#     print(chance_cor)
#     histogram_plot(variable, statistic='var', method='m1g', tau_values=tau_values, tau_original=tau, chance_cor=chance_cor)
#     print(variable.name, tau, p, np.nansum(np.array(p) > 0.05))

def kendall_tau_valhist(path='./1/'):
    print("Variables present in the current run are:")
    for name in names:
        print(name)

    print("Enter the short name for the state variable:")
    state_variable_input = input()
    state_variable = [variable for variable in variables if variable.name == state_variable_input][0]

    if state_variable.temporal:
        print("EW signals present are:", ews_temporal_signals)
    elif state_variable.spatial:
        print("EW signals present are:", ews_spatial_signals)

    print("Enter the short name for the EW signal:")
    summary_statistic = input()
    tau_m, p_m = kendalltau_stats(state_variable=state_variable, sum_stat=summary_statistic)

    print("Do you want to make a histogram comparison graph? [Y/n]")
    histogram = input()
    if histogram == 'Y' or histogram == 'y':
        print("For which null model do you want to test? [m1g, m2g, m3g]")
        method = input()
        tau_values, p_values = kendalltau_stats_dummy(state_variable=state_variable, sum_stat=summary_statistic)  # TODO - Add method
        chance_cor = chance_value(tau_values, tau_m)

        histogram_plot(variable=state_variable, statistic=summary_statistic, method=method, tau_values=tau_values, tau_original=tau_m, chance_cor=chance_cor, path=path)
    else:
        print("Print Kendall tau and p values? [Y/n]")
        print_on = input()
        if print_on == 'Y' or print_on == 'y':
            print(tau_m, p_m)

    print("Would you like to make another plot? [Y/n]")
    answer = input()
    if answer == 'Y' or answer == 'y':
        kendall_tau_valhist(path=path)
    else:
        print("Terminated histogram maker. Goodbye.")


#kendall_tau_valhist(path='./1/')

# TODO - Variables sorted as in weekly_plots, move other things to cfg
# TODO - Change method names


def window(timeseries, window_size, window_overlap):
    actual_window_overlap = window_size - window_overlap
    sh = (timeseries.size - window_size + 1, window_size)
    st = timeseries.strides * 2
    if window_overlap != 0:
        view = np.lib.stride_tricks.as_strided(timeseries, strides=st, shape=sh)[0::actual_window_overlap]
    elif window_overlap == 0:
        view = np.lib.stride_tricks.as_strided(timeseries, strides=st, shape=sh)[0::window_size]
    return view



def test_windowsize(path='./1/'):
    fname = ews.file_name_str('bioA', cfg.number_of_timesteps_weekly)
    fpath = os.path.join(path + fname)
    biomass_timeseries = np.loadtxt(fpath + '.numpy.txt')

    # window_sizes = ews.divisor_generator(500, cfg.number_of_timesteps_weekly)[:-1]
    # window_overlaps = [0, 10]
    window_sizes = np.arange(1000, cfg.number_of_timesteps_weekly // 2 + 1, 1000)
    window_overlaps = np.arange(0, 1000, 10)
    x, y = np.meshgrid(window_overlaps, window_sizes)

    tau_arr = np.zeros((len(window_sizes), len(window_overlaps)))
    p_arr = np.zeros((len(window_sizes), len(window_overlaps)))

    for i, window_size in enumerate(window_sizes):
        print(f"Moved to windowsize {window_size}")
        for j, window_overlap in enumerate(window_overlaps):
            # stack_of_windows = np.array([biomass_timeseries[i:i + window_size] for i in range(0, len(biomass_timeseries), window_size - window_overlap)])

            stack_of_windows = window(biomass_timeseries, window_size, window_overlap)

            mean = np.nanmean(stack_of_windows, axis=1)
            stat = ews.temporal_AR1(stack_of_windows)

            tau, p = scipy.stats.kendalltau(stat, mean, nan_policy='propagate')
            tau_arr[i, j] = tau
            p_arr[i, j] = p

    z_array = [tau_arr, p_arr]
    fig, axs = plt.subplots(ncols=2)
    for ax, z in zip(axs, z_array):
        # Contour plot
        cs = ax.contourf(x, y, z)
        ax.contour(cs, colors='k')
        if z.all() == tau_arr.all():
            ax.set_title('Biomass variance windowtest \u03C4-value')
            min_value = np.min(z)
            local_min_index = np.where(z == min_value)
            min_x, min_y = x[local_min_index[0], local_min_index[1]], y[local_min_index[0], local_min_index[1]]
            ax.plot(min_x, min_y, color='blue', marker="v", zorder=10, markersize=10, clip_on=False)
        elif z.all() == p_arr.all():
            ax.set_title('Biomass variance windowtest p-value')
        ax.set_ylabel('Window size')
        ax.set_xlabel('Window overlap')
        ax.grid(c='k', alpha=0.3)
        fig.colorbar(cs, ax=ax)

        # Max & Min value
        max_value = np.max(z)
        local_max_index = np.where(z == max_value)
        max_x, max_y = x[local_max_index[0], local_max_index[1]], y[local_max_index[0], local_max_index[1]]
        ax.plot(max_x, max_y, color='red', marker="v", zorder=10, markersize=10, clip_on=False)

    plt.show()


#test_windowsize()
kendall_tau_valhist()

# a = np.arange(0, 100, 1)
# print(window(a, 9, 5))
