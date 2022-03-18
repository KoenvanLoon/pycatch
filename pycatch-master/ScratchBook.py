import numpy as np
import math
import os
from cycler import cycler
import matplotlib.pyplot as plt
from datetime import datetime

import EWS_configuration as cfg
import EWSPy as ews
import EWS_StateVariables as ews_sv

## State variables for EWS ##
# State variables present in EWS_StateVariables can be added through EWS_main_configuration.py
variables = ews_sv.variables_weekly
# variables = ews_sv.variables_hourly

timeseries = ews_sv.variables_weekly

names = []
for variable in variables:
    names.append([f'{variable.full_name} as {variable.name}'])

names_timeseries = []
for ts in timeseries:
    names_timeseries.append([f'{ts.full_name} as {ts.name}'])

## Number of timesteps over which EWS can be calculated ##
# This number can be different for the weekly and hourly model
number_of_timesteps = cfg.number_of_timesteps_weekly
# number_of_timesteps = cfg.number_of_timesteps_hourly

## Statistical EWS ##
ews_temporal_signals = {'mn': "mean", 'std': "standard deviation", 'var': "variance",
                        'cv': "coefficient of variation", 'skw': "skewness", 'krt': "kurtosis",
                        'dfa': "detrended fluctuation analysis", 'acr': "autocorrelation", 'AR1': "AR1",
                        'rr': "return rate", 'coh': "conditional heteroskedasticity", 'timeseries': "timeseries",
                        'gauss': "gauss", 'linear': "linear"}
ews_spatial_signals = {'mn': "mean", 'std': "standard deviation", 'var': "variance", 'skw': "skewness",
                       'krt': "kurtosis", 'mI': "Moran's I"}


def user_input_weekly_hourly_coupled(path='./1/'):
    timeseries_list = []
    variables_list = []
    signals_list = []

    print("Timeseries present in the current run are:")
    for name in names_timeseries:
        print(name)

    nr_of_variables = 0
    cont = True
    while cont:
        nr_of_variables += 1

        print("Enter the short name for the timeseries:")
        timeseries_input = input()

        ts_weekly = [ts for ts in timeseries if ts.name == timeseries_input][0]
        if ts_weekly.temporal:
            timeseries_list.append(ts_weekly)
        elif ts_weekly.spatial:
            print("Spatial data can not be used for timeseries.")

        if nr_of_variables == 9:
            cont = False
        elif nr_of_variables < 9:
            print("Include another timeseries? [Y/n]")
            another_input = input()
            if another_input == 'Y' or another_input == 'y':
                cont = True
            else:
                cont = False

    if nr_of_variables < 9:
        cont = True

        print("Variables present in the current run are:")
        for name in names:
            print(name)

    while cont:
        nr_of_variables += 1

        print(f"Enter the short name for state variable {nr_of_variables}")
        variable_input = input()
        variable_name = [var for var in variables if var.name == variable_input][0]
        variables_list.append(variable_name)

        if variable_name.temporal:
            print("EW signals present are:", ews_temporal_signals)
            dim = '.t.'
        elif variable_name.spatial:
            print("EW signals present are:", ews_spatial_signals)
            dim = '.s.'

        print(f"Enter the short name for the signal for variable {nr_of_variables}")
        signal_input = input()
        signals_list.append(signal_input)

        if nr_of_variables == 9:
            cont = False
        elif nr_of_variables < 9:
            print("Include another variable? [Y/n]")
            another_input = input()
            if another_input == 'Y' or another_input == 'y':
                cont = True
            else:
                cont = False

    print("Add a legend to the plot? [Y/n]")
    legend_input = input()
    if legend_input == 'Y' or legend_input == 'y':
        legend = True
    else:
        legend = False

    print("Save the plot as a .pdf? [Y/n]")
    save_plot = input()
    if save_plot == 'Y' or save_plot == 'y':
        save = True
    else:
        save = False

    print("Show the plot when finished? [Y/n]")
    print("Note that the program is still running if the plot stays open.")
    show_plot = input()
    if show_plot == 'Y' or show_plot == 'y':
        show = True
    else:
        show = False

    plot_maker_weekly_hourly_coupled(timeseries=timeseries_list, variables=variables_list, signals=signals_list,
                                     path=path, legend=legend, save=save, show=show)


def plot_maker_weekly_hourly_coupled(timeseries, variables, signals, path='./1/', trendline_on=False, numbers_on=False,
                                     legend=False, save=False, show=False):
    fig, ax1 = plt.subplots()

    snapshot_timeseries = np.loadtxt('./snapshot_times.numpy.txt')
    timeseries_x_axis = np.arange(0, cfg.number_of_timesteps_weekly, 1)

    nr_of_variables = len(timeseries) + len(variables)
    axes = [ax1]
    plots = []
    offset = 60
    for i in np.arange(nr_of_variables):
        if nr_of_variables > i + 1:
            ax = ax1.twinx()
            ax.spines["right"].set_position(("outward", offset * i))
            axes.append(ax)

    # Grid
    ax1.grid(which='minor', linestyle=':', alpha=0.2)
    ax1.grid(which='major', linestyle='--', alpha=0.5)

    # X axis label (uses only 1 for the whole plot)
    ax1.set_xlabel("time (weeks)")

    # Linestyles
    linestyles = [
        (0, ()),  # solid
        (0, (1, 1)),  # dotted
        (0, (5, 5)),  # dashed
        (0, (3, 1, 1, 1)),  # densely dashdotted
        (0, (3, 1, 1, 1, 1, 1)),  # densely dashdotdotted
        (0, (1, 1)),  # dotted
        (0, (5, 5)),  # dashed
        (0, (3, 1, 1, 1)),  # densely dashdotted
        (0, (3, 1, 1, 1, 1, 1))  # densely dashdotdotted
    ]

    # Colours list
    colours_list = [plt.cm.tab10(i) for i in np.linspace(0, 1, 9)]

    for i in np.arange(nr_of_variables):
        ax = axes[i]

        # Cycle colours and linestyles
        even = i % 2  # True if even number, False if odd number
        idx = math.floor(i/2)
        if even:
            colours = np.concatenate((np.asarray(colours_list)[idx:], np.asarray(colours_list)[:idx]))
        if not even:
            colours = np.concatenate((np.asarray(colours_list)[idx:], np.asarray(colours_list)[:idx]))
            colours = colours[::-1]
        ax.set_prop_cycle(cycler(color=colours, linestyle=linestyles))

        # Ticks
        ax.minorticks_on()
        ax.tick_params(axis='y', which='both', colors=colours[i])

        if i <= len(timeseries) - 1:
            fname_ts_weekly = ews.file_name_str(timeseries[i].name, cfg.number_of_timesteps_weekly)
            fpath_ts_weekly = f"./inputs_from_weekly/{fname_ts_weekly}"
            timeseries_y_axis = np.loadtxt(fpath_ts_weekly + '.numpy.txt')
            plot = ax.plot(timeseries_x_axis, timeseries_y_axis, label=f"Timeseries of {timeseries[i].full_name}", color=colours[i])

            for p in plot:
                plots.append(p)

        else:
            snapshot_y_axis = []
            snapshot_x_axis = snapshot_timeseries

            if variables[i - len(timeseries)].temporal:
                dim = '.t.'
                slabel = f'{variables[i - len(timeseries)].full_name} {ews_temporal_signals[signals[i - len(timeseries)]]}'
                ax.set_ylabel(f'{ews_temporal_signals[signals[i - len(timeseries)]]}', c=plt.cm.tab10(i))
            elif variables[i - len(timeseries)].spatial:
                dim = '.s.'
                slabel = f'{variables[i - len(timeseries)].full_name} {ews_spatial_signals[signals[i - len(timeseries)]]}'
                ax.set_ylabel(f'{ews_spatial_signals[signals[i - len(timeseries)]]}', c=plt.cm.tab10(i))

            fname_signal_hourly = variables[i - len(timeseries)].name + dim + signals[i - len(timeseries)] + '.numpy.txt'
            print(fname_signal_hourly)

            for nr, _ in enumerate(snapshot_timeseries):
                fpath_signal_hourly = os.path.join("./h" + str(nr).zfill(2) + "/" + fname_signal_hourly)
                statistic = np.loadtxt(fpath_signal_hourly)
                if statistic.ndim == 0:
                    snapshot_y_axis.append(statistic)
                elif statistic.ndim == 1:
                    snapshot_y_axis.append(statistic[-1])

            plot = ax.scatter(snapshot_x_axis, snapshot_y_axis, label=slabel, c=plt.cm.tab10(i))

            if trendline_on:
                z = np.polyfit(snapshot_x_axis, snapshot_y_axis, 1)
                p = np.poly1d(z)
                ax.plot(snapshot_x_axis, p(snapshot_x_axis), "r--")

            if numbers_on:
                for i, nr in enumerate(snapshot_timeseries):
                    ax.annotate(int(i), (snapshot_x_axis[i], snapshot_y_axis[i]))

            plots.append(plot)

    # Legend
    if legend:
        labs = [p.get_label() for p in plots]
        ax1.legend(plots, labs)

    # Saving file
    fig.tight_layout()
    if save:
        format = "pdf"  # either "pdf" or "png"
        timestamp = datetime.now()

        dir_name = './plots/'
        if os.path.isdir(dir_name) == False:
            os.makedirs(dir_name)

        if len(variables) > 1:
            fig.savefig(dir_name + f"{variables[0].name}_{signals[0]}_{signals[1]}_with_{nr_of_variables}_variables"
                        f"_{timestamp.hour}.{timestamp.minute}.{timestamp.second}.{format}", dpi=300, format=format)
        else:
            fig.savefig(dir_name + f"{variables[0].name}_{signals[0]}_{timestamp.hour}.{timestamp.minute}."
                        f"{timestamp.second}.{format}", dpi=300, format=format)


    # Showing plot
    if show:
        plt.draw()
        plt.show()
    elif not show:
        plt.close()


def plot_maker(variables_input, signals, path='/1/', legend=False, save=False, show=False):
    fig, ax1 = plt.subplots()

    nr_of_variables = len(variables_input)
    axes = [ax1]
    plots = []
    offset = 60
    for i in range(10):
        if nr_of_variables > i + 1:
            ax = ax1.twinx()
            ax.spines["right"].set_position(("outward", offset * i))
            axes.append(ax)

    # Grid
    ax1.grid(which='minor', linestyle=':', alpha=0.2)
    ax1.grid(which='major', linestyle='--', alpha=0.5)

    # X axis label (uses only 1 for the whole plot)
    ax1.set_xlabel("time (weeks)")

    # Linestyles
    linestyles = [
        (0, ()),  # solid
        (0, (1, 1)),  # dotted
        (0, (5, 5)),  # dashed
        (0, (3, 1, 1, 1)),  # densely dashdotted
        (0, (3, 1, 1, 1, 1, 1)),  # densely dashdotdotted
        (0, (1, 1)),  # dotted
        (0, (5, 5)),  # dashed
        (0, (3, 1, 1, 1)),  # densely dashdotted
        (0, (3, 1, 1, 1, 1, 1))  # densely dashdotdotted
    ]

    # Colours list
    colours_list = [plt.cm.tab10(i) for i in np.linspace(0, 1, 9)]

    for i in np.arange(nr_of_variables):
        ax = axes[i]
        variable = variables_input[i]
        signal = signals[i]

        # Cycle colours and linestyles
        even = i % 2  # True if even number, False if odd number
        idx = math.floor(i/2)
        if even:
            colours = np.concatenate((np.asarray(colours_list)[idx:], np.asarray(colours_list)[:idx]))
        if not even:
            colours = np.concatenate((np.asarray(colours_list)[idx:], np.asarray(colours_list)[:idx]))
            colours = colours[::-1]
        ax.set_prop_cycle(cycler(color=colours, linestyle=linestyles))

        # Ticks
        ax.minorticks_on()
        ax.tick_params(axis='y', which='both', colors=colours[i])

        # Dimension and x axis
        if variable.spatial:
            dim = '.s.'
            start = cfg.interval_map_snapshots
            step = start
        elif variable.temporal:
            dim = '.t.'
            start = variable.window_size
            step = variable.window_size - variable.window_overlap

        if cfg.cutoff:
            x_axis = np.arange(start, cfg.cutoff_point + 1, step)
        else:
            x_axis = np.arange(start, number_of_timesteps + 1, step)

        # Signal
        if signal == "timeseries":
            fname = ews.file_name_str(variable.name, number_of_timesteps)
            fpath = os.path.join(path + fname)
            timeseries_y_axis = np.loadtxt(fpath + '.numpy.txt')
            timeseries_x_axis = np.arange(0, number_of_timesteps, 1)
            plot = ax.plot(timeseries_x_axis, timeseries_y_axis, label=f"Timeseries of {variable.full_name}", color=colours[i])
            ax.set_ylabel(f"{ews_temporal_signals[signal]} ({variable.unit})", color=colours[i])

        elif signal == "gauss":
            fname = ews.file_name_str(variable.name + 'g', number_of_timesteps)
            fpath = os.path.join(path + fname)
            timeseries_y_axis = np.loadtxt(fpath + '.numpy.txt')
            timeseries_x_axis = np.arange(0, number_of_timesteps, 1)
            plot = ax.plot(timeseries_x_axis, timeseries_y_axis, label=f"Gaussian filter of {variable.full_name}", color=colours[i])
            ax.set_ylabel(f"{ews_temporal_signals[signal]} ({variable.unit})", color=colours[i])

        elif signal == "linear":
            fname = ews.file_name_str(variable.name + 'l', number_of_timesteps)
            fpath = os.path.join(path + fname)
            timeseries_y_axis = np.loadtxt(fpath + '.numpy.txt')
            timeseries_x_axis = np.arange(0, number_of_timesteps, 1)
            plot = ax.plot(timeseries_x_axis, timeseries_y_axis, label=f"Linear detrending of {variable.full_name}", color=colours[i])
            ax.set_ylabel(f"{ews_temporal_signals[signal]} ({variable.unit})", color=colours[i])

        elif signal != 'None':
            fpath = os.path.join(path + variable.name + dim + signal)
            signal_array = np.loadtxt(fpath + '.numpy.txt')
            if signal_array.ndim > 1 and variable.temporal:
                signal_array = signal_array.T
                plot = ax.plot(x_axis, signal_array)
                lines = ax.get_lines()
                ax.set_ylabel(f"{ews_temporal_signals[signal]}", color=colours[i])
                for loc, line in enumerate(lines):
                    line.set_label(f"{variable.full_name} {loc + 1} - {ews_temporal_signals[signal]}")
            elif variable.spatial:
                plot = ax.plot(x_axis, signal_array, label=f"{variable.full_name} {ews_spatial_signals[signal]}", color=colours[i])
                ax.set_ylabel(f"{ews_spatial_signals[signal]}", color=colours[i])
            elif variable.temporal:
                plot = ax.plot(x_axis, signal_array, label=f"{variable.full_name} {ews_temporal_signals[signal]}", color=colours[i])
                ax.set_ylabel(f"{ews_temporal_signals[signal]}", color=colours[i])

        for p in plot:
            plots.append(p)

    # Legend
    if legend:
        labs = [p.get_label() for p in plots]
        ax1.legend(plots, labs)

    # Saving file
    fig.tight_layout()
    if save:
        format = "pdf"  # either "pdf" or "png"
        timestamp = datetime.now()

        if nr_of_variables > 1:
            fig.savefig(path + f"{variables_input[0].name}_{signals[0]}_{signals[1]}_with_{nr_of_variables}_variables_"
                               f"{timestamp.hour}.{timestamp.minute}.{timestamp.second}.{format}", dpi=300, format=format)
        else:
            fig.savefig(path + f"{variables_input[0].name}_{signals[0]}_{timestamp.hour}.{timestamp.minute}."
                               f"{timestamp.second}.{format}", dpi=300, format=format)

    # Showing plot
    if show:
        plt.draw()
        plt.show()
    elif not show:
        plt.close()


def user_input_plotmaker(path='./1/'):
    variables_list = []
    signals_list = []

    print("Variables present in the current run are:")
    for name in names:
        print(name)

    nr_of_variables = 0
    cont = True
    while cont:
        nr_of_variables += 1

        print(f"Enter the short name for state variable {nr_of_variables}")
        variable_input = input()
        variable_name = [var for var in variables if var.name == variable_input][0]
        variables_list.append(variable_name)

        if variable_name.temporal:
            print("EW signals present are:", ews_temporal_signals)
        elif variable_name.spatial:
            print("EW signals present are:", ews_spatial_signals)
        print(f"Enter the short name for the signal for variable {nr_of_variables}")
        signal_input = input()
        signals_list.append(signal_input)

        if nr_of_variables == 9:
            cont = False
        elif nr_of_variables < 9:
            print("Include another variable? [Y/n]")
            another_input = input()
            if another_input == 'Y' or another_input == 'y':
                cont = True
            else:
                cont = False

    print("Add a legend to the plot? [Y/n]")
    legend_input = input()
    if legend_input == 'Y' or legend_input == 'y':
        legend = True
    else:
        legend = False

    print("Save the plot as a .pdf? [Y/n]")
    save_plot = input()
    if save_plot == 'Y' or save_plot == 'y':
        save = True
    else:
        save = False

    print("Show the plot when finished? [Y/n]")
    print("Note that the program is still running if the plot stays open.")
    show_plot = input()
    if show_plot == 'Y' or show_plot == 'y':
        show = True
    else:
        show = False

    plot_maker(variables_input=variables_list, signals=signals_list, path=path, legend=legend, save=save, show=show)


def user_input_plotmaker_looper(path='./1/'):
    user_input_plotmaker(path=path)
    print("Would you like to make another plot? [Y/n]")
    answer = input()
    if answer == 'Y' or answer == 'y':
        user_input_plotmaker_looper(path=path)
    elif answer == 'N' or answer == 'n':
        print("Terminated plotmaker. Goodbye.")
    else:
        print("Invalid input, terminated plotmaker. Goodbye.")


user_input_plotmaker_looper()
