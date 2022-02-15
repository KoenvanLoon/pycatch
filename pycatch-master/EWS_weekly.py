import EWSPy as ews
from pcraster import *
import numpy as np
import os
import time
from scipy import ndimage

import EWS_configuration as cfg
import NULL_models_timeseries_weekly as temp_NULL
import NULL_models_spatial_weekly as spat_NULL
import EWS_StateVariables as ews_sv

### User input ###

## State variables for EWS ##
variables = ews_sv.variables_weekly  # State variables present in EWS_StateVariables can be added through configuration

### End user input ###

### Information from configuration weekly ###
## Realizations/MC runs ##
realizations = cfg.nrOfSamples
# realizations = 1 # for test cases

## Timesteps, intervals ##
spatial_ews_present = cfg.map_data

if not cfg.cutoff:
    spatial_ews_interval = np.arange(cfg.interval_map_snapshots, cfg.number_of_timesteps_weekly +
                                     cfg.interval_map_snapshots, cfg.interval_map_snapshots)
elif cfg.cutoff:
    spatial_ews_interval = np.arange(cfg.interval_map_snapshots, cfg.cutoff_point + cfg.interval_map_snapshots,
                                     cfg.interval_map_snapshots)

temporal_ews_present = cfg.mean_timeseries_data
temporal_ews_interval = cfg.number_of_timesteps_weekly  # the interval is defined as t=0 : the last timestep.


### Functions ###

def time_series2time_windows(timeseries, window_size=100, window_overlap=0):
    actual_window_overlap = window_size - window_overlap
    sh = (timeseries.size - window_size + 1, window_size)
    st = timeseries.strides * 2
    if window_overlap != 0:
        view = np.lib.stride_tricks.as_strided(timeseries, strides=st, shape=sh)[0::actual_window_overlap]
    elif window_overlap == 0:
        view = np.lib.stride_tricks.as_strided(timeseries, strides=st, shape=sh)[0::window_size]
    return view

def generate_datasets(variable, path='./1/', nr_realizations=1, detrending_temp='None', sigma=50,
                      method1=False, method2=False, method3=False):
    if variable.temporal:
        state_variable_timeseries, files_present = temporal_data_file_loading(variable, path=path)

        if files_present:
            if state_variable_timeseries.ndim == 1:
                # Detrending: 'None', 'Gaussian'
                if detrending_temp == 'None':
                    temp_NULL.detrend_(state_variable_timeseries, realizations=nr_realizations, path=path,
                                       file_name=variable.name)
                elif detrending_temp == 'Gaussian':  # TODO - Multiple sigmas?
                    gaussian_filter = ndimage.gaussian_filter1d(state_variable_timeseries, sigma)
                    state_variable_timeseries = state_variable_timeseries - gaussian_filter
                    if cfg.save_detrended_data:
                        # Only 1 realization, as the detrending method parameters do not change.
                        temp_NULL.detrend_(state_variable_timeseries, gaussian_filter, realizations=nr_realizations,
                                           path=path, file_name=variable.name)

                # Generate dummy datasets
                if method1:
                    temp_NULL.method1_(state_variable_timeseries, realizations=nr_realizations, path=path,
                                       file_name=variable.name)
                if method2:
                    temp_NULL.method2_(state_variable_timeseries, realizations=nr_realizations, path=path,
                                       file_name=variable.name)
                if method3:
                    temp_NULL.method3_(state_variable_timeseries, realizations=nr_realizations, path=path,
                                       file_name=variable.name)
            else:
                print(f"Multiple dimensions are currently not supported for generated datasets, so no datasets are being "
                      f"generated for {variable.name}.")

    if variable.spatial:
        state_variable_snapshots, files_present = spatial_data_file_loading(variable, path=path)
        if files_present:
            state_variable_snapshots = np.asarray(state_variable_snapshots)

            # Generate dummy datasets
            if method1:
                spat_NULL.method1_(state_variable_snapshots, realizations=nr_realizations, path=path,
                                   file_name=variable.name)
            if method2:
                spat_NULL.method2_(state_variable_snapshots, realizations=nr_realizations, path=path,
                                   file_name=variable.name)
            if method3:
                spat_NULL.method3_(state_variable_snapshots, realizations=nr_realizations, path=path,
                                   file_name=variable.name)


def ews_calculations_generated_datasets_init(variable, path='./1/', nr_realizations=1, timer_on=False, method1=False,
                                             method2=False, method3=False):
    generated_number_length = ews.generated_number_length(nr_realizations)

    if cfg.save_detrended_data and variable.temporal:
        ews_calculations_generated_datasets_main(variable, 'dtr', gen_nr_len=generated_number_length, path=path,
                                                 nr_realizations=1, timer_on=timer_on)
    if method1:
        ews_calculations_generated_datasets_main(variable, 'm1g', gen_nr_len=generated_number_length, path=path,
                                                 nr_realizations=nr_realizations, timer_on=timer_on)
    if method2:
        ews_calculations_generated_datasets_main(variable, 'm2g', gen_nr_len=generated_number_length, path=path,
                                                 nr_realizations=nr_realizations, timer_on=timer_on)
    if method3:
        ews_calculations_generated_datasets_main(variable, 'm3g', gen_nr_len=generated_number_length, path=path,
                                                 nr_realizations=nr_realizations, timer_on=timer_on)


def ews_calculations_generated_datasets_main(variable, method, gen_nr_len=4, path='./1/', nr_realizations=1, timer_on=False):
    for realization in range(nr_realizations):
        generated_number_string = method + str(realization).zfill(gen_nr_len) + '/'
        dir_name = os.path.join(path + generated_number_string)
        ews_calculations_init(variable, path=dir_name, timer_on=timer_on)


def temporal_data_file_loading(variable, path='./1/'):
    state_variable_timeseries = []
    EWS_calculations = True
    if variable.datatype == 'numpy':
        file_name = ews.file_name_str(variable.name, temporal_ews_interval)
        if os.path.exists(path + file_name + ".numpy.txt"):
            state_variable_timeseries = np.loadtxt(path + file_name + ".numpy.txt")
        else:
            print(f"{file_name}.numpy.txt not found in {path}")
            EWS_calculations = False
    else:
        print(f"Datatype for {variable.name} currently not supported.")
        EWS_calculations = False

    return state_variable_timeseries, EWS_calculations


def window_stacker(variable, state_variable_timeseries):
    nr_dim = state_variable_timeseries.ndim
    if nr_dim == 1:
        if cfg.cutoff:
            state_variable_timeseries = state_variable_timeseries[:cfg.cutoff_point]
        stack_of_windows = time_series2time_windows(state_variable_timeseries, variable.window_size,
                                                    variable.window_overlap)
        stack_x, stack_y = np.asarray(stack_of_windows).shape
    else:
        stack_of_windows = [0.0] * np.asarray(state_variable_timeseries).shape[1]
        for k, timeseries in enumerate(state_variable_timeseries.T):
            if cfg.cutoff:
                stack_of_windows[k] = time_series2time_windows(timeseries[:cfg.cutoff_point],
                                                               variable.window_size, variable.window_overlap)
            elif not cfg.cutoff:
                stack_of_windows[k] = time_series2time_windows(timeseries, variable.window_size,
                                                               variable.window_overlap)
        stack_x, stack_y, stack_z = np.asarray(stack_of_windows).shape
        stack_of_windows = np.asarray(stack_of_windows).reshape(-1, stack_z)
    return stack_of_windows, nr_dim, stack_x, stack_y


def spatial_data_file_loading(variable, path='./1/'):
    state_variable_snapshots = [0.0] * len(spatial_ews_interval)
    EWS_calculations = True
    if variable.datatype == 'numpy':
        for k, timestep in enumerate(spatial_ews_interval):
            file_name = ews.file_name_str(variable.name, timestep)
            if os.path.exists(path + file_name + ".numpy.txt"):
                state_variable_snapshots[k] = np.loadtxt(path + file_name + 'numpy.txt')
            else:
                print(f"{file_name}.numpy.txt not found in {path}.")
                EWS_calculations = False

    if variable.datatype == 'map':
        for k, timestep in enumerate(spatial_ews_interval):
            file_name = ews.file_name_str(variable.name, timestep)
            if os.path.exists(path + file_name):
                state_variable_snapshots[k] = pcr2numpy(readmap(path + file_name), np.NaN)
            else:
                print(f"{file_name} not found in {path}.")
                EWS_calculations = False
    else:
        print(f"Datatype for {variable.name} currently not supported.")
        EWS_calculations = False

    return state_variable_snapshots, EWS_calculations


def ews_calc_and_save(variable, data, method_name, method_function, path='./1/', nr_dim=1, stack_x=1, stack_y=1):
    fpath = os.path.join(path + variable.name + method_name)
    signal = method_function(data)
    if nr_dim > 1:
        signal = signal.reshape(stack_x, stack_y)
    np.savetxt(fpath + '.numpy.txt', signal)


def temporal_ews_calculations(variable, state_variable_timeseries, path='./1/'):
    stack_of_windows, nr_dim, stack_x, stack_y = window_stacker(variable, state_variable_timeseries)

    ews_calc_and_save(variable, stack_of_windows, '.t.mn', ews.temporal_mean, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.std', ews.temporal_std, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.var', ews.temporal_var, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.cv', ews.temporal_cv, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.skw', ews.temporal_skw, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.krt', ews.temporal_krt, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.acr', ews.temporal_autocorrelation, path=path, nr_dim=nr_dim,
                      stack_x=stack_x, stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.AR1', ews.temporal_AR1, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.rr', ews.temporal_returnrate, path=path, nr_dim=nr_dim,
                      stack_x=stack_x, stack_y=stack_y)

    # Temporal dfa TODO - returns 3-4 values, save only 1?
    fpath = os.path.join(path + variable.name + '.t.dfa')
    _, _, _, temporal_statistic = ews.temporal_dfa(stack_of_windows, window_size=variable.window_size)
    # scales, fluct, coeff, propagator (== temporal statistic)
    if nr_dim > 1:
        temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
    np.savetxt(fpath + '.numpy.txt', temporal_statistic)

    # Temporal cond. het. TODO - returns 2 values, save only 1?
    fpath = os.path.join(path + variable.name + '.t.coh')
    save_p = True
    if save_p and nr_dim == 1:
        temporal_statistic = [[0.0], [0.0]]
        statistic, p_val = ews.temporal_cond_het(stack_of_windows)
        temporal_statistic[0] = statistic
        temporal_statistic[1] = p_val
    else:
        temporal_statistic, _ = ews.temporal_cond_het(stack_of_windows)  # _ is the p-value of the test, not saved
        if nr_dim > 1:
            temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
    np.savetxt(fpath + '.numpy.txt', temporal_statistic)


def spatial_ews_calculations(variable, state_variable_maps, path='./1/'):
    ews_calc_and_save(variable, state_variable_maps, '.s.mn', ews.spatial_mean, path=path)
    ews_calc_and_save(variable, state_variable_maps, '.s.std', ews.spatial_std, path=path)
    ews_calc_and_save(variable, state_variable_maps, '.s.var', ews.spatial_var, path=path)
    ews_calc_and_save(variable, state_variable_maps, '.s.skw', ews.spatial_skw, path=path)
    ews_calc_and_save(variable, state_variable_maps, '.s.krt', ews.spatial_krt, path=path)
    ews_calc_and_save(variable, state_variable_maps, '.s.mI', ews.spatial_corr, path=path)
    # ews_calc_and_save(variable, state_variable_maps, '.s.dft', ews.spatial_DFT, path=path)


def ews_calculations_init(variable, path='./1/', timer_on=False):
    if variable.temporal:
        if temporal_ews_present:
            ews_calculations_main(variable, temporal_data_file_loading, temporal_ews_calculations, path=path,
                                  timer_on=timer_on)
        elif not temporal_ews_present:
            print(f"Mean timeseries data == False in configuration_weekly.py, could not calculate EWS for "
                  f"{variable.name}.")
    elif variable.spatial:
        if spatial_ews_present:
            ews_calculations_main(variable, spatial_data_file_loading, spatial_ews_calculations, path=path,
                                  timer_on=timer_on)
        elif not spatial_ews_present:
            print(f"Map data == False in configuration, could not calculate EWS for {variable.name}.")


def ews_calculations_main(variable, loading_function, calculation_function, path='./1/', timer_on=False):
    state_variable, files_present = loading_function(variable, path=path)

    if files_present:
        print(f"Started EWS calculations for {variable.name} in {path}")
        if timer_on:
            start_time = time.time()

        calculation_function(variable, state_variable, path=path)

        if timer_on:
            end_time = time.time()
            print(f"Elapsed time for EWS calculations for {variable.name} in {path} equals:", end_time - start_time,
                  '\n')


def EWS_weekly_calculations():
    start_time = time.time()
    for realization in range(1, realizations + 1):
        for variable in variables:
            ews_calculations_init(variable, path=f'./{realization}/', timer_on=True)
            if cfg.generate_dummy_datasets:
                generate_datasets(variable, path=f'./{realization}/', nr_realizations=cfg.nr_generated_datasets,
                                  detrending_temp='None', sigma=100, method1=cfg.method_1, method2=cfg.method_2,
                                  method3=cfg.method_3)
                ews_calculations_generated_datasets_init(variable, path=f'./{realization}/',
                                                         nr_realizations=cfg.nr_generated_datasets,
                                                         timer_on=True, method1=cfg.method_1, method2=cfg.method_2,
                                                         method3=cfg.method_3)
    end_time = time.time() - start_time
    print(f"Total elapsed time equals: {end_time} seconds")


EWS_weekly_calculations()


# def time_series2time_windows(time_series, window_size=100, window_overlap=0):
#     return np.array([time_series[i:i + window_size] for i in range(0, len(time_series), window_size - window_overlap)])


# def ews_calculations(variable, path='./1/', timer_on=False):
#     ## Temporal EWS calculations ##
#     EWS_calculations = True  # If files are not found, set to False; no further calculations done on this state variable
#
#     if variable.temporal:
#         if temporal_ews_present:
#             print(f"Started temporal EWS calculations for {variable.name}")
#
#             ## Start timer if set to True ##
#             if timer_on:
#                 start_time = time.time()
#
#             ## Timeseries file loading ##
#             state_variable_timeseries = []
#             if variable.datatype == 'numpy':
#                 file_name = ews.file_name_str(variable.name, temporal_ews_interval)
#                 if os.path.exists(path + file_name + ".numpy.txt"):
#                     state_variable_timeseries = np.loadtxt(path + file_name + ".numpy.txt")
#                 else:
#                     print(f"{file_name}.numpy.txt not found in {path}")
#                     EWS_calculations = False
#             else:
#                 print(f"Datatype for {variable.name} currently not supported.")
#
#             if EWS_calculations:
#                 ## Splitting timeseries into (overlapping) windows ##
#                 # if cutoff:
#                 #     state_variable_timeseries = state_variable_timeseries[:cutoff_point]
#
#                 if state_variable_timeseries.ndim == 1:
#                     if cfg.cutoff:
#                         state_variable_timeseries = state_variable_timeseries[:cfg.cutoff_point]
#                     stack_of_windows = time_series2time_windows(state_variable_timeseries, variable.window_size,
#                                                                 variable.window_overlap)
#                 else:
#                     stack_of_windows = [0.0] * np.asarray(state_variable_timeseries).shape[1]
#                     for k, timeseries in enumerate(state_variable_timeseries.T):
#                         if cfg.cutoff:
#                             stack_of_windows[k] = time_series2time_windows(timeseries[:cfg.cutoff_point],
#                                                                            variable.window_size, variable.window_overlap)
#                         else:
#                             stack_of_windows[k] = time_series2time_windows(timeseries, variable.window_size,
#                                                                            variable.window_overlap)
#                     stack_x, stack_y, stack_z = np.asarray(stack_of_windows).shape
#                     stack_of_windows = np.asarray(stack_of_windows).reshape(-1, stack_z)
#
#                 ## EWS calculations ###
#                 # Temporal mean #
#                 fpath = os.path.join(path + variable.name + '.t.mn')
#                 temporal_statistic = ews.temporal_mean(stack_of_windows)
#                 if state_variable_timeseries.ndim > 1:
#                     temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#                 # Temporal std #
#                 fpath = os.path.join(path + variable.name + '.t.std')
#                 temporal_statistic = ews.temporal_std(stack_of_windows)
#                 if state_variable_timeseries.ndim > 1:
#                     temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#                 # Temporal var #
#                 fpath = os.path.join(path + variable.name + '.t.var')
#                 temporal_statistic = ews.temporal_var(stack_of_windows)
#                 if state_variable_timeseries.ndim > 1:
#                     temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#                 # Temporal cv #
#                 fpath = os.path.join(path + variable.name + '.t.cv')
#                 temporal_statistic = ews.temporal_cv(stack_of_windows)
#                 if state_variable_timeseries.ndim > 1:
#                     temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#                 # Temporal skw #
#                 fpath = os.path.join(path + variable.name + '.t.skw')
#                 temporal_statistic = ews.temporal_skw(stack_of_windows)
#                 if state_variable_timeseries.ndim > 1:
#                     temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#                 # Temporal krt #
#                 fpath = os.path.join(path + variable.name + '.t.krt')
#                 temporal_statistic = ews.temporal_krt(stack_of_windows)
#                 if state_variable_timeseries.ndim > 1:
#                     temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#                 # Temporal dfa #
#                 fpath = os.path.join(path + variable.name + '.t.dfa')
#                 _, _, _, temporal_statistic = ews.temporal_dfa(stack_of_windows, window_size=variable.window_size)  # scales, fluct, coeff, propagator
#                 if state_variable_timeseries.ndim > 1:
#                     temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#                 # Temporal autocorr. #
#                 fpath = os.path.join(path + variable.name + '.t.acr')
#                 temporal_statistic = ews.temporal_autocorrelation(stack_of_windows)
#                 if state_variable_timeseries.ndim > 1:
#                     temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#                 # Temporal AR1 #
#                 fpath = os.path.join(path + variable.name + '.t.AR1')
#                 temporal_statistic = ews.temporal_AR1(stack_of_windows)
#                 if state_variable_timeseries.ndim > 1:
#                     temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#                 # Temporal return rate #
#                 fpath = os.path.join(path + variable.name + '.t.rr')
#                 temporal_statistic = ews.temporal_returnrate(stack_of_windows)
#                 if state_variable_timeseries.ndim > 1:
#                     temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#                 # Temporal cond. het. #
#                 fpath = os.path.join(path + variable.name + '.t.coh')
#                 save_p = True
#                 if save_p and state_variable_timeseries.ndim == 1:
#                     temporal_statistic = [[0.0], [0.0]]
#                     statistic, p_val = ews.temporal_cond_het(stack_of_windows)
#                     temporal_statistic[0] = statistic
#                     temporal_statistic[1] = p_val
#                 else:
#                     temporal_statistic, _ = ews.temporal_cond_het(stack_of_windows)  # _ is the p-value of the test, not saved
#                     if state_variable_timeseries.ndim > 1:
#                         temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
#                 np.savetxt(fpath + '.numpy.txt', temporal_statistic)
#
#             ## End timer if set to True##
#             if timer_on:
#                 end_time = time.time()
#                 print("Elapsed time for temporal EWS calculations equals:", end_time - start_time, '\n')
#
#         elif temporal_ews_present == False:
#             print(
#                 f"Mean timeseries data == False in configuration_weekly.py, could not calculate EWS for {variable.name}.")
#
#     ## Spatial EWS calculations ##
#     EWS_calculations = True  # If files are not found, set to False; no further calculations done on this state variable
#
#     if variable.spatial:
#         if spatial_ews_present:
#             print(f"Started spatial EWS calculations for {variable.name}")
#
#             ## Start timer if set to True ##
#             if timer_on:
#                 start_time = time.time()
#
#             ## Spatial maps file loading ##
#             state_variable_snapshots = [0.0] * len(spatial_ews_interval)
#             if variable.datatype == 'numpy':
#                 for k, timestep in enumerate(spatial_ews_interval):
#                     file_name = ews.file_name_str(variable.name, timestep)
#                     if os.path.exists(path + file_name + ".numpy.txt"):
#                         state_variable_snapshots[k] = np.loadtxt(path + file_name + 'numpy.txt')
#                     else:
#                         print(f"{file_name}.numpy.txt not found in {path}.")
#                         EWS_calculations = False
#
#             if variable.datatype == 'map':
#                 for k, timestep in enumerate(spatial_ews_interval):
#                     file_name = ews.file_name_str(variable.name, timestep)
#                     if os.path.exists(path + file_name):
#                         state_variable_snapshots[k] = pcr2numpy(readmap(path + file_name), np.NaN)
#                     else:
#                         print(f"{file_name} not found in {path}.")
#                         EWS_calculations = False
#             else:
#                 print(f"Datatype for {variable.name} currently not supported.")
#
#             ## EWS calculations ##
#             if EWS_calculations:
#                 state_variable_snapshots = np.asarray(state_variable_snapshots)
#
#                 # Spatial mean #
#                 fpath = os.path.join(path + variable.name + '.s.mn')
#                 np.savetxt(fpath + '.numpy.txt', ews.spatial_mean(state_variable_snapshots))
#
#                 # Spatial std #
#                 fpath = os.path.join(path + variable.name + '.s.std')
#                 np.savetxt(fpath + '.numpy.txt', ews.spatial_std(state_variable_snapshots))
#
#                 # Spatial var #
#                 fpath = os.path.join(path + variable.name + '.s.var')
#                 np.savetxt(fpath + '.numpy.txt', ews.spatial_var(state_variable_snapshots))
#
#                 # Spatial skw #
#                 fpath = os.path.join(path + variable.name + '.s.skw')
#                 np.savetxt(fpath + '.numpy.txt', ews.spatial_skw(state_variable_snapshots))
#
#                 # Spatial krt #
#                 fpath = os.path.join(path + variable.name + '.s.krt')
#                 np.savetxt(fpath + '.numpy.txt', ews.spatial_krt(state_variable_snapshots))
#
#                 # Spatial correlation (Moran's I) #
#                 fpath = os.path.join(path + variable.name + '.s.mI')
#                 np.savetxt(fpath + '.numpy.txt', ews.spatial_corr(state_variable_snapshots))
#
#                 # # spatial DFT #
#                 # fpath = os.path.join(path + variable.name + '.s.dft')
#                 # np.savetxt(fpath + '.numpy.txt', ews.spatial_DFT(state_variable_snapshots))
#
#             ## End timer if set to True##
#             if timer_on:
#                 end_time = time.time()
#                 print("Elapsed time for temporal EWS calculations equals:", end_time - start_time, '\n')
#
#         elif spatial_ews_present == False:
#             print(f"Map data == False in configuration, could not calculate EWS for {variable.name}.")


# def spatial_ews_calculations(variable, state_variable_maps, path='./1/'):
#     state_variable_snapshots = np.asarray(state_variable_maps)
#
#     # Spatial mean
#     fpath = os.path.join(path + variable.name + '.s.mn')
#     np.savetxt(fpath + '.numpy.txt', ews.spatial_mean(state_variable_snapshots))
#
#     # Spatial std
#     fpath = os.path.join(path + variable.name + '.s.std')
#     np.savetxt(fpath + '.numpy.txt', ews.spatial_std(state_variable_snapshots))
#
#     # Spatial var
#     fpath = os.path.join(path + variable.name + '.s.var')
#     np.savetxt(fpath + '.numpy.txt', ews.spatial_var(state_variable_snapshots))
#
#     # Spatial skw
#     fpath = os.path.join(path + variable.name + '.s.skw')
#     np.savetxt(fpath + '.numpy.txt', ews.spatial_skw(state_variable_snapshots))
#
#     # Spatial krt
#     fpath = os.path.join(path + variable.name + '.s.krt')
#     np.savetxt(fpath + '.numpy.txt', ews.spatial_krt(state_variable_snapshots))
#
#     # Spatial correlation (Moran's I)
#     fpath = os.path.join(path + variable.name + '.s.mI')
#     np.savetxt(fpath + '.numpy.txt', ews.spatial_corr(state_variable_snapshots))
#
#     # # spatial DFT
#     # fpath = os.path.join(path + variable.name + '.s.dft')
#     # np.savetxt(fpath + '.numpy.txt', ews.spatial_DFT(state_variable_snapshots))


# def ews_calculations_generated_datasets_init(variable, path='./1/', nr_realizations=1, timer_on=False, method1=False,
#                                              method2=False, method3=False):
#     generated_number_length = 4
#     if len(str(realizations)) > 4:
#         generated_number_length = len(str(realizations))
#
#     if cfg.save_detrended_data and variable.temporal:
#         generated_number_string = 'dtr' + str(0).zfill(generated_number_length) + '/'
#         dir_name = os.path.join(path + generated_number_string)
#         ews_calculations(variable, path=dir_name, timer_on=timer_on)
#
#     if method1:
#         for realization in range(nr_realizations):
#             generated_number_string = 'm1g' + str(realization).zfill(generated_number_length) + '/'
#             dir_name = os.path.join(path + generated_number_string)
#             ews_calculations(variable, path=dir_name, timer_on=timer_on)
#     if method2:
#         for realization in range(nr_realizations):
#             generated_number_string = 'm2g' + str(realization).zfill(generated_number_length) + '/'
#             dir_name = os.path.join(path + generated_number_string)
#             ews_calculations(variable, path=dir_name, timer_on=timer_on)
#     if method3:
#         for realization in range(nr_realizations):
#             generated_number_string = 'm3g' + str(realization).zfill(generated_number_length) + '/'
#             dir_name = os.path.join(path + generated_number_string)
#             ews_calculations(variable, path=dir_name, timer_on=timer_on)


# def ews_calculations(variable, path='./1/', timer_on=False):
#     # Temporal EWS calculations
#
#     if variable.temporal:
#         if temporal_ews_present:
#             print(f"Started temporal EWS calculations for {variable.name}")
#
#             # Timer start
#             if timer_on:
#                 start_time = time.time()
#
#             # Timeseries file loading
#             state_variable_timeseries, files_present = temporal_data_file_loading(variable, path=path)
#
#             if files_present:
#                 # Temporal EWS calculations
#                 temporal_ews_calculations(variable, state_variable_timeseries, path=path)
#
#             # Timer end
#             if timer_on:
#                 end_time = time.time()
#                 print("Elapsed time for temporal EWS calculations equals:", end_time - start_time, '\n')
#
#         elif not temporal_ews_present:
#             print(f"Mean timeseries data == False in configuration_weekly.py, could not calculate EWS for {variable.name}.")
#
#     elif variable.spatial:
#         if spatial_ews_present:
#             print(f"Started spatial EWS calculations for {variable.name}")
#
#             # Timer start
#             if timer_on:
#                 start_time = time.time()
#
#             # Spatial data file loading
#             state_variable_maps, files_present = spatial_data_file_loading(variable, path=path)
#
#             if files_present:
#                 # Spatial EWS calculations
#                 spatial_ews_calculations(variable, state_variable_maps, path=path)
#
#             # Timer end
#             if timer_on:
#                 end_time = time.time()
#                 print("Elapsed time for temporal EWS calculations equals:", end_time - start_time, '\n')
#
#         elif not spatial_ews_present:
#             print(f"Map data == False in configuration, could not calculate EWS for {variable.name}.")
