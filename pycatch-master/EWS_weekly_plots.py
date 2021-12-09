from pcraster import *
import numpy as np
import os
import time
from scipy import ndimage

import configuration_weekly as cfg
import NULL_models_timeseries as temp_NULL
import NULL_models_spatial as spat_NULL
import EWSPy as ews
import EWS_StateVariables as ews_sv

import matplotlib.pyplot as plt

## State variables for EWS ##
variables = ews_sv.variables # State variables present in EWS_StateVariables can be added through configuration_weekly

## Statistical EWS ##
ews_temporal_signals = ['t.mn', 't.std', 't.var', 't.cv', 't.skw', 't.krt', 't.dfa', 't.acr', 't.rr', 't.coh']
ews_spatial_signals = ['s.mn', 's.std', 's.var', 's.skw', 's.krt', 's.mI']

## Functions ## - TODO: Maybe plot for 2 variables? (e.g. bioA timeseries & bioM var)
def plot(variable, variable_signal1, variable_signal2='None', path='./1/', save=False, show=False):
    if variable.spatial:
        x_axis = np.arange(cfg.interval_map_snapshots, cfg.numberOfTimeSteps + cfg.interval_map_snapshots,
                           cfg.interval_map_snapshots)
    if variable.temporal:
        x_axis = np.arange(0, cfg.numberOfTimeSteps, variable.window_size - variable.window_overlap)

    fpath = os.path.join(path + variable.name + '.' + variable_signal1)
    variable_signal1_array = np.loadtxt(fpath + '.numpy.txt')
    plt.plot(x_axis, variable_signal1_array, label=f'{variable.name} {variable_signal1}')

    if variable_signal2 == 'timeseries':
        fname = ews.file_name_str(variable.name, cfg.numberOfTimeSteps)
        fpath = os.path.join(path + fname)
        variable_signal2 = np.loadtxt(fpath + '.numpy.txt')
        timeseries_x_axis = np.arange(0, cfg.numberOfTimeSteps, 1)
        plt.plot(timeseries_x_axis, variable_signal2, label=f'Continues measurement of {variable.name}')
    elif variable_signal2 != 'None':
        fpath = os.path.join(path + variable.name + '.' + variable_signal2)
        variable_signal2_array = np.loadtxt(fpath + '.numpy.txt')
        plt.plot(x_axis, variable_signal2_array, label=f'{variable.name} {variable_signal2}')

    plt.xlabel('x - axis')
    plt.ylabel('y - axis')
    plt.title('Title')
    plt.legend()

    if save:
        if variable_signal2 != 'None':
            plt.savefig(path + f"{variable.name}_{variable_signal1}_{variable_signal2}.pdf", format="pdf")
        else:
            plt.savefig(path + f"{variable.name}_{variable_signal1}.pdf", format="pdf")
    if show:
        plt.show()

for variable in variables:
    plot(variable, 't.var', save=False, show=True)
