"""
EWS - Early Warning Signals
EWS Configuration (cfg)

@authors: KoenvanLoon & TijmenJanssen
"""

import pathlib

# Pycatch configuration
"""
Settings for hour/week pycatch models.

Args:
-----

inputFolder : str, folder where inputs for the hourly model are stored. Inputs for the weekly model is 'weekly_inputs'
    by default.
    
cloneString : str, .map file location of 'clone.map' for hourly model.

locations : str, .map file location of 'locs.map' for hourly model.

number_of_timesteps_hourly/weekly : The number of steps in time for which the respective model (hourly/weekly) 
    calculates model outputs.

map_data : Selects whether spatial data is saved in map files. Usually True.

interval_map_snapshots : The number of timesteps between the optional saving of spatial data in map files.

mean_timeseries_data : Selects whether temporal data is saved in numpy.txt files. The saved value is the mean of the
    modelled spatial data. Usually True.

loc_timeseries_data : Selects whether temporal data is saved in numpy.txt files. The saved value is the value found at
    specified locations stored in ./inputs_weekly/mlocs.map . Usually True.

setOfVariablesToReport : Selects which set of variables are reported, either 'full' or 'filtering'. These are passed to 
    the class of a component where the variables that can be reported can be defined. Can be found at the end of this
    file.

timesteps_to_report_all_hourly/weekly : Defines timesteps at which all variables are reported.

timesteps_to_report_some_hourly/weekly : Defines timesteps at which some variables are reported.

doReportComponentsDynamicAsNumpy : Switch to report for locations in the hourly model as "small" numpy files, mainly 
    used for particle filtering.

nrOfSamples : The number of Monte Carlo (MC) samples or particles, realizations are written in folder(s) 1, 2, ...

readDistributionOfParametersFromDisk : When True, one can read a set of parameters for all Monte Carlo realizations from
    disk (e.g. representing probability distributions from a calibration). First time users should use False.

filtering : When True, a particle filtering run is done. Usually False for first time users.

createRealizations : Selects whether a single, given value is used for a number of parameters, or whether a realization
    for that parameter is drawn randomly. Usually False for first time users.

calculateUpstreamTotals : Selects whether upstream totals are calculated (accuflux) in the subsurfacewateronelayer and
    interceptionuptomaxstore modules. May be needed for some reports and possibly budget checks (if one needs these).
    For normal use, this is set to False.

fixedStates : Option to fix both the regolith and the vegetation states (week model). Usually False.

changeGeomorphology : Option to call on the methods that change the geomorphology (week model). Usually True.

swapCatchments : Switch to swap parameter values between two catchments, first time users will need to set this to 
    False.

rainstorm_probability : Chance of rainstorm per week, e.g. if set to 0.4, there is a 40% chance of rain per week.

rainstorm_duration : The time in hours for modelled rainfall events.

rainstorm_expected_intensity : Expected value of the intensity of a given rainstorm.

rainstorm_gamma_shape_param : Gamma shape parameter used for rainstorm intensity.

maxSurfaceStoreValue : int/flt, the maximum store of water at the surface level (in m).

saturatedConductivityMetrePerDayValue : int/flt, the saturated conductivity in meter per day.

limitingPointFractionValue : float, the lower limit of soil moisture fraction.

mergeWiltingPointFractionFSValue : float, the lower limit of soil moisture fraction at which plants start to wilt.

fieldCapacityFractionValue : float, the equilibrium of soil moisture fraction in the field.

initialSoilMoistureFractionCFG : float, the initial soil moisture fraction at the start of the model run.

soilPorosityFractionValue : float, the porosity of the soil.

rel_start_grazing : Ratio between 0 and 1 which sets the starting point over the number of timesteps

tot_increase_grazing : The total increase of grazing over the grazing period.

return_ini_grazing : Selects whether the grazing rates return to the initial value after the halfway point of the
    grazing period. Usually False. (Note that this halfway point occurs on (1 - rel_start) * total time / 2) .
    
startTimeYearValue , startTimeMonthValue , startTimeDayValue : int, starting time of model run specified by YYYY/MM/DD

withShading : bool, selects whether shading is used or not in the hourly model. This is not used for the
    evaporationsimple module.

-----
"""

# Hourly model input folder/file(s)
inputFolder = "inputs_from_weekly"
cloneString = str(pathlib.Path(inputFolder, "clone.map"))
locations = str(pathlib.Path(inputFolder, "clone.map"))

# Duration of model runs
number_of_timesteps_hourly = 8760  # ~1y in hours (24 hours per day)
# number_of_timesteps_hourly = 2000  # smaller test-run
number_of_timesteps_weekly = 104000  # ~2000y in weeks (52 weeks per year)
# number_of_timesteps_weekly = 5200  # smaller test-run of ~100y in weeks

# Saving spatial data
map_data = True
interval_map_snapshots = 100

# Saving temporal data
mean_timeseries_data = True
loc_timeseries_data = True

# Reporting of variables
# setOfVariablesToReport = 'full'
setOfVariablesToReport = 'filtering'
# setOfVariablesToReport = 'None'

# Timesteps of reporting variables
timesteps_to_report_all_hourly = list(range(1, number_of_timesteps_hourly + 1, 1))  # TODO - set step-size to map interval?
timesteps_to_report_all_weekly = list(range(0, number_of_timesteps_weekly + 1, interval_map_snapshots))

timesteps_to_report_some_hourly = list(range(100, number_of_timesteps_hourly + 1, 100))  # TODO - set step-size to map interval?
timesteps_to_report_some_weekly = list(range(0, number_of_timesteps_weekly + 1, interval_map_snapshots))

# Hourly model report as numpy
doReportComponentsDynamicAsNumpy = False

# Monte Carlo (MC) realizations
nrOfSamples = 1

# Read set of parameters for all MC realizations from disk
readDistributionOfParametersFromDisk = False

# Particle filtering
filtering = False

# Create realizations
createRealizations = False

# Calculate upstream totals
calculateUpstreamTotals = False

# Fixed regolith and vegetation states
fixedStates = False

# Change geomorphology
changeGeomorphology = True

# Swap catchments
swapCatchments = False

# Rainstorm parameters
# # scenario: original
rainstorm_probability = 0.4
rainstorm_duration = 2.0
rainstorm_expected_intensity = 0.002
rainstorm_gamma_shape_param = 100

# Surface store
maxSurfaceStoreValue = 0.0001

# Groundwater (saturated flow) and soil (moisture) parameters
saturatedConductivityMetrePerDayValue = 12.5
limitingPointFractionValue = 0.05
mergeWiltingPointFractionFSValue = 0.019
fieldCapacityFractionValue = 0.22

initialSoilMoistureFractionCFG = 0.22
soilPorosityFractionValue = 0.43

# Grazing
rel_start_grazing = 0
tot_increase_grazing = 0.00025
return_ini_grazing = False

# Real-time of first time step of model run(s)
# ! - Note: This is now UTC time almost certainly at least for shading.
startTimeYearValue = 2005
startTimeMonthValue = 7
startTimeDayValue = 1

# Shading
with_shading = True  # TODO - Check if this can be deleted
if with_shading is False:
    fractionReceivedValue = 1.0
    fractionReceivedFlatSurfaceValue = 1.0

# EWS configuration
"""
Settings for EWS calculations.

   ! - Note that some information, such as timesteps and intervals, used in the pycatch models is also used in the EWS 
   calculations. Because of this, if for example map_data is set to False, it is assumed that no spatial data is 
   present for spatial EWS calculations.

Args:
-----

state_variables_for_ews_hourly/weekly : State variables for which early warning signals/statistics are calculated.
    Either list of strings or 'full'.

stepsInShift : int, number of snapshots for which the hourly model is going to be run in a specified period.

stepsTotal : int, number of snapshots for which the hourly model is going to be run over the whole period as specified
    for the weekly model. Evenly distributed before and after the stepsInShift period if possible.

generate_dummy_datasets : Selects whether null model realizations are generated or not.

nr_generated_datasets : Selects how many null model realizations are generated.

method_1 : If True, the given number of null model realizations are created by shuffling data (similar mean and 
    variance).

method_2 : If True, the given number of null model realizations are created with the same Fourier spectrum and 
    amplitudes (similar autocorrelations, mean and variance).

method_3 : If True, the given number of null model realizations are created with an AR(1) model fitted to the data
    (similar autocorrelations, mean and variance).

detrended_method : Either 'None' or 'Gaussian', selects whether none detrending occurs or Gaussian filtering detrending 
    using scipy.ndimage.gaussian_filter1d().
    
detrended_sigma : If detrended_method is 'Gaussian', selects the sigma used in scipy.ndimage.gaussian_filter1d().

save_detrended_data : Selects whether detrended temporal data used in the generation of null models is saved. Only
    relevent when detrending in EWS_weekly.py is not set to 'None'.

cutoff : Selects whether a cutoff point is used at a defined point for null model realizations, calculations and plots.
    Usually True, speeds up calculation time and makes null model realizations more accurate (see null models).

cutoff_point : Time at which states shift. Retrieved from plotting the biomass timeseries.

-----
"""

# State variables
state_variables_for_ews_hourly = ['Gs']
# state_variables_for_ews_hourly = 'full'  # - TODO check the 'full' list in EWS_StateVariables.py
# state_variables_for_ews_weekly = ['bioA', 'bioM']
state_variables_for_ews_weekly = 'full'  # - TODO check the 'full' list in EWS_StateVariables.py

# Loop over snapshots in hourly model
stepsInShift = 20
stepsTotal = 30

# Generate null models
generate_dummy_datasets = True
nr_generated_datasets = 1

# Methods for generated null models
method_1 = True
method_2 = True
method_3 = True

# Temporal data detrending
detrended_method = 'Gaussian'
detrended_sigma = 5000
save_detrended_data = True

# Cutoff transition
cutoff = False
cutoff_point = 96000  # TODO - Implement a way to cutoff data before/at (?) CT --> if tuple, elif int/scalar.

# Reporting for the model components (both hourly and weekly)
"""
Reporting model components of the pycatch model (see time_steps_to_report_all/some_hourly/weekly and 
timeStepsToReportRqs)

"""

if setOfVariablesToReport == 'full':
    interception_report_rasters = ["Vo", "Vi", "Vgf", "Vms"]
    #   reports of totals (Vot) only make sense if calculateUpstreamTotals is True
    infiltration_report_rasters_weekly = ["Ii", "Is", "Iks"]
    infiltration_report_rasters = ["Ii", "Ij", "Is", "Iks"]  # TODO - might want to rename this to ""_hourly, as above
    runoff_report_rasters = ["Rq", "Rqs"]
    subsurface_report_rasters = ["Gs", "Go"]
    #   reports of totals (Gxt, Got) only make sense if calculateUpstreamTotals is True
    shading_report_rasters = ["Mfs", "Msc", "Msh"]
    surfacestore_report_rasters = ["Ss", "Sc"]
    rainfalleventsfromgammadistribution_report_rasters = ["Pf"]
    exchange_report_rasters = ["Xrc"]
    soilwashMMF_report_rasters = ["Wde", "Wdm", "Wfl"]
    regolith_report_rasters = ["Ast"]
    bedrockweathering_report_rasters = ["Cwe"]
    evapotrans_report_rasters = ["Ep", "Epc"]
    evapotranspirationsimple_report_rasters = ["Ep", "Ea"]
    biomassmodifiedmay_report_rasters = ["Xs"]
    baselevel_report_rasters = ["Ll"]
    creep_report_rasters = ["Ds"]
    randomparameters_report_rasters = ["RPic", "RPks", "RPrt", "RPsc", "RPmm"]
elif setOfVariablesToReport == 'filtering':
    interception_report_rasters = []
    #   reports of totals (Vot) only make sense if calculateUpstreamTotals is True
    infiltration_report_rasters_weekly = ["Iks"]
    infiltration_report_rasters = []
    runoff_report_rasters = []
    subsurface_report_rasters = []
    #   reports of totals (Gxt, Got) only make sense if calculateUpstreamTotals is True
    shading_report_rasters = []
    surfacestore_report_rasters = []
    rainfalleventsfromgammadistribution_report_rasters = []
    exchange_report_rasters = []
    soilwashMMF_report_rasters = []
    regolith_report_rasters = []
    bedrockweathering_report_rasters = []
    evapotrans_report_rasters = []
    evapotranspirationsimple_report_rasters = []
    biomassmodifiedmay_report_rasters = []
    baselevel_report_rasters = []
    creep_report_rasters = []
    randomparameters_report_rasters = []
elif setOfVariablesToReport == 'None':
    interception_report_rasters = []
    #   reports of totals (Vot) only make sense if calculateUpstreamTotals is True
    infiltration_report_rasters_weekly = []
    infiltration_report_rasters = []
    runoff_report_rasters = []
    subsurface_report_rasters = []
    #   reports of totals (Gxt, Got) only make sense if calculateUpstreamTotals is True
    shading_report_rasters = []
    surfacestore_report_rasters = []
    rainfalleventsfromgammadistribution_report_rasters = []
    exchange_report_rasters = []
    soilwashMMF_report_rasters = []
    regolith_report_rasters = []
    bedrockweathering_report_rasters = []
    evapotrans_report_rasters = []
    evapotranspirationsimple_report_rasters = []
    biomassmodifiedmay_report_rasters = []
    baselevel_report_rasters = []
    creep_report_rasters = []
    randomparameters_report_rasters = []
