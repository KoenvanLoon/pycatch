import pathlib

### BEGINNING ADDITIONS ###
filtering = False
### END ADDITIONS ###

# use for other runs
# numberOfTimeSteps=5200   # test run
numberOfTimeSteps = 104000  # Y2K long run

# option to fix both the regolith and the vegetation, not typically used
# in normal simulations
fixedStates = False

# option to call the methods that change the geomorphology
# this is typically on
changeGeomorphology = True

# number of realizations
nrOfSamples = 1

# Saving map files which can be used for spatial EWS
map_data = True
mean_timeseries_data = True
interval_map_snapshots = 100

# option to define some parameters as realizations from random
# functions, typically False
createRealizations = False

theDurationOfRainstorm = 2.0

# option to define the variables to report to disc
# these are defined in the modules, it is typically either 'full' or 'filtering'
# setOfVariablesToReport = 'filtering'
setOfVariablesToReport = 'full'

# state_variables_for_ews = 'full'
state_variables_for_ews = ['bioM', 'bioA', 'moiM', 'moiA']

# calculate upstream totals (with accuflux) in subsurfacewateronelayer module
# and interceptionuptomaxstore module
# may be needed for some reports
# and possibly for budget checks (if one needs these)
# normal use: set to False
calculateUpstreamTotals = False

#
# Reporting for the model components
#
if setOfVariablesToReport == 'full':
    interception_report_rasters = ["Vo", "Vi", "Vgf",
                                   "Vms"]  # reports of totals (Vot) only make sense if calculateUpstreamTotals is True
    infiltration_report_rasters = ["Ii", "Is", "Iks"]
    runoff_report_rasters = ["Rq", "Rqs"]
    subsurface_report_rasters = ["Gs",
                                 "Go"]  # reports of totals (Gxt, Got) only make sense if calculateUpstreamTotals is True
    surfacestore_report_rasters = ["Ss", "Sc"]
    rainfalleventsfromgammadistribution_report_rasters = ["Pf"]
    exchange_report_rasters = ["Xrc"]
    soilwashMMF_report_rasters = ["Wde", "Wdm", "Wfl"]
    regolith_report_rasters = ["Ast"]
    bedrockweathering_report_rasters = ["Cwe"]
    evapotranspirationsimple_report_rasters = ["Ep", "Ea"]
    biomassmodifiedmay_report_rasters = ["Xs"]
    baselevel_report_rasters = ["Ll"]
    creep_report_rasters = ["Ds"]
elif setOfVariablesToReport == 'filtering':
    interception_report_rasters = []  # reports of totals (Vot) only make sense if calculateUpstreamTotals is True
    infiltration_report_rasters = []
    runoff_report_rasters = []
    subsurface_report_rasters = []  # reports of totals (Gxt, Got) only make sense if calculateUpstreamTotals is True
    surfacestore_report_rasters = []
    rainfalleventsfromgammadistribution_report_rasters = []
    soilwashMMF_report_rasters = []
    exchange_report_rasters = ["Xrc", "Xra"]
    soilwashMMF_report_rasters = []
    regolith_report_rasters = []
    bedrockweathering_report_rasters = []
    evapotranspirationsimple_report_rasters = []
    biomassmodifiedmay_report_rasters = []
    baselevel_report_rasters = []
    creep_report_rasters = []
