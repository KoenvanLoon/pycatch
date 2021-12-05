# R package loader ### ADDED - KL
import rpy2
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)   # select the first mirror in the list for R packages
packageNamesR = ('gstat', 'automap', 'e1071', 'tseries')

# check rpy2 version and R packages for variogram calculations
if rpy2.__version__ != '2.9.4':
    print("Tested for rpy2 version 2.9.4, current version is", rpy2.__version__)
    print("Please make sure you use the correct version.")
names_to_install = [x for x in packageNamesR if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    print(f"Installing the following R packages: {names_to_install}")
    utils.install_packages(StrVector(names_to_install))
### END OF ADDITIONS ###

# BIOMASS
variable = self.d_biomassModifiedMay.biomass
variableSampled = ifthen(self.someLocs, variable)

self.historyOfBiomass = generalfunctions_test01.keepHistoryOfMaps(self.historyOfBiomass,
                                                                  variableSampled,
                                                                  self.durationHistory)
stackOfMapsAsListVariable = list(self.historyOfBiomass)

if save_maps:
    generalfunctions_test01.report_as_map(variable, 'bioM', self.currentSampleNumber(), self.currentTimeStep())

    # if cfg.variances:
    #     # Test case
    #     # bins, gamma = generalfunctions_test01.variogramValuesKoen(stackOfMapsAsListVariable, boundVector)
    #     # numpy.savetxt(generateNameST('biTS', self.currentSampleNumber(), self.currentTimeStep()),
    #     #               numpy.array(gamma))
    #
    #     # temporal
    #     # dist, gamma = generalfunctions_test01.experimentalVariogramValues(stackOfMapsAsListVariable,
    #                                                                       boundVector, 0, 1,
    #                                                                       'test', 2.0)
    #     #dist, gamma = generalfunctions_test01.experimentalVariogramValuesInTime(stackOfMapsAsListVariable,
    #     #                                                                        list(boundVector))
    #     numpy.savetxt(generateNameST('bioT', self.currentSampleNumber(), self.currentTimeStep()),
    #                   # Added + '.numpy.txt'
    #                   numpy.array(gamma))
    #
    #     # spatial
    #     dist, gamma = generalfunctions_test01.experimentalVariogramValues(stackOfMapsAsListVariable,
    #                                                                       boundVector, 1, 1,
    #                                                                       'test', 2.0)
    #     numpy.savetxt(generateNameST('bioS', self.currentSampleNumber(), self.currentTimeStep()),
    #                   # Added + '.numpy.txt'
    #                   numpy.array(gamma))
    #
    #     # mean and var ### ADDITION - KL ###
    #     # MeanVarVariable = generalfunctions_test01.descriptiveStatistics(stackOfMapsAsListVariable)
    #     # numpy.savetxt(generateNameST('biMV', self.currentSampleNumber(), self.currentTimeStep()),
    #     #               numpy.array(MeanVarVariable))
    #     #
    #     # # lag-1 autocorrelation
    #     # lag1Variable = generalfunctions_test01.autocor1(stackOfMapsAsListVariable)
    #     # numpy.savetxt(generateNameST('biLO', self.currentSampleNumber(), self.currentTimeStep()),
    #     #               numpy.array(lag1Variable))
    #     # END OF ADDITION ###

    # # mean
    # meanVariable = areaaverage(variable, self.zoneMap)
    # generalfunctions_test01.reportLocationsAsNumpyArray(self.aLocation, meanVariable, 'bioA',
    #                                                     self.currentSampleNumber(), self.currentTimeStep())
    #
    # generalfunctions_test01.reportLocationsAsNumpyArray(self.aLocation, variable, 'bioF',
    #                                                     self.currentSampleNumber(), self.currentTimeStep())

# SOIL MOISTURE
self.d_subsurfaceWaterOneLayer.calculateSoilMoistureFraction()
variable = self.d_subsurfaceWaterOneLayer.soilMoistureFraction
variableSampled = ifthen(self.someLocs, variable)

self.historyOfSoilMoistureFraction = generalfunctions_test01.keepHistoryOfMaps(
    self.historyOfSoilMoistureFraction,
    variableSampled,
    self.durationHistory)
stackOfMapsAsListVariable = list(self.historyOfSoilMoistureFraction)

# # net weathering # - Testcase
# variable = self.d_bedrockweathering.weatheringMetrePerYear
# if save_np_temporal_mean == True:
#     generalfunctions_test01.report_locations_as_mean_np(
#         variable, 'weaA', self.currentSampleNumber(), self.currentTimeStep())
# if save_np_spatial_snapshots == True:
#     generalfunctions_test01.report_locations_as_np_arr(
#         variable, 'weaN', self.currentSampleNumber(), self.currentTimeStep())
# if save_maps == True:
#     generalfunctions_test01.report_as_map(
#         variable, 'weaM', self.currentSampleNumber(), self.currentTimeStep())

def postmcloop(self):
    # import generalfunctions # not sure why this needs to be imported again
    # names = [] # ['gA', 'bioA', 'bioM', 'demA', 'regA', 'sfA', 'qA', 'gpA', 'grA', 'grNA', 'depA', 'weaA', 'creA']
    # for name in names:
    #     aVariable = generalfunctions_test01.openSamplesAndTimestepsAsNumpyArraysAsNumpyArray(
    #         name, range(1, cfg.nrOfSamples + 1), timeStepsWithStatsCalculated)
    #     numpy.save(name, aVariable)
    pass