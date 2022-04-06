import datetime
from collections import deque
import sys
import numpy as np
import EWS_configuration as cfg
import time as timeit

sys.path.append("./pcrasterModules/")

from pcrasterModules import generalfunctions_test01
from pcrasterModules import datetimePCRasterPython
from pcrasterModules import interceptionuptomaxstore
from pcrasterModules import surfacestore
from pcrasterModules import infiltrationonlyksat
from pcrasterModules import subsurfacewateronelayer
from pcrasterModules import runoffaccuthreshold
from pcrasterModules import rainfalleventsfromgammadistribution
from pcrasterModules import exchangevariables_weekly
from pcrasterModules import soilwashMMF
from pcrasterModules import regolith
from pcrasterModules import bedrockweathering
from pcrasterModules import evapotranspirationsimple
from pcrasterModules import biomassmodifiedmay
from pcrasterModules import baselevel
from pcrasterModules import creep

from pcraster.framework import *

if cfg.fixedStates:
    fixedStatesReg = spatial(scalar(float(sys.argv[1])))
    fixedStatesBio = spatial(scalar(float(sys.argv[2])))


def calculateGapFractionAndMaxIntStoreFromLAI(leafAreaIndex):
    maximumInterceptionCapacityPerLAI = scalar(0.001)
    gapFraction = exp(-0.5 * leafAreaIndex)  # equation 40 in Brolsma et al 2010a
    maximumInterceptionStore = maximumInterceptionCapacityPerLAI * leafAreaIndex
    return gapFraction, maximumInterceptionStore


class CatchmentModel(DynamicModel, MonteCarloModel):

    def __init__(self):
        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)
        setclone('inputs_weekly/clone.map')

        # fix the seed for random functions (pcraster)
        setrandomseed(101)

        if cfg.filtering:
            ParticleFilterModel.__init__(self)

    def premcloop(self):

        self.clone = boolean("inputs_weekly/clone.map")
        self.numberOfCellsOnMap = maptotal(ifthenelse(self.clone, scalar(1), scalar(1)))

        # Locations where values are reported as a numpy array to disk.
        # In this model iteration, all locations are reported. Hence, zoning and sample locations are also omitted.
        self.all_locs = nominal("inputs_weekly/clone.map")
        self.locations2report = pcr2numpy(readmap('./inputs_weekly/mlocs.map'), 0).astype(bool)

        self.createInstancesPremcloop()

        # Duration history is used to store the data of previous timesteps, in this model iteration the number of
        # timesteps reported in the final timeseries is also dependant on this duration history.
        self.durationHistory = cfg.number_of_timesteps_weekly

        # time step duration in hours, typically (and only tested) one week, i.e. 7.0*24.0
        self.timeStepDuration = 7.0 * 24.0

    def initial(self):

        self.initializeTime(cfg.startTimeYearValue, cfg.startTimeMonthValue, cfg.startTimeDayValue,
                            self.timeStepDuration)

        self.createInstancesInitial()

        self.d_exchangevariables.upwardSeepageFlux = scalar(0)
        self.d_exchangevariables.evapFromSoilMultiplier = scalar(1)

        self.timeStepDurationYears = self.timeStepDuration / (365.0 * 24.0)

        self.actualAbstractionFluxFromSubsurface = 0.0

        # functions and settings for saving timeseries
        self.historyOfSoilMoistureFractionMean = deque([])
        self.historyOfSoilMoistureFractionLoc = deque([])
        self.historyOfBiomassMean = deque([])
        self.historyOfBiomassLoc = deque([])
        self.historyOfRegolithThicknessMean = deque([])
        self.historyOfRegolithThicknessLoc = deque([])
        self.historyOfDemMean = deque([])
        self.historyOfDemLoc = deque([])
        self.historyOfTotQ = deque([])
        self.history_of_grazing_rate = deque([])
        self.history_of_growth_part = deque([])
        self.history_of_grazing_part = deque([])
        self.history_of_net_growth = deque([])
        self.history_of_net_deposition_mean = deque([])
        self.history_of_net_deposition_loc = deque([])
        self.history_of_net_weathering_mean = deque([])
        self.history_of_net_weathering_loc = deque([])
        self.history_of_net_creep_deposition_mean = deque([])
        self.history_of_net_creep_deposition_loc = deque([])
        self.history_of_max_int_mean = deque([])
        self.history_of_max_int_loc = deque([])
        self.history_of_LAI_mean = deque([])
        self.history_of_LAI_loc = deque([])

        # budgets
        self.d_exchangevariables.cumulativePrecipitation = scalar(0)

        # initial values
        self.grazingRate = cfg.initial_grazing
        self.runoffMetreWaterDepthPerHour = scalar(0.0)
        self.creepDeposition = spatial(scalar(0.0))

    def dynamic(self):
        # biomass
        if cfg.fixedStates:
            self.d_regolithdemandbedrock.setNewRegolith(spatial(scalar(fixedStatesReg)))
            self.d_biomassModifiedMay.setNewBiomass(spatial(scalar(fixedStatesBio)))

        # grazing pressure driver
        # # Static at first to set a baseline for the first state, increase and decrease afterwards (same timespan)
        # relative_start_of_grazing = cfg.rel_start_grazing
        # grazingRateIncreaseTotal = cfg.tot_increase_grazing
        # grazingRateIncrease = \
        #     grazingRateIncreaseTotal / ((cfg.number_of_timesteps_weekly * (1 - relative_start_of_grazing)) / 2)
        #
        # if cfg.return_ini_grazing:
        #     if self.currentTimeStep() < (cfg.number_of_timesteps_weekly * relative_start_of_grazing):
        #         self.grazingRate = self.grazingRate  # Equal to the initial grazing rate (see above)
        #     elif self.currentTimeStep() < (cfg.number_of_timesteps_weekly / ((1 - relative_start_of_grazing) / 2)):
        #         self.grazingRate = self.grazingRate + grazingRateIncrease
        #     else:
        #         self.grazingRate = self.grazingRate - grazingRateIncrease
        #
        # if not cfg.return_ini_grazing:
        #     if self.currentTimeStep() < (cfg.number_of_timesteps_weekly * relative_start_of_grazing):
        #         self.grazingRate = self.grazingRate  # Equal to the initial grazing rate (see above)
        #     else:
        #         self.grazingRate = self.grazingRate + (grazingRateIncrease / 2)

        relative_start_of_grazing = cfg.rel_start_grazing
        relative_end_of_grazing = cfg.rel_end_grazing
        grazingRateTime = (relative_end_of_grazing - relative_start_of_grazing)
        grazingRateIncrease = cfg.tot_increase_grazing / ((cfg.number_of_timesteps_weekly * grazingRateTime) / 2)

        if self.currentTimeStep() >= (cfg.number_of_timesteps_weekly * relative_start_of_grazing):
            if self.currentTimeStep() < (cfg.number_of_timesteps_weekly * relative_end_of_grazing):
                if not cfg.return_ini_grazing:
                    self.grazingRate += (grazingRateIncrease / 2)
                if cfg.return_ini_grazing:
                    if self.currentTimeStep < ((cfg.number_of_timesteps_weekly * relative_start_of_grazing) + ((cfg.number_of_timesteps_weekly * grazingRateTime) / 2) ):
                        self.grazingRate += grazingRateIncrease
                    else:
                        self.grazingRate -= grazingRateIncrease

        runoffMetreWaterDepthPerWeek = self.runoffMetreWaterDepthPerHour * cfg.rainstorm_duration
        self.biomass, self.LAI = self.d_biomassModifiedMay.update(self.actualAbstractionFluxFromSubsurface,
                                                                  runoffMetreWaterDepthPerWeek, self.grazingRate)

        # update gap fraction and maximum interception store
        gapFraction, maxIntStore = calculateGapFractionAndMaxIntStoreFromLAI(self.LAI)
        self.d_interceptionuptomaxstore.setGapFraction(gapFraction)
        self.d_interceptionuptomaxstore.setMaximumStore(maxIntStore)

        # update stone cover for erosion
        fractionOfVegetationAboveSoil = 0.7
        vegetationCoverForErosion = (1.0 - gapFraction) * fractionOfVegetationAboveSoil
        self.d_soilwashMMF.updateStoneOrVegetationCover(vegetationCoverForErosion)

        # precipitation
        isRaining, rainfallFlux, rainfallAmount = self.d_rainfalleventsfromgammadistribution.getRainstorm()

        self.d_exchangevariables.cumulativePrecipitation = \
            self.d_exchangevariables.cumulativePrecipitation + rainfallFlux * self.timeStepDuration

        if isRaining:
            # interception store
            actualAdditionFluxToInterceptionStore = self.d_interceptionuptomaxstore.addWater(rainfallFlux)
            throughfallFlux = rainfallFlux - actualAdditionFluxToInterceptionStore

            # surface store
            # Currently no upward seepage
            # totalToSurfaceFlux = throughfallFlux + self.d_exchangevariables.upwardSeepageFlux
            totalToSurfaceFlux = throughfallFlux
            potentialToSurfaceStoreFlux = self.d_surfaceStore.potentialToFlux()

            # potential infiltration
            self.d_infiltrationonlyksat.setSaturatedConductivityFluxAsFunctionOfBiomass(self.biomass)
            potentialHortonianInfiltrationFlux = self.d_infiltrationonlyksat.potentialInfiltrationFluxFunction()
            maximumSaturatedOverlandFlowInfiltrationFlux = self.d_subsurfaceWaterOneLayer.getMaximumAdditionFlux()
            potentialInfiltrationFlux = min(potentialHortonianInfiltrationFlux,
                                            maximumSaturatedOverlandFlowInfiltrationFlux)

            # abstraction from surface water
            potentialAbstractionFromSurfaceWaterFlux = potentialToSurfaceStoreFlux + potentialInfiltrationFlux
            actualAbstractionFromSurfaceWaterFlux, runoffCubicMetrePerHour = self.d_runoffAccuthreshold.update(
                totalToSurfaceFlux, potentialAbstractionFromSurfaceWaterFlux)
            potentialOutSurfaceStoreFlux = self.d_surfaceStore.potentialOutFlux()

            # infiltration
            availableForInfiltrationFlux = potentialOutSurfaceStoreFlux + actualAbstractionFromSurfaceWaterFlux
            availableForInfiltrationNotExceedingMaximumSaturatedOverlandFlowFlux = min(
                availableForInfiltrationFlux, maximumSaturatedOverlandFlowInfiltrationFlux)
            actualInfiltrationFlux = self.d_infiltrationonlyksat.update(
                availableForInfiltrationNotExceedingMaximumSaturatedOverlandFlowFlux)

            # surface store
            surfaceStoreChange = actualAbstractionFromSurfaceWaterFlux - actualInfiltrationFlux
            self.d_surfaceStore.update(surfaceStoreChange)
            actualAdditionFlux = self.d_subsurfaceWaterOneLayer.addWater(actualInfiltrationFlux)
            # empty surface store again
            self.d_surfaceStore.emptyIt()

            # surface wash
            self.runoffMetreWaterDepthPerHour = runoffCubicMetrePerHour / cellarea()
            netDeposition, netDepositionMetre, lateralFluxKg, totalDetachKgPerCell, transportCapacityKgPerCell = \
                self.d_soilwashMMF.calculateWash(self.runoffMetreWaterDepthPerHour, rainfallFlux, throughfallFlux)

            # TODO - Work out the statement below:
            # LET OP: dit is metre flux, maar dat zou het zijn als er slechts 1 regenbui is per jaar
            # het is dus de ene week uitgemiddeld over een jaar
            # om echt iets te krijgen met een eenheid m/jaar (dwz wat 'zou' de depositie zijn
            # als dit event elke week zou optreden), moet dit keer 52 weken
            # hetzelfde geldt voor actual deposition flux hieronder
            netDepositionMetreFlux = netDepositionMetre / self.timeStepDurationRegolithInYears

            # LDD, surface
            if cfg.changeGeomorphology:
                actualDepositionFlux = self.d_regolithdemandbedrock.updateWithDeposition(netDepositionMetreFlux)
                regolithThickness, demOfBedrock, dem, bedrockLdd, surfaceLdd = \
                    self.d_regolithdemandbedrock.getRegolithProperties()
                amountOfMoistureThickNetAdded = \
                    self.d_subsurfaceWaterOneLayer.updateRegolithThickness(regolithThickness)
                self.d_soilwashMMF.setSurfaceProperties(surfaceLdd, dem)
                self.d_runoffAccuthreshold.setSurfaceProperties(surfaceLdd)

        else:
            # surface wash
            netDeposition, netDepositionMetre, lateralFluxKg, totalDetachKgPerCell, transportCapacityKgPerCell = \
                self.d_soilwashMMF.noWash()
            actualDepositionFlux = spatial(scalar(0))
            self.runoffMetreWaterDepthPerHour = scalar(0)

        if cfg.changeGeomorphology:
            # random noise
            netDepositionMetreNoiseFlux = normal(1) / 5000
            # LDD, surface
            actualDepositionNoiseFlux = self.d_regolithdemandbedrock.updateWithDeposition(netDepositionMetreNoiseFlux)
            regolithThickness, demOfBedrock, dem, bedrockLdd, surfaceLdd = self.d_regolithdemandbedrock.getRegolithProperties()

        # potential evapotranspiration, m/hour
        fWaterPotential = self.d_subsurfaceWaterOneLayer.getFWaterPotential()
        potentialEvapotranspirationFlux = self.d_evapotranspirationSimple.potentialEvapotranspiration(fWaterPotential,
                                                                                                      self.biomass)
        # evapotranspirate first from interception store
        # assume this does not depend on vegetation, and does not influence transpiration
        # assume it immediately empties (ie, within a week)
        potentialEvaporationFromInterceptionStore = 99999.9
        actualAbstractionFluxFromInterceptionStore = self.d_interceptionuptomaxstore.abstractWater(
            potentialEvaporationFromInterceptionStore)

        # evapotranspirate from subsurface store
        potentialEvapotranspirationFluxFromSubsurface = \
            max(0.0, potentialEvapotranspirationFlux)

        self.actualAbstractionFluxFromSubsurface = \
            self.d_subsurfaceWaterOneLayer.abstractWater(potentialEvapotranspirationFluxFromSubsurface)

        # lateral flow in subsurface and upward seepage from subsurface storage
        # typically switched off and never tested
        # # self.d_exchangevariables.upwardSeepageFlux=self.d_subsurfaceWaterOneLayer.lateralFlow()

        # self.checkBudgets(self.currentSampleNumber(), self.currentTimeStep())
        # self.printMemberVariables()

        ###############
        # geomorphology
        ###############
        if cfg.changeGeomorphology and (self.currentTimeStep() % 52 == 0):
            # bedrock weathering
            regolithThickness, demOfBedrock, dem, bedrockLdd, surfaceLdd = self.d_regolithdemandbedrock.getRegolithProperties()
            bedrockWeatheringFlux = self.d_bedrockweathering.weatheringRate(regolithThickness)
            # LDD, bedrock
            self.d_regolithdemandbedrock.updateWithBedrockWeathering(bedrockWeatheringFlux)

            # creep
            regolithThickness, demOfBedrock, dem, bedrockLdd, surfaceLdd = self.d_regolithdemandbedrock.getRegolithProperties()
            newRegolithThickness, outflow, flowOverBoundaries, correctedFactor, amountX, amountY, inflowX, inflowY = \
                self.d_creep.diffuse(regolithThickness, dem, 1)
            self.creepDeposition = newRegolithThickness - regolithThickness

            # LDD, surface
            self.d_regolithdemandbedrock.setNewRegolith(newRegolithThickness)

            # update bedrock with baselevel change
            baselevel = self.d_baselevel.getBaselevel(self.currentTimeStep())
            # LDD, surface, bedrock
            self.d_regolithdemandbedrock.setBaselevel(baselevel)

            # update subsurface store with new regolith thickness
            regolithThickness, demOfBedrock, dem, bedrockLdd, surfaceLdd = self.d_regolithdemandbedrock.getRegolithProperties()
            amountOfMoistureThickNetAdded = self.d_subsurfaceWaterOneLayer.updateRegolithThickness(regolithThickness)
            # no lateral flow, so bedrock does not need to be updated
            # self.d_subsurfaceWaterOneLayer=updateBedrock(self,bedRockLdd,demOfBedrock)

            # update soil wash and runoff with new surface properties
            self.d_soilwashMMF.setSurfaceProperties(surfaceLdd, dem)
            self.d_runoffAccuthreshold.setSurfaceProperties(surfaceLdd)

        # reports
        self.reportComponentsDynamic()
        save_maps = (self.currentTimeStep() % cfg.interval_map_snapshots) == 0 and cfg.map_data
        save_mean_timeseries = self.currentTimeStep() == cfg.number_of_timesteps_weekly and cfg.mean_timeseries_data
        save_loc_timeseries = self.currentTimeStep() == cfg.number_of_timesteps_weekly and cfg.loc_timeseries_data

        ###########
        # Variable reporting
        ###########

        # MAXIMUM INTERCEPTION STORE
        variable = maxIntStore
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))
        variable_loc = pcr2numpy(variable, np.NaN)[self.locations2report]

        self.history_of_max_int_mean = generalfunctions_test01.keepHistoryOfMaps(self.history_of_max_int_mean,
                                                                                 variable_mean, self.durationHistory)
        self.history_of_max_int_loc = generalfunctions_test01.keepHistoryOfMaps(self.history_of_max_int_loc,
                                                                                variable_loc, self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'micM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.history_of_max_int_mean)
            generalfunctions_test01.report_as_array(variable_mean_array, 'micA', self.currentSampleNumber(),
                                                    self.currentTimeStep())
        if save_loc_timeseries:
            variable_loc_array = np.array(self.history_of_max_int_loc)
            generalfunctions_test01.report_as_array(variable_loc_array, 'micL', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # LAI
        variable = self.LAI
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))
        variable_loc = pcr2numpy(variable, np.NaN)[self.locations2report]

        self.history_of_LAI_mean = generalfunctions_test01.keepHistoryOfMaps(self.history_of_LAI_mean, variable_mean,
                                                                             self.durationHistory)
        self.history_of_LAI_loc = generalfunctions_test01.keepHistoryOfMaps(self.history_of_LAI_loc, variable_loc,
                                                                            self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'laiM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.history_of_LAI_mean)
            generalfunctions_test01.report_as_array(variable_mean_array, 'laiA', self.currentSampleNumber(),
                                                    self.currentTimeStep())
        if save_loc_timeseries:
            variable_loc_array = np.array(self.history_of_LAI_loc)
            generalfunctions_test01.report_as_array(variable_loc_array, 'laiL', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # SOIL MOISTURE
        self.d_subsurfaceWaterOneLayer.calculateSoilMoistureFraction()
        variable = self.d_subsurfaceWaterOneLayer.soilMoistureFraction
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))
        variable_loc = pcr2numpy(variable, np.NaN)[self.locations2report]

        self.historyOfSoilMoistureFractionMean = generalfunctions_test01.keepHistoryOfMaps(
            self.historyOfSoilMoistureFractionMean, variable_mean, self.durationHistory)
        self.historyOfSoilMoistureFractionLoc = generalfunctions_test01.keepHistoryOfMaps(
            self.historyOfSoilMoistureFractionLoc, variable_loc, self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'moiM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.historyOfSoilMoistureFractionMean)
            generalfunctions_test01.report_as_array(variable_mean_array, 'moiA', self.currentSampleNumber(),
                                                    self.currentTimeStep())
        if save_loc_timeseries:
            variable_loc_array = np.array(self.historyOfSoilMoistureFractionLoc)
            generalfunctions_test01.report_as_array(variable_loc_array, 'moiL', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # BIOMASS
        variable = self.d_biomassModifiedMay.biomass
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))
        variable_loc = pcr2numpy(variable, np.NaN)[self.locations2report]

        self.historyOfBiomassMean = generalfunctions_test01.keepHistoryOfMaps(self.historyOfBiomassMean, variable_mean,
                                                                          self.durationHistory)
        self.historyOfBiomassLoc = generalfunctions_test01.keepHistoryOfMaps(self.historyOfBiomassLoc, variable_loc,
                                                                          self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'bioM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.historyOfBiomassMean)
            generalfunctions_test01.report_as_array(variable_mean_array, 'bioA', self.currentSampleNumber(),
                                                    self.currentTimeStep())
        if save_loc_timeseries:
            variable_loc_array = np.array(self.historyOfBiomassLoc)
            generalfunctions_test01.report_as_array(variable_loc_array, 'bioL', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # REGOLITH THICKNESS
        variable = self.d_regolithdemandbedrock.regolithThickness
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))
        variable_loc = pcr2numpy(variable, np.NaN)[self.locations2report]

        self.historyOfRegolithThicknessMean = generalfunctions_test01.keepHistoryOfMaps(
            self.historyOfRegolithThicknessMean, variable_mean, self.durationHistory)
        self.historyOfRegolithThicknessLoc = generalfunctions_test01.keepHistoryOfMaps(
            self.historyOfRegolithThicknessLoc, variable_loc, self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'regM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.historyOfRegolithThicknessMean)
            generalfunctions_test01.report_as_array(variable_mean_array, 'regA', self.currentSampleNumber(),
                                                    self.currentTimeStep())
        if save_loc_timeseries:
            variable_loc_array = np.array(self.historyOfRegolithThicknessLoc)
            generalfunctions_test01.report_as_array(variable_loc_array, 'regL', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # DEM
        variable = self.d_regolithdemandbedrock.dem
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))
        variable_loc = pcr2numpy(variable, np.NaN)[self.locations2report]

        self.historyOfDemMean = generalfunctions_test01.keepHistoryOfMaps(self.historyOfDemMean, variable_mean,
                                                                          self.durationHistory)
        self.historyOfDemLoc = generalfunctions_test01.keepHistoryOfMaps(self.historyOfDemLoc, variable_loc,
                                                                         self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'demM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.historyOfDemMean)
            generalfunctions_test01.report_as_array(variable_mean_array, 'demA', self.currentSampleNumber(),
                                                    self.currentTimeStep())
        if save_loc_timeseries:
            variable_loc_array = np.array(self.historyOfDemLoc)
            generalfunctions_test01.report_as_array(variable_loc_array, 'demL', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # Discharge
        # ! - Note that although RunoffCubicMetrePerHour is used, it returns the runoff per time step.
        downstreamEdge = generalfunctions_test01.edge(self.clone, 2, 0)
        pits = pcrne(pit(self.d_runoffAccuthreshold.ldd), 0)
        outflowPoints = pcrand(downstreamEdge, pits)
        totQ = ifthen(self.clone, maptotal(ifthenelse(outflowPoints, self.d_runoffAccuthreshold.RunoffCubicMetrePerHour, scalar(0))))

        variable = totQ
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))

        self.historyOfTotQ = generalfunctions_test01.keepHistoryOfMaps(self.historyOfTotQ, variable_mean,
                                                                       self.durationHistory)

        if save_mean_timeseries:
            variable_mean_array = np.array(self.historyOfTotQ)
            generalfunctions_test01.report_as_array(variable_mean_array, 'qA', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # grazing rate
        variable = self.grazingRate

        self.history_of_grazing_rate = generalfunctions_test01.keepHistoryOfMaps(self.history_of_grazing_rate,
                                                                                 variable,
                                                                                 self.durationHistory)

        if save_mean_timeseries:
            variable_mean_array = np.array(self.history_of_grazing_rate)
            generalfunctions_test01.report_as_array(variable_mean_array, 'gA', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # growth part
        variable = self.d_biomassModifiedMay.growthPart
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))

        self.history_of_growth_part = generalfunctions_test01.keepHistoryOfMaps(self.history_of_growth_part,
                                                                                variable_mean,
                                                                                self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'gpM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.history_of_growth_part)
            generalfunctions_test01.report_as_array(variable_mean_array, 'gpA', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # grazing part
        variable = 0.0 - self.d_biomassModifiedMay.grazing
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))

        self.history_of_grazing_part = generalfunctions_test01.keepHistoryOfMaps(self.history_of_grazing_part,
                                                                                 variable_mean, self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'grM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.history_of_grazing_part)
            generalfunctions_test01.report_as_array(variable_mean_array, 'grA', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # net growth
        variable = self.d_biomassModifiedMay.netGrowth
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))

        self.history_of_net_growth = generalfunctions_test01.keepHistoryOfMaps(self.history_of_net_growth,
                                                                               variable_mean, self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'grNM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.history_of_net_growth)
            generalfunctions_test01.report_as_array(variable_mean_array, 'grNA', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # net deposition
        variable = actualDepositionFlux
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))
        variable_loc = pcr2numpy(variable, np.NaN)[self.locations2report]

        self.history_of_net_deposition_mean = generalfunctions_test01.keepHistoryOfMaps(
            self.history_of_net_deposition_mean, variable_mean, self.durationHistory)
        self.history_of_net_deposition_loc = generalfunctions_test01.keepHistoryOfMaps(
            self.history_of_net_deposition_loc, variable_loc, self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'depM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.history_of_net_deposition_mean)
            generalfunctions_test01.report_as_array(variable_mean_array, 'depA', self.currentSampleNumber(),
                                                    self.currentTimeStep())
        if save_loc_timeseries:
            variable_loc_array = np.array(self.history_of_net_deposition_loc)
            generalfunctions_test01.report_as_array(variable_loc_array, 'depL', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # net weathering
        variable = self.d_bedrockweathering.weatheringMetrePerYear
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))
        variable_loc = pcr2numpy(variable, np.NaN)[self.locations2report]

        self.history_of_net_weathering_mean = generalfunctions_test01.keepHistoryOfMaps(
            self.history_of_net_weathering_mean, variable_mean, self.durationHistory)
        self.history_of_net_weathering_loc = generalfunctions_test01.keepHistoryOfMaps(
            self.history_of_net_weathering_loc, variable_loc, self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'weaM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.history_of_net_weathering_mean)
            generalfunctions_test01.report_as_array(variable_mean_array, 'weaA', self.currentSampleNumber(),
                                                    self.currentTimeStep())
        if save_loc_timeseries:
            variable_loc_array = np.array(self.history_of_net_weathering_loc)
            generalfunctions_test01.report_as_array(variable_loc_array, 'weaL', self.currentSampleNumber(),
                                                    self.currentTimeStep())

        # net creep deposition
        variable = self.creepDeposition
        variable_mean = np.nanmean(pcr2numpy(variable, np.NaN))
        variable_loc = pcr2numpy(variable, np.NaN)[self.locations2report]

        self.history_of_net_creep_deposition_mean = generalfunctions_test01.keepHistoryOfMaps(
            self.history_of_net_creep_deposition_mean, variable_mean, self.durationHistory)
        self.history_of_net_creep_deposition_loc = generalfunctions_test01.keepHistoryOfMaps(
            self.history_of_net_creep_deposition_loc, variable_loc, self.durationHistory)

        if save_maps:
            generalfunctions_test01.report_as_map(variable, 'creM', self.currentSampleNumber(), self.currentTimeStep())
        if save_mean_timeseries:
            variable_mean_array = np.array(self.history_of_net_creep_deposition_mean)
            generalfunctions_test01.report_as_array(variable_mean_array, 'creA', self.currentSampleNumber(),
                                                    self.currentTimeStep())
        if save_loc_timeseries:
            variable_loc_array = np.array(self.history_of_net_weathering_loc)
            generalfunctions_test01.report_as_array(variable_loc_array, 'creL', self.currentSampleNumber(),
                                                    self.currentTimeStep())

    def postmcloop(self):
        pass

    def createInstancesPremcloop(self):
        pass

    def createInstancesInitial(self):
        timeStepsToReportAll = cfg.timesteps_to_report_all_weekly
        timeStepsToReportSome = cfg.timesteps_to_report_some_weekly

        # Class for exchange variables in initial and dynamic, introduced to make filtering possible
        self.d_exchangevariables = exchangevariables_weekly.ExchangeVariables(timeStepsToReportSome,
                                                                              cfg.exchange_report_rasters)

        # base level
        deterministicDem = scalar('inputs_weekly/demini.map')
        dem = deterministicDem
        baselevelRise = -0.0001
        self.d_baselevel = baselevel.Baselevel(
            generalfunctions_test01.bottom(self.clone),
            deterministicDem,
            baselevelRise,
            self.timeStepDuration / (365.0 * 24.0),
            timeStepsToReportAll,
            cfg.baselevel_report_rasters)

        weatheringRateBareBedrock = 0.0005
        weatheringExponentParameter = 4.0
        self.d_bedrockweathering = bedrockweathering.BedrockWeathering(
            weatheringRateBareBedrock,
            weatheringExponentParameter,
            timeStepsToReportAll,
            cfg.bedrockweathering_report_rasters)
        steadyStateSoilDepth = self.d_bedrockweathering.steadyStateSoilDepth(0 - baselevelRise)
        self.report(steadyStateSoilDepth, 'sssd')

        # regolith
        regolithThickness = spatial(steadyStateSoilDepth)

        self.timeStepDurationRegolithInYears = 1.0
        self.d_regolithdemandbedrock = regolith.RegolithDemAndBedrock(
            dem,
            regolithThickness,
            self.timeStepDurationRegolithInYears,
            timeStepsToReportAll,
            cfg.regolith_report_rasters)

        regolithThickness, demOfBedrock, dem, bedrockLdd, surfaceLdd = self.d_regolithdemandbedrock.getRegolithProperties()
        self.report(regolithThickness, 'regIni')  # Added self., removed .map

        # biomass
        initialBiomass = 2.0
        waterUseEfficiency = cfg.waterUseEfficiency
        maintenanceRate = cfg.maintenanceRate

        # original
        gamma = 0.004  # runoff
        # # gamma zero
        # gamma=0.0001 # runoff
        # # gamma low
        # gamma=0.001 # runoff
        # # gamma medium
        # gamma=0.002 # runoff

        alpha = 0.4  # grazing
        dispersion = 0.01 / (365.0 * 24)
        runoff = 0.0
        sdOfNoise = 0.000000000001
        LAIPerBiomass = 2.5
        self.d_biomassModifiedMay = biomassmodifiedmay.BiomassModifiedMay(
            initialBiomass,
            waterUseEfficiency,
            maintenanceRate,
            gamma,
            alpha,
            dispersion,
            runoff,
            sdOfNoise,
            LAIPerBiomass,
            self.timeStepDuration,
            timeStepsToReportAll,
            cfg.biomassmodifiedmay_report_rasters)

        # precipitation
        probabilityOfARainstorm = cfg.rainstorm_probability
        durationOfRainstorm = cfg.rainstorm_duration
        expectedRainfallIntensity = cfg.rainstorm_expected_intensity
        gammaShapeParameter = cfg.rainstorm_gamma_shape_param

        self.d_rainfalleventsfromgammadistribution = \
            rainfalleventsfromgammadistribution.RainfallEventsFromGammaDistribution(
                probabilityOfARainstorm,
                durationOfRainstorm,
                expectedRainfallIntensity,
                gammaShapeParameter,
                timeStepsToReportAll,
                cfg.rainfalleventsfromgammadistribution_report_rasters)

        # interception
        initialLeafAreaIndex = initialBiomass * LAIPerBiomass
        initialInterceptionStore = scalar(0.000001)
        gapFraction, maximumInterceptionStore = calculateGapFractionAndMaxIntStoreFromLAI(initialLeafAreaIndex)

        self.d_interceptionuptomaxstore = interceptionuptomaxstore.InterceptionUpToMaxStore(
            spatial(ldd(5)),
            initialInterceptionStore,
            maximumInterceptionStore,
            gapFraction,
            cfg.calculateUpstreamTotals,
            durationOfRainstorm,
            timeStepsToReportAll,
            cfg.interception_report_rasters)

        # surface store
        initialSurfaceStore = scalar(0.0)
        maxSurfaceStore = scalar(cfg.maxSurfaceStoreValue)
        self.d_surfaceStore = surfacestore.SurfaceStore(
            initialSurfaceStore,
            maxSurfaceStore,
            durationOfRainstorm,
            timeStepsToReportAll,
            cfg.surfacestore_report_rasters)

        # infiltration
        bareSoilSaturatedConductivityFlux = scalar(0.0001)
        # maxSaturatedConductivityFluxFromVegetation=scalar(0.01)
        maxSaturatedConductivityFluxFromVegetation = scalar(0.1)
        biomassHalfSaturation = scalar(1.0)
        ksat = bareSoilSaturatedConductivityFlux
        self.d_infiltrationonlyksat = infiltrationonlyksat.InfiltrationOnlyKsat(
            ksat,
            bareSoilSaturatedConductivityFlux,
            maxSaturatedConductivityFluxFromVegetation,
            biomassHalfSaturation,
            durationOfRainstorm,
            timeStepsToReportAll,
            cfg.infiltration_report_rasters_weekly)

        # subsurface water
        initialSoilMoistureFraction = scalar(cfg.initialSoilMoistureFractionCFG)

        soilPorosityFraction = scalar(cfg.soilPorosityFractionValue)
        fieldCapacityFraction = scalar(cfg.fieldCapacityFractionValue)
        limitingPointFraction = scalar(cfg.limitingPointFractionValue)
        wiltingPointFraction = scalar(cfg.mergeWiltingPointFractionFSValue)

        saturatedConductivityMetrePerDay = generalfunctions_test01.mapuniformBounds(
            2, 8, scalar(12.5), cfg.createRealizations)

        self.d_subsurfaceWaterOneLayer = subsurfacewateronelayer.SubsurfaceWaterOneLayer(
            bedrockLdd,
            demOfBedrock,
            regolithThickness,
            initialSoilMoistureFraction,
            soilPorosityFraction,
            wiltingPointFraction,
            fieldCapacityFraction,
            limitingPointFraction,
            saturatedConductivityMetrePerDay,
            cfg.calculateUpstreamTotals,
            self.timeStepDurationHours,
            timeStepsToReportAll,
            cfg.subsurface_report_rasters)

        # evapotranspiration
        beta = 1.0
        maximumEvapotranspirationFlux = 0.8 / (365.0 * 24.0)
        self.d_evapotranspirationSimple = evapotranspirationsimple.EvapotranspirationSimple(
            self.timeStepDuration,
            beta,
            maximumEvapotranspirationFlux,
            timeStepsToReportAll,
            cfg.evapotranspirationsimple_report_rasters)

        # runoff
        self.d_runoffAccuthreshold = runoffaccuthreshold.RunoffAccuthreshold(
            surfaceLdd,
            durationOfRainstorm,
            timeStepsToReportAll,
            cfg.runoff_report_rasters)

        # soilwash
        plantHeightMetres = 5.0
        stoneCoverFraction = 0.1
        vegetationCoverOfSoilFraction = 0.1
        manningsN = 0.03  # 'original'

        # standard erosion scenario
        detachabilityOfSoilRaindrops = 1.6  # 'original'  (used for all scenarios)
        detachabilityOfSoilRunoff = 6.4  # 'original'

        # # more erosion scenario
        # detachabilityOfSoilRaindrops=16
        # detachabilityOfSoilRunoff=64

        self.d_soilwashMMF = soilwashMMF.SoilWashMMF(
            surfaceLdd,
            dem,
            durationOfRainstorm,
            plantHeightMetres,
            detachabilityOfSoilRaindrops,
            stoneCoverFraction,
            detachabilityOfSoilRunoff,
            vegetationCoverOfSoilFraction,
            manningsN,
            soilPorosityFraction,
            timeStepsToReportAll,
            cfg.soilwashMMF_report_rasters)

        # creep
        diffusion = 0.01
        self.d_creep = creep.Creep(
            dem,
            self.timeStepDurationRegolithInYears,
            diffusion,
            timeStepsToReportAll,
            cfg.creep_report_rasters)

    def reportComponentsDynamic(self):
        components = [
            self.d_exchangevariables,
            self.d_evapotranspirationSimple,
            self.d_regolithdemandbedrock,
            self.d_bedrockweathering,
            self.d_baselevel,
            self.d_rainfalleventsfromgammadistribution,
            self.d_interceptionuptomaxstore,
            self.d_surfaceStore,
            self.d_infiltrationonlyksat,
            self.d_runoffAccuthreshold,
            self.d_subsurfaceWaterOneLayer,
            self.d_soilwashMMF,
            self.d_creep,
            self.d_biomassModifiedMay
        ]

        for component in components:
            component.reportAsMaps(self.currentSampleNumber(), self.currentTimeStep())

    def printMemberVariables(self):
        components = [
            self.d_exchangevariables,
            self.d_interceptionuptomaxstore,
            self.d_surfaceStore,
            self.d_infiltrationonlyksat,
            self.d_runoffAccuthreshold,
            self.d_subsurfaceWaterOneLayer
        ]

        for component in components:
            generalfunctions_test01.printMemberVariables(component)

    def printComponentsDynamic(self):
        self.d_dateTimePCRasterPython.printit()

    def initializeTime(self, startTimeYear, startTimeMonth, startTimeDay, timeStepDurationHours):
        startTime = datetime.datetime(year=startTimeYear, month=startTimeMonth, day=startTimeDay)
        self.timeStepDurationHours = timeStepDurationHours
        self.timeStepDatetimeFormat = datetime.timedelta(hours=self.timeStepDurationHours)
        self.d_dateTimePCRasterPython = datetimePCRasterPython.DatetimePCRasterPython(startTime,
                                                                                      self.timeStepDatetimeFormat)

    def checkBudgets(self, currentSampleNumber, currentTimeStep):

        increaseInPrecipitationStore = 0.0 - self.d_exchangevariables.cumulativePrecipitation
        report(increaseInPrecipitationStore, generateNameST('incP', currentSampleNumber, currentTimeStep))

        increaseInInterceptionStore = self.d_interceptionuptomaxstore.budgetCheck(currentSampleNumber, currentTimeStep)
        report(increaseInInterceptionStore, generateNameST('incI', currentSampleNumber, currentTimeStep))

        increaseInSurfaceStore = self.d_surfaceStore.budgetCheck(currentSampleNumber, currentTimeStep)
        report(increaseInSurfaceStore, generateNameST('incS', currentSampleNumber, currentTimeStep))
        increaseInSurfaceStoreQM = catchmenttotal(increaseInSurfaceStore, self.ldd) * cellarea()
        report(increaseInSurfaceStoreQM, generateNameST('testb', currentSampleNumber, currentTimeStep))

        # let op: infiltration store is directly passed to subsurface store, thus is not a real store
        increaseInInfiltrationStore = self.d_infiltrationonlyksat.budgetCheck(currentSampleNumber, currentTimeStep)

        increaseInSubSurfaceWaterStore, lateralFlowInSubsurfaceStore, abstractionFromSubSurfaceWaterStore = \
            self.d_subsurfaceWaterOneLayer.budgetCheck(currentSampleNumber, currentTimeStep)
        increaseInSubSurfaceStoreQM = catchmenttotal(increaseInSubSurfaceWaterStore, self.ldd) * cellarea()

        increaseInRunoffStoreCubicMetresInUpstreamArea = self.d_runoffAccuthreshold.budgetCheck()

        totalIncreaseInStoresCubicMetresInUpstreamArea = 0.0
        stores = [increaseInPrecipitationStore, increaseInInterceptionStore, increaseInSurfaceStore,
                  increaseInSubSurfaceWaterStore]
        for store in stores:
            increaseInStoreCubicMetresInUpstreamArea = catchmenttotal(store, self.ldd) * cellarea()
            totalIncreaseInStoresCubicMetresInUpstreamArea = totalIncreaseInStoresCubicMetresInUpstreamArea \
                + increaseInStoreCubicMetresInUpstreamArea

        report(totalIncreaseInStoresCubicMetresInUpstreamArea,
               generateNameST('inSt', currentSampleNumber, currentTimeStep))
        report(increaseInRunoffStoreCubicMetresInUpstreamArea,
               generateNameST('inRu', currentSampleNumber, currentTimeStep))
        report(catchmenttotal(self.d_exchangevariables.upwardSeepageFlux, self.ldd) * cellarea(),
               generateNameST('inSe', currentSampleNumber, currentTimeStep))
        # total budget is total increase in stores plus the upward seepage flux for each ts that is passed to the next
        # timestep and thus not taken into account in the current timestep budgets
        budget = totalIncreaseInStoresCubicMetresInUpstreamArea + increaseInRunoffStoreCubicMetresInUpstreamArea \
            + lateralFlowInSubsurfaceStore * cellarea() + catchmenttotal(abstractionFromSubSurfaceWaterStore, self.ldd) \
            * cellarea() + catchmenttotal(self.d_exchangevariables.upwardSeepageFlux, self.ldd) * cellarea()
        report(budget, generateNameST('B-tot', currentSampleNumber, currentTimeStep))
        budgetRel = budget / increaseInRunoffStoreCubicMetresInUpstreamArea
        report(budgetRel, generateNameST('B-rel', currentSampleNumber, currentTimeStep))

start_time = timeit.time()

myModel = CatchmentModel()
dynamicModel = DynamicFramework(myModel, cfg.number_of_timesteps_weekly)
mcModel = MonteCarloFramework(dynamicModel, cfg.nrOfSamples)
mcModel.setForkSamples(True, 10)
mcModel.run()

end_time = timeit.time() - start_time
print(f"Total elapsed time equals: {end_time} seconds")
