import numpy as np
import math
import os
from cycler import cycler
import matplotlib.pyplot as plt
from datetime import datetime

import EWS_configuration as cfg
import EWSPy as ews
import EWS_StateVariables as ews_sv

relative_start_of_grazing = cfg.rel_start_grazing
relative_end_of_grazing = 1

grazingRateIncreaseTotal = cfg.tot_increase_grazing
grazingRateIncrease = \
grazingRateIncreaseTotal / ((cfg.number_of_timesteps_weekly * (1 - relative_start_of_grazing)) / 2)

if cfg.return_ini_grazing:
    if self.currentTimeStep() < (cfg.number_of_timesteps_weekly * relative_start_of_grazing):
        self.grazingRate = self.grazingRate  # Equal to the initial grazing rate (see above)
    elif self.currentTimeStep() < (cfg.number_of_timesteps_weekly / ((1 - relative_start_of_grazing) / 2)):
        self.grazingRate = self.grazingRate + grazingRateIncrease
    else:
        self.grazingRate = self.grazingRate - grazingRateIncrease

if not cfg.return_ini_grazing:
    if self.currentTimeStep() < (cfg.number_of_timesteps_weekly * relative_start_of_grazing):
        self.grazingRate = self.grazingRate  # Equal to the initial grazing rate (see above)
    else:
        self.grazingRate = self.grazingRate + (grazingRateIncrease / 2)



relative_start_of_grazing = 0.1
relative_end_of_grazing = 0.9
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
