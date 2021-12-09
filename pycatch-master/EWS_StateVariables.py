import configuration_weekly as cfg

variables = []


# Class StateVariable for Variable objects #


class StateVariable:
    def __init__(self, name, spatial=False, temporal=False, snapshot_interval=cfg.interval_map_snapshots,
                 window_size=100, window_overlap=0, datatype='map', full_name='', unit='unit'):
        self.name = name
        self.spatial = spatial
        self.temporal = temporal
        self.snapshot_interval = snapshot_interval
        self.window_size = window_size
        self.window_overlap = window_overlap
        self.datatype = datatype
        self.full_name = full_name
        self.unit = unit


# State variables for EWS #
# Maximum interception store
micM = StateVariable('micM', spatial=True, full_name='Maximum interception storage spatial')
micA = StateVariable('micA', temporal=True, datatype='numpy', full_name='Maximum interception storage temporal')

# LAI
laiM = StateVariable('laiM', spatial=True, full_name='LAI spatial')
laiA = StateVariable('laiA', temporal=True, datatype='numpy', full_name='LAI temporal')

# Soil moisture
moiM = StateVariable('moiM', spatial=True, full_name='Soil moisture spatial')
moiA = StateVariable('moiA', temporal=True, datatype='numpy', full_name='Soil moisture temporal')

# Biomass
bioM = StateVariable('bioM', spatial=True, full_name='Biomass spatial')
bioA = StateVariable('bioA', temporal=True, datatype='numpy', full_name='Biomass temporal')

# Regolith thickness
regM = StateVariable('regM', spatial=True, full_name='Regolith thickness spatial')
regA = StateVariable('regA', temporal=True, datatype='numpy', full_name='Regolith thickness temporal')

# DEM
demM = StateVariable('demM', spatial=True, full_name='DEM spatial')
demA = StateVariable('demA', temporal=True, datatype='numpy', full_name='DEM temporal')

# Discharge
qA = StateVariable('qA', temporal=True, datatype='numpy', full_name='Discharge temporal')

# Grazing rate
gM = StateVariable('gM', spatial=True, full_name='Grazing rate spatial')
gA = StateVariable('gA', temporal=True, datatype='numpy', full_name='Grazing rate temporal')

# Growth part
gpM = StateVariable('gpM', spatial=True, full_name='Growth part spatial')
gpA = StateVariable('gpA', temporal=True, datatype='numpy', full_name='Growth part temporal')

# Grazing part
grM = StateVariable('grM', spatial=True, full_name='Grazing part spatial')
grA = StateVariable('grA', temporal=True, datatype='numpy', full_name='Grazing part temporal')

# Net growth
grnM = StateVariable('grnM', spatial=True, full_name='Net growth spatial')
grnA = StateVariable('grnA', temporal=True, datatype='numpy', full_name='Net growth temporal')

# Net deposition
depM = StateVariable('depM', spatial=True, full_name='Net deposition spatial')
depA = StateVariable('depA', temporal=True, datatype='numpy', full_name='Net deposition temporal')

# Net weathering
weaM = StateVariable('weaM', spatial=True, full_name='Net weathering spatial')
weaA = StateVariable('weaA', temporal=True, datatype='numpy', full_name='Net weathering temporal')

# Net creep deposition
creM = StateVariable('creM', spatial=True, full_name='Net creep deposition spatial')
creA = StateVariable('creA', temporal=True, datatype='numpy', full_name='Net creep deposition temporal')

# Grazing pressure array
grazing = StateVariable('grazing', temporal=True, full_name='Grazing pressure temporal')

# Ksat value
Iks = StateVariable('Iks', spatial=True, full_name='Ksat value spatial')

# Check which variables are present in the configuration and append these to the list of variables #

full_set_of_variables = [micM, micA, laiM, laiA, moiM, moiA, bioM, bioA, regM, regA, demM, demA, qA, gM, gA, gpM, gpA,
                         grM, grA, grnM, grnA, depM, depA, weaM, weaA, creM, creA, grazing, Iks]

if cfg.state_variables_for_ews == 'full':
    variables = full_set_of_variables
else:
    for state_variable in cfg.state_variables_for_ews:
        for variable in full_set_of_variables:
            if variable.name == state_variable:
                variables.append(variable)
