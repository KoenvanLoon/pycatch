import pandas as pd
import numpy as np

from statsmodels.tsa.arima_process import ArmaProcess
from statsmodels.tsa.arima_model import ARMA

from statsmodels.tsa.ar_model import AutoReg, ar_select_order
from statsmodels.tsa.arima.model import ARIMA

from statsmodels.stats.diagnostic import het_arch

import EWSPy as ews

ar1 = np.array([1, -0.9])
AR_object1 = ArmaProcess(ar1)
simulated_data_1 = AR_object1.generate_sample(nsample=1000)
simulated_data_2 = AR_object1.generate_sample(nsample=1000)

mod = AutoReg(simulated_data_1, 1).fit()
print(mod.summary())
print("When the true phi=0.9, the estimate of the constant and phi are:", mod.params)
print(np.reciprocal(mod.params))

residuals = mod.resid
AutoRegCondHet = het_arch(residuals, nlags=4, ddof=1)
print(AutoRegCondHet)
