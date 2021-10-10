
def mean(data):
    number = len(data)
    mean = sum(data) / number
    return mean

def variance(data):
    number = len(data)
    mean = sum(data) / number
    variance = sum([(x - mean) ** 2 for x in data]) / number
    return variance

def varianceSample(data):
    number = len(data)
    mean = sum(data) / number
    variance = sum([(x - mean) ** 2 for x in data]) / (number - 1)
    error = (2 * (variance ** 2) / (number - 1)) ** (1/2)
    return variance, error

def skewness(data):
    number = len(data)
    mean = sum(data) / number
    numerator = sum([(x - mean) ** 3 for x in data]) / number
    denominator = sum([(x - mean) ** 2 for x in data]) / number
    skewness = numerator / (denominator ** (3/2))
    return skewness

def skewnessSample(data):
    number = len(data)
    mean = sum(data) / number
    numerator = sum([(x - mean) ** 3 for x in data]) / number
    denominator = sum([(x - mean) ** 2 for x in data]) / (number - 1)
    skewness = numerator / (denominator ** (3/2))
    error = (6 * (number - 2) / ((number + 1) * (number + 3))) ** (1/2)
    return skewness, error

def kurtosis(data):
    number = len(data)
    mean = sum(data) / number
    numerator = sum([(x - mean) ** 4 for x in data]) / number
    denominator = sum([(x - mean) ** 2 for x in data]) / number
    kurtosis = numerator / (denominator ** 2) - 3
    return kurtosis

def kurtosisSample(data):
    number = len(data)
    mean = sum(data) / number
    numerator = sum([(x - mean) ** 4 for x in data]) / number
    denominator = sum([(x - mean) ** 2 for x in data]) / number
    kurtosis = numerator / (denominator ** 2) - 3
    error = (
            (24 * (number - 2) * (number - 3))
            / (((number + 1) ** 2) * (number + 3) * (number + 5))
    ) ** (1/2)
    return kurtosis, error

def autocovariance(data, lag=1):
    number = len(data)
    mean = sum(data) / number
    autocovariance = sum([(data[i] - mean) * (data[i+lag] - mean) for i in range(number - lag)]) / number
    return autocovariance

def autocorrelation(data, lag=1):
    return autocovariance(data, lag=lag) / variance(data)

#########################################


import rpy2
import rpy2.robjects.packages as rpackages
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
import rpy2.robjects as robjects
from rpy2.robjects.vectors import StrVector

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)   # select the first mirror in the list for R packages
packageNamesR = ('gstat', 'automap', 'e1071', 'tseries','earlywarnings')

# check rpy2 version and R packages
if rpy2.__version__ != '2.9.4':
    print("Tested for rpy2 version 2.9.4, current version is", rpy2.__version__)
    print("Please make sure you use the correct version.")
names_to_install = [x for x in packageNamesR if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    print(f"Installing the following R packages: {names_to_install}")
    utils.install_packages(StrVector(names_to_install))

def generic_ews_timeseries(timeseries):
    robjects.r('''
  generic_ews <- function(timeseries) {
        library(earlywarnings)
        return = ch_ews(timeseries)
        return
        }
        ''')
    generic_ews_timeseries = robjects.r['generic_ews']
    timeseries_R = numpy_to_R(timeseries)
    result = generic_ews_timeseries(timeseries_R)
    return result

def numpy_to_R(numpy_object):
    R_object = 0
    if len(numpy_object.shape) not in (1, 2):
        raise ValueError("Zero, or more than 3 dimensions to numpy object not supported")
        print("Empty numpy object")
    elif len(numpy_object.shape) == 1:
        R_object = robjects.r.matrix(numpy_object, nrow=numpy_object.shape[0])
        # R_object = robjects.FloatVector(numpy_object)
    elif len(numpy_object.shape) == 2:
        R_object = robjects.r.matrix(numpy_object, nrow=numpy_object.shape[0],
                                     ncol=numpy_object.shape[1])
    return robjects.r.assign("Numpy object", R_object)

#########################################

import numpy as np
np.random.seed(101)
data = np.random.randint(0, 100, 1000)
