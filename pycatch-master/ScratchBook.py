import numpy as np
import os

a = np.loadtxt('C:/Users/koenv/Desktop/Thesis/pycatch-master/pycatch-master/1/m2g0000/bioA0005.200.numpy.txt')
b = np.loadtxt('C:/Users/koenv/Desktop/Thesis/pycatch-master/pycatch-master/1/bioA0005.200.numpy.txt')

print(np.nanmean(a))
print(np.nanmean(b))
