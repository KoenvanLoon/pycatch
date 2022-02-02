import numpy as np
import math
from pcraster import *
import matplotlib.pyplot as plt

Fs = 8000
sample = 8000
x = np.arange(sample)
y = np.sin(2 * np.pi * 5 * x / Fs) + np.sin(2 * np.pi * 5 * x / Fs)

plt.plot(x, y)
plt.show()
