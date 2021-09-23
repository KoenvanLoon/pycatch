import matplotlib.pyplot as plt
import numpy

# script to create plots for timeseries output of main_weekly.py

# first index time steps, second index samples, third row number, fourth column number 
gA = numpy.load('gA.npy')        # GRAZING RATE #
bioA = numpy.load('bioA.npy')    # BIOMASS #
regA = numpy.load('regA.npy')    # REGOLITH THICKNESS #
sfA = numpy.load('sfA.npy')      # moisture content
qA = numpy.load('qA.npy')        # DISCHARGE #
gpA = numpy.load('gpA.npy')      # growth part
grA = numpy.load('grA.npy')      # grazing
grNA = numpy.load('grNA.npy')    # net growth
depA = numpy.load('depA.npy')    # net deposition
weaA = numpy.load('weaA.npy')    # net weathering
creA = numpy.load('creA.npy')    # net deposition due to creep

demA = numpy.load('demA.npy')    # DEM #

fig=plt.figure()
sample = 0
row = 0

location = 0
#bioAMA=bioA[:,sample,row,location]
test01 = demA[:,sample,row,location]
one = fig.add_subplot(231)
#one.plot(bioAMA)
one.plot(test01)

location = 7
#bioAMA=bioA[:,sample,row,location]
test01 = demA[:,sample,row,location]
one = fig.add_subplot(232)
#one.plot(bioAMA)
one.plot(test01)

plt.show()

fig.savefig("plotTimeseries_weekly.pdf",format="pdf")
