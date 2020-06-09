import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.pyplot import figure
from matplotlib import ticker, cm
from matplotlib import rc

if len(sys.argv) != 3:
    print "USAGE: %s <input file> <output file>" %sys.argv[0]
    sys.exit(1)

inFileName = sys.argv[1]
outFileName = sys.argv[2]
print "Reading from", inFileName, "and writing to", outFileName

data = np.load(inFileName)

tau = data[0]
mass = data[1]

xbin = np.logspace(2.45, 5.05, num=27)
ybin = np.linspace(950, 3050, 22)
w = np.repeat(0.001, len(mass))

font = {'family' : 'normal',
        'size'   : 15}

rc('font', **font)

counts, _, _ = np.histogram2d(tau, mass, bins=(xbin, ybin), weights = w)

fig, ax = plt.subplots(figsize = (12, 8))
cs = ax.pcolormesh(xbin, ybin, counts.T, norm = mpl.colors.LogNorm())
cbar = fig.colorbar(cs)

ax.set_xscale('log')

plt.title(outFileName)
plt.xlabel('Lifetime (ctau)')
plt.ylabel('Mass (GeV)')

plt.savefig("%s.png" %outFileName)
