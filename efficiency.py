import ROOT
import sys
import numpy as np
from numpy import sqrt, arctanh
from numpy.linalg import norm
import matplotlib.pyplot as plt


def beta(e,m):
    return sqrt(1-(m/e)**2)

def unit_vector(vector):
    return vector / norm(vector)

def unit_proj(vector):
    vector = unit_vector(vector)
    return sqrt(vector[0]**2+vector[1]**2)

def eta(p):     #Returns pseudorapidity
    return arctanh(1-unit_proj(p))

def getDist(vtx, dvtx):
    return norm(dvtx-vtx)

def getHit(vtx, p, r):      #Returns the distance btw vtx and hit point at barrel 

    m = p[1]/p[0]
    k = vtx[1]-m*vtx[0]     

    if p[0]<0:
        x = (-2*m*k-sqrt((2*m*k)**2-4*(m**2+1)*(k**2-r**2)))/(2*(m**2+1))
        y = m*x+k
        return sqrt((vtx[0]-x)**2+(vtx[1]-y)**2)

    else:
        x = (-2*m*k+sqrt((2*m*k)**2-4*(m**2+1)*(k**2-r**2)))/(2*(m**2+1))
        y = m*x+k
        return sqrt((vtx[0]-x)**2+(vtx[1]-y)**2)


if len(sys.argv) != 4:
    print "USAGE: %s <barrel radius> <input file> <output file>" %sys.argv[0]
    sys.exit(1)

R = float(sys.argv[1])
inFileName = sys.argv[2]
outFileName = sys.argv[3]
print "Reading from", inFileName, "and writing to", outFileName

inFile = ROOT.TFile.Open(inFileName, "READ")
tree = inFile.Get("tree")

mass = []
tau = []


for entryNum in range(0, tree.GetEntries()):

    detect = False

    if entryNum % 1000 == 0: print entryNum

    tree.GetEntry(entryNum)
    id = np.array(getattr(tree, "id"))
    tau0 = np.array(getattr(tree, "tau0"))
    e = np.array(getattr(tree, "e"))
    m = np.array(getattr(tree, "m"))
    vtx = np.array([getattr(tree, "vtx_x"), getattr(tree, "vtx_y"), getattr(tree, "vtx_z")]).T
    dvtx = np.array([getattr(tree, "dvtx_x"), getattr(tree, "dvtx_y"), getattr(tree, "dvtx_z")]).T
    p = np.array([getattr(tree, "px"), getattr(tree, "py"), getattr(tree, "pz")]).T
    r = sqrt(dvtx[:,0]**2+dvtx[:,1]**2)

    for i in range(0, id.size):
        if r[i]<R and eta(p[i])<1.48:

            g_mass = m[i]
            g_tau0 = tau0[i]

            t = r[i]/(beta(e[i], m[i])*unit_proj(p[i]))         #Momentum change during intermediate steps are negligible
            t_0 = r[i]/unit_proj(p[i])

            if t-t_0>910: detect = True

    if(detect):
        mass.append(g_mass) 
        tau.append(g_tau0)

data = np.array([tau, mass])
np.save(outFileName, data)

'''
xbin = np.logspace(2.5, 5, num=6)
ybin = np.linspace(1000, 3500, 6)
w = np.repeat(1/200, len(mass))

plt.hist2d(tau, mass, bins = [xbin, ybin])
plt.xscale('log')
plt.savefig(outFileName)
'''

