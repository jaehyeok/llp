import ROOT
import sys
import numpy as np
from numpy import sqrt, arctanh
from numpy.linalg import norm


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


R = 1200    #Radius of barrel in mm

gluino_id = 1000021
gluon_id = 21
gravitino_id = 1000039

if len(sys.argv) != 3:
    print "USAGE: %s <input file> <output file>" %sys.argv[0]
    sys.exit(1)

inFileName = sys.argv[1]
outFileName = sys.argv[2]
print "Reading from", inFileName, "and writing to", outFileName

inFile = ROOT.TFile.Open(inFileName, "READ")
tree = inFile.Get("tree")

hist = ROOT.TH1D("hist", "Hit time(in ns)", 100, 0, 20)


for entryNum in range(0, tree.GetEntries()):

    if entryNum % 100 == 0: print entryNum

    tree.GetEntry(entryNum)
    id = np.array(getattr(tree, "id"))
    status = np.array(getattr(tree, "status"))
    dtr1 = np.array(getattr(tree, "dtr1"))
    dtr2 = np.array(getattr(tree, "dtr2"))
    e = np.array(getattr(tree, "e"))
    m = np.array(getattr(tree, "m"))
    vtx = np.array([getattr(tree, "vtx_x"), getattr(tree, "vtx_y"), getattr(tree, "vtx_z")]).T
    dvtx = np.array([getattr(tree, "dvtx_x"), getattr(tree, "dvtx_y"), getattr(tree, "dvtx_z")]).T
    p = np.array([getattr(tree, "px"), getattr(tree, "py"), getattr(tree, "pz")]).T
    r = sqrt(dvtx[:,0]**2+dvtx[:,1]**2)

    for i in range(0, id.size):
        if r[i]<R and id[i] == gluino_id and abs(status[i])<30 and eta(p[i])<1.48:

            t = r[i]/(beta(e[i], m[i])*unit_proj(p[i]))         #Momentum change during intermediate steps are negligible
            
            trueDecay = False       #True when ~g -> g ~G

            while trueDecay == False:       #Tracks the decay using that there is only one susy particle in a single step
                for j in range(dtr1[i], dtr2[i]+1):
                    if abs(id[j]) > 1000000 and abs(id[j]) != gravitino_id:      #Susy particles except gravitinos
                        i = j
                        break
                    else: trueDecay = True

            for j in range(dtr1[i], dtr2[i]+1):
                    if id[j] == gluon_id:
                        i = j
                        t = t + getHit(vtx[i],p[i],R)/unit_proj(p[i])
                        hist.Fill(0.0033 * t)       #mm/c to ns




'''
hist.SetDirectory(0)    #save as root file
inFile.Close()
outHistFile = ROOT.TFile.Open(outFileName, "RECREATE")
outHistFile.cd()
hist.Write()
hist.Draw()
hist.Print()
outHistFile.Close()
'''

can = ROOT.TCanvas("can", "hist")       #Save as image
can.SetLogy()
hist.Draw()
can.SaveAs(outFileName)
