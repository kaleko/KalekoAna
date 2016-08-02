from ROOT import larlite as fmwk
import numpy as np
import matplotlib.pyplot as plt

a = fmwk.MuTrackMomentumSpline()

xes = np.linspace(1,1000,100)
ys = np.array([a.GetMuMomentum(poop) for poop in xes])
plt.plot(xes,ys,'ro--')
plt.title('Muon Energy based on Range from Spline',fontsize=16)
plt.xlabel('Muon Range [cm]',fontsize=16)
plt.ylabel('Muon Kinetic Energy [GeV]',fontsize=16)
#plt.xscale('log')
plt.grid(True,'both')
plt.show()

