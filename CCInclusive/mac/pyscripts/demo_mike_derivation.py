from ROOT import *
from ROOT import larlite as fmwk
calc = fmwk.NuEnergyCalc()
En = 938.+25
Ep = 1263.
thetamu = 0.79453
thetap = 0.987699
pp = 845.769
print "Beginning with the following: (En = %0.2f, Ep = %0.2f, thetamu = %0.2f, thetap = %0.2f, pp = %0.2f)" \
% ( En, Ep, thetamu, thetap, pp )

Emu = calc.ComputeEmu1mu1pQuadratic(En,Ep,thetamu,thetap,pp) * 1000.
delta = Emu - max(106.,pow(Emu*Emu - 106.*106.,0.5))
print "Iteration 0: Emu = %0.2f, delta = %0.5f" % (Emu, delta)

Emu = calc.ComputeEmu1mu1pQuadratic(En,Ep,thetamu,thetap,pp,delta) * 1000.
delta = Emu - max(106.,pow(Emu*Emu - 106.*106.,0.5))
print "Iteration 1: Emu = %0.2f, delta = %0.5f" % (Emu, delta)

Emu = calc.ComputeEmu1mu1pQuadratic(En,Ep,thetamu,thetap,pp,delta) * 1000.
delta = Emu - max(106.,pow(Emu*Emu - 106.*106.,0.5))
print "Iteration 2: Emu = %0.2f, delta = %0.5f" % (Emu, delta)

Emu = calc.ComputeEmu1mu1pQuadratic(En,Ep,thetamu,thetap,pp,delta) * 1000.
delta = Emu - max(106.,pow(Emu*Emu - 106.*106.,0.5))
print "Iteration 3: Emu = %0.2f, delta = %0.5f" % (Emu, delta)

Emu = calc.ComputeEmu1mu1pQuadratic(En,Ep,thetamu,thetap,pp,delta) * 1000.
delta = Emu - max(106.,pow(Emu*Emu - 106.*106.,0.5))
print "Iteration 4: Emu = %0.2f, delta = %0.5f" % (Emu, delta)

Emu = calc.ComputeEmu1mu1pQuadratic(En,Ep,thetamu,thetap,pp,delta) * 1000.
delta = Emu - max(106.,pow(Emu*Emu - 106.*106.,0.5))
print "Iteration 5: Emu = %0.2f, delta = %0.5f" % (Emu, delta)

Emu = calc.ComputeEmu1mu1pQuadratic(En,Ep,thetamu,thetap,pp,delta) * 1000.
delta = Emu - max(106.,pow(Emu*Emu - 106.*106.,0.5))
print "Iteration 6: Emu = %0.2f, delta = %0.5f" % (Emu, delta)

print "Now I will use the self-contained C++ compiled iterative technique:"
print "Beginning with the following: (En = %0.2f, Ep = %0.2f, thetamu = %0.2f, thetap = %0.2f, pp = %0.2f)" \
% ( En, Ep, thetamu, thetap, pp )
result = calc.ComputeEmu1mu1pQuadraticIterative(En,Ep,thetamu,thetap,pp)
print "Iterative technique has found this result for Emu:",result