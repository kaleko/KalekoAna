from ROOT import *

f = TFile("/Users/davidkaleko/Data/larlite/mcc7_022316/scan_prodgenie_bnb_nu_uboone_mcc7_detsim_v1/scan_prodgenie_bnb_nu_uboone_mcc7_detsim_v1_all.root","READ")
t = f.mctruth_generator_tree
blah = t.GetEntry(0)
b = t.mctruth_generator_branch
h = TH1F("h","Angle b/t Neutrino and z-direction [degrees]",1000,0,1)
zdir = TVector3(0,0,1)
for x in xrange(5000):
    blah = t.GetEntry(x)
    if b.size() != 1: continue
    angle = b.at(0).GetNeutrino().Nu().Trajectory().front().Momentum().Vect().Angle(zdir)
    angle *= (180./3.14159265)
    dummy = h.Fill(angle)

c1 = TCanvas()
h.Draw()
c1.Modified()
c1.Update()
gStyle.SetOptStat(1111111)
raw_input("")
