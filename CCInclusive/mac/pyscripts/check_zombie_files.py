import os
from ROOT import *

basedir = '/Users/davidkaleko/Data/larlite/080116_selection_output/'

allfiles = os.listdir(basedir)
rootfiles = [ x for x in allfiles if '.root' in x ]

print "Checking for zombie files in base directory => %s" % basedir
counter = 0
for rootfile in rootfiles:
	counter += 1
	f = TFile(basedir+rootfile,'READ')
	if f.IsZombie():
		print "FOUND ZOMBIE FILE! => %s" % rootfile
	f.Close()

print "Done checking for zombie files (checked %d total files)." % counter