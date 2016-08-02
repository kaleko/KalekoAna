import sys

if len(sys.argv) != 2:
    print "Usage: python %s input_file.root" % sys.argv[0]

infile = sys.argv[1]

issignalsample = 'signal' in infile

import pandas as pd
from root_numpy import root2array

df = pd.DataFrame( root2array( infile, 'evttree' ) )

print "Total number of events analyzed: ",len(df)

cuts = []

if issignalsample:
    cuts.append('is_truth_fiducial')
    cuts.append('is_numuCC')
    

cuts.append('flash_in_bgw')
cuts.append('is_areco_vtx_in_fidvol')
cuts.append('is_atrack_fromvtx')
cuts.append('longest_track_nearflash_z')
cuts.append('longest_trk_contained')
cuts.append('longest_trk_range_longenough')


if issignalsample:
	cuts.append('longest_trk_start_near_truth_nu')

nocuts = "is_truth_fiducial and is_numuCC"

mycut = ''
for x in xrange(len(cuts)):
    if not x:
        mycut += cuts[x]
    else:
        mycut += ' and %s' % cuts[x]
    print 'After applying cut on %s, %d events remain.' % (cuts[x],len(df.query(mycut)))


num = float(len(df.query(mycut)))

denom = len(df) if not issignalsample else len(df.query(nocuts))

print "Total Efficiency: %d/%d = %0.2f%%" % (num,denom,100*(num/denom))
