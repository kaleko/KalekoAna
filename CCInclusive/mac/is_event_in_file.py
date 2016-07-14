import sys
from ROOT import *

event = int(sys.argv[1])
infile = sys.argv[2]

f = TFile(infile,"READ")
t = f.larlite_id_tree
entries = t.GetEntries()
found_event = False
print "# entries in this file is ",entries

for x in xrange(entries):
    blah = t.GetEntry(x)
    if t._event_id == event:
        print "EVENT %d WAS INDEED IN THIS FILE! TTREE INDEX %d" % (event,x)
        found_event = True
        break

if not found_event:
    print "Event %d NOT FOUND in this file." % event

quit()
