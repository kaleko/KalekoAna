import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE RUN SUBRUN EVENTID\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem
from ROOT import larlite as fmwk

#gSystem.Load('libKalekoAna_scratch_ana.so')

# Create ana_processor instance
my_proc = fmwk.ana_processor()
my_proc.enable_filter(True)

# Set input root file
for x in xrange(len(sys.argv)-3):
    if x == 0: continue
    my_proc.add_input_file(sys.argv[x])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD) #kBOTH if you want to save that event

# Specify output root file name
my_proc.set_ana_output_file("");
#my_proc.set_output_file("FoundEvent.root") #uncomment if you want to save that event

# Attach a template process
mod = fmwk.FindEventByRunSubrunEvt()
mod.setRunNo(int(sys.argv[-3]))
mod.setSubRunNo(int(sys.argv[-2]))
mod.setEvtNo(int(sys.argv[-1]))
my_proc.add_process(mod)

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run();

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

