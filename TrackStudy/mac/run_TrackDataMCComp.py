import sys

if len(sys.argv) < 3:
    msg  = '\n'
    msg += "Usage 1: %s mc_1_data_0 $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem
from ROOT import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()
my_proc.enable_filter(True)
# Set input root files
for x in xrange(len(sys.argv)):
    if x < 2: continue
    my_proc.add_input_file(sys.argv[x])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

mc_1_data_0 = int(sys.argv[1])

# Specify output root file name
my_proc.set_ana_output_file("TrackDataMCComp_MC=%d_ana_out.root"%mc_1_data_0)

# Attach plot maker
mymod = fmwk.TrackDataMCComp()
mymod.setTrackProducer("pandoraNuPMA")
mymod.setRunningOnData(False)
my_proc.add_process(mymod)

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run(0,100)

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

