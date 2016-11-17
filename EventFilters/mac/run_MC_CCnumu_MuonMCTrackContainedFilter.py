import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem
from ROOT import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
my_proc.add_input_file(sys.argv[1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)
my_proc.enable_filter(True)

# Specify output root file name
my_proc.set_ana_output_file("");
my_proc.set_output_file('MC_CCnumu_MuonMCTrackContainedFilter_savedevents.root')

# Attach a template process
my_proc.add_process(fmwk.MC_CCnumu_MuonMCTrackContainedFilter());

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run(0,1000)

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

