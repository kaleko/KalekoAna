import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)
my_proc.enable_filter(True)

# Specify output root file name
my_proc.set_ana_output_file("TestMCTruth_ana_out.root");

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.

my_proc.add_process(fmwk.TestMCTruth())

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
#the event discussed in nevis mtg 061016 is (47,1)
#from ~/Data/larlite/mcc7_022316/scan_prodgenie_bnb_nu_uboone_mcc7_detsim_v1/scan_prodgenie_bnb_nu_uboone_mcc7_detsim_v1_all.root
my_proc.run(47,1) 
# my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
