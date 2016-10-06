import sys

if len(sys.argv) < 5:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s) TRACK_PRODUCER SPS1Track0 SEGLEN\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-4):
    my_proc.add_input_file(sys.argv[x+1])

track_producer = sys.argv[-3]
sps1track0 = int(sys.argv[-2])
seglen = int(sys.argv[-1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("spstestout/TestTrkSPSAssn_ana_out_track%s_SpsOrTrack_%d_seglen%d.root" % (track_producer, sps1track0,seglen));

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
mymod = fmwk.TestTrkSPSAssn()
mymod.SetTrackProducer(track_producer)
mymod.SetSPS1Track0(sps1track0)
mymod.SetSeglen(seglen)

my_proc.add_process(mymod)

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
