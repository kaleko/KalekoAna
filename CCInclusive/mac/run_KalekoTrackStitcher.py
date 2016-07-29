
import sys

if len(sys.argv) != 5:
    msg = '\n'
    msg += "Usage 1: %s base_producer match_producer input_root_file output_root_file\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


base_producer = sys.argv[1]
match_producer = sys.argv[2]
out_producer = "Kaleko%sPlus%s" % (base_producer,match_producer)


from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
my_proc.add_input_file(sys.argv[3])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)
my_proc.enable_filter(True)

# Specify output root file name
my_proc.set_ana_output_file("dummy_output")
my_proc.set_output_file(sys.argv[4])

my_proc.set_data_to_read(fmwk.data.kTrack,base_producer)
my_proc.set_data_to_read(fmwk.data.kTrack,match_producer)
my_proc.set_data_to_write(fmwk.data.kTrack,out_producer)


mystitch = fmwk.KalekoTrackStitcher()
mystitch.setBaseProducer(base_producer)
mystitch.setMatchProducer(match_producer)
mystitch.setOutputProducer(out_producer)
my_proc.add_process(mystitch)

print
print "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
