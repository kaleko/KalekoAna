import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE OUTPUT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem
from ROOT import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

#multiple input files
for x in xrange(len(sys.argv)-1):
    if x == 0:
        continue
    my_proc.add_input_file(sys.argv[x])

# Specify output root file name
my_proc.set_ana_output_file(sys.argv[-1])

# Attach a template process
myunit = fmwk.MCSBiasStudyDriver()
#myunit.SetAnalysisType(fmwk.MCSBiasStudy.kSingleMuonMCTrack)
#myunit.SetAnalysisType(fmwk.MCSBiasStudy.kSingleMuonRecoTrack)
#myunit.SetAnalysisType(fmwk.MCSBiasStudy.kMCBNBSelectedRecoTrack)
myunit.SetAnalysisType(fmwk.MCSBiasStudy.kMCBNBRecoTrack)
#myunit.SetAnalysisType(fmwk.MCSBiasStudy.kDataBNBSelectedRecoTrack)

my_proc.add_process(myunit)



print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

