import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE OUTPUT_FILE ANALYSIS_TYPE\n" % sys.argv[0]
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
for x in xrange(len(sys.argv)-2):
    if x == 0:
        continue
    my_proc.add_input_file(sys.argv[x])

# Specify output root file name
my_proc.set_ana_output_file(sys.argv[-2])

ana_type_str = sys.argv[-1]
if ana_type_str == 'SingleMuonMCTrack':
	ana_type = fmwk.MCSBiasStudy.kSingleMuonMCTrack
elif ana_type_str == 'SingleMuonRecoTrack':
	ana_type = fmwk.MCSBiasStudy.kSingleMuonRecoTrack
elif ana_type_str == 'MCBNBSelectedRecoTrack':
	ana_type = fmwk.MCSBiasStudy.kMCBNBSelectedRecoTrack
elif ana_type_str == 'MCBNBRecoTrack':
	ana_type = fmwk.MCSBiasStudy.kMCBNBRecoTrack
elif ana_type_str == 'MCBNBMCTrack':
    ana_type = fmwk.MCSBiasStudy.kMCBNBMCTrack
elif ana_type_str == 'DataBNBSelectedRecoTrack':
	ana_type = fmwk.MCSBiasStudy.kDataBNBSelectedRecoTrack
elif ana_type_str == 'MCBNBMCTrackExiting':
    ana_type = fmwk.MCSBiasStudy.kMCBNBMCTrackExiting
else:
	print "ERROR ANA TYPE NOT CORRECTLY SPECIFIED"
	quit()

# Attach a template process
myunit = fmwk.MCSBiasStudyDriver()
myunit.SetAnalysisType(ana_type)
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

