import sys

if len(sys.argv) < 4:
    msg = '\n'
    msg += "Usage 1: %s data_1_mc_0 mcbnbcosmic_0_bnbdata_1_bnbextdata_2 $INPUT_ROOT_FILE(s)\n" % sys.argv[
        0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk


# Create ana_processor instance
my_proc = fmwk.ana_processor()
data_1_mc_0 = int(sys.argv[1])
mcbnbcosmic_0_bnbdata_1_bnbextdata_2 = int(sys.argv[2])
print "DATA_1_MC_0 IS %d, mcbnbcosmic_0_bnbdata_1_bnbextdata_2 IS %d." % (data_1_mc_0, mcbnbcosmic_0_bnbdata_1_bnbextdata_2)
# Set input root file
for x in xrange(len(sys.argv) - 3):
    my_proc.add_input_file(sys.argv[x + 3])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)
my_proc.set_output_file("XiaoEventAna_andStitching_out_data1mc0_%d_type_%d.root"
                        % (data_1_mc_0, mcbnbcosmic_0_bnbdata_1_bnbextdata_2))
# my_proc.set_io_mode(fmwk.storage_manager.kREAD)
my_proc.enable_filter(True)

# Specify output root file name
my_proc.set_ana_output_file("XiaoEventAna_andStitching_ana_out_data1mc0_%d_type_%d.root"
                            % (data_1_mc_0, mcbnbcosmic_0_bnbdata_1_bnbextdata_2))

base_producer = "pandoraNuPMA"
match_producer = "trackkalmanhit"
out_producer = "Kaleko%sPlus%s" % (base_producer, match_producer)


mystitch = fmwk.KalekoTrackStitcher()
mystitch.setBaseProducer(base_producer)
mystitch.setMatchProducer(match_producer)
mystitch.setOutputProducer(out_producer)
my_proc.add_process(mystitch)


if not data_1_mc_0:
    my_proc.set_data_to_read(fmwk.data.kMCFlux, "generator")
    my_proc.set_data_to_read(fmwk.data.kMCTruth, "generator")
    my_proc.set_data_to_read(fmwk.data.kMCShower, "mcreco")
    my_proc.set_data_to_read(fmwk.data.kMCTrack, "mcreco")

# my_proc.set_data_to_read(fmwk.data.kOpFlash, "opflashSat")
my_proc.set_data_to_read(fmwk.data.kTrack, base_producer)
my_proc.set_data_to_read(fmwk.data.kTrack, match_producer)
# my_proc.set_data_to_read(fmwk.data.kVertex, "pmtrack")
# my_proc.set_data_to_read(fmwk.data.kCalorimetry, "pandoraNuPMAcalo")
# my_proc.set_data_to_read(fmwk.data.kAssociation, "pandoraNuPMAcalo")

if not data_1_mc_0:
    my_proc.set_data_to_write(fmwk.data.kMCFlux, "generator")
    my_proc.set_data_to_write(fmwk.data.kMCTruth, "generator")
    my_proc.set_data_to_write(fmwk.data.kMCShower, "mcreco")
    my_proc.set_data_to_write(fmwk.data.kMCTrack, "mcreco")

# my_proc.set_data_to_write(fmwk.data.kOpFlash, "opflashSat")
# my_proc.set_data_to_write(fmwk.data.kTrack, base_producer)
# my_proc.set_data_to_write(fmwk.data.kTrack, match_producer)
my_proc.set_data_to_write(fmwk.data.kTrack, out_producer)
# my_proc.set_data_to_write(fmwk.data.kVertex, "pmtrack")
# my_proc.set_data_to_write(fmwk.data.kCalorimetry, "pandoraNuPMAcalo")
# my_proc.set_data_to_write(fmwk.data.kAssociation, "pandoraNuPMAcalo")


myxiao = fmwk.XiaoEventAna()
myxiao.setRunningOnData(data_1_mc_0)
myxiao.setVtxSphereRadius(4.0)
mytype = fmwk.kBNBCosmic
if mcbnbcosmic_0_bnbdata_1_bnbextdata_2 == 1:
	mytype = fmwk.kOnBeam
if mcbnbcosmic_0_bnbdata_1_bnbextdata_2 == 2:
	mytype = fmwk.kOffBeam
myxiao.setInputType(mytype)  # fmwk.kBNBCosmic fmwk.kCorsikaInTime
# print fmwk.kBNBCosmic
my_proc.add_process(myxiao)

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
