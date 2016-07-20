import sys

if len(sys.argv) < 2:
    msg = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk


# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv) - 1):
    my_proc.add_input_file(sys.argv[x + 1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)
my_proc.set_output_file("aho.root")
#my_proc.set_io_mode(fmwk.storage_manager.kREAD)
my_proc.enable_filter(True)

# Specify output root file name
my_proc.set_ana_output_file("temp.root")

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
# my_proc.add_process(fmwk.MC_1mu1pNn0else_Filter_MCTracks())#fmwk.MC_1mu1pNn0else_Filter())
#myfilter = fmwk.NuMuCCFilter()
#myfilter.FlipFilter(False)
# myfilter.SetNuMuFromKaonOnly(False)
# myfilter.SetMinNuEnergy(2.4);
#my_proc.add_process(myfilter)
myxiao = fmwk.XiaoEventAna()
myxiao.setRunningOnData(False)
myxiao.setVtxSphereRadius(4.0)
myxiao.setInputType(fmwk.kBNBCosmic)#fmwk.kBNBCosmic fmwk.kCorsikaInTime
#print fmwk.kBNBCosmic
my_proc.add_process(myxiao)

my_proc.set_data_to_read(fmwk.data.kMCTruth,"generator")
my_proc.set_data_to_read(fmwk.data.kMCShower,"mcreco")
my_proc.set_data_to_read(fmwk.data.kOpFlash,"opflashSat")
my_proc.set_data_to_read(fmwk.data.kTrack,"pandoraNuPMA")
my_proc.set_data_to_read(fmwk.data.kVertex,"pmtrack")
my_proc.set_data_to_read(fmwk.data.kCalorimetry,"pandoraNuPMAcalo")
my_proc.set_data_to_read(fmwk.data.kMCTrack,"mcreco")
my_proc.set_data_to_read(fmwk.data.kAssociation,"pandoraNuPMAcalo")

my_proc.set_data_to_write(fmwk.data.kMCTruth,"generator")
my_proc.set_data_to_write(fmwk.data.kMCShower,"mcreco")
my_proc.set_data_to_write(fmwk.data.kOpFlash,"opflashSat")
my_proc.set_data_to_write(fmwk.data.kTrack,"pandoraNuPMA")
my_proc.set_data_to_write(fmwk.data.kVertex,"pmtrack")
my_proc.set_data_to_write(fmwk.data.kCalorimetry,"pandoraNuPMAcalo")
my_proc.set_data_to_write(fmwk.data.kMCTrack,"mcreco")
my_proc.set_data_to_write(fmwk.data.kAssociation,"pandoraNuPMAcalo")
"""
    [NORMAL]  <close> TTree "gtruth_generator_tree" for gtruth written with 15 events...
    [NORMAL]  <close> TTree "mctruth_corsika_tree" for mctruth written with 15 events...
    [NORMAL]  <close> TTree "mctruth_generator_tree" for mctruth written with 15 events...
    [NORMAL]  <close> TTree "mcflux_generator_tree" for mcflux written with 15 events...
    [NORMAL]  <close> TTree "mcshower_mcreco_tree" for mcshower written with 15 events...
    [NORMAL]  <close> TTree "opflash_opflashSat_tree" for opflash written with 15 events...
    [NORMAL]  <close> TTree "track_pandoraNu_tree" for track written with 15 events...
    [NORMAL]  <close> TTree "track_pandoraNuPMA_tree" for track written with 15 events...
    [NORMAL]  <close> TTree "vertex_pandoraNu_tree" for vertex written with 15 events...
    [NORMAL]  <close> TTree "vertex_pmtrack_tree" for vertex written with 15 events...
    [NORMAL]  <close> TTree "calo_pandoraNuPMAcalo_tree" for calo written with 15 events...
    [NORMAL]  <close> TTree "pfpart_pandoraNu_tree" for pfpart written with 15 events...
    [NORMAL]  <close> TTree "mctrack_mcreco_tree" for mctrack written with 15 events...
    [NORMAL]  <close> TTree "ass_opflashSat_tree" for ass written with 15 events...
    [NORMAL]  <close> TTree "ass_pandoraNu_tree" for ass written with 15 events...
    [NORMAL]  <close> TTree "ass_pandoraNuPMA_tree" for ass written with 15 events...
    [NORMAL]  <close> TTree "ass_pandoraNuPMAcalo_tree" for ass written with 15 events...
    [NORMAL]  <close> TTree "ass_pmtrack_tree" for ass written with 15 events...
    [NORMAL]  <close> TTree "ass_showerrecopandora_tree" for ass written with 15 events...
"""

print
print "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run(0,2000)

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
