import sys

# python run_XiaoEventAna_batch.py vtx_sphere_radius flip_filter input_files

if len(sys.argv) < 3:
    msg = '\n'
    msg += "Usage 1: %s vtx_sphere_radius flip_filter $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk

vtx_sphere_radius = float(sys.argv[1])
flip_filter = bool(int(sys.argv[2]))

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv) - 3):
    my_proc.add_input_file(sys.argv[x + 3])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)
my_proc.enable_filter(True)

# Specify output root file name
my_proc.set_ana_output_file("batch_out/XiaoEventAna_out_filterflip%d_radius%0.2f.root"%(int(flip_filter),vtx_sphere_radius))

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
# my_proc.add_process(fmwk.MC_1mu1pNn0else_Filter_MCTracks())#fmwk.MC_1mu1pNn0else_Filter())
myfilter = fmwk.NuMuCCFilter()
myfilter.FlipFilter(flip_filter)
# myfilter.SetNuMuFromKaonOnly(False)
# myfilter.SetMinNuEnergy(2.4);
my_proc.add_process(myfilter)
myxiao = fmwk.XiaoEventAna()
myxiao.setRunningOnData(False)
myxiao.setInputType(fmwk.kBNBCosmic)#fmwk.kBNBCosmic fmwk.kCorsikaInTime
myxiao.setVtxSphereRadius(vtx_sphere_radius)
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
