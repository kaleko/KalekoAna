import sys
import os

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem
from ROOT import larlite as fmwk
from ROOT import ertool

from seltool.primarycosmicDef import GetPrimaryCosmicFinderInstance

gSystem.Load('libKalekoAna_scratch_ana.so')

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)):
    if x == 0: continue
    my_proc.add_input_file(sys.argv[x])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("ana_out_MCShowerStudy_MCTrackStudy_CosmicValidation.root");

# Attach a template process
mcsstudy = fmwk.MCShowerStudy()
mcsstudy.setMCShowerProducer('mcreco')
my_proc.add_process(mcsstudy)

mctstudy = fmwk.MCTrackStudy()
mctstudy.setMCShowerProducer('mcreco')
my_proc.add_process(mctstudy)

cosmicprimary_algo = GetPrimaryCosmicFinderInstance()
cosmicsecondary_algo = ertool.ERAlgoCRSecondary()
cosmicorphanalgo = ertool.ERAlgoCROrphan()

Ecut = 50 # in MeV

cos_ana = ertool.ERAnaCosmicValidation()
cos_ana.SetECut(Ecut)

main_anaunit = fmwk.ExampleERSelection()
main_anaunit._mgr.ClearCfgFile()
main_anaunit._mgr.AddCfgFile(os.environ['LARLITE_USERDEVDIR']+'/SelectionTool/ERTool/dat/ertool_default.cfg')

main_anaunit.SetShowerProducer(True,'mcreco')
main_anaunit.SetTrackProducer(True,'mcreco')

# from seltool.trackDresserDef import GetTrackDresserInstance
# dresser = GetTrackDresserInstance()
# main_anaunit._mgr.AddAlgo(dresser)

main_anaunit._mgr.AddAlgo(cosmicprimary_algo)
main_anaunit._mgr.AddAlgo(cosmicsecondary_algo)
main_anaunit._mgr.AddAlgo(cosmicorphanalgo)

main_anaunit._mgr.AddAna(cos_ana)
main_anaunit._mgr._profile_mode = True

main_anaunit.SetMinEDep(Ecut)
main_anaunit._mgr._mc_for_ana = True

my_proc.add_process(main_anaunit)


# Let's run it.
my_proc.run();

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

