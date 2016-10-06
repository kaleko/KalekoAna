#track_producer = KalekopandoraNuPlustrackkalmanhit #KalekopandoraNuPMAPlustrackkalmanhit #pandoraNuPMA #pandoraNu
#vtx_producer = pmtrack #pandoraNu

import sys

if len(sys.argv) < 8:
    msg = '\n'
    msg += "Usage 1: %s track_producer vtx_producer calo_producer mcc71_ext2_bnb3 mcs_minlen $OUTPUT_DIR $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk

track_producer = sys.argv[1]
vtx_producer = sys.argv[2]
calo_producer = sys.argv[3]
mcc71_ext2_bnb3 = int(sys.argv[4])
mcs_minlen = float(sys.argv[5])
outdir = sys.argv[6]

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv) - 7):
    my_proc.add_input_file(sys.argv[x + 7])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)
my_proc.enable_filter(True)

# Specify output root file name
anaoutname = 'ana_out_trk%s_vtx%s_mcc71ext2bnb3_%d_MCSminlen_%0.2f.root'%(track_producer,vtx_producer,mcc71_ext2_bnb3,mcs_minlen)
my_proc.set_ana_output_file(outdir+'/'+anaoutname)

myxiao = fmwk.XiaoEventAna()
myxiao.setMCSMinLen(mcs_minlen)

run_on_data = True
if mcc71_ext2_bnb3 == 1: run_on_data = False
myxiao.setRunningOnData(run_on_data)
myxiao.setVtxSphereRadius(4.0)
mytype = 4#fmwk.kBNBCosmic
if mcc71_ext2_bnb3 == 2: mytype = 1#fmwk.kOffBeam
if mcc71_ext2_bnb3 == 3: mytype = 0#fmwk.kOnBeam

myxiao.setInputType(mytype)
myxiao.setTrackProducer(track_producer)
myxiao.setVtxProducer(vtx_producer)
myxiao.setCaloProducer(calo_producer)
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
