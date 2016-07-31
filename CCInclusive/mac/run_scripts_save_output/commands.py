#Usage 1: run_Xiao_write_output_batch.py track_producer vtx_producer calo_producer mcc71_ext2_bnb3 $OUTPUT_DIR $INPUT_ROOT_FILE(s)

import sys,os

if len(sys.argv) != 2:
    print "USAGE: python commands.py NUMBER_OF_COMMAND_TO_RUN"

num = int(sys.argv[1])

outdir = '/uboone/data/users/kaleko/073116_selection_output/'

trackprod, vtxprod, caloprod = 'fuck', 'fuck', 'fuck'
mcc71_ext2_bnb3 = 1234

#mcc7 commands
if num == 1:
    trackprod, vtxprod, caloprod = 'pandoraNuPMA', 'pmtrack', 'pandoraNuPMAcalo'
    mcc71_ext2_bnb3 = 1
if num == 2:
    trackprod, vtxprod, caloprod = 'pandoraNu', 'pandoraNu', 'pandoraNucalo'
    mcc71_ext2_bnb3 = 1
if num == 3:
    trackprod, vtxprod, caloprod = 'KalekopandoraNuPMAPlustrackkalmanhit', 'pmtrack', 'pandoraNuPMAcalo'
    mcc71_ext2_bnb3 = 1
if num == 4:
    trackprod, vtxprod, caloprod = 'KalekopandoraNuPlustrackkalmanhit', 'pandoraNu', 'pandoraNucalo'
    mcc71_ext2_bnb3 = 1
#ext commands
if num == 5:
    trackprod, vtxprod, caloprod = 'pandoraNuPMA', 'pmtrack', 'pandoraNuPMAcalo'
    mcc71_ext2_bnb3 = 2
if num == 6:
    trackprod, vtxprod, caloprod = 'pandoraNu', 'pandoraNu', 'pandoraNucalo'
    mcc71_ext2_bnb3 = 2
if num == 7:
    trackprod, vtxprod, caloprod = 'KalekopandoraNuPMAPlustrackkalmanhit', 'pmtrack', 'pandoraNuPMAcalo'
    mcc71_ext2_bnb3 = 2
if num == 8:
    trackprod, vtxprod, caloprod = 'KalekopandoraNuPlustrackkalmanhit', 'pandoraNu', 'pandoraNucalo'
    mcc71_ext2_bnb3 = 2
#bnb commands
if num == 9:
    trackprod, vtxprod, caloprod = 'pandoraNuPMA', 'pmtrack', 'pandoraNuPMAcalo'
    mcc71_ext2_bnb3 = 3
if num == 10:
    trackprod, vtxprod, caloprod = 'pandoraNu', 'pandoraNu', 'pandoraNucalo'
    mcc71_ext2_bnb3 = 3
if num == 11:
    trackprod, vtxprod, caloprod = 'KalekopandoraNuPMAPlustrackkalmanhit', 'pmtrack', 'pandoraNuPMAcalo'
    mcc71_ext2_bnb3 = 3
if num == 12:
    trackprod, vtxprod, caloprod = 'KalekopandoraNuPlustrackkalmanhit', 'pandoraNu', 'pandoraNucalo'
    mcc71_ext2_bnb3 = 3

infiles = 'fuck'
if mcc71_ext2_bnb3 == 1:
    infiles = '/uboone/data/users/kaleko/kaleko_mcc7_bnbcosmic_samdef_fullsample_moretrackproducers_larlite_out/*mcinfo* /uboone/data/users/kaleko/kaleko_mcc7_bnbcosmic_samdef_fullsample_moretrackproducers_larlite_out/*reco* /uboone/data/users/kaleko/kaleko_mcc7_bnbcosmic_samdef_fullsample_moretrackproducers_larlite_out/stitched_tracks_KalekopandoraNuPMAPlus* /uboone/data/users/kaleko/kaleko_mcc7_bnbcosmic_samdef_fullsample_moretrackproducers_larlite_out/stitched_tracks_KalekopandoraNuPlus*'
if mcc71_ext2_bnb3 == 2:
    infiles = '/uboone/data/users/kaleko/kaleko_extbnb_reco_neutrion2016_samdef_fullsample_moretrackproducers_larlite_out/*reco* /uboone/data/users/kaleko/kaleko_extbnb_reco_neutrion2016_samdef_fullsample_moretrackproducers_larlite_out/stitched_tracks_KalekopandoraNuPMAPlus* /uboone/data/users/kaleko/kaleko_extbnb_reco_neutrion2016_samdef_fullsample_moretrackproducers_larlite_out/stitched_tracks_KalekopandoraNuPlus*'
if mcc71_ext2_bnb3 == 3:
    infiles = '/uboone/data/users/kaleko/kaleko_beamfilter_reco_neutrion2016_samdef_fullsample_moretrackproducers_larlite_out/*reco* /uboone/data/users/kaleko/kaleko_beamfilter_reco_neutrion2016_samdef_fullsample_moretrackproducers_larlite_out/fuckme_KalekopandoraNuPMAPlus* /uboone/data/users/kaleko/kaleko_beamfilter_reco_neutrion2016_samdef_fullsample_moretrackproducers_larlite_out/fuckme_KalekopandoraNuPlus*'


mycommand = "python run_Xiao_write_output_batch.py %s %s %s %d %s %s"%(trackprod,vtxprod,caloprod,mcc71_ext2_bnb3,outdir,infiles)
print "RUNNING COMMAND:"
print mycommand
sys.stdout.flush()
os.system(mycommand)
