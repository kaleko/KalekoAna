import os, sys

data_basedir = '/Users/davidkaleko/Data/larlite/'

anatype_dict = {
#'SingleMuonMCTrack' : data_basedir + 'kaleko_prod_muminus_100616_larlite_out/larlite_*.root',
#'SingleMuonRecoTrack' : data_basedir + 'kaleko_prod_muminus_100616_larlite_out/larlite_*.root',

#'MCBNBSelectedRecoTrack' : data_basedir + '082416_mcc7_savedlongesttracks.root',
#'MCBNBRecoTrack' : data_basedir + 'kaleko_mcc7_bnb_minimalproducts/ccnumu_longmuontrack_fullycontained_savedevents.root',
#'MCBNBMCTrack' : data_basedir + 'kaleko_mcc7_bnb_minimalproducts/ccnumu_longmuontrack_fullycontained_savedevents.root',
#'DataBNBSelectedRecoTrack' : data_basedir + 'data_selected_contained_1m/savedlongesttracks.root'
'MCBNBMCTrackExiting' : data_basedir + 'kaleko_mcc7_bnb_minimalproducts/ccnumu_longmuontrack_exiting_savedevents.root'
}

outdir = '/Users/davidkaleko/Desktop/MCS_Technote/anafiles/'

runscript = '/Users/davidkaleko/larlite/UserDev/KalekoAna/scratch_ana/mac/run_MCSBiasStudyDriver.py'

for myanatype, myinfiles in anatype_dict.iteritems():
	print "RUNNING TYPE:",myanatype
	outfile = outdir + 'MCSBiasStudy_%s_anaout_10cmseg_2res_bothscatters_nonrelfix.root' % myanatype
	os.system('python %s %s %s %s' % (runscript, myinfiles, outfile, myanatype))

print "Done!"