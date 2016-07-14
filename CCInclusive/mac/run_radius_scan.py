import os

radii = [ 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10., 20. ]

infiletext = '/Users/davidkaleko/Data/larlite/kaleko_mcc7_bnbcosmic_samdef_fullsample_larlite_out/'
infiletext += 'kaleko_mcc7_bnbcosmic_samdef_fullsample_mcinfo_file* '
infiletext += '/Users/davidkaleko/Data/larlite/kaleko_mcc7_bnbcosmic_samdef_fullsample_larlite_out/'
infiletext += 'kaleko_mcc7_bnbcosmic_samdef_fullsample_reco_file*'

for radius in radii:
	print " \n\n\n\n\n RUNNING FLIP TRUE WITH RADIUS %f \n\n\n\n\n" % radius
	runcommand = 'python run_XiaoEventAna_batch.py %f 1 ' % radius
	runcommand += infiletext
	os.system(runcommand)
	print " \n\n\n\n\n RUNNING FLIP FALSE WITH RADIUS %f \n\n\n\n\n" % radius
	runcommand = 'python run_XiaoEventAna_batch.py %f 0 ' % radius
	runcommand += infiletext
	os.system(runcommand)