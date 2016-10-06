import os, sys

mysamples = [1, 2, 3]
#mylens = [10,20,30,40,50,60,70,80,90,100,150,200]
mylens = [125, 175, 225, 250, 275, 300]

outdir = '/Users/davidkaleko/Data/larlite/081716_mcs_len_scan/'
infiledir = '/Users/davidkaleko/Data/larlite/080116_selection_output/'

counter = 0
countermax = len(mysamples)*len(mylens)
for mylen in mylens:
	for mysample in mysamples:
		counter += 1
		infile = infiledir + 'saved_evts_trkpandoraNuPMA_vtxpmtrack_mcc71ext2bnb3_%d.root' % mysample
		command = 'python run_XiaoEventAna_MCSlenscan.py pandoraNuPMA pmtrack pandoraNuPMAcalo %d %0.1f %s %s'%(mysample, mylen, outdir, infile)
		print "Counter is at %d of %d..." % (counter, countermax)
		print "Command is: %s" % command
		sys.stdout.flush()
		os.system(command)

