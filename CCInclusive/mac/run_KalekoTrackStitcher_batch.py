import os, sys


#base_producer = "pandoraNuPMA"
#match_producer = "trackkalmanhit"

if len(sys.argv) != 5:
	print
	print "USAGE: python %s base_producer match_producer output_directory input_directory" % sys.argv[0]
	print
	quit()

#Make sure larlite is sourced correctly
if os.environ['LARLITE_BASEDIR'] != '/uboone/app/users/kaleko/larlite':
	print "Wrong larlite sourced, or no larlite sourced. Quitting."
	quit()

base_producer = sys.argv[1]
match_producer = sys.argv[2]
#outdir = '/uboone/data/users/kaleko/kaleko_mcc7_bnbcosmic_samdef_fullsample_moretrackproducers_larlite_out/'
outdir = sys.argv[3]
out_producer = "Kaleko%sPlus%s" % (base_producer,match_producer)

if not os.path.exists(outdir):
    os.makedirs(outdir)

#basedir = '/pnfs/uboone/scratch/users/kaleko/'
#basedir += 'kaleko_mcc7_bnbcosmic_samdef_fullsample_moretrackproducers_larlite_out/'
#basedir += 'v05_09_00/'
basedir = sys.argv[4]

infiles = os.listdir(basedir)
infiles = [ x for x in infiles if 'opreco' in x ]
n_infiles= len(infiles)
print "%d input files found." % n_infiles
#infiles = infiles[:2]

# for right now each opreco file is 600M, each mcinfo file is 200M
# goal is to have 10 output files for reco, 10 for mcinfo
n_final_files = 10

if n_infiles % n_final_files != 0:
	print "Number of desired final files is not an integer multiple of input files. Quitting."
	quit()

n_in_per_final = n_infiles / n_final_files


for ifile in xrange(n_final_files):
	infilelist = []
	these_infiles = infiles[ifile*n_in_per_final:(ifile*n_in_per_final+n_in_per_final)]
	for infile in these_infiles:
		infilelist.append(basedir+infile)
		
	# i ahve no fucking idea why but the name stitcked_tracks_%s blah blah segfaults but fuckme doesn't
	outfilename = 'fuckme_%s_file%d.root'%(out_producer,ifile)
	print
	print "Attempting to run over %d files to file %s..."%(len(infilelist),outfilename)
	print
	sys.stdout.flush()
	command = 'python '+os.environ['LARLITE_USERDEVDIR']+'/KalekoAna/CCInclusive/mac/run_KalekoTrackStitcher.py ' 
	command += base_producer + ' ' + match_producer + ' ' + ' '.join(infilelist) + ' ' + outdir
	command += outfilename
	print "about to do command: %s"%command
	os.system(command)

print "DONE!"
