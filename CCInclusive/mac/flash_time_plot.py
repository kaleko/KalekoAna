from ROOT import *
import sys

if len(sys.argv) != 4:
	print "USAGE: python %s input_file.root mc_true_data_false on_beam_true_off_beam_false" % sys.argv[0]
	quit()


f = TFile(sys.argv[1],'READ')
mc_true_data_false = int(sys.argv[2])
on_beam_true_off_beam_false = int(sys.argv[3])
t = f.opflash_opflashSat_tree
blah = t.GetEntry(0)
b = t.opflash_opflashSat_branch
entries = t.GetEntriesFast()

h = TH1F("h_flashtimes","Flash Times: Events w/ >50 PE Flash in BSW. Only >50 PE Flashes Drawn",1000,-10,30)
hall = TH1F("h_flashtimes_all","Flash Times: All Flashes, All Events.",1000,-10,30)
hbig = TH1F("h_flashtimes_big","Flash Times: Flashes > 50 PE, All Events.",1000,-10,30)

BGW_mintime = 3.3
BGW_maxtime = 4.9
if not on_beam_true_off_beam_false:
    BGW_mintime = 3.65
    BGW_maxtime = 5.25
if not on_beam_true_off_beam_false and mc_true_data_false:
    BGW_mintime = 3.2
    BGW_maxtime = 4.8

kept_evts, rej_evts = 0, 0
if entries > 20000: entries = 20000
print 'Looping over %d entries ... '%entries
for x in xrange(entries):
     blah = t.GetEntry(x)

     keep_evt = False
     for i in xrange(len(b)):
        blah = hall.Fill(b.at(i).Time())
        if b.at(i).TotalPE() > 50: blah = hbig.Fill(b.at(i).Time())
        if b.at(i).TotalPE() > 50 and b.at(i).Time() > BGW_mintime and b.at(i).Time() < BGW_maxtime: keep_evt = True

     if not keep_evt: 
     	rej_evts += 1
     	continue

     kept_evts += 1
     for i in xrange(len(b)):
         if b.at(i).TotalPE() > 50: poop = h.Fill(b.at(i).Time())

fout = TFile('flash_time_plot_out_mc=%d_onbeam=%d.root'%(int(mc_true_data_false), int(on_beam_true_off_beam_false)),'RECREATE')
fout.cd()
h.Write()
hall.Write()
hbig.Write()
fout.Close()
f.Close()

print "Number of total events analyzed: %d" % entries
print "Number of events with >50 flash in BSW: %d" % kept_evts
print "Number of events without >50 flash in BSW: %d" % rej_evts