import numpy as np
import ROOT
from common_funcs import *
from collections import defaultdict
import os
path = "/mnt/sdb/mdaq_dec2022_data/ROOT_file_data/"
self_trigger_fn = "self_trigger_run0.root"
#files = [self_trigger_fn]
#LENGTH_CUTOFF = 19990
#files = [ "run70.root", "run71.root", "run72.root", "run73.root", "run74.root", "run75.root", "run76.root", "run77.root", "run78.root", "run79.root", "run80.root", "run81.root", "run82.root", "run83.root", "run84.root"]
files = sorted([x for x in os.listdir(path) if x[:3] == "run" and  x[-5:] == ".root"])
LENGTH_CUTOFF = 1990
HIST_TOP = 8700
HIST_BOT = 0
count = 0
center = 100
holding_list = []
def new_th2d():
    global count
    print("NEW KEY ", count)
    count += 1
    h2 = ROOT.TH2D("h2_%i" %count ,"" , 800, 0, 800, HIST_TOP-HIST_BOT, HIST_BOT, HIST_TOP)
    holding_list.append(h2)
    return h2

histograms = defaultdict(new_th2d)
for fn in files:
    for i, wfs in enumerate(read_waveforms(path+fn)):
        wfs = wfs[:96, :LENGTH_CUTOFF]
        if(i%100 == 0):
            print(i)
        locs = np.argwhere((wfs-8192 < -50))
        diff_locs = np.diff(locs, axis=0)
        start_locs = np.where(np.logical_or(diff_locs[:,-1] > 50, diff_locs[:,0]))[0]
        start_locs = locs[start_locs]
        for chan, sample in start_locs:
            wf = wfs[chan]
            start = max(0, sample-50)
            end = min(LENGTH_CUTOFF, sample+100)
            min_loc = start + np.argmin(wf[start:end])
            min_val = wf[min_loc]

            window_len = 20
            mean_wf = np.mean(rolling_window(wf[min_loc:], window_len), axis=-1)
            cond1 = np.logical_and(mean_wf > 8191, mean_wf < 8193)
            mean_wf2 = np.roll(mean_wf, -window_len) 
            cond2 = np.logical_and(mean_wf2 > 8191, mean_wf2 < 8193)

            actual_end = wf.shape[-1]
            end_zone = np.where(np.logical_and(cond1, cond2))[0]
            if(end_zone.shape[-1]):
                actual_end = min_loc + end_zone[0]

            template = wf[start:actual_end]
            min_loc = np.argmin(template)

            if (template[min_loc] != min_val):
                continue

            k = np.round(min_val,decimals=-1) 
            for i,samp in enumerate(template):
                histograms[k].Fill(center + (i-min_loc), samp)

#graphs = defaultdict(ROOT.TGraph)
#for k,v in histograms.items():
#    count = 0
#    pf = v.ProfileX()
#    for i in range(v.GetNbinsX()):
#        bin_num = i+1
#        x = v.GetXaxis().GetBinCenter(bin_num)
#        if(x < 55 or x > 500):
#            continue
#        graphs[k].SetPoint(count, x, pf.GetBinContent(i))
#        #vals_y = [v.GetBinContent(bin_num, j) for j in range(1, v.GetNbinsY()+1)]
#        #max_bin = np.argmax(vals_y)+1
#        #graphs[k].SetPoint(count, x, v.GetYaxis().GetBinCenter(max_bin))
#        count+=1
print("Done graphing")

#c = ROOT.TCanvas()
#tline= ROOT.TLine(50, 8191.8, 500, 8191.8)
#tline.SetLineWidth(2)
#tline.SetLineColor(2)
#for k in reversed(sorted(graphs.keys())):
#    graphs[k].SetTitle("Peak = %i" % (8192-k))
#    graphs[k].Draw("aline")
#    graphs[k].GetYaxis().SetRangeUser(8160, 8350)
#    tline.Draw('same')
#
#    c.Update()
#    c.Print("overshoot_graphs_led3.gif+")

#c2 = ROOT.TCanvas()
#g3 = ROOT.TGraph()
#for cnt, (k, g) in enumerate(sorted(graphs.items(), key=lambda x:x[0], reverse=True)):
#    ys = np.array([g.GetY()[i] for i in range(g.GetN())])
#    xs = np.array([g.GetX()[i] for i in range(g.GetN())])
#    region = np.linspace(180,220)
#    dx = np.diff(region[:10])[0]
#    val = np.trapz(np.interp(region, xs, ys-8192),dx=dx)
#    g3.SetPoint(cnt, 8192-k, val)
#g3.Draw("aline")

#line = ROOT.TLine(0, 8192,600, 8192)
#line.SetLineColor(2)
#for k in reversed(sorted(means.keys())):
#    g = plot_it(means[k])
#    g.Draw("aline")
#    line.Draw()
#    c.Update()
#    c.Print("waveformtemplates.gif+")
