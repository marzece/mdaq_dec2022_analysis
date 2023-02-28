# coding: utf-8
import numpy as np
import ROOT
from collections import defaultdict
from common_funcs import *

beam_spot_channels = [4, 1, 26,25, 0, 28, 27, 6, 18, 67, 24, 75, 5, 2, 19,46]
files = [ "run70.root", "run71.root", "run72.root", "run73.root", "run74.root", "run75.root", "run76.root", "run77.root", "run78.root", "run79.root", "run80.root", "run81.root", "run82.root", "run83.root", "run84.root"]
files = ["run55.root"]
path = "/mnt/sdb/mdaq_dec2022_data/ROOT_file_data/"
charge_all = {}
var_all = {}

#CHANNELS = beam_spot_channels
CHANNELS = list(range(96))
#WINDOWS = list(range(2,40))
start = 990
end = 1030
h = ROOT.TH1D("", "", 500, -1, 3)
h.SetStats(0)
h.SetLineColor(ROOT.kBlack)
h.SetLineWidth(2)
h2 = ROOT.TH1D("", "", 100, -1, 3)
h2.SetLineWidth(2)
for fn in files[:1]:
    runid = fn[:5]
    data = []
    spe_wfs = []
    charges = defaultdict(list)
    variances = defaultdict(list)
    sat_charges= defaultdict(list)
    for i, wfs in enumerate(read_waveforms(path+fn)):
        if((i%500) == 0):
            print(i)
        for chan in CHANNELS:
            wf = wfs[chan]
            bline, _ = find_baseline(wf)
            #varis = np.array([np.var(rolling_window(wf, width), axis=-1) for width in WINDOWS])
            #wf = np.clip(wfs[chan], a_min=None, a_max=bline+2)
            # Find the pulse (if there is one) in the window
            pulse_len = 15
            pulse_start = start + np.argmin(wf[start:end]) - 3
            #variance = np.max(varis[:, pulse_start:pulse_start+pulse_len], axis=-1)
            if(np.min(wf[start:end]) == 0):
               sat_charges[chan].append(calculate_charge(wf[pulse_start: pulse_start+pulse_len], bline))
            else:
                charges[chan].append(calculate_charge(wf[pulse_start: pulse_start+pulse_len], bline))
                #variances[chan].append(variance)
    charge_all[runid] = (charges, sat_charges)
    #var_all[runid] = variances

    for chan in CHANNELS:
        if(chan%4 in [0,3]):
                continue
        for v in charges[chan]:
            h.Fill(v)
        for v in sat_charges[chan]:
            h2.Fill(v)
        h.GetXaxis().SetTitle("Charge [pC]")
        h.GetXaxis().SetTitleFont(132)
        h.GetXaxis().SetLabelFont(132)
        h.GetXaxis().SetTitleSize(0.05)
        h.GetYaxis().SetTitle("Counts")
        h.GetYaxis().SetTitleFont(132)
        h.GetYaxis().SetLabelFont(132)
        h.GetYaxis().SetTitleSize(0.05)
#c = ROOT.TCanvas()
#h.Draw()

#    fout = open("SPE_CHARGE_VAR_INFO_%s.dat" % runid, 'w')
#    fout.write(str(dict(charge_all[runid][0])))
#    fout.write("\n")
#    fout.write(str(dict(var_all[runid])))
#    fout.close()
#    c2 = ROOT.TCanvas("can1", "can1", 1500, 1200)
#    c2.Divide(2,2)
#    for idx, wind in enumerate(WINDOWS):
#        h2_noisey = ROOT.TH2D("h22_noisey", "Noisey-%i;Charge [pC];Variance" % wind, 100, -1, 3, 250, 0, 500)
#        h2_clean = ROOT.TH2D("h22_clean", "Clean-%i;Charge [pC];Variance" % wind, 100, -1, 3, 250, 0, 500)
#        h2_clean.SetStats(0)
#        h2_noisey.SetStats(0)
#        for chan in CHANNELS:
#            h2 = h2_noisey
#            if(chan %4  in [1,2]):
#                h2 = h2_clean
#
#            varis = np.array(var_all[runid][chan])[:,idx]
#            
#            for _c, v in zip(charge_all[runid][0][chan], varis):
#                h2.Fill(_c, v)
#
#        for h2 in [h2_clean, h2_noisey]:
#            for ax in [h2.GetXaxis(), h2.GetYaxis()]:
#                ax.SetTitleFont(132)
#                ax.SetLabelFont(132)
#                ax.SetTitleSize(0.05)
#
#        pad = c2.cd(1)
#        pad.SetLogz()
#        h2_noisey.Draw("colz")
#        pad = c2.cd(2)
#        pad.SetLogz()
#        h2_clean.Draw("colz")
#
#        pad = c2.cd(3)
#        pad.SetLogz()
#        h2_noisey.GetYaxis().SetRangeUser(0,40)
#        h2_noisey.Draw("colz")
#
#        pad = c2.cd(4)
#        pad.SetLogz()
#        h2_clean.GetYaxis().SetRangeUser(0,40)
#        h2_clean.Draw("colz")
#
#        c2.Print("VAR_SPE_SCAN.gif+")
#        c2.Update()
#        del h2_clean
#        del h2_noisey
#
#    c = ROOT.TCanvas()
#    c.SetLogy()
#    h2.SetLineColor(2)
#    h2.SetStats(0)
#    if(sat_charges):
#        if(len(sat_charges) >  len(charges)):
#            h2.Draw()
#            h.Draw("same")
#        else:
#            h.Draw()
#            h2.Draw("same")
#        leg = ROOT.TLegend(0.3, 0.6, 0.5,0.9)
#        leg.SetTextFont(132)
#        leg.AddEntry(h, "Unsaturated", "l")
#        leg.AddEntry(h2, "Saturated", "l")
#        leg.SetLineColor(0)
#        leg.Draw("same")
#    else:
#        h.Draw()
#    c.Print(str(runid) + "charge_dist.pdf")

#ks = sorted(charge_all.keys())
#for chan in CHANNELS:
#    g = ROOT.TGraphErrors()
#    g2 = ROOT.TGraphErrors()
#    count1 = 0
#    count2 = 0
#    for i, k in enumerate(ks):
#        chgs = charge_all[k][0][chan]
#        sat_chgs = charge_all[k][1][chan]
#        m1 = np.mean(chgs) if len(chgs)>100 else None
#        std1 = np.std(chgs) if len(chgs)>100 else None
#        m2 = np.mean(sat_chgs) if len(sat_chgs)>100 else None
#        std2 = np.std(sat_chgs) if len(sat_chgs)>100 else None
#        x = 9160 + 200*i if i<=6 else 9160+200*(i-1)
#        if(m1):
#            g.SetPointError(count1, 0, std1/np.sqrt(len(chgs)))
#            g.SetPoint(count1, x, m1)
#            count1+=1
#        if(m2):
#            g2.SetPointError(count2, 0, std2/np.sqrt(len(sat_chgs)))
#            g2.SetPoint(count2, x, m2)
#            count2+=1    
#
#    c = ROOT.TCanvas()
#    mg = ROOT.TMultiGraph()
#    mg.Add(g)
#    mg.Add(g2)
#    mg.Draw("aline*")
#    g2.SetMarkerColor(ROOT.kRed)
#    mg.GetXaxis().SetTitle("LED Intensity [Arb.]")
#    mg.GetYaxis().SetTitle("Observed Charge [pC]")
#    mg.GetYaxis().SetTitleFont(132)
#    mg.GetYaxis().SetLabelFont(132)
#    mg.GetXaxis().SetLabelFont(132)
#    mg.GetXaxis().SetTitleFont(132)
#    mg.GetXaxis().SetTitleSize(0.05)
#    mg.GetYaxis().SetTitleSize(0.05)
#    leg = ROOT.TLegend(0.1, 0.6, 0.4, 0.9)
#    leg.AddEntry(g, "Unsaturated", "lp")
#    leg.AddEntry(g2, "Saturated", "lp")
#    leg.SetLineColorAlpha(0,0)
#    leg.SetFillColorAlpha(0,0)
#    if(count2 > 0):
#        leg.Draw()
#    c.Print("intensity_scan_chan%i.pdf" % chan)
#    yax = mg.GetYaxis()
#    yax.SetRangeUser(0.1, yax.GetXmax())
#    c.SetLogy()
#    c.Print("intensity_scan_chan%i_log.pdf" % chan)
