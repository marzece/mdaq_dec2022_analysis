# coding: utf-8
import numpy as np
import ROOT
from collections import defaultdict
from common_funcs import *

#files = [ "run70.root", "run71.root", "run72.root", "run73.root", "run74.root", "run75.root", "run76.root", "run77.root", "run78.root", "run79.root", "run80.root", "run81.root", "run82.root", "run83.root", "run84.root"]
files = ["run0.root", "run1.root", "run2.root", "run15.root", "run16.root", "run17.root", "run28.root", "run29.root", "run30.root", "run70.root", "run71.root", "run72.root"]
path = "/mnt/sdb/mdaq_dec2022_data/ROOT_file_data/"
charge_all = {}
var_all = {}
conv_all = {}

inner_channels = range(0, 96)
veto_channels = range(96, 128)

# Uncomment below to look at inner or veto PMT channels
#CHANNELS = inner_channels
CHANNELS = veto_channels

# Helper function to create a TH1D and apply my preferred style
def my_new_th1d(name, title, binsx, low, high):
    h = ROOT.TH1D(name, title, binsx,low,high)
    h.SetStats(0)
    h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(2)
    for ax in [h.GetXaxis(), h.GetYaxis()]:
        ax.SetTitleFont(132)
        ax.SetLabelFont(132)
        ax.SetTitleSize(0.05)
    return h

h_noisey = my_new_th1d("h_noisey", "Noisey;Charge [pC];Counts", 500, -1, 3)
h_clean = my_new_th1d("h_clean", "Clean;Charge [pC];Counts", 500, -1, 3)
h_uncut_noisey = my_new_th1d("h_uncut_noisey", "Noisey;Charge [pC];Counts", 500, -1, 3)
h_uncut_clean = my_new_th1d("h_uncut_clean", "Clean;Charge [pC];Counts", 500, -1, 3)
h_uncut_clean.SetLineColor(2)
h_uncut_noisey.SetLineColor(2)

end_limit = 950 # only look at the first 950 samples to avoid any LED induced noise
VARIANCE_THRESHOLD = 23 # Variance threshold for where a SPE pulse is found
for fn in files:
    runid = fn[:fn.find(".")]
    charges = defaultdict(list)
    uncut_charges = defaultdict(list)
    variances = defaultdict(list)
    convs = defaultdict(list)

    for i, wfs in enumerate(read_waveforms(path+fn)):
        if((i%500) == 0):
            print(i)
        for chan in CHANNELS:
            wf = wfs[chan, :end_limit]
            bline, _ = find_baseline(wf)
            conv = np.convolve(kernels[0], wf-bline, mode="same") 
            conv2 = np.convolve(kernel2, wf-bline, mode="same")
            varis = np.var(rolling_window(wf, 6), axis=-1)
            # Find locations where the variance goes above threshold
            locs = np.where(varis > VARIANCE_THRESHOLD)[0]
            # Find only the start of the PMT waveforms, i.e where ever the threshold is crossed first
            start_locs = np.concatenate((locs[0:1], locs[1:][np.diff(locs)>1]))
            for loc in start_locs:
                # Find the start and end of waveform. Unless it's at the start/end of the
                # waveform use a fixed 15 sample window
                start = max(loc-3, 0)
                end = min(len(wf), start+15)
                # Find the largest values for the variance & convolution variables within the pulse window
                c = np.max(conv[start:end])
                c2 = np.max(conv2[start:end])
                v = np.max(varis[start:end])

                # Get the PMT pulse charge
                this_charge = calculate_charge(wf[start:end], bline)

                # Store values before applying cuts/filters
                uncut_charges[chan].append(this_charge)
                variances[chan].append(v)
                convs[chan].append((c,c2))
                # Apply cuts
                cut = cut_value(v)
                if(v < 180 and (c>cut or c2 < 3)):
                    continue
                charges[chan].append(this_charge)

    # Store pulse values, indexed by run, for later inspection
    charge_all[runid] = (charges, uncut_charges)
    var_all[runid] = variances
    conv_all[runid] = convs

    # Fill histograms
    for chan in CHANNELS:
        h = h_clean
        h_uncut = h_uncut_clean
        if(chan%4 in [0,3]):
            h = h_noisey
            h_uncut = h_uncut_noisey
        for v in charges[chan]:
            h.Fill(v)
        for v in uncut_charges[chan]:
            h_uncut.Fill(v)

c = ROOT.TCanvas("can1", "can1", 1500, 500)
c.Divide(2,1)
pad = c.cd(1)
pad.SetLogy()
h_noisey.Draw()
h_uncut_noisey.Draw("same")
pad = c.cd(2)
pad.SetLogy()
h_clean.Draw()
h_uncut_clean.Draw("same")
c2 = ROOT.TCanvas()
c2.SetLogy()
hsum = my_new_th1d("hsum", "Veto PMTs;Charge [pC];Counts", 500, -1,3)
hsum.Add(h_clean)
hsum.Add(h_noisey)
hsum.Draw();
