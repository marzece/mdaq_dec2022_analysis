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

CHANNELS = inner_channels

#def cut_value(var):
#    return 8.425/44.314*var - 0.901

def cut_value(var):
    return 0.040559*var +0.516458

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

kernels = np.array([eval(x) for x in open("kernels.txt").readlines()])
kernels = [k/np.sum(np.abs(k)) for k in kernels]

kernel2 = [-0.11557788944719505, 0.025125628140813205, -0.4045226130656374, 0.22613065326549986, 0.15326633165750536, -0.04648241205995873, -2.1469849246232116, -8.247487437185555, -12.341708542713604, -7.153266331658415, -3.4874055415621115, -3.1368221941993397, -2.892676767676676, -1.4266750948163462, -0.9816687737038592, -1.2249683143218135, -0.840966921119616, -0.3613231552162688, -0.5318471337577648, -0.5650510204077364]
#kernel2 = np.array([-3.18, -13.18, -40.18, -33.18, -16.18,  -9.18, -10.18, -10.18, -3.18], dtype=np.float64)
kernel2 -= np.sum(kernel2)/float(len(kernel2))
kernel2 = kernel2/np.sum(np.abs(kernel2))

#end_limit = -1
end_limit = 950
saved_wfs = defaultdict(list)
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
            locs = np.where(varis > 23)[0]
            start_locs = np.concatenate((locs[0:1], locs[1:][np.diff(locs)>1]))
            for loc in start_locs:
                start = max(loc-3, 0)
                end = min(len(wf), start+15)
                c = np.max(conv[start:end])
                c2 = np.max(conv2[start:end])
                v = np.max(varis[start:end])
                this_charge = calculate_charge(wf[start:end], bline)
                uncut_charges[chan].append(this_charge)
                variances[chan].append(v)
                convs[chan].append((c,c2))
                cut = cut_value(v)
                if(v < 180 and (c>cut or c2 < 3)):
                    continue
                charges[chan].append(this_charge)
                if(this_charge < 0.5 and this_charge > 0.3 and ((chan%4) in [1,2])):
                    saved_wfs[chan].append(wf[ max(start-15, 0):min(end+15, len(wf))])

    charge_all[runid] = (charges, uncut_charges)
    var_all[runid] = variances
    conv_all[runid] = convs

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
