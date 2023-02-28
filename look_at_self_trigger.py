import numpy as np
import ROOT
from common_funcs import *


if __name__ == "__main__":
    if(True):
        fn = "self_trigger_run0.root"
        print(fn)
        f = ROOT.TFile(fn, "READ")
        t= f.Get("tree")
        #hists_v = [ROOT.TH1D("hv%i" %i, "Chan%i;Variance;Counts"%i,500,0,50) for i in range(128)]
        #hists_avg = [ROOT.TH1D("hmean%i" %i, "Chan%i;Mean;Counts"%i,200,-10,10) for i in range(128)]
        #hists2d = [ROOT.TH2D("h2%i" % i, "Chan%i;Variance;Mean" %i, 500, 0,50, 200,-10,10) for i in range(128)]
    #    data = []
    #    hists_v = np.array([np.histogram([], 500, (0,50))[0] for _ in range(128)])
    #    hists_avg = np.array([np.histogram([], 500, (-50,50))[0] for _ in range(128)])
        #hists2d = np.array([np.histogram2d([], [], (500, 200), ((0, 50),(-10,10)))[0] for _ in range(128)])
        bincounts = [200, 200, 30, 30]
        ranges = [(0,50), (-10,10), (0,30), (0,30)]
        histsdd = np.array([np.histogramdd(np.zeros((0,4)), bincounts, ranges)[0] for _ in range(128)])
        #ops = [np.var, np.mean, np.min, np.max]
        dset = []
        for i, ev in enumerate(t):
            print(i)
            wfs = np.reshape(list(ev.FADC), (128, ev.nsample))
            wfs = wfs[:,4000:19900]
            wfs = wfs - 8192
#            data = np.array([np.var(rolling_window(x,25), axis=-1) ,
#                             np.mean(rolling_window(x,25), axis=-1),
#                             np.min(rolling_window(x,25), axis=-1),
#                             np.max(rolling_window(x,25), axis=-1)] for x in wfs]
            for i, x in enumerate(wfs):
                varwf = np.var(rolling_window(x,25), axis=-1)
                meanwf = np.mean(rolling_window(x,25), axis=-1)
                minwf = np.min(rolling_window(x,25), axis=-1)
                maxwf = np.max(rolling_window(x,25), axis=-1)
                histsdd[i] += np.histogramdd([varwf, meanwf, minwf, maxwf], bincounts, ranges)[0]

        np.save("multid_wf_data.npy", np.array(histsdd))

    #        varwfs = np.array([ for x in wfs])
    #        meanwf = np.array([np.mean(rolling_window(x,25)-8192, axis=-1) for x in wfs])
#            histsdd += np.array([np.histogramdd( [np.var(rolling_window(x,25), axis=-1) ,
#                                                  np.mean(rolling_window(x,25), axis=-1) ,
#                                                  np.min(rolling_window(x,25), axis=-1) ,
#                                                  np.max(rolling_window(x,25), axis=-1)],
#                                                  bincounts, ranges)[0] for x in wfs])
#            histsdd += np.array([np.histogramdd([varwf,
#                                                meanwf,
#                                                minwf,
#                                                maxwf],
#                                                (500, 200, 50, 50), ((0, 50),(-10,10), (-50,0), (0, 50)))[0] for x in wfs])
    #        for vwf, mwf, hv, hm, h2 in zip(varwfs, meanwf, hists_v, hists_avg, hists2d):
    #            for v, m in zip(vwf, mwf):
    #                h2.Fill(v, m)
    else:

        hist2d = np.load("2d_var_mean_histogram.npy")
        _, binsx, binsy= np.histogram2d([], [], (500, 200), ((0,50),(-10,10)))
        binsx
        binsx = np.mean(np.dstack((binsx[0:500], binsx[1:])),axis=-1)[0]
        binsy = np.mean(np.dstack((binsy[0:200], binsy[1:])),axis=-1)[0]
        root_hist2d = [ROOT.TH2D("h2%i" %i, "Chan%i;Variance;Mean" %i, 500, 0, 50, 200, -10, 10) for i in range(128)]
        for chan, h in enumerate(root_hist2d):
            for i, xb in enumerate(binsx):
                for j, yb in enumerate(binsy):
                    h.Fill(xb,yb, hist2d[chan, i,j])


        histsum1 = ROOT.TH2D("hs_1", "Sum1;Variance;Mean", 500, 0, 50, 200, -10, 10)
        histsum2 = ROOT.TH2D("hs_2", "Sum2;Variance;Mean", 500, 0, 50, 200, -10, 10)
            
        for chan, h in enumerate(root_hist2d):
            if(chan%4 in [0,4]):
                continue
            histsum1.Add(h)


        c = ROOT.TCanvas()
        histsum1.Draw("colz")
        for chan, h in enumerate(root_hist2d):
            if(chan%4 not in [0,4]):
                continue
            histsum2.Add(h)

        c2 = ROOT.TCanvas()
        histsum2.Draw("colz")
