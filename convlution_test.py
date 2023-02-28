import numpy as np
import ROOT
from common_funcs import *

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("rn", type=int)
    args = parser.parse_args()
    rn = args.rn
    fn = "run%i.root" % rn
    if(rn == 0):
        fn = "self_trigger_run0.root"

    #files = [ "run70.root", "run71.root", "run72.root", "run73.root", "run74.root", "run75.root", "run76.root", "run77.root", "run78.root", "run79.root", "run80.root", "run81.root", "run82.root", "run83.root", "run84.root"]
    kernels = np.array([eval(x) for x in open("kernels.txt").readlines()])
    kernels = [k/np.sum(np.abs(k)) for k in kernels]
    kernel2 = np.array([-3.18, -13.18, -40.18, -33.18, -16.18,  -9.18, -10.18, -10.18, -3.18], dtype=np.float64)
    kernel2 -= np.sum(kernel2)/float(len(kernel2))
    kernel2 = kernel2/np.sum(np.abs(kernel2))
    noisey_channels = np.array([ i for i in range(128) if i%4 in [0, 3]])
    clean_channels = np.array([ i for i in range(128) if i%4 in [1, 2]])

    binsx = 625
    binsy = 500
    rangey = (-2, 10)
    rangey2 = (-10, 50)
    rangex = (0,1000)
    h2_noisey = ROOT.TH2D("h_noisey", "Noisey;variance;conv", binsx,rangex[0], rangex[1], binsy, rangey[0], rangey[1])
    h2_clean = ROOT.TH2D("h_clean", "Clean;variance;conv", binsx, rangex[0], rangex[1], binsy,rangey[0], rangey[1])
    h22_noisey = ROOT.TH2D("h2_noisey", "Noisey;variance;conv2", binsx,rangex[0], rangex[1], binsy, rangey2[0], rangey2[1])
    h22_clean = ROOT.TH2D("h2_clean", "Clean;variance;conv2", binsx, rangex[0], rangex[1], binsy,rangey2[0], rangey2[1])

    h23_noisey = ROOT.TH2D("h3_noisey", "Noisey;conv1;conv2", binsx,rangey[0], rangey[1], binsy, rangey2[0], rangey2[1])
    h23_clean = ROOT.TH2D("h3_clean", "Clean;conv1;conv2", binsy, rangey[0], rangey[1], binsy, rangey2[0], rangey2[1])

    for h in [h2_noisey, h2_clean, h22_noisey, h22_clean, h23_noisey, h23_clean]:
        h.SetStats(0)
        for ax in [h.GetXaxis(), h.GetYaxis()]:
            ax.SetTitleFont(132)
            ax.SetLabelFont(132)
            ax.SetTitleSize(0.05)
            ax.SetTitleOffset(0.93)

    h2_noisey.SetStats(0)
    h2_clean.SetStats(0)
    end_limit = 950
    if(rn == 0):
        end_limit = -1
#    h_clean, binsx, binsy = np.histogram2d([], [], bins, ranges)
#    h_noisey, _, _ = np.histogram2d([], [], bins, ranges)
    print(fn)
    stash = []
    for i, wfs in enumerate(read_waveforms("/mnt/sdb/mdaq_dec2022_data/ROOT_file_data/"+fn)):
        if not i%100:
            print(i)
        wfs = wfs[:, :end_limit]
        blines = np.array([find_baseline(wf)[0] for wf in wfs])
        wfs = wfs - blines[:, np.newaxis]
        clean_wfs = wfs[clean_channels]
        noisey_wfs = wfs[noisey_channels]
        for h2, h22, h23, wfs in zip([h2_clean, h2_noisey],[h22_clean, h22_noisey], [h23_clean, h23_noisey],[clean_wfs, noisey_wfs]):
            if(i==0):
                h2.SetStats(0)
            convs = np.array([np.convolve(kernels[0], wf, mode="same") for wf in wfs])
            convs2 = np.array([np.convolve(kernel2, wf, mode="same") for wf in wfs])
            #convs = np.array([np.mean(rolling_window(wf,25), axis=-1) for wf in wfs])
            varis = np.array([np.var(rolling_window(wf, 6), axis=-1) for wf in wfs])
            #varis = np.append(np.zeros((len(varis), 1)), np.diff(varis, axis=1), axis=1)


            start_locs = np.argwhere(varis > 23)
            if(not np.any(start_locs)):
                continue
            starts = np.logical_or(np.abs(np.diff(start_locs,axis=0)[:,-1]) > 2, np.any(np.diff(start_locs,axis=0)[:,:-1], axis=1))
            starts = np.append([True], starts)
            start_locs = start_locs[starts]

            ## Find the END of the event
#            for event, chan, sample in start_locs:
#                end = np.argwhere(vardset[event, chan, sample:] < 7)[:,0]
#                if not len(end) or end[0]+10 > 900:
#                    continue
#                length = end[0] + 20 # Pad with 10 samples at start and end (so 20)
#                sample = sample -34 # Move the start back by 24, then 10
#                final.append((event, chan, sample, sample+length))
#            final = np.array(final)

            start_locs[:, -1] = np.clip(start_locs[:, -1]-20, a_min=0, a_max=None)

            # Remove events that are at the very start
            start_locs = start_locs[start_locs[:,-1] != 0]

            end_locs = np.copy(start_locs)
            end_locs[:,-1] = np.clip(end_locs[:, -1]+25, a_max= wfs.shape[-1], a_min=None)

            xvals = [np.max(varis[chan, sample:end_sample]) for (chan, sample), (_, end_sample) in zip(start_locs, end_locs)]
            yvals = [np.max(convs[chan, sample:end_sample]) for x, (chan, sample), (_, end_sample) in zip(xvals, start_locs, end_locs)]
            yvals2 = [np.max(convs2[chan, sample:end_sample]) for x, (chan, sample), (_, end_sample) in zip(xvals, start_locs, end_locs)]
#
#            temp = np.array(list(zip(xvals,yvals)))
#            cond = np.logical_and(np.logical_and(np.logical_and(temp[:,0] > 60, temp[:, 0] <80), temp[:,1] > 0.09), temp[:, 1] < 0.13)
#            if(rn == 0 and h2 == h2_noisey and np.any(cond)):
#                import ipdb; ipdb.set_trace()
#                stash.append((wfs, convs, varis, start_locs))
#

            for v1, v2, v22 in zip(xvals, yvals, yvals2):
                h2.Fill(v1, v2)
                h22.Fill(v1, v22)
                if(v1 <150):
                    h23.Fill(v2, v22)

#            varis2 = np.array([np.max(rolling_window(x, 5),axis=-1) for x in varis])
#            convs2 = np.array([np.max(rolling_window(x, 5),axis=-1) for x in convs])
#
#            varis2 = np.reshape(varis2, np.prod(varis2.shape))
#            convs2 = np.reshape(convs2, np.prod(convs2.shape))
#
            #h2 += np.histogram2d(varis, convs, bins, ranges)[0]

#            for v,c in zip(varis2, convs2):
#                h2.Fill(v,c)
#                    for v2, v1 in zip(c,v):
#                        h2.Fill(v1,v2)
#            #dset.append(np.array([wfs, convs, varis]))

        #[ [[h2.Fill(v1,v2) for v2, v1 in zip(c,v)] for c,v in zip(convs, varis)] for _, convs, varis in dset]
#
#    file_out = ROOT.TFile("histogram_out2_run%i.root" % rn, "RECREATE")
#    h2_noisey.Write()
#    h2_clean.Write()
#    file_out.Close()

#
    c = ROOT.TCanvas("can1", "can1", 1500, 500)
    c.Divide(2, 1)
    c.cd(1)
    h2_noisey.Draw("colz")
    c.cd(1).SetLogz()
    c.cd(1).Update()

    c.cd(2)
    h2_clean.Draw("colz")
    c.cd(2).SetLogz()
    c.cd(2).Update()

    c2 = ROOT.TCanvas("can2", "can2", 1500, 500)
    c2.Divide(2, 1)
    c2.cd(1)
    h22_noisey.Draw("colz")
    c2.cd(1).SetLogz()
    c2.cd(1).Update()

    c2.cd(2)
    h22_clean.Draw("colz")
    c2.cd(2).SetLogz()
    c2.cd(2).Update()

    c3 = ROOT.TCanvas("can3", "can3", 1500, 500)
    c3.Divide(2, 1)
    c3.cd(1)
    h23_noisey.Draw("colz")
    c3.cd(1).SetLogz()
    c3.cd(1).Update()

    c3.cd(2)
    h23_clean.Draw("colz")
    c3.cd(2).SetLogz()
    c3.cd(2).Update()
