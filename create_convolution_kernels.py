# This script is for finding the convolution kernels for noise waveforms
import numpy as np 
import ROOT
from scipy import signal
from common_funcs import *

def get_correlation_lag(wf1, wf2, mode):
    corr = signal.correlate(wf1, wf2, mode)
    xax = correlation_lags(len(wf1), len(wf2), mode)
    return xax[np.argmax(corr)]

if __name__ == "__main__":
    #fn = "run70.root"
    files = ["run0.root", "run70.root", "run71.root", "run72.root", "run73.root", "run74.root", "run75.root", "run76.root", "run77.root", "run78.root", "run79.root", "run80.root", "run81.root", "run82.root", "run83.root", "run84.root"]
    library = []
    library_locs = []
    for fn in files:
        print(fn)
        kernel = [0,-1,2,-2,3,-3,2,2,-3,3,-2,2,-1,0]
        pos_sum = sum([x for x in kernel if x>0])
        neg_sum = sum([-x for x in kernel if x < 0])
        kernel = [(x/float(pos_sum+neg_sum))*(neg_sum if x >0 else pos_sum) for x in kernel]
        #kernel =[0.5,-1,0.5]
        noisey_channels = np.array([ i for i in range(128) if i%4 in [0, 3]])
        clean_channels = np.array([ i for i in range(128) if i%4 in [1, 2]])
        dset = []
        for i, wfs in enumerate(read_waveforms(fn)):
            wflen = wfs.shape[-1]
            noisey_dset = wfs[noisey_channels, :900]
            dset.append(noisey_dset)
        dset = np.array(dset)
        vardset = np.array([[np.var(rolling_window(wf, 25),axis=-1) for wf in event] for event in dset])
        maxdset = np.array([[np.max(rolling_window(wf-8192, 25),axis=-1) for wf in event] for event in dset])
        # Box of interest, max in window > 7 varianec > 15, ends when variance < 7
        start_locs = np.argwhere(np.logical_and(vardset>=15, maxdset >= 6))
        starts = np.logical_or(np.abs(np.diff(start_locs,axis=0)[:,-1]) > 2, np.any(np.diff(start_locs,axis=0)[:,:2], axis=1))
        starts = np.append([True], starts)

        start_locs = start_locs[starts]
        start_locs = start_locs[start_locs[:,-1] >= 34]


        # I can't think of a clever numpy way of doing this, so I'm just gonna use a for loop
        final = []
        for event, chan, sample in start_locs:
            end = np.argwhere(vardset[event, chan, sample:] < 7)[:,0]
            if not len(end) or end[0]+10 > 900:
                continue
            length = end[0] + 20 # Pad with 10 samples at start and end (so 20)
            sample = sample -34 # Move the start back by 24, then 10
            final.append((event, chan, sample, sample+length))
        final = np.array(final)

        for event, chan, start, end in final:
            library_locs.append((event, chan, start, end))
            library.append(dset[event, chan, start:end])

    library = [x-np.mean(x) for x in library]
    standard = library[0]
    mode = "same"

    lags = []
    for wf in library:
        lags.append(get_correlation_lag(standard, wf, mode))
    graphs = []
    for lag, wf in zip(lags, library):
        g=  ROOT.TGraph()
        for i,v in enumerate(wf):
            g.SetPoint(i,i+lag, v)
        graphs.append(g)

    points = []
    for lag, wf in zip(lags,library):
        for i,v in enumerate(wf):
            points.append((i+lag,v))
    points = np.array(points)
    h2, binsx, binsy = np.histogram2d(points[:,0], points[:,1], [60, 100], [(0,60), (-50,50)])
    best_wf = binsy[np.argmax(h2, axis=1)]
    mean_wf = [np.sum([binsy[j]*h2[i,j]/np.sum(h2[i,:]) for j in range(100)]) for i in range(60)]

    # Cut the fat (values determined from ocular inspection)
#    best_wf = best_wf[20:45]
#    mean_wf = mean_wf[20:45]

    #"normalize" so it's DC balanced
    best_wf -= np.sum(best_wf)/len(best_wf)
    mean_wf -= np.sum(mean_wf)/len(mean_wf)
