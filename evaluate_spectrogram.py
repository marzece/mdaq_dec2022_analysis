# This script is for looking at the spectrogram for various waveforms
import numpy as np 
import ROOT
from scipy import signal
from common_funcs import *

if __name__ == "__main__":
    #fn = "run70.root"
    files = ["run0.root", "run70.root", "run71.root", "run72.root", "run73.root", "run74.root", "run75.root", "run76.root", "run77.root", "run78.root", "run79.root", "run80.root", "run81.root", "run82.root", "run83.root", "run84.root"]
    fn = files[0]
    print(fn)
    clean_led_dset = []
    noisey_led_dset = []
    clean_dset = []
    noisey_dset = []
    NPERSEG = 50
    end_limit = -1
    FFT_MODE = "magnitude"
    binsx = 625
    binsy = 500
    rangex = (0,1000)
    rangey = (0, 10)
    h2_clean = ROOT.TH2D("h_clean", "clean;variance;conv", binsx, rangex[0], rangex[1], binsy,rangey[0], rangey[1])
    h2_noisey = ROOT.TH2D("h_clean", "clean;variance;conv", binsx, rangex[0], rangex[1], binsy,rangey[0], rangey[1])
    for i, wfs in enumerate(read_waveforms(fn)):
        #spect = signal.spectrogram(wf, nperseg=NPERSEG, mode=FFT_MODE)[-1]
        #ratios = np.sum(spect[17:25, :], axis=0)/ np.sum(spect, axis=0)

        clean_wfs = wfs[clean_channels, :end_limit]
        clean_wfs = clean_wfs - 8192
        noisey_wfs = wfs[noisey_channels, :end_limit]
        noisey_wfs = noisey_wfs - 8192
        for h2, wfs in zip([h2_clean, h2_noisey],[clean_wfs, noisey_wfs]):
            #convs = np.array([signal.spectogram(kernels[0], wf, mode="same") for wf in wfs])
            varis = np.array([np.var(rolling_window(wf, 25), axis=-1) for wf in wfs])


            start_locs = np.argwhere(varis > 10)
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
            end_locs[:,-1] = np.clip(end_locs[:, -1]+50, a_max= wfs.shape[-1], a_min=None)
            import ipdb;ipdb.set_trace()

            for (chan, sample), (_, end_sample) in zip(start_locs, end_locs):
                if((sample - end_sample) != 50):
                    continue
                fft = np.fft(wfs[chan, sample:end_sample]

