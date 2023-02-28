# This script is for looking at the spectrogram for various waveforms
import numpy as np 
import ROOT
from scipy import signal
from common_funcs import *
canvas_index = 0
def gramit(wf):
    spectrogram = signal.spectrogram(wf,nperseg=50)
    t = spectrogram[1]
    f = spectrogram[0]
    xax, yax = np.meshgrid(t,f)
    xax = xax.reshape(np.prod(xax.shape))
    yax = yax.reshape(np.prod(yax.shape))
    gram = spectrogram[2].reshape(np.prod(spectrogram[2].shape))
    g = ROOT.TGraph2D()
    for i, (x, y, v) in enumerate(zip(xax,yax, gram)):
        g.SetPoint(i, x, y, v)
    return g

def showme(wf):
    global canvas_index
    g = gramit(wf)
    c = ROOT.TCanvas("canvas %i" % canvas_index, "canvas", 1500, 500)
    canvas_index += 1
    c.Divide(2,1)
    c.cd(1)
    g.Draw("col")

    c.cd(2)
    g2 = plot_it(wf)
    g2.Draw("aline")

    return c, g, g2

if __name__ == "__main__":
    #fn = "run70.root"
    files = ["run0.root", "run70.root", "run71.root", "run72.root", "run73.root", "run74.root", "run75.root", "run76.root", "run77.root", "run78.root", "run79.root", "run80.root", "run81.root", "run82.root", "run83.root", "run84.root"]
    fn = files[0]
    print(fn)
    noisey_channels = np.array([ i for i in range(128) if i%4 in [0, 3]])
    clean_channels = np.array([ i for i in range(128) if i%4 in [1, 2]])
    clean_led_dset = []
    noisey_led_dset = []
    clean_dset = []
    noisey_dset = []
    NPERSEG = 50
    FFT_MODE = "magnitude"
    for i, wfs in enumerate(read_waveforms(fn)):
        wflen = wfs.shape[-1]
        noisey_dset.append(wfs[noisey_channels, :900])
        clean_dset.append(wfs[clean_channels, :900])
        clean_led_dset.append(wfs[clean_channels,:1900])
        noisey_led_dset.append(wfs[noisey_channels,:1900])

    hist = np.zeros(26)
    hratio1 = ROOT.TH1D("hr1", "hr1", 500, 0, 1)
    hratio2 = ROOT.TH1D("hr2", "hr2", 500, 0, 1)
    hargmax1 = ROOT.TH1D("hargmax1", "hargmax1", 26, 0, 26)
    hargmax2 = ROOT.TH1D("hargmax2", "hargmax2", 26, 0, 26)
    f_ax,  t_ax, _ = signal.spectrogram(noisey_dset[0][0], nperseg=NPERSEG, mode=FFT_MODE)
    for event in noisey_dset:
        for wf in event:
            spect = signal.spectrogram(wf, nperseg=NPERSEG, mode=FFT_MODE)[-1]
            ratios = np.sum(spect[17:25, :], axis=0)/ np.sum(spect, axis=0)
            argmaxes = np.argmax(spect,axis=0)
            hist+= np.sum(spect, axis=1)
            for v in argmaxes:
                hargmax1.Fill(v)
            for v in ratios:
                hratio1.Fill(v)

    hist2 = np.zeros(26)
    for event in clean_dset:
        for wf in event:
            spect = signal.spectrogram(wf, nperseg=NPERSEG, mode=FFT_MODE)[-1]
            hist2+= np.sum(spect, axis=1)
            ratios = np.sum(spect[17:25, :], axis=0)/ np.sum(spect, axis=0)
            hist2+= np.sum(spect, axis=1)
            argmaxes = np.argmax(spect,axis=0)
            for v in argmaxes:
                hargmax2.Fill(v)
            for v in ratios:
                hratio2.Fill(v)

    hratio3 = ROOT.TH1D("hr3", "hr3", 500, 0, 1)
    for event in clean_led_dset:
        for wf in event:
            spect = signal.spectrogram(wf, nperseg=NPERSEG,mode=FFT_MODE )[-1]
            ratios = np.sum(spect[17:25, :], axis=0)/ np.sum(spect, axis=0)
            for v in ratios:
                hratio3.Fill(v)

    hratio4 = ROOT.TH1D("hr4", "hr4", 500, 0, 1)
    for event in noisey_led_dset:
        for wf in event:
            spect = signal.spectrogram(wf, nperseg=NPERSEG, mode=FFT_MODE)[-1]
            ratios = np.sum(spect[17:25, :], axis=0)/ np.sum(spect, axis=0)
            for v in ratios:
                hratio4.Fill(v)

    h1 = ROOT.TH1D("h1", "h1", 26,0,26)
    h2 = ROOT.TH1D("h2", "h2", 26,0,26)
    for i,v in enumerate(hist):
        h1.SetBinContent(i+1,v)
    for i,v in enumerate(hist2):
        h2.SetBinContent(i+1,v)
