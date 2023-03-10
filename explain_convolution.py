# coding: utf-8
import os

from common_funcs import *
gen = read_waveforms("self_trigger_run0.root")
wfs = next(gen)
#wf = wfs[0,4050:]-8192

wf = wfs[3,1950:]-8192
#k = kernels[0]
k = kernel2
klen = len(k)

try:
    os.remove("convolution_explain.gif")
except OSError:
    pass

c = ROOT.TCanvas("can", "", 1400, 1000)
c.GetFrame().SetLineColor(0)
c.GetFrame().SetLineWidth(0)
c.Update()
c.Modify()
MAX = 200
integral = []
for i in range(MAX-klen):
    print(i, MAX-klen)
    temp = wf[i:i+klen]
    g = plot_it(wf[:MAX+50])
    g.SetLineColor(ROOT.kGray)
    prod = k*temp
    g2 = ROOT.TGraph();
    for j in range(klen):
        g2.SetPoint(j, i+j, 5*prod[j]+11);
    this_integral = np.trapz(prod)

    g5 = ROOT.TGraph()
    g5.SetPoint(0, i+klen/2, this_integral+15)

    integral.append(this_integral)

    g3 = ROOT.TGraph()
    for j in range(klen):
        g3.SetPoint(j,i+j, 5*k[j]+7)

    g4 = ROOT.TGraph()
    for j, v in enumerate(integral):
        g4.SetPoint(j, j+klen/2, v+15)

    g6=  ROOT.TGraph()
    for j, v in enumerate(temp):
        g6.SetPoint(j, i+j, v)




    mg= ROOT.TMultiGraph()

    mg.Add(g)
    mg.Add(g2)
    mg.Add(g3)
    mg.Add(g4)
    mg.Add(g6)
    g3.SetLineColor(ROOT.kBlue)
    g4.SetLineColor(ROOT.kGreen)
    g2.SetLineColor(ROOT.kRed)
    g.SetLineWidth(2)
    g2.SetLineWidth(2)
    g3.SetLineWidth(2)
    g4.SetLineWidth(2)
    g5.SetMarkerSize(3)
    mg.Draw('aline')
    mg.GetYaxis().SetNdivisions(0)
    #mg.GetYaxis().SetAxisColor(0)
    g5.Draw("same*")
    mg.GetXaxis().SetRangeUser(0,200)
    mg.GetYaxis().SetRangeUser(-25, 20)
    mg.GetXaxis().SetTitle("Time [ticks]")
    l1 = ROOT.TLine(i, -6, i, 13);
    l2 = ROOT.TLine(i+klen, -6, i+klen, 13);
    l1.Draw()
    l2.Draw()

    leg = ROOT.TLegend(0.7,0.67,1.0, 0.97)
    leg.SetLineColor(0)
    leg.AddEntry(g, "Waveform", "l")
    leg.AddEntry(g3, "Template", "l")
    leg.AddEntry(g2, "Waveform * Template", "l")
    leg.AddEntry(g5, "#int Waveform * Template", "p")
    leg.AddEntry(g4, "Convolution", "l")
    leg.Draw()


    c.Update()
    c.Print("convolution_explain.gif+")

