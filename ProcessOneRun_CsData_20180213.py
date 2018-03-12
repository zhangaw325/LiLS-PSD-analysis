import h5py
import numpy as np
#import matplotlib.pyplot as plt
#import graphUtils
import ROOT
from ROOT import TH1F, TH1D, TTree, TFile, TCanvas
import datetime
import time
from array import array
import zipfile
import os
import sys
from operator import itemgetter

#ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(1111)

T=[200, 300, 400, 500, 600, 700, 800, 900, 1000]
hTotalCharge=[]
for i in range(len(T)):
    name = "hTotCharge_T" + str(T[i])
    hist = TH1D(name,"",600,-100,200)
    hist.SetXTitle("Charge (V/50#Ohm*0.2ns)")
    hist.SetYTitle("Counts")
    hTotalCharge.append(hist)

hBaseline = TH1D("hBaseline","",100,-20,-10)
hBaseline.SetXTitle("Baseline (mV)");
hBaseline.SetYTitle("Counts");

hAmp = TH1D("hAmplitude","",800,0,800)
hAmp.SetXTitle("Amplitude (mV)")
hAmp.SetYTitle("Counts")

hAmpIndex = TH1D("hAmpIndex","",100,2480,2580)
hAmpIndex.SetXTitle("Samping point where pulse reaches maximum")
hAmpIndex.SetYTitle("Counts")
#hIntChargeInBaseRegion = TH1F("hIntChargeInBaseRegion","",5000,0.0,10000.0)

#hIntChargeInBaseRegion.SetXTitle("Charge (mV)");
#hIntChargeInBaseRegion.SetYTitle("Counts");

filename = "/media/sf_H_DRIVE/PROSPECT_LiLS_PSD_data/LiLS_batch21_sample1/run3576493547.h5"
print filename
rtname = str.replace(filename,".h5","_runCs137.root")
file = h5py.File(filename,'r')

rtfile = TFile(rtname,"recreate")

ph, ct = [], [] #array of scope reaadout, clock time
subwave = []
cnt = 0
cWaveCheck = TCanvas()
for w in file['Waveforms']:
    a, b = w[0]
    if len(a)>0:
        cnt = cnt+1
        amp, idx = min((amp, idx) for (idx, amp) in enumerate(a))
        base = np.mean(a[idx-300:idx-100:1])
        for i in range(len(T)):
            charge = base*(50+T[i]) - np.sum(a[idx-50:idx+T[i]:1]) #35, 2490, 2525
            hTotalCharge[i].Fill(charge) # charge in mC
        hAmp.Fill(-1000.0*amp)
        hAmpIndex.Fill(idx)
        hBaseline.Fill(base*1000.0)

#c = TCanvas()
#hAmp.Draw()
hAmp.Write()
#c.SaveAs("testAmp.png")
#c2 = TCanvas()
#hAmpIndex.Draw()
hAmpIndex.Write()
#c2.SaveAs("testAmpIndex.png")
#cBase = TCanvas()
#hBaseline.Draw()
hBaseline.Write()
#cBase.SaveAs("testBaseline.png")
#cCharge = TCanvas()
#hCharge.Draw()
for i in range(len(T)):
    hTotalCharge[i].Write()
#cCharge.SaveAs("testCharge.png")
#cCharge.SaveAs("testCharge.root")
#cWaveCheck.SaveAs("testWaveforms.png")
#cWaveCheck.Write()
#hIntChargeInBaseRegion.Write()
#cWaveCheck.SaveAs("testWaveforms.root")

#print("Keys: %s" % file.keys())
#waveforms = list(file[file.keys()[1]])
#timestamp = list(file[file.keys()[2]])#

#nEvents = len(waveforms)
#print nEvents
#for ievt in range(0,nEvents,1):
#    print timestamp[ievt]

#print waveforms[1]

file.close()
