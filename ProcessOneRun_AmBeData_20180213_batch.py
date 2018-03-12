import h5py
import numpy as np
#import matplotlib.pyplot as plt
#import graphUtils
import ROOT
from ROOT import TH1F, TH1D, TH2F, TH2D, TTree, TFile, TCanvas, TF1, TGraphErrors
import datetime
import time
from array import array
import zipfile
import os
import sys
from operator import itemgetter

#ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(1111)

F=[60, 70, 80, 90, 100, 110, 120, 130, 140] #9 values
T=[200, 300, 400, 500, 600, 700, 800, 900, 1000] #9 values
chargerange1 = [5, 15, 25, 35, 45, 55, 65]
hTotalCharge=[]
hFastCharge=[]
hPSD = [ [], [], [], [], [], [], [], [], [] ]
hPSDvsQtot = [ [], [], [], [], [], [], [], [], [] ]
hFOMvsQtot = [ [ TGraphErrors() for j in range(len(F))] for i in range(len(T)) ]
hPSD_subregion =[ [ [TH1F() for k in range(len(chargerange1))] for j in range(len(F))] for i in range(len(T))] # 3d array of histogram

for i in range(len(F)):
    #fast charge histogram
    name = "hFastCharge_F" + str(F[i])
    hist = TH1D(name,"",100,-100,200)
    hist.SetXTitle("Charge (Volts/50 #Omega*0.2 ns)")
    hist.SetYTitle("Counts")
    hFastCharge.append(hist)
for i in range(len(T)):
    #total charge histogram
    name = "hTotCharge_T" + str(T[i])
    hist = TH1D(name,"",100,-100,200)
    hist.SetXTitle("Charge (Volts/50 #Omega*0.2 ns)")
    hist.SetYTitle("Counts")
    hTotalCharge.append(hist)
for i in range(len(T)):
    for j in range(len(F)):
        #PSD distributions, depend on both T and F
        name = "hPSD_T" + str(T[i]) + "_F" + str(F[j])
        hist = TH1D(name,"",100,-0.1,0.5)
        hist.SetXTitle("PSD")
        hist.SetYTitle("Counts")
        hPSD[i].append(hist)
        #PSD vs Qtot scattering plot, depend on both T and F
        name = "hPSDvsQtot_T" + str(T[i]) + "_F" + str(F[j])
        hist = TH2D(name,"",100,-100,200, 100, -0.1, 0.5)
        hist.SetXTitle("Total Charge (Volts/50 #Omega*0.2 ns)")
        hist.SetYTitle("PSD")
        hPSDvsQtot[i].append(hist)
        #FOM vs Qtot scattering plot, depend on both T and F
        name = "hFOMvsQtot_T" + str(T[i]) + "_F" + str(F[j])
        hFOMvsQtot[i][j].SetName(name)
        #hFOMvsQtot[i][j] = TH2F(name,"",600,-100,200, 50, 0, 5)
        #hFOMvsQtot[i][j].SetXTitle("Total Charge (Volts/50 #Omega*0.2 ns)")
        #hFOMvsQtot[i][j].SetYTitle("FOM")
        for k in range(len(chargerange1)):
            name = "hPSDvsQtot_T" + str(T[i]) + "_F" + str(F[j]) + "_QtotRange_" + str(chargerange1[k]) + "_" + str(chargerange1[k]+10)
            hPSD_subregion[i][j][k] = TH1D(name,"",100,-0.1,0.5)
            hPSD_subregion[i][j][k].SetXTitle("PSD")
            hPSD_subregion[i][j][k].SetYTitle("Counts")

#exit()

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

filename = "/media/sf_H_DRIVE/PROSPECT_LiLS_PSD_data/" + sys.argv[1] #LiLS_batch6_sample1/run3572274897.h5"
print filename
rtname = str.replace(filename,".h5","_runAmBe.root")
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
        if cnt%500==0:
            print "evt: ", cnt
        amp, idx = min((amp, idx) for (idx, amp) in enumerate(a))
        base = np.mean(a[idx-300:idx-100:1])
        #calculate fast charge and fill histograms
        for i in range(len(T)):
            totalcharge = base*(50+T[i]) - np.sum(a[idx-50:idx+T[i]:1]) #35, 2490, 2525
            hTotalCharge[i].Fill(totalcharge) # charge in mC
            for j in range(len(F)):
                fastcharge = base*(50+F[j]) - np.sum(a[idx-50:idx+F[j]:1]) #35, 2490, 2525
                if i==0: #only need to fill FastCharge histogram once
                    hFastCharge[j].Fill(fastcharge) # charge in mC
                psd = (totalcharge-fastcharge)/totalcharge
                hPSD[i][j].Fill(psd)
                hPSDvsQtot[i][j].Fill(totalcharge,psd)
                for k in range(len(chargerange1)):
                    if totalcharge>=chargerange1[k] and totalcharge<chargerange1[k]+10:
                        hPSD_subregion[i][j][k].Fill(psd)
        hAmp.Fill(-1000.0*amp)
        hAmpIndex.Fill(idx)
        hBaseline.Fill(base*1000.0)

#calculate FOM, Z, and save plots
# read the waveforms again
#for w in file['Waveforms']:
#    a, b = w[0]
#    if len(a)>0:
#        amp, idx = min((amp, idx) for (idx, amp) in enumerate(a))
#        base = np.mean(a[idx-300:idx-100:1])
for i in range(len(T)):
#    totalcharge = base*(50+T[i]) - np.sum(a[idx-50:idx+T[i]:1]) #35, 2490, 2525
    totalcharge = 0.0
    for j in range(len(F)):
#        fastcharge = base*(50+F[j]) - np.sum(a[idx-50:idx+F[j]:1]) #35, 2490, 2525
        fom = 0.0
        fomerr = 0.0
        Z = 0.0
        for k in range(len(chargerange1)):
            #fit the psd distribution with two Gaussians
            fitfunc = TF1("twogaus","gaus(0)+gaus(1)",-0.1,0.5)
            mean = hPSD_subregion[i][j][k].GetMean()
            hPSD_subregion[i][j][k].Fit("gaus","RQ","",0,mean)
            c1 = hPSD_subregion[i][j][k].GetFunction("gaus").GetParameter(0)
            mu1 = hPSD_subregion[i][j][k].GetFunction("gaus").GetParameter(1)
            width1 = hPSD_subregion[i][j][k].GetFunction("gaus").GetParameter(2)
            hPSD_subregion[i][j][k].Fit("gaus","RQ","",mean,0.4)
            c2 = hPSD_subregion[i][j][k].GetFunction("gaus").GetParameter(0)
            mu2 = hPSD_subregion[i][j][k].GetFunction("gaus").GetParameter(1)
            width2 = hPSD_subregion[i][j][k].GetFunction("gaus").GetParameter(2)
            fitfunc.SetParameters(c1, mu1, width1, c2, mu2, width2)
            hPSD_subregion[i][j][k].Fit("twogaus","QR", "", 0, 0.4)
            mu1 = fitfunc.GetParameter(1)
            mu1err = fitfunc.GetParError(1)
            width1 = fitfunc.GetParameter(2)
            width1err = fitfunc.GetParError(2)
            mu2 = fitfunc.GetParameter(4)
            mu2err = fitfunc.GetParError(4)
            width2 = fitfunc.GetParameter(5)
            width2err = fitfunc.GetParError(5)
#            print T[i], F[j], chargerange1[k], mu1, mu1err, width1, width1err, mu2, mu2err, width2, width2err
            fom = (mu2-mu1)/(width1**2.0+width2**2.0)**0.5
#            err1 = (mu1err**2+mu2err**2)**0.5
#            err2 = ((2.0*width1*width1err)**2+(2.0*width2*width2err)**2)**0.5
#            fomerr = fom*((err1/(mu2-mu1))**2+(err2/(width1**2+width2**2))**2)**0.5
            Z = mu2 - mu1
#        hFOMvsQtot[i][j].Fill(totalcharge, fom)
            hFOMvsQtot[i][j].SetPoint(k,chargerange1[k]+5, fom)
#            hFOMvsQtot[i][j].SetPointError(k,5, fomerr)
                

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
for i in range(len(F)):
    hFastCharge[i].Write()
for i in range(len(T)):
    for j in range(len(F)):
        hPSD[i][j].Write()
        hPSDvsQtot[i][j].Write()
        cFOMvsQtot = TCanvas()
        name = hFOMvsQtot[i][j].GetName()
        name = "c_" + name
        cFOMvsQtot.SetName(name)
        hFOMvsQtot[i][j].Draw("AP")
        cFOMvsQtot.Write()
        cFOMvsQtot.Close()
        for k in range(len(chargerange1)):
            hPSD_subregion[i][j][k].Write()

file.close()
