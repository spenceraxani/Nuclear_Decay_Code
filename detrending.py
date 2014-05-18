import os
import glob
import math
from math import log10, floor
import time
import sys
import scipy.signal
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy
from ROOT import gROOT, TCanvas, TF1, TGraph, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt
from scipy.stats import *
from scipy.stats.stats import pearsonr


c1= TCanvas('c1', '',600,1800)
c1.Draw()
c1.cd()
c1p2 = TPad('Bin 1','p',0.0,0.0,1,0.5)
c1p2.SetGrid()
c1p2.Draw()
c1p1 = TPad('Counts versus Bin','p',0.0,0.5,1,1)
c1p1.SetGrid()
c1p1.Draw()
c1p1.cd()
'''
displacement_10 = [3.688, 7.376,11.06,14.75]
displacement_09 = [5.506, 11.01,16.52,22.02]
pressure = [1,2,3,4]
pressure = array("d",pressure)
displacement_10 = array("d",displacement_10)
displacement_09 = array("d",displacement_09)

dis_gr = TGraph(len(pressure), pressure, displacement_10)
dis_gr.Draw("ap")
dis_gr.SetTitle(";Pressure [kPa]; Displacement [\mu m]")
fsignal = TF1("FSignal","pol1", pressure[0] ,pressure[-1])
dis_gr.Fit(fsignal,'R')
dis_gr.SetMarkerColor(4)
dis_gr.SetMarkerSize(1)
dis_gr.SetMarkerStyle(21)
c1p1.Modified()
c1p1.Update()
'''
date , counts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_total_counts.txt', unpack=True)
date , pressure = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_pressure.txt', unpack=True)
date , temperature = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_temperature.txt', unpack=True)
mean_pressure = np.mean(pressure)
mean_temperature = np.mean(temperature)
try:
	os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended.txt')
	os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended_b.txt')
except OSError:
	pass
outfile= open('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended.txt','a')
outfileb= open('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended_b.txt','a')
x1 = array("d",date)
y1 = array("d",temperature)
print(mean_pressure)
print(mean_temperature)
det_counts = []
print(counts)
dis_to_pressure = 3.687
dis_to_temperature = 22.2*0.142000
eff_to_dis = 0.000051#508857
for i in range(len(date)):
	det_counts.append(counts[i]*(1.0  - (pressure[i]-mean_pressure)*dis_to_pressure*eff_to_dis + (temperature[i]-mean_temperature)*dis_to_temperature*eff_to_dis))
	outfile.write(str(date[i]) + " \t" +str(counts[i]*(1.0  - (pressure[i]-mean_pressure)*dis_to_pressure*eff_to_dis + (temperature[i]-mean_temperature)*dis_to_temperature*eff_to_dis))+"\n" )
new_date = []
print(len(det_counts))
for j in range(len(date)):
	new_date.append(date[j]-date[0])
for k in range(len(new_date)):
	outfileb.write(str(new_date[k]) + " \t" +str(det_counts[k])+"\n" )
x = array("d",new_date)
y = array("d",det_counts)
dis_gr = TGraph(len(x), x, y)
dis_gr.Draw("ap")
dis_gr.SetTitle(";date; counts")
fsignal = TF1("FSignal","expo", new_date[0] ,new_date[-1])
dis_gr.Fit(fsignal,'R')
dis_gr.SetMarkerColor(4)
dis_gr.SetMarkerSize(1)
#dis_gr.SetMarkerStyle(2)
c1p1.Modified()
c1p1.Update()

c1p2.cd()
sec_gr = TGraph(len(x1), x1, y1)
sec_gr.Draw("apl")
c1p2.Modified()
c1p2.Update()
raw_input("done")
