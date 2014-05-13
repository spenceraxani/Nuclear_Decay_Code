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

ccb= TCanvas('a', 'Binned Data Canvas',600,1800)
ccb.Draw()
ccb.cd()
pcb1 = TPad('Bin 1','p',0.0,0.0,1,0.5)
pcb1.SetGrid()
pcb1.Draw()
pcb2 = TPad('Counts versus Bin','p',0.0,0.5,1,1)
pcb2.SetGrid()
pcb2.Draw()
pcb2.cd()

first_day = 56547.8942593
last_day = 56673.0953472
juliandate , pressure = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/binned_counts.txt', unpack=True)
y = array("d",pressure)
t = array("d",juliandate)



juliandate , counts = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/binned_counts.txt', unpack=True)
x = array("d",counts)
graph = TGraph(len(x), x, y)
graph.SetTitle(";Counts; Pressure ")
fsignal = TF1("FSignal","pol1", 0.9985 ,1.002)
graph.Fit(fsignal,'R')
kpa_per_counts = fsignal.GetParameter(1)
graph.Draw("Ap")
graph.SetMarkerStyle(6)
print("A Least squares gives: y = "+ str(fsignal.GetParameter(1))+"x " + str(fsignal.GetParameter(0)))
pcb2.Modified()
pcb2.Update()

pcb1.cd()
graph2 = TGraph(len(t), t, y)
graph2.SetTitle(";Counts; Pressure ")
graph2.Draw('al')
fsignal2 = TF1("FSignal","pol1", juliandate[0] ,juliandate[-1])
graph2.Fit(fsignal2,'R')


pcb1.Modified()
pcb1.Update()

#######################
#Now Detrend it
#######################
def div_pressure(time):
	best = fsignal2.GetParameter(1) * time + fsignal2.GetParameter(0)
	return best

detrend_data = []
time_data = []

for i in range(len(juliandate)):
	time_data.append(juliandate[i])
	detrend_data.append(counts[i]-(pressure[i]-div_pressure(juliandate[i]))/5000)
	print(pressure[i]-div_pressure(juliandate[i]))
new_data = array("d",detrend_data)
new_time = array("d",time_data)
detrend_graph = TGraph(len(new_time), new_time, new_data)
detrend_graph.SetTitle(";Normalized Count Rate; Modified Julian Date MST [days]")
fsignal3 = TF1("FSignal","pol1", juliandate[0] ,juliandate[-1])
detrend_graph.Fit(fsignal3,'R')

c2= TCanvas('ad', 'Binned Data Canvas',600,1800)
c2.Draw()
c2.cd()
p1a = TPad('Bin 1d','p',0.0,0.0,1,0.5)
p1a.SetGrid()
p1a.Draw()
p2a = TPad('Bin 1d','p',0.0,0.5,1,1)
p2a.SetGrid()
p2a.Draw()
p2a.cd()
detrend_graph.Draw("Al")
p2a.Modified()
p2a.Update()

p1a.cd()
normal_graph = TGraph(len(t), t, x)
normal_graph.Draw("Al")
fsignal4 = TF1("FSignal","pol1", juliandate[0] ,juliandate[-1])
normal_graph.Fit(fsignal4,'R')
p1a.Modified()
p1a.Update()
raw_input("done")