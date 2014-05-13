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
from ROOT import gROOT, TGaxis,TCanvas, TF1, TGraph, TLine, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt
cscar = TCanvas('cscar','scar',600,1800)
cscar.Draw()
cscar.cd()
pscar1 = TPad('pscar1','p',0.0,0.0,1,.5)
pscar1.SetGrid()
pscar1.Draw()
pscar2 = TPad('pscar2','p',0.0,0.5,1,1)
pscar2.SetGrid()
pscar2.Draw()
pscar2.cd()

first_day = 56506.8942593
last_day = 56673.0953472

time , counts, errs = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/hist_mean_fileout.txt', unpack=True)
time_list = []
sample_x = []
mean = np.mean(counts)
for i in range(len(time)):
	time_list.append(time[i]-time[0])
	sample_x.append((counts[i]-mean)*1000)
sample_mean = np.mean(sample_x)*1.0
x = array("d",time_list)
time_series = np.array(x, np.float64)
y = array("d",sample_x)
measurement = np.array(y, np.float64)
f = np.linspace(0.01, 1, 200)
x1 = array("d",f)
norm = time_series.shape[0]
lombs = sp.signal.lombscargle(time_series , measurement,f)
y1 = array("d",lombs)
gr_1 = TGraph(len(x1), x1, y1)
gr_1.SetTitle("Real Data;frequency [rads/day]; Amplitude [Arbitrary Units]")

gr_1.GetXaxis().SetLimits(0,1)
gr_1.Draw('al')
gr_1.SetLineColor(1)

pscar2.Modified()
pscar2.Update()

pscar1.cd()
variance = 0.0
for k in range(len(sample_x)):
	variance += (sample_x[k] - sample_mean)**2
sample_deviation = np.sqrt((1.0/len(sample_x))*variance)
print(sample_deviation)
time_test = []
sample_test = []
for i in range(len(time)):
	time_test.append(time[i]-time[0])
	sample_test.append(np.random.normal(0, sample_deviation-0.005)+0.089*np.sin(2*np.pi/7.99 *time[i])+0.098*np.sin(2*np.pi/38.37 *time[i])+0.08*np.sin(2*np.pi/28.72 *time[i])+0.087*np.sin(2*np.pi/128.89 *time[i]))
test_mean = np.mean(sample_test)
variance_test = 0.0
for j in range(len(sample_test)):
	variance_test += (sample_test[j] - test_mean)**2
test_deviation = np.sqrt((1.0/len(sample_test))*variance_test)
print(test_deviation)
x2 = array("d",time_test)
t3 = np.array(x2, np.float64)
y2 = array("d",sample_test)
x3 = np.array(y2, np.float64)
f1 = np.linspace(0.01, 1, 200)
lombs1 = sp.signal.lombscargle(t3, x3,f1)
x4 = array("d",f)
y4 = array("d",lombs1)
gr_2 = TGraph(len(x4), x4, y4)
gr_2.SetTitle("Sim Data; Frequency [rads/day]; Amplitude [Arbitrary Units]")
gr_2.GetXaxis().SetLimits(0,1)
gr_2.Draw('ALP')
gr_2.SetLineColor(1)
pscar1.Modified()
pscar1.Update()
'''
cscar1 = TCanvas('cscar2','scar1',600,1800)
cscar1.Draw()
cscar1.cd()
pscar11 = TPad('MYDATA','MYDATA',0.0,0.0,1,.5)
pscar11.SetGrid()
pscar11.Draw()
pscar21 = TPad('SIMDATA','SIMDATA',0.0,0.5,1,1)
pscar21.SetGrid()
pscar21.Draw()

pscar21.cd()
gr_2 = TGraph(len(x2), x2, y2)
gr_2.SetTitle("Simulated Data; Frequency; Amplitude")
gr_2.Set
gr_2.Draw('ap')
pscar21.Modified()
pscar21.Update()

pscar11.cd()
gr_3 = TGraph(len(x), x, y)
gr_3.SetTitle("Real Data; Arbitrary Units; frequency  [rads/day]")
gr_3.Draw('ap')
pscar11.Modified()
pscar11.Update()
'''
raw_input("done")
