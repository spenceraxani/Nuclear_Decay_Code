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

def exp_counts(A,C,t):
	counts_exp = np.exp(A)*numpy.exp((t)*C)
	return(counts_exp)

juliandate , net = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended.txt', unpack=True)

juliandate1 , total_counts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_total_counts.txt', unpack=True)


first_day = juliandate[0]
last_day = juliandate[-1]
time = []
time2 = []
for i in range(len(juliandate)):
	time.append(juliandate[i] - first_day)
	time2.append(juliandate[i] - first_day)
x = array("d",time)
y = array("d",net)

x1 = array("d",time2)
y1 = array("d",total_counts)

c = TCanvas('canvas','can1',600,1800)
c.Draw()
c.cd()
p1 = TPad('MYDATA','MYDATA',0.0,0.0,1,.5)
p1.SetGrid()
p1.Draw()
p2 = TPad('SIMDATA','SIMDATA',0.0,0.5,1,1)
p2.SetGrid()
p2.Draw()

p1.cd()

gr = TGraph(len(x), x, y)
gr.Draw('al')
fitsignal = TF1("FSignal","expo", first_day-first_day, last_day-first_day) # this gaussian can then be integrated to get the two tailed p-values. Input these values into correlation_gaussian.np code to integrate.
gr.Fit(fitsignal,'R')

p2.cd()
gr2 = gr = TGraph(len(x1), x1, y1)
gr2.Draw('al')
fitsignal2 = TF1("FSignal2","expo", first_day-first_day, last_day-first_day) # this gaussian can then be integrated to get the two tailed p-values. Input these values into correlation_gaussian.np code to integrate.
gr2.Fit(fitsignal2,'R')


residual = []
old_residual = []
for i in range(len(x)):
	residual.append(y[i] - exp_counts(fitsignal.GetParameter(0),fitsignal.GetParameter(1),x[i]))
	old_residual.append(y1[i] - exp_counts(fitsignal2.GetParameter(0),fitsignal2.GetParameter(1),x[i]))
print(residual)

fileout = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended_residual.txt'
fileout2 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_residual.txt'
try:
	os.remove(fileout)
	os.remove(fileout2)
except OSError:
	pass

outfile = open(fileout,'a')
outfile2 = open(fileout2,'a')
for j in range(len(juliandate)):
	outfile.write(str(juliandate[j]) + " \t" +str(residual[j])+"\n" )
	outfile2.write(str(juliandate[j]) + " \t" +str(old_residual[j])+"\n" )

raw_input("done")
