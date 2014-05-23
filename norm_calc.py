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
from ROOT import gROOT, TGaxis, TCanvas, TF1, TGraph, TLine, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt
import time
gStyle.SetOptStat("nei")
def exp_counts(t):
	counts_exp = numpy.exp(6.93518037533)*numpy.exp(-t*0.00036002491846)
	return(counts_exp)

def fit_counts(B,A,t):
	counts_exp = A * numpy.exp( B * t)
	return(counts_exp)

juliandate , net, err = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_data_detrended_with_errors_copy.txt', unpack=True)
#juliandate , net  = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_total_counts.txt', unpack=True)


first_day = juliandate[0]
last_day = juliandate[-1]

counts = net.tolist()
date = []
for i in range(len(juliandate)):
	date.append(juliandate[i]-first_day)

x = array("d",date)
y = array("d",counts)
c1 = TCanvas('c1', 'Peak Information',600,1800)
c1.Draw()
c1.cd()
p11 = TPad('p11','p',0.0,0.5,1,1)
p11.SetGrid()
p11.Draw()
p12 = TPad('p12','p',0.0,0.0,1,0.5)
p12.Draw()
p11.cd()

real_gr = TGraph(len(x), x, y)
real_gr.Draw('ap')
real_gr.SetMarkerStyle(2)
fsignal = TF1("FSignal","expo", first_day -first_day,last_day-first_day)
fsignal.SetLineColor(2)
real_gr.Fit(fsignal,'R')

normalization = []

try:
	os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_data_detrended_norm2.txt')
except OSError:
	pass
outfile = open('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_data_detrended_norm2.txt','a')

print(fsignal.GetParameter(0))
print(fsignal.GetParameter(1))

p12.cd()

juliandate2 , zeroes , er117, c133, er133, net, ernet, exp, realcount, livetime, fwhm, excesscounts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/real_data_livetime.txt', unpack=True)
first_day2 = juliandate2[0]
last_day2 = juliandate2[-1]
date1 = []
for i in range(len(juliandate2)):
	date1.append(juliandate2[i]-first_day2)

x2 = array("d",date1)
y2 = array("d",net)

real_gr2 = TGraph(len(x2), x2, y2)
real_gr2.Draw('ap')
real_gr2.SetMarkerStyle(2)
fsignal2 = TF1("FSignal","expo",first_day2 - first_day2,last_day2 - first_day2)
fsignal2.SetLineColor(2)
real_gr2.Fit(fsignal2,'R')
B = fsignal2.GetParameter(1)
A = math.exp(fsignal2.GetParameter(0))
print(A)
print(B)
p12.Clear()

fit_dev = TH1F('Residual','Residual',100,-5,5)
fit_dev.Draw('al')
for k in range(len(date1)):
	#print(fit_counts(fsignal2.GetParameter(1),math.exp(fsignal2.GetParameter(0)),date1[k]))

	diff = fit_counts(fsignal2.GetParameter(1),math.exp(fsignal2.GetParameter(0)),date1[k]) - y2[k]
	fit_dev.Fill(diff,1.0/len(juliandate2))
fit_dev.GetYaxis().SetTitleOffset(10.4)
fit_dev.GetXaxis().SetTitleOffset(10.4)
fit_dev.SetLineColor(1)
fit_dev.GetYaxis().SetTitle("Probability")
fit_dev.GetXaxis().SetTitle("Deviation from mean [cps]")

fit_dev.Fit("gaus",'','', -3, 3 )
p12.Update()
'''
for i in range(len(date)):
	outfile.write(str(date[i]+first_day) + "\t" +str(counts[i]/exp_counts(date[i]))+ "\t" +str(err[i]/exp_counts(date[i]))+"\n" )
'''
raw_input("done")