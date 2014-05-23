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
from ROOT import gROOT, TCanvas, TF1, TGraph, TGaxis, TH1, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum, TLine
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt
from scipy.stats import *
from scipy.stats.stats import pearsonr
import random

juliandate , net  = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_data_residual_detrended.txt', unpack=True)

first_day = juliandate[0]
last_day = juliandate[-1]

counts = net.tolist()
date = juliandate.tolist()

x = array("d",date)
y = array("d",counts)


date1 = first_day
date2 = 56601.6
counter = 0
for k in range(len(date)):
	if date[k]> first_day and date[k]<date2:
		counter+=1
total_counts = 0

for i in range(len(date)):
	if date[i]>date1 and date[i]<date2:
		total_counts += counts[i]

n_points = 200
c = TCanvas('canvas','can1',600,400)
c.Draw()
c.cd()
p1 = TPad('MYDATA','MYDATA',0.0,0.0,1,1)
p1.SetGrid()
p1.Draw()
p1.cd()
gr_2 = TH1F('h1','',100, -2, 2)
last_residual_count = 0
summation = []
maximum = 0
increment = (date2 - date1)/(n_points)
a = int(math.floor((last_day-first_day)/increment))
print("The size of the final interval is: " +str(increment))
for i in range(a):
	if i == n_points:
		for j in range(len(date)):
			if date[j] >= (increment * i + date1) and date[j] <= (increment * (i + 1.0) +date1):
				summation.append(counts[j])
		start_date = (increment * i + date1)
		last_residual_count = np.mean(summation)
		gr_2.Fill(np.mean(summation))
		summation = []
	else:
		for j in range(len(date)):
			if date[j] >= (increment * i + date1) and date[j] <= (increment * (i + 1.0) +date1):
				summation.append(counts[j])
		gr_2.Fill(np.mean(summation))
		summation = []
gr_2.Draw()
print("double check the start date: " + str(start_date))
print("The residual counts in the increment before date2 is: " + str(last_residual_count))
print("The probability of observing this result is: " + str(gr_2.Integral(0,gr_2.GetXaxis().FindBin(last_residual_count))/a))
gr_2.Fit("gaus",'','', -3, 3) 

l1 = TLine(last_residual_count,0,last_residual_count,1000)
l1.SetLineStyle(2)
l1.SetLineColor(40)
l1.SetLineWidth(3)
l1.Draw()
p1.Modified()
p1.Update()
gr_2.Clear()

raw_input("done")