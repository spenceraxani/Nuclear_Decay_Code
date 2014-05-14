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
from ROOT import gROOT, TCanvas, TF1, TGraph, TLine, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt

def drange(start, stop, step):
	r = start
	while r < stop:
		yield r
		r += step
LOMB = 1

clomb = TCanvas('clomb', 'Lomb Scargle Discrete Transform',800,800)
clomb.Draw()
clomb.cd()
p1lomb = TPad('p1lomb','pl',0.05,0.05,0.5,.5)
p1lomb.Draw()
p2lomb = TPad('p1lomb','p2',0.05,0.5,0.95,.98)
p2lomb.Draw()
p3lomb = TPad('p3lomb','p3',0.55,0.05,0.95,.5)
p3lomb.Draw()
p1lomb.cd()

time , counts, errs= numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/hist_mean_fileout.txt', unpack=True)
time_list = []
sample_x = []
mean = np.mean(counts)
for i in range(len(time)):
	time_list.append(time[i]-time[0])
	sample_x.append((counts[i]-mean)*1000)
	#sample_x.append(np.random.normal(0, 2))# + 0.1*np.sin(6.28/30 * time[i]) ) #+ np.sin(2 * 3.1415 / 30 * i))
x = array("d",time_list)
time_series = np.array(x, np.float64)
y = array("d",sample_x)
measurement = np.array(y, np.float64)
f = np.linspace(0.01, 1, 100)
norm = time_series.shape[0]
lombs = sp.signal.lombscargle(time_series , measurement,f)

sample_mean = np.mean(sample_x)*1.0
variance = 0.0
for k in range(len(sample_x)):
	variance += (sample_x[k] - sample_mean)**2
sample_deviation = np.sqrt((1.0/len(sample_x))*variance)
print(sample_deviation)
print(sample_mean)

time_test = []
sample_test = []
for i in range(len(time)):
	time_test.append(time[i]-time[0])
	sample_test.append(np.random.normal(0, sample_deviation))


lombs_gr = TGraph(len(f), f, np.sqrt(4*lombs/norm))
lombs_gr.GetXaxis().SetTitle("Angular Frequency")
lombs_gr.GetYaxis().SetTitle("Some Sort of Probability")
lombs_gr.SetTitle("Lomb Scargle of Averaged data points")
lombs_gr.GetXaxis().SetLimits(0,1);      
lombs_gr.GetHistogram().SetMinimum(0)
lombs_gr.Draw('al')
p1lomb.Update()

N = len(sample_x)*1.0
sumh = 0 
sigma2 = 0 
for j in range(len(sample_x)-1):
	sumh += sample_x[j]
hbar = (1/N)*sumh
hhbar2 = 0
for k in range(len(sample_x)-1):
	hhbar2 += (sample_x[k] - hbar)**2
sigma2 = (1/(N-1))*hhbar2

ang_freq_domain = []
periodogram = []
sum_sin_2wtau = 0
sum_cos_2wtau = 0
sum_x_cos_wttau = 0
sum_cos2_wttau = 0
sum_x_sin_wttau = 0
sum_sin2_wttau = 0

i0=drange(0.01, 1.00, 0.1)
for w in i0:
	print(w)
	for t in time_list:
		sum_sin_2wtau += np.sin(2.0*w*t)
		sum_cos_2wtau += np.cos(2.0*w*t)
	tau = np.arctan(sum_sin_2wtau/sum_cos_2wtau)/(2.0*w)

	for t in range(len(time_list)):
		sum_x_cos_wttau +=  (sample_x[t]-hbar)*np.cos(w*(time_list[t]-tau))
		sum_x_sin_wttau +=  (sample_x[t]-hbar)*np.sin(w*(time_list[t]-tau))
		sum_sin2_wttau += np.sin(w*(time_list[t]-tau))**2
		sum_cos2_wttau += np.cos(w*(time_list[t]-tau))**2
	P = (0.5/sigma2)*(sum_x_cos_wttau**2)/sum_cos2_wttau + (0.5/sigma2)*(sum_x_sin_wttau**2)/sum_sin2_wttau
	periodogram.append(P)
	ang_freq_domain.append(w)
	sum_sin_2wtau = 0
	sum_cos_2wtau = 0
	sum_x_cos_wttau = 0
	sum_cos2_wttau = 0
	sum_x_sin_wttau = 0
	sum_sin2_wttau = 0
p2lomb.cd()

x1 = array("d",ang_freq_domain)
y1 = array("d",periodogram)
sample_gr = TGraph(len(x1), x1, y1)
sample_gr.GetXaxis().SetLimits(0,1);
sample_gr.Draw('al')
p2lomb.Update()

p3lomb.cd()

x2 = array("d",time_list)
y2 = array("d",sample_x)
sample1_gr = TGraph(len(x2), x2, y2)
sample1_gr.GetXaxis().SetLimits(0,1);
sample1_gr.Draw('al')
sample1_gr.GetXaxis().SetTitle("time")
sample1_gr.GetYaxis().SetTitle("normalized counts")
sample1_gr.SetTitle("Lomb Scargle of Averaged data points")
p3lomb.Update()

raw_input("done")

'''

if LOMB == True:
	try:
		os.remove(lomb_fileout)
	except OSError:
		pass
	time , counts, errs= numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/away_from_mean.txt', unpack=True)
	time_series = array("d",time)
	mean = np.mean(counts)
	time_list = []
	sample_x = []
	for i in range(len(time_series)):
		time_list.append(time_series[i]-time_series[0])
		sample_x.append((counts[i]-mean)*1000)


	#for i in time_list:
		#sample_x.append(np.random.normal(0, 2) + 4.1*np.sin(6.28/30 * i) ) #+ np.sin(2 * 3.1415 / 30 * i))
	

	x = np.array(time_list, np.float64)

	N = len(sample_x)*1.0
	sumh = 0 
	sigma2 = 0 
	for j in range(len(sample_x)-1):
		sumh += sample_x[j]
	hbar = (1/N)*sumh
	hhbar2 = 0
	for k in range(len(sample_x)-1):
		hhbar2 += (sample_x[k] - hbar)**2
	sigma2 = (1/(N-1))*hhbar2

	ang_freq_domain = []
	periodogram = []
	ang_freq = 1
	sum_sin_2wtau = 0
	sum_cos_2wtau = 0
	sum_x_cos_wttau = 0
	sum_cos2_wttau = 0
	sum_x_sin_wttau = 0
	sum_sin2_wttau = 0

	i0=drange(0.0, 1.0, 0.005)
	wnum = 0.0
	for i in i0:
		wnum += 1.0

	po=0.01
	z = - np.log(1-np.power(1-po,1/wnum))

	print("the number of frequencies: " +str(wnum))
	print('the z value is :' + str(z))
	print("the variance is: " +str(sigma2))
	print("the number of points is N: " +str(N))
	print("the mean is: " +str(hbar)+ ' ' +str(np.mean(sample_x)))

	i0=drange(0.01, 1.01, 0.005)
	for w in i0:
		print(w)
		for t in time_list:
			sum_sin_2wtau += np.sin(2.0*w*t)
			sum_cos_2wtau += np.cos(2.0*w*t)
		tau = np.arctan(sum_sin_2wtau/sum_cos_2wtau)/(2.0*w)

		for t in range(len(time_list)):
			sum_x_cos_wttau +=  (sample_x[t]-hbar)*np.cos(w*(time_list[t]-tau))
			sum_x_sin_wttau +=  (sample_x[t]-hbar)*np.sin(w*(time_list[t]-tau))
			sum_sin2_wttau += np.sin(w*(time_list[t]-tau))**2
			sum_cos2_wttau += np.cos(w*(time_list[t]-tau))**2
		P = (0.5/sigma2)*(sum_x_cos_wttau**2)/sum_cos2_wttau + (0.5/sigma2)*(sum_x_sin_wttau**2)/sum_sin2_wttau
		periodogram.append(P)
		ang_freq_domain.append(w)
		sum_sin_2wtau = 0
		sum_cos_2wtau = 0
		sum_x_cos_wttau = 0
		sum_cos2_wttau = 0
		sum_x_sin_wttau = 0
		sum_sin2_wttau = 0

	for i in range(len(periodogram)):
		lomboutput = open(lomb_fileout,'a')
		lomboutput.write(str(ang_freq_domain[i])+ "\t" + str(periodogram[i])  + "\n" )
		'''
