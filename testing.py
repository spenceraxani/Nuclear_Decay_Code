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
##############################################
##
##############################################
def drange(start, stop, step):
	r = start
	while r < stop:
		yield r
		r += step
LOMB = 1

clomb = TCanvas('clomb', 'Lomb Scargle Discrete Transform',800,800)
clomb.Draw()
clomb.cd()
p1lomb = TPad('p1lomb','pl',0.05,0.05,0.95,.5)
p1lomb.Draw()
p2lomb = TPad('p1lomb','pl',0.05,0.5,0.95,.98)
p2lomb.Draw()
p1lomb.cd()
lomb_fileout ='/Users/spenceraxani/Documents/Nuclear_Decay/Data/lomb_fileout_test.txt'

if LOMB == True:
	try:
		os.remove(lomb_fileout)
	except OSError:
		pass
	time , counts, errs= numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/away_from_mean.txt', unpack=True)
	time_series = array("d",time)
	mean = np.mean(counts)
	time_list = []
	sample_x = []

	for i in range(len(time_series)):
		time_list.append(time_series[i]-time_series[0])
		#sample_x.append((counts[i]-mean)*1000)
		sample_x.append(np.random.normal(0, 2) + 0.1*np.sin(6.28/30 * time_series[i]) ) #+ np.sin(2 * 3.1415 / 30 * i))

	x = np.array(time_list, np.float64)
	measurement = array("d",counts)
	y = np.array(measurement, np.float64)
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

	i0=drange(0.005, 2.005, 0.005)
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

	i0=drange(0.01, 1.00, 0.01)
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


p1sigma = 1-0.68269
p2sigma = 1-0.9545
p3sigma = 1-0.9973
p5sigma = 1-0.999999426
z1 = - np.log(1-np.power(1-p1sigma,1/200.0))
z2 = - np.log(1-np.power(1-p2sigma,1/200.0))
z3 = - np.log(1-np.power(1-p3sigma,1/200.0))
z5 = - np.log(1-np.power(1-p5sigma,1/200.0))
print("the most obvious period is: " +str(2*3.14159265/0.164147))
print("the second obvious period is: " +str(2*3.14159265/0.0527322))
print("the third obvious period is: " +str(2*3.14159265/0.221163))
p1 = TLine(0,z1,1,z1);
p1.SetLineStyle(2);
p1.SetLineColor(40)
p1.SetLineWidth(3)
p2 = TLine(0,z2,1,z2);
p2.SetLineStyle(2);
p3 = TLine(0,z3,1,z3);
p3.SetLineStyle(2);
p3.SetLineColor(30)
p3.SetLineWidth(3)
p5 = TLine(0,z5,1,z5);
p5.SetLineStyle(2);
p5.SetLineColor(45)
p5.SetLineWidth(3)

freq, amp  = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/lomb_fileout_test.txt', unpack=True)
time , counts, errs= numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/away_from_mean.txt', unpack=True)
x1 = array("d",freq)
y1 = array("d",amp)

x = array("d",time_list)
y = array("d",sample_x)

lombs_gr = TGraph(len(x1), x1, y1)
lombs_gr.GetXaxis().SetTitle("Angular Frequency")
lombs_gr.GetYaxis().SetTitle("Some Sort of Probability")
lombs_gr.SetTitle("Lomb Scargle of Averaged data points")
lombs_gr.GetXaxis().SetLimits(0,1);
lombs_gr.GetHistogram().SetMaximum(40)          
lombs_gr.GetHistogram().SetMinimum(0)

lombs_gr.Draw('al')
p1.Draw('Same')
p3.Draw('Same')
p5.Draw('Same')
p1lomb.Update()

p2lomb.cd()
measurement = array("d",sample_x)
y = np.array(measurement, np.float64)
sample_gr = TGraph(len(x), x, y)
sample_gr.GetXaxis().SetTitle("Angular Frequency")
sample_gr.GetYaxis().SetTitle("Some Sort of Probability")
sample_gr.SetTitle("Lomb Scargle of Averaged data points")
sample_gr.Draw('al')
p2lomb.Update()
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