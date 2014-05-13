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

def drange(start, stop, step):
	r = start
	while r < stop:
		yield r
		r += step

time , counts, errs= numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/away_from_mean.txt', unpack=True)
time_series = array("d",time)
mean = np.mean(counts)
print(mean)
time_list = []
sample_x = []
for i in range(len(time_series)):
	time_list.append(time_series[i]-time_series[0])
	#sample_x.append((counts[i]-mean)*1000)
x = np.array(time_list, np.float64)
measurement = array("d",counts)
y = np.array(measurement, np.float64)

for i in time_list:
	sample_x.append(np.random.normal(0, 2) + 0.7*np.sin(6.28/30 * i)+1) #+ np.sin(2 * 3.1415 / 30 * i))
	#print(sample_x)

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
print("the most obvious period is: " +str(2*3.14159265/0.165))
print("the second obvious period is: " +str(2*3.14159265/0.22))
print("the third obvious period is: " +str(2*3.14159265/0.788))
i0=drange(0.0, 1.0, 0.005)
for w in i0:
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
print("There is a period at: " + str(2*3.141592/0.057))
plt.figure(1)
plt.subplot(211)
plt.plot(ang_freq_domain, periodogram, 'bo-')
plt.axis([0, 1, 0,20])
plt.subplot(212)
plt.plot(time_list, sample_x, 'bo-')
plt.axis([0, 166, -20.999,20.001])
plt.show()















'''



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
print("The angular frequency is " + str(2 * 3.1415 / 30))
for w in i0:
	for t in time_list:
		sum_sin_2wtau += np.sin(2.0*w*t)
		sum_cos_2wtau += np.cos(2.0*w*t)
	tau = np.arctan(sum_sin_2wtau/sum_cos_2wtau)/(2.0*w)

	for t in range(len(time_list)):
		sum_x_cos_wttau +=  sample_x[t]*np.cos(w*(time_list[t]-tau))
		sum_x_sin_wttau +=  sample_x[t]*np.sin(w*(time_list[t]-tau))
		sum_sin2_wttau += np.sin(w*(time_list[t]-tau))**2
		sum_cos2_wttau += np.cos(w*(time_list[t]-tau))**2
	P = 0.5*(sum_x_cos_wttau**2)/sum_cos2_wttau + 0.5*(sum_x_sin_wttau**2)/sum_sin2_wttau
	periodogram.append(P)
	ang_freq_domain.append(w)
	sum_sin_2wtau = 0
	sum_cos_2wtau = 0
	sum_x_cos_wttau = 0
	sum_cos2_wttau = 0
	sum_x_sin_wttau = 0
	sum_sin2_wttau = 0
print("There is a period at: " + str(2*3.141592/0.057))
plt.figure(1)
plt.subplot(211)
plt.plot(ang_freq_domain, periodogram, 'bo-')
plt.axis([0, 1, 0,180])
plt.subplot(212)
plt.plot(time_list, sample_x, 'bo-')
plt.axis([0, 166, 0.999,1.001])
plt.show()

'''












'''



ang_freq_domain = []
periodogram = []

ang_freq = 100
sum_sin_2wtau = 0
sum_cos_2wtau = 0
sum_x_cos_wttau = 0
sum_cos2_wttau = 0
sum_x_sin_wttau = 0
sum_sin2_wttau = 0

for w in range(ang_freq):
	for t in time_series:
		sum_sin_2wtau += np.sin(2.0*w*t)
		sum_cos_2wtau += np.cos(2.0*w*t)
	tau = np.arctan(sum_sin_2wtau/sum_cos_2wtau)/(2.0*w)

	for t in range(len(time_series)):
		sum_x_cos_wttau +=  measurement[t]*np.cos(w*(time_series[t]-tau))
		sum_x_sin_wttau +=  measurement[t]*np.sin(w*(time_series[t]-tau))
		sum_sin2_wttau += np.sin(w*(time_series[t]-tau))**2
		sum_cos2_wttau += np.cos(w*(time_series[t]-tau))**2
	P = 0.5*(sum_x_cos_wttau**2)/sum_cos2_wttau + 0.5*(sum_x_sin_wttau**2)/sum_sin2_wttau
	periodogram.append(P)
	ang_freq_domain.append(w)
	sum_sin_2wtau = 0
	sum_cos_2wtau = 0
	sum_x_cos_wttau = 0
	sum_cos2_wttau = 0
	sum_x_sin_wttau = 0
	sum_sin2_wttau = 0
print(len(periodogram))
print(len(ang_freq_domain))

'''
