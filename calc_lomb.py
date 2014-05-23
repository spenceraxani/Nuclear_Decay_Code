import os
import glob
import math
from math import log10, floor
import time
import sys
import scipy.signal
import scipy.stats as st
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
###########################################
#So, this will first manually calculates a Lomb Scargle periodogram for the first graph, and then
###########################################


def exp_counts(t): #this just calculates the expected count rate given the values of the best fit
	half_life = 1926.14 	#days	
	initial_counts = 1027.57125
	counts_exp = (initial_counts)*numpy.exp(-numpy.log(2)*(t-56546.9071712)/half_life)
	return(counts_exp)

print(exp_counts(56547.9071712))

time , counts= numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/nitrogen_data.txt', unpack=True) #Calc Lomb Scargle from this file, and give confidence intervals
time_series = array("d",time)
mean = np.mean(counts)
time_list = []
sample_x = []
for i in range(len(time_series)):
	time_list.append(time_series[i]-time_series[0])
	#print(time_series[i])
	sample_x.append((counts[i]-exp_counts(time_series[i]))*1000.0)
print(sample_x)
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

i0 = np.linspace(0.0172142, 2* np.pi, 600)
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

i0 = np.linspace(0.0172142, 2* np.pi, 600)
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
x1 = array("d",ang_freq_domain)
x2 = np.array(x1, np.float64)
y1 = array("d",periodogram)
y2 = np.array(y1, np.float64)

c1 = TCanvas('cscar2','scar1',600,1800)
c1.Draw()
c1.cd()
c1p1 = TPad('MYDATA','MYDATA',0.0,0.0,1,.5)
c1p1.SetGrid()
c1p1.Draw()
c1p2 = TPad('SIMDATA','SIMDATA',0.0,0.5,1,1)
c1p2.SetGrid()
c1p2.Draw()
c1p2.cd()
gr_2 = TGraph(len(x2), x2, y2)
gr_2.GetHistogram().SetMaximum(30)          
gr_2.GetHistogram().SetMinimum(0)
gr_2.GetXaxis().SetLimits(0,2* np.pi)
gr_2.SetTitle("Simulated Data; Frequency [day^{-1}]; Periodogram Amplitude [Arbitrary Units]")
gr_2.SetLineColor(1)
gr_2.SetLineWidth(1)
gr_2.Draw('al')
z1 = 0.317311
z3 = 0.0027
z5 = float(0.0000006004)

v5 = - np.log(1-np.power(1.0-z5,1.0/600.0))
v3 = - np.log(1-np.power(1.0-z3,1.0/600.0))
v1 = - np.log(1-np.power(1.0-z1,1.0/600.0))

p1 = TLine(0,v1,2* np.pi,v1)
p1.SetLineStyle(2)
p1.SetLineColor(40)
p1.SetLineWidth(3)
p3 = TLine(0,v3,2* np.pi,v3)
p3.SetLineStyle(2);
p3.SetLineColor(30)
p3.SetLineWidth(3)
p5 = TLine(0,v5,2* np.pi,v5);
p5.SetLineStyle(2);
p5.SetLineColor(45)
p5.SetLineWidth(3)

legMC = TLegend(0.44,0.71,0.89,0.89)
legMC.SetFillColor(0)
legMC.AddEntry(gr_2,"Lomb Scargle Periodogram","l")
legMC.AddEntry(p1,"1 #sigma","l")
legMC.AddEntry(p3,"3 #sigma","l")
legMC.AddEntry(p5,"5 #sigma","l")
legMC.Draw()

p1.Draw("same")
p3.Draw("same")
p5.Draw("same")

c1p2.Modified()
c1p2.Update()

c1p1.cd()
'''
time , counts = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/binned_total_counts.txt', unpack=True)
time_series = array("d",time)
mean_counts = np.mean(counts)
time_list = []
sample_x = []
for i in range(len(time_series)):
	time_list.append(time_series[i]-time_series[0])
	sample_x.append((counts[i]-mean_counts)*1000.0)
print(sample_x)
x = np.array(time_list, np.float64)
measurement = array("d",counts)
y = np.array(measurement, np.float64)
'''
time , counts= numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_total_counts.txt', unpack=True) #use this file to calculate the second lombs scargle
time_series = array("d",time)
mean = np.mean(counts)
time_list = []
sample_x = []
for i in range(len(time_series)):
	time_list.append(time_series[i]-time_series[0])
	#print(time_series[i])
	sample_x.append((counts[i]-exp_counts(time_series[i]))*1000.0)
print(sample_x)
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

i0 = np.linspace(0.0172142, 2* np.pi, 600)
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

i0 = np.linspace(0.0172142, 2* np.pi, 600)
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
x3 = array("d",ang_freq_domain)
x4 = np.array(x3, np.float64)
y3 = array("d",periodogram)
y4 = np.array(y3, np.float64)
gr_3 = TGraph(len(x4), x4, y4)
gr_3.GetHistogram().SetMaximum(30)          
gr_3.GetHistogram().SetMinimum(0)
gr_3.GetXaxis().SetLimits(0,2*np.pi);
gr_3.SetTitle("Normal Data; Frequency [day^{-1}]; Periodogram Amplitude [Arbitrary Units]")
gr_3.SetLineColor(1)
gr_3.SetLineWidth(1)
gr_3.Draw('al')
legMC.Draw()
#gr_2.Draw("same")
p1.Draw("same")
p3.Draw("same")
p5.Draw("same")

c1p1.Modified()
c1p1.Update()
raw_input("done")

