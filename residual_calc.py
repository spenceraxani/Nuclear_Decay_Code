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

def exp_counts(t):
	counts_exp = 1027.5*numpy.exp(-(t-56546.9071712)*0.000360576)
	return(counts_exp)

binned_residual_dtemp = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_detrended_residual.txt"
binned_residual_norm = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_detrended_normalization.txt"
residual_fileout = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/residual_dtemp.txt"
residual_org = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/residual.txt"
try:
	os.remove(residual_fileout)
	os.remove(binned_residual_dtemp)
	os.remove(binned_residual_norm)
except OSError:
	pass

juliandate , net = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_data_detrended.txt', unpack=True)
x = array("d",juliandate)
y = array("d",net)

date_1 , net_1 = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_data_detrended.txt', unpack=True)
x = array("d",date_1)
y = array("d",net_1)

date = juliandate.tolist()
counts = net.tolist()

first_day = date[0]
last_day = date[-1]

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
real_gr.Draw('al')
fsignal = TF1("FSignal","expo", first_day ,last_day)
fsignal.SetLineColor(2)
real_gr.Fit(fsignal,'R')


time , temperature = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_temperature.txt', unpack=True)
time2 , pressure = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_pressure.txt', unpack=True)
mean_pressure = np.mean(pressure)
mean_temperature = np.mean(temperature)

detrend_temp = []
detrend_pressure = []
detrend_all = []
try:
	os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended_with_temp.txt')
	os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended_with_pressure.txt')
	os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_detrended.txt')

except OSError:
	pass
outfile = open('/Users/spenceraxani/Documents/Nuclear_Decay/Data/detrended_with_temp.txt','a')
pressure_detrended_outfile = open('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended_with_pressure.txt','a')
temperature_detrended_outfile = open('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended_with_temp.txt','a')
dis_to_temperature = 22.2*0.142000
eff_to_dis = 0.000051#508857
for i in range(len(time)):
	detrend_temp.append(net[i]*(1.0  + (temperature[i]-mean_temperature)*dis_to_temperature*eff_to_dis))
	temperature_detrended_outfile.write(str(time[i]) + " \t" +str(net[i]*(1.0  + (temperature[i]-mean_temperature)*dis_to_temperature*eff_to_dis))+"\n")
	outfile.write(str(date[i]) + " \t" +str(counts[i]*(1.0  - (pressure[i]-mean_pressure)*dis_to_pressure*eff_to_dis + (temperature[i]-mean_temperature)*dis_to_temperature*eff_to_dis))+"\n" )
	pressure_detrended_outfile.write(str(date[i]) + " \t" +str(counts[i]*(1.0  - (pressure[i]-mean_pressure)*dis_to_pressure*eff_to_dis)+"\n" )

print(len(time))
print(len(date))

for i in range(len(counts)):
	residual_file = open(residual_fileout,'a')
	#print(exp_counts(date[i]))
	residual_file.write(str(date[i]) + "\t" + str(detrend_temp[i]-exp_counts(date[i])) + "\n")
raw_input("done")