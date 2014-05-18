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
from ROOT import gROOT, TGaxis, TPaveText,TCanvas, TF1, TGraph, TLine, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt
import time
from scipy.stats import *
from scipy.stats.stats import pearsonr
import random


def correlation(x,y):
	#input two sets of arrays. like x = array("d",counts) and y = array("d",electron_energy)	
	sum_diff_x_square = 0
	sum_diff_y_square = 0
	mean_x = np.mean(x)
	mean_y = np.mean(y)
	sum_diff_x_times_diff_y =0
	for i in range(len(x)):
		sum_diff_x_times_diff_y += (x[i]-mean_x)*(y[i]-mean_y)
		sum_diff_x_square += (x[i]-mean_x)**2
		sum_diff_y_square += (y[i]-mean_y)**2
	r = sum_diff_x_times_diff_y/(np.sqrt(sum_diff_x_square)*np.sqrt(sum_diff_y_square))
	return(r)


#now permute test that mother!
trials = 50000 #number of bootstraping data sets to generate. 50000 generatlly
#pick the file to bootstrap.

file1 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_residual.txt' #comparing file
file2 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_temperature.txt' #this one, bootstrap
date1, value1 = numpy.loadtxt(file1, unpack=True)
date2, value2 = numpy.loadtxt(file2, unpack=True)

p_corr_list = []
new_temp_dict = {}
temp_values_dict = {}
another_dict = {}
temp_another_dict = {}
print("Calculating the permutation distribution...")
for i in range(trials):
	#print(i)
	for i in range(len(date1)):
		new_temp_dict[date1[i]] = [value1[i]] #here, column[0] is the date. So element [date] in dictionary is [value(column[1])]
	for j in range(len(date2)):
		if date2[j] in new_temp_dict:
			new_temp_dict[date2[j]].append(value2[j]) 
	for k in new_temp_dict:
		if len(new_temp_dict[k]) == 2:
			temp_values_dict[k]=new_temp_dict[k] # now all the time stamps and data sets line up.
	list1 = []
	list2 = []
	for i in temp_values_dict:
		list1.append(float(temp_values_dict[i][0]))
		list2.append(float(temp_values_dict[i][1]))
	x = array("d",list1)
	y = array("d",sorted(list2, key=lambda *args: random.random())) #this randomly sorts the list, this is the "permutation test"
	p_corr_list.append(correlation(x,y)) #this list now holds all the various values for the pearsons correlation
	#print(correlation(x,y))
	#print(pearsonr(x, y))
	new_temp_dict = {}
	temp_values_dict = {}

#Calc the actual correlation of the data
for i in range(len(date1)):
	another_dict[date1[i]] = [value1[i]] #here, column[0] is the date. So element [date] in dictionary is [value(column[1])]
for j in range(len(date2)):
	if date2[j] in another_dict:
		another_dict[date2[j]].append(value2[j])
for k in another_dict:
	if len(another_dict[k]) == 2: # two because you are looking at only 2 variables
		temp_another_dict[k] = another_dict[k]
print(len(temp_another_dict))
list1 = []
list2 = []
for i in temp_another_dict:
	list1.append(float(temp_another_dict[i][0]))
	list2.append(float(temp_another_dict[i][1]))
x = array("d",list1)
y = array("d",list2)
corr = correlation(x,y)



c = TCanvas('canvas','can1',600,600)
c.Draw()
c.cd()
tp1 = TPad('MYDATA','MYDATA',0.0,0.0,1,1)
tp1.SetGrid()
tp1.Draw()
tp1.cd()

maximum = 1
bins =500
bin_to_value = (bins/2)/maximum
gr_2 = TH1F('h1','myhist',bins, -maximum, maximum)
for i in range(len(p_corr_list)):
	gr_2.Fill(p_corr_list[i],1.0/trials)
gr_2.Draw()
fitsignal = TF1("FSignal","gaus", -maximum, maximum) # this gaussian can then be integrated to get the two tailed p-values. Input these values into correlation_gaussian.np code to integrate.
gr_2.Fit(fitsignal,'R')

norm_gaus = fitsignal.Integral(-1,1)

if corr < 0:
	integrals = 2*fitsignal.Integral(-1,corr)/norm_gaus
	print("The two tailed integral is : "+ str(2*fitsignal.Integral(-1, corr)))
if corr > 0:
	integrals = 2*fitsignal.Integral(corr,1)/norm_gaus
	print("The two tailed integral is : "+ str(2*fitsignal.Integral(corr,1)/norm_gaus))
	#print("The full gaussian integral is : "+ str(2*fitsignal.Integral(-1,1)))
	#print("The integral of the histogram is : "+ str(gr_2.ComputeIntegral()))
p1 = TLine(corr,0,corr,1)
p1.SetLineStyle(2)
p1.SetLineColor(40)
p1.SetLineWidth(3)

legMC = TLegend(0.44,0.71,0.89,0.89)
legMC.SetFillColor(0)
legMC.AddEntry(gr_2,"Permutation Histogram","l")
legMC.AddEntry(p1,"The PCC","l")
#legMC.Draw()
p1.Draw("same")

print("My calculation yeilds a [r, p] = " +str(corr)+ ", "+ str(integrals))
print("Pythons calculatino yeilds a [r, p] = " +str(pearsonr(x, y)))

tp1.Modified()
tp1.Update()
raw_input("done")