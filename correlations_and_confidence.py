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
from ROOT import gROOT, TCanvas, TF1, TGraph, TH1, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum, TLine
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt
from scipy.stats import *
from scipy.stats.stats import pearsonr
import random

##############################
#Input data sets that are of the same time series (all those that have been binned, found in the binned folder)
###############################

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

dictionary = {}
counter = 0
for file in glob.glob('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/*.txt'): #these files are the binned files to calculate the Pearsons correlations using my calc and the Scipy pearsonsr calc
	print(file)
	fh = open(file)
	lines = fh.readlines()
	if counter == 0: #counter==0 takes the first file in the folder. All the other files (ie. counter =! 0) will then be lined up with the same time stamps
		for element in lines:
			columns = element.split('\t')
			dictionary[columns[0]] = [columns[1]] #here, column[0] is the date. So element [date] in dictionary is [value(column[1])]
			counter += 1
	else:
		for element in lines:
			columns = element.split('\t')
			if columns[0] in dictionary: #
				dictionary[columns[0]].append(columns[1]) #so now the dictionary is dictionary[date]=[counts,deadtime,electron,...]
#Now calculate the correlation or whatever from the elements in the dictionary.
print("Calculating Pearsons Correlations ...")
new_dict = {}
number_of_variables = 16 #basically number of files in the folder
for element in dictionary: #This loop gets rid of missing data from any of the data sets, and only uses the dates at which there is data for every variable.
	#print(len(dictionary[element]))
	if len(dictionary[element]) == number_of_variables:
		new_dict[element] = dictionary[element] #this new_dict now always has length number_of_variables
#print(len(new_dict))
for j in range(number_of_variables): #Now loop through all entries to calculate correlation stuff. like counts vs counts, then counts vs dead_time... 
	for k in range(number_of_variables):
		list1 = []
		list2 = []
		for i in new_dict:
			list1.append(float(new_dict[i][j]))
			list2.append(float(new_dict[i][k]))
		#print(list1)
		x = array("d",list1)
		y = array("d",list2)	
		#print(correlation(x,y))
		#print(pearsonr(x, y))#Check with scipy version of pearsons correlation.
		#print(probabilityOfResult(x,y)[0])

#now bootstrap that mother!
trials = 50 #number of bootstraping data sets to generate. 50000 generatlly
#pick the file to bootstrap.

file1 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_data_detrended_residual.txt' #comparing file
file2 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_proton_100MeV.txt' #this one, bootstrap
date1, value1 = numpy.loadtxt(file1, unpack=True)
date2, value2 = numpy.loadtxt(file2, unpack=True)

p_corr_list = []
new_temp_dict = {}
temp_values_dict = {}
another_dict = {}
temp_another_dict = {}
print("Calculating the bootstrapping distribution ...")
for i in range(trials):
	for i in range(len(date1)):
		new_temp_dict[date1[i]] = [value1[i]] #here, column[0] is the date. So element [date] in dictionary is [value(column[1])]
	for j in range(len(date2)):
		if date2[j] in new_temp_dict:
			new_temp_dict[date2[j]].append(value2[random.randint(0,len(value2)-1)])
			#print(random.randint(0,len(value2)))
	for k in new_temp_dict:
		if len(new_temp_dict[k]) == 2:
			temp_values_dict[k]=new_temp_dict[k]
	list1 = []
	list2 = []
	for i in temp_values_dict:
		list1.append(float(temp_values_dict[i][0]))
		list2.append(float(temp_values_dict[i][1]))
	x = array("d",list1)
	y = array("d",list2)
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
list1 = []
list2 = []
for i in temp_another_dict:
	list1.append(float(temp_another_dict[i][0]))
	list2.append(float(temp_another_dict[i][1]))
x = array("d",list1)
y = array("d",list2)
corr = correlation(x,y)

print(corr)
print(pearsonr(x, y))
zr = 0.5*np.log((1+0.176401467636)/(1-0.176401467636))
za = 1.96
cl = zr -za*(np.sqrt(1/(len(x)-3)))
rl = (np.exp(2*cl)-1)/(np.exp(2*cl)+1)
print(rl)
print(corr*np.sqrt(len(x)-2)/np.sqrt(1-corr**2))

c = TCanvas('canvas','can1',600,1800)
c.Draw()
c.cd()
p1 = TPad('MYDATA','MYDATA',0.0,0.0,1,.5)
p1.SetGrid()
p1.Draw()
p2 = TPad('SIMDATA','SIMDATA',0.0,0.5,1,1)
p2.SetGrid()
p2.Draw()
p2.cd()

maximum = 0.1
bins =100
bin_to_value = (bins/2)/maximum

gr_2 = TH1F('h1','myhist',bins, -maximum, maximum)
for i in range(len(p_corr_list)):
	gr_2.Fill(p_corr_list[i])
gr_2.Draw()

gr_2.Fit("gaus",'','', -maximum, maximum) # this gaussian can then be integrated to get the two tailed p-values. Input these values into correlation_gaussian.np code to integrate.

p1 = TLine(corr,0,corr,1000)
p1.SetLineStyle(2)
p1.SetLineColor(40)
p1.SetLineWidth(3)

legMC = TLegend(0.44,0.71,0.89,0.89)
legMC.SetFillColor(0)
legMC.AddEntry(gr_2,"boostrap graph","l")
legMC.AddEntry(p1,"The Pearsons value","l")
legMC.Draw()

p1.Draw("same")

if corr >= 0:
	#print(gr_2.Integral(int(gr_2.GetXaxis().FindBin(corr)),bins/2))
	print(2*gr_2.Integral(int(gr_2.GetXaxis().FindBin(corr)),bins/2)/trials)
else:
	print(2*gr_2.Integral(-bins/2,int(gr_2.GetXaxis().FindBin(corr)))/trials)

p2.Modified()
p2.Update()
raw_input("done")

#print(p_corr_list)

#Check with scipy version of pearsons correlation.
#print(new_temp_dict[date2[j]])