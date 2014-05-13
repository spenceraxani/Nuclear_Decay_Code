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

def probabilityOfResult(X,Y,direction=0):
	x = len(X)
	if x != len(Y):
	    raise ValueError("variables not same len: " + str(x) + ", and " + str(len(Y)))
	if x < 6:
	    raise ValueError("must have at least 6 samples, but have " + str(x))
	(corr, prb_2_tail) = stats.pearsonr(X, Y)
	
	if not direction:
	    #print(str(corr)+ " " + str(1 - prb_2_tail))
	    return (corr, prb_2_tail)
	    
	prb_1_tail = prb_2_tail / 2
	if corr * direction > 0:
		#print(str(corr)+ " " + str(1 - prb_1_tail))
		return (corr, prb_1_tail)
	return (corr, 1 - prb_1_tail)

def correlation(x,y):
	#input needs to be the entire path location and formated as a text file of the form x \t y
	#x, y = numpy.loadtxt(input, unpack=True)
	mean_x = 0.0
	mean_y = 0.0
	N = float(len(x))
	sum_x = 0.0
	sum_y = 0.0
	sum_x_squared = 0.0
	sum_y_squared = 0.0
	xy = 0.0
	diff_x = 0.0
	diff_y = 0.0
	diff_x_squared = 0.0
	diff_y_squared = 0.0
	diff_x_diff_y = 0.0


	for j in range(len(x)):
		sum_x = x[j]
		sum_y = y[j]
	mean_x = sum_x / N
	mean_y = sum_y / N
	for i in range(len(x)):
		sum_x_squared += x[i]**2
		sum_y_squared += y[i]**2
		xy += x[i]*y[i]
		diff_x += x[i] - mean_x
		diff_y += y[i] - mean_y
		diff_x_squared += (x[i] - mean_x)**2
		diff_y_squared += (y[i] - mean_y)**2
		diff_x_diff_y += (x[i] - mean_x)*(y[i] - mean_y)
	#Sx = (diff_x_squared/(N-1))**(0.5)
	Sx = (diff_x_squared)**(0.5)
	#Sy = (diff_y_squared/(N-1))**(0.5)
	Sy = (diff_y_squared)**(0.5)
	Sxy = diff_x_diff_y / (N-1)
	Cov = (xy - N*mean_y*mean_x)/N-1





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
	#print("The covariance is: " + str(Sxy))
	#print("The correlation coefficient is: " + str(r)) 
	return(r)

dictionary = {}
counter = 0
for file in glob.glob('/Users/spenceraxani/Documents/499_Thesis/data/datapack/binned/*.txt'):
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
new_dict = {}
number_of_variables = 15 #basically number of files in the folder
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
		print(correlation(x,y))
		print(pearsonr(x, y))#Check with scipy version of pearsons correlation.
		#print(probabilityOfResult(x,y)[0])

raw_input("done")