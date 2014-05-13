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
	    print(str(corr)+ " " + str(1 - prb_2_tail))
	    return (corr, prb_2_tail)
	    
	prb_1_tail = prb_2_tail / 2
	if corr * direction > 0:
		print(str(corr)+ " " + str(1 - prb_1_tail))
		return (corr, prb_1_tail)
	return (corr, 1 - prb_1_tail)

def correlation(input):
	#input needs to be the entire path location and formated as a text file of the form x \t y
	x, y = numpy.loadtxt(input, unpack=True)
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
		diff_x_squared = (x[i] - mean_x)**2
		diff_y_squared = (y[i] - mean_y)**2
		diff_x_diff_y = (x[i] - mean_x)*(y[i] - mean_y)
	Sx = (diff_x_squared/(N-1))**(0.5)
	Sy = (diff_y_squared/(N-1))**(0.5)
	Sxy = diff_x_diff_y / (N-1)
	Cov = (xy - N*mean_y*mean_x)/N-1

	r = (xy - N*mean_x*mean_y)/((sum_x_squared - sum_x**2/N)*(sum_y_squared - sum_y**2/N))**(0.5)
	print("The covariance is: " + str(Sxy))
	print("The correlation coefficient is: " + str(r)) 
	return(r)


COR = 1
if COR == True:
	for file in glob.glob('/Users/spenceraxani/Documents/499_Thesis/data/datapack/binned/*.txt'):
		fh = open(file)
		lines = fh.readlines()
		print(lines)






'''
	try:
		os.remove('/Users/spenceraxani/Desktop/cor_file.txt')	
	except OSError:
		pass
	for i in range(1000):
		deadoutput = open('/Users/spenceraxani/Desktop/cor_file.txt','a')
		deadoutput.write(str(i) + "\t" + str(i*10*numpy.sin(i)+i) + "\n" )










x , y = numpy.loadtxt('/Users/spenceraxani/Desktop/cor_file.txt', unpack=True)
x = array("d",x)
y = array("d",y)

correlation('/Users/spenceraxani/Desktop/cor_file.txt')
probabilityOfResult(x,y,0)

cfwhm = TCanvas('c1cfwhm', 'FWHM and Baseline',800,800)
cfwhm.Draw()
cfwhm.cd()
p1fwhm = TPad('p1fwhm','p',0.05,0.05,0.95,.95)
p1fwhm.SetGrid()
p1fwhm.Draw()

fwhm_graph = TGraph(len(x), x, y)
fwhm_graph.Draw("Ap")
'''
raw_input("done")