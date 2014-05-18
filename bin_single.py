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
BINXRAY = 1
###############################
#Used to bin data sets so that all the time series line up.
#When calculating the correlation coefficients, the data needs all the time series to be lined up
###############################

first_day = 56506.8942593
last_day = 56673.0953472
ave_pressure = 94.0170575183
ave_temperature = 15.0531177202
c7x = TCanvas('c7x','Peakx',600,1800)
c7x.Draw()
c7x.cd()
p70x = TPad('p70x','px',0.0,0.5,1,1)
p70x.SetGrid()
p70x.Draw()
p71x = TPad('p71x','px',0.0,0.0,1,0.5)
p71x.SetGrid()
p71x.Draw()
p71x.cd()

xray_bin = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_single.txt" #These are the output file names and locations
n_points = 2000 #all the others were binned into 2000 bins

if BINXRAY == True:
	try:
		os.remove(xray_bin)
	except OSError:
		pass
	#xray_date , zeroes , er117, c133, er133, net, ernet, xray_long, realcount, livetime, fwhm, excesscounts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/real_data.txt', unpack=True)

	xray_date  , xray_long = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/residual.txt', unpack=True) #these are the input files
	increment = (last_day - first_day)/n_points
	date = xray_date.tolist()
	value = xray_long.tolist()
	xray_bin_out= open(xray_bin,'a')
	x_ray_counts = []
	x_ray_sums = 0.0

	for i in range(n_points):
		print(i)
		for j in range(len(date)):
			if date[j] >= (increment * i + first_day) and date[j] <= (increment * (i + 1.0) +first_day):
				x_ray_counts.append(0.00001)
				#x_ray_counts.append(value[j]*1.0)
				#print(value[j]*1.0)
		if x_ray_counts:
			x_ray_sums = np.average(x_ray_counts)
			xray_bin_out.write(str(((increment * i + first_day)+(increment * (i + 1) +first_day))/2.0)+ '\t' +str(x_ray_sums*1.0)+"\n" )
		else:
			xray_bin_out.write(str(((increment * i + first_day)+(increment * (i + 1) +first_day))/2.0)+ '\t' +str(0.0)+"\n" )#WATCH THIS ONE! If no data, what do you want to put in the bin??

		x_ray_sums = 0.0
		x_ray_counts = []