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
BINXRAY = 0

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

xray_bin = "/Users/spenceraxani/Documents/499_Thesis/data/datapack/binned_humidity.txt"
counts_bin = "/Users/spenceraxani/Documents/499_Thesis/data/datapack/binned_pressure.txt"
n_points = 2000

if BINXRAY == True:
	try:
		os.remove(xray_bin)
		#os.remove(counts_bin)
	except OSError:
		pass
	#count_date , zeroes , er117, c133, er133, net, ernet, exp, realcount, livetime,counts, excesscounts = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/real_data.txt', unpack=True)
	xray_date  , xray_long = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/temperture_file.txt', unpack=True)
	#xray_date , xray_long , peak2 = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/peak_data.txt', unpack=True)
	count_date  , counts = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/pressure_data.txt', unpack=True)
	increment = (last_day - first_day)/n_points
	x_ray_sums = 0
	counts_sums = 0
	x_ray_counts = []
	counts_counts = []

	xray_bin_out= open(xray_bin,'a')
	counts_bin_out= open(counts_bin,'a')

	for i in range(n_points):
		print(i)

		for j in range(len(xray_date)):
			if xray_date[j] >= (increment * i + first_day) and xray_date[j] <= (increment * (i + 1) +first_day):
				x_ray_counts.append(xray_long[j])
		#print(len(x_ray_counts))
		if x_ray_counts:
			#print(str(x_ray_counts) + "\n")
			#print(str(np.average(x_ray_counts)) + "\n")
			x_ray_sums = np.average(x_ray_counts)
			xray_bin_out.write(str(((increment * i + first_day)+(increment * (i + 1) +first_day))/2)+ '\t' +str(x_ray_sums)+"\n" )
		else:
			xray_bin_out.write(str(((increment * i + first_day)+(increment * (i + 1) +first_day))/2)+ '\t' +str(ave_temperature)+"\n" )

		x_ray_sums = 0
		x_ray_counts = []

		for k in range(len(count_date)):
			if count_date[k] >= (increment * i + first_day) and count_date[k] <= (increment * (i + 1) +first_day):
				counts_counts.append(counts[k])
		'''
		if counts_counts:
			#print(str(np.average(counts_counts)) + "\n")
			counts_sums = np.average(counts_counts)
			counts_bin_out.write(str(((increment * i + first_day)+(increment * (i + 1) +first_day))/2)+ '\t' +str(counts_sums)+"\n" )
		else:			
			counts_bin_out.write(str(((increment * i + first_day)+(increment * (i + 1) +first_day))/2)+ '\t' +str(ave_pressure)+"\n" )
		'''
		counts_counts = []
		counts_sums = 0
	
xray_dates , ave_xray = numpy.loadtxt(xray_bin, unpack=True)
xray_hist = TH1F('X-Ray Histrogram','xh',n_points,first_day,last_day)
xray_hist.SetTitle("; Modified Julian Date MST [days];W")
for i in range(len(xray_dates)):
	#print(np.amax(ave_xray))
	if ave_xray[i] != "nan":
		#print(ave_xray[i])
		xray_hist.Fill(xray_dates[i],ave_xray[i])
xray_hist.Draw()
p71x.Modified()
p71x.Update()

p70x.cd()
counts_dates , ave_counts = numpy.loadtxt(counts_bin, unpack=True)
counts_hist = TH1F('Counts Histrogram','xh',n_points,first_day,last_day)
counts_hist.SetTitle("; Modified Julian Date MST [days];Count Rate [cps]")
for j in range(len(counts_dates)):	
	counts_hist.Fill(counts_dates[j],ave_counts[j])
counts_hist.Draw()
p70x.Modified()
p70x.Update()
raw_input("done")