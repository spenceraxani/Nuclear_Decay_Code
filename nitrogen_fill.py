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
getcontext().prec = 12

try:
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/nitrogen_data.txt')
except OSError:
	pass
temp = open('/Users/spenceraxani/Documents/499_Thesis/data/datapack/nitrogen.txt', 'r')
counter1 = 0
for line in temp:
	counter1 += 1;
	column = line.split('\t')
	columns = [col.strip() for col in column]
	if columns[-1] == "1":
		on = columns[-1]
		clock1 = columns[0]
		(day, month, year) = clock1.split('/')
		year = int(year)
		month = int(month)
		day = int(day)
		clock2 = columns[1]
		(watch, noon) = clock2.split(' ')
		(hours, mins) = watch.split(":")
		if noon == "AM":
			if int(hours) == int(12):
				hour = 0
			else:
				hour = int(hours)
		else:
			if int(hours) == int(12):
				hour = 12
				#print(hour)
			else:
				hour = int(hours) + 12
			#print(hour)
		actual_time  = float((int(hour) * 3600 + int(mins)*60))/86400
		current_time = time.mktime((year, month, day, 0, 0, 0, 0, 0, 0))
		a =time.gmtime(current_time) 
		proper_time= float(time.gmtime(current_time)[7]) + float(actual_time)
		if year == 2013:
			temp_out = open('/Users/spenceraxani/Documents/499_Thesis/data/datapack/nitrogen_data.txt', 'a')
			temp_out.write(str(56293.50 + proper_time) + "\t" + str(1.0007) + "\n")
		if year == 2014:
			temp_out = open('/Users/spenceraxani/Documents/499_Thesis/data/datapack/nitrogen_data.txt', 'a')
			temp_out.write(str(56293.50 + 365+ proper_time) + "\t" + str(1.0007) + "\n")

raw_input("done")

