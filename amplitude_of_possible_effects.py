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

file1 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_data_residual_detrended.txt' 
file2 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_xray_short.txt' 
date1, value1 = numpy.loadtxt(file1, unpack=True)
date2, value2 = numpy.loadtxt(file2, unpack=True)
date1=date1.tolist()
date2=date2.tolist()
value1=value1.tolist()
value2=value2.tolist()

new_values1 = []
#amplitude = 0.0000000168 # cps/m2-s-sr Good for E 0.8 MeV
#amplitude = 0.00000067 # cps/m2-s-sr Good for E = 2 MeV
#amplitude =  0.00000063 cps/m2-s-sr #Good for P = 1 MeV
#amplitude = 0.0000034 cps/m2-s-sr #Good for P = 10 MeV
#amplitude = 0.00099 cps/m2-s-sr # Good for P = 100 MeV
# amplitude = 16200 #cps - W/s good for Xlong
amplitude = 80290 # good for Xshort
conversion = 1.0

for i in range(len(date1)):
	try:
		new_values1.append(value1[i]+amplitude*conversion*value2[i])
	except IndexError:
		pass

x = array("d",value2)
y = array("d",new_values1)

print(pearsonr(x, y))



