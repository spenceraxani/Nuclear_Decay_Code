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

residual_fileout = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_residual.txt" #These are the output file names and locations
try:
	os.remove(residual_fileout)
except OSError:
	pass
date1  , counts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_total_counts.txt', unpack=True) #these are the input files
date2  , expected = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_expected.txt', unpack=True) #these are the input files
date = date1.tolist()
net = counts.tolist()
exp = expected.tolist()
 
for i in range(len(date)):
	residual_file = open(residual_fileout,'a')
	residual_file.write(str(date[i]) + "\t" + str(net[i]-exp[i]) + "\n")