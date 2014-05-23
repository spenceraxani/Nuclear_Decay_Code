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

file1 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/XRAY_5m_Long1.txt' #comparing file
date1, value1 = numpy.loadtxt(file1, unpack=True)

maximum = 0.0
max_time = 0.0
for i in range(len(value1)):
	if value1[i] > 0.0001:
		maximum = value1[i]
		max_time = date1[i]
		print[maximum]
		print[max_time]


