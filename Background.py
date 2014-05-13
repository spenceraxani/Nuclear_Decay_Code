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

#################
#Simple graph for background
#################

c = TCanvas('canvas','can',600,600)
c.Draw()
c.cd()
p = TPad('p','pad',0.0,0.0,1,1)
p.SetGrid()
p.Draw()

file_name = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/BACKGROUND_JULY_29_2013.Spe'
#file_name = '/Users/spenceraxani/Documents/Nuclear_Decay/data/datapack/data/2013_Aug_17_L2-008_germanium/Co60_1__800LIVE_000002.Spe'
h1 = TH1F('Background','BACKGROUND_JULY_29_2013',16383,0,16383)
h1.SetTitle("L2-008 Germanium Detector ^{60}CO Spectrum;Energy [keV]; Count Rate [cps]")
h1.Draw("C")
bin = 0
for line in open(file_name):
	if bin>12 and bin < (16383):
		line = int(line)
		h1.SetBinContent(bin, line)
	bin += 1

c.Modified()
c.Update()
print(h1.Integral(12280,12400))

raw_input("done")
