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
from ROOT import gROOT, TCanvas, TF1, TGraph, TLine, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt

cscar = TCanvas('cscar','scar',800,800)
cscar.Draw()
cscar.cd()
pscar1 = TPad('pscar1','p',0.05,0.59,0.95,.98)
pscar1.SetGrid()
pscar1.Draw()

pscar1.cd()
pscar1.SetLogy()
h1 = TH1F('Scargle Amplitude','scargy',80,0,40)
h1.SetTitle("Scargle Amplitude;Amplitude; Counts")
h1.Draw()
spectrum = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/MC_data.txt', unpack=True)
for i in range(len(spectrum)):
	h1.Fill(spectrum[i],1)


pscar1.Modified();
pscar1.Update()
raw_input("done")
