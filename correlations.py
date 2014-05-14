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

cave = TCanvas('cave', 'ave',600,1800)
cave.Draw()
cave.cd()
p1ave = TPad('p1ave1','p',0.0,0.5,1,1)
p1ave.SetGrid()
p1ave.Draw()
p2ave = TPad('p2ave2','p',0.0,0.0,1,0.5)
p2ave.SetGrid()
p2ave.Draw()
p1ave.cd()

juliandate1 , away1  = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_counts.txt', unpack=True)
x1 = array("d",juliandate1)
y1 = array("d",away1)

juliandate2 , away2  = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_pressure.txt', unpack=True)
x2 = array("d",juliandate2)
y2 = array("d",away2)

plot_away = TGraph(len(y1), y1, y2)
plot_away.Draw("ap")
p1ave.Update()
raw_input("done")