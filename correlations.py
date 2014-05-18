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

juliandate1 , away1  = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/residual.txt', unpack=True)
x1 = array("d",juliandate1)
y1 = array("d",away1)

juliandate2 , away2  = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_pressure.txt', unpack=True)
x2 = array("d",juliandate2)
y2 = array("d",away2)

juliandate3 , away3  = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/residual.txt', unpack=True)
x3 = array("d",juliandate3)
y3 = array("d",away3)

juliandate4 , away4  = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_temperature.txt', unpack=True)
x4 = array("d",juliandate4)
y4 = array("d",away4)

plot_away = TGraph(len(y2), y2, y1)
plot_away.SetTitle(";Pressure [kPa]; Residual Count Rate [cps]")
plot_away.GetXaxis().SetLimits(90,97)
plot_away.SetMarkerStyle(25)
plot_away.Fit("pol1",'','', 90, 97) 
pt = TPaveText(90.5,1,93,1.5);
pt.AddText("dc/dP = 0.187 \pm 0.011~kPa^{-1}")
pt.SetFillColor(2)
pt.SetFillStyle(0)

plot_away.Draw("ap")
pt.Draw()
p1ave.Update()

p2ave.cd()

plot_away2 = TGraph(len(y4), y4, y3)
plot_away2.SetTitle(";Temperature [Celcius]; Residual Count Rate [cps]")
plot_away2.GetXaxis().SetLimits(14,17)
plot_away2.SetMarkerStyle(25)
plot_away2.Fit("pol1",'','', 14, 17) 
pt1 = TPaveText(15.5,1,16.75,1.5);
pt1.AddText("dc/T = - 0.100 \pm 0.035~T^{-1}")
pt1.SetFillColor(1)
pt1.SetFillStyle(0)

plot_away2.Draw("ap")
pt1.Draw()

p2ave.Update()

raw_input("done")