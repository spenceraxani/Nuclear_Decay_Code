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

c23proton = TCanvas('c23pro1ton', '1Proton Flux Maximum Activity',600,1800)
c23proton.Draw()
c23proton.cd()
p23proton1 = TPad('p23prot1on1','p21',0.0,0.5,1,1)
p23proton1.SetGrid()
p23proton1.Draw()
p23proton2 = TPad('p23prot1on2','p21',0.0,0.0,1,0.5)
p23proton2.SetGrid()
p23proton2.Draw()
p23proton2.cd()

juliandate , zeroes , er117, c133, er133, net, ernet, exp, realcount, livetime,fwhm, excesscounts = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/real_data.txt', unpack=True)
temp_date = []
for i in range(len(juliandate)):
	temp_date.append(juliandate[i]-juliandate[0])
xreal = array("d",temp_date)
#xreal = array("d",juliandate)
yreal = array("d",net)
real_gr = TGraph(len(xreal), xreal, yreal)

real_gr.GetHistogram().SetMaximum(1050)
real_gr.GetHistogram().SetMinimum(980)

real_gr.GetXaxis().SetTitle("Modified Julian Date MST [Days]")
real_gr.GetYaxis().SetTitle("Count Rate [cps]")
real_gr.SetTitle("Cobalt 60 decay count rate")
#real_gr.GetXaxis().SetLimits(first_day,last_day);
#real_gr.GetYaxis().SetLimits(980,1050);

real_gr.Draw("ap")
fsignal = TF1("FSignal","expo", 0 ,temp_date[-1])
#real_gr.Fit("fsignal","","",first_day ,last_day)
fsignal.SetLineColor(2)
real_gr.Fit(fsignal,'R')
p23proton2.Modified();
p23proton2.Update()
raw_input("done")