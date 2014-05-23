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
from ROOT import gROOT, TCanvas, TF1, TGraph, TGaxis, TH1, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum, TLine
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt
from scipy.stats import *
from scipy.stats.stats import pearsonr
import random

Red = TColor(2000,1,0,0)
Redish = TColor(2001,1,.4,0)
Redishish = TColor(2002,1,.8,0)	
Yellow = TColor(2003,0.5,1,0)	
Yellowish = TColor(2004,0,1,1)		
Orange = TColor(2005,0,.5,1)

first_day = 56506.8942593
last_day = 56673.0953472

fileout_xray1  	= '/Users/spenceraxani/Documents/Nuclear_Decay/Data/XRAY_5m_Long.txt'
fileout_xray2 	= '/Users/spenceraxani/Documents/Nuclear_Decay/Data/XRAY_5m_Short.txt'

c = TCanvas('canvas','solars',600,500)
c.Draw()
c.cd()
p2 = TPad('p2','p2',0.0,0.0,1,1)
p2.SetGrid()
p2.Draw()
p2.cd()

juliandate , away, erry = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_data_detrended_norm_with_errors.txt', unpack=True)
x5 = array("d",juliandate)
y5 = array("d",away)

errx = []
for i in range(len(erry)):
	errx.append(0)

ey = array("d",erry)
ex = array("d",errx)
detrend_gr = TGraphErrors(len(x5), x5, y5,ex,ey)
detrend_gr.SetMarkerStyle(4)
detrend_gr.SetMarkerSize(1)
#####################
#This code only works in sections, you need to comment out 2 of the 3 graphs:
#####################

#Graph 1
xl = []
xld = []
xs = []
xsd = []

solar_date_xray_l  , solar_flux_xray_l = numpy.loadtxt(fileout_xray1, unpack=True)
solar_date_xray_s 	, solar_flux_xray_s = numpy.loadtxt(fileout_xray2, unpack=True)

for i in range(len(solar_flux_xray_l)):
	xl.append(solar_flux_xray_l[i]*20+0.994)
	xld.append(solar_date_xray_l[i])
	xs.append(solar_flux_xray_s[i]*20+0.994)
	xsd.append(solar_date_xray_s[i])
	#print(solar_flux_xray_s[i])
x_xrayl2 = array("d",xl)
x_xrays2 = array("d",xs)
x_xrayld2 = array("d",xld)
x_xraysd2 = array("d",xsd)
xldg1 = TGraphErrors(len(x_xrayld2), x_xrayld2, x_xrayl2)
xldg1.SetMarkerColor(2005)
xldg1.SetLineColor(2005)
xsdg1 = TGraph(len(x_xraysd2), x_xraysd2, x_xrays2)
xsdg1.SetMarkerColor(2004)
xsdg1.SetLineColor(2004)

mg_norm1 = TMultiGraph()
#mg_norm1.SetTitle(";Modified Julian Date MST [days]; Normalized Count Rate [cps]")
mg_norm1.SetTitle(";Modified Julian Date MST [days]; ")

mg_norm1.Add(xldg1)
mg_norm1.Add(xsdg1)
mg_norm1.Add(detrend_gr)
mg_norm1.Draw('ap')

#mg_norm1.GetYaxis().SetTitleOffset(1.5)
mg_norm1.GetYaxis().SetTitleOffset(1.)
a = mg_norm1.GetYaxis().GetLabelFont()
print(a)
mg_norm1.GetHistogram().SetMaximum(1.003)      
mg_norm1.GetHistogram().SetMinimum(0.994)
#mg_norm1.GetXaxis().SetLimits(56587,56597)
#mg_norm1.GetXaxis().SetLimits(56597,56603)
mg_norm1.GetXaxis().SetLimits(56661,56667)

axis1 = TGaxis(56667,0.994,56667,1.003,0,0.000450675,510,"+L")
axis1.SetName("axis1")
axis1.SetLabelColor(1)
axis1.SetTitle("X-Ray Flux [W/m^{2}]")
axis1.SetTitleOffset(1.2)
axis1.SetLabelFont(42)
axis1.SetTitleFont(42)
axis1.Draw()
'''
leg1 = TLegend(0.48, 0.74, 0.89, 0.89)
leg1.SetFillColor(0)
leg1.AddEntry(detrend_gr, "^{60}Co Count Rate", "lp")
leg1.AddEntry(xldg1, "Long Wavelength X-Rays: 0.1 - 0.8 nm", "lp")
leg1.AddEntry(xsdg1, "Short Wavelength X-Rays: 0.05 - 0.4 nm", "lp")
leg1.Draw()
'''
p2.Modified()
p2.Update()

raw_input("done")