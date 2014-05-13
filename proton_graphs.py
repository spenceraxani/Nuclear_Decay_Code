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

cproton = TCanvas('cproton', 'Proton',600,1800)
cproton.Draw()
cproton.cd()
Proton_pad1 = TPad('x10','p',0.0,0.5,1,1)
Proton_pad1.SetGrid()
Proton_pad1.Draw()
Proton_pad2 = TPad('x11','p',0.0,0.0,1,0.5)
Proton_pad2.SetGrid()
Proton_pad2.Draw()
Proton_pad2.cd()
print("Calculating Solar Proton Data ...")
xl = []
xld = []
xs = []
xsd = []
solar_date1, solar_flux1 = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/GP_5m_proton_1MeV.txt', unpack=True)
solar_date2, solar_flux2 = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/GP_5m_proton_100MeV.txt', unpack=True)
fill_date, on = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/nitrogen_data.txt', unpack=True)
juliandate , away, errors = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/away_from_mean.txt', unpack=True)

first_day = juliandate[0]
last_day = juliandate[-1]

x5 = array("d",juliandate)
y5 = array("d",away)

xf = array("d",fill_date)
yf = array("d",on)

detrend_gr = TGraph(len(x5), x5, y5)

for i in range(len(solar_flux1)):
	xl.append(solar_flux1[i]/5000000+0.998)
	xld.append(solar_date1[i])
for j in range(len(solar_flux2)):
	xs.append(solar_flux2[j]/5000000+0.998)
	xsd.append(solar_date2[j])
x_xrayl2 = array("d",xl)
x_xrays2 = array("d",xs)
x_xrayld2 = array("d",xld)
x_xraysd2 = array("d",xsd)
#print(x_xrayl2)
xldg1 = TGraph(len(x_xrayld2), x_xrayld2, x_xrayl2)
xldg1.SetMarkerColor(100)
xldg1.SetLineColor(100)
xsdg1 = TGraph(len(x_xraysd2), x_xraysd2, x_xrays2)
xsdg1.SetMarkerColor(93)
xsdg1.SetLineColor(93)


gr_fill = TGraph(len(xf), xf, yf)
gr_fill.SetMarkerStyle(23)
gr_fill.SetMarkerColor(2)
gr_fill.SetDrawOption("AP")

mg_norm1 = TMultiGraph()
mg_norm1.SetTitle(";Modified Julian Date MST [days]; Total Counts")
mg_norm1.Add(xldg1)
mg_norm1.Add(xsdg1)
mg_norm1.Add(detrend_gr)
#mg_norm1.Add(gr_fill)
mg_norm1.Draw("alp")
mg_norm1.GetYaxis().SetTitleOffset(1.5)
mg_norm1.GetHistogram().SetMaximum(1.0015)      
mg_norm1.GetHistogram().SetMinimum(0.998)
mg_norm1.GetXaxis().SetLimits(first_day,last_day)

legx1 = TLegend(0.65, 0.65, 0.89, 0.89)
legx1.SetFillColor(0)
legx1.AddEntry(detrend_gr, "^{60}Co Count Rate", "lp")
legx1.AddEntry(xldg1, "Proton Flux: 1.0 MeV", "lp")
legx1.AddEntry(xsdg1, "Proton Flux: 100.0 MeV", "lp")
legx1.Draw()

Proton_pad2.Modified()
Proton_pad2.Update()
Proton_pad1.cd()
raw_input("done")