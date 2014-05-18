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

Red = TColor(2000,1,0,0)
Redish = TColor(2001,1,.4,0)
Redishish = TColor(2002,1,.8,0)	
Yellow = TColor(2003,0.5,1,0)	
Yellowish = TColor(2004,0,1,1)		
Orange = TColor(2005,0,.5,1)

first_day = 56506.8942593
last_day = 56673.0953472

fileout_xray1  	= '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_xray_long.txt'
fileout_xray2 	= '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_xray_short.txt'
fileout_proton1 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_proton_1MeV.txt'
fileout_proton2 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_proton_10MeV.txt'

c = TCanvas('canvas','solars',600,1800)
c.Draw()
c.cd()
p1 = TPad('p1','p1',0.0,0.5,1,1)
p1.SetGrid()
p1.Draw()
p2 = TPad('p2','p2',0.0,0.0,1,0.5)
p2.SetGrid()
p2.Draw()
p2.cd()

juliandate , away, errors = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/away_from_mean.txt', unpack=True)
x5 = array("d",juliandate)
y5 = array("d",away)
detrend_gr = TGraph(len(x5), x5, y5)
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
	xl.append(solar_flux_xray_l[i]*70+0.998)
	xld.append(solar_date_xray_l[i])
	xs.append(solar_flux_xray_s[i]*70+0.998)
	xsd.append(solar_date_xray_s[i])
	print(solar_flux_xray_s[i])
x_xrayl2 = array("d",xl)
x_xrays2 = array("d",xs)
x_xrayld2 = array("d",xld)
x_xraysd2 = array("d",xsd)
xldg1 = TGraph(len(x_xrayld2), x_xrayld2, x_xrayl2)
xldg1.SetMarkerColor(2005)
xldg1.SetLineColor(2005)
xsdg1 = TGraph(len(x_xraysd2), x_xraysd2, x_xrays2)
xsdg1.SetMarkerColor(2004)
xsdg1.SetLineColor(2004)

mg_norm1 = TMultiGraph()
mg_norm1.SetTitle(";Modified Julian Date MST [days]; Normalized Count Rate [cps]")
mg_norm1.Add(xldg1)
mg_norm1.Add(xsdg1)
mg_norm1.Add(detrend_gr)
mg_norm1.Draw("al")
mg_norm1.GetYaxis().SetTitleOffset(1.5)
a = mg_norm1.GetYaxis().GetLabelFont()
print(a)
mg_norm1.GetHistogram().SetMaximum(1.0015)      
mg_norm1.GetHistogram().SetMinimum(0.998)
mg_norm1.GetXaxis().SetLimits(first_day,last_day)

axis1 = TGaxis(last_day,0.998,last_day,1.0015,0,0.0000499206,510,"+L")
axis1.SetName("axis1")
axis1.SetLabelColor(1)
axis1.SetTitle("X-Ray Flux [W/m^{2}]")
axis1.SetTitleOffset(1.2)
axis1.SetLabelFont(42)
axis1.SetTitleFont(42)
axis1.Draw()



leg1 = TLegend(0.4, 0.73, 0.89, 0.89)
leg1.SetFillColor(0)
leg1.AddEntry(detrend_gr, "^{60}Co Count Rate", "lp")
leg1.AddEntry(xldg1, "Long Wavelength X-Rays: 0.1 -- 0.8 nm", "lp")
leg1.AddEntry(xsdg1, "Short Wavelength X-Rays: 0.05 -- 0.4 nm", "lp")
leg1.Draw()

p2.Modified()
p2.Update()

'''
#Graph 2
p1.cd()
pl = []
pld = []
ps = []
psd = []

solar_date_proton_l, solar_flux_proton_l = numpy.loadtxt(fileout_proton1, unpack=True)
solar_date_proton_s, solar_flux_proton_s = numpy.loadtxt(fileout_proton2, unpack=True)

for i in range(len(solar_flux_proton_l)):
	pl.append(solar_flux_proton_l[i]/7000000+0.998)
	pld.append(solar_date_proton_l[i])
	ps.append(solar_flux_proton_s[i]/7000000+0.998)
	psd.append(solar_date_proton_s[i])
print(max(pl))
x_protonl2 = array("d",pl)
x_protons2 = array("d",ps)
x_protonld2 = array("d",pld)
x_protonsd2 = array("d",psd)

pldg1 = TGraph(len(x_protonld2), x_protonld2, x_protonl2)
pldg1.SetMarkerColor(2000)
pldg1.SetLineColor(2000)
psdg1 = TGraph(len(x_protonsd2), x_protonsd2, x_protons2)
psdg1.SetMarkerColor(2001)
psdg1.SetLineColor(2001)

mg_norm2 = TMultiGraph()
mg_norm2.SetTitle("new;Modified Julian Date MST [days]; Normalized Count Rate [cps]")
mg_norm2.Add(pldg1)
mg_norm2.Add(psdg1)
mg_norm2.Add(detrend_gr)
mg_norm2.Draw("al")

mg_norm2.GetYaxis().SetTitleOffset(1.5)
mg_norm2.GetHistogram().SetMaximum(1.0015)      
mg_norm2.GetHistogram().SetMinimum(0.998)
mg_norm2.GetXaxis().SetLimits(first_day,last_day)

axis2 = TGaxis(last_day,0.998,last_day,1.0015,0,24490.1*10000,510,"+L")
axis2.SetName("axis2")
axis2.SetLabelColor(1)
axis2.SetTitle("Proton Flux [p^{+}/m^{2}-s-sr]")
axis2.SetTitleOffset(1.2)
axis2.SetLabelFont(42)
axis2.SetTitleFont(42)
axis2.Draw()

leg2 = TLegend(0.45, 0.73, 0.89, 0.89)
leg2.SetFillColor(0)
leg2.AddEntry(detrend_gr, "^{60}Co Count Rate", "lp")
leg2.AddEntry(pldg1, "Low Energy Protons: 1 MeV", "lp")
leg2.AddEntry(psdg1, "High Energy Protons: 10 MeV", "lp")
leg2.Draw()

p1.Modified()
p1.Update()

#Graph 3
fileout_electron1 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_electron_0.8MeV.txt'
fileout_electron2 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_electron_2MeV.txt'

c2 = TCanvas('canvas2','solars',600,1800)
c2.Draw()
c2.cd()
p12 = TPad('p12','p1',0.0,0.5,1,1)
p12.SetGrid()
p12.Draw()
p22 = TPad('p22','p2',0.0,0.0,1,0.5)
p22.SetGrid()
p22.Draw()
p22.cd()

el = []
eld = []
es = []
esd = []

solar_date_electron_l  , solar_flux_electron_l = numpy.loadtxt(fileout_electron1, unpack=True)
solar_date_electron_s 	, solar_flux_electron_s = numpy.loadtxt(fileout_electron2, unpack=True)


for k in range(len(solar_flux_electron_l)):
	el.append(solar_flux_electron_l[k]/50000000+0.998)
	eld.append(solar_date_electron_l[k])
	es.append(solar_flux_electron_s[k]/50000000+0.998)
	esd.append(solar_date_electron_s[k])
	print(solar_flux_electron_s[k])

x_electronl2 = array("d",el)
x_electrons2 = array("d",es)
x_electronld2 = array("d",eld)
x_electronsd2 = array("d",esd)
eldg1 = TGraph(len(x_electronld2), x_electronld2, x_electronl2)
eldg1.SetMarkerColor(2003)
eldg1.SetLineColor(2003)
esdg1 = TGraph(len(x_electronsd2), x_electronsd2, x_electrons2)
esdg1.SetMarkerColor(2002)
esdg1.SetLineColor(2002)

mg_norm3 = TMultiGraph()
mg_norm3.SetTitle(";Modified Julian Date MST [days]; Normalized Count Rate [cps]")
mg_norm3.Add(eldg1)
mg_norm3.Add(esdg1)
mg_norm3.Add(detrend_gr)
mg_norm3.Draw("al")
mg_norm3.GetYaxis().SetTitleOffset(1.5)
mg_norm3.GetHistogram().SetMaximum(1.0015)      
mg_norm3.GetHistogram().SetMinimum(0.998)
mg_norm3.GetXaxis().SetLimits(first_day,last_day)

axis3 = TGaxis(last_day,0.998,last_day,1.0015,0,174300*10000,510,"+L")
axis3.SetName("axis3")
axis3.SetLabelColor(1)
axis3.SetTitle("Electron Flux [e^{-}/m^{2}-s-sr ]")
axis3.SetLabelFont(42)
axis3.SetTitleFont(42)
axis3.SetTitleOffset(1.2)
axis3.Draw()


leg3 = TLegend(0.45, 0.73, 0.89, 0.89)
leg3.SetFillColor(0)
leg3.AddEntry(detrend_gr, "^{60}Co Count Rate", "lp")
leg3.AddEntry(eldg1, "Low Energy Electrons: 0.8 MeV", "lp")
leg3.AddEntry(esdg1, "High Energy Electrons: 2.0 MeV", "lp")
leg3.Draw()

p22.Modified()
p22.Update()
'''
raw_input("done")