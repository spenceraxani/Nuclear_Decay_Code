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
first_day = 56506.8942593
last_day = 56673.0953472
def fit_counts(B,A,t):
	counts_exp = A * numpy.exp( B * t)
	return(counts_exp)

SOLAR_DATA = 0
AVE=0

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
juliandate , c117 , zeros, c133, er133, net, ernet, exp , realcount, livetime, fwhm, excesscounts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/real_data.txt', unpack=True)
ave_out_file = 'ave_detector_counts.txt'
if AVE == True:
	try:
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + ave_out_file)
	except OSError:
		pass
	num_files = 0
	list_net = []
	list_date = []
	list_error = []
	ave_sum = 0
	ave_error = 0
	ave_date = 0
	for k in range(len(juliandate)):
		num_files += 1
		list_net.append(realcount[k])
		list_date.append(juliandate[k])	
		if num_files == 50:
			for i in range(len(list_net)):
				ave_sum += list_net[i]
				ave_date += list_date[i]
			ave_error = math.sqrt(ave_sum)
			ave_error = ave_error/num_files
			ave_sum = ave_sum/num_files
			ave_date = ave_date/num_files
			ave_outfile = open(ave_out_file,'a')
			ave_outfile.write(str(ave_date) + "\t" + str(round_sig(ave_sum/livetime[k]))+ "\t" + str(round_sig(ave_error/livetime[k]))+"\n")
			num_files = 0
			ave_sum = 0
			ave_date = 0
			list_net = []
			list_date = []			
juliandate , counts, error = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/ave_detector_counts.txt', unpack=True)

for k in range(len(juliandate)):
	diff = 1 - (fit_counts(B,A,juliandate[k]) - counts[k])/fit_counts(B,A,juliandate[k])
	errdifftop = 1 + error[k]/fit_counts(B,A,juliandate[k])
	errdiffbottom = 1 - error[k]/fit_counts(B,A,juliandate[k])
	away_fileout = open('/Users/spenceraxani/Documents/Nuclear_Decay/Data/away_from_mean.txt','a')
	away_fileout.write(str(juliandate[k]) + "\t" + str(diff)+"\t" + str(1-errdiffbottom)+ "\n" )
for k in range(10000):
	away_error = open('/Users/spenceraxani/Documents/Nuclear_Decay/Data/away_error_fileout.txt','a')
	away_error.write(str(first_day+0.1*k) + "\t" +str(1)+ "\t" + str(0) +"\t" + str(math.fabs(1-errdiffbottom))  +"\t" +"\t" +  str(math.fabs(2*(1-errdiffbottom)))+"\n" )
plot_ave = TGraph()
plot_ave.SetLineColor(2)
for k in range(len(juliandate)):
	plot_ave.SetPoint(k,juliandate[k],counts[k])
plot_ave.GetYaxis().SetTitle("Normalized Count Rate [cps]")
plot_ave.GetXaxis().SetTitle("Modified Julian Date MST [Days]")
plot_ave.GetXaxis().SetLimits(first_day,last_day); 
plot_ave.SetDrawOption("l")
plot_ave.SetFillStyle(3013)
plot_ave.SetLineColor(1)
plot_ave.SetLineWidth(2)
plot_ave.Draw("al")
plot_ave.GetYaxis().SetTitleOffset(1.6);
p1ave.Update()


cdetre = TCanvas('cdetrend', 'detrend',600,1800)
cdetre.Draw()
cdetre.cd()
pd1 = TPad('x10','p',0.0,0.5,1,1)
pd1.SetGrid()
pd1.Draw()
pd2 = TPad('x11','p',0.0,0.0,1,0.5)
pd2.SetGrid()
pd2.Draw()
pd2.cd()

date_1 , data_detrended = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended.txt', unpack=True)

plot_detrend = TGraph()
plot_detrend.SetLineColor(2)
counter = 0
ave_detrend = []
ave_date = []
detrended_average = []
detrended_time = []
for i in range(len(date_1)):
	counter += 1
	ave_detrend.append(data_detrended[i])
	ave_date.append(date_1[i])
	if counter == 10:
		count = np.mean(ave_detrend)
		detrended_average.append(count)
		time = np.mean(ave_date)
		detrended_time.append(time)
		ave_detrend = []
		ave_date = []
		counter = 0
x = np.array(detrended_time, np.float64)
y = np.array(detrended_average, np.float64)

for k in range(len(x)):
	plot_detrend.SetPoint(k,x[k],y[k])
plot_detrend.GetYaxis().SetTitle("Count Rate [cps]")
plot_detrend.GetXaxis().SetTitle("Modified Julian Date MST [Days]")
plot_detrend.GetXaxis().SetLimits(first_day,last_day); 
plot_detrend.SetDrawOption("l")
plot_detrend.SetFillStyle(3013)
plot_detrend.SetLineColor(1)
plot_detrend.SetLineWidth(2)
plot_detrend.Draw("ap")
plot_detrend.GetYaxis().SetTitleOffset(1.6);
pd2.Update()

cxray = TCanvas('cxray', 'X-Ray',600,1800)
cxray.Draw()
cxray.cd()
px1 = TPad('x10','p',0.0,0.5,1,1)
px1.SetGrid()
px1.Draw()
px2 = TPad('x11','p',0.0,0.0,1,0.5)
px2.SetGrid()
px2.Draw()
print("Calculating Solar X-Ray Data ...")

fileout_xray1 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_xray_long.txt'
fileout_xray2 = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned/binned_xray_short.txt'
if SOLAR_DATA == True:
	try:
		os.remove(fileout_xray1)
		os.remove(fileout_xray2)
	except OSError:
		pass
		
	for file in glob.glob('/Users/spenceraxani/Documents/Nuclear_Decay/Data/xray/*_Gp_xr_5m.txt'):
		#print(file)
		list = file
		fh = open(file)
		lines = fh.readlines()
		for line in lines:
			columns = line.split(' ')
			columns = [col.strip() for col in columns]
			if columns[0] == "2013":
				#print(str(columns))
				output = open(fileout_xray1,'a')
				if float(columns[-1]) > 0 and float(columns[-5]) > 0 : #Xray long
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000- 0.2916667) + "\t" +  str(columns[-5]) + "\n")
				else:
					continue
				
				output = open(fileout_xray2,'a') #Xray Short
				if float(columns[-1]) > 0 and float(columns[-9]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000 - 0.2916667)+ "\t" +  str(columns[-9]) + "\n")
				else:
					continue
					
			if columns[0] == "2014":
				#print(str(columns))
				output = open(fileout_xray1,'a')
				if float(columns[-1]) > 0 and float(columns[-5]) > 0 : #Xray long
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000- 0.2916667) + "\t" +  str(columns[-5]) + "\n")
				else:
					continue
				
				output = open(fileout_xray2,'a') #Xray Short
				if float(columns[-1]) > 0 and float(columns[-9]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000- 0.2916667) + "\t" +  str(columns[-9]) + "\n")
				else:
					continue

solar_date_xray_l  , solar_flux_xray_l = numpy.loadtxt(fileout_xray1, unpack=True)
solar_date_xray_s 	, solar_flux_xray_s = numpy.loadtxt(fileout_xray2, unpack=True)

x_xray_l = array("d",solar_date_xray_l)
y_xray_l = array("d",solar_flux_xray_l)
x_xray_s = array("d",solar_date_xray_s)
y_xray_s = array("d",solar_flux_xray_s)
solar_graphs_xray_l= TGraph(len(x_xray_l), x_xray_l, y_xray_l)
solar_graphs_xray_l.SetMarkerColor(2005)
solar_graphs_xray_l.SetLineColor(2005)
solar_graphs_xray_s = TGraph(len(x_xray_s), x_xray_s, y_xray_s)
solar_graphs_xray_s.SetMarkerColor(2004)
solar_graphs_xray_s.SetLineColor(2004)
px2.cd()
mg_xray = TMultiGraph()
mg_xray.SetTitle("Solar Flare Xray Flux Data GOES15 Satelite;Modified Julian Date MST [days]; X-Ray Flux [Watts/m^{2}]")
mg_xray.Add(solar_graphs_xray_l)
mg_xray.Add(solar_graphs_xray_s)

#mg_xray.SetMaximum(0.00001);           
mg_xray.SetMinimum(0);

mg_xray.Draw("ALp")

mg_xray.GetXaxis().SetLimits(56546.9071712,56672.8875958) 
px2.Modified();
px2.Update()
xl = []
xld = []
xs = []
xsd = []
print(len(solar_flux_xray_l))
print(len(solar_flux_xray_s))
for i in range(len(solar_flux_xray_s)):
	xl.append(solar_flux_xray_l[i]*250000.0 +980)
	xld.append(solar_date_xray_l[i])
	xs.append(solar_flux_xray_s[i]*250000.0 +980)
	xsd.append(solar_date_xray_s[i])
x_xrayl = array("d",xl)
x_xrays = array("d",xs)
x_xrayld = array("d",xld)
x_xraysd = array("d",xsd)
xldg = TGraph(len(x_xrayld), x_xrayld, x_xrayl)
xldg.SetMarkerColor(2005)
xldg.SetLineColor(2005)
xsdg = TGraph(len(x_xraysd), x_xraysd, x_xrays)
xsdg.SetMarkerColor(2004)
xsdg.SetLineColor(2004)
px1.cd()
mg_norm = TMultiGraph()
mg_norm.SetTitle(";Modified Julian Date MST [days]; Count Rate [cps]")
fxraysignal = TF1("FSignal","expo", first_day ,last_day)
fxraysignal.SetLineWidth(1)
plot_detrend.Fit(fxraysignal,'R')
mg_norm.Add(xldg)
mg_norm.Add(xsdg)
mg_norm.Add(plot_detrend)
mg_norm.Draw("AL")
mg_norm.GetYaxis().SetTitleOffset(1.4)
mg_norm.GetHistogram().SetMaximum(1030);        
mg_norm.GetHistogram().SetMinimum(980); 
mg_norm.GetXaxis().SetLimits(56546.9071712,56672.8875958)

legx = TLegend(0.65, 0.65, 0.89, 0.89)
legx.SetFillColor(0)
legx.AddEntry(plot_ave, "^{60}Co Count Rate", "lp")
legx.AddEntry(fxraysignal, "Exponential Fit", "lp")
legx.AddEntry(xldg, "Long Wavelength X-Ray: 0.1 - 0.8 nm", "lp")
legx.AddEntry(xsdg, "Short Wavelength X-Ray: 0.05 - 0.4 nm", "lp")
legx.Draw()
px1.Modified();
px1.Update()

p2ave.Modified();
p2ave.Update()

c7x = TCanvas('c7x','Peakx',600,1800)
c7x.Draw()
c7x.cd()
p70x = TPad('p70x','px',0.0,0.5,1,1)
p70x.SetGrid()
p70x.Draw()
p71x = TPad('p71x','px',0.0,0.0,1,0.5)
p71x.SetGrid()
p71x.Draw()
p71x.cd()
xl = []
xld = []
xs = []
xsd = []

solar_date_xray_l  , solar_flux_xray_l = numpy.loadtxt(fileout_xray1, unpack=True)
solar_date_xray_s 	, solar_flux_xray_s = numpy.loadtxt(fileout_xray2, unpack=True)

juliandate , away, errors = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/away_from_mean.txt', unpack=True)
x5 = array("d",juliandate)
y5 = array("d",away)
detrend_gr = TGraph(len(x5), x5, y5)

for i in range(len(solar_flux_xray_l)):
	xl.append(solar_flux_xray_l[i]*10+0.998)
	xld.append(solar_date_xray_l[i])
	xs.append(solar_flux_xray_s[i]*10+0.998)
	xsd.append(solar_date_xray_s[i])
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
mg_norm1.SetTitle(";Modified Julian Date MST [days]; Count Rate [cps]")
mg_norm1.Add(xldg1)
mg_norm1.Add(xsdg1)
mg_norm1.Add(detrend_gr)
mg_norm1.Draw("al")
mg_norm1.GetYaxis().SetTitleOffset(1.5)
mg_norm1.GetHistogram().SetMaximum(1.0015)      
mg_norm1.GetHistogram().SetMinimum(0.998)
mg_norm1.GetXaxis().SetLimits(first_day,last_day)

legx1 = TLegend(0.65, 0.65, 0.89, 0.89)
legx1.SetFillColor(0)
legx1.AddEntry(detrend_gr, "^{60}Co Count Rate", "lp")
legx1.AddEntry(xldg1, "Long Wavelength X-Ray: 0.1 - 0.8 nm", "lp")
legx1.AddEntry(xsdg1, "Short Wavelength X-Ray: 0.05 - 0.4 nm", "lp")
legx1.Draw()

p71x.Modified()
p71x.Update()
p70x.cd()
raw_input("done")