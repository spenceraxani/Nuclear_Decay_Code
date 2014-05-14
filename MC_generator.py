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
from ROOT import gROOT, TGaxis,TCanvas, TF1, TGraph, TLine, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt

#####################################




MC = 0
ang_freq = 200
trials = 10000
MC_fileout = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/MC_data.txt'
if MC == True:
	try:
	    os.remove(MC_fileout)	
	except OSError:
		pass
	time , counts,errs = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/hist_mean_fileout.txt', unpack=True)
	time_list = []
	sample_x = []
	mean = np.mean(counts)
	for i in range(len(time)):
		time_list.append(time[i]-time[0])
		sample_x.append((counts[i]-mean)*1000)

	sample_mean = np.mean(sample_x)*1.0
	variance = 0.0
	for k in range(len(sample_x)):
		variance += (sample_x[k] - sample_mean)**2
	sample_deviation = np.sqrt((1.0/len(sample_x))*variance)
	print("The standard deviation of input file is: " + str(sample_deviation))
	print("The input file mean is: " + str(sample_mean))
	time_test = []
	sample_test = []
	counter = 0 
	MCfileout = open(MC_fileout,'a')
	maximum = 0
	for j in range(trials):
		print(counter)
		counter +=1
		for i in range(len(time)):
			time_test.append(time[i]-time[0])
			sample_test.append(np.random.normal(0, sample_deviation))
		x = array("d",time_test)
		t1 = np.array(x, np.float64)
		y = array("d",sample_test)
		x1 = np.array(y, np.float64)
		f = np.linspace(0.01, 1, ang_freq)
		lombs = sp.signal.lombscargle(t1, x1,f)
			#MCfileout.write(str(f.item(l))+ "\t" +str(np.round(lombs.item(l),3)) + "\n" )
		MCfileout.write(str(np.amax(lombs))+"\n" )
		time_test = []
		sample_test = []
	MCfileout.close()
cscar = TCanvas('cscar','scar',600,1800)
cscar.Draw()
cscar.cd()
pscar1 = TPad('pscar1','p',0.0,0.0,1,.5)
pscar1.SetGrid()
pscar1.Draw()
pscar2 = TPad('pscar2','p',0.0,0.5,1,1)
pscar2.SetGrid()
pscar2.Draw()

bins = 60
max_amp = 30

pscar1.cd()
z1 = 0.317311
z2 = 0.045501
z3 = 0.0027
z5 = float(0.0000006004)
v5a = - np.log(1-np.power(1.0-z5,1.0/200.0))

pscar1.cd()

pscar1.SetLogy()
h1 = TH1F('Maximums','scargy',bins,0,max_amp)
h1.SetTitle(";Maximum Periodogram Amplitude Detected; Counts")
h1.SetLineColor(1)
h1.Draw()
spectrum = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/MC_total.txt', unpack=True)
obs_max = np.amax(spectrum)
h1.SetStats(False)
for i in range(len(spectrum)):
	h1.Fill(spectrum[i],1.0)
print("The maximum peak out of the " +str(trials)+" trials was " +str(np.amax(spectrum)))

bins_to_amp = 1.0*bins/max_amp

total_bins = 0
for j in range(bins):
	total_bins += h1.GetBinContent(j)
	if total_bins <= len(spectrum)*(1.0-z1):
		v1 = j/bins_to_amp
		height1 = h1.GetBinContent(j)
	
	if total_bins <= len(spectrum)*(1.0-z3):
		v3 = j/bins_to_amp
		height3 = h1.GetBinContent(j)
	if total_bins <= len(spectrum)*(1.0-z5):
		v5 = j/bins_to_amp
		height5 = h1.GetBinContent(j)
		#print("tata"+str(v5))
print(str(v1))
#v5 = obs_max
#height5 = 1.0/len(spectrum)

p1 = TLine(v1,0,v1,height1)
p1.SetLineStyle(2)
p1.SetLineColor(40)
p1.SetLineWidth(3)
p3 = TLine(v3,0,v3,height3)
p3.SetLineStyle(2);
p3.SetLineColor(30)
p3.SetLineWidth(3)
p5 = TLine(v5,0,v5,height5);
p5.SetLineStyle(2);
p5.SetLineColor(45)
p5.SetLineWidth(3)
p5a = TLine(v5a,0,v5a,height5);
p5a.SetLineStyle(2);
p5a.SetLineColor(1)
p5a.SetLineWidth(3)
pblack = TLine()
legMC = TLegend(0.59,0.66,0.88,0.88)
legMC.SetFillColor(0)
legMC.AddEntry(pblack,"Lomb Scargle Maximum Amplitude Data","l")
legMC.AddEntry(p1,"1 #sigma","l")
legMC.AddEntry(p3,"3 #sigma","l")
legMC.AddEntry(p5,"5 #sigma","l")
legMC.Draw()

p1.Draw("same")
p3.Draw("same")
p5.Draw("same")

pscar1.Modified();
pscar1.Update()

pscar2.cd()
#time , zeroes , er117, c133, er133, net, ernet, exp, realcount, livetime, fwhm, counts = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/real_data.txt', unpack=True)

time , counts, errs = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/hist_mean_fileout.txt', unpack=True)
time_list = []
sample_x = []
mean = np.mean(counts)
for i in range(len(time)):
	time_list.append(time[i]-time[0])
	sample_x.append((counts[i]-mean)*1000)
x = array("d",time_list)
time_series = np.array(x, np.float64)
y = array("d",sample_x)
measurement = np.array(y, np.float64)
f = np.linspace(0.01, 1, ang_freq)
norm = time_series.shape[0]
lombs = sp.signal.lombscargle(time_series , measurement,f)

h2 = TH1F('Scargle Amplitude1','scargy1',200,0,1)
h2.SetTitle("Scargle Amplitude;Angular Frequency [rads/day]; Arbitrary Units ")

line_1s = TLine(0,v1,1,v1);
line_1s.SetLineStyle(2)
line_1s.SetLineColor(40)
line_1s.SetLineWidth(4)
line_3s = TLine(0,v3,1,v3);
line_3s.SetLineStyle(2);
line_3s.SetLineColor(30)
line_3s.SetLineWidth(4)
line_5s = TLine(0,v5,1,v5);
line_5s.SetLineStyle(2);
line_5s.SetLineColor(45)
line_5s.SetLineWidth(3)


print(v3)

for i in range(len(lombs)):
	h2.Fill(f[i],lombs[i])
s = TSpectrum(4,1)
nfound = s.Search(h2,1,"new")
peaks = s.GetPositionX()
height = s.GetPositionY()


h2.Draw('l')
h2.SetLineColor(1)

line_3s.Draw("Same")
line_1s.Draw("Same")
line_5s.Draw("Same")

print("Found candidate peaks = " + str(peaks[0]) + 'keV with ' +str(height[0]) +' counts, and '+ str(peaks[1]) + 'keV with ' +str(height[1])+ ' counts')
h2.SetStats(False)
axis6 = TGaxis(2,0,2,0.4,0,2,50510,"+L")
axis6.Draw()

leghist = TLegend(0.59,0.66,0.88,0.88)
leghist.SetFillColor(0)
leghist.AddEntry(pblack,"Lomb Scargle Periodogram","l")
leghist.AddEntry(line_1s,"1 #sigma","l")
leghist.AddEntry(line_3s,"3 #sigma","l")
leghist.AddEntry(line_5s,"5 #sigma","l")
leghist.Draw()

pscar2.Modified();
pscar2.Update()
raw_input("done")
