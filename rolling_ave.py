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
from ROOT import gROOT, TGaxis, TCanvas, TF1, gPad, TGraph, TGaxis, TLine, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum
from array import array
from random import random
from decimal import *
from random import gauss
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt
import time
getcontext().prec = 12
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
first_day = 56506.8942593
last_day = 56673.0953472
juliandate , counts, error = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/ave_detector_counts.txt', unpack=True)

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
new_graph = TF1("sin","sin((2*3.1415*x/38.37)+3.9)",first_day,last_day)
new_graph.Draw("al")
new_graph.GetYaxis().SetTitleOffset(1.6);
new_graph.GetHistogram().SetMaximum(2)      
new_graph.GetHistogram().SetMinimum(-2)
#plot_ave.Draw("al")
#plot_ave.GetYaxis().SetTitleOffset(1.6);
p1ave.Update()




p2ave.cd()
juliandate , away, errors = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/away_from_mean.txt', unpack=True)
x = array("d",juliandate)
y = array("d",away)
plot_away = TGraph(len(x), x, y)
plot_away.SetMarkerStyle(7)

inc , one, errorx, errory , doubleerror = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/away_error_fileout.txt', unpack=True)
j = array("d",inc)
r = array("d",one)
q = array("d",errorx)
e = array("d",errory)
w = array("d",doubleerror)
#print(type(j))
plot_awaymin = TGraphErrors(len(j), j, r , q , e)
plot_awaymin.SetLineColor(40);
plot_awaymin.SetMarkerColor(40);
plot_awaymin.SetLineWidth(-802);
plot_awaymax = TGraphErrors(len(j), j, r , q , w)
plot_awaymax.SetLineColor(30);
mg_away = TMultiGraph()
mg_away.SetTitle("Normalized Counts;Modified Julian Date MST [days]; Normalized Counts")
#mg_away.SetTitle("Probability;Residual [cps]; Probability")

mg_away.Add(plot_awaymax,"p")
mg_away.Add(plot_awaymin,"p")
mg_away.Add(plot_away,"lp")
mg_away.Draw("a")

mg_away.GetYaxis().SetTitleOffset(1.5);
mg_away.GetXaxis().SetLimits(first_day,last_day);



axis7 = TGaxis(last_day,0,last_day,1,0,2000,50510,"+L")
axis7.SetName("axis7")
#axis7.SetLabelOffset(0.01)
axis7.Draw()

p2ave.Modified()
p2ave.Update()

c1 = TCanvas("c1","Examples of Gaxis",10,10,700,500);

c1.Range(-10,-1,10,1)
axis1 = TGaxis(-4.5,-0.2,5.5,-0.2,-6,8,510,"")
axis1.SetName("axis1")
axis1.Draw()

axis2 = TGaxis(-4.5,0.2,5.5,0.2,0.001,10000,510,"G")
axis2.SetName("axis2")
axis2.Draw()



raw_input("done")