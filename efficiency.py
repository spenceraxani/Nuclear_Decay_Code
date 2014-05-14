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

########################################
#This was just a quick piece of code to calculate the change in efficiency with repect to a change in distance of the HPGe, 
########################################
cproton = TCanvas('cproton', 'Proton',600,900)
cproton.Draw()
cproton.cd()
Proton_pad2 = TPad('x11','p',0.0,0.0,1,0.5)
Proton_pad2.SetGrid()
Proton_pad2.Draw()
Proton_pad2.cd()

position, net = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/efficiency.txt', unpack=True) #were the data is found.
eff_117 = []
eff_117_error = []
eff_133 = []
eff_133_error = []
eff_total = []
eff_total_error = []
pos = []
pos_error = []
for i in range(len(position)):
	#eff_117.append(net_117[i]/30000000.0*100)
	#eff_117_error.append(np.sqrt(net_117[i])/30000000.0*100)
	#eff_133.append(net_133[i]/30000000.0*100)
	#eff_133_error.append(np.sqrt(net_133[i])/30000000.0*100)
	eff_total.append(net[i]/50000000.0*1)
	eff_total_error.append(np.sqrt(net[i])/50000000.0*1)
	pos.append(position[i]*1.0)
	pos_error.append(0.0)
print(eff_117_error)
[36.88, 73.76,110.6,147.5]
#y2 = array("d",eff_133)
#erry2 = array("d",eff_133_error)
y3 = array("d",eff_total)
erry3 = array("d",eff_total_error)
x = array("d",pos)
errx = array("d",pos_error)

gr = TGraphErrors(len(x),x,y3,errx,erry3)
gr.SetTitle("TGraphErrors Example")
gr.SetMarkerColor(4)
gr.SetMarkerStyle(21)
gr.Draw("AP")
gr.SetMarkerStyle(24)
gr.SetMarkerSize(1.7)
fsignal = TF1("FSignal","pol1", -800 ,800)
gr.Fit(fsignal,'R')
gr.SetTitle("; Distance from Central Position [#mum]; Combined Efficiency [%]")
gr.GetYaxis().SetTitleOffset(1.4)
Proton_pad2.Modified()
Proton_pad2.Update()
raw_input("done")

