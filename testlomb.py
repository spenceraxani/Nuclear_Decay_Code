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
from ROOT import gROOT, TCanvas, TF1, TGraph, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad,THStack, TAxis,TRandom3, TH1F, gStyle, TH1D,TMath, TSpectrum
from array import array
from random import random
from random import gauss
from decimal import *
from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt

#time , zero, a,b,v,c,x,s,counts, livetime, fwhm, excess = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/real_data.txt', unpack=True)
time , counts, errs= numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/away_from_mean.txt', unpack=True)
x = array("d",time)
time_series = np.array(x, np.float64)
y = array("d",counts)
measurement = np.array(y, np.float64)
f = np.linspace(0.0001, 2, 10000)

print("interesting peak at w = 0.34, which is equivalent to a peariod of: T = "+ str(2 * np.pi / 0.35))

lombs = sp.signal.lombscargle(time_series , measurement,f)
lombs_array = array("d",lombs)
frequency_array = array("d",f)
print(type(lombs))
plt.subplot(2, 1, 1)
plt.plot(f,lombs)
plt.subplot(2, 1, 2)
plt.plot(time, counts)
plt.show()






'''

list_time = []
list_counts = []
t = 30
w = 2 * np.pi/t
phi = 0.5 * np.pi

print(w)
for i in range(1000):
	list_time.append(i)
	list_counts.append(np.sin(w * i))

time_series = np.array(list_time, np.float64)
measurement = np.array(list_counts, np.float64)

f = np.linspace(0.01, 2, 10000)
lombs = sp.signal.lombscargle(time_series , measurement,f)
plt.subplot(2, 1, 2)
plt.plot(f,lombs)
plt.subplot(2, 1, 1)
plt.plot(list_time, list_counts)
plt.show()


list_time = []
list_counts = []
t = 30
w = 2 * np.pi/t
phi = 0.5 * np.pi

print(w)
for i in range(1000):
	list_time.append(i)
	list_counts.append(gauss(1.17, 0.018)+0.02*np.sin(w*i))

time_series = np.array(list_time, np.float64)
measurement = np.array(list_counts, np.float64)

f = np.linspace(0.01, 2, 10000)
lombs = sp.signal.lombscargle(time_series , measurement,f)
plt.subplot(2, 1, 2)
plt.plot(f,lombs)
plt.yscale('log')
plt.subplot(2, 1, 1)
plt.plot(list_time, list_counts)
plt.show()

'''









'''
#b = array(list1)
for i in range(10):
	gauss(1.17, 0.018)
list2 = []
list1 = []
b = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/CO60_1_.txt',dtype=numpy.ndarray)
print(b)
a = ([])
for i in range(len(b)):
	a = np.append(i)
#print(a)
'''
'''
#xdata , zero, a,b,v,c,x,s,ydata, livetime, fwhm, excess = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/real_data.txt', unpack=True)
x1data , y1data, errs = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/away_from_mean.txt', unpack=True)

another_list = [] 
list_count = []
for i in range(1000):
	another_list.append(numpy.sin(i*6.28/30.0)*2*numpy.sin(i*6.28/60.0))
	#another_list.append(1000*(gauss(1, 0.002)+0.01*numpy.sin(i*6.28/30)))
	list_count.append(i/100.0)

x = array("d",another_list)
y = array("d",list_count)
x_array = np.array(x, np.float64)
y_array = np.array(y, np.float64)

f = np.linspace(0.001, 100, 10000)
pgram = sp.signal.lombscargle(x_array, y_array, f)

plt.subplot(2, 1, 1)
plt.plot(list_count, another_list)
plt.subplot(2, 1, 2)

plt.plot(f, pgram)

plt.show()






#x, y, err = np.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/away_from_mean.txt', usecols = (0,1), unpack=True)
y = y - y.mean()

W = fftfreq(y.size, d=(x[1]-x[0])*86400)
plt.subplot(2,1,1)
plt.plot(x,y)
plt.xlabel('Time (days)')

f_signal = fft(y)
plt.subplot(2,1,2)
plt.plot(W, abs(f_signal)**2)
plt.xlabel('Frequency (Hz)')

plt.xscale('log')
plt.xlim(10**(-6), 10**(-5))
plt.show()

A = 2.
w = 1.
phi = 0.5 * np.pi
nin = 1000
nout = 100000
frac_points = 0.9 # Fraction of points to select
r = np.random.rand(nin)

x = np.linspace(0.01, 10*np.pi, nin)
x = x[r >= frac_points]

#print(x)
normval = x.shape[0] # For normalization of the periodogram
y = A * np.sin(w*x+phi)

f = np.linspace(0.01, 10, nout)
pgram = sp.signal.lombscargle(x, y, f)
plt.subplot(2, 1, 1)

plt.plot(x, y, 'b+')
plt.subplot(2, 1, 2)

plt.plot(f, np.sqrt(4*(pgram/normval)))

plt.show()
'''