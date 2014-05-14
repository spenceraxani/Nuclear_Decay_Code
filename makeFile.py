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

##########################################################
# This is just a mess of code to generate a tonne of graphs.  It also takes all the solar data and count rate data to make it usefull.  You'll spend most of your time in this file. 
###########################################################



SOLAR_DATA = 0
TEMP_DATA = 0
SYNC = 1
CHI = 0
LIKELIHOOD = 0
FLAT = 0
REAL = 0
PRESSURE = 0
FWHM = 0
CHIgraph = 0
AVE = 0
THREED = 0
TEMPERATURE = 0
LINE = 0
MODTEMP = 0
SUN_DISTANCE=0
LOMB = 0
HIST = 0
MEAN =0
MC = 0
BUN = 1

gStyle.SetEndErrorSize(0)
gStyle.SetOptStat("nei")

##########################################################################################
### Goodness of Fit
##########################################################################################		
def drange(start, stop, step):
	r = start
	while r < stop:
		yield r
		r += step

def exp_counts(t):
	half_life = 1925.28 	#days	
	initial_counts = 1043.57
	counts_exp = (initial_counts)*numpy.exp(-numpy.log(2)*(t-56506.8942593)/half_life)
	return(counts_exp)

def solar_distance(t):
	perihelion = numpy.float128(147166462000.0)	
	aphelion =  numpy.float128(152171522000.0)
	eccentricity = numpy.float128(1 - 2/((aphelion/perihelion)+1))
	#print(eccentricity)
	pi = 3.14159265397
	distance_to_sun = numpy.float128(aphelion*(1-math.pow(eccentricity,2))/(1+eccentricity*math.cos((2*pi/365)*(t - 56478.5))))
	return(distance_to_sun)
	
def fit_counts(B,A,t):
	counts_exp = A * numpy.exp( B * t)
	return(counts_exp)

def toy_counts(t): 
	half_life = 1925.28 	#days	
	initial_counts = 1036
	counts_toy = (initial_counts)*numpy.exp(-numpy.log(2)*(t)/half_life) + .1*numpy.sin(100*t/numpy.pi + numpy.pi/4)
	return(counts_toy)
	
def iround(x):
    return int(round(x) - .5) + (x > 0)

def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line
            
def round_sig(x, sig=6):
	return round(x, sig-int(floor(log10(x)))-1)

def correlation(input):
	#input needs to be the entire path location and formated as a text file of the form x \t y
	x, y = numpy.loadtxt(input, unpack=True)
	mean_x = 0.0
	mean_y = 0.0
	N = float(len(x))
	sum_x = 0.0
	sum_y = 0.0
	sum_x_squared = 0.0
	sum_y_squared = 0.0
	xy = 0.0
	diff_x = 0.0
	diff_y = 0.0
	diff_x_squared = 0.0
	diff_y_squared = 0.0
	diff_x_diff_y = 0.0

	for j in range(len(x)):
		sum_x = x[j]
		sum_y = y[j]
	mean_x = sum_x / N
	mean_y = sum_y / N
	print(90+1.0)
	for i in range(len(x)):
		sum_x_squared += x[i]**2
		sum_y_squared += y[i]**2
		xy += x[i]*y[i]
		diff_x += x[i] - mean_x
		diff_y += y[i] - mean_y
		diff_x_squared = (x[i] - mean_x)**2
		diff_y_squared = (y[i] - mean_y)**2
		diff_x_diff_y = (x[i] - mean_x)*(y[i] - mean_y)
	Sx = (diff_x_squared/(N-1))**(0.5)
	Sy = (diff_y_squared/(N-1))**(0.5)
	Sxy = diff_x_diff_y / (N-1)
	Cov = (xy - N*mean_y*mean_x)/N-1
	r = (xy - N*mean_x*mean_y)/((sum_x_squared - sum_x**2/N)*(sum_y_squared - sum_y**2/N))**(0.5)
	print("The covariance is: " + str(Sxy))
	print("The correlation coefficient is: " + str(r))
	print("The coefficient of determination is: " + str(r**2)) 

	return(r)

##########################################################################################
###	Data statistics
##########################################################################################
fileout = "L2-008_germanium_data.txt"
peak_counts_correlation = "peak_counts_correlation.txt"
deadtimefileout = "L2-008_germanium_dead_time.txt"
chi_fileout = "chi_squared_data.txt"
expected_fileout = "expected_data.txt"
toy_fileout = "toy_data.txt"
chi_fileout_toy = "toy_chi_squared_data.txt"
likelihood_fileout = "likelihood_data.txt"
peak_fileout = "peak_data.txt"
FWHM_fileout = "FWHM_data.txt"
flat_fileout = "flat_data.txt"
real_fileout = "real_data.txt"
pressure_fileout = "pressure_data.txt"
away_from_mean_fileout = "away_from_mean.txt"
away_error_fileout = "away_error_fileout.txt"
lomb_fileout = "lomb_fileout.txt"
try:
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + away_error_fileout)	
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + away_from_mean_fileout)
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout)
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + deadtimefileout)
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + chi_fileout)
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + expected_fileout)
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + toy_fileout)
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + chi_fileout_toy)
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + likelihood_fileout)
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + peak_fileout)
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + FWHM_fileout)
    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + peak_counts_correlation)
except OSError:
	pass

##########################################################################################
### Real Counts
##########################################################################################
number = 0
counter=0
error117 = 0
error133 = 0
last_day = 0
list_117 = []
list_133 = []
list_time =[]
list_deat_time = []
list_error117 = []
list_error133 = []
live_time = 0

fileout = "L2-008_germanium_data.txt"
peak_counts_correlation = "peak_counts_correlation.txt"
deadtimefileout = "L2-008_germanium_dead_time.txt"
chi_fileout = "chi_squared_data.txt"
expected_fileout = "expected_data.txt"
toy_fileout = "toy_data.txt"
chi_fileout_toy = "toy_chi_squared_data.txt"
likelihood_fileout = "likelihood_data.txt"
peak_fileout = "peak_data.txt"
FWHM_fileout = "FWHM_data.txt"
sample_fileout = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/CO60_800LIVE_000001.Spe"
flat_fileout = "flat_data.txt"
#real_fileout = "ave_out_file"
pressure_fileout = "pressure_data.txt"
leastsquares_fileout = "leastsquares.txt"
temperature_file = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/Temperature.txt"
temperature_out_file = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/detector_temperature.txt"
humidity_out_file = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/detector_humitidy.txt"
x_out_file = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/x.txt"
y_out_file = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/y.txt"

if REAL == True:
	try:
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + real_fileout)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + deadtimefileout)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + peak_fileout)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + FWHM_fileout)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + peak_counts_correlation)
		os.remove(x_out_file)
		os.remove(y_out_file)
	except OSError:
		pass
	
	for file in glob.glob('/Users/spenceraxani/Documents/Nuclear_Decay/Data//data/2013_Aug_17_L2-008_germanium/*.Spe'):
		h1 = TH1F('^{60}CO Spectrum','^{60}CO Spectrum',16384,0,16384)
		h1.Draw()
		fh = open(file)
		lines = fh.readlines()
		number = number +1
		net117 = 0
		net133 = 0
		deadtime = 0
		counter = counter + 1
		line_number  = 0 
		bin_number = 0
		total_counts_117 =0
		total_counts_133 = 0
		start_of_file = 12
		end_of_file = 16396
		total_counts = 0
		for line in lines:
			line_number += 1
			if line_number == 10:
				columns = line.split(' ')
				live_time = float(columns[0])
				print(live_time)
				deadtime = 100*(float(columns[1])-float(columns[0]))/float(columns[1])
				print(deadtime)
			if line_number == 8:
				columns = line.split(' ')
				print(columns)
				date = columns[0]
				clock = columns[1]
				(hour, min, sec) = clock.split(':')
				(month, day, year) = date.split('/')
				actual_time  = float((int(hour) * 3600 + int(min) * 60 + int(sec)))/86400
				year = int(year)
				month = int(month)
				day = int(day)
				t = time.mktime((year, month, day, 0, 0, 0, 0, 0, 0))
				time.gmtime(t) 
				if year == 2013:
					proper_time= float(time.gmtime(t)[7]) + float(actual_time) + 56293.50
				elif year == 2014:
					proper_time= float(time.gmtime(t)[7]) + float(actual_time) + 56658.500000
				elif year == 2015:
					proper_time= float(time.gmtime(t)[7]) + float(actual_time) + 56658.500000 +365
				elif year == 2015:
					proper_time= float(time.gmtime(t)[7]) + float(actual_time) + 56658.500000 + 2 * 365
				print(proper_time)
				if counter == 1:
					first_day = proper_time
				if last_day < proper_time:
					last_day = proper_time	

			if line_number > start_of_file and line_number < end_of_file:
				line = ' '.join(line.split())
				bin_number += 1
				h1.AddBinContent(bin_number,int(line))

		s = TSpectrum(2)
		found = s.Search(h1,1,"new")
		peaks = s.GetPositionX()
		height = s.GetPositionY()
		print("Found candidate peaks = " + str(peaks[0]) + 'bin with ' +str(height[0]) +' counts, and '+ str(peaks[1]) + 'bin with ' +str(height[1])+ ' counts')
		bin = 0
		low117 = 0
		high117 = 0 
		binlow117 = 0
		binhigh117 = 0
		bin_difference = peaks[1]-peaks[0]
		energy_difference = 1332.492 - 1173.228
		bins_per_energy = bin_difference/energy_difference
		for i in range(100):
			low117 = h1.GetBinContent(int(peaks[0])-i)
			binlow117 = int(peaks[0])-i
			if low117 < h1.GetBinContent(int(peaks[0]))/2:
				break
		for i in range(100):
			high117 = h1.GetBinContent(int(peaks[0])+i)
			binhigh117 = int(peaks[0])+i
			if high117 < h1.GetBinContent(int(peaks[0]))/2:
				break
		for i in range(500):
			excess_counts = h1.Integral(int(peaks[1])+100 , int(peaks[1])+600)
		FWHM_energy = (binhigh117 - binlow117)/bins_per_energy
		FWHM_bin = binhigh117 - binlow117
		print("The FWHM is : "+str(FWHM_energy)+" KeV")
		for i in range(int(peaks[0])-100,int(peaks[0])+100):
			bin_counts_117 = h1.GetBinContent(i)
			total_counts_117 += bin_counts_117
		for j in range(int(peaks[1])-100,int(peaks[1])+100):
			bin_counts_133 = h1.GetBinContent(j)
			total_counts_133 += bin_counts_133
		net_counts = total_counts_117 + total_counts_133		
		expected_counts_poisson = exp_counts(proper_time)
		error_net_counts = math.sqrt(total_counts_117+ total_counts_133)
		real_output = open(real_fileout,'a')
		real_output.write(str(proper_time) + "\t" + str(0) + "\t" + "\t" + str(round_sig(math.sqrt(total_counts_117/live_time))) + "\t"+  "\t"+str(round_sig(total_counts_133/live_time)) + "\t"+"\t"+ str(round_sig(math.sqrt(total_counts_133/live_time))) + "\t"+"\t"+str(round_sig(net_counts/live_time))+"\t"+"\t" + str(round_sig(error_net_counts/live_time))+"\t" + str(expected_counts_poisson)+"\t"+str(round_sig(net_counts))+"\t"+str(round_sig(live_time))+ '\t'+ str(FWHM_energy) +"\t"+str(round_sig(excess_counts))+ "\n")
	 	deadoutput = open(deadtimefileout,'a')
	 	deadoutput.write(str(proper_time) + "\t" + str(deadtime) + "\n" )
	 	x_outfile= open(x_out_file,'a')
	 	x_outfile.write(str(proper_time) +"\n" )
	 	y_outfile= open(y_out_file,'a')
	 	y_outfile.write(str(round_sig(net_counts/live_time))+ "\n" )
	 	peak_output = open(peak_fileout,'a')
	 	peak_output.write(str(proper_time) + "\t" + str(peaks[0]) + '\t' + str(peaks[1]) +"\n")
	 	fh.close()	 	

juliandate , zeroes , er117, c133, er133, net, ernet, exp, realcount, livetime, fwhm, excesscounts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/'+real_fileout, unpack=True)

first_day = 0
last_day = 0
last_count = 0
first_count = 0
number_of_files = 0
total_counts_ever = 0
for k in range(len(juliandate)):
	#exp_gr.SetPoint(k,juliandate[k],exp_counts(juliandate[k]))
	first_day = juliandate[0]
	total_counts_ever += realcount[k]
	first_count = net[0]
	number_of_files += 1
	if juliandate[k] > last_day:
		last_day = juliandate[k]
		last_count = net[k]
print("The experiment started on: " + str(first_day))
print("The experiment ended on: " + str(last_day))
print("The experiment ran for: " + str(last_day- first_day))
print("And we took : " + str(number_of_files) + " spectra")
print("And we saw a total number of decays : " + "%.4g" % total_counts_ever)
c1 = TCanvas('c1', 'Peak Information',600,1800)
c1.Draw()
c1.cd()
p11 = TPad('p11','p',0.0,0.5,1,1)
#p11.SetGrid()
p11.Draw()
p12 = TPad('p12','p',0.0,0.0,1,0.5)
p12.Draw()
p11.cd()


##########################################################################################
### Distance to sun
##########################################################################################	
sun_file = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/solar_distance.txt"

if SUN_DISTANCE == True:
	try:
		os.remove(sun_file)
	except OSError:
		pass
	
	date = numpy.linspace(first_day-100, last_day+100, 20000)
	for i in range(len(date)):
		sunfile = open(sun_file,'a')
	 	sunfile.write(str(date[i]) + "\t" + str((solar_distance(date[i]) - 147166462000.0)*(first_count - last_count)/10000506000.0+980) + "\n" )
date_2 , height = numpy.loadtxt("/Users/spenceraxani/Documents/Nuclear_Decay/Data/solar_distance.txt", unpack=True)
date2 = array("d",date_2)
height2 = array("d",height)

solar_dis = TGraph(len(date2), date2, height2)
solar_dis.GetXaxis().SetLimits(first_day,last_day);
solar_dis.SetLineColor(2)

##########################################################################################
### end
##########################################################################################	

xreal = array("d",juliandate)
yreal = array("d",net)
errorreal = array("d",ernet)
zeros = array("d",zeroes)

real_gr = TGraph(len(xreal), xreal, yreal)
real_gr.GetHistogram().SetMaximum(1050)
real_gr.GetHistogram().SetMinimum(980)
real_gr.GetXaxis().SetTitle("Modified Julian Date MST [Days]")
real_gr.GetYaxis().SetTitle("Count Rate [cps]")
real_gr.SetTitle("Cobalt 60 decay count rate")
real_gr.GetXaxis().SetLimits(first_day,last_day);
real_gr.GetYaxis().SetLimits(980,1050);
p11.Update()
#exp_gr.SetDrawOption("l")
real_gr.SetMarkerColor(1)
real_gr.Draw("ap")

fsignal = TF1("FSignal","expo", first_day ,last_day)
#real_gr.Fit("fsignal","","",first_day ,last_day)
fsignal.SetLineColor(2)
real_gr.Fit(fsignal,'R')
print("            " +str(fsignal.GetParameter(1)))
print("So the actual Chi2 Fit give a half life of " + str(-numpy.log(2)/fsignal.GetParameter(1))+ " +/- "+str((fsignal.GetParError(1)/fsignal.GetParameter(1))*numpy.log(2)/fsignal.GetParameter(1))+ "  "+ str(fsignal.GetParameter(0)))
B = fsignal.GetParameter(1)
A = math.exp(fsignal.GetParameter(0))
fit_dev = TH1F('Deviation from Mean','Deviation from mean',200,-5,5)
for k in range(len(juliandate)):
	diff = fit_counts(fsignal.GetParameter(1),math.exp(fsignal.GetParameter(0)),juliandate[k]) - net[k]
	fit_dev.Fill(diff,1.0/len(juliandate))
fit_dev.GetYaxis().SetTitleOffset(10.4)
fit_dev.GetXaxis().SetTitleOffset(10.4)
fit_dev.SetLineColor(4)
fit_dev.GetYaxis().SetTitle("Probability")
fit_dev.GetXaxis().SetTitle("Deviation from mean [cps]")
fit_dev.Fit("gaus",'','', -3, 3 )

#fit_gr.Draw("ELsame")
solar_dis.Draw("ELsame")

#exp_gr.Draw("Elsame")
p12.cd()
juldate , peak1 , peak2 = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/'+ peak_fileout, unpack=True)
x = array("d",juldate)
y = array("d",peak1)
z = array("d",peak2)
peak1_graph = TGraph(len(x), x, y)
peak1_graph.GetXaxis().SetLimits(first_day,last_day);

peak2_graph = TGraph(len(x), x, z)
peak1_graph.SetLineColor(9)
peak2_graph.SetLineColor(3)
multi_peak = TMultiGraph()
multi_peak.SetTitle("Peak Location;Modified Julian Date MST [days];Bin Number")
multi_peak.Add(peak1_graph)
multi_peak.Draw("AP")
multi_peak.GetYaxis().SetTitleOffset(1.4);
multi_peak.GetXaxis().SetLimits(first_day,last_day)
p12.Modified();
p12.Update()

c2 = TCanvas('c2', 'Deadtime',600,1800)
c2.Draw()
c2.cd()
p21 = TPad('p21','p',0.0,0.5,1,1)
p21.SetGrid()
p21.Draw()
p22 = TPad('p22','p',0.0,0.0,1,0.5)
p22.Draw()
p21.cd()
real_gr.Draw("ap")
#fit_gr.Draw("Elsame")
real_gr.GetXaxis().SetLimits(first_day,last_day);
p21.Update()

print("Calculating Deadtime ...")
p22.cd()
juldate , deadtime = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/'+ deadtimefileout, unpack=True)
x = array("d",juldate)
y = array("d",deadtime)
deadtime_graph = TGraph(len(x), x, y)
deadtime_graph.GetXaxis().SetTitle("Modified Julian Date MST [Days]")
deadtime_graph.GetYaxis().SetTitle("Deadtime [%]")
deadtime_graph.SetTitle("Deadtime")
deadtime_graph.GetXaxis().SetLimits(first_day,last_day);
deadtime_graph.Draw("AP")
p22.Update()

######################################
error_for2 = "error_for.txt"
if LINE == True:
    try:
        os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + error_for2)	
    except OSError:
    	pass
    for u in range(10000):
	    error_for = open(error_for2,'a')
	    error_for.write(str(56565+0.002*u) + "\t" +str(fit_counts(B,A,56565+0.002*u))+ "\t" +str(0)+ "\t"+ str(math.sqrt(fit_counts(B,A,56565+0.002*u)*800)/800)+ "\n" )


##########################################################################################
### FWHM and baseline
##########################################################################################	

cfwhm = TCanvas('c1cfwhm', 'FWHM and Baseline',600,1800)
cfwhm.Draw()
cfwhm.cd()
p1fwhm = TPad('p1fwhm','p',0.0,0.5,1,1)
p1fwhm.SetGrid()
p1fwhm.Draw()
p1base = TPad('p1base','pa',0.0,0.0,1,0.5)
p1base.SetGrid()
p1base.Draw()
p1base.cd()

juliandate , zeroes , er117, c133, er133, net, ernet, exp, realcount, livetime, fwhm, excesscounts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/'+real_fileout, unpack=True)

x = array("d",juliandate)
y1 = array("d",fwhm)
y2 = array("d",excesscounts)

fwhm_graph = TGraph(len(x), x, y1)
fwhm_graph.SetMarkerStyle(7)
fwhm_graph.GetXaxis().SetTitle("Modified Julian Date MST [Days]")
fwhm_graph.GetYaxis().SetTitle("Full Width Half Maximum [keV]")
fwhm_graph.SetTitle("FWHM")
fwhm_graph.GetXaxis().SetLimits(first_day,last_day);
fwhm_graph.Draw("AP")
p1base.Update()
p1fwhm.cd()
baseline_graph = TGraph(len(x), x, y2)
baseline_graph.GetXaxis().SetTitle("Modified Julian Date MST [Days]")
baseline_graph.GetYaxis().SetTitle("Sum counts above 1.33 MeV Peak")
baseline_graph.SetTitle("Baseline")
baseline_graph.GetXaxis().SetLimits(first_day,last_day);
baseline_graph.Draw("AP")
p1fwhm.Update()
##########################################################################################
### Temperature
##########################################################################################	
print("Calculating Temperature and Humidity ...")

if TEMP_DATA == True:
	try:
		os.remove(temperature_out_file)
		os.remove(humidity_out_file)
	except OSError:
		pass
	temp = open(temperature_file, 'r')
	counter1 = 0
	for line in temp:
		counter1 += 1;
		if counter1 == 89:
			counter1 = 0
			column = line.split('\t')
			columns = [col.strip() for col in column]
			#print(columns)
			if columns[-1] == "9":
				if not ":" in columns[-3]: 
					temperature = columns[2]
					humidity = columns[-2]
				clock1 = columns[0]
				(day, month, year) = clock1.split('/')
				year = int(year)
				month = int(month)
				day = int(day)
				clock2 = columns[1]
				#print(clock2)
				(watch, noon) = clock2.split(' ')
				(hours, mins) = watch.split(":")
			
				if noon == "AM":
					#print(hours)
					if int(hours) == int(12):
						hour = 0
						#print(hour)
					else:
						hour = int(hours)
					#print(hour)
				else:
					if int(hours) == int(12):
						hour = 12
						#print(hour)
					else:
						hour = int(hours) + 12
					#print(hour)
				actual_time  = float((int(hour) * 3600 + int(mins)*60))/86400
				current_time = time.mktime((year, month, day, 0, 0, 0, 0, 0, 0))
				a =time.gmtime(current_time) 
				proper_time= float(time.gmtime(current_time)[7]) + float(actual_time)
				if year == 2013:
					temp_out = open(temperature_out_file, 'a')
					temp_out.write(str(56293.50 + proper_time) + "\t" + str(temperature) + "\n")
					hum_out = open(humidity_out_file, 'a')
					hum_out.write(str(56293.50 + proper_time) + "\t" + str(humidity) + "\n")
				if year == 2014:
					temp_out = open(temperature_out_file, 'a')
					temp_out.write(str(56293.50 + 365+ proper_time) + "\t" + str(temperature) + "\n")
					hum_out = open(humidity_out_file, 'a')
					hum_out.write(str(56293.50 +365+ proper_time) + "\t" + str(humidity) + "\n")

t1 = TCanvas('temperature', 'Temperature',600,1800)
t1.Draw()
t1.cd()
t11 = TPad('t11','p',0.0,0.5,1,1)
t11.Draw()
t11.SetGrid()
t12 = TPad('t12','p',0.0,0.0,1,0.5)
t12.Draw()
t12.SetGrid()

temp_date  , temp_det = numpy.loadtxt(temperature_out_file, unpack=True)
julian_date  , hum_det 	= numpy.loadtxt(humidity_out_file, unpack=True)

x_temp_date 	= array("d",temp_date)
y_temp 	= array("d",temp_det)
x_date 	= array("d",julian_date)
y_hum 	= array("d",hum_det)
gedet_temp 	= TGraph(len(x_temp_date), x_temp_date, y_temp)
gedet_temp.SetMarkerColor(2)
gedet_temp.SetLineColor(2)
gedet_hum 	= TGraph(len(x_date), x_date, y_hum)
gedet_hum.SetMarkerColor(1)
gedet_hum .SetLineColor(1)

t11.cd()
mg_temp = TMultiGraph()
mg_temp.SetTitle("Aluminum cap temperature;Modified Julian Date MST [days]; Temperature [#circC]")
mg_temp.Add(gedet_temp)
mg_temp.Draw("AL")
if SYNC == True:
	mg_temp.GetXaxis().SetLimits(first_day,last_day); 
t11.Modified();
t11.Update()
t12.cd()
real_gr.Draw("Ap")
#fit_gr.Draw("Elsame")
real_gr.GetXaxis().SetLimits(first_day,last_day);
t12.Update()

h1 = TCanvas('h11', 'humidity',600,1800)
h1.Draw()
h1.cd()
h11 = TPad('p301','p',0.0,0.5,1,1)
h11.SetGrid()
h11.Draw()
h12 = TPad('p311','p',0.0,0.0,1,0.5)
h12.SetGrid()
h12.Draw()

h11.cd()
mg_hum = TMultiGraph()
mg_hum.SetTitle("L2-008 relative humidity;Modified Julian Date MST [days]; Relative Humidity")
mg_hum.Add(gedet_hum)
mg_hum.Draw("AL")
if SYNC == True:
	mg_hum.GetXaxis().SetLimits(first_day,last_day); 
h11.Modified();
h11.Update()

h12.cd()
real_gr.Draw("ap")
#fit_gr.Draw("Elsame")
real_gr.GetXaxis().SetLimits(first_day,last_day)
h12.Update()

##########################################################################################
### Temperature of cap at every data point
##########################################################################################
modifieddata = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/modified_with_temp'
if MODTEMP == True:
	try:
	    os.remove(modifieddata)
	except OSError:
	    pass
	the_time, the_temp = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/detector_temperature.txt', unpack=True)
	time = array("d",the_time)
	temp = array("d",the_temp)

	the_julian, non, err11, errs, ewrwer,  coustsps, err_countsps, exp_countsps, total_counts_seen, live_times = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/real_data.txt', unpack=True)
	time_counts = array("d",the_julian)
	counts_ps =  array("d",coustsps)
	err_counts_ps = array("d",err_countsps)
	tot_counts = array("d",total_counts_seen)

	temperature_at_point =[]
	das_date = []
	temp_at_data = 0
	modified_data = open(modifieddata,'a')
	for i in range(len(time_counts)):
		for j in range(len(time)):
			if time_counts[i] > time[j]:
				temps = temp[j]	
				dasdate = time[j]
		temperature_at_point.append(temps)
		das_date.append(dasdate)
				#print(temperature_at_point)
		modified_data.write(str(time_counts[i]) + "\t" + str(counts_ps[i]) + "\t" + str(tot_counts[i]) + "\t"+ str(err_counts_ps[i]) + "\t" + str(temps) + "\t" + str(tot_counts[i] + (temps - 15)*(22.2*math.pow(10,-6))*0.101*100000) + "\n" )

##########################################################################################
### TEMP versus HUM
##########################################################################################
c6 = TCanvas('c6', 'Standard Deviation',600,1800)
c6.Draw()
c6.cd()
p61 = TPad('p61','p',0.0,0.5,1,1)
p61.SetGrid()
p61.Draw()
p62 = TPad('p62','p',0.0,0.0,1,0.5)
p62.Draw()
p61.cd()
real_gr.Draw("ap")
real_gr.GetXaxis().SetLimits(first_day,last_day);
p61.Update()

p62.cd()
temp_date  , temp_det = numpy.loadtxt(temperature_out_file, unpack=True)
julian_date  , hum_det 	= numpy.loadtxt(humidity_out_file, unpack=True)
y_temp 	= array("d",temp_det)
x_hum 	= array("d",hum_det)
temp_hum = TGraph(len(x_hum), x_hum, y_temp)
temp_hum.GetXaxis().SetTitle("Relativity Humidity")
temp_hum.GetYaxis().SetTitle("Temperature [#circC]")
temp_hum.SetTitle("Temperature Humidity Correspondance")

temp_hum.Draw("ap")
p62.Update()

##########################################################################################
### Averaging over multiple data points
##########################################################################################
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
juliandate , c117 , zeros, c133, er133, net, ernet, exp , realcount, livetime, fwhm, excesscounts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/'+real_fileout, unpack=True)
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


mg_away.Add(plot_awaymax,"p")
mg_away.Add(plot_awaymin,"p")
mg_away.Add(plot_away,"lp")
mg_away.Draw("a")

mg_away.GetYaxis().SetTitleOffset(1.5);
mg_away.GetXaxis().SetLimits(first_day,last_day);

p2ave.Modified();
p2ave.Update()

##########################################################################################
### Histogram dates and such
##########################################################################################	
print('Calculating Daily Histogram...')
juliandate , c117 , zeros, c133, er133, net, ernet, exp , realcount, livetime, fwhm, excesscounts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/'+real_fileout, unpack=True)
hist_file = 'hist_fileout.txt'


if HIST == True:
	try:
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + hist_file)
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
		if num_files == 1:
			for i in range(len(list_net)):
				ave_sum += list_net[i]
				ave_date += list_date[i]
			ave_error = math.sqrt(ave_sum)
			ave_error = ave_error/num_files
			ave_sum = ave_sum/num_files
			ave_date = ave_date/num_files
			histoutfile = open(hist_file,'a')
			histoutfile.write(str(ave_date) + "\t" + str(round_sig(ave_sum/livetime[k]))+ "\t" + str(round_sig(ave_error/livetime[k]))+"\n")
			num_files = 0
			ave_sum = 0
			ave_date = 0
			list_net = []
			list_date = []	
mean_hist_fileout = "hist_mean_fileout.txt"
if MEAN:
	try:
	    os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + mean_hist_fileout)
	except OSError:
		pass
	juliandate, counts, error = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/hist_fileout.txt', unpack=True)
	for k in range(len(juliandate)):
		diff = 1 - (fit_counts(B,A,juliandate[k]) - counts[k])/fit_counts(B,A,juliandate[k])
		errdifftop = 1 + error[k]/fit_counts(B,A,juliandate[k])
		errdiffbottom = 1 - error[k]/fit_counts(B,A,juliandate[k])
		mean_histfileout = open(mean_hist_fileout,'a')
		mean_histfileout.write(str(juliandate[k]) + "\t" + str(diff) + "\t" + str(1-errdiffbottom) + "\n" )

chist = TCanvas('chist', 'Daily Cycle',600,1800)
chist.Draw()
chist.cd()
p1hist = TPad('p1hist','p',0.0,0.5,1,1)
p1hist.SetGrid()
p1hist.Draw()
p2hist = TPad('p2hist','p',0.0,0.0,1,0.5)
p2hist.SetGrid()
p2hist.Draw()
p1hist.cd()

juliandate , away, errors = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/hist_mean_fileout.txt', unpack=True)
x = array("d",juliandate)
y = array("d",away)
plot_hist_away = TGraph(len(x), x, y)
plot_hist_away.Draw("ap")
plot_hist_away.SetMarkerStyle(6)
plot_hist_away.GetXaxis().SetTitle("Modified Julian Date [days]")
plot_hist_away.GetYaxis().SetTitle("Count Rate [cps]")
plot_hist_away.GetYaxis().SetTitleOffset(1.3)
plot_hist_away.SetTitle("")

plot_hist_away.GetXaxis().SetLimits(first_day,last_day)
p1hist.Modified();
p1hist.Update()

juliandate , away, errors = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/hist_fileout.txt', unpack=True)
x = array("d",juliandate)
y = array("d",away)

histogram = "/Users/spenceraxani/Documents/Nuclear_Decay/Data/daily_counts.txt"
try:
	os.remove(histogram)
except OSError:
	pass

p2hist.cd()
daily = TH1F("hist","histo",24,0,24)
#daily.Sumw2()

separation = 1000/24
print(separation)
total_sum = 0
count_sum = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
ave_count_sum = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
counting = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
hour = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
for i in range(len(juliandate)):
	total_sum += away[i]
	t = float(juliandate[i])
	number = "%.3f" % round(t,4)
	number = str(number)
	columns = number.split('.')
	for j in range(1,25):
		if int(columns[1]) < j * separation and int(columns[1]) > (j-1) * separation:
			counting[j-1]+= 1
			#print(columns[1])
			count_sum[j-1] += away[i]

for j in range(len(count_sum)):
	ave_count_sum[j] = count_sum[j]/counting[j] * 620
	#print(str(hour[j]) + "\t" + str(count_sum[j]))
	daily.Fill(hour[j],ave_count_sum[j])
daily.GetXaxis().SetTitle("Time of Day [hours]")
daily.GetYaxis().SetTitle("Sum of Count Rate [cps]")
daily.SetTitle("Hourly Time Variations")
daily.Draw()
daily.Draw('Esame')
daily.SetStats(False)
daily.SetMaximum(629000)
daily.SetMinimum(625000)

p2hist.Modified();
p2hist.Update()

##########################################################################################
### Lomb Scargle
##########################################################################################	
print("Calculating Lomb Scargle Spectra...")
clomb = TCanvas('clomb', 'Lomb Scargle Discrete Transform',600,1800)
clomb.Draw()
clomb.cd()
p1lomb = TPad('p1lomb','pl',0.0,0.5,1,1)
#p1lomb.SetGrid()
p1lomb.Draw()
p2lomb = TPad('p2lomb','pla',0.0,0.0,1,0.5)
p2lomb.Draw()
p2lomb.SetGrid()
p2lomb.cd()

mg_away.Draw("a")
p2lomb.Modified();
p2lomb.Update()
p1lomb.cd()

if LOMB == True:
	try:
		os.remove(lomb_fileout)
	except OSError:
		pass
	time , counts= numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/data_detrended.txt', unpack=True)
	time_series = array("d",time)
	mean = np.mean(counts)
	time_list = []
	sample_x = []
	for i in range(len(time_series)):
		time_list.append(time_series[i]-time_series[0])
		sample_x.append((counts[i]-mean)*1000)
	x = np.array(time_list, np.float64)
	measurement = array("d",counts)
	y = np.array(measurement, np.float64)
	N = len(sample_x)*1.0
	sumh = 0 
	sigma2 = 0 
	for j in range(len(sample_x)-1):
		sumh += sample_x[j]
	hbar = (1/N)*sumh
	hhbar2 = 0
	for k in range(len(sample_x)-1):
		hhbar2 += (sample_x[k] - hbar)**2
	sigma2 = (1/(N-1))*hhbar2

	ang_freq_domain = []
	periodogram = []
	ang_freq = 1
	sum_sin_2wtau = 0
	sum_cos_2wtau = 0
	sum_x_cos_wttau = 0
	sum_cos2_wttau = 0
	sum_x_sin_wttau = 0
	sum_sin2_wttau = 0

	i0=drange(0.005, 2.005, 0.005)
	wnum = 0.0
	for i in i0:
		wnum += 1.0
	po=0.01
	z = - np.log(1-np.power(1-po,1/wnum))

	print("the number of frequencies: " +str(wnum))
	print('the z value is :' + str(z))
	print("the variance is: " +str(sigma2))
	print("the number of points is N: " +str(N))
	print("the mean is: " +str(hbar)+ ' ' +str(np.mean(sample_x)))

	i0=drange(0.005, 2.00, 0.005)
	for w in i0:
		print(w)
		for t in time_list:
			sum_sin_2wtau += np.sin(2.0*w*t)
			sum_cos_2wtau += np.cos(2.0*w*t)
		tau = np.arctan(sum_sin_2wtau/sum_cos_2wtau)/(2.0*w)

		for t in range(len(time_list)):
			sum_x_cos_wttau +=  (sample_x[t]-hbar)*np.cos(w*(time_list[t]-tau))
			sum_x_sin_wttau +=  (sample_x[t]-hbar)*np.sin(w*(time_list[t]-tau))
			sum_sin2_wttau += np.sin(w*(time_list[t]-tau))**2
			sum_cos2_wttau += np.cos(w*(time_list[t]-tau))**2
		P = (0.5/sigma2)*(sum_x_cos_wttau**2)/sum_cos2_wttau + (0.5/sigma2)*(sum_x_sin_wttau**2)/sum_sin2_wttau
		periodogram.append(P)
		ang_freq_domain.append(w)
		sum_sin_2wtau = 0
		sum_cos_2wtau = 0
		sum_x_cos_wttau = 0
		sum_cos2_wttau = 0
		sum_x_sin_wttau = 0
		sum_sin2_wttau = 0

	for i in range(len(periodogram)):
		lomboutput = open(lomb_fileout,'a')
		lomboutput.write(str(ang_freq_domain[i])+ "\t" + str(periodogram[i])  + "\n" )
p1sigma = 0.317311
p2sigma = 0.045501
p3sigma = 0.0027
p5sigma = 0.0000006004

z1 = - np.log(1-np.power(1-p1sigma,1/200.0))
z2 = - np.log(1-np.power(1-p2sigma,1/200.0))
z3 = - np.log(1-np.power(1-p3sigma,1/200.0))
z5 = - np.log(1-np.power(1-p5sigma,1/200.0))
print("the most obvious period is: " +str(2*3.14159265/0.164147))
print("the second obvious period is: " +str(2*3.14159265/0.0527322))
print("the third obvious period is: " +str(2*3.14159265/0.221163))
p1 = TLine(0,z1,2,z1);
p1.SetLineStyle(2);
p1.SetLineColor(40)
p1.SetLineWidth(3)
p2 = TLine(0,z2,2,z2);
p2.SetLineStyle(2);
p3 = TLine(0,z3,2,z3);
p3.SetLineStyle(2);
p3.SetLineColor(30)
p3.SetLineWidth(3)
p5 = TLine(0,z5,2,z5);
p5.SetLineStyle(2);
p5.SetLineColor(45)
p5.SetLineWidth(3)

oneweek = TLine(2*math.pi/7,0,2.0*math.pi/7,40);
oneweek.SetLineStyle(2);
oneweek.SetLineColor(1)
oneweek.SetLineWidth(3)

onemonth = TLine(2*math.pi/30.4,0,2*math.pi/30.4,40);
onemonth.SetLineStyle(2);
onemonth.SetLineColor(1)
onemonth.SetLineWidth(3)

oneyear = TLine(2*math.pi/365,0,2*math.pi/365.4,40);
oneyear.SetLineStyle(2);
oneyear.SetLineColor(1)
oneyear.SetLineWidth(3)

freq, amp  = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/lomb_fileout.txt', unpack=True)
x1 = array("d",freq)
y1 = array("d",amp)

lombs_gr = TGraph(len(x1), x1, y1)
lombs_gr.GetXaxis().SetTitle("Angular Frequency [rads/day]")
lombs_gr.GetYaxis().SetTitle("Arbitrary Units")
lombs_gr.SetTitle("Lomb Scargle of Averaged data points")
lombs_gr.GetXaxis().SetLimits(0,2);
lombs_gr.GetHistogram().SetMaximum(40)          
lombs_gr.GetHistogram().SetMinimum(0)
lombs_gr.Draw('al')

p1.Draw('Same')
p3.Draw('Same')
p5.Draw('Same')

p = TLine()

leg = TLegend(0.65,0.65,0.85,0.85)
leg.SetFillColor(0)
leg.AddEntry(p,"Lomb Scargle Periodogram","l")
leg.AddEntry(p1,"1 #sigma","l")
leg.AddEntry(p3,"3 #sigma","l")
leg.AddEntry(p5,"5 #sigma","l")
leg.Draw()
#oneweek.Draw('same')
#onemonth.Draw('same')
#oneyear.Draw('same')
p1lomb.Update()
p1lomb.Update()

##########################################################################################
### Lomb Scargle
##########################################################################################	
ang_freq = 200
trials = 10
MC_fileout = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/MC_data.txt'
if MC == True:
	try:
	    os.remove(MC_fileout)	
	except OSError:
		pass
	time , counts, errs= numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/hist_mean_fileout.txt', unpack=True)
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
		#print(counter)
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
z5 = 0.0000006004
v5a = - np.log(1-np.power(1-z5,1/200.0))

pscar1.cd()
pscar1.SetLogy()
h1hist = TH1F('Maximums','scargy',bins,0,max_amp)
h1hist.SetTitle("Maximums;Arbitrary Units; Angular Frequency [rads/day]")
h1hist.Draw()
spectrum = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/MC_data.txt', unpack=True)
for i in range(len(spectrum)):
	h1hist.Fill(spectrum[i],1)
print("The maximum peak out of the " +str(trials)+" trials was " +str(np.amax(spectrum)))

bins_to_amp = 1.0*bins/max_amp

total_bins = 0
for j in range(bins):
	total_bins += h1hist.GetBinContent(j)
	if total_bins < trials * (1-z1):
		v1 = j/bins_to_amp
		height1 = h1hist.GetBinContent(j)
	if total_bins < trials * (1-z3):
		v3 = j/bins_to_amp
		height3 = h1hist.GetBinContent(j)
	if total_bins < trials * (1-z5):
		v5 = j/bins_to_amp
		height5 = h1hist.GetBinContent(j)

p1line = TLine(v1,0,v1,height1)
p1line.SetLineStyle(2)
p1line.SetLineColor(40)
p1line.SetLineWidth(4)
p3line = TLine(v3,0,v3,height3)
p3line.SetLineStyle(2);
p3line.SetLineColor(30)
p3line.SetLineWidth(4)
p5line = TLine(v5,0,v5,height5);
p5line.SetLineStyle(2);
p5line.SetLineColor(45)
p5line.SetLineWidth(4)
p5aline = TLine(v5a,0,v5a,height5);
p5aline.SetLineStyle(2);
p5aline.SetLineColor(1)
p5aline.SetLineWidth(4)

p1line.Draw("same")
p3line.Draw("same")
p5line.Draw("same")
pscar1.Modified();
pscar1.Update()

pscar2.cd()

#time , counts, errs= numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/hist_mean_fileout.txt', unpack=True)
time , counts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/binned_counts_copy.txt', unpack=True)

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

h2hist = TH1F('Scargle Amplitude1','scargy1',200,0,1)
h2hist.SetTitle("Scargle Amplitude; Angular Frequency [rads/day];Arbitrary Units")

for i in range(len(lombs)):
	h2hist.Fill(f[i],lombs[i])	

s = TSpectrum(4,1)
nfound = s.Search(h2hist,1,"new")
peaks = s.GetPositionX()
height = s.GetPositionY()
print("Found candidate peaks = " + str(2*np.pi/peaks[0]) + ' '+ str(2*np.pi/peaks[1])+ ' '+ str(2*np.pi/peaks[2])+ ' '+ str(2*np.pi/peaks[3]))

h2hist.Draw('l')
p1sigma = 0.317311
p3sigma = 0.0027
p5sigma = 0.0000006004

z1 = - np.log(1-np.power(1-p1sigma,1/200.0))
z3 = - np.log(1-np.power(1-p3sigma,1/200.0))
z5 = - np.log(1-np.power(1-p5sigma,1/200.0))

p1hist = TLine(0,z1,1,z1);
p1hist.SetLineStyle(2);
p1hist.SetLineColor(40)
p1hist.SetLineWidth(4)
p3hist = TLine(0,z3,1,z3);
p3hist.SetLineStyle(2);
p3hist.SetLineColor(30)
p3hist.SetLineWidth(4)
p5hist = TLine(0,z5,1,z5);
p5hist.SetLineStyle(2);
p5hist.SetLineColor(45)
p5hist.SetLineWidth(4)
h2hist.SetStats(False)

p1hist.Draw('Same')
p3hist.Draw('Same')
p5hist.Draw('Same')

leghist = TLegend(0.65,0.65,0.85,0.85)
leghist.SetFillColor(0)
leghist.AddEntry(p,"Lomb Scargle Periodogram","l")
leghist.AddEntry(p1,"1 #sigma","l")
leghist.AddEntry(p3,"3 #sigma","l")
leghist.AddEntry(p5,"5 #sigma","l")
leghist.Draw()

pscar2.Modified();
pscar2.Update()



##########################################################################################
### Pressure
##########################################################################################
cpressure = TCanvas('cpressure', 'Pressure',600,1800)
cpressure.Draw()
cpressure.cd()
pp1 = TPad('pp1','p',0.0,0.5,1,1)
pp1.SetGrid()
pp1.Draw()
pp2 = TPad('pp2','p',0.0,0.0,1,0.5)
pp2.Draw()
pp2.SetGrid()
pp2.cd()
print("Calculating Pressure ...")
if PRESSURE == True:
	try:
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + pressure_fileout)
	except OSError:
		pass
	input = '/Users/spenceraxani/Documents/Nuclear_Decay/Data/pressure.txt'
	#h4 = TH1F('Pressure','Pressure',16384,0,16384)
	#h4.Draw()
	fh = open(input)
	lines = fh.readlines()
	for line in lines:
		columns = line.split('-')
		columns = [col.strip() for col in columns]
		if columns[0] == "2013":
			year = columns[0]
			month = columns[1]
			col = columns[2].split('\t')
			day = col[0]
			pressure = col[-1]
			times = col[1].split(':')
			hour = times[0]
			#print(year + " " + month + " " + day)
			min = times[1]
			sec = times[2]
			actual_time  = float((int(hour) * 3600 + int(min) * 60 + int(sec)))/86400
			year = int(year)
			month = int(month)
			day = int(day)
			t = time.mktime((year, month, day, 0, 0, 0, 0, 0, 0))
			time.gmtime(t) 
			proper_time= float(time.gmtime(t)[7]) + float(actual_time) +56293.50
			#print(proper_time)
			pressureoutput = open(pressure_fileout,'a')
		 	pressureoutput.write(str(proper_time) + "\t" + str(pressure) + "\n" )
		if columns[0]=="2014":
			year = columns[0]
			month = columns[1]
			col = columns[2].split('\t')
			day = col[0]
			pressure = col[-1]
			times = col[1].split(':')
			hour = times[0]
			#print(year + " " + month + " " + day)
			min = times[1]
			sec = times[2]
			actual_time  = float((int(hour) * 3600 + int(min) * 60 + int(sec)))/86400
			year = int(year)
			month = int(month)
			day = int(day)
			t = time.mktime((year, month, day, 0, 0, 0, 0, 0, 0))
			time.gmtime(t) 
			proper_time= float(time.gmtime(t)[7]) + float(actual_time) + 56658.500000
			#print(proper_time)
			pressureoutput = open(pressure_fileout,'a')
		 	pressureoutput.write(str(proper_time) + "\t" + str(pressure) + "\n" )
		
x1 , y1 = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + pressure_fileout, unpack=True)
x = array("d", x1)
y = array("d", y1)
pressure_gr = TGraph(len(x), x, y)
pressure_gr.GetXaxis().SetLimits(first_day,last_day);
pressure_gr.GetXaxis().SetTitle("Modified Julian Date MST [Days]")
pressure_gr.GetYaxis().SetTitle("Pressure [kPa]")
pressure_gr.SetTitle("Atmospheric Pressure in L2-004")
pressure_gr.SetLineColor(9)
pressure_gr.Draw("AL")
pp2.Update()

pp1.cd()
real_gr.Draw("ap")
#fit_gr.Draw("Elsame")
real_gr.GetXaxis().SetLimits(first_day,last_day);
pp1.Update()



##########################################################################################
### Temperature pressure and relative humidity
##########################################################################################	
cevery = TCanvas('cevery', 'cevery',600,1800)
cevery.Draw()
cevery.cd()
pp1e = TPad('pp1','p',0.0,0.7,1,1)
pp1e.SetGrid()
pp1e.Draw()
pp2e = TPad('pp2','p',0.0,0.0,1,0.3)
pp2e.Draw()
pp2e.SetGrid()
pp3e = TPad('pp3','p',0.0,0.35,1,0.65)
pp3e.Draw()
pp3e.SetGrid()
pp3e.cd()
pressure_gr.Draw("AL")
pp2e.cd()
mg_temp.Draw("AL")
pp1e.cd()
mg_hum.Draw("AL")

##########################################################################################
### Solar Data Proton and Electron
##########################################################################################	

#ELECTRON
fileout  = "GP_5m_electron_2MeV.txt"
fileout1 = "GP_5m_electron_0.8MeV.txt"
fileout2 = "GP_5m_proton_1MeV.txt"
fileout3 = "GP_5m_proton_5MeV.txt"
fileout4 = "GP_5m_proton_10MeV.txt"
fileout5 = "GP_5m_proton_30MeV.txt"
fileout6 = "GP_5m_proton_50MeV.txt"
fileout7 = "GP_5m_proton_100MeV.txt"
print("Calculating Solar Proton Data ...")
if SOLAR_DATA == True:
	try:
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout1)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout2)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout3)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout4)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout5)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout6)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout7)

	except OSError:
		pass
	for file in glob.glob('/Users/spenceraxani/Documents/Nuclear_Decay/Data/particle/*_Gp_part_5m.txt'):
		list = file
		fh = open(file)
		lines = fh.readlines()
		for line in lines:
			columns = line.split(' ')
			columns = [col.strip() for col in columns]
			if columns[0] == "2013":
				#print(str(columns))
				output = open(fileout,'a')
				if float(columns[-2]) > 0: # pull 2 MeV electron Data
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000 - 0.2916667) + "\t" +  str(columns[-2]) + "\n")
				else:
					continue
				
				output = open(fileout1,'a') #pull 0.8 MeV electron Data
				if float(columns[-4]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000 - 0.2916667) + "\t" +  str(columns[-4]) + "\n")
				else:
					continue
				
				output = open(fileout7,'a') #pull 100 MeV proton Data
				if float(columns[-6]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000 - 0.2916667) + "\t" +  str(columns[-6]) + "\n")
				else:
					continue
					
				output = open(fileout6,'a') #pull 50 MeV proton Data
				if float(columns[-8]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000 - 0.2916667) + "\t" +  str(columns[-8]) + "\n")
				else:
					continue
					
				output = open(fileout5,'a') #pull 30 MeV proton Data
				if float(columns[-10]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000 - 0.2916667) + "\t" +  str(columns[-10]) + "\n")
				else:
					continue
					
				output = open(fileout4,'a') #pull 10 MeV proton Data
				if float(columns[-12]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000 - 0.2916667) + "\t" +  str(columns[-12]) + "\n")
				else:
					continue	
					
				output = open(fileout3,'a') #pull 5 MeV proton Data
				if float(columns[-14]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000 - 0.2916667) + "\t" +  str(columns[-14]) + "\n")
				else:
					continue		
					
				output = open(fileout2,'a') #pull 1 MeV proton Data
				if float(columns[-16]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000 - 0.2916667) + "\t" +  str(columns[-16]) + "\n")
				else:
					continue	
					
			
			if columns[0] == "2014":
				#print(str(columns[-2]) + str(columns[-2]) + str(columns[-2]))
				output = open(fileout,'a')
				if float(columns[-2]) > 0: # pull 2 MeV electron Data
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000- 0.2916667) + "\t" +  str(columns[-2]) + "\n")
				else:
					continue
				
				output = open(fileout1,'a') #pull 0.8 MeV electron Data
				if float(columns[-4]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000- 0.2916667) + "\t" +  str(columns[-4]) + "\n")
				else:
					continue
				
				output = open(fileout7,'a') #pull 100 MeV proton Data
				if float(columns[-6]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000- 0.2916667) + "\t" +  str(columns[-6]) + "\n")
				else:
					continue
					
				output = open(fileout6,'a') #pull 50 MeV proton Data
				if float(columns[-8]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000- 0.2916667) + "\t" +  str(columns[-8]) + "\n")
				else:
					continue
					
				output = open(fileout5,'a') #pull 30 MeV proton Data
				if float(columns[-10]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000- 0.2916667) + "\t" +  str(columns[-10]) + "\n")
				else:
					continue
					
				output = open(fileout4,'a') #pull 10 MeV proton Data
				if float(columns[-12]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000- 0.2916667) + "\t" +  str(columns[-12]) + "\n")
				else:
					continue	
					
				output = open(fileout3,'a') #pull 5 MeV proton Data
				if float(columns[-14]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000- 0.2916667) + "\t" +  str(columns[-14]) + "\n")
				else:
					continue		
					
				output = open(fileout2,'a') #pull 1 MeV proton Data
				if float(columns[-16]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000 - 0.2916667)+ "\t" +  str(columns[-16]) + "\n")
				else:
					continue						
					

Red = TColor(2000,1,0,0)
Redish = TColor(2001,1,.4,0)
Redishish = TColor(2002,1,.8,0)	
Yellow = TColor(2003,0.5,1,0)	
Yellowish = TColor(2004,0,1,1)		
Orange = TColor(2005,0,.5,1)

cproton = TCanvas('cproton', 'Proton Flux',600,1800)
cproton.Draw()
cproton.cd()
pproton1 = TPad('pproton1','p',0.0,0.5,1,1)
pproton1.SetGrid()
pproton1.Draw()
pproton2 = TPad('pproton2','p',0.0,0.0,1,0.5)
pproton2.SetGrid()
pproton2.Draw()
pproton2.cd()
	
solar_date_p100 , solar_flux_p100 = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout7, unpack=True)
solar_date_p50 	, solar_flux_p50 = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout6, unpack=True)
solar_date_p30 	, solar_flux_p30 = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout5, unpack=True)
solar_date_p10	, solar_flux_p10 = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout4, unpack=True)
solar_date_p5 	, solar_flux_p5 = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout3, unpack=True)
solar_date_p1	, solar_flux_p1 = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout2, unpack=True)

x_p100 = array("d",solar_date_p100)
y_p100 = array("d",solar_flux_p100)
x_p50 = array("d",solar_date_p50)
y_p50 = array("d",solar_flux_p50)
x_p30 = array("d",solar_date_p30)
y_p30 = array("d",solar_flux_p30)
x_p10 = array("d",solar_date_p10)
y_p10 = array("d",solar_flux_p10)
x_p5 = array("d",solar_date_p5)
y_p5 = array("d",solar_flux_p5)
x_p1 = array("d",solar_date_p1)
y_p1 = array("d",solar_flux_p1)

solar_graphs_p100 = TGraph(len(x_p100), x_p100, y_p100)
solar_graphs_p100.SetMarkerColor(2005)
solar_graphs_p100.SetLineColor(2005)
solar_graphs_p50 = TGraph(len(x_p50), x_p50, y_p50)
solar_graphs_p50.SetMarkerColor(2004)
solar_graphs_p50.SetLineColor(2004)
solar_graphs_p30 = TGraph(len(x_p30), x_p30, y_p30)
solar_graphs_p30.SetMarkerColor(2003)
solar_graphs_p30.SetLineColor(2003)
solar_graphs_p10 = TGraph(len(x_p10), x_p10, y_p10)
solar_graphs_p10.SetMarkerColor(2002)
solar_graphs_p10.SetLineColor(2002)
solar_graphs_p5 = TGraph(len(x_p5), x_p5, y_p5)
solar_graphs_p5.SetMarkerColor(2001)
solar_graphs_p5.SetLineColor(2001)
solar_graphs_p1 = TGraph(len(x_p1), x_p1, y_p1)
solar_graphs_p1.SetMarkerColor(2000)
solar_graphs_p1.SetLineColor(2000)


mg_p = TMultiGraph()
mg_p.SetTitle("Solar Flare Proton Flux Data GOES15 Satelite;Modified Julian Date MST [days]; Proton Flux [P^{+}/cm#bullets#bulletsr]")

leg = TLegend(0.1, 0.7, 0.3, 0.9)
leg.SetFillColor(0)
leg.SetHeader("test legend")
leg.AddEntry(solar_graphs_p1, "graph 1", "lp")
leg.AddEntry(solar_graphs_p5, "graph 2", "lp")
leg.AddEntry(solar_graphs_p10, "graph 3", "lp")
leg.Draw()

mg_p.Add(solar_graphs_p100)
mg_p.Add(solar_graphs_p50)
mg_p.Add(solar_graphs_p30)
mg_p.Add(solar_graphs_p10)
mg_p.Add(solar_graphs_p5)
mg_p.Add(solar_graphs_p1)
mg_p.SetMinimum(0);  
mg_p.Draw("AL")
mg_p.GetXaxis().SetLimits(first_day,last_day); 
pproton2.Modified();
pproton2.Update()
mg_p.GetYaxis().SetTitleOffset(1.4);
mg_p.GetXaxis().SetTitleOffset(1.4);

pproton1.cd()
real_gr.Draw("AP")
#fit_gr.Draw("Elsame")
real_gr.GetYaxis().SetTitleOffset(1.4);
real_gr.GetXaxis().SetTitleOffset(1.4);
real_gr.GetXaxis().SetLimits(first_day,last_day);
pproton1.Update()

##########################################################################################
### 2nd proton data
##########################################################################################	
c2proton = TCanvas('c2pro1ton', '1Proton Flux Maximum Activity',600,1800)
c2proton.Draw()
c2proton.cd()
p2proton1 = TPad('p2prot1on1','p21',0.0,0.5,1,1)
p2proton1.SetGrid()
p2proton1.Draw()
p2proton2 = TPad('p2prot1on2','p21',0.0,0.0,1,0.5)
p2proton2.SetGrid()
p2proton2.Draw()
p2proton2.cd()

juls, countsnstuff ,zilch,errorforcounts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/'+error_for2, unpack=True)
xya = array("d",juls)
yya = array("d",countsnstuff)
errorya = array("d",errorforcounts)
zeroya = array("d",zilch)
real_gr_fit = TGraphErrors(len(xya), xya, yya , zeroya , errorya)
real_gr_fit.SetLineColor(40);
real_gr_fit.SetMarkerColor(40);

juliandate , zeroes , er117, c133, er133, net, ernet, exp, realcount, livetime,fwhm, excesscounts = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/'+real_fileout, unpack=True)
xreal = array("d",juliandate)
yreal = array("d",net)
errorreal = array("d",ernet)
zeros = array("d",zeroes)
real_gr_error = TGraphErrors(len(xreal), xreal, yreal , zeros , errorreal)
real_gr_error.SetMarkerStyle(7)

mg_p.Draw("AL")
mg_p.SetMinimum(0);
mg_p.GetXaxis().SetLimits(56565,56568); 
mg_p.GetYaxis().SetTitleOffset(1.4);
mg_p.GetXaxis().SetTitleOffset(1.4);
p2proton2.Modified();
p2proton2.Update()

p2proton1.cd()
error_for = TMultiGraph()
error_for.SetTitle("Normalized Counts;Modified Julian Date MST [days];Normalized Counts")
error_for.Add(real_gr_fit,"p")
error_for.Add(real_gr_error,"p")
error_for.Draw("a")
error_for.GetYaxis().SetTitleOffset(1.4);
error_for.GetYaxis().SetTitleOffset(1.4);
error_for.GetXaxis().SetTitleOffset(1.4);
error_for.GetHistogram().SetMaximum(1025);        
error_for.GetHistogram().SetMinimum(1015); 
error_for.GetXaxis().SetLimits(56565,56568);
p2proton1.Modified();
p2proton1.Update()

##########################################################################################
### 3rd proton data
##########################################################################################	

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

mg_p.Draw("AL")
mg_p.SetMinimum(0);
mg_p.GetXaxis().SetLimits(56664,56668); 
mg_p.GetYaxis().SetTitleOffset(1.4);
mg_p.GetXaxis().SetTitleOffset(1.4);
p23proton2.Modified();
p23proton2.Update()

p23proton1.cd()
error_for.Draw("a")
error_for.GetYaxis().SetTitleOffset(1.4);
error_for.GetYaxis().SetTitleOffset(1.4);
error_for.GetXaxis().SetTitleOffset(1.4);
error_for.GetHistogram().SetMaximum(990);        
error_for.GetHistogram().SetMinimum(980); 
error_for.GetXaxis().SetLimits(56664,56668);
p23proton1.Modified();
p23proton1.Update()

##########################################################################################
### end
##########################################################################################	


ce = TCanvas('ce', 'Electron Flux',600,1800)
ce.Draw()
ce.cd()
pe1 = TPad('pe1','p',0.0,0.5,1,1)
pe1.SetGrid()
pe1.Draw()
pe2 = TPad('pe2','p',0.0,0.0,1,0.5)
pe2.SetGrid()
pe2.Draw()
pe2.cd()
print("Calculating Solar Electron Data ...")
solar_date_e0pt8 	, solar_flux_e0pt8 	= numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' 	+ fileout, unpack=True)
solar_date_e2 		, solar_flux_e2 	= numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' 	+ fileout1, unpack=True)
x_e0pt8 = array("d",solar_date_e0pt8)
y_e0pt8 = array("d",solar_flux_e0pt8)
x_e2 = array("d",solar_date_e2)
y_e2 = array("d",solar_flux_e2)
solar_graphs_e0pt8 = TGraph(len(x_e0pt8), x_e0pt8, y_e0pt8)
solar_graphs_e0pt8.SetMarkerColor(2005)
solar_graphs_e0pt8.SetLineColor(2005)
solar_graphs_e2 = TGraph(len(x_e2), x_e2, y_e2)
solar_graphs_e2.SetMarkerColor(2004)
solar_graphs_e2.SetLineColor(2004)
mg_e = TMultiGraph()
mg_e.SetTitle("Solar Flare Electron Flux Data GOES15 Satelite; Modified Julian Date MST [days]; Electron Flux [P/cm sr]")
mg_e.Add(solar_graphs_e0pt8)
mg_e.Add(solar_graphs_e2)
#mg_e.SetMaximum(10000);           
mg_e.SetMinimum(0);  
mg_e.Draw("AL")
if SYNC == True:
	mg_e.GetXaxis().SetLimits(first_day,last_day); 
pe2.Modified();
pe2.Update()

pe1.cd()
real_gr.Draw("ap")
#fit_gr.Draw("Elsame")
real_gr.GetXaxis().SetLimits(first_day,last_day);
pe1.Update()


##########################################################################################
### Get a nice ave
##########################################################################################	
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

##########################################################################################
### Solar Data XRAY
##########################################################################################	
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

fileout_xray1  	= "XRAY_5m_Long.txt"
fileout_xray2 	= "XRAY_5m_Short.txt"

if SOLAR_DATA == True:
	try:
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout_xray1)
		os.remove('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout_xray2)
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
					
solar_date_xray_l  , solar_flux_xray_l = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout_xray1, unpack=True)
solar_date_xray_s 	, solar_flux_xray_s = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout_xray2, unpack=True)

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

if SYNC == True:
	mg_xray.GetXaxis().SetLimits(56546.9071712,56672.8875958) 
px2.Modified();
px2.Update()
xl = []
xld = []
xs = []
xsd = []
for i in range(len(solar_flux_xray_l)):
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

##########################################################################################
### xray on top
##########################################################################################

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
solar_date_xray_l  , solar_flux_xray_l = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout_xray1, unpack=True)
solar_date_xray_s 	, solar_flux_xray_s = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/' + fileout_xray2, unpack=True)

juliandate , away, errors = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/away_from_mean.txt', unpack=True)
x5 = array("d",juliandate)
y5 = array("d",away)
detrend_gr = TGraph(len(x5), x5, y5)

for i in range(len(solar_flux_xray_l)):
	xl.append(solar_flux_xray_l[i]*10+0.998)
	xld.append(solar_date_xray_l[i])
	xs.append(solar_flux_xray_s[i]*10+0.998)
	xsd.append(solar_date_xray_s[i])
print(max(xl))
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

axis7 = TGaxis(last_day,0.998,last_day,1.0015,0,1.0015*10,510,"+L")
axis7.SetName("axis7")
axis7.SetLabelColor(1)
axis7.SetTitle("X-Ray Flux [Watts/m^{2}]")
axis7.Draw()

legx1 = TLegend(0.45, 0.73, 0.89, 0.89)
legx1.SetFillColor(0)
legx1.AddEntry(detrend_gr, "^{60}Co Count Rate", "lp")
legx1.AddEntry(xldg1, "Long Wavelength X-Rays: 0.1 - 0.8 nm", "lp")
legx1.AddEntry(xsdg1, "Short Wavelength X-Rays: 0.05 - 0.4 nm", "lp")
legx1.Draw()

p71x.Modified()
p71x.Update()

p70x.cd()
pl = []
pld = []
ps = []
psd = []
solar_date_proton_l, solar_flux_proton_l = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/GP_5m_proton_1MeV.txt', unpack=True)
solar_date_proton_s, solar_flux_proton_s = numpy.loadtxt('/Users/spenceraxani/Documents/Nuclear_Decay/Data/GP_5m_proton_10MeV.txt', unpack=True)

for i in range(len(solar_flux_proton_l)):
	pl.append(solar_flux_proton_l[i]/3000000+0.998)
	pld.append(solar_date_proton_l[i])
	ps.append(solar_flux_proton_s[i]/1000000+0.998)
	psd.append(solar_date_proton_s[i])
print(max(pl))
x_protonl2 = array("d",pl)
x_protons2 = array("d",ps)
x_protonld2 = array("d",pld)
x_protonsd2 = array("d",psd)
pldg1 = TGraph(len(x_protonld2), x_protonld2, x_protonl2)
pldg1.SetMarkerColor(2005)
pldg1.SetLineColor(2005)
psdg1 = TGraph(len(x_protonsd2), x_protonsd2, x_protons2)
psdg1.SetMarkerColor(2004)
psdg1.SetLineColor(2004)

mg_norm2 = TMultiGraph()
mg_norm2.SetTitle("s;Modified Julian Date MST [days]; Count Rate [cps]")
mg_norm2.Add(pldg1)
mg_norm2.Add(psdg1)
mg_norm2.Add(detrend_gr)
mg_norm2.Draw("al")
mg_norm2.GetYaxis().SetTitleOffset(1.5)
mg_norm2.GetHistogram().SetMaximum(1.0015)      
mg_norm2.GetHistogram().SetMinimum(0.998)
mg_norm2.GetXaxis().SetLimits(first_day,last_day)

axis2 = TGaxis(last_day,0.998,last_day,1.0015,0,1.0015/3000000,510,"+L")
axis2.SetName("axis2")
axis2.SetLabelColor(1)
axis2.SetTitle("Proton Flux [Watts/m^{2}]")
axis2.Draw()

legp2 = TLegend(0.45, 0.73, 0.89, 0.89)
legp2.SetFillColor(0)
legp2.AddEntry(detrend_gr, "^{60}Co Count Rate", "lp")
legp2.AddEntry(pldg1, "Low Energy Protons: 1 MeV", "lp")
legp2.AddEntry(psdg1, "High Energy Protons: 10 MeV", "lp")
legp2.Draw()

p70x.Modified()
p70x.Update()

'''
BinCanvas = TCanvas('BinCanvas','Peakx',600,1800)
BinCanvas.Draw()
BinCanvas.cd()
BinPad2 = TPad('binpad2','px',0.0,0.5,1,1)
BinPad2.SetGrid()
BinPad2.Draw()
BinPad1 = TPad('binpad2','px',0.0,0.0,1,0.5)
BinPad1.SetGrid()
BinPad1.Draw()
BinPad1.cd()

xray_bin = "/Users/spenceraxani/Documents/499_Thesis/data/datapack/x_ray_binned.txt"
n_points = 500.0

if BUN == True:
	try:
		os.remove(xray_bin)
	except OSError:
		pass
	xray_date  , xray_long = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/XRAY_5m_Long.txt', unpack=True)
	increment = (last_day - first_day)/n_points
	counter = 0
	x_ray_sums = 0
	x_ray_counts = []
	x_ray_date = []
	for i in range(n_points):
		if xray_date[i] >= (increment * i + first_day) and xray_date[i] >= (increment * (i + 1) +first_day):
			counter += 1
			x_ray_sums += xray_date[i]
		x_ray_sums = xray_sums/counter
		x_ray_counts.append(x_ray_sums)
		x_ray_date
		counter = 0 

	count_date  , counts, error = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/away_from_mean.txt', unpack=True)
	xray_hist = TH1F('X-Ray Histrogram','xh',n_points,first_day,last_day)
	xray_hist_counts = TH1F('X-Ray counts','xh',n_points,first_day,last_day)
	xray_hist.SetTitle("; Modified Julian Date MST [days];counts")
	for i in range(len(xray_date)):
		#xray_hist.Fill(xray_date[i],xray_long[i])	
		xray_hist_counts.Fill(xray_date[i])


	for j in range(n_points):
	 	xraybin = open(xray_bin,'a')
	 	xraybin.write(str(xray_hist.GetBinContent[j]) + "\t" + str(deadtime) + "\n" )


xray_hist.Draw()
BinPad2.cd()
xray_hist_counts.Draw()
'''

##########################################################################################
### Calculate Chi Squared
##########################################################################################
'''
if CHI == True:
	chi_fileout = '/Users/spenceraxani/Documents/499_Thesis/data/datapack/chi_squared_data.txt'
	chi = open(chi_fileout,'r')
	chi_squared_value = 0
	toy_fileout = '/Users/spenceraxani/Documents/499_Thesis/data/datapack/toy_chi_squared_data.txt'
	toy = open(toy_fileout,'r')
	toy_chi_squared_value = 0
	for line in chi:
		columns = line.split('\t')
		columns = [col.strip() for col in columns]
		num_exp = float(columns[2])
		error_obs = float(columns[1])
		num_obs = float(columns[0])
		chi_squared_value += (num_obs - num_exp)**2/num_exp
	print("Exponential Decay Chi Square = " + str(chi_squared_value))
	
	for line in toy:
		columns = line.split('\t')
		columns = [col.strip() for col in columns]
		num_exp = float(columns[2])
		error_obs = float(columns[1])
		num_obs = float(columns[0])
		toy_chi_squared_value += (num_obs - num_exp)**2/num_exp
	print("Toy Chi Square = " + str(toy_chi_squared_value))
	'''
##########################################################################################
### Sample graphs
##########################################################################################		

c7 = TCanvas('c7','Peak',600,1800)
c7.Draw()
c7.cd()
p70 = TPad('p70','p',0.0,0.5,1,1)
p70.SetGrid()
p70.Draw()
p71 = TPad('p71','p',0.0,0.0,1,0.5)
p71.SetGrid()
p71.Draw()

p71.cd()
p71.SetLogy()
h1 = TH1F('Cobalt 60','Cobalt 60',16384,0,16384)
h1.SetTitle("L2-008 Germanium Detector ^{60}CO Spectrum;Energy [keV]; Count Rate [cps]")
h1.Draw("C")
bin = 0
spectrum = open(sample_fileout, 'r')
lines = spectrum.readlines()
for line in lines:
	bin += 1
	line = int(line)
	h1.SetBinContent(bin, line)
	
s = TSpectrum(2)
nfound = s.Search(h1,1,"new")
peaks = s.GetPositionX()
print("Found candidate peaks = " + str(peaks[0]) + ' ' + str(peaks[1]))

bin_difference = peaks[1]-peaks[0]
energy_difference = 1332.492 - 1173.228

bins_per_energy = bin_difference/energy_difference

low117 = 0
high117 = 0 
binlow117 = 0
binhigh117 = 0
for i in range(100):
	low117 = h1.GetBinContent(int(peaks[0])-i)
	binlow117 = int(peaks[0])-i
	if low117 < h1.GetBinContent(int(peaks[0]))/2:
		break
for i in range(100):
	high117 = h1.GetBinContent(int(peaks[0])+i)
	binhigh117 = int(peaks[0])+i
	if high117 < h1.GetBinContent(int(peaks[0]))/2:
		break
FWHM_energy = (binhigh117 - binlow117)/bins_per_energy
FWHM_bin = binhigh117 - binlow117
print("The FWHM is : "+str(FWHM_energy)+" KeV")

#h1.Fit("gaus",'','', peaks[0]-20,peaks[0]+20 )
c7.Update()
c7.cd(2)
	
p71.Modified();
p71.Update()

p70.cd()
p70.SetLogy()
h2 = TH1F('Cobalt 60 spectrum','Cobalt 60 spectrum',16384,0,int(16384/bins_per_energy))
h2.SetTitle("L2-008 Germanium Detector ^{60}CO Spectrum; Energy [keV]; Count Rate [cps]")
h2.Draw("C")
legend=TLegend(0.1,0.6,0.8,0.8)
legend.AddEntry(h2, h2.GetName(), "l")
legend.Draw()
p70.Update()

bin = 0
spectrum = open(sample_fileout, 'r')
lines = spectrum.readlines()
for line in lines:
	bin += 1
	line = int(line)
	h2.SetBinContent(bin, line)

s = TSpectrum(2)
nfound = s.Search(h2,1,"new")
peaks = s.GetPositionX()
height = s.GetPositionY()
print("Found candidate peaks = " + str(peaks[0]) + 'keV with ' +str(height[0]) +' counts, and '+ str(peaks[1]) + 'keV with ' +str(height[1])+ ' counts')

#g1 = TF1('m1', "gaus", int(binlow117/bins_per_energy),int(binhigh117/bins_per_energy))
#g2 = TF1('m2', "expo", peaks[0]-30,peaks[0]-6 )
#g3 = TF1('m2', "expo", peaks[0]+4,peaks[0]+30 )
#h2.Fit(g1,"R") 
#h2.Fit(g2,"R+") 
#h2.Fit(g3,"R+") 

p70.Modified();
p70.Update()
##########################################################################################
### X^2
##########################################################################################
'''
if CHIgraph == True:
	CHI_data = "CHI_data.txt"
	try:
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + CHI_data)
	except OSError:
		pass
		
	juliandate , c117 , er117, c133, er133, net, ernet, exp, counts, livetime = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+ real_fileout, unpack=True)
	x = array("d",juliandate)
	y = array("d",net)
	chi_square = 0
	for xo in range(0,50):
		No = 1040.0 + xo/10.0
		for yo in range(0,40):
			tau = 1/(0.00034 + yo/10000)
			chi_square = 0
			for i in range(len(x)):
		 		chi_square += (x[i] - No*math.exp(-y[i]/tau))/math.sqrt(x[i])
			chi2 = math.pow(chi_square,2)
			CHI_fileout = open(CHI_data,'a')
			CHI_fileout.write(str(No) + "\t" + str(tau) + "\t" + str(chi2) + "\n" )
			CHI_fileout.close()
'''

##########################################################################################
### Sun rotation
##########################################################################################	
'''
c1test = TCanvas('c1test', 'Ptest',600,1800)
c1test.Draw()
c1test.cd()
p11test = TPad('p11teast','p11taest',0,0,1,1)
p11test.SetGrid()
p11test.Draw()
p11test.cd()
p11test.GetFrame().SetFillColor(21);
p11test.GetFrame().SetBorderSize(12);
mg_away.GetXaxis().SetTitleOffset(1.14);
mg_away.GetYaxis().SetTitleOffset(1.6);
mg_away.Draw("a")

overlay = TPad("a","a",0.,0.,1,1);
overlay.SetFillStyle(4000);
overlay.SetFillColor(0);
overlay.SetFrameFillStyle(4000);
overlay.Draw();
overlay.cd();

fa1 = TF1("fa1","0.5*sin(x*0.35)-0.3",0,90)
fa1.SetTitle('')
hframe = overlay.DrawFrame(0,0,30,1)
hframe.GetXaxis().SetTitleOffset(99)
hframe.GetYaxis().SetTitleOffset(99)
#mg_xray.Draw("AL")
fa1.Draw("al")
'''
'''
c9 = TCanvas('c9', 'Peak',800,800)
c9.Draw()
c9.cd()
p90 = TPad('p80','p',0.05,0.59,0.95,.98)
p90.Draw()
p91 = TPad('p81','p',0.05,0.05,0.95,0.55)
p91.Draw()

p91.cd()

count =0 
for file in glob.glob('/Users/spenceraxani/Documents/499_Thesis/data/datapack/data/2013_Aug_17_L2-008_germanium/*.Spe'):
	count += 1
integral_spectrum = TH1F('Compton Valley Integral','Compton Valley Integral', count,0,count)
integral_spectrum.SetTitle("Compton Valley Integral; Bin; Integrated Count Rate [cps]")
integral_spectrum.Draw()

low_compton = 4600
high_compton = 4800
bins_compton = high_compton - low_compton

if FLAT == True:
	try:
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + flat_fileout)
	except OSError:
		pass
	flat_output = open(flat_fileout,'a')
	file_number = 0
	flat = [] 
	integral_spectrum.Draw()

	for file in glob.glob('/Users/spenceraxani/Documents/499_Thesis/data/datapack/data/2013_Aug_17_L2-008_germanium/*.Spe'):
		file_number += 1
		sum = 0
		list = file
		fh = open(file)
		lines = fh.readlines()
		counter = 0
		for line in lines:
			counter += 1
			columns = line.split('\t')
			columns = [col.strip() for col in columns]
			if counter > low_compton and counter < high_compton:
				sum =+ int(columns[0])
			else:
			 	continue
		flat_output.write(str(file_number) + "\t" + str(sum) +"\n")
		integral_spectrum.SetBinContent(file_number,sum)
p91.Modified();
p91.Update()



fh = open('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' +flat_fileout)
lines = fh.readlines()
for line in lines:
	columns = line.split('\t')
	columns = [col.strip() for col in columns]
	integral_spectrum.Fill(int(columns[0]),int(columns[1]))

#if SYNC == True:
#	integral_spectrum.GetXaxis().SetLimits(first_day,last_day); 
p91.Modified()
p91.Update()	

p90.cd()

gaus_int= TH1F('Compton Valley Bins','Compton Valley Bins', 200,200,300)

fh = open('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' +flat_fileout)
lines = fh.readlines()
for line in lines:
	columns = line.split('\t')
	columns = [col.strip() for col in columns]
	gaus_int.Fill(int(columns[1]),1)
gaus_int.Draw()
gaus_int.SetLineColor(2)
c9.Update()
g10 = TF1('m10', "gaus", 200,300)
gaus_int.Fit(g10,"R") 

p90.Modified()
p90.Update()

gaus_int_theory= TH1F('Compton Valley Bins Theory','Compton Valley Binsa', 200,200,300)

average_number_of_counts = 0
rpois = TRandom3(1242)

for j in range(count):
	for i in range(bins_compton):
		n_obs = rpois.Poisson(gaus_int.GetMean());
		average_number_of_counts += n_obs
	ave_num = average_number_of_counts/bins_compton
	gaus_int_theory.Fill(int(ave_num),1)
	average_number_of_counts = 0	
gaus_int_theory.Draw('same')
gaus_int_theory.SetLineColor(3)

p90.Modified()
p90.Update()


'''
raw_input("done")
