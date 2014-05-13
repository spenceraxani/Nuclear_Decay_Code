import os
import glob
import math
import time
import sys
import numpy
from ROOT import gROOT, TCanvas, TF1, TGraph, TGraphErrors, gApplication, TMultiGraph, TColor, TLegend, TPad, TAxis, TH1F, gStyle, TH1D,TMath
from array import array
from random import random

SOLAR_DATA = False
TEMP_DATA = False
SYNC = False
CHI = False
LIKELIHOOD = False
FLAT = False

c0 = TCanvas('c0', 'proton',800,800)
c0.Draw()
c0.cd()
p0 = TPad('p0','p',0.05,0.59,0.95,.98)
p0.Draw()
p1 = TPad('p1','p',0.05,0.05,0.95,0.55)
p1.Draw()

gStyle.SetEndErrorSize(0)

##########################################################################################
### Goodness of Fit
##########################################################################################		
def exp_counts(t):
	half_life = 1925.28 	#days	
	initial_counts = 780000/live_time	
	counts_exp = (initial_counts)*numpy.exp(-numpy.log(2)*(t)/half_life)
	return(counts_exp)
	
def toy_counts(t): 
	half_life = 1925.28 	#days	
	initial_counts = 780000/live_time
	counts_toy = (initial_counts)*numpy.exp(-numpy.log(2)*(t)/half_life) + 10*numpy.sin(3*t/numpy.pi + numpy.pi/4)
	return(counts_toy)
	
def iround(x):
    return int(round(x) - .5) + (x > 0)

##########################################################################################
###	Data statistics
##########################################################################################
fileout = "L2-008_germanium_data.txt"
deadtimefileout = "L2-008_germanium_dead_time.txt"
chi_fileout = "chi_squared_data.txt"
expected_fileout = "expected_data.txt"
toy_fileout = "toy_data.txt"
chi_fileout_toy = "toy_chi_squared_data.txt"
likelihood_fileout = "likelihood_data.txt"
peak_fileout = "peak_data.txt"
FWHM_fileout = "FWHM_data.txt"
sample_fileout = "CO60_800LIVE_000000.Spe"
flat_fileout = "flat_data.txt"
try:
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + deadtimefileout)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + chi_fileout)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + expected_fileout)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + toy_fileout)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + chi_fileout_toy)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + likelihood_fileout)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + peak_fileout)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + FWHM_fileout)
	
except OSError:
	pass
	
number = 0

num_files= 1

counter=0
error117 = 0
error133 = 0

list_117 = []
list_133 = []
list_time =[]
list_deat_time = []
list_error117 = []
list_error133 = []
first_day = 0
last_day = 0
live_time = 0

for file in glob.glob('/Users/spenceraxani/Documents/499_Thesis/data/datapack/L2-008_germanium/*.txt'):
	list = file
	print(list)
	fh = open(file)
	lines = fh.readlines()
	number = number +1
	net117 = 0
	net133 = 0
	deadtime = 0
	counter = counter + 1
	#print(file)
	
	for line in lines:
		line = ' '.join(line.split())
		#print(line)
		columns = line.split(' ')
		#print(columns)
		columns = [col.strip() for col in columns]
		if len(columns) > 2:
			if columns[0] == "Detector":
				live_time = float(columns[-1])
				deadtime = 100*(float(columns[-4])-float(columns[-1]))/float(columns[-4])
				date = columns[3]
				clock = columns[5]
				(hour, min, sec) = clock.split(':')
				(day, month, year) = date.split('-')
				actual_time  = float((int(hour) * 3600 + int(min) * 60 + int(sec)))/86400

				if month == "Jan":
					month = 1
				elif month == "Feb":
					month = 2
				elif month == "Mar":
					month = 3
				elif month == "Apr":
					month = 4
				elif month == "May":
					month = 5
				elif month == "Jun":
					month = 6
				elif month == "Jul":
					month = 7
				elif month == "Aug":
					month = 8
				elif month == "Sep":
					month = 9
				elif month == "Oct":
					month = 10
				elif month == "Nov":
					month = 11
				elif month == "Dec":
					month = 12
					
				year = int(year)
				month = int(month)
				day = int(day)
				t = time.mktime((year, month, day, 0, 0, 0, 0, 0, 0))
				time.gmtime(t) 
				proper_time= float(time.gmtime(t)[7]) + float(actual_time)
				#print(proper_time)
				list_time.append(proper_time)
				if number == 1:
					first_day = proper_time + 56293.50
				
			elif columns[0] == "1":
				gross117 = float(columns[3])/live_time
				net117 = float(columns[4])/live_time
				centroid117 = float(columns[6])
				FWHM117 = float(columns[8])
				error117 = float(columns[6])/live_time
				list_117.append(int(net117))
				#print(columns)
			elif columns[0] == "2":
				gross133 =float(columns[3])/live_time
				net133 = float(columns[4])/live_time
				#print(columns)
				centroid133 = float(columns[6])
				FWHM133 = float(columns[8])
				error133 = float(columns[6])/live_time
				list_133.append(int(net133))	
			error = math.sqrt(error117**2 + error133**2)
	
			
	if counter == num_files:
		average_117 = sum(list_117) / float(len(list_117))
		average_133 = sum(list_133) / float(len(list_133))	
		average_time = sum(list_time) / float(len(list_time))	
		
		if number == 1:
			start_time = 56293.50 + average_time
		last_day = list_time[-1] + 56293.50
		
		expected_counts_poisson = exp_counts(56293.50 + average_time - start_time)
		toy_counts_poisson 		= toy_counts(56293.50 + average_time - start_time)
		
		del list_117[:]
		del list_time[:]
		del list_133[:]
		counter = 0
		average_net = (average_117 + average_133)
		output = open(fileout,'a')
	 	output.write(str(56293.50 + average_time) + "\t" + str(average_net) + "\t" + str(0) + "\t" + str(error) +"\n")
	 	deadoutput = open(deadtimefileout,'a')
	 	deadoutput.write(str(56293.50 + average_time) + "\t" + str(deadtime) + "\n" )
	 	chioutput = open(chi_fileout,'a')
	 	chioutput.write(str(average_net) + "\t" + str(error) + "\t" + str(expected_counts_poisson) + "\n" )
	 	toychioutput = open(chi_fileout_toy,'a')
	 	toychioutput.write(str(average_net) + "\t" + str(error) + "\t" + str(toy_counts_poisson) + "\n" )
	 	expected_output = open(expected_fileout,'a')
	 	expected_output.write(str(56293.50 + average_time) + "\t" + str(expected_counts_poisson) + "\n")
	 	toy_output = open(toy_fileout,'a')
	 	toy_output.write(str(56293.50 + average_time) + "\t" + str(toy_counts_poisson) +"\n")
	 	likelihood_output = open(likelihood_fileout,'a')
	 	likelihood_output.write(str(iround(expected_counts_poisson)) + "\t" + str(iround(toy_counts_poisson)) + '\t' + str(iround(average_net)) +"\n")
		peak_output = open(peak_fileout,'a')
	 	peak_output.write(str(56293.50 + average_time) + "\t" + str(centroid117) + '\t' + str(centroid133) +"\n")
	 	FWHM_output = open(FWHM_fileout,'a')
	 	FWHM_output.write(str(56293.50 + average_time) + "\t" + str(FWHM117) + '\t' + str(FWHM133) +"\n")
##########################################################################################
### Make Graphs
##########################################################################################		
x1 , y1 , erx1, ery1 = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+fileout, unpack=True)
x = array("d",x1)
er_x = array("d",erx1)
y = array("d",y1)
er_y = array("d",ery1)
gr = TGraphErrors(len(x), x, y, er_x, er_y)
#gr.SetTitle("Germanium Detector L2-008 ^{60}Co Count Rate;Julian Date [days];Sum Counts [cps]")
gr.SetLineColor(9)

x_exp , y_exp = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + expected_fileout, unpack=True)
x_e = array("d",x_exp)
y_e = array("d",y_exp)
exp = TGraph(len(x_e),x_e,y_e)
#exp.SetTitle("Expected Decay of ^{60}Co ;Julian Date [days];Sum Counts [cps]")
exp.SetMarkerColor(2)
exp.SetMarkerSize(2)
exp.SetLineColor(2)
exp.SetLineWidth(4)

x_toy , y_toy = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + toy_fileout, unpack=True)
x_toy = array("d",x_toy)
y_toy = array("d",y_toy)
toy = TGraph(len(x_toy), x_toy, y_toy)
#toy.SetTitle("Expected Decay of ^{60}Co ;Julian Date [days];Sum Counts [cps]")
toy.SetMarkerColor(2)
toy.SetMarkerSize(2)
toy.SetLineColor(2)
toy.SetLineWidth(4)

multi = TMultiGraph()
multi.SetTitle("Germanium Detector L2-008 ^{60}Co Count Rate;Julian Date [days];Sum Counts [cps]")
multi.Add(gr)
multi.Add(exp)
multi.Add(toy)

#myfit = TF1("myfit", "sin(x)",0, 100)
#myfit.SetLineColor(2)
#myfit.SetLineWidth(10)
#gr.Fit("myfit")
print("the first day of the data set is " + str(first_day))
print("the last day of the data set is " + str(last_day))

c0.cd()
p0.cd()
multi.Draw("AP")
multi.GetXaxis().SetLimits(first_day,last_day);
p0.Update()

##########################################################################################
###DEADTIME FILE
##########################################################################################	
c01 = TCanvas('c01', 'DeadTime',800,800)
c01.Draw()
c01.cd()
p010 = TPad('p110','p',0.05,0.59,0.95,.98)
p010.Draw()
p011 = TPad('p111','p',0.05,0.05,0.95,0.55)
p011.Draw()
p011.cd()
Julian_time		, deadtime	= numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' 	+ deadtimefileout, unpack=True)
x = array("d",Julian_time)
y = array("d",deadtime)
deadtime_graph = TGraph(len(x), x, y)
deadtime_graph.SetTitle("L2-008 Germanium Detector Deadtime; Julian Date [days]; Deadtime [%]")
deadtime_graph.SetMarkerColor(1)
deadtime_graph.SetLineColor(1)
deadtime_graph.Draw("AL")
if SYNC == True:
	deadtime_graph.GetXaxis().SetLimits(first_day,last_day)
p011.Modified()
p011.Update()

p010.cd()
multi.Draw("AC")
multi.GetXaxis().SetLimits(first_day,last_day);
p010.Update()

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

if SOLAR_DATA == True:
	try:
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout)
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout1)
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout2)
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout3)
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout4)
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout5)
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout6)
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout7)

	except OSError:
		pass
	for file in glob.glob('/Users/spenceraxani/Documents/499_Thesis/data/datapack/particle/*_Gp_part_5m.txt'):
		list = file
		fh = open(file)
		lines = fh.readlines()
		for line in lines:
			columns = line.split(' ')
			columns = [col.strip() for col in columns]
			if columns[0] == "2013":
				#print(str(columns[-2]) + str(columns[-2]) + str(columns[-2]))
				output = open(fileout,'a')
				if float(columns[-2]) > 0: # pull 2 MeV electron Data
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "\t" +  str(columns[-2]) + "\n")
				else:
					continue
				
				output = open(fileout1,'a') #pull 0.8 MeV electron Data
				if float(columns[-4]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "\t" +  str(columns[-4]) + "\n")
				else:
					continue
				
				output = open(fileout7,'a') #pull 100 MeV proton Data
				if float(columns[-6]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "\t" +  str(columns[-6]) + "\n")
				else:
					continue
					
				output = open(fileout6,'a') #pull 50 MeV proton Data
				if float(columns[-8]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "\t" +  str(columns[-8]) + "\n")
				else:
					continue
					
				output = open(fileout5,'a') #pull 30 MeV proton Data
				if float(columns[-10]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "\t" +  str(columns[-10]) + "\n")
				else:
					continue
					
				output = open(fileout4,'a') #pull 10 MeV proton Data
				if float(columns[-12]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "\t" +  str(columns[-12]) + "\n")
				else:
					continue	
					
				output = open(fileout3,'a') #pull 5 MeV proton Data
				if float(columns[-14]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "\t" +  str(columns[-14]) + "\n")
				else:
					continue		
					
				output = open(fileout2,'a') #pull 1 MeV proton Data
				if float(columns[-16]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "\t" +  str(columns[-16]) + "\n")
				else:
					continue						
					
Red = TColor(2000,1,0,0)
Redish = TColor(2001,1,.4,0)
Redishish = TColor(2002,1,.8,0)	
Yellow = TColor(2003,0.5,1,0)	
Yellowish = TColor(2004,0,1,1)		
Orange = TColor(2005,0,.5,1)

c0.cd()
p1.cd()
	
solar_date_p100 , solar_flux_p100 = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout7, unpack=True)
solar_date_p50 	, solar_flux_p50 = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout6, unpack=True)
solar_date_p30 	, solar_flux_p30 = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout5, unpack=True)
solar_date_p10	, solar_flux_p10 = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout4, unpack=True)
solar_date_p5 	, solar_flux_p5 = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout3, unpack=True)
solar_date_p1	, solar_flux_p1 = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout2, unpack=True)

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
mg_p.SetTitle("Solar Flare Proton Flux Data GOES15 Satelite; Julian Date [days]; Proton Flux [P^{+}/cm#bullets#bulletsr]")

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

mg_p.Draw("AL")
if SYNC == True:
	mg_p.GetXaxis().SetLimits(first_day,last_day); 
p1.Modified();
p1.Update()

c11 = TCanvas('c11', 'electron',800,800)
c11.Draw()
c11.cd()
p110 = TPad('p110','p',0.05,0.59,0.95,.98)
p110.Draw()
p111 = TPad('p111','p',0.05,0.05,0.95,0.55)
p111.Draw()

p111.cd()

solar_date_e0pt8 	, solar_flux_e0pt8 	= numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' 	+ fileout, unpack=True)
solar_date_e2 		, solar_flux_e2 	= numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' 	+ fileout1, unpack=True)
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
mg_e.SetTitle("Solar Flare Electron Flux Data GOES15 Satelite; Julian Date [days]; Electron Flux [P^{+}/cm#bullets#bulletsr]")
mg_e.Add(solar_graphs_e0pt8)
mg_e.Add(solar_graphs_e2)
mg_e.Draw("AP")
if SYNC == True:
	mg_e.GetXaxis().SetLimits(first_day,last_day); 
p111.Modified();
p111.Update()

p110.cd()
multi.Draw("AC")
multi.GetXaxis().SetLimits(first_day,last_day);
p110.Update()

##########################################################################################
### Solar Data XRAY
##########################################################################################	

#XRAY
c1 = TCanvas('c1', 'Xray',800,800)
c1.Draw()
c1.cd()
p10 = TPad('p10','p',0.05,0.59,0.95,.98)
p10.Draw()
p11 = TPad('p11','p',0.05,0.05,0.95,0.55)
p11.Draw()


fileout_xray1  	= "XRAY_5m_Long.txt"
fileout_xray2 	= "XRAY_5m_Short.txt"

if SOLAR_DATA == True:
	try:
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout_xray1)
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout_xray2)
	except OSError:
		pass
		
	for file in glob.glob('/Users/spenceraxani/Documents/499_Thesis/data/datapack/xray/*_Gp_xr_5m.txt'):
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
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "\t" +  str(columns[-5]) + "\n")
				else:
					continue
				
				output = open(fileout_xray2,'a') #Xray Short
				if float(columns[-1]) > 0 and float(columns[-9]) > 0:
					output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "\t" +  str(columns[-9]) + "\n")
				else:
					continue
					
solar_date_xray_l  , solar_flux_xray_l = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout_xray1, unpack=True)
solar_date_xray_s 	, solar_flux_xray_s = numpy.loadtxt('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + fileout_xray2, unpack=True)

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
p11.cd()
mg_xray = TMultiGraph()
mg_xray.SetTitle("Solar Flare Xray Flux Data GOES15 Satelite; Julian Date [days]; Electron Flux [P^{+}/cm#bullets#bulletsr]")
mg_xray.Add(solar_graphs_xray_l)
mg_xray.Add(solar_graphs_xray_s)
mg_xray.Draw("AL")
if SYNC == True:
	mg_xray.GetXaxis().SetLimits(first_day,last_day) 
p11.Modified();
p11.Update()

p10.cd()
multi.Draw("AC")
multi.GetXaxis().SetLimits(first_day,last_day);
p10.Update()



##########################################################################################
### Temperature
##########################################################################################	
temperature_file = "/Users/spenceraxani/Documents/499_Thesis/data/datapack/Temperature.txt"
temperature_out_file = "/Users/spenceraxani/Documents/499_Thesis/data/datapack/detector_temperature.txt"
humidity_out_file = "/Users/spenceraxani/Documents/499_Thesis/data/datapack/detector_humitidy.txt"

if TEMP_DATA == True:


	try:
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+temperature_out_file)
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+humidity_out_file)
	except OSError:
		pass
	temp = open(temperature_file, 'r')
	
	for line in temp:
		columns = line.split('\t')
		columns = [col.strip() for col in columns]
		#print(columns)
		if columns[-1] == "9":
			temperature = columns[2]
			humidity = columns[3]
			clock1 = columns[0]
			(day, month, year) = clock1.split('/')
			year = int(year)
			month = int(month)
			day = int(day)
			clock2 = columns[1]
			(watch, noon) = clock2.split(' ')
			(hours, min) = watch.split(":")
			if noon == "AM":
				if hours == "12":
					hour = 0
				else:
					hour = int(hours)
			else:
				hour = int(hours) + 12
			#print(str(year) + ' ' + str(month) +' ' + str(day) +' '+ str(hour) +' ' + str(min))
			actual_time  = float((int(hour) * 3600 + int(min)*60))/86400
			current_time = time.mktime((year, month, day, 0, 0, 0, 0, 0, 0))
			a =time.gmtime(current_time) 
			proper_time= float(time.gmtime(current_time)[7]) + float(actual_time)
			#print(proper_time)
			#list_time.append(proper_time)
			temp_out = open(temperature_out_file, 'a')
			temp_out.write(str(56293.50 + proper_time) + "\t" + str(temperature) + "\n")
			hum_out = open(humidity_out_file, 'a')
			hum_out.write(str(56293.50 + proper_time) + "\t" + str(humidity) + "\n")
	
c3 = TCanvas('c3', 'temp',800,800)
c3.Draw()
c3.cd()
p30 = TPad('p30','p',0.05,0.59,0.95,.98)
p30.Draw()
p31 = TPad('p31','p',0.05,0.05,0.95,0.55)
p31.Draw()

julian_date  , temp_det = numpy.loadtxt(temperature_out_file, unpack=True)
julian_date  , hum_det 	= numpy.loadtxt(humidity_out_file, unpack=True)
x_date 	= array("d",julian_date)
y_temp 	= array("d",temp_det)
x_date 	= array("d",julian_date)
y_hum 	= array("d",hum_det)
gedet_temp 	= TGraph(len(x_date), x_date, y_temp)
gedet_temp.SetMarkerColor(2)
gedet_temp.SetLineColor(2)
gedet_hum 	= TGraph(len(x_date), x_date, y_hum)
gedet_hum.SetMarkerColor(1)
gedet_hum .SetLineColor(1)


p31.cd()
mg_temp = TMultiGraph()
mg_temp.SetTitle("L2-008 Germanium Detector Aluminum Temperature; Julian Date [days]; Temperature [#circC]")
mg_temp.Add(gedet_temp)
mg_temp.Draw("Ap")
if SYNC == True:
	mg_temp.GetXaxis().SetLimits(first_day,last_day); 
p31.Modified();
p31.Update()

p30.cd()
multi.Draw("AC")
multi.GetXaxis().SetLimits(first_day,last_day);
p30.Update()

c31 = TCanvas('c31', 'humidity',800,800)
c31.Draw()
c31.cd()
p301 = TPad('p301','p',0.05,0.59,0.95,.98)
p301.Draw()
p311 = TPad('p311','p',0.05,0.05,0.95,0.55)
p311.Draw()

p311.cd()
mg_hum = TMultiGraph()
mg_hum.SetTitle("L2-008 Germanium Detector Aluminum Humidity; Julian Date [days]; Relative Humidity [%]")
mg_hum.Add(gedet_hum)
mg_hum.Draw("Ap")
if SYNC == True:
	mg_hum.GetXaxis().SetLimits(first_day,last_day); 
p311.Modified();
p311.Update()

p301.cd()
multi.Draw("AC")
multi.GetXaxis().SetLimits(first_day,last_day);
p301.Update()


##########################################################################################
### Calculate Chi Squared
##########################################################################################
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
	
	
##########################################################################################
### Peak location
##########################################################################################	
c5 = TCanvas('c5', 'Peak',800,800)
c5.Draw()
c5.cd()
p50 = TPad('p50','p',0.05,0.59,0.95,.98)
p50.Draw()
p51 = TPad('p51','p',0.05,0.05,0.95,0.55)
p51.Draw()
	
julian_date  , peak117, peak133 	= numpy.loadtxt(peak_fileout, unpack=True)
julian_date  , FWHM117, FWHM133 	= numpy.loadtxt(FWHM_fileout, unpack=True)
x_date 	= array("d",julian_date)
y_peak117 	= array("d",peak117)
y_peak133	= array("d",peak133)

x_date 	= array("d",julian_date)
y_FWHM117 = array("d",FWHM117)
y_FWHM133 = array("d",FWHM133)

gedet_peak117 	= TGraph(len(x_date), x_date, y_peak117)
gedet_peak133 	= TGraph(len(x_date), x_date, y_peak133)
gedet_peak117.SetMarkerColor(2)
gedet_peak133.SetMarkerColor(2)
gedet_peak117.SetLineColor(2)
gedet_peak133.SetLineColor(2)
gedet_FWHM117 	= TGraph(len(x_date), x_date, y_FWHM117)
gedet_FWHM133 	= TGraph(len(x_date), x_date, y_FWHM133)
gedet_FWHM117.SetMarkerColor(1)
gedet_FWHM133.SetMarkerColor(1)
gedet_FWHM117 .SetLineColor(1)
gedet_FWHM133 .SetLineColor(1)

p51.cd()
mg_peak = TMultiGraph()
mg_peak.SetTitle("L2-008 Germanium Detector Peak117; Julian Date [days]; Peak Energy [keV]")
mg_peak.Add(gedet_peak117)
mg_peak.Add(gedet_peak133)

mg_peak.Draw("Ap")

#mg_peak.Add(gedet_FWHM117)
#mg_peak.Add(gedet_FWHM133)
if SYNC == True:
	mg_peak.GetXaxis().SetLimits(first_day,last_day); 
p51.Modified();
p51.Update()

p50.cd()
multi.Draw("AC")
multi.GetXaxis().SetLimits(first_day,last_day);
p50.Update()

c6 = TCanvas('c6', 'FWHM',800,800)
c6.Draw()
c6.cd()
p60 = TPad('p60','p',0.05,0.59,0.95,.98)
p60.Draw()
p61 = TPad('p61','p',0.05,0.05,0.95,0.55)
p61.Draw()

p61.cd()
mg_FWHM = TMultiGraph()
mg_FWHM.SetTitle("L2-008 Germanium Detector FWHM; Julian Date [days]; FWHM [keV]")
mg_FWHM.Add(gedet_FWHM117)
mg_FWHM.Add(gedet_FWHM133)

mg_FWHM.Draw("AL")

if SYNC == True:
	mg_FWHM.GetXaxis().SetLimits(first_day,last_day); 
p61.Modified();
p61.Update()

p60.cd()
multi.Draw("AC")
multi.GetXaxis().SetLimits(first_day,last_day);
p60.Update()


	
##########################################################################################
### Sample graphs
##########################################################################################		
c7 = TCanvas('c7', 'Peak',800,800)
c7.Draw()
c7.cd()
p70 = TPad('p70','p',0.05,0.59,0.95,.98)
p70.Draw()
p71 = TPad('p71','p',0.05,0.05,0.95,0.55)
p71.Draw()

p71.cd()
p71.SetLogy()
h1 = TH1F('^{60}CO Spectrum','^{60}CO Spectrum',16384,0,16384)
h1.SetTitle("L2-008 Germanium Detector ^{60}CO Spectrum; Julian Date [days]; Count Rate [cps]")
h1.Draw()
bin = 0
spectrum = open(sample_fileout, 'r')
lines = spectrum.readlines()
for line in lines:
	bin += 1
	line = int(line)
	h1.SetBinContent(bin, line)
	
p71.Modified();
p71.Update()


p70.cd()

count =0 
for file in glob.glob('/Users/spenceraxani/Documents/499_Thesis/data/datapack/L2-008_germanium/*.Spe'):
	count += 1
integral_spectrum = TH1F('Compton Valley Integral','Compton Valley Integral', count,0,count)
integral_spectrum.SetTitle("Compton Valley Integral; Julian Date [days]; Integrated Count Rate [cps]")
integral_spectrum.Draw()

if FLAT == True:
	try:
		os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' + flat_fileout)
	except OSError:
		pass
	flat_output = open(flat_fileout,'a')
	file_number = 0
	flat = [] 
	integral_spectrum.Draw()

	for file in glob.glob('/Users/spenceraxani/Documents/499_Thesis/data/datapack/L2-008_germanium/*.Spe'):
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
			if counter > 4600 and counter < 4800:
				sum =+ int(columns[0])
			else:
			 	continue
		flat_output.write(str(file_number) + "\t" + str(sum) +"\n")
		integral_spectrum.SetBinContent(file_number,sum)
		p70.Modified();
		p70.Update()

fh = open('/Users/spenceraxani/Documents/499_Thesis/data/datapack/' +flat_fileout)
lines = fh.readlines()
for line in lines:
	try:
		columns = line.split('\t')
		columns = [col.strip() for col in columns]
		integral_spectrum.SetBinContent(int(columns[0]),int(columns[1]))
	except IndexError:
		pass
		
if SYNC == True:
	integral_spectrum.GetXaxis().SetLimits(first_day,last_day); 
p70.Modified();
p70.Update()	
raw_input("done")
'''	


##########################################################################################
### Likelihood Calculation
##########################################################################################
if LIKELIHOOD == True:
	
	likelihood = open(likelihood_fileout, 'r')    
	likelihood_significance = 0
 	lines = likelihood.readlines()
	
	log_ratio1_hist = TH1D("log_ratio1_hist","log_ratio1_hist",10000,-500,500);
	log_ratio2_hist = TH1D("log_ratio2_hist","log_ratio2_hist",10000,-500,500);
	l_obs_null_exp_null = 0
	l_obs_other_exp_null = 0
	for line in lines:
	 	columns = line.split('\t')
		columns = [col.strip() for col in columns]
		#NOTE: I'm using stirlings approximation, only valid for high count rates!!!!!!
		n_exp 		= round(float(columns[0]))
		n_toy_exp 	= round(float(columns[1]))
		n_obs 		= round(float(columns[2]))
	
		l_obs_other_exp_null += - n_toy_exp + n_obs*TMath.Log(n_toy_exp) + n_obs*TMath.Log((n_obs!=0)*n_obs - n_obs)

		l_obs_null_exp_null += - n_exp + n_obs*TMath.Log(n_exp) + n_obs*TMath.Log((n_obs!=0)*n_obs - n_obs);

		l_obs_other_exp_other += - n_toy_exp + n_obs*TMath.Log(n_toy_exp) + n_obs*TMath.Log((n_obs!=0)*n_obs - n_obs);
		  
		l_obs_null_exp_other += - n_exp + n_obs*TMath::Log(n_exp) - (TMath::Log(TMath::Factorial((n_obs<160)*n_obs)) + (n_obs>=160)*(n_obs*TMath::Log((n_obs!=0)*n_obs + (n_obs==0)) - n_obs));
	
	  
	  log_ratio1 = l_obs_other_exp_other - l_obs_other_exp_null;
	  log_ratio1_hist.Fill(log_ratio1);
	  l_obs_other_exp_other = 0;
	  l_obs_other_exp_null = 0;
	  
	  log_ratio2 = l_obs_null_exp_other - l_obs_null_exp_null;
	  log_ratio2_hist.Fill(log_ratio2);
	  l_obs_null_exp_other = 0;
	  l_obs_null_exp_null = 0;
    

      numBins = log_ratio1_hist->GetXaxis()->GetNbins();
     
      //integral for p value
      integral = 0;
      
      //integral to find 90% point in curve or 50% or whatever is needed
      percent_integral = 0;
      
      for (int i = 1; i < log_ratio_bins; i++)
	{
	  if(percent_integral < percent_entries)
	    percent_integral += log_ratio1_hist->GetBinContent(i);
	  else
	    {
	      percent_pt = log_ratio1_hist->GetBinCenter(i);
	      //cout << percent_pt << endl;
	      break;
	    }
	}
      
      
      //normalize first
      log_ratio1_hist->Scale(1./log_ratio1_hist->Integral());
      log_ratio2_hist->Scale(1./log_ratio2_hist->Integral());

      //integrate over region
      for (int i = 1; i<log_ratio_bins; i++)
	{
	  if (log_ratio2_hist->GetXaxis()->GetBinCenter(i) >= percent_pt)
	    {
	      integral += TMath::Min(log_ratio2_hist->GetBinContent(i),log_ratio1_hist->GetBinContent(i));
	    }
	  
	}


      //significance value
      nsigma = TMath::Sqrt(2)*TMath::ErfInverse(1. - integral);
      //in case of 0sigma , set it so that we have the maximum for
      //the number of trials used
      if (nsigma == 0.)
	nsigma = TMath::Sqrt(2)*TMath::ErfInverse(1. - (1./float(num_trials)));

      //Write to file, or print to log file

      theta_string.replace(2,2,".");  //put ".' back in
      delta_string.replace(1,2,".");

      cout << theta_string << "    "  << delta_string << "    " << nsigma  << endl;
      //prob_file << pts_theta << "    "  << pts_delta << "    " << prob  << endl;
      
      delete log_ratio1_hist;
      delete log_ratio2_hist;
      delete otherH_N_smear;
      
//    }
  //prob_file.close();
  //close input file
  //pts_to_do.close();

	
	
	
##########################################################################################
### CO60 Plot
##########################################################################################	

c4 = TCanvas('c4', 'D',800,800)
c4.Draw()
c4.cd()

h0 = TH1F("Co60","Co60",16384,0,16384*.082756 + 0.188941);

bin = 0

source_file = "/Users/spenceraxani/Documents/499_Thesis/data/datapack/CO60.txt"
f = open(source_file, 'r')

for line in f:
	#print(line)
	bin = bin +1
	h0.Fill(line, bin)
'''