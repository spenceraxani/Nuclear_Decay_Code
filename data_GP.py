import os
import glob
import math

#ELECTRON
fileout = "GP_5m_electron_2MeV.txt"
fileout1 = "GP_5m_electron_0.8MeV.txt"
fileout2= "GP_5m_proton_1MeV.txt"
fileout3 = "GP_5m_proton_5MeV.txt"
fileout4 = "GP_5m_proton_10MeV.txt"
fileout5 = "GP_5m_proton_30MeV.txt"
fileout6 = "GP_5m_proton_50MeV.txt"
fileout7 = "GP_5m_proton_100MeV.txt"

try:
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+fileout)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+fileout1)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+fileout2)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+fileout3)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+fileout4)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+fileout5)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+fileout6)
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+fileout7)

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
		
			print(str(columns[-2]) + str(columns[-2]) + str(columns[-2]))
			output = open(fileout,'a')
			if float(columns[-2]) > 0: #pull 2 MeV electron Data
				output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-2]) + "\n")
			else:
				continue
			
			output = open(fileout1,'a') #pull 0.8 MeV electron Data
			if columns[-4] > 0:
				output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-4]) + "\n")
			else:
				continue
			
			output = open(fileout7,'a') #pull 100 MeV proton Data
			if columns[-6] > 0:
				output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-6]) + "\n")
			else:
				continue
				
			output = open(fileout6,'a') #pull 50 MeV proton Data
			if columns[-8] > 0:
				output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-8]) + "\n")
			else:
				continue
				
			output = open(fileout5,'a') #pull 30 MeV proton Data
			if columns[-10] > 0:
				output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-10]) + "\n")
			else:
				continue
				
			output = open(fileout4,'a') #pull 10 MeV proton Data
			if columns[-12] > 0:
				output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-12]) + "\n")
			else:
				continue	
				
			output = open(fileout3,'a') #pull 5 MeV proton Data
			if columns[-14] > 0:
				output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-14]) + "\n")
			else:
				continue		
				
			output = open(fileout2,'a') #pull 5 MeV proton Data
			if columns[-16] > 0:
				output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-16]) + "\n")
			else:
				continue	
				
			output = open(fileout1,'a') #pull 1 MeV proton Data
			if columns[-18] > 0:
				output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-18]) + "\n")
			else:
				continue			
