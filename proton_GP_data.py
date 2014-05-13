import os
import glob
import math

fileout = "GP_5m_proton_1MeV.txt"
fileout2 = "GP_5m_proton_5MeV.txt"
fileout3 = "GP_5m_proton_10MeV.txt"
fileout4 = "GP_5m_proton_30MeV.txt"
fileout5 = "GP_5m_proton_50MeV.txt"
fileout6 = "GP_5m_proton_100MeV.txt"

try:
	os.remove('/Users/spenceraxani/Desktop/499/data/datapack/'+fileout)
	os.remove('/Users/spenceraxani/Desktop/499/data/datapack/'+fileout2)
	os.remove('/Users/spenceraxani/Desktop/499/data/datapack/'+fileout3)
	os.remove('/Users/spenceraxani/Desktop/499/data/datapack/'+fileout4)
	os.remove('/Users/spenceraxani/Desktop/499/data/datapack/'+fileout5)
	os.remove('/Users/spenceraxani/Desktop/499/data/datapack/'+fileout6)

except OSError:
	pass

for file in glob.glob('/Users/spenceraxani/Desktop/499/data/datapack/particle/*_Gp_part_5m.txt'):
	list = file
	fh = open(file)
	lines = fh.readlines()
	for line in lines:
		columns = line.split(' ')
		columns = [col.strip() for col in columns]
		if columns[0] == "2013":
			print(columns[-6])
			output = open(fileout,'a')
			output.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-16]) + "\n")
			output2 = open(fileout2,'a')
			output2.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-14]) + "\n")
			output3 = open(fileout3,'a')
			output3.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-12]) + "\n")
			output4 = open(fileout4,'a')
			output4.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-10]) + "\n")
			output5 = open(fileout5,'a')
			output5.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-8]) + "\n")
			output6 = open(fileout6,'a')
			output6.write(str(float(columns[7]) + 4.246284501061571125265392*float(columns[4])/10000) + "  " +  str(columns[-6]) + "\n")
