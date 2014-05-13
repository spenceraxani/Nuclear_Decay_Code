import os
import glob
import math

fileout = "L2-008_germanium_dead_time.txt"


try:
	os.remove('/Users/spenceraxani/Documents/499_Thesis/data/datapack/'+fileout)

except OSError:
	pass
	
sum = 0


for file in glob.glob('/Users/spenceraxani/Documents/499_Thesis/data/datapack/L2-008_germanium/*.txt'):
	list = file
	fh = open(file)
	lines = fh.readlines()
	sum = sum +1
	deadtime = 0
	for line in lines:
		columns = line.split(' ')
		columns = [col.strip() for col in columns]
		if len(columns) > 2:
			if columns[0] == "Detector":
				print(columns[-6])
				print((float(columns[-6])-800)/float(columns[-6])*100)
				deadtime = (float(columns[-6])-800)/float(columns[-6])*100
				output = open(fileout,'a')
				output.write(str(sum) + "\t" + str(deadtime) + "\n" )

	