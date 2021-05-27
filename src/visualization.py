'''
FILE: visualization.py
PURPOSE: visualizes
INPUT: none
OUTPUT: visualization as pie graph
CREATED BY: Anne-Sophie Fratzscher
'''
import config
import pandas as pd
import numpy as np
import math
from pathlib import Path
import time
import matplotlib.pyplot as plt
import os

def getVisualization(df, name, cscopy, fullFlag):
	if df.empty: # if nothing found
		plt.title("No %s".center(80) %(name))
		if fullFlag == 0: # NOT full length haplotype
				graphname = name + ".png"
		else:
			graphname = "full_" + name + ".png"
		plt.savefig(visualizationFolder + graphname)
		plt.clf()

	else:
		graphnum = 0

		#graph coloring scheme
		cmap = plt.get_cmap("magma")

		for (index, x) in df.iterrows():
			graphnum += 1
			pops = []
			values = []
			labeling = []
			cs = []
			for i in range(0, len(x)):
				if i == (len(x) - 4): # if length
					length = str(int(x[i]))
				elif i == (len(x) - 3): # if start
					start = str(int(x[i]))
				elif i == (len(x) - 2): # if end
					end = str(int(x[i]))
				elif i == (len(x) - 1): # if haplotypeID
					haploID = str(x[i])
				elif x[i] != 0:
					values.append(int(x[i]))
					pops.append(df.columns[i])
					labeling.append(str(df.columns[i]) + "," + str(int(x[i])))
					cs.append(cscopy[i])

			csfull = cmap(cs)

			# visualize (save as file)
			data = values

			fig, ax = plt.subplots()
			ax.axis('equal')
			pie = ax.pie(data, radius = 1, labeldistance=1.005, labels=labeling, 
	              rotatelabels =True, startangle=90,counterclock=False, 
	              colors = csfull, textprops={'fontsize': 8})

			graphType = ''
			if (fullFlag != 0): # 0 if not full, not 0 if full
				graphType = 'Full length'
			if name == 'most_frequent':
				graphType = 'Most frequent'
			else:
				graphType += ' multiple ' + name

			plt.title("%s haplotype (ID: %s): \nlength %s (chr %s: nt %s - %s)" %(graphType, haploID, length, config.__CHR__, start, end))
			if fullFlag == 0: # NOT full length haplotype
				graphname = name + str(graphnum) + ".png"
			else:
				graphname = "full_" + name + str(graphnum) + ".png"
			plt.savefig(visualizationFolder + graphname, dpi=300)
			plt.clf()

def visualizeFull():
	df = pd.read_csv(fullPopFile, sep="\t")
	lastSNPpos = df.columns[df.columns.get_loc("start")-2]
	info = df.head(3) #stores id, ref, alt (not used but important)
	df = df.drop(df.head(3).index)
	df = df[df['start'] <= config.__START__]
	df = df[df['end'] >= int(float(lastSNPpos))]
	df = pd.concat([info, df])
	df.to_csv(fullFile, sep="\t", mode='a', index = False)
	return(df)

def visualizeMultiplePopOnly(df):
	df = df.drop(df.head(3).index) # remove ID, REF, ALT
	df = df[df['numberOfPops'] > 1]
	return df

# drop rows (so only have number of rows = num)
def visualizeNumber(df, num):
	if not (num == 'all'):
		df = df.drop(df.index[num:len(df)])
	return df

def visualizePop(df, num, fullFlag):
	# get only columns for pop
	cols = [c for c in df.columns if c in config.__POPS__]
	cols.append('length')
	cols.append('start')
	cols.append('end')
	cols.append('haplotypeID')
	df = df[cols]
	df = visualizeNumber(df, num)
	getVisualization(df, "pop", config.__CSPOP__, fullFlag)
	return
	
def visualizeSuper(df, num, fullFlag):
	cols = [c for c in df.columns if c in config.__SUPERPOPS__]
	cols.append('length')
	cols.append('start')
	cols.append('end')
	cols.append('haplotypeID')
	df = df[cols]
	df = visualizeNumber(df, num)
	getVisualization(df, "super", config.__CSSUPER__, fullFlag)
	return

def visualizeMostFrequent(df):
	df = df.drop(df.head(3).index) # remove ID, REF, ALT
	df = df.sort_values(['counts', 'length'], ascending=[False, False]) 
		# picks longest most frequent haplotype
	df = df.drop(df.index[1:len(df)])
	cols = [c for c in df.columns if c in config.__POPS__]
	cols.append('length')
	cols.append('start')
	cols.append('end')
	cols.append('haplotypeID')
	df = df[cols]
	fullFlag = 0
	getVisualization(df, "most_frequent", config.__CSPOP__, fullFlag)
	return

def getInput(inputI):
	vistype = inputI
	
	if vistype == 'most frequent':
		full = 'n'
	else:
		full = 'y'

	if full == 'y' or full == 'yes':
		fileCheck = Path(fullFile)
		if fileCheck.is_file():
			df = pd.read_csv(fullFile, sep="\t")
		else:
			df = visualizeFull()
	elif full == 'n' or full == 'no':
		df = pd.read_csv(mostFreqFile, sep="\t")

	if vistype == 'most frequent':
		visualizeMostFrequent(df)
		return

	df = visualizeMultiplePopOnly(df) # only if in multiple populations
	num = 'all'
	
	fullFlag = 1
	if vistype == 'pop':
		visualizePop(df, num, fullFlag)
	elif vistype == 'super':
		visualizeSuper(df, num, fullFlag)
	# elif vistype == 'gender':
	# 	visualizeGender(df, num, fullFlag)		

def main():
	print("*****STARTING VISUALIZATION*****")
	global fullFile
	global fullPopFile
	global mostFreqFile
	global visualizationFolder
	# mostFreqFile = config.__FILEPATH__ + config.__FOLDERNAME__ + "/mostfreq_" + config.__FILENAME__
	# fullPopFile = config.__FILEPATH__ + config.__FOLDERNAME__ + "/identical_" + config.__FILENAME__
	# fullFile = config.__FILEPATH__ + config.__FOLDERNAME__ + "/full_length_haplotypes_" + config.__FILENAME__

	mostFreqFile = config.__FILEPATH__ + "/mostfreq_" + config.__FILENAME__
	fullPopFile = config.__FILEPATH__  + "/identical_" + config.__FILENAME__
	fullFile = config.__FILEPATH__ + "/full_length_haplotypes_" + config.__FILENAME__

	visualizationFolder = (config.__FOLDERPATH__ + "/results/" 
		+ config.__FOLDERNAME__ + "/visualization/")
	# visualizationFolder = (config.__FOLDERPATH__ + "/results/" 
	# 	+ config.__GENENAME__ + "/" + config.__FOLDERNAME__ + "/visualization/")
	if not (os.path.isdir(visualizationFolder)):
		try:
		    Path(visualizationFolder).mkdir(parents=True, exist_ok=True)
		except OSError:
		    print ("Creation of the directory %s failed" % visualizationFolder)
	
	# call getInput using grouping and mostfreq
	if config.__MOSTFREQ__ == 1:
		getInput('most frequent')
	for i in config.__GROUPING__:
		getInput(i)
