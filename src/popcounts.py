'''
FILE: counts.py
PURPOSE: obtain populations from distinct haplotypes
INPUT: none
OUTPUT: 2 csv files with population counts called 	
	1) mostfreq_chr1_159203314-159283887.vcf
	2) identical_chr1_159203314-159283887.vcf
CREATED BY: Anne-Sophie Fratzscher
'''
import config
import pandas as pd
import numpy as np
import math
from pathlib import Path
import time

'''
FUNCTION: getPopsMostFreq(filename)
PURPOSE: gets population/ superpopulation count based on 
	subsamples found (used for most frequent)
INPUT: filename
OUTPUT: csv called "mostfreq_chr1_159203314-159283887.vcf"
'''
def getPopsMostFreq(filename):
	popFile = config.__FILEPATH__ + config.__FOLDERNAME__ + "/mostfreq_" + config.__FILENAME__ 
	# if already have populations
	fileCheck = Path(popFile)
	if fileCheck.is_file():
		return

	df = pd.read_csv(filename, sep="\t")
	df = replacePop(df, 'subsamples')
	df = df.reset_index(drop=True)
	info = df.head(3) #stores id, ref, alt (not used but important)
	df = df.drop(df.head(3).index)

	df = df.sort_values(by='counts', ascending=False)
	df = pd.concat([info, df])
	df.to_csv(popFile, sep="\t", mode='a', index = False)
	
	
'''
FUNCTION: getPopsIdentical(filename)
PURPOSE: gets population/ superpopulation count based on 
	identical samples found (used for multiple populations)
INPUT: filename
OUTPUT: csv called "identical_chr1_159203314-159283887.vcf"
'''
def getPopsIdentical(filename):
	popFile = config.__FILEPATH__ + config.__FOLDERNAME__ + "/identical_" + config.__FILENAME__ 
	
	# if already have populations
	fileCheck = Path(popFile)
	if fileCheck.is_file():
		return

	df = pd.read_csv(filename, sep="\t")
	df = replacePop(df, 'identical')
	df.to_csv(popFile, sep="\t", mode='a', index = False)
	
'''
FUNCTION: replacePop(df, replaceType)
PURPOSE: algortihm to count population/superpopulation
INPUT: dataframe and replaceType (identical or subsample)
OUTPUT: dataframe
'''
def replacePop(df, replaceType):
	tf = pd.read_csv('integrated_call_samples_v3.20130502.ALL.panel.txt', sep="\t")

	# add columns for these pops
	for x in config.__POPS__:
		df[x] = 0
	df['numberOfPops'] = 0

	# add columns for super pops
	for x in config.__SUPERPOPS__:
		df[x] = 0
	df['numberOfSuperpops'] = 0

	# remove ID, ALT, REF
	info = df.head(3) #stores id, ref, alt (not used but important)
	df = df.drop(df.head(3).index)
	df = df.reset_index(drop=True)

	# get pops (using info from 2504 gp file)
	xidx = 0
	for x in df[replaceType]:
		dup = [] # to account for 2 copies of DNA
		dup_pop = []
		dup_super_pop = []
		if not (pd.isnull(x)):
			x = x.split(", ")
			for i in range(0, len(x)):
				# look for sample in "samples" df and get pop, super_pop, gender
				if x[i][2:] not in dup:
					idx = tf[tf['sample'] == x[i][2:]].index.tolist()
					sample_pop = tf.iloc[idx[0]]['pop']
					sample_super_pop = tf.iloc[idx[0]]['super_pop']

					# increment count if new pop
					if sample_pop not in dup_pop:
						df.at[xidx, 'numberOfPops']+=1
					# increment count if new super pop
					if sample_super_pop not in dup_super_pop:
						df.at[xidx, 'numberOfSuperpops']+=1

					df.at[xidx, sample_pop] += 1 #increment pop
					df.at[xidx, sample_super_pop] +=1 # increment superpop

					#append gender, pop, superpop
					dup_pop.append(sample_pop)
					dup_super_pop.append(sample_super_pop)

				dup.append(x[i][2:])
		xidx+=1
	df = pd.concat([info, df])
	return df

def main():
	print("*****STARTING POP COUNTER*****")
	distinctFile = config.__FILEPATH__ + config.__FOLDERNAME__ + "/distinct_" + config.__FILENAME__

	getPopsMostFreq(distinctFile)
	getPopsIdentical(distinctFile)
