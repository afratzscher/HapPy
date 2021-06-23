'''
FILE: distinct.py
PURPOSE: obtain distinct haplotypes
INPUT: none
OUTPUT: 2 csv files -> 1 with counts, 1 with distinct haplotypes
CREATED BY: Anne-Sophie Fratzscher
'''
import config
from pathlib import Path
import pandas as pd
import numpy as np
import math
import itertools

'''
FUNCTION: getCounts(filename)
PURPOSE: gets count for each haplotype (number of times either as duplicate or 
	as subcomponent)
INPUT: filename
OUTPUT: csv with counts, sorted by number of counts (descending)
'''
def getCounts(filename):
	# NOTE: longest haplotype first
	df = pd.read_csv(filename, sep="\t")
	info = df.head(3) #stores id, ref, alt (not used but important)
	df = df.drop(df.head(3).index)
	df = df.loc[(df.iloc[:, 1:-5] != 'N').any(1)] # removes haplotypes where only N
	start = df['start']
	end = df['end']
	indel_num = df['indel_num']
	del df['start']
	del df['end']
	del df['indel_num']
	lengthdf = df['length'] # append later
	posdf = df['POS']
	del df['POS']
	del df['length']
	del df['numHeteroSNPs']

	d = df.loc[:, ~df.columns.str.contains('^Unnamed')] # remove conserved parts

	# get list of samples and list of haplotypes
	x = d.to_string(header=False, index=False,
								index_names=False).split('\n')
	haplotypes = [''.join(ele.split()) for ele in x]

	y = posdf.to_string(header=False, index=False).split('\n')
	samples = [''.join(ele.split()) for ele in y]
	
	identical = samples.copy()
	subSample = samples.copy()
	counts = [0] * len(identical) # DOESNT include sample
	longest = [0] * len(identical) 
		# for a (longer haplotype being compared), 0 intially, 1 if already found longest subsequence,
												# 2 if no longest subsequence found 0 after comparing to all other haplotypes
	found = [0] * len(identical)
		# for b (shorter haplotype being compared), 0 initially, 1 if subsamples already found for haplo

	for ((aidx,a),(bidx,b)) in itertools.combinations(enumerate(haplotypes), 2):
		if haplotypes[aidx] != -1:
			if not (haplotypes[aidx] is None or haplotypes[bidx] is None):
				if (a == b) and (a is not None) and (b is not None): # identical haplotypes
					haplotypes[bidx] = None
					identical[aidx] = identical[aidx] + ", " + identical[bidx]
					identical[bidx] = None
					if subSample[aidx] == None:
						subSample[bidx] = samples[aidx]
					elif subSample[bidx] == None:
						subSample[bidx] = subSample[aidx]
					else:
						subSample[aidx] = subSample[aidx] + ", " + subSample[bidx]
				else: # not identical haplotypes				
					if ((b == (len(samples) - 1)) and (longest[aidx] == 0)): # if didnt find subseqeuence, set longest to 2
						longest[aidx] = 2
					if (longest[aidx] != 2 and found[aidx] != 1):
					# compare SNP by SNP
						occ = [b.find('A'), b.find('G'), b.find('C'),
										b.find('T'), b.find('-')] #, len(b)+1] # find first non N
						# assumes have non N -> earlier case is if ALL N
						startidx = min(i for i in occ if i >= 0) # if one of them is 0, means no N
						endidx = max(b.rfind('A'), b.rfind('G'), b.rfind('C'),
										b.rfind('T'), b.rfind('-')) #len(b)-1)
						if b[startidx:endidx+1] == a[startidx:endidx+1]:
							if (longest[aidx] == 0):
								if subSample[aidx] == None:
									subSample[bidx] = samples[aidx]
								elif subSample[bidx] == None:
									subSample[bidx] = subSample[aidx]
								else:
									subSampleList = set(list(subSample[bidx].split(", "))).union(set(list(subSample[aidx].split(", "))))
									for i in identical[aidx].split(", "):
										subSampleList.remove(i)
									subSample[bidx] = ', '.join(subSampleList)
									subSample[aidx] = None
								longest[aidx] = 1
								found[bidx] = 1

	# get counts
	for i in range(0, len(counts)):
		if subSample[i] is None:
			subSample[i] = identical[i]
		if subSample[i] == samples[i] and (identical[i] is not None):
			subSample[i] = identical[i]
		counts[i] = len(list(subSample[i].split(", ")))

	df.insert(loc=0, column='sampleID', value = samples)
	df.insert(loc=1, column='subsamples', value=subSample)
	df.insert(loc=2, column='identical', value = identical)
	df['start'] = start
	df['end'] = end
	df['length'] = lengthdf + indel_num
	df['counts'] = counts
	
	# format info
	info.rename(columns={'POS':'sampleID'}, inplace=True)
	del info['indel_num']
	info.insert(loc=1, column='haplotypeID', value = '-')
	info.insert(loc=2, column='subsamples', value = '-')
	info.insert(loc=3, column='identical', value = '-')
	info['counts'] = 0
	del info['numHeteroSNPs']

	# remove identical haplotypes
	df = df.dropna(axis=0, subset=['identical'])
	df.insert(loc=1, column='haplotypeID', value=["HAP" + str(i+1) for i in range(len(df.index))])
	# sort by count number (descending)
	countsorted = df.sort_values(by='counts', ascending=False)
	df2 = pd.concat([info, countsorted])
	print(df2)
	df2.to_csv(countFile, sep="\t", mode='a', index = False)
	
	getDistinct(df, info)

'''
FUNCTION: getDistinct(df, samples)
PURPOSE: creates csv with only distinct haplotypes, with samples
	 with that haplotype stored
INPUT: df, info (from count method)
OUTPUT: csv of distinct haplotypes, sorted by length (descending)
'''
def getDistinct(df, info):
	#sort by length
	df = df.sort_values(by='length', ascending=False)
	df = pd.concat([info, df])
	
	df.to_csv(distinctFile, sep="\t", mode='a', index = False)
	
def main():
	print('here')
	print('*****STARTING DISTINCT*****')
	global countFile
	global distinctFile
	filename = config.__FILEPATH__ + "haplotypes_" + config.__FILENAME__
	countFile = config.__FILEPATH__ + "count_" + config.__FILENAME__
	distinctFile = config.__FILEPATH__ +  "distinct_" + config.__FILENAME__

	# if already have counts
	fileCheck = Path(countFile)
	if fileCheck.is_file():
		return

	getCounts(filename)
