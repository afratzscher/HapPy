'''
FILE: output.py
PURPOSE: creates clean output for user
INPUT: most frequent file
OUTPUT: sequence and meta files
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

def getMeta(df, metaFile):
	idx = df.columns.get_loc('start')
	meta = df.iloc[3:, idx:]
	haploID = df['haplotypeID']
	meta.insert(0, 'haplotypeID', haploID)
	meta.to_csv(metaFile, sep="\t", mode='a', index = False)

def getSequence(df, seqFile):
	idx = df.columns.get_loc(config.__POPS__[0])
	df = df.iloc[:, 0:idx]
	del df['subsamples']
	del df['identical']
	del df['sampleID']

	list = df.columns.tolist()
	list.insert(1, list.pop(list.index('length'))) 
	list.insert(2, list.pop(list.index('start'))) 
	list.insert(3, list.pop(list.index('end'))) 
	list.insert(4, list.pop(list.index('counts'))) 
	df = df.reindex(columns= list)

	# set length, start, end for REF, ALT (clear for ID)
	df.iloc[0, df.columns.get_loc("haplotypeID")] = 'SNV_ID'
	df.iloc[1, df.columns.get_loc("haplotypeID")] = 'REF'
	df.iloc[2, df.columns.get_loc("haplotypeID")] = 'ALT'
	df.iloc[0, df.columns.get_loc("length")] = '-'
	df.iloc[1, df.columns.get_loc("length")] = config.__END__ - config.__START__ + 1
	df.iloc[2, df.columns.get_loc("length")] = config.__END__ - config.__START__ + 1
	df.iloc[0, df.columns.get_loc("start")] = '-'
	df.iloc[1, df.columns.get_loc("start")] = config.__START__
	df.iloc[2, df.columns.get_loc("start")] = config.__START__
	df.iloc[0, df.columns.get_loc("end")] = '-'
	df.iloc[1, df.columns.get_loc("end")] = config.__END__
	df.iloc[2, df.columns.get_loc("end")] = config.__END__ 
	df.iloc[0, df.columns.get_loc("counts")] = '-'

	unnamed = [col for col in df.columns if 'Unnamed' in col]
	newnames = [None] * len(unnamed)
	df.rename(columns=dict(zip(unnamed, newnames)),inplace=True)

	haploID = df['haplotypeID']
	del df['haplotypeID']
	df = df.replace('N', ' ', regex=True)
	df.insert(0, "haploID", haploID)

	for i in df.columns:
		if i not in ['haplotypeID', 'length', 'start', 'end', 'counts', None, "", 'sampleID']:
			ref = df.loc[1][i]
			df[i] = df[i].replace('-', ref, regex=True) # replace - with nucleotide

	df.to_csv(seqFile, sep="\t", mode='a', index = False)

def outputFile(filename):
	seqFile = config.__FILEPATH__ + config.__FOLDERNAME__ + "/sequence_" + config.__FILENAME__
	seqFile = seqFile[:-3] + 'vcf'
	metaFile = config.__FILEPATH__ + config.__FOLDERNAME__ + "/meta_" + config.__FILENAME__
	metaFile = metaFile[:-3] + 'vcf'

	# if already have output
	fileCheck = Path(seqFile)
	if fileCheck.is_file():
		return

	df = pd.read_csv(filename, sep="\t")

	getMeta(df, metaFile)
	getSequence(df, seqFile)

def main():
	print('*****STARTING OUTPUT FORMATTING*****')
	countFile = config.__FILEPATH__ + config.__FOLDERNAME__ + "mostfreq_" + config.__FILENAME__
	outputFile(countFile)
