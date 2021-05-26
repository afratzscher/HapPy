'''
FILE: cleaner.py
PURPOSE: removes SNPs that are out of specified range from file 
		and removes samples that have >1 heterozygous SNP
INPUT: none
OUTPUT: cleaned csv file named "cleaned_1000G_chr1_159203314-159283887.vcf"
CREATED BY: Anne-Sophie Fratzscher
'''
import config
import read_vcf
from pathlib import Path
import numpy as np
import pandas as pd

'''
FUNCTION: getIndividuals()
PURPOSE: gets ID for 2504 individuals from phase 3 of 1000GP
INPUT: none
OUTPUT: list of individuals from phase 3 of 1000GP
'''
def getIndividuals():
	# this file holds ID for the 2504 individuals 
	indFile = 'integrated_call_samples_v3.20130502.ALL.panel.txt'
	data = pd.read_csv(indFile, sep="\t")
	indi = []
	for i in data['sample']:
		indi.append(i)
	return indi

def manualclean(df):
	df = df[df['ID'] != 'rs529460769']
	return df
'''
FUNCTION: clean()
PURPOSE: removes samples that are not part the 2504 
		individuals for phase 3 AND removes samples
		with > 1 heterozygous SNP
INPUT: list of individuals in 1000GP phase 3
OUTPUT: cleaned csv file named "cleaned_1000G_chr1_159203314-159283887.vcf"
'''
def clean(individuals):
	df = read_vcf.read_vcf(filename)
	
	copy = rangeSelection(df)
	
	columns = copy.columns
	infoNames = {'CHROM','POS','ID','REF','ALT',
		'QUAL','FILTER','INFO','FORMAT'}
	
	# get all sample names
	initsamples = []
	for i in copy.columns:
		if not (i in infoNames):
			initsamples.append(i)

	cleanedsamples = []
	dups = []

	# get samples that are in indFile AND with at most 1 hetero SNP within range
	for i in initsamples:
		if i in individuals and i not in dups:
			# gets counts for different values (0|0, 0|1, 1|0, ...)
			count = copy[i].value_counts()
			totalcount = df[i].value_counts()

			# if homozygous (ie. 0|0, 1|1, 2|2, ... 6|6), dont count
			totalSNPS = 0
			rangeSNPS = 0

			# get range SNPS first
			for c in range(0, len(count)):
				if (count.index[c][0] != count.index[c][-1]):
						rangeSNPS += count[c]

			if rangeSNPS <= 1:
				cleanedsamples.append(i)
				# get total SNPS if part of clean sample
				for c in range(0, len(totalcount)):
					if (totalcount.index[c][0] != totalcount.index[c][-1]):
							totalSNPS += totalcount[c]
				df.at[-1,i] = totalSNPS
			dups.append(i)

	# remove samples from file that are not part of cleanedsamples list
	for i in initsamples:
		if not(i in cleanedsamples):
			del df[i]

	df = manualclean(df)
	df.to_csv(cleanedName, sep="\t", mode='a', index=False)
	
'''CLEAN WITHOUT FILTERING INDIVIDUALS'''
def clean():
	df = read_vcf.read_vcf(filename)
	
	copy = rangeSelection(df)
	
	columns = copy.columns
	infoNames = {'CHROM','POS','ID','REF','ALT',
		'QUAL','FILTER','INFO','FORMAT'}
	
	# get all sample names
	initsamples = []
	for i in copy.columns:
		if not (i in infoNames):
			initsamples.append(i)

	cleanedsamples = []
	dups = []

	# get samples that are in indFile AND with at most 1 hetero SNP within range
	for i in initsamples:
		if i not in dups:
			# gets counts for different values (0|0, 0|1, 1|0, ...)
			count = copy[i].value_counts()
			totalcount = df[i].value_counts()

			# if homozygous (ie. 0|0, 1|1, 2|2, ... 6|6), dont count
			totalSNPS = 0
			rangeSNPS = 0

			# get range SNPS first
			for c in range(0, len(count)):
				if (count.index[c][0] != count.index[c][-1]):
						rangeSNPS += count[c]

			if rangeSNPS <= 1:
				cleanedsamples.append(i)
				# get total SNPS if part of clean sample
				for c in range(0, len(totalcount)):
					if (totalcount.index[c][0] != totalcount.index[c][-1]):
							totalSNPS += totalcount[c]
				df.at[-1,i] = totalSNPS
			dups.append(i)

	# remove samples from file that are not part of cleanedsamples list
	for i in initsamples:
		if not(i in cleanedsamples):
			del df[i]

	df = manualclean(df)
	df.to_csv(cleanedName, sep="\t", mode='a', index=False)

'''
FUNCTION: rangeSelection(df)
PURPOSE: removes SNPs that are outside of range
INPUT: df
OUTPUT: returns ranged df
'''
def rangeSelection(df):
	mask = (df['POS'] >= config.__GENESTART__) & (df['POS'] <= config.__GENEEND__)
	df = df.loc[mask]
	return df
	
def main():
	print('*****STARTING CLEANER*****')
	global filename
	global cleanedName
	filename = config.__FILEPATH__ + config.__FILENAME__
	cleanedName = config.__FILEPATH__ + "cleaned_" + config.__FILENAME__
	
	# if already cleaned
	fileCheck = Path(cleanedName)
	if fileCheck.is_file():
		return

	# individuals = getIndividuals()
	# clean(individuals)
	clean()
