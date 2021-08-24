'''
FILE: haplotype.py
PURPOSE: obtain haplotypes
INPUT: none
OUTPUT: csv file named "haplotypes_1000G_chr1_159203314-159283887.vcf"
CREATED BY: Anne-Sophie Fratzscher
'''
import config
import sequence
import pandas as pd
import numpy as np
import math
from pathlib import Path
import time

'''
FUNCTION: getHaplotypes(filename)
PURPOSE: gets haplotypes
INPUT: filename
OUTPUT: haplotype csv
'''
def getHaplotypes(filename):
	df = pd.read_csv(filename, sep="\t")
	
	#if no samples left (i.e. all samples have >1 hetero SNP), warn user and quit
	if len(df.columns) == 9: # means no samples
		print("*****WARNING: NO SAMPLES REMAINING AFTER CLEANING, ENDING PROGRAM*****")
		print("*****This means there are no haplotypes that can be found using this approach*****")
		quit()

	# only keep samples, SNP, and ref/alt nucleotide
	infoNames = {'CHROM', 'QUAL','FILTER','INFO','FORMAT'}
	for i in infoNames:
		del df[i]

	# call sequence file to get FASTA sequence and add
	df = sequence.main(df) 
	df = df.dropna(axis=0, how="all").reset_index(drop=True) # delete rows with all nan

	print('start unambi')
	# get unambiguous haplotypes
	df = getUnambiguous(df)
	
	# split into maternal and paternal haplotype
	df = splitHaplotypes(df)

	# replace with N all outside of start/end
	for i in df.columns:
		if i not in ['POS', 'ID', 'REF', 'ALT']:
			startval = int(df.loc[len(df.index)-3][i])
			endval = int(df.loc[len(df.index)-2][i])
			if startval > config.__START__:
				startidx = df.index[df['POS'] == startval-1].tolist()
				if startidx:
					idx = startidx[0]
					df.loc[0:idx, i] = 'N'
			if endval < config.__END__:
				endidx = df.index[df['POS'] == endval+1].tolist()
				if endidx:
					idx = endidx[0]
					df.loc[idx:len(df.index) - 4, i] = 'N'

	## TO DO: SPEED THIS UP! -> currently 350ish seconds for this (majority of time in this file)
	# replace values with nucleotides (only if SNP -> if ref, keep as '-')
	# also count indel 
	for (index,row) in df.iterrows():
		if not (index > len(df.index) - 5): # ignore last 5 rows (start, end, #SNPS...)
			if not df.loc[index]['ID'] is None:
				ref = ['-']
				reflen = len(row['REF'])
				alt = row['ALT'].split(',')
				# list of nucleotides (pos 0 = ref, pos 1 - 6 = alts)
				nt = ref + alt
				for i in range(4,len(row)):
					if not df.iloc[index, i] == 'N':
						start = df.iat[-3, i]
						end = df.iat[-2, i]
						indel = df.iat[-1, i]
						# if last SNP is ALT and is full length haplotype, change end site to position of LAST SNP
						if (index == len(df.index)-6): # if last row
							if (end == config.__END__): # if is full length haplotype
								if not (df.iat[index, i] == '0|0'): # if last SNP is alt, alter end site
									val = row[i]
									if "p_" in row.index[i]:
										if (val[0] != '0'): 
											df.iat[-2, i] = int(df.loc[index]['POS']) # if alt, set POS = last SNP pos
											df.iat[index+1, i] = "N"
									else:
										if (val[-1] != '0'):
											df.iat[-2, i] = int(df.loc[index]['POS']) # if alt, set POS = last SNP pos
											df.iat[index+1, i] = "N"
											
						if not (df.iat[index, 0] < start) or (df.iat[index, 0] > end):
							val = row[i]
							if (df.iat[index,i] == '.|.' or df.iat[index,i] == '.'): # if cannot make call for locus, set to 'Q'	
								df.iat[index, i] = 'Q'
							else:		
								if "p_" in row.index[i]: # first occurence
									df.iat[index,i] = nt[int(val[0])]
									if (nt[int(val[0])] != '-'):
										df.iat[-1, i] += (len(nt[int(val[0])]) - reflen)
								else:
									df.iat[index,i] = nt[int(val[-1])]
									if (nt[int(val[-1])] != '-'):
										df.iat[-1, i] += (len(nt[int(val[-1])]) - reflen)						
	# sort by length
	df = df.rename(columns = {'POS': 'idx'})
	df = df.set_index('idx')
	df = sortHaplotypes(df)

	df.to_csv(haplotypeFile, sep="\t", mode='a', index = False)
	return df

'''
FUNCTION: getUnambiguous(df)
PURPOSE: gets part of haplotype that is unambiguous
INPUT: dataframe
OUTPUT: dataframe of unambiguous haplotypes with start and end row
'''
def getUnambiguous(df):
	homozy = ["0|0", "1|1", "2|2", "3|3", "4|4", "5|5", "6|6"]
	
	SNPS = df.tail(1)

	# add rows to store start and end for haplotype
	df.loc[len(df)] = "-" # start
	df.loc[len(df)] = "-" # end

	d = df.dropna()
	d = d.mask(d.isin(homozy)).drop(d.tail(3).index)
	d = d.reset_index(drop=True)

	gene = d[(d['POS'] > config.__GENESTART__) & (d['POS'] <= config.__GENEEND__)] # SNPS within gene
	gene = gene.reset_index(drop=True)
	before = d.drop(d[d['POS'] > config.__GENESTART__].index) # SNPs before gene
	after = d.drop(d[d['POS'] < config.__GENEEND__].index) # SNPs after gene
	after = after.reset_index(drop=True)

	try:
		start = before.notna()[::-1].idxmax()
	except: # no SNPs before gene, then have NONE flag
		start = None
	try:
		end = after.notna().idxmax()
	except: # if no SNP after gene, then have NONE flag
		end = None
	if end is not None:
		endsecond = after.notna().cumsum().eq(2).idxmax() #index of second SNP after gene
	else:
		endsecond = None
	within = gene.notna().sum(axis = 0) # number of hetero SNP in gene (0 or 1)

	for i in df.columns:
		if i in ['POS', 'ID' 'REF', 'ALT']:
			if start is not None:
				start[i] = "-"
			if end is not None:
				end[i] = "-"
		else:
			if within[i] == 1: # if 1 hetero SNP in gene (dont add another hetero)
				# get start site
				if start is None: # if no SNP before gene, set to end of range
					df.iloc[-2][i] = config.__START__
				elif(pd.isna(before.iloc[start[i]][i])): # if no SNP before, set to start 
					df.iloc[-2][i] = config.__START__
				else:
					if (start[i]+1 < len(before)):
						df.iloc[-2][i] = (before.iloc[start[i]]['POS']) + 1
					else:
						df.iloc[-2][i] = config.__START__
				# get end site
				if end is None: # if no SNP after, set to end of range
					df.iloc[-1][i] = config.__END__
				elif(pd.isna(after.iloc[end[i]][i])):
					df.iloc[-1][i] = config.__END__
				else:
					df.iloc[-1][i] = (after.iloc[end[i]]['POS']) - 1 
			else:
				# get start site
				if start is None: # if no SNP before gene, set to end of range
					df.iloc[-2][i] = config.__START__
				elif(pd.isna(before.iloc[start[i]][i])): # if no SNP before, set to start 
					df.iloc[-2][i] = config.__START__
				else:
					if (start[i]+1 < len(before)):
						df.iloc[-2][i] = (before.iloc[start[i]]['POS']) + 1
					else:
						df.iloc[-2][i] = config.__START__
				# get end site -> includes 1 hetero SNP
				if end is None: # if no SNP after, set to end of range
					df.iloc[-1][i] = config.__END__
				elif(pd.isna(after.iloc[endsecond[i]][i])):
					df.iloc[-1][i] = config.__END__
				else:
					df.iloc[-1][i] = (after.iloc[endsecond[i]]['POS']) - 1 

	df.iloc[-1]["POS"] = 'end'
	df.iloc[-2]["POS"] = 'start'
	df.iloc[-3]["POS"] = 'numHeteroSNPs'
	df.iloc[-3]["ID"] = 0
	df.iloc[-3]["REF"] = 0
	df.iloc[-3]["ALT"] = 0
	df.iloc[-1]["ID"] = 0
	df.iloc[-1]["REF"] = 0
	df.iloc[-1]["ALT"] = 0
	df.iloc[-2]["REF"] = 0
	df.iloc[-2]["ALT"] = 0
	df.iloc[-2]["ID"] = 0

	# save added/removed nucleotide length
	df.loc[len(df)] = 0 # indel length
	df.iloc[-1]["POS"] = 'indel_num'

	return df	

'''
FUNCTION: splitHaplotypes(df)
PURPOSE: splits into maternal and paternal haplotypes
INPUT: dataframe
OUTPUT: dataframe with maternal and paternal haplotypes
'''
def splitHaplotypes(df):
	for i in df.columns:
		if not (i == 'POS' or i == 'ID' or i == 'REF' or i == 'ALT'):
			df['m_'+i] = df[i]
			df.rename(columns={i: 'p_'+i}, inplace=True)
	return df

'''
FUNCTION: sortHaplotypes(df)
PURPOSE: sorts haplotypes by length
INPUT: dataframe
OUTPUT: sorted dataframe by haplotype length
'''
def sortHaplotypes(df):
	df = df.transpose()
	df['length'] = (df['end'] - df['start']) + 1
	excluded = df.head(3)
	included = df.iloc[3:]
	included = included.sort_values(by ='length', ascending=[False]) #longest at top
	df = pd.concat([excluded, included])
	df = df.reset_index()
	df = df.rename(columns = {'index': 'POS'})
	return df

# def main(df):
def main():
	print('*****STARTING HAPLOTYPES*****')
	global haplotypeFile
	filename = config.__FILEPATH__ + "cleaned_" + config.__FILENAME__
	distinctFile = config.__FILEPATH__ +  "distinct_" + config.__FILENAME__
	haplotypeFile = config.__FILEPATH__ +  "distinct_" + config.__FILENAME__

	# if already have counts
	fileCheck = Path(distinctFile)
	if fileCheck.is_file():
		return

	getHaplotypes(filename)
	return
