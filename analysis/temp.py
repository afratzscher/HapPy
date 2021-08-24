# tries to build reference from full length haplos
import os
import pandas as pd
from itertools import compress
import itertools
import math
import numpy as np

def fetchDF(gene):
	print('here')
	file = direct + gene + "/distinct_" + gene + ".vcf"
			
	df = pd.read_csv(file, sep="\t")
	print('done reading')
	df = df[df.columns[~df.isnull().any()]] # removes columns that dont have SNPs (keep columns with rs ID)
	newCols = df.iloc[0][:].tolist() # save ID row (to rename columns later)
	
	#save positions of SNPs (in case that have no overlapping SNPs)
	position = df.columns[df.columns.get_loc("identical")+1:df.columns.get_loc("start")]

	df = df[3:][:].reset_index() # removes REF and ALT and ID rows

	# df['counts'] = df['SAS'] # POPULATION
	df = df.sort_values(by='counts', ascending=False)
	df = df[df['counts'] != 0] #POPULATION
	
	#ranks haplos based on counts -> if multiple have same, do average
	df['rank'] = df['counts'].rank(ascending=False) 
	rank = df.iloc[:, df.columns.get_loc("rank")]
	# rank= np.sqrt(df['rank'])
	#FREQUENCY INSTEAD OF RANK
	# rank = df['counts'] /  df['counts'].sum()

	beginidx = df.columns.get_loc("identical")
	idx = df.columns.get_loc("start")
	haploID = df.iloc[:, df.columns.get_loc("haplotypeID")]
	df = df.iloc[:, beginidx+1:idx] # keeps only columns with SNV data
	if len(df) == 0:
		print('length issue')
		
	
	print("DONE")
	print(df)
	exit()

def run(genes, direction):
	fetchDF(genes[0])
	print('out')

def main():
	global direct
	global offset
	global maxval
	global threshold
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/results/DONE/"
	genes = ['PPIAL4E', 'NBPF15']
	run(genes, 'downstream')

if __name__ == '__main__':
	main()