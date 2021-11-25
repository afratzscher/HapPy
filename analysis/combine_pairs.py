# tries to build reference from full length haplos
import os
import pandas as pd
from itertools import compress
import itertools
import math
import numpy as np

def getDF(file, numberrows = None):
	if numberrows is not None:
		print('start reading')
		df = pd.read_csv(file, sep="\t")
		# df = pd.read_csv(file, sep='\t', nrows = numberrows)
		print('done reading')
	else:
		df = pd.read_csv(file, sep="\t")
	df = df[df.columns[~df.isnull().any()]] # removes columns that dont have SNPs (keep columns with rs ID)
	newCols = df.iloc[0][:].tolist() # save ID row (to rename columns later)
	
	#save positions of SNPs (in case that have no overlapping SNPs)
	position = df.columns[df.columns.get_loc("identical")+1:df.columns.get_loc("start")]

	df = df[3:][:].reset_index() # removes REF and ALT and ID rows

	# df['counts'] = df['SAS'] # BY POPULATION
	if numberrows is None:
		df = df.sort_values(by='counts', ascending=False)
		df = df[df['counts'] != 0] 
	
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
	df.columns = newCols[beginidx:idx-1] # have to use SNP IDs instead of vals b/c sometimes stored as 1.0, others as 1
		
	if df.empty:
		return (None, None, None, None)
	return df, haploID, rank, position

def fetchDF(gene):
	file = direct + gene + "/full_length_haplotypes_" + gene + ".vcf"
	try:
		(df, haploID, rank, position) = getDF(file)
		if df is None:
			raise Exception
	except:
		print('checking other file')
		file = direct + gene + "/distinct_" + gene + ".vcf"
		(df, haploID, rank, position) = getDF(file) #, numberrows = 1000)
		print(len(df), 'length of df')
	return df, haploID, rank, position

def toFile(genes, res):
	idx = pd.DataFrame(res, columns = genes)

	df = pd.DataFrame()
	for gene in genes:
		print('gene fetch: ', gene)
		if (idx[gene] == '-').any():
			print(gene, 'no LCSH')
			break
		genedf, haploID = fetchDF(gene)
		temp = genedf.iloc[idx[gene]]
		newcols = [col for col in temp.columns if col not in df.columns]

		if gene == genes[0]:
			df = temp
		else:
			df = df.join(temp[newcols].reset_index(drop = True))

	df.reset_index(inplace = True, drop = True)
	print(df)
	
	# df.to_csv('LCSH_158781205-159712288.csv')

def run(genes, direction):
	if direction == 'upstream':
		genes = list(reversed(genes))
	pairs = list(zip(genes, genes[1:]))

	overlap = []
	for gene, secondgene in pairs:
		print(gene, secondgene)
		if gene == genes[0]:
			df, firstID, rank, position = fetchDF(gene)
		# get df for next gene
		nextDF, nextID, nextRANK, nextPOSITION = fetchDF(secondgene)

		# combine twice
		tmp = []

		mutuals = [col for col in df.columns if col in nextDF.columns]
	
		# add 10000 to index (b/c number of full length haplos < 10000 ) -> can differentiate between haplos from df and nextDF
		nextDF.index+=offset

		# combines df for both genes with same SNPs
		allSNPs = df[mutuals].append(nextDF[mutuals]) 
			#NOTE: df[mutuals] CAN have duplicates (b/c different outside this range, but only looking at specific range)
		
		#case where no overlap
		if not mutuals:
			if direction == 'downstream':
				if position[-1] < nextPOSITION[0]:
					print("NO MUTUALS BC NO SNPS BETWEEN, TRY ALL COMBINATIONS")
					idx = [tuple(df.index.tolist() + nextDF.index.tolist())]
			if direction == 'upstream':
				if position[-1] > nextPOSITION[0]:
					print("*****NO MUTUALS BC NO SNPS BETWEEN, TRY ALL COMBINATIONS*****")
					idx = [tuple(df.index.tolist() + nextDF.index.tolist())]	
		else:
			allSNPs = allSNPs[allSNPs.duplicated(keep=False)] # keeps only identical rows
			idx = allSNPs.groupby(list(allSNPs)).apply(lambda x: tuple(x.index)).tolist() # gets indices of identical rows
	
		# get combinations (e.g. make [1,2,3] to [1,2], [2,3], [1,3])
		lst = []
		for i in idx:
			val = list(itertools.combinations(i, 2))
			lst.extend(val)

		# only keep if overlapping match for DIFFERENT gene (NOT same gene)
		idx = []
		for i in lst:
			length = len(i)
			if (i[0] < offset) and (i[1] < offset): # both in gene
				continue
			elif (offset <= i[0]) and (offset <= i[1]): # both in next gene
				continue
			else:
				if i[0] >= offset:
					idx.append(((i[0]-offset, i[1]), (nextRANK[i[0]-offset], rank[i[1]]))) # do - so that get index in df for next gene
				if i[1] >= offset:
					idx.append(((i[0], i[1]-offset), (rank[i[0]], nextRANK[i[1]-offset]))) # do - so that get index in df for next gene
		overlap.append(idx)

		#remove offset
		nextDF.index-=offset
	
		# at end, go to next 
		df, firstID, rank, position = nextDF, nextID, nextRANK, nextPOSITION

	res = []
	for i in range(0, len(overlap)-1):
		temp = res
		print(i,i+1, 'ADDING ', genes[i+2])
		if i == 0:
			res = overlap[0]
			print(len(res), ' initial')
			res.sort(key = lambda x: x[0][-1])
			next_overlap = overlap[i+1]
			next_overlap.sort(key = lambda x: x[0][0])
			res = [[[u[0][0], u[0][1], v[0][1]], [u[1][0] + u[1][1] + v[1][1]]] 
				for u in res for v in next_overlap
				if u[0][-1] == v[0][0]]
		else:
			res.sort(key = lambda x: x[0][-1])
			next_overlap = overlap[i+1]
			next_overlap.sort(key = lambda x: x[0][0])
			res = [[u[0] + [v[0][1]], [u[1][0]+v[1][1]]] 
				for u in res for v in next_overlap
				if u[0][-1] == v[0][0]]
		if not res: # if res empty
			print('empty when adding: ', genes[i+2])
			res = temp
			print(len(res), " LCSH")
			break

		#sort res by rank
		res.sort(key = lambda x: x[1])

		#if > 10,000, keep only first 10,000
		if len(res) > threshold:
			if res[threshold-1][1] == res[threshold][1]: 
				tmp = [i for i in res if i[1] < res[threshold][1]]
			else: # otherwise, take first 10,000 instances
				tmp = res[:threshold]
				# if none left, keep all
			if len(tmp) == 0:
				print('cannot reduce, otherwise none left, NOTHING DONE')

			print('reduced from ', len(res), ' to ', len(tmp), ' for ', genes[i+2])
			res = tmp
		else:
			print('length: ', len(res))
	print(str(len(res)) + ' LCSHs\n\n')

def main():
	global direct
	global offset
	global maxval
	global threshold
	# direct = '/'.join(os.getcwd().split('/')[:-1]) + "/results/DONE/"
	direct = '/Volumes/AF_SSD/ACKR1-Algorithm/results/'
	offset = 10000 # to distinguish 2 df
	threshold = 10000 # 1 million takes too long, 500000 and 100,000 take too long as well
	numGenes = 2500
	
	
	genes = ['FAM72C', 'PPIAL4E', 'NBPF15', 'PPIAL4F']
	print(len(genes))
	run(genes, 'downstream')


if __name__ == '__main__':
	main()