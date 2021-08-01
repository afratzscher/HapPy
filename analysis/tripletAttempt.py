# tries to build reference from full length haplos
import os
import pandas as pd
from itertools import compress
import itertools
import math
import numpy as np

def fetchDF(gene):
	file = direct + gene + "/full_length_haplotypes_" + gene + ".vcf"
	df = pd.read_csv(file, sep="\t")
	df = df[df.columns[~df.isnull().any()]] # removes columns that dont have SNPs (keep columns with rs ID)
	newCols = df.iloc[0][:].tolist() # save ID row (to rename columns later)
	df = df[3:][:].reset_index() # removes REF and ALT and ID rows

	df = df.sort_values(by='counts', ascending=False)
	# get 25 most frequent full length haplotypes
	# df = df.iloc[0:25, :]	
	#ranks haplos based on counts -> if multiple have same, do average
	df['rank'] = df['counts'].rank(ascending=False) 

	idx = df.columns.get_loc("start")
	beginidx = df.columns.get_loc("identical")
	# rank = df.iloc[:, df.columns.get_loc("rank")]
	rank= np.sqrt(df['rank'])
	haploID = df.iloc[:, df.columns.get_loc("haplotypeID")]
	df = df.iloc[:, beginidx+1:idx] # keeps only columns with SNV data
	df.columns = newCols[beginidx:idx-1] # have to use SNP IDs instead of vals b/c sometimes stored as 1.0, others as 1

	return df, rank

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
	pairs = list(zip(genes, genes[1:]))
	pairs = list(zip(genes, genes[1:], genes[2:]))
		
	overlap = []
	# get for first gene
	for first, second, third in pairs:
		print(first, second, third)
		if first == genes[0]:
			firstDF, firstRANK = fetchDF(first)
			secondDF, secondRANK = fetchDF(second)
		# get df for next gene (third)
		thirdDF, thirdRANK = fetchDF(third)

		tmp = []
		for i in range(0,2):
			if i == 0:
				df = firstDF
				nextDF = secondDF
				rank = firstRANK
				nextRANK = secondRANK
			else:
				df = secondDF
				nextDF = thirdDF
				rank = secondRANK
				nextRANK = thirdRANK
			mutuals = [col for col in df.columns if col in nextDF.columns]
			
			# add 10000 to index (b/c number of full length haplos < 10000 ) -> can differentiate between haplos from df and nextDF
			nextDF.index+=offset

			# combines df for both genes with same SNPs
			allSNPs = df[mutuals].append(nextDF[mutuals]) 
				#NOTE: df[mutuals] CAN have duplicates (b/c different outside this range, but only looking at specific range)
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
			tmp.append(idx)

			#remove offset
			nextDF.index-=offset
		
		#combine into triplets
		tmp2 = [((u[0][0], u[0][1], v[0][1]), (u[1][0], u[1][1], v[1][1]) ) for u in tmp[0] for v in tmp[1] if u[0][1] == v[0][0]]
		overlap.append(tmp2)

		# at end, set to
		firstDF, firstRANK = secondDF, secondRANK
		secondDF, secondRANK = thirdDF, thirdRANK

	if direction == 'downstream':
		# FOR DOWNSTREAM GENE ADDED
		res = []
		for i in range(0, len(overlap)-1):
			temp = res
			print(i,i+1, 'ADDING ', genes[i+2])
			if i == 0:
				res = overlap[0]
				res.sort(key = lambda x: x[0][-2])
				next_overlap = overlap[i+1]
				next_overlap.sort(key = lambda x: x[0][0])
				res = [[[u[0][0], u[0][1], u[0][2], v[0][2]], [u[1][0] + u[1][1] + u[1][2] + v[1][2]]] for u in res for v in next_overlap
					if u[0][1] == v[0][0] and u[0][2] == v[0][1]]
			else:
				res.sort(key = lambda x: x[0][-2])
				next_overlap = overlap[i+1]
				next_overlap.sort(key = lambda x: x[0][0])
				res = [[u[0] + [v[0][2]], [u[1][0]+v[1][2]]] for u in res for v in next_overlap
					if u[0][-2] == v[0][0] and u[0][-1] == v[0][1]]
			if not res: # if res empty
				print('empty when adding: ', genes[i+3])
				res = temp
				print(len(res), " LCSH")
				break

			#sort res by rank
			res.sort(key = lambda x: x[1])

			#if > 10,000, keep only first 10,000
			if len(res) > threshold:
				tmp = res[:threshold]
				print('reduced from ', len(res), ' to ', len(tmp), ' for ', genes[i+2])
				res = tmp
			else:
				print('length: ', len(res))
		print(str(len(res)) + ' LCSHs\n\n')
		# print(res)

	if direction == 'upstream':
		# FOR UPSTREAM GENE ADDED -> FOCUS ON THIS FOR NOW
		res = []
		for i in reversed(range(1, len(overlap))):
			temp = res
			print(i,i-1, 'ADDING ', genes[i-1])
			if i == (len(overlap)-1):
				res = overlap[-1]
				res = [[[v[0][0], u[0][0], u[0][1]], [v[1][0] + u[1][0] + u[1][1]]] for u in res for v in overlap[i-1]
					if u[0][0] == v[0][1]]
			else:
				res = [[[v[0][0]] + u[0], [u[1][0]+v[1][0]]] for u in res for v in overlap[i-1]
					if u[0][0] == v[0][1]]
			if not res: # if res empty
				print('empty when adding: ', genes[i-1])
				res = temp
				print(len(res), " LCSH")
				break

			#sort res by rank
			res.sort(key = lambda x: x[1])

			#if > 10,000, keep only first 10,000
			if len(res) > threshold:
				#if #10,000 has rank of 2 and so does 9999, remove BOTH
				if res[threshold-1][1] == res[threshold][1]: 
					tmp = [i for i in res if i[1] < res[threshold][1]]
				else: # otherwise, take first 10,000 instances
					tmp = res[:threshold]
				tmp = res[:threshold]
				print('reduced from ', len(res), ' to ', len(tmp), ' for ', genes[i-1])
				res = tmp
			else:
				print('length: ', len(res))
		print(str(len(res)) + " LCSH")
		# print(res)

	# toFile(genes, res)

def main():
	global direct
	global offset
	global maxval
	global threshold
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/results/"
	offset = 10000 # to distinguish 2 df
	threshold = 10000 # 1 million takes too long, 500000 and 100,000 take too long as well

	# ALL GENES UPSTREAM OF ACKR1 + ACKR1 + 4 GENES DOWNSTREAM
	genes = ['ATP8B2', 'IL6R', 'UBE2Q1', 'ADAR', 'KCNN3', 
		 'PMVK', 'PBXIP1', 'PYGO2', 'SHC1', 'CKS1B', 
		 'FLAD1', 'LENEP', 'ZBTB7B', 'DCST1', 'ADAM15', 
		 'EFNA4', 'EFNA3', 'EFNA1', 'SLC50A1', 'DPM3', 
		 'KRTCAP2', 'TRIM46', 'MUC1', 'THBS3', 'GBA', 
		 'FAM189B', 'SCAMP3', 'CLK2', 'HCN3', 'FDPS', 
		 'RUSC1', 'ASH1L', 'MSTO1', 'DAP3', 'GON4L', 
		 'SYT11', 'RIT1', 'KHDC4', 'RXFP4', 'ARHGEF2', 
		 'SSR2', 'UBQLN4', 'LAMTOR2', 'RAB25', 'MEX3A', 
		 'LMNA', 'SEMA4A', 'SLC25A44', 'PMF1-BGLAP', 'GLMP', 
		 'VHLL', 'CCT3', 'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 
		 'NAXE', 'HAPLN2', 'BCAN', 'NES', 'CRABP2', 
		 'METTL25B', 'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 
		 'PEAR1', 'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 
		 'FCRL4', 'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 
		 'KIRREL1', 'CD1D', 'CD1A', 'CD1B', 'CD1E', 
		 'OR10T2', 'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 
		 'OR6P1', 'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 
		 'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2',
	'CADM3', 'ACKR1', 'FCER1A', 'OR10J1', 'OR10J5', 'APCS', 'CRP'] 

	run(genes, 'downstream')


if __name__ == '__main__':
	main()