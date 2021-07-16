# tries to build reference from full length haplos
# TO DO: decide on EITHER PYDC5 OR IFI16 
# CHECK WITH WF AND KS ON WHICH 11 GENES TO USE!
import os
import pandas as pd
from itertools import compress
import itertools

def fetchDF(gene):
	file = direct + gene + "/full_length_haplotypes_" + gene + ".vcf"
	df = pd.read_csv(file, sep="\t")
	df = df[df.columns[~df.isnull().any()]] # removes columns that dont have SNPs (keep columns with rs ID)
	df = df[3:][:].reset_index()
	del df['index']
	idx = df.columns.get_loc("start")
	beginidx = df.columns.get_loc("identical")
	haploID = df.iloc[:, df.columns.get_loc("haplotypeID")]
	df = df.iloc[:, beginidx+1:idx] # keeps only columns with SNV data
	return df, haploID

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


def main():
	global direct
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/results/"
	# genes = ['MNDA', 'PYHIN1', 'AIM2', 'CADM3', 'ACKR1', 'FCER1A', 'OR10J1', 
	#  'OR10J5', 'APCS', 'CRP']

	# genes = ['OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 'OR10X1', 'SPTA1',
		# 'OR6K2', 'OR6K3', 'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3',
		#  'ACKR1', 'FCER1A', 'OR10J1', 'OR10J5', 'APCS', 'CRP']

	# genes = ['OR10K1', 'OR10R2', 'OR6Y1', 

	genes = ['OR10K1', 'OR10R2', 'OR6Y1',
	'OR6P1', 'OR10X1', 'SPTA1','OR6K2', 'OR6K3',
	'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2'] # 'CADM3',
	 # 'ACKR1', 'FCER1A', 'OR10J1', 'OR10J5', 'APCS', 'CRP']



	pairs = list(zip(genes, genes[1:]))

	overlap = []
	# get for first gene
	for gene, nextgene in pairs:
		print(gene, nextgene)
		if gene == genes[0]:
			df, geneID = fetchDF(gene)
		# get df for next gene
		nextDF, nextID = fetchDF(nextgene)

		''' 
		block of code below shows that, for a gene, there are NO duplciated full length haplotypes
		this means that all duplicates in combined df have same conserved region for both genes
		'''
		# df = df[df.duplicated(keep=False)] # proves NO duplicated full-length haplotypes
		# nextDF = nextDF[nextDF.duplicated(keep=False)]

		mutuals = [col for col in df.columns if col in nextDF.columns]
		dfLen = df[mutuals].shape[0]
		nextLen = nextDF[mutuals].shape[0]
		maxdf = df[mutuals].shape[0] - 1
		maxnextdf = nextDF[mutuals].shape[0] + maxdf
		allSNPs = df[mutuals].append(nextDF[mutuals]).reset_index(drop = True) # combines df for both genes with same SNPs

		#NOTE: df[mutuals] CAN have duplicates (b/c different outside this range, but only looking at specific range)
		#-> NEED TO ACCOUNT FOR THIS 

		allSNPs = allSNPs[allSNPs.duplicated(keep=False)] # keeps only identical rows
		idx = allSNPs.groupby(list(allSNPs)).apply(lambda x: tuple(x.index)).tolist() # gets indices of identical rows

		# get combinations (e.g. make [1,2,3] to [1,2], [2,3], [1,3])
		lst = []
		for i in idx:
			val = list(itertools.combinations(i, 2))
			lst.extend(val)

		# remove if overlapping match for same gene
		idx = []
		for i in lst:
			length = len(i)
			if (i[0] <= maxdf) and (i[1] <= maxdf): # both in gene
				continue
			elif (maxdf < i[0] <= maxnextdf) and (maxdf < i[1] <= maxnextdf): # both in next gene
				continue
			else:
				idx.append((i[0], i[1]-maxdf-1)) # do - so that get index in df for next gene

		overlap.append(idx)
		# at end, set df to " nextgene"
		df, geneID = nextDF, nextID

	
	print('done getting vals')
	res = []
	for i in range(0, len(overlap)-1):
		temp = res
		print(i,i+1)
		# print(res)
		if i == 0:
			res = overlap[0]
			# print(res)
			res = [[u[0], u[1], v[1]] for u in res for v in overlap[i+1]
				if u[1] == v[0]]
			# print(type(overlap[0][1][1]))
			test = (res[0] + [(overlap[0][1][1])])
			# print(test[-1])
			# print(test)
			# print(overlap[0][1])
			# break
		else:
			res = [u + [v[1]] for u in res for v in overlap[i+1]
			if u[-1] == v[0]]
		# print(res)	
		if not res: # if res empty
			res = [u +  ["-"] for u in temp]
			print('empty')

	# for i in reversed(range(1, len(overlap))):
	# 	temp = res
	# 	print(i,i-1)
	# 	if i == (len(overlap)-1):
	# 		res = overlap[-1]
	# 	res = [(v[0], u) for u in res for v in overlap[i-1]
	# 		if u[0] == v[1]]
	# 	if not res: # if res empty
	# 		res = [("-", u) for u in temp]
	# 		print('empty')
	
	# print(res)
	# # print(len(res))
	print('starting to file')
	print(len(res))
	toFile(genes, res)

if __name__ == '__main__':
	main()