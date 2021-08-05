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
	
	#save positions of SNPs (in case that have no overlapping SNPs)
	position = df.columns[df.columns.get_loc("identical")+1:df.columns.get_loc("start")]

	df = df[3:][:].reset_index() # removes REF and ALT and ID rows

	df = df.sort_values(by='counts', ascending=False)
	
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
		print("LENGTH ISSUE")
	df.columns = newCols[beginidx:idx-1] # have to use SNP IDs instead of vals b/c sometimes stored as 1.0, others as 1

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
	if direction == 'downstream':
		triplets = list(tuple(genes[i:i+3]) for i in range(len(genes)-1) if i%2==0)
	if direction == 'upstream':
		genes = list(reversed(genes))
		triplets = list(tuple(genes[i:i+3]) for i in range(len(genes)-1) if i%2==0)
	
	#in case DONT have multiple of 3, add last gene
	last = None
	if len(triplets[-1]) == 2:
		last = list(triplets[-1])
		triplets = triplets[:-1]

	overlap = []
	for gene, secondgene, thirdgene in triplets:
		print(gene, secondgene, thirdgene)
		if gene == genes[0]:
			firstDF, firstID, firstRANK, firstPOSITION = fetchDF(gene)
		# get df for next gene
		secondDF, secondID, secondRANK, secondPOSITION = fetchDF(secondgene)
		thirdDF, thirdID, thirdRANK, thirdPOSITION = fetchDF(thirdgene)

		# combine twice
		tmp = []
		for i in range(0,2):
			if i == 0:
				df = firstDF
				nextDF = secondDF
				rank = firstRANK
				nextRANK = secondRANK
				position = firstPOSITION
				nextPOSITION = secondPOSITION
			else:
				df = secondDF
				nextDF = thirdDF
				rank = secondRANK
				nextRANK = thirdRANK
				position = secondPOSITION
				nextPOSITION = thirdPOSITION

			mutuals = [col for col in df.columns if col in nextDF.columns]
		
			# add 10000 to index (b/c number of full length haplos < 10000 ) -> can differentiate between haplos from df and nextDF
			nextDF.index+=offset

			# combines df for both genes with same SNPs
			allSNPs = df[mutuals].append(nextDF[mutuals]) 
				#NOTE: df[mutuals] CAN have duplicates (b/c different outside this range, but only looking at specific range)
			
			#case where no overlap
			if not mutuals:
				if position[-1] < nextPOSITION[0]:
					print("NO MUTUALS BC NO SNPS BETWEEN, TRY ALL COMBINATIONS")
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
			tmp.append(idx)

			#remove offset
			nextDF.index-=offset
		
		#combine duplets into triplets
		overlap.append([((u[0][0], u[0][1], v[0][1]), (u[1][0], u[1][1],v[1][1])) for u in tmp[0] for v in tmp[1]
				if u[0][1] == v[0][0]])
	
		# at end, go to next 
		firstDF, firstID, firstRANK, firstPOSITION = thirdDF, thirdID, thirdRANK, thirdPOSITION

	# run on single double (if not odd number of genes)
	if last:
		print(last[0], last[1], 'last')
		df, haploID, rank, position = fetchDF(last[0])
		nextDF, nextID, nextRANK, nextPOSITION = fetchDF(last[1])
		mutuals = [col for col in df.columns if col in nextDF.columns]
		
		# add 10000 to index (b/c number of full length haplos < 10000 ) -> can differentiate between haplos from df and nextDF
		nextDF.index+=offset
		# combines df for both genes with same SNPs
		allSNPs = df[mutuals].append(nextDF[mutuals]) 
			#NOTE: df[mutuals] CAN have duplicates (b/c different outside this range, but only looking at specific range)
			
		#case where no overlap
		if not mutuals:
			if position[-1] < nextPOSITION[0]:
				print("NO MUTUALS BC NO SNPS BETWEEN, TRY ALL COMBINATIONS")
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
		#add to end of list 
		overlap.append(idx)

		#remove offset
		nextDF.index-=offset

	res = []
	for i in range(0, len(overlap)-1):
		val = 2*i
		temp = res
		if (len(overlap[i+1][0][0])) == 3:
			print(val,val+1,val+2, 'ADDING ', genes[val+3], genes[val+4])
		else:
			print(val,val+1, 'ADDING ', genes[val+3])
		if i == 0:
			res = overlap[0]
			res.sort(key = lambda x: x[0][-1])
			next_overlap = overlap[i+1]
			next_overlap.sort(key = lambda x: x[0][0])
			res = [[[u[0][0], u[0][1], v[0][0], v[0][1], v[0][2]], [u[1][0] + u[1][1] + u[1][2] + v[1][1] + v[1][2]]] for u in res for v in next_overlap
				if u[0][-1] == v[0][0]]
		else:
			res.sort(key = lambda x: x[0][-1])
			next_overlap = overlap[i+1]
			next_overlap.sort(key = lambda x: x[0][0])
			if (len(overlap[i+1][0][0])) == 3:
				res = [[u[0] + [v[0][1], v[0][2]], [u[1][0]+v[1][1]+v[1][2]]] for u in res for v in next_overlap
					if u[0][-1] == v[0][0]]
			else:
				res = [[u[0] + [v[0][1]], [u[1][0]+v[1][1]]] for u in res for v in next_overlap
					if u[0][-1] == v[0][0]]
		if not res: # if res empty
			print('empty when adding: ', genes[val+3])
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

			if (len(overlap[i+1][0][0])) == 3:
				print('reduced from ', len(res), ' to ', len(tmp), ' for ', genes[val+3], genes[val+4])
			else:
				print('reduced from ', len(res), ' to ', len(tmp), ' for ', genes[val+3])
			res = tmp
		else:
			print('length: ', len(res))
	print(str(len(res)) + ' LCSHs\n\n')
	# print(res)


	# if direction == 'upstream':
	# 	# FOR UPSTREAM GENE ADDED -> FOCUS ON THIS FOR NOW
	# 	res = []
	# 	for i in reversed(range(1, len(overlap))):
	# 		temp = res
	# 		val = 2*i
	# 		if (i == 1) and len(overlap[-1][0][0]) == 2:
	# 			print(val-1,val-2, 'ADDING ', genes[val-2])
	# 		else:
	# 			print(val-1,val-2,val-3, 'ADDING ', genes[val-2], genes[val-3])
		
	# 		if i == (len(overlap)-1):
	# 			res = overlap[-1]
	# 			res.sort(key = lambda x: x[0][-1])
	# 			next_overlap = overlap[i-1]
	# 			next_overlap.sort(key = lambda x: x[0][0])
	# 			res = [[[v[0][0], v[0][1], u[0][0], u[0][1], u[0][2]], 
	# 				[u[1][0] + u[1][1] + u[1][2] + v[1][1] + v[1][0]]] 
	# 				for u in res for v in next_overlap
	# 				if u[0][0] == v[0][-1]]
	# 		else:
	# 			res.sort(key = lambda x: x[0][-1])
	# 			next_overlap = overlap[i-1]
	# 			next_overlap.sort(key = lambda x: x[0][0])
	# 			res = [[[v[0][0]] + u[0], [u[1][0]+v[1][0]]] for u in res for v in next_overlap
	# 				if u[0][0] == v[0][1]]
	# 			#triplet at end
	# 			if (i == 1) and len(overlap[-1][0][0]) == 3:
	# 				res = [[[v[0][0], v[0][1]] + u[0], [u[1][0]+v[1][1]+v[1][0]]] 
	# 					for u in res for v in next_overlap
	# 					if u[0][0] == v[0][-1]]
	# 			else: # pair at end
	# 				res = [[[v[0][0]] + u[0], [u[1][0]+v[1][0]]] 
	# 					for u in res for v in next_overlap
	# 					if u[0][0] == v[0][-1]]

	# 		if not res: # if res empty
	# 			print('empty when adding: ', genes[i-1])
	# 			res = temp
	# 			print(len(res), " LCSH")
	# 			break

	# 		#sort res by rank
	# 		res.sort(key = lambda x: x[1])

	# 		#if > 10,000, keep only first 10,000
	# 		if len(res) > threshold:
	# 			#if #10,000 has rank of 2 and so does 9999, remove BOTH
	# 			if res[threshold-1][1] == res[threshold][1]: 
	# 				tmp = [i for i in res if i[1] < res[threshold][1]]
	# 			else: # otherwise, take first 10,000 instances
	# 				tmp = res[:threshold]
	# 			if len(tmp) == 0:
	# 				print('cannot reduce, otherwise none left')
	# 			if (len(overlap[i+1][0][0])) == 3:
	# 				print('reduced from ', len(res), ' to ', len(tmp), ' for ', genes[i-1], genes[i-2])
	# 			else:
	# 				print('reduced from ', len(res), ' to ', len(tmp), ' for ', genes[i-1])
	# 			res = tmp
	# 		else:
	# 			print('length: ', len(res))
	# 	print(str(len(res)) + " LCSH")
	# 	print(res)

	# toFile(genes, res)

def main():
	global direct
	global offset
	global maxval
	global threshold
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/results/"
	offset = 10000 # to distinguish 2 df
	threshold = 10000 # 1 million takes too long, 500000 and 100,000 take too long as well
	
	#TRIPLETS
	#101 genes
	genes  = ['SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1', 'CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 'ACKR1', 'FCER1A', 'OR10J1', 'OR10J5', 'APCS', 'CRP', 'DUSP23', 'FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1',
	'CD48', 'SLAMF7', 'LY9', 'CD244', 'ITLN1', 'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1']
	run(genes, 'downstream')

	#DOUBLES
	#101 genes
	# genes  = ['SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1', 'CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 'ACKR1', 'FCER1A']
	# genes = ['OR10J1', 'OR10J5', 'APCS', 'CRP', 'DUSP23', 'FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1']
	# genes = ['CD48', 'SLAMF7', 'LY9', 'CD244', 'ITLN1', 'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1']
	# run(genes, 'downstream')

	# genes  = ['SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1']
	# genes = ['CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 'ACKR1', 'FCER1A']
	# genes = ['OR10J1', 'OR10J5', 'APCS', 'CRP', 'DUSP23']
	# genes = ['FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1', 'CD48']
	# genes = ['SLAMF7', 'LY9', 'CD244', 'ITLN1', 'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1']
	# run(genes, 'upstream')

	# test set = 150 genes (75 + ACKR1 + 74)
	# genes = ['MUC1', 'THBS3', 'GBA', 'FAM189B', 'SCAMP3', 
	# 		'CLK2', 'HCN3', 'FDPS', 'RUSC1', 'ASH1L', 
	# 		'MSTO1', 'DAP3', 'GON4L', 'SYT11', 'RIT1', 
	# 		'KHDC4', 'RXFP4', 'ARHGEF2']
	# genes = ['SSR2', 'UBQLN4', 
	# 		'LAMTOR2', 'RAB25', 'MEX3A', 'LMNA', 'SEMA4A', 
	# 		'SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 
	# 		'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 
	# 		'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 
	# 		'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 
	# 		'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 
	# 		'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1', 
	# 		'CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 
	# 		'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 
	# 		'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 
	# 		'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3',
	# 		'ACKR1',
	# 		'FCER1A']
	# genes = ['OR10J1', 'OR10J5', 'APCS', 'CRP', 
	# 		'DUSP23', 'FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 
	# 		'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 
	# 		'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 
	# 		'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 
	# 		'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1']
	# genes = ['CD48', 'SLAMF7', 'LY9', 'CD244', 'ITLN1', 
	# 		'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 
	# 		'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 
	# 		'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1', 
	# 		'MPZ', 'SDHC', 'FCGR2A', 'HSPA6', 'FCGR3A', 
	# 		'FCGR3B', 'FCGR2B']
	# genes = ['FCRLA', 'FCRLB', 'DUSP12', 
	# 		'ATF6', 'OLFML2B', 'NOS1AP', 'SPATA46', 'C1orf226', 
	# 		'SH2D1B', 'UHMK1', 'UAP1', 'DDR2', 'HSD17B7', 
	# 		'CCDC190', 'RGS4', 'RGS5', 'NUF2']
	# run(genes, 'downstream')

	# genes = ['MUC1', 'THBS3', 'GBA', 'FAM189B', 'SCAMP3', 
	# 		'CLK2', 'HCN3', 'FDPS', 'RUSC1', 'ASH1L', 
	# 		'MSTO1', 'DAP3', 'GON4L', 'SYT11', 'RIT1', 
	# 		'KHDC4', 'RXFP4', 'ARHGEF2'] 
	# genes = ['SSR2', 'UBQLN4', 
	# 		'LAMTOR2', 'RAB25', 'MEX3A', 'LMNA', 'SEMA4A', 
	# 		'SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 
	# 		'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 
	# 		'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 
	# 		'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 
	# 		'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 
	# 		'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1'] 
	# genes = ['CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 
	# 		'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 
	# 		'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 
	# 		'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3',
	# 		'ACKR1',
	# 		'FCER1A']
	# genes = ['OR10J1', 'OR10J5', 'APCS', 'CRP', 'DUSP23']
	# genes = ['FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 
	# 		'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 
	# 		'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 
	# 		'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 
	# 		'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1',
	# 		'CD48']
	# genes = ['SLAMF7', 'LY9', 'CD244', 'ITLN1', 
	# 		'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 
	# 		'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 
	# 		'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1', 
	# 		'MPZ', 'SDHC', 'FCGR2A', 'HSPA6', 'FCGR3A', 
	# 		'FCGR3B', 'FCGR2B']
	# genes = ['FCRLA', 'FCRLB', 'DUSP12', 
	# 		'ATF6', 'OLFML2B', 'NOS1AP', 'SPATA46', 'C1orf226', 
	# 		'SH2D1B', 'UHMK1', 'UAP1', 'DDR2', 'HSD17B7', 
	# 		'CCDC190', 'RGS4', 'RGS5']
	# run(genes, 'upstream')


if __name__ == '__main__':
	main()