#USING THIS ONE CURRENTLY
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

	# df['counts'] = df['EAS'] # BY POPULATION
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
	
	# RESULTS: #
	# 1.ANKRD35 - TXNIP (empty when adding HJV), 
	# 2.NBPF10-PPIAL4H (ends with NBPF12)
	# HAVE A LONG BREAK
	# 3.NBPF12 - PPIAL4H (empty when adding NBPF12)
	# ANOTHER LONG BREAK
	# 4.GJA5 - GPR89B (empty when adding NBPF11)
	# 5.NOTCH2NLC - H3C14 (empty when adding:  H2BC18)
	# SHORT BREAK (maybe because missing H3C15)
#NOTE: in this region where dont have breaks (but not multiple genes after), usually reading from full length haplos
# (instead of distinct haplos, which we do if no full length haplos)
	# 6.H4C15 - TUFT1 (empty when adding SNX27)
	# 7.SNX27 - NPR1 (empty when adding:  INTS3)
	# 8.INTS3 - TPM3 (empty because missing CFAP141)
	# 9.UBAP2L - ARHGEF2 (empty when adding:  SSR2)
	# 10.SSR2 - ARHGEF11 (empty when adding:  ETV3L)
	# 11.ETV3L - FCER1A (empty when adding:  OR10J1)
	# OR10J1 - SLAMF1 (empty when adding:  CD48)
	# CD48 - FCGR2B (empty when adding FCRLA)
	# FCRLA - NUF2 (empty when adding PBX1)
	# break PBX1, LMX1A (reading from distinct)
	# RXRG - TBX19 (empty when adding:  XCL2)
	# break, XCL2, XCL1, DPT
	# DPT - NMT2 (empty when adding: GORAB)
	# break GORAB
	# PRRX1 - FASLG (empty when adding: TNFSF18, reading from distinct)
	#break TNFSF18
	# TNFSF4 - RC3H1 (break because no RABGAP1L file)
	# break because no RABGAP1L file
	# CACYBP - KIAA0040 (empty when adding: TNR)
	# break - tnr
	# 20. COP1 - BRINP2 (break because no SEC16B file)
	# break - CRYZL2P-SEC16B, RASAL2 (no file),TEX35, RALGPS2
	# FAM20B - TDRD5 ( empty when adding:  FAM163A)
	# FAM163A - QSOX1 ( empty when adding:  ACBD6)
	#break ACBD6
	# XPR1 - IER5 (empty when adding:  CACNA1E)
	#break CACNA1E
	# ZNF648- TEDDM1 (empty when adding:  RGSL1)
	# RGSL1 - EDEM3 (empty when adding:  NIBAN1)
	# NIBAN1 - IVNS1ABP (empty when adding:  HMCN1)
	# HMCN1 - PTGS2 (too large: 'PLA2G4A')
	#break BRINP3, 'BRINP3', 'RGS18',
	# RGS21 - RGS13 empty when adding:  RGS2
	# RGS2 - CDC73 (empty when adding:  KCNT2)
	#break KCNT2
	# 30. CFH - NEK7 empty when adding:  ATP6V1G3
	# break ATP6V1G3, PTPRC, NR5A2, ZNF281
	# KIF14 - GPR25 empty when adding:  INAVA
	# INAVA - LAD1 (empty when adding:  TNNI1)
	# TNNI1 - GPR37L1 (empty when adding:  PTPN7)
	# PTPN7 - UBE2T (empty when adding:  PPP1R12B)
	# break PPP1R12B, KDM5B
	# RABIF - OPTC (empty when adding:  ATP2B4)
	# break ATP2B4 , LAX1, ZC3H11A, SNRPE
	# SOX13 - SLC45A3 (empty when adding:  NUCKS1)
	# NUCKS1 - CTSE (RHEX missing)
	# AVPR1B - PLXNA2 empty when adding:  CAMK1G
	# CAMK1G - SYT14 empty when adding:  SERTAD4
	# break - 'SERTAD4', 'HHAT', KCNH1(KCNH1 does not exist)
	# 40. RCOR3 - DTL (empty when adding:  PPP2R5A)
	# PPP2R5A - SMYD2 (CENPF error?)
	#BREAK KCNK2, KCTD3, USH2A, ESRRG, GPATCH2, SPATA17, RRP15, TGFB2
	# LYPLAL1 - HLX empty when adding:  DUSP10
	# break DUSP10
	#  HHIPL2- CAPN8 (empty when adding:  CAPN2)
	# break CAPN2,TP53BP2
	# FBXO28 - CNIH4 (empty when adding:  CNIH3)
	# BREAK CNIH3, DNAH14, LBR, ENAH
	# SRP9 - PARP1 (empty when adding:  STUM)
	# break STUM,ITPKB 
	# PSEN2 - CDC42BPA empty when adding:  ZNF678
	# ZNF678 - WNT9A (empty when adding:  WNT3A)
	# WNT3A - RHOU empty when adding:  RAB4A
	# RAB4A - GALNT2 empty when adding:  PGBD5
	# PGBD5 - TSNAX empty when adding:  DISC1
	# break DISC1
	# SIPA1L2 - PCNX2 empty when adding:  MAP3K21
	# MAP3K21 - IRF2BP2 empty when adding:  TOMM20
	# TOMM20 - ARID4B empty when adding:  GNG4
	# BREAK = GNG4,LYST,NID1
	# GPR137B - HEATR1 empty when adding:  ACTN2
	# ACTN2 - MT1HL1 empty when adding:  RYR2
	#  break RYR2', 'ZP4', 'CHRM3', FMN2 error, GREM2
	# RGS7-EXO1 empty when adding:  BECN2
	#### break BECN2, MAP1LC3C, PLD5, AKT3, ZBTB18
	# C1orf100 - HNRNPU (empty when adding:  EFCAB2)
	# break EFCAB2, KIF26B, SMYD3 (missing)
	# TFB2M - SCCPDH (empty when adding:  AHCTF1)
	# AHCTF1 - OR2T6 (empty OR2T1
	 # break or2t1, missing OR2T7, ORT2T2
	 # ORT2T2 - OR2T27 (empty when adding:  OR14I1)
	 # break OR14I1,LYPD8
	 # SH3BP5L - PGBD2


	genes = ['SH3BP5L', 'ZNF672', 'ZNF692', 'PGBD2']
	print(len(genes))
	genes = genes[0:60]
	print(len(genes))
	run(genes, 'downstream')


if __name__ == '__main__':
	main()