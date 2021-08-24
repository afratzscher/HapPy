# tries to build reference from full length haplos
import os
import pandas as pd
from itertools import compress
import itertools
import math
import numpy as np

def fetchDF(gene):
	full = True
	while True:
		if full:
			file = direct + gene + "/full_length_haplotypes_" + gene + ".vcf"
		else:
			file = direct + gene + "/haplotypes_" + gene + ".vcf"
			print('here')
		df = pd.read_csv(file, sep="\t")
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
			if not full:
				print("STILL LENGTH ISSUE")
				exit()
			full = False
		else:
			df.columns = newCols[beginidx:idx-1] # have to use SNP IDs instead of vals b/c sometimes stored as 1.0, others as 1
			break
	print("DONE")
	print(df)
	exit()
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
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/results/DONE/"
	offset = 10000 # to distinguish 2 df
	threshold = 10000 # 1 million takes too long, 500000 and 100,000 take too long as well
	
	#PAIRS
	#101 genes 
	#DOWNSTREAM


	genes = ['FAM72C', 'PPIAL4E', 'NBPF15', 'PPIAL4F', 'SRGAP2B', 'LOC101929805', 'PPIAL4D', 'NBPF20', 'GPR89A', 'ANKRD35', 'ITGA10', 'PEX11B', 'RBM8A', 'LIX1L', 'POLR3GL', 'TXNIP', 'HJV', 'NBPF10', 'NOTCH2NLA', 'PPIAL4H', 'NBPF12', 'CHD1L', 'GJA5', 'GJA8', 'GPR89B', 'NBPF11', 'PPIAL4G', 'NBPF14', 'NOTCH2NLB', 'NUDT4B', 'PDE4DIP', 'NBPF9', 'NOTCH2NLC', 'NBPF19', 'PPIAL4C', 'H2BC18', 'H3C13', 'H4C14', 'H3C14', 'H2AC18', 'H2AC19', 'H3C15', 'H4C15', 'H2BC21', 'H2AC20', 'H2AC21', 'BOLA1', 'SV2A', 'SF3B4', 'OTUD7B', 'VPS45', 'PLEKHO1', 'ANP32E', 'APH1A', 'CIART', 'MRPS21', 'PRPF3', 'RPRD2', 'TARS2', 'ECM1', 'ADAMTSL4', 'MCL1', 'ENSA', 'GOLPH3L', 'HORMAD1', 'CTSS', 'CTSK', 'ARNT', 'CTXND2', 'SETDB1', 'CERS2', 'ANXA9', 'MINDY1', 'PRUNE1', 'BNIPL', 'CDC42SE1', 'MLLT11', 'GABPB2', 'SEMA6C', 'LYSMD1', 'VPS72', 'PIP5K1A', 'PSMD4', 'PI4KB', 'RFX5', 'SELENBP1', 'PSMB4', 'POGZ', 'CGN', 'TUFT1', 'SNX27', 'CELF3', 'RIIAD1', 'TDRKH', 'LINGO4', 'RORC', 'C2CD4D', 'THEM5', 'THEM4', 'S100A10', 'S100A11', 'LOC100131107', 'TCHHL1', 'TCHH', 'RPTN', 'HRNR', 'FLG', 'FLG2', 'CRNN', 'LCE5A', 'CRCT1', 'LCE3E', 'LCE3D', 'LCE3C', 'LCE3B', 'LCE3A', 'LCE2D', 'LCE2C', 'LCE2B', 'LCE2A', 'LCE4A', 'C1orf68', 'KPRP', 'LCE1F', 'LCE1E', 'LCE1D', 'LCE1C', 'LCE1B', 'LCE1A', 'LCE6A', 'SMCP', 'IVL', 'SPRR4', 'SPRR1A', 'SPRR3', 'SPRR1B', 'SPRR2D', 'SPRR2A', 'SPRR2B', 'SPRR2E', 'SPRR2F', 'SPRR2G', 'LELP1', 'PRR9', 'LORICRIN', 'PGLYRP3', 'PGLYRP4', 'S100A9', 'S100A12', 'S100A8', 'S100A7', 'S100A6', 'S100A5', 'S100A4', 'S100A3', 'S100A2', 'S100A16', 'S100A14', 'S100A13', 'NPR1', 'INTS3', 'SLC27A3', 'GATAD2B', 'DENND4B', 'CRTC2', 'SLC39A1', 'RAB13', 'RPS27', 'NUP210L', 'TPM3', 'C1orf189', 'UBAP2L', 'HAX1', 'AQP10', 'ATP8B2', 'IL6R', 'UBE2Q1', 'ADAR', 'KCNN3', 'PMVK', 'PBXIP1', 'PYGO2', 'SHC1', 'CKS1B', 'FLAD1', 'LENEP', 'ZBTB7B', 'DCST1', 'ADAM15', 'EFNA4', 'EFNA3', 'EFNA1', 'SLC50A1', 'DPM3', 'KRTCAP2', 'TRIM46', 'MUC1', 'THBS3', 'GBA', 'FAM189B', 'SCAMP3', 'CLK2', 'HCN3', 'FDPS', 'RUSC1', 'ASH1L', 'MSTO1', 'DAP3', 'GON4L', 'SYT11', 'RIT1', 'KHDC4', 'RXFP4', 'ARHGEF2', 'SSR2', 'UBQLN4', 'LAMTOR2', 'RAB25', 'MEX3A', 'LMNA', 'SEMA4A', 'SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1', 'CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 'ACKR1', 'FCER1A', 'OR10J1', 'OR10J5', 'APCS', 'CRP', 'DUSP23', 'FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1', 'CD48', 'SLAMF7', 'LY9', 'CD244', 'ITLN1', 'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1']
	genes = ['PPIAL4E', 'NBPF15']
	run(genes, 'downstream')

	# genes = [ 'S100A8', 'S100A7', 'S100A6', 'S100A5', 'S100A4', 'S100A3', 'S100A2', 'S100A16', 'S100A14', 'S100A13', 'NPR1']
	
	# #TESTING
	# genes = ['S100A14', 'S100A13', 'NPR1','INTS3', 'SLC27A3', 'GATAD2B'] ## WORKS
	# #END TESTING

	# genes = ['INTS3', 'SLC27A3', 'GATAD2B', 'DENND4B', 'CRTC2', 'SLC39A1', 'RAB13', 'RPS27', 'NUP210L', 'TPM3', 'C1orf189', 'UBAP2L', 'HAX1', 'AQP10', 'ATP8B2', 'IL6R', 'UBE2Q1', 'ADAR']
	
	# #TESTING
	# genes = ['IL6R', 'UBE2Q1', 'ADAR','KCNN3', 'PMVK', 'PBXIP1'] ## WORKS
	# #END TESTING

	# genes = ['KCNN3', 'PMVK', 'PBXIP1', 'PYGO2', 'SHC1', 'CKS1B', 'FLAD1', 'LENEP', 'ZBTB7B', 'DCST1', 'ADAM15', 'EFNA4', 'EFNA3', 'EFNA1', 'SLC50A1', 'DPM3', 'KRTCAP2', 'TRIM46', 'MUC1', 'THBS3', 'GBA', 'FAM189B', 'SCAMP3', 'CLK2', 'HCN3', 'FDPS', 'RUSC1', 'ASH1L', 'MSTO1', 'DAP3', 'GON4L', 'SYT11', 'RIT1', 'KHDC4', 'RXFP4', 'ARHGEF2']
	
	# #TESTING
	# genes = ['KHDC4', 'RXFP4', 'ARHGEF2','SSR2', 'UBQLN4', 'LAMTOR2'] ## WORKS
	# #END TESTING

	# genes = ['SSR2', 'UBQLN4', 'LAMTOR2', 'RAB25', 'MEX3A', 'LMNA', 'SEMA4A', 
	# 		'SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 'ARHGEF11', 'ETV3L']
	
	# #TESTING
	# genes = ['PEAR1', 'ARHGEF11', 'ETV3L','ETV3L', 'ETV3', 'FCRL5'] ## WORKS
	# #END TESTING

	# genes = ['ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1', 'CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 'OR10K2',
	# 		'OR10K2', 'OR10K1','OR10R2', 'OR6Y1','OR6P1', 'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 'OR6N1',
	# 		'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 'ACKR1', 'FCER1A', 'OR10J1']
	
	# #TESTING
	# genes = ['ACKR1', 'FCER1A', 'OR10J1','OR10J1', 'OR10J5', 'APCS'] ## WORKS
	# #END TESTING

	# genes = ['OR10J1', 'OR10J5', 'APCS', 'CRP', 'DUSP23', 'FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1', 'CD48']
	
	# #TESTING
	# genes = ['CD84', 'SLAMF1', 'CD48','CD48', 'SLAMF7', 'LY9'] ## WORKS
	# #END TESTING

	# genes = ['CD48', 'SLAMF7', 'LY9', 'CD244', 'ITLN1', 'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1']
	# run(genes, 'downstream')

	#CONITNUE ADDING

	#CAN KEEP GOING
	genes = ['CRCT1', 'LCE3E', 'LCE3D', 'LCE3C', 'LCE3B', 'LCE3A', 'LCE2D', 'LCE2C', 'LCE2B', 'LCE2A', 'LCE4A', 'C1orf68', 'KPRP', 'LCE1F', 'LCE1E', 'LCE1D', 'LCE1C', 'LCE1B', 'LCE1A', 'LCE6A', 'SMCP', 'IVL', 'SPRR4', 'SPRR1A', 'SPRR3', 'SPRR1B', 'SPRR2D', 'SPRR2A', 'SPRR2B', 'SPRR2E', 'SPRR2F', 'SPRR2G', 'LELP1', 'PRR9', 'LORICRIN', 'PGLYRP3', 'PGLYRP4', 'S100A9', 'S100A12', 
	'S100A8', 'S100A7', 'S100A6', 'S100A5', 'S100A4', 'S100A3', 'S100A2', 'S100A16', 'S100A14', 'S100A13', 'NPR1', 
	'INTS3', 'SLC27A3', 'GATAD2B', 'DENND4B', 'CRTC2', 'SLC39A1', 'RAB13', 'RPS27', 'NUP210L', 'TPM3', 'C1orf189', 'UBAP2L', 'HAX1', 'AQP10', 'ATP8B2', 'IL6R', 'UBE2Q1', 'ADAR', 'KCNN3', 'PMVK', 'PBXIP1', 'PYGO2', 'SHC1', 'CKS1B', 'FLAD1', 'LENEP', 'ZBTB7B', 'DCST1', 'ADAM15', 'EFNA4', 'EFNA3', 'EFNA1', 'SLC50A1', 'DPM3', 'KRTCAP2', 'TRIM46', 'MUC1', 'THBS3', 'GBA', 'FAM189B', 'SCAMP3', 'CLK2', 'HCN3', 'FDPS', 'RUSC1','ASH1L', 'MSTO1', 'DAP3', 'GON4L', 'SYT11', 'RIT1', 'KHDC4', 'RXFP4', 'ARHGEF2']
	# genes = ['ARHGEF2', 'SSR2', 'UBQLN4', 'LAMTOR2', 'RAB25', 'MEX3A', 'LMNA', 'SEMA4A', 
			# 'SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1'] 
	# genes = ['KIRREL1', 'CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 'ACKR1', 'FCER1A']
	# genes = ['FCER1A', 'OR10J1', 'OR10J5', 'APCS', 'CRP', 'DUSP23']
	# genes = ['DUSP23', 'FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1']
	# genes = ['CD48', 'SLAMF7', 'LY9', 'CD244', 'ITLN1', 'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1']
	# run(genes, 'upstream')

if __name__ == '__main__':
	main()