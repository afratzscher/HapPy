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
	direct = os.getcwd()[:-(len('src/getHaplotype/'))] + "results/"
	offset = 10000 # to distinguish 2 df
	threshold = 10000 # 1 million takes too long, 500000 and 100,000 take too long as well
	
	genes = ['H3-2', 'FAM72C', 'PPIAL4E', 'NBPF15', 'PPIAL4F', 'SRGAP2B', 'LOC101929805', 'PPIAL4D', 'NBPF20', 'GPR89A', 'ANKRD35', 'ITGA10', 'PEX11B', 'RBM8A', 'LIX1L', 'POLR3GL', 'TXNIP', 'HJV', 'NBPF10', 'NOTCH2NLA', 'PPIAL4H', 'NBPF12', 'CHD1L', 'GJA5', 'GJA8', 'GPR89B', 'NBPF11', 'PPIAL4G', 'NBPF14', 'NOTCH2NLB', 'NUDT4B', 'PDE4DIP', 'NBPF9', 'NOTCH2NLC', 'NBPF19', 'PPIAL4C', 'H2BC18', 'H3C13', 'H4C14', 'H3C14', 'H2AC18', 'H2AC19', 'H3C15', 'H4C15', 'H2BC21', 'H2AC20', 'H2AC21', 'BOLA1', 'SV2A', 'SF3B4', 'OTUD7B', 'VPS45', 'PLEKHO1', 'ANP32E', 'APH1A', 'CIART', 'MRPS21', 'PRPF3', 'RPRD2', 'TARS2', 'ECM1', 'ADAMTSL4', 'MCL1', 'ENSA', 'GOLPH3L', 'HORMAD1', 'CTSS', 'CTSK', 'ARNT', 'CTXND2', 'SETDB1', 'CERS2', 'ANXA9', 'MINDY1', 'PRUNE1', 'BNIPL', 'CDC42SE1', 'MLLT11', 'GABPB2', 'SEMA6C', 'LYSMD1', 'VPS72', 'PIP5K1A', 'PSMD4', 'PI4KB', 'RFX5', 'SELENBP1', 'PSMB4', 'POGZ', 'CGN', 'TUFT1', 'SNX27', 'CELF3', 'RIIAD1', 'TDRKH', 'LINGO4', 'RORC', 'C2CD4D', 'THEM5', 'THEM4', 'S100A10', 'S100A11', 'LOC100131107', 'TCHHL1', 'TCHH', 'RPTN', 'HRNR', 'FLG', 'FLG2', 'CRNN', 'LCE5A', 'CRCT1', 'LCE3E', 'LCE3D', 'LCE3C', 'LCE3B', 'LCE3A', 'LCE2D', 'LCE2C', 'LCE2B', 'LCE2A', 'LCE4A', 'C1orf68', 'KPRP', 'LCE1F', 'LCE1E', 'LCE1D', 'LCE1C', 'LCE1B', 'LCE1A', 'LCE6A', 'SMCP', 'IVL', 'SPRR4', 'SPRR1A', 'SPRR3', 'SPRR1B', 'SPRR2D', 'SPRR2A', 'SPRR2B', 'SPRR2E', 'SPRR2F', 'SPRR2G', 'LELP1', 'PRR9', 'LORICRIN', 'PGLYRP3', 'PGLYRP4', 'S100A9', 'S100A12', 'S100A8', 'S100A7', 'S100A6', 'S100A5', 'S100A4', 'S100A3', 'S100A2', 'S100A16', 'S100A14', 'S100A13', 'NPR1', 'INTS3', 'SLC27A3', 'GATAD2B', 'DENND4B', 'CRTC2', 'SLC39A1', 'RAB13', 'RPS27', 'NUP210L', 'TPM3', 'C1orf189', 'UBAP2L', 'HAX1', 'AQP10', 'ATP8B2', 'IL6R', 'UBE2Q1', 'ADAR', 'KCNN3', 'PMVK', 'PBXIP1', 'PYGO2', 'SHC1', 'CKS1B', 'FLAD1', 'LENEP', 'ZBTB7B', 'DCST1', 'ADAM15', 'EFNA4', 'EFNA3', 'EFNA1', 'SLC50A1', 'DPM3', 'KRTCAP2', 'TRIM46', 'MUC1', 'THBS3', 'GBA', 'FAM189B', 'SCAMP3', 'CLK2', 'HCN3', 'FDPS', 'RUSC1', 'ASH1L', 'MSTO1', 'DAP3', 'GON4L', 'SYT11', 'RIT1', 'KHDC4', 'RXFP4', 'ARHGEF2', 'SSR2', 'UBQLN4', 'LAMTOR2', 'RAB25', 'MEX3A', 'LMNA', 'SEMA4A', 'SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1', 'CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 'ACKR1', 'FCER1A', 'OR10J1', 'OR10J5', 'APCS', 'CRP', 'DUSP23', 'FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1', 'CD48', 'SLAMF7', 'LY9', 'CD244', 'ITLN1', 'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1', 'MPZ', 'SDHC', 'FCGR2A', 'HSPA6', 'FCGR3A', 'FCGR3B', 'FCGR2B', 'FCRLA', 'FCRLB', 'DUSP12', 'ATF6', 'OLFML2B', 'NOS1AP', 'SPATA46', 'C1orf226', 'SH2D1B', 'UHMK1', 'UAP1', 'DDR2', 'HSD17B7', 'CCDC190', 'RGS4', 'RGS5', 'NUF2', 'PBX1', 'LMX1A', 'RXRG', 'LRRC52', 'MGST3', 'ALDH9A1', 'TMCO1', 'UCK2', 'FAM78B', 'POGK', 'TADA1', 'MAEL', 'GPA33', 'STYXL2', 'POU2F1', 'CD247', 'CREG1', 'RCSD1', 'MPZL1', 'ADCY10', 'DCAF6', 'GPR161', 'TIPRL', 'SFT2D2', 'TBX19', 'XCL2', 'XCL1', 'DPT', 'NME7', 'SLC19A2', 'F5', 'SELP', 'SELL', 'SELE', 'C1orf112', 'KIFAP3', 'NTMT2', 'GORAB', 'PRRX1', 'MROH9', 'FMO3', 'FMO2', 'FMO1', 'FMO4', 'PRRC2C', 'MYOCOS', 'MYOC', 'VAMP4', 'METTL13', 'DNM3', 'C1orf105', 'SUCO', 'FASLG', 'TNFSF18', 'TNFSF4', 'PRDX6', 'SLC9C2', 'KLHL20', 'DARS2', 'ZBTB37', 'SERPINC1', 'RC3H1', 'RABGAP1L', 'CACYBP', 'MRPS14', 'TNN', 'KIAA0040', 'TNR', 'COP1', 'PAPPA2', 'ASTN1', 'BRINP2', 'CRYZL2P-SEC16B', 'TEX35', 'RALGPS2', 'FAM20B', 'TOR3A', 'ABL2', 'SOAT1', 'AXDND1', 'TDRD5', 'FAM163A', 'TOR1AIP2', 'TOR1AIP1', 'CEP350', 'QSOX1', 'ACBD6', 'XPR1', 'KIAA1614', 'STX6', 'MR1', 'IER5', 'CACNA1E', 'ZNF648', 'GLUL', 'TEDDM1', 'RGSL1', 'RNASEL', 'RGS16', 'RGS8', 'NPL', 'DHX9', 'SHCBP1L', 'LAMC1', 'LAMC2', 'NMNAT2', 'SMG7', 'NCF2', 'ARPC5', 'RGL1', 'COLGALT2', 'TSEN15', 'C1orf21', 'EDEM3', 'NIBAN1', 'RNF2', 'SWT1', 'IVNS1ABP', 'HMCN1', 'TPR', 'ODR4', 'PDC', 'PTGS2', 'PLA2G4A', 'BRINP3', 'RGS18', 'RGS21', 'RGS1', 'RGS13', 'RGS2', 'UCHL5', 'GLRX2', 'CDC73', 'KCNT2', 'CFH', 'CFHR3', 'CFHR1', 'CFHR4', 'CFHR2', 'CFHR5', 'F13B', 'ASPM', 'ZBTB41', 'CRB1', 'DENND1B', 'C1orf53', 'LHX9', 'NEK7', 'ATP6V1G3', 'PTPRC', 'NR5A2', 'ZNF281', 'KIF14', 'DDX59', 'CAMSAP2', 'GPR25', 'INAVA', 'KIF21B', 'CACNA1S', 'ASCL5', 'TMEM9', 'IGFN1', 'PKP1', 'TNNT2', 'LAD1', 'TNNI1', 'PHLDA3', 'CSRP1', 'NAV1', 'IPO9', 'SHISA4', 'LMOD1', 'TIMM17A', 'RNPEP', 'ELF3', 'GPR37L1', 'PTPN7', 'LGR6', 'UBE2T', 'PPP1R12B', 'KDM5B', 'RABIF', 'KLHL12', 'ADIPOR1', 'CYB5R1', 'TMEM183A', 'PPFIA4', 'MYOG', 'ADORA1', 'MYBPH', 'CHI3L1', 'CHIT1', 'BTG2', 'FMOD', 'PRELP', 'OPTC', 'ATP2B4', 'LAX1', 'ZC3H11A', 'SNRPE', 'SOX13', 'ETNK2', 'REN', 'KISS1', 'GOLT1A', 'PLEKHA6', 'PPP1R15B', 'PIK3C2B', 'MDM4', 'LRRN2', 'NFASC', 'CNTN2', 'TMEM81', 'RBBP5', 'DSTYK', 'TMCC2', 'NUAK2', 'KLHDC8A', 'LEMD1', 'CDK18', 'MFSD4A', 'ELK4', 'SLC45A3', 'NUCKS1', 'RAB29', 'SLC41A1', 'PM20D1', 'SLC26A9', 'RAB7B', 'CTSE', 'RHEX', 'AVPR1B', 'SRGAP2', 'IKBKE', 'RASSF5', 'DYRK3', 'MAPKAPK2', 'IL19', 'IL20', 'FCMR', 'PIGR', 'FCAMR', 'C1orf116', 'PFKFB2', 'C4BPB', 'C4BPA', 'CD55', 'CR2', 'CR1', 'CR1L', 'CD46', 'CD34', 'PLXNA2', 'CAMK1G', 'LAMB3', 'G0S2', 'HSD11B1', 'TRAF3IP3', 'IRF6', 'UTP25', 'SYT14', 'SERTAD4', 'HHAT', 'KCNH1', 'RCOR3', 'TRAF5', 'RD3', 'SLC30A1', 'NEK2', 'LPGAT1', 'INTS7', 'DTL', 'PPP2R5A', 'PACC1', 'NENF', 'ATF3', 'FAM71A', 'BATF3', 'NSL1', 'SPATA45', 'FLVCR1', 'VASH2', 'ANGEL2', 'RPS6KC1', 'PROX1', 'SMYD2', 'PTPN14', 'CENPF', 'KCNK2', 'KCTD3', 'USH2A', 'ESRRG', 'GPATCH2', 'SPATA17', 'RRP15', 'TGFB2', 'LYPLAL1', 'ZC3H11B', 'SLC30A10', 'EPRS1', 'BPNT1', 'IARS2', 'RAB3GAP2', 'MARK1', 'C1orf115', 'MTARC2', 'MTARC1', 'HLX', 'DUSP10', 'HHIPL2', 'TAF1A', 'MIA3', 'AIDA', 'FAM177B', 'DISP1', 'TLR5', 'SUSD4', 'CCDC185', 'CAPN8', 'CAPN2', 'TP53BP2', 'FBXO28', 'DEGS1', 'NVL', 'CNIH4', 'CNIH3', 'DNAH14', 'LBR', 'ENAH', 'SRP9', 'TMEM63A', 'LEFTY1', 'PYCR2', 'LEFTY2', 'SDE2', 'H3-3A', 'ACBD3', 'MIXL1', 'LIN9', 'PARP1', 'STUM', 'ITPKB', 'PSEN2', 'COQ8A', 'CDC42BPA', 'ZNF678', 'SNAP47', 'PRSS38', 'WNT9A', 'WNT3A', 'ARF1', 'C1orf35', 'MRPL55', 'GUK1', 'GJC2', 'IBA57', 'OBSCN', 'TRIM11', 'TRIM17', 'H3-4', 'H2AW', 'H2BU1', 'RNF187', 'RHOU', 'RAB4A', 'CCSAP', 'ACTA1', 'NUP133', 'ABCB10', 'TAF5L', 'URB2', 'GALNT2', 'PGBD5', 'COG2', 'AGT', 'CAPN9', 'C1orf198', 'TTC13', 'ARV1', 'FAM89A', 'TRIM67', 'C1orf131', 'GNPAT', 'EXOC8', 'SPRTN', 'EGLN1', 'TSNAX', 'DISC1', 'SIPA1L2', 'MAP10', 'PCNX2', 'MAP3K21', 'KCNK1', 'SLC35F3', 'COA6', 'TARBP1', 'IRF2BP2', 'TOMM20', 'RBM34', 'ARID4B', 'GNG4', 'LYST', 'NID1', 'GPR137B', 'EDARADD', 'HEATR1', 'ACTN2', 'MTR', 'MT1HL1', 'RYR2', 'ZP4', 'CHRM3', 'FMN2', 'GREM2', 'RGS7', 'FH', 'KMO', 'WDR64', 'EXO1', 'BECN2', 'MAP1LC3C', 'PLD5', 'AKT3', 'ZBTB18', 'C1orf100', 'CATSPERE', 'DESI2', 'COX20', 'HNRNPU', 'EFCAB2', 'KIF26B', 'SMYD3', 'TFB2M', 'CNST', 'SCCPDH', 'AHCTF1', 'ZNF695', 'ZNF670', 'ZNF669', 'ZNF124', 'ZNF496', 'NLRP3', 'OR2B11', 'GCSAML', 'OR2G2', 'OR2G3', 'OR13G1', 'OR6F1', 'OR14A2', 'OR1C1', 'OR14A16', 'OR11L1', 'TRIM58', 'OR2W3', 'OR2T8', 'OR2AJ1', 'OR2L13', 'OR2M5', 'OR2M2', 'OR2M3', 'OR2M4', 'OR2T33', 'OR2T12', 'OR2M7', 'OR14C36', 'OR2T4', 'OR2T6', 'OR2T1', 'OR2T7', 'OR2T2', 'OR2T3', 'OR2T5', 'OR2G6', 'OR2T29', 'OR2T34', 'OR2T10', 'OR2T11', 'OR2T35', 'OR2T27', 'OR14I1', 'LYPD8', 'SH3BP5L', 'ZNF672', 'ZNF692', 'PGBD2']
	run(genes, 'downstream')


if __name__ == '__main__':
	main()