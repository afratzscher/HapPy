# tries to build reference from full length haplos
import os
import pandas as pd
from itertools import compress
import itertools
import math
import numpy as np

def fetchDF(gene):
	print(gene)
	file = direct + gene + "/full_length_haplotypes_" + gene + ".vcf"
	file = direct + gene + "/distinct_" + gene + ".vcf"
	try:
		df = pd.read_csv(file, sep="\t")
		if len(df) < 4:
			short.append(gene) 
	except:
		empty.append(gene)

def run(genes):
	for gene in genes:
		fetchDF(gene)
	print(short, ' no full length')
	print(empty, ' MISSING')

def main():
	global direct
	global short
	global empty
	short = []
	empty = []
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/results/DONE/"

	genes = ['FAM72C', 'PPIAL4E', 'NBPF15', 'PPIAL4F', 'SRGAP2B', 'LOC101929805', 'PPIAL4D', 'NBPF20', 'GPR89A', 'ANKRD35', 'ITGA10', 'PEX11B', 'RBM8A', 'LIX1L', 'POLR3GL', 'TXNIP', 'HJV', 'NBPF10', 'NOTCH2NLA', 'PPIAL4H', 'NBPF12', 'CHD1L', 'GJA5', 'GJA8', 'GPR89B', 'NBPF11', 'PPIAL4G', 'NBPF14', 'NOTCH2NLB', 'NUDT4B', 'PDE4DIP', 'NBPF9', 'NOTCH2NLC', 'NBPF19', 'PPIAL4C', 'H2BC18', 'H3C13', 'H4C14', 'H3C14', 'H2AC18', 'H2AC19', 'H3C15', 'H4C15', 'H2BC21', 'H2AC20', 'H2AC21', 'BOLA1', 'SV2A', 'SF3B4', 'OTUD7B', 'VPS45', 'PLEKHO1', 'ANP32E', 'APH1A', 'CIART', 'MRPS21', 'PRPF3', 'RPRD2', 'TARS2', 'ECM1', 'ADAMTSL4', 'MCL1', 'ENSA', 'GOLPH3L', 'HORMAD1', 'CTSS', 'CTSK', 'ARNT', 'CTXND2', 'SETDB1', 'CERS2', 'ANXA9', 'MINDY1', 'PRUNE1', 'BNIPL', 'CDC42SE1', 'MLLT11', 'GABPB2', 'SEMA6C', 'LYSMD1', 'VPS72', 'PIP5K1A', 'PSMD4', 'PI4KB', 'RFX5', 'SELENBP1', 'PSMB4', 'POGZ', 'CGN', 'TUFT1', 'SNX27', 'CELF3', 'RIIAD1', 'TDRKH', 'LINGO4', 'RORC', 'C2CD4D', 'THEM5', 'THEM4', 'S100A10', 'S100A11', 'LOC100131107', 'TCHHL1', 'TCHH', 'RPTN', 'HRNR', 'FLG', 'FLG2', 'CRNN', 'LCE5A', 'CRCT1', 'LCE3E', 'LCE3D', 'LCE3C', 'LCE3B', 'LCE3A', 'LCE2D', 'LCE2C', 'LCE2B', 'LCE2A', 'LCE4A', 'C1orf68', 'KPRP', 'LCE1F', 'LCE1E', 'LCE1D', 'LCE1C', 'LCE1B', 'LCE1A', 'LCE6A', 'SMCP', 'IVL', 'SPRR4', 'SPRR1A', 'SPRR3', 'SPRR1B', 'SPRR2D', 'SPRR2A', 'SPRR2B', 'SPRR2E', 'SPRR2F', 'SPRR2G', 'LELP1', 'PRR9', 'LORICRIN', 'PGLYRP3', 'PGLYRP4', 'S100A9', 'S100A12', 'S100A8', 'S100A7', 'S100A6', 'S100A5', 'S100A4', 'S100A3', 'S100A2', 'S100A16', 'S100A14', 'S100A13', 'NPR1', 'INTS3', 'SLC27A3', 'GATAD2B', 'DENND4B', 'CRTC2', 'SLC39A1', 'RAB13', 'RPS27', 'NUP210L', 'TPM3', 'C1orf189', 'UBAP2L', 'HAX1', 'AQP10', 'ATP8B2', 'IL6R', 'UBE2Q1', 'ADAR', 'KCNN3', 'PMVK', 'PBXIP1', 'PYGO2', 'SHC1', 'CKS1B', 'FLAD1', 'LENEP', 'ZBTB7B', 'DCST1', 'ADAM15', 'EFNA4', 'EFNA3', 'EFNA1', 'SLC50A1', 'DPM3', 'KRTCAP2', 'TRIM46', 'MUC1', 'THBS3', 'GBA', 'FAM189B', 'SCAMP3', 'CLK2', 'HCN3', 'FDPS', 'RUSC1', 'ASH1L', 'MSTO1', 'DAP3', 'GON4L', 'SYT11', 'RIT1', 'KHDC4', 'RXFP4', 'ARHGEF2', 'SSR2', 'UBQLN4', 'LAMTOR2', 'RAB25', 'MEX3A', 'LMNA', 'SEMA4A', 'SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1', 'CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 'ACKR1', 'FCER1A', 'OR10J1', 'OR10J5', 'APCS', 'CRP', 'DUSP23', 'FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1', 'CD48', 'SLAMF7', 'LY9', 'CD244', 'ITLN1', 'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1']
	run(genes)
if __name__ == '__main__':
	main()