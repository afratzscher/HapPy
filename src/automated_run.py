import config
import os
import pandas as pd
from pathlib import Path
import subprocess


def run_commands(*commands):
	os.system(' ; '.join(commands))


# for P3 grch38 b154
def makeCommands(names):
	cmds = []

	for name in names:
		cmd = 'python3 main.py -g ' + name
		cmds.append(cmd)

	return cmds
	
def main():
	names = ['UBAP2L', 'HAX1', 'AQP10', 'ATP8B2', 'IL6R', 
		'UBE2Q1', 'ADAR', 'KCNN3', 'PMVK', 'PBXIP1', 
		'PYGO2', 'SHC1', 'CKS1B', 'FLAD1', 'LENEP', 
		'ZBTB7B', 'DCST1', 'ADAM15', 'EFNA4', 'EFNA3', 
		'EFNA1', 'SLC50A1', 'DPM3', 'KRTCAP2', 'TRIM46', 
		'MUC1', 'THBS3', 'GBA', 'FAM189B', 'SCAMP3', 
		'CLK2', 'HCN3', 'FDPS', 'RUSC1', 'ASH1L', 
		'MSTO1', 'DAP3', 'GON4L', 'SYT11', 'RIT1', 
		'KHDC4', 'RXFP4', 'ARHGEF2', 'SSR2', 'UBQLN4', 
		'LAMTOR2', 'RAB25', 'MEX3A', 'LMNA', 'SEMA4A', 
		'SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 
		'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 
		'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 
		'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 
		'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 
		'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1', 
		'CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 
		'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 
		'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 
		'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 
		'ACKR1', 
		'FCER1A', 'OR10J1', 'OR10J5', 'APCS', 'CRP']

	# TO DO: FIX AQP10
	#ADDED ALL GENES UPSTREM OF ACKR1 ON LONG ARM
	names = ['RPTN', 'TCHH', 'TCHHL1', 'LOC100131107', 'S100A11', 'S100A10', 'THEM4', 'THEM5', 'C2CD4D', 'RORC', 'LINGO4', 'TDRKH', 'RIIAD1', 'CELF3', 'SNX27', 'TUFT1', 'CGN', 'POGZ', 'PSMB4', 'SELENBP1', 'RFX5', 'PI4KB', 'PSMD4', 'PIP5K1A', 'VPS72', 'LYSMD1', 'SEMA6C', 'GABPB2', 'MLLT11', 'CDC42SE1', 'BNIPL', 'PRUNE1', 'MINDY1', 'ANXA9', 'CERS2', 'SETDB1', 'CTXND2', 'ARNT', 'CTSK', 'CTSS', 'HORMAD1', 'GOLPH3L', 'ENSA', 'MCL1', 'ADAMTSL4', 'ECM1', 'TARS2', 'RPRD2', 'PRPF3', 'MRPS21', 'CIART', 'APH1A', 'ANP32E', 'PLEKHO1', 'VPS45', 'OTUD7B', 'SF3B4', 'SV2A', 'BOLA1', 'H2AC21', 'H2AC20', 'H2BC21', 'H4C15', 'H3C15', 'H2AC19', 'H2AC18', 'H3C14', 'H4C14', 'H3C13', 'H2BC18', 'PPIAL4C', 'NBPF19', 'NOTCH2NLC', 'NBPF9', 'PDE4DIP', 'NUDT4B', 'NOTCH2NLB', 'NBPF14', 'PPIAL4G', 'NBPF11', 'GPR89B', 'GJA8', 'GJA5', 'CHD1L', 'NBPF12', 'PPIAL4H', 'NOTCH2NLA', 'NBPF10', 'HJV', 'TXNIP', 'POLR3GL', 'LIX1L', 'RBM8A', 'PEX11B', 'ITGA10', 'ANKRD35', 'GPR89A', 'NBPF20', 'PPIAL4D', 'LOC101929805', 'SRGAP2B', 'PPIAL4F', 'NBPF15', 'PPIAL4E', 'FAM72C', 'H3-2']
	# AQP10
	cmd = makeCommands(names)
	run_commands(*cmd)

if __name__ == '__main__':
	main()
	