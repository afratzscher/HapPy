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

	# TO DO: FIX AQP10, H4C15, H3C15, H2AC19, 
	#ADDED ALL GENES UPSTREM OF ACKR1 ON LONG ARM

	# TO CHECK: H3C15, CTXND2, TUFT1
	names = ['XCL2', 'TIPRL', 'GPR161', 'DCAF6', 'ADCY10', 'MPZL1', 'RCSD1', 'CREG1', 'CD247', 'POU2F1', 'STYXL2', 'GPA33', 'MAEL', 'TADA1', 'POGK', 'FAM78B', 'UCK2', 'TMCO1', 'ALDH9A1', 'MGST3', 'LRRC52', 'RXRG', 'LMX1A', 'PBX1', 'NUF2', 'RGS5', 'RGS4', 'CCDC190', 'HSD17B7', 'DDR2', 'UAP1', 'UHMK1', 'SH2D1B', 'C1orf226', 'SPATA46', 'NOS1AP', 'OLFML2B', 'ATF6', 'DUSP12', 'FCRLB', 'FCRLA', 'FCGR2B', 'FCGR3B', 'FCGR3A', 'HSPA6', 'FCGR2A', 'SDHC', 'MPZ', 'PCP4L1', 'NR1I3', 'APOA2', 'FCER1G', 'ADAMTS4', 'USP21', 'UFC1', 'PFDN2', 'KLHDC9', 'NECTIN4', 'ARHGAP30', 'USF1', 'TSTD1', 'F11R', 'ITLN2', 'ITLN1', 'CD244', 'LY9', 'SLAMF7', 'CD48', 'SLAMF1', 'CD84', 'SLAMF6', 'VANGL2', 'NHLH1', 'NCSTN', 'COPA', 'PEX19', 'DCAF8', 'PEA15', 'CASQ1', 'ATP1A4', 'ATP1A2', 'IGSF8', 'KCNJ9', 'KCNJ10', 'PIGM', 'SLAMF9', 'IGSF9', 'TAGLN2', 'CFAP45', 'VSIG8', 'SLAMF8', 'FCRL6', 'DUSP23']
	names = ['ADCY10']
	cmd = makeCommands(names)
	run_commands(*cmd)

if __name__ == '__main__':
	main()
	