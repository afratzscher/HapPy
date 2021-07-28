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
	names = ['H3-2']
	cmd = makeCommands(names)
	run_commands(*cmd)

if __name__ == '__main__':
	main()
	