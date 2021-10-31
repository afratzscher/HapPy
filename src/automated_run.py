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
	# ISSUES TO RUN
	names = ['NBPF14', 'PDE4DIP', 'H3C15', 'RABGAP1L', 'PDC', 'CDC73', 'RHEX']
	
	#CURRENT BLOCK
	names = ['PCNX2', 'KCNK1', 'SLC35F3', 'COA6', 'TARBP1', 'IRF2BP2', 'TOMM20', 'RBM34', 'ARID4B', 'GNG4', 'LYST', 'NID1', 'GPR137B', 'EDARADD', 'HEATR1', 'ACTN2', 'MTR', 'MT1HL1', 'RYR2', 'ZP4', 'CHRM3', 'FMN2', 'GREM2', 'RGS7', 'FH', 'KMO', 'WDR64', 'EXO1', 'BECN2', 'MAP1LC3C', 'PLD5', 'AKT3', 'ZBTB18', 'C1orf100', 'CATSPERE', 'DESI2', 'COX20', 'HNRNPU', 'EFCAB2', 'KIF26B', 'SMYD3', 'TFB2M', 'CNST', 'SCCPDH', 'AHCTF1', 'ZNF695', 'ZNF670', 'ZNF669', 'ZNF124', 'ZNF496', 'NLRP3', 'OR2B11', 'GCSAML', 'OR2G2', 'OR2G3', 'OR13G1', 'OR6F1', 'OR14A2', 'OR1C1', 'OR14A16', 'OR11L1', 'TRIM58', 'OR2W3', 'OR2T8', 'OR2AJ1', 'OR2L13', 'OR2M5', 'OR2M2', 'OR2M3', 'OR2M4', 'OR2T33', 'OR2T12', 'OR2M7', 'OR14C36', 'OR2T4', 'OR2T6', 'OR2T1', 'OR2T7', 'OR2T2', 'OR2T3', 'OR2T5', 'OR2G6', 'OR2T29', 'OR2T34', 'OR2T10', 'OR2T11', 'OR2T35', 'OR2T27', 'OR14I1', 'LYPD8', 'SH3BP5L', 'ZNF672', 'ZNF692', 'PGBD2', 'RABGAP1L', 'PDC', 'CDC73',]
	print(len(names))
	cmd = makeCommands(names)
	run_commands(*cmd)

if __name__ == '__main__':
	main()
	