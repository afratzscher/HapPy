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
	# TO DO: FIX AQP10, H4C15, H3C15, H2AC19, 
	#ADDED ALL GENES UPSTREM OF ACKR1 ON LONG ARM

	# TO CHECK: H3C15, CTXND2, TUFT1
	# names = ['S100A2', 'S100A4', 'S100A5']
	
	names = ['NBPF12', 'NBPF20', 
	'NOTCH2NLA', 'PEX11B', 'POLR3GL',
	'PPIAL4D', 'PPIAL4H', 'RBM8A',
	'SRGAP2B', 'TXNIP']

	names = ['PPIAL4E']

	print(len(names))
	cmd = makeCommands(names)
	run_commands(*cmd)

if __name__ == '__main__':
	main()
	