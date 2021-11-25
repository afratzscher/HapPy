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
	names = ['RABGAP1L', 'PDC', 'CDC73', 'RHEX']
	
	#CURRENT BLOCK
	names = ['RABGAP1L', 'PDC', 'CDC73',]
	print(len(names))
	cmd = makeCommands(names)
	run_commands(*cmd)

if __name__ == '__main__':
	main()
	