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
	names = ['ACKR1', 'AIM2', 'APCS', 'CADM3',
		'FCER1A', 'IFI16', 'OR10J1', 'CRP', 
		'OR10J5', 'PYDC5', 'PYHIN1']
	cmd = makeCommands(names)
	run_commands(*cmd)

if __name__ == '__main__':
	main()
	