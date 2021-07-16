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
	names = ['OR10R2', 'OR6Y1', 
		'OR6P1', 'OR10X1', 'OR6K2', 'OR6K3', 'OR6K6']
	cmd = makeCommands(names)
	run_commands(*cmd)

if __name__ == '__main__':
	main()
	