'''
FILE: main.py
PURPOSE: main entry point to application
INPUT: none
OUTPUT: none
'''
import config
import options
import run
import os
import time

def main():	
	opts = options.main()
	run.main(opts)

if __name__ == '__main__':
	main()
