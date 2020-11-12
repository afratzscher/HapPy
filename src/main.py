'''
FILE: main.py
PURPOSE: main entry point to application
INPUT: none
OUTPUT: none
'''
import config
import selection
import run
import os
import time

def main():
	# INSERT EMAIL HERE
	config.__EMAIL__ = "example@example.com"
	
	start_time = time.time()
	selection.main()
	run.main()

	print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
	main()
