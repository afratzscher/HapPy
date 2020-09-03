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
	# INSERT EMAIL AND REFERENCE GENOME VERSION HERE
	config.__EMAIL__ = 'test@example.com'
	config.__REFVER__ = '38' # INPUT '37' OR '38'
	
	start_time = time.time()
	selection.main()
	run.main()

	print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
	main()
