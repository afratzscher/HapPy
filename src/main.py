'''
FILE: main.py
PURPOSE: main entry point to application
INPUT: none
OUTPUT: none
'''
import config
import settings
import selection
import run
import os
import time

def main():
	# INSERT EMAIL HERE
	# config.__EMAIL__ = "example@example.com"
	config.__EMAIL__ = "afratzscher@yahoo.com"
	
	start_time = time.time()
	settings.main()

	print('data notes')
	print(config.__CHR__)
	print(config.__START__)
	print(config.__END__)
	print(config.__GENESTART__)
	print(config.__GENEEND__)
	print(config.__GENENAME__)


	selection.main()
	run.main()

	print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
	main()
