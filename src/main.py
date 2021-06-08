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
	start_time = time.time()
	settings.main()
	print("--- settings %s seconds ---" % (time.time() - start_time))

	# print('data notes')
	# print(config.__CHR__)
	# print(config.__START__)
	# print(config.__END__)
	# print(config.__GENESTART__)
	# print(config.__GENEEND__)
	# print(config.__GENENAME__)
	# print(config.__CHRVERSION__)

	# config.__START__ = int('159203314')
	# config.__END__ = int('159283574')
	# config.__GENESTART__ = int('159204875')
	# config.__GENEEND__ = int('159206500')
	# print(config.__FILTERED__)


	selection.main()
	print("--- selection %s seconds ---" % (time.time() - start_time))
	run.main()

	print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
	main()
