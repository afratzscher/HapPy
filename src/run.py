'''
FILE: run.py
PURPOSE: runs computations (after data has been retrieved)
INPUT: none
OUTPUT: none
CREATED BY: Anne-Sophie Fratzscher
'''
import config
import time
import cleaner
import haplotypes
import distinct
import visualization
import popcounts
import userfile
import time

'''
FUNCTION: computation()
PURPOSE: runs algorithm
INPUT: none
OUTPUT: none
'''
def computation():
	start_time = time.time()
	cleaner.main()
	print("--- cleaner %s seconds ---" % (time.time() - start_time))
	start_time = time.time()
	haplotypes.main()
	print("--- haplo %s seconds ---" % (time.time() - start_time))
	start_time = time.time()
	distinct.main()
	print("--- distinct %s seconds ---" % (time.time() - start_time))
	start_time = time.time()
	popcounts.main()
	print("--- popocount %s seconds ---" % (time.time() - start_time))
	start_time = time.time()
	# visualization.main()
	userfile.main()
	print("--- sequence %s seconds ---" % (time.time() - start_time))

def main():
	computation()
	