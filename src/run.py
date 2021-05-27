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

'''
FUNCTION: computation()
PURPOSE: runs algorithm
INPUT: none
OUTPUT: none
'''
def computation():
	cleaner.main()
	haplotypes.main()
	distinct.main()
	popcounts.main()
	# visualization.main()
	userfile.main()

def main():
	computation()
	