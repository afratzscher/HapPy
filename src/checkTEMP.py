'''
FILE: haplotype.py
PURPOSE: obtain haplotypes
INPUT: none
OUTPUT: csv file named "haplotypes_1000G_chr1_159203314-159283887.vcf"
CREATED BY: Anne-Sophie Fratzscher
'''

import pandas as pd
import numpy as np
import math
from pathlib import Path
import time


def main():
	readingfile =  '/Users/afratzscher/Documents/GitHub/ACKR1-Algorithm/results/ACKR1/distinct_1000G_ACKR1.vcf'
	newfile =  '/Users/afratzscher/Documents/GitHub/ACKR1-Algorithm/results/ACKR1/new_distinct_ACKR1.vcf'
	olddf = pd.read_csv(readingfile, sep="\t")
	newdf = pd.read_csv(newfile, sep="\t")
	print(olddf.equals(newdf))
	print(olddf.columns.equals(newdf.columns))

if __name__ == '__main__':
	main()
