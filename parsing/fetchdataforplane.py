'''
FILE: getTelomeres.py
PURPOSE: get telomere and centromere sites
NOTE: FASTA starts with index 0, saves assuming index 0 for first nt
INPUT: df
OUTPUT: df with sequence added
CREATED BY: Anne-Sophie Fratzscher
'''

import os
from Bio import Entrez, SeqIO
import pandas as pd


def eFetch(chrom,start,end):
	Entrez.email = "afratzscher@yahoo.com"

	try:
		handle = Entrez.efetch(db="nucleotide", id=chrom,
			seq_start=start , seq_stop= end,
			rettype="fasta",retmode="text")
	except:
		return(None, -1)
	
	record = SeqIO.read(handle, "fasta")
	handle.close()
	return record, 0

def getFASTA(chrom):
	lst = []
	start = 0
	# start = 123400000
	# end = 123500000
	end = 500000 # go in increments of 1 million
	flag = True
	file = open("chr1.txt", 'a')
	while flag:
		# print(start, '-', end, 'vales')
		record, errCode= eFetch(chrom,start,end)
		if errCode == -1: # stops when done
			flag = False
			break
		seq = str(record.seq)
		file.write(seq+"\n")

	
def main():
	print('*****STARTING TELOMERE/CENTROMERE CALCULATION*****')
	name = os.getcwd().rsplit("/",1)[0] + "/src/GRCh38_chr_versions.txt"
	df = pd.read_csv(name, sep = '\t')

	for index,row in df.iterrows():
		chrom = row['version']
		print(chrom)
		seq = getFASTA(chrom)
		break

if __name__ == '__main__':
	main()
