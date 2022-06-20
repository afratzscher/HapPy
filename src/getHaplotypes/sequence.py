'''
FILE: sequence.py
PURPOSE: get sequence around SNP 
INPUT: df
OUTPUT: df with sequence added
CREATED BY: Anne-Sophie Fratzscher
'''

import sys
from Bio import Entrez, SeqIO
import xmltodict
import config
import pandas as pd
from pathlib import Path
from toolz import interleave

def eFetch():
	Entrez.email = config.__EMAIL__
	term = str(config.__CHRVERSION__)

	handle = Entrez.efetch(db="nucleotide", id=term,
		seq_start=config.__START__, seq_stop=config.__END__,
		strand=1, rettype="gb",retmode="text")
	record = SeqIO.read(handle, "genbank")
	handle.close()
	return record

def getFASTA():
	record = eFetch()
	seq = str(record.seq)
	return seq

def replace(df, seq):
	pos = []
	location = []
	other = ['POS', 'ID', 'REF', 'ALT']
	for (index, row) in df.iterrows():
		if index < len(df.index)-1:
			pos.append(int(float(df.iloc[index]['POS'])) - config.__START__)
	pos.append(config.__END__ - config.__START__)

	sequence = []
	for i in range(0, len(pos)):
		if i == 0:
			seqtoadd = seq[0:pos[i]]
			loc = config.__START__
		elif i == (len(pos)-1):
			seqtoadd = seq[(pos[i-1]+1):(pos[i]+1)]
			loc = config.__START__ + pos[i-1]+1
		else:
			seqtoadd = seq[(pos[i-1]+1):pos[i]]
			loc = config.__START__ + pos[i-1]+1
		if seqtoadd == "":
			sequence.append(None)
			location.append(None)
		else:
			sequence.append(seqtoadd)
			location.append(loc)
	
	rowseq = []
	for i in range(0, len(sequence)):
		rowseq.append([sequence[i]] * (len(df.columns)-2))
	
	new = pd.DataFrame(rowseq)
	new.insert(0, 'POS', location)
	new.insert(1, 'ID', None)
	
	new.columns = df.columns
	df = pd.DataFrame(interleave([new.values, df.values]))
	df.columns = new.columns
	return df

def main(df):
	seq = getFASTA()
	df = replace(df, seq)
	return df
