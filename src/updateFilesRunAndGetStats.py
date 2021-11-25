import pandas as pd
import os
from pathlib import Path
import math
import numpy as np

file = 'chr1_long'
direct = os.getcwd() + '/data/'
newfile = direct + 'run_log_' + file + '.json'


#if no file, create
fileCheck = Path(newfile)
if not fileCheck.is_file():
	df = pd.read_json(direct + file + '.json')
	df = df[df['skip']!= True]
	new = df[['gene', 'start', 'end', 'range_start', 'range_end']].copy()
	new['gene_length'] = new['end'] - new['start']

	#change types from object to int64
	new['range_start'] = new['range_start'].astype(str).astype(int)
	new['range_end'] = new['range_end'].astype(str).astype(int)

	new['upstream_inter'] = new['start'] - new['range_start']
	new['upstream_inter'] = new['start'] - new['range_start']
	new.insert(0, 'error', '')
	new.insert(0, 'done', '')
	new['full_haplos'] = ''
	new.to_csv(newfile, sep="\t", mode='a', index = False)

# check if have folder for gene
df = pd.read_csv(newfile, sep="\t")
todo = df[df['done'].isnull()]
# todo = df[df['full_haplos'].isnull()]
# todo = todo.append(df[df['error'] == True])

for index, row in todo.iterrows():
	gene = row['gene']
	print(gene)
	fileCheck = Path('/Volumes/AF_SSD/ACKR1-Algorithm/results/' + gene) 
	if fileCheck.is_dir():
		fileCheck = Path('/Volumes/AF_SSD/ACKR1-Algorithm/results/' + gene + "/meta_" + gene + ".vcf") 
		if (fileCheck.is_file()): # if have meta file, do this block
			#check if have full length haplos
			temp = pd.read_csv('/Volumes/AF_SSD/ACKR1-Algorithm/results/' + gene + "/full_length_haplotypes_" + gene + ".vcf", sep="\t")
			if len(temp) == 3: #only ID, REF, ALT rows = no full length haplos
				df.loc[df['gene'] == gene, 'full_haplos'] = False
			df.loc[df['gene'] == gene, 'done'] = True
			df.loc[df['gene'] == gene, 'error'] = np.nan
		else:
			df.loc[df['gene'] == gene, 'done'] = "error"
			df.loc[df['gene'] == gene, 'error'] = True
		continue

df.to_csv(newfile, sep="\t", index = False)

# check stats left
df = pd.read_csv(newfile, sep="\t")
numTotal = df.shape
print('Total number of genes: ', numTotal[0])
notRun = df[df['done'].isnull()]

torun = df[df['done'].isnull()]
error = df[df['error'] == True]

print('Number done running (no error): ', (numTotal[0] - len(error) - len(torun)))
print('Number to run: ' + str(len(torun)))
print(torun['gene'].tolist())
print('Number of error: ' + str(len(error)))
print(error['gene'].tolist())

print(list(df['gene']))
