import config
import pandas as pd
import json
import os

def getVersion():
	versions = pd.read_csv(os.getcwd() + '/data/GRCh38_chr_versions.txt', sep = '\t')
	config.__CHRVERSION__ = versions.loc[versions['chr'] == str(config.__CHR__)]['version'].values[0]


def getRange():
	for i in range(1,23):
		chrom = i
		df = pd.read_json(os.getcwd() + '/data/chr' + str(chrom) + '.json')

		# swap start/end for minus
		if config.__GENENAME__ in df.gene.values:
			idx = df[df['gene'] == config.__GENENAME__].index.tolist()[0]
			gene = df.loc[idx]

			config.__CHR__ = int(gene['chr'])
			config.__GENESTART__ = int(gene['start'])
			config.__GENEEND__ = int(gene['end'])
			config.__START__ = int(gene['range_start'])
			config.__END__ = int(gene['range_end'])

			getVersion()

			return(0)
	else: # cant find gene
		print('cant find gene')
		return(-1)

def main():

	print("*****STARTING GENE SEARCH*****")
	errCode = getRange() #gets range of gene)
	if (errCode == 0):
		return(0)
	else:
		print('Gene not found. Please try again')
		return(-1)
