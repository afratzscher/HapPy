import config
import pandas as pd
import json

def getVersion():
	versions = pd.read_csv('GRCh38_chr_versions.txt', sep = '\t')
	config.__CHRVERSION__ = versions.loc[versions['chr'] == str(config.__CHR__)]['version'].values[0]


def getRange():
	df = pd.read_json('data.json')
	if config.__GENENAME__ in df.gene.values:
		idx = df.index[df['gene'] == config.__GENENAME__]
		if (idx+1 >= len(df)) or (idx-1 < 0): # last or first gene
			print('1st or last gene')
			return(-1)

		gene = df.loc[idx]
		nextGene = df.loc[idx+1]
		prevGene = df.loc[idx-1]

		config.__CHR__ = int(gene['chr'])
		config.__GENESTART__ = int(gene['start'])
		config.__GENEEND__ = int(gene['end'])
		config.__START__ = int(prevGene['end']) + 1
		config.__END__ = int(nextGene['start'])- 1

		getVersion()

		return(0)
	else: # cant find gene
		print('cant find gene')
		return(-1)

def main():
	print(config.__GENENAME__)

	print("*****STARTING GENE SEARCH*****")
	errCode = getRange() #gets range of gene)
	if (errCode == 0):
		return(0)
	else:
		print('Gene not found. Please try again')
		return(-1)
