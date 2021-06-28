import config
import pandas as pd
import json

def getVersion():
	versions = pd.read_csv('GRCh38_chr_versions.txt', sep = '\t')
	config.__CHRVERSION__ = versions.loc[versions['chr'] == str(config.__CHR__)]['version'].values[0]


def getRange():
	df = pd.read_json('data.json')
	minus = df[df['strand'] == 'minus'].index.tolist()
	df.loc[minus,['start','end']] = df.loc[minus,['end','start']].values
			# swaps start/end for genes on minus strand (so that can find region on pos strand and not have error)

	# swap start/end for minus
	if config.__GENENAME__ in df.gene.values:
		idx = df[df['gene'] == config.__GENENAME__].index.tolist()[0]
		if (idx+1 >= len(df)) or (idx-1 < 0): # last or first gene
			print('TO DO: 1st or last gene')
			return(-1)
		
		print(df.loc[idx-5:idx+5])
		print("HERE")
		gene = df.loc[idx]
		df = df.sort_values('start').reset_index()
		nextGene = df.loc[idx+1]
			#sort by start (since doesnt work completely in getGenes.py b/c negative genes have flipped start/end)
		# check if gene overlaps with next gene -> if so, increase idx
		i = 0
		while nextGene['start'] < gene['end']:
			i+=1
			nextGene = df.loc[idx+i]

		df = df.sort_values('end').reset_index()
		prevGene = df.loc[idx-1]
		# check if gene overlaps with previous gene -> if so, decrease idx
		i = 0
		while prevGene['end'] > gene['start']:
			i+=1
			prevGene = df.loc[idx-i]

		geneFlag = (gene['start'] < gene['end'])
		nextFlag = (nextGene['start'] < nextGene['end'])
		prevFlag = (prevGene['start'] < prevGene['end'])
		print(prevFlag, geneFlag, nextFlag, 'start < end', (gene['start'] > prevGene['end']), 'before', (gene['end'] < nextGene['end']), 'after')
		print(prevGene['gene'], prevGene['start'], prevGene['end'])
		print(gene['gene'], gene['start'], gene['end'])
		print(nextGene['gene'], nextGene['start'], nextGene['end'])


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
