import pandas as pd
import os
import numpy as np
import ujson
from datetime import date, time 


def getRange(df, first, last, longFlag):
	startlist = []
	endlist = []

	# go through all genes and get range
	for idx in range(0, len(df.index)):
		rangeEnd = 0
		rangeStart = 0

		gene = df.loc[idx]
		df = df.sort_values('start').reset_index(drop=True)
			#sort by start (since doesnt work completely in getGenes.py b/c negative genes have flipped start/end)
		if idx == (len(df.index)-1): 
			rangeEnd = last  #CHANGE LATER!

		else: # if not last gene
			nextGene = df.loc[idx+1]
			# check if gene overlaps with next gene -> if so, increase idx
			i = idx
			while nextGene['start'] < gene['end']:
				i+=1
				if i > len(df.index)-1:
					endFlag = True
					break
				nextGene = df.loc[i]
			rangeEnd = nextGene['start'] -1
		endlist.append(rangeEnd)

		df = df.sort_values('end').reset_index(drop=True)
		if idx == 0: # if not first gene
			rangeStart = first  
		else:
			prevGene = df.loc[idx-1]
			# check if gene overlaps with previous gene -> if so, decrease idx
			i = idx
			while prevGene['end'] > gene['start']:
				i-=1
				if i < 0:
					beginFlag = True
					break
				prevGene = df.loc[i]
			rangeStart = prevGene['end'] +1
		startlist.append(rangeStart)

	df['range_start'] = startlist
	df['range_end'] = endlist
	return df

def toCSV(df, file):
	today = date.today().strftime("%B %d %Y")

	df = df.replace(False, "")
	name = file[:-5] + '.csv'
	f = open(name, "w", newline='')
	f.write("# Updated on: " + today + "\n")
	f.write("# NOTE: short arm = nt 10000-125184587; long arm = nt 143184588-248946422 \n")

	# remove false and replace with None for simplicity when reading for overlapping genes
	df.to_csv(f, index = None)
	f.close()

def main():
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/src"
	file = direct +  "/data.json"
	df = pd.read_json(file)

	# swap minus start/end so that larger number is end
	minus = df[df['strand'] == 'minus'].index.tolist()
	df.loc[minus,['start','end']] = df.loc[minus,['end','start']].values
			# swaps start/end for genes on minus strand (so that can find region on pos strand and not have error)

	#sort by start b/c inverted minus strand genes
	df = df.sort_values('start').reset_index(drop=True)

	# do for short arm first
	shortarm =  df[df['arm'] == 'short']
	shortarm = getRange(shortarm, 10000, 125184587, False)
	file = direct + "/chr1_short.json"
	shortarm.insert(loc=0, column='gene_num', value=range(1, len(shortarm)+1))
	shortarm.to_json(file, orient='records', indent=4)
	toCSV(shortarm, file)

	# then do forst long arm
	longarm =  df[df['arm'] == 'long'].reset_index(drop=True)
	longarm = getRange(longarm, 143184588, 248946422, True)
	file = direct + "/chr1_long.json"
	longarm.insert(loc=0, column='gene_num', value=range(1, len(longarm)+1))
	longarm.to_json(file, orient='records', indent=4)
	toCSV(longarm, file)

	#have file with both arms together
	del longarm['gene_num']
	del shortarm['gene_num']
	df = shortarm.append(longarm)
	df = df.sort_values('start').reset_index(drop=True)

	# check for overlap
	df['start_overlap'] = df['range_start'].duplicated(keep=False) # True if have duplicate, False if no duplicate
	df['end_overlap'] = df['range_end'].duplicated(keep=False)
	# if True for start and end, then is overlapping gene
	df['overlap'] = (df['start_overlap']) | (df['end_overlap'])
	del df['start_overlap']
	del df['end_overlap']
	df.insert(loc=0, column='gene_num', value=range(1, len(df)+1))

	file = direct + "/updated_data.json"
	df.to_json(file, orient='records', indent=4)
	toCSV(df, file)

if __name__ == '__main__':
	main()



